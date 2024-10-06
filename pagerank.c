#include "xerrori.h"

#define BUF_SIZE 10
#define QUI __LINE__,__FILE__

// dati thread gestione segnali
int max_node = -1;
double max_pr = 0.0;

typedef struct {
    bool *prstarted;
    int *numiter;
    int *max_node;
    double *max_pr;
    pthread_mutex_t *mutex;
} signal_data;

void *tsignal(void *arg) {
    signal_data *d = (signal_data *) arg;
    pthread_mutex_t *m = d->mutex;
    
    sigset_t mask;
    sigemptyset(&mask);
    sigaddset(&mask,SIGUSR1);
    sigaddset(&mask,SIGUSR2);

    int s;
    while (true) {

        int e = sigwait(&mask,&s);
        if(e!=0) xtermina("Errore sigwait", QUI);

        if (s==SIGUSR2) break;

        if (!*(d->prstarted)) fprintf(stderr, "Calcolo Pagerank non ancora iniziato.\n");
        else {
            xpthread_mutex_lock(m,QUI);
            fprintf(stderr, "Iterazione corrente: %d\n", *(d->numiter));
            fprintf(stderr, "Nodo con massimo Pagerank: %d (%f)\n", *(d->max_node), *(d->max_pr));
            xpthread_mutex_unlock(m,QUI);
        }
    }
    return (void *) 0;
}

// definizione tipo inmap
// insieme di nodi j con archi di destinazione i
typedef struct inmap{
    int *in;
    int size;
    int messi;
} inmap;

// aggiunge nodo j all'inmap i
inmap add_node(inmap i, int j) {

    if (i.messi == i.size) {
        i.size *= 2;
        i.in = realloc(i.in, sizeof(int)*i.size);
        if (i.in==NULL) xtermina("realloc failed",QUI);
    }

    if (i.messi == 0 || i.in[i.messi-1] < j) {
        i.in[i.messi] = j;
        
        i.messi++;
    }
    else {
        int n = 0;
        while (j>i.in[n]) n++;
        memmove(&i.in[n+1], &i.in[n], (i.messi - n)*sizeof(int));
        i.in[n] = j;
        i.messi++;
    }

    return i;
}

// crea inmap
inmap create_inmap () {
    inmap i;

    i.size = 8;
    i.messi = 0;
    i.in = malloc(sizeof(int)*i.size);
    if (i.in == NULL) xtermina("insufficient memory",QUI);

    return i;
}

// distrugge inmap
void destroy_inmap(inmap i) {
    free(i.in);
    return;
}

// funzione di confronto per bsearch
int compare(const void *a, const void *b) {
    return (*(int *)a-*(int *)b);
}

// trova gli n nodi con il pr più alto
double *find_max(int n, double *pr, int size, int *nodes) {
    double *max = calloc(n, sizeof(double));
    for (int i=0;i<n;i++) {
        for (int j=0; j<size; j++) {
            if (pr[j]>max[i]) {
                if (i==0) {
                    max[i] = pr[j];
                    nodes[i] = j;
                }
                else if (pr[j]<=max[i-1] && j!=nodes[i-1]) {
                    max[i] = pr[j];
                    nodes[i] = j;
                }
            }
        }
    }
    return max;
}

// struct grafo
typedef struct {
    int N; // numero nodi grafo
    int *out; // array con il numero di archi uscenti da ogni nodo
    inmap *in; // array con gli insiemi di archi entranti in ogni nodo
} grafo;

// struct per memorizzare archi letti da file
typedef struct {
    int out;
    int in;
} arco;

// struct dati thread lettura
typedef struct {
    grafo *grafo;
    arco *buffer; 
    int *pcindex;
    pthread_mutex_t *mutex;
    pthread_mutex_t *mutex_graph;
    sem_t *sem_free_slots;
    sem_t *sem_data_items;
} tdati;

// funzione thread lettura dati
void *tread(void *arg) {
    // legge archi dal buffer e memorizza i dati in un grafo'''

    tdati *a = (tdati *)arg;
    pthread_mutex_t *m = a->mutex;
    pthread_mutex_t *mg = a->mutex_graph;
    grafo *g = a->grafo;
    arco arc;

    do {
        // legge arco dal buffer
        xsem_wait(a->sem_data_items,QUI);
        xpthread_mutex_lock(m,QUI);
    
        arc = a->buffer[*(a->pcindex) % BUF_SIZE];
        *(a->pcindex)+=1;
        xpthread_mutex_unlock(m,QUI);
        xsem_post(a->sem_free_slots,QUI);

        if (arc.out == -1 && arc.in == -1) break;

        int out = arc.out - 1;
        int in = arc.in - 1;

        if (out < 0 || out > g->N-1 || in < 0 || in > g->N-1)
            xtermina("Arco non valido.", QUI);

        // controlla se i=j e in caso ignora l'arco 
        if (out != in) {
            xpthread_mutex_lock(mg,QUI);
            // controlla se i è in g.in[j] e in caso contrario lo aggiunge
            if (!bsearch(&out,g->in[in].in,g->in[in].messi,sizeof(int),compare)) {
                g->in[in] = add_node(g->in[in],out);
                // fprintf(stderr,"added arc %d -> %d\n", out+1, in+1);
                // fprintf(stderr, "archi aggiunti a inmap di %d: %d\n", in+1, g->in[in].messi);
                //incremenet g.out
                g->out[out]++;
            }
            xpthread_mutex_unlock(mg,QUI);
        }

    } while (true);

    return (void *)0;
}

// struct dati thread calcolo
typedef struct {
    int start;
    int end;
    grafo *grafo;
    double tp;
    double d;
    double eps;
    int maxiter;
    int *numiter;
    double *x;
    double *y;
    double *xnext;
    double *s;
    double *e;
    int *x_done;
    int *y_done;
    int *xnext_done;
    bool *stop;
    double *max;
    int *node;
    pthread_cond_t *cv_x;
    pthread_cond_t *cv_y;
    pthread_cond_t *cv_xnext;
    pthread_mutex_t *mutex;
} tcalcolo;

// funzione thread calcolo
void *tpagerank(void *arg) {
    tcalcolo *a = (tcalcolo *)arg;
    grafo *g = a->grafo;
    double n = (double)g->N;
    double *max = a->max;
    int *node = a->node;

    while(true) {

        if (*(a->numiter) == a->maxiter) break;

        if (*(a->numiter) == 0) {

            for(int i=a->start;i<a->end;i++) { // capire come dividere l'array
                a->x[i] = 1/n;
                
                // fprintf(stderr,"PR nodo %d all'iterazione %d: %.10f\n", i, *(a->numiter), a->x[i]);

                xpthread_mutex_lock(a->mutex,QUI);
                *(a->x_done) += 1;

                 if (a->x[i]>*max) {
                    *max = a->x[i];
                    *node = i;
                }

                if ((*(a->x_done) % g->N) == 0) {
                    *(a->numiter) += 1;
                    // fprintf(stderr, "==%d== Array filled.\n", gettid());
                    max_pr = *max;
                    max_node = *node;
                    *max = 0.0;
                    *node = -1;
                    xpthread_cond_broadcast(a->cv_x,QUI);
                }
                xpthread_mutex_unlock(a->mutex,QUI);
            }

            // attendo che gli altri thread finiscano di riempire l'array
            xpthread_mutex_lock(a->mutex,QUI);
            while ((*(a->x_done) % g->N) > 0) {
                // fprintf(stderr, "==%d== I'm waiting for the other threads to finish filling the array.\n", gettid());
                xpthread_cond_wait(a->cv_x,a->mutex,QUI);
            }
            xpthread_mutex_unlock(a->mutex,QUI);
 
        }
        else {

            for (int i=a->start;i<a->end;i++) {
                // calcolo y e s
                if (g->out[i] > 0){
                    a->y[i] = a->x[i] / (double)g->out[i];
                    // fprintf(stderr, "==%d== Y(t) nodo %d: %.4f\n", gettid(), i+1, a->y[i]);
                    }
                else {
                    // fprintf(stderr, "X(t) di %d: %.4f\n", i+1, a->x[i]);
                    xpthread_mutex_lock(a->mutex,QUI);
                    *(a->s) += a->x[i];
                    xpthread_mutex_unlock(a->mutex,QUI);
                }
                xpthread_mutex_lock(a->mutex,QUI);
                *(a->y_done) += 1;

                if ((*(a->y_done) % g->N) == 0) {
                    // fprintf(stderr, "==%d== Yt and St calculated for every node.\n", gettid());
                    // fprintf(stderr, "St = %.4f\n", *(a->s));

                    xpthread_cond_broadcast(a->cv_y,QUI);
                }
                xpthread_mutex_unlock(a->mutex,QUI);
            }   
            
            // attendo che gli altri thread calcolino y e s
            xpthread_mutex_lock(a->mutex,QUI);
            while ((*(a->y_done) % g->N) > 0) {
                // fprintf(stderr, "==%d== I'm waiting for the other threads to finish calculating Yt and St.\n", gettid());
                xpthread_cond_wait(a->cv_y,a->mutex,QUI);
            }
            xpthread_mutex_unlock(a->mutex,QUI);

            // fprintf(stderr, "==%d== I restarted. Calculating X(t+1)...\n", gettid());

            for (int i=a->start;i<a->end;i++) {
                double sumy = 0;
                for (int j=0; j<g->in[i].messi; j++) {
                        sumy += a->y[g->in[i].in[j]];
                }
                //fprintf(stderr, "Sum of Y(t) for every node in IN(%d): %.4f\n", i+1, sumy);

                a->xnext[i] = a->tp + a->d * sumy + (a->d/n) * (*(a->s));
                
                // fprintf(stderr,"==%d== PR nodo %d all'iterazione %d: %.4f\n", gettid(), i+1, *(a->numiter), a->xnext[i]);
                xpthread_mutex_lock(a->mutex,QUI);
                *(a->e) += fabs(a->xnext[i] - a->x[i]);
                xpthread_mutex_unlock(a->mutex,QUI);

                a->x[i] = a->xnext[i];

                xpthread_mutex_lock(a->mutex,QUI);
                *(a->xnext_done) += 1;

                if (a->xnext[i]>*max) {
                    *max = a->xnext[i];
                    *node = i;
                }

                if ((*(a->xnext_done) % g->N) == 0) {
                    // fprintf(stderr, "Errore all'iterazione %d: %.8f\n", *(a->numiter), *(a->e));
                    if (*(a->e) < a->eps) {
                        *(a->stop) = true;
                        xpthread_cond_broadcast(a->cv_xnext,QUI);
                        xpthread_mutex_unlock(a->mutex,QUI);
                        break;
                    }
                    *(a->numiter) += 1;

                    // resetta i valori di St ed e per la prossima iterazione
                    *(a->s) = 0;
                    *(a->e) = 0;
                    max_pr = *max;
                    max_node = *node;
                    *max = 0.0;
                    *node = -1;

                    // fprintf(stderr, "==%d== X(t+1) calculated for every node.\n", gettid());
                    xpthread_cond_broadcast(a->cv_xnext,QUI);
                }
                xpthread_mutex_unlock(a->mutex,QUI);               
            }

            xpthread_mutex_lock(a->mutex, QUI);

                while ((*(a->xnext_done) % g->N) > 0 && !(*(a->stop))){
                    // fprintf(stderr, "==%d== I'm waiting for the other threads to finish calculating X(t+1).\n", gettid());
                    xpthread_cond_wait(a->cv_xnext,a->mutex,QUI);
                }

                if (*(a->stop)) {
                    xpthread_mutex_unlock(a->mutex,QUI);
                    break;
                }

            // fprintf(stderr, "==%d== I restarted. Moving to the next iteration...\n", gettid());
            xpthread_mutex_unlock(a->mutex,QUI);
        }
    }

    return (void *) 0;
}

// funzione calcolo pagerank
double *pagerank(grafo *g, double d, double eps, int maxiter, int taux, int *numiter) {

    // alloca 3 vettori di g.N double
    double *x = calloc(g->N,sizeof(double));
    if (x == NULL) xtermina ("insufficient memory", QUI);
    double *y = calloc(g->N,sizeof(double));
    if (y == NULL) xtermina ("insufficient memory", QUI);
    double *xnext = calloc(g->N,sizeof(double));
    if (xnext == NULL) xtermina ("insufficient memory", QUI);

    double tp = (1-d)/(double)g->N; // teleporting

    // inizializza dati thread calcolo
    pthread_t t[taux];
    tcalcolo a[taux];
    double s = 0, e = 0;
    int x_done = 0, y_done = 0, xnext_done = 0;
    bool stop = false;
    double max = 0.0;
    int node = -1;
    pthread_cond_t cv_x = PTHREAD_COND_INITIALIZER;
    pthread_cond_t cv_y = PTHREAD_COND_INITIALIZER;
    pthread_cond_t cv_xnext = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    
    // inizializza thread calcolo
    for (int i=0; i<taux; i++) {
        int n = g->N/taux;  
        a[i].start = n*i; 
        a[i].end = (i==taux-1) ? g->N : n*(i+1);
        a[i].grafo = g;
        a[i].tp = tp;
        a[i].d = d;
        a[i].eps = eps;
        a[i].maxiter = maxiter;
        a[i].numiter = numiter;
        a[i].x = x;
        a[i].y = y;
        a[i].xnext = xnext;
        a[i].s = &s;
        a[i].e = &e;
        a[i].x_done = &x_done;
        a[i].y_done = &y_done;
        a[i].xnext_done = &xnext_done;
        a[i].stop = &stop;
        a[i].max = &max;
        a[i].node = &node;
        a[i].cv_x = &cv_x;
        a[i].cv_y = &cv_y;
        a[i].cv_xnext = &cv_xnext;
        a[i].mutex = &m;
        xpthread_create(&t[i],NULL,&tpagerank,&a[i],QUI);
    }

    // join thread
    for (int i=0; i<taux; i++) {
        xpthread_join(t[i],NULL,QUI);
    }

    free(x);
    free(y);
    xpthread_cond_destroy(&cv_x,QUI);
    xpthread_cond_destroy(&cv_y,QUI);
    xpthread_cond_destroy(&cv_xnext,QUI);
    xpthread_mutex_destroy(&m,QUI);

    return xnext; // restituisce vettore di tutti i pr
}

int main(int argc, char *argv[]) {
    time_t s1;
    time(&s1);

    // dati gestore segnali
    pthread_mutex_t sigmutex = PTHREAD_MUTEX_INITIALIZER;
    int numiter = 0;
    bool prstarted = false;

    signal_data data = {
        .prstarted = &prstarted,
        .numiter = &numiter,
        .max_node = &max_node,
        .max_pr = &max_pr,
        .mutex = &sigmutex
    };

    sigset_t mask;
    sigfillset(&mask);
    sigdelset(&mask,SIGQUIT);
    pthread_sigmask(SIG_BLOCK,&mask,NULL);

    pthread_t handler;
    xpthread_create(&handler,NULL,&tsignal,&data,QUI);


    // leggi input
    if (argc < 2) {
        fprintf(stderr, "Use\n\t%s infile", argv[0]);
        exit(1);
    }

    // leggi argomenti opzionali da linea di comando
    extern char *optarg;
    int k = 3;
    int maxiter = 100;
    int taux = 3;
    double d = 0.9;
    double eps = 1.0e-7;

    int c = 0;
    while((c = getopt(argc, argv, ":k:m:d:e:t:"))!=-1) {
        switch (c) {
            case 'k':
                k = atoi(optarg);
                break;
            case 'm':
                maxiter = atoi(optarg);
                break;
            case 'd':
                d = atof(optarg);
                break;
            case 'e':
                eps = atof(optarg);
                break;
            case 't':
                taux = atoi(optarg);
                break;
            case ':':
                break;
            case '?':
                fprintf(stderr, "Unknown option: -%c.\n", optopt);
                break;
        }
    }

    arco arc;
    grafo g;

    // apre file e legge prima linea 
    FILE *f = fopen(argv[argc-1], "r");
    if (f==NULL) {
        xtermina("File opening error", QUI);
    }

    // ignora i commenti
    while (true) {
        if (getc(f)=='%') {
            while (getc(f) != '\n') {
                continue;
            }
        }
        else break;
    }

    fseek(f, -1, SEEK_CUR);
    int n, tmp1, tmp2;
    int e = fscanf(f, "%d %d %d", &n, &tmp1, &tmp2);
    if (e!=3) xtermina("fscanf error", QUI);

    fprintf(stderr, "Reading mtx file and generating graph...\n");

    g.N = n;

    g.out = calloc(n,sizeof(int));
    if (g.out == NULL) xtermina("insufficient memory",QUI);

    g.in = malloc(sizeof(inmap)*n);
    if (g.in == NULL) xtermina("insufficient memory",QUI);
    for (int i=0;i<n;i++) {
        g.in[i] = create_inmap();
    }

    // inizializza thread lettura dati
    pthread_t t[taux];
    tdati a[taux];
    arco buffer[BUF_SIZE];
    int pindex=0, cindex=0;
    pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t mg = PTHREAD_MUTEX_INITIALIZER;
    sem_t sem_free_slots, sem_data_items;
    xsem_init(&sem_free_slots,0,BUF_SIZE,QUI);
    xsem_init(&sem_data_items,0,0,QUI);

    for (int i=0; i<taux;i++) {
        a[i].grafo = &g;
        a[i].buffer = buffer;
        a[i].pcindex = &cindex;
        a[i].mutex = &m;
        a[i].mutex_graph = &mg;
        a[i].sem_data_items = &sem_data_items;
        a[i].sem_free_slots = &sem_free_slots;
        xpthread_create(&t[i],NULL,&tread,&a[i],QUI);
    }

    // scrive sul buffer gli archi letti dal file
    while (true) {
        e = fscanf(f, "%d %d", &arc.out, &arc.in);
        if (e!=2) {
            break;
        }
        xsem_wait(&sem_free_slots,QUI);
        xpthread_mutex_lock(&m,QUI);
        buffer[pindex % BUF_SIZE] = arc;
        // fprintf(stderr, "arc in buffer: %d -> %d\n", buffer[pindex % BUF_SIZE].out,buffer[pindex % BUF_SIZE].in);
        pindex += 1;
        xpthread_mutex_unlock(&m,QUI);
        xsem_post(&sem_data_items,QUI);
    }
    if(fclose(f)!=0) xtermina("Input file closing error",QUI);

    // inserisce valori di terminazione nel buffer
    for (int i=0; i<taux; i++) {
        arc.out = -1;
        arc.in = -1;
        xsem_wait(&sem_free_slots,QUI);
        xpthread_mutex_lock(&m,QUI);
        buffer[pindex % BUF_SIZE] = arc;
        // fprintf(stderr, "arc in buffer: %d -> %d\n", buffer[pindex % BUF_SIZE].out,buffer[pindex % BUF_SIZE].in);
        pindex += 1;
        xpthread_mutex_unlock(&m,QUI);
        xsem_post(&sem_data_items,QUI);
    }

    // join thread lettura dati
    for (int i=0; i<taux; i++) {
        xpthread_join(t[i],NULL,QUI);
    }

    fprintf(stderr, "Finished reading and generating graph.\n");

    // dealloca
    xsem_destroy(&sem_data_items,QUI);
    xsem_destroy(&sem_free_slots,QUI);
    xpthread_mutex_destroy(&m,QUI);
    xpthread_mutex_destroy(&mg,QUI);

    // count dead-end nodes and valid arcs
    int deadend = 0;
    int valid = 0;
    for (int i=0; i<g.N; i++) {
        if (g.out[i] == 0) deadend++;
        if (g.in[i].size > g.in[i].messi) {
            g.in[i].in = realloc(g.in[i].in, sizeof(int)*g.in[i].messi);
            if (g.in[i].messi != 0 && g.in[i].in == NULL) xtermina("realloc fallita", QUI);
        }
        valid += g.out[i];
    }

    fprintf(stderr, "Calculating Pagerank...\n");
    // chiama pagerank()
    prstarted = true;
    double *pr = pagerank(&g,d,eps,maxiter,taux,&numiter);
    fprintf(stderr,"Pagerank calculated.\n");

    // calcola somma dei pr
    double prsum = 0;
    for (int i=0; i<g.N; i++) {
        prsum += pr[i];
    }


    // trova i k nodi con il pr più alto
    double *max = malloc(sizeof(double)*k);
    int *nodes = malloc(sizeof(int)*k);
    max = find_max(k, pr, g.N, nodes);
    

    // print risultati
    printf("Number of nodes: %d\n", g.N);
    printf("Number of dead-end nodes: %d\n", deadend);
    printf("Number of valid arcs: %d\n", valid);
    if (numiter < maxiter)
        printf("Converged after %d iterations\n", numiter);
    else
        printf("Did not converge after %d iterations\n", maxiter);
    printf("Sum of ranks: %.4f (should be 1)\n", prsum);
    printf("Top %d nodes:\n", k);
    for (int i=0; i<k; i++) {
        printf("\t%d %.6f\n", nodes[i], max[i]);
    }

    // dice all'handler di fermarsi e fa la join
    pthread_kill(handler, SIGUSR2);
    xpthread_join(handler,NULL,QUI);

    // dealloca
    free(g.out);
    for (int i=0; i<g.N; i++) {
        destroy_inmap(g.in[i]);
    }
    free(g.in);
    free(pr);
    free(max);
    free(nodes);
    xpthread_mutex_destroy(&sigmutex,QUI);

    time_t s2;
    time(&s2);

    fprintf(stderr, "Run time: %ld m %ld s\n", (s2-s1)/60, (s2-s1)%60);

    return 0;
}
