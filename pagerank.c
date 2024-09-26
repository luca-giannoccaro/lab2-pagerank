#include "xerrori.h"
#include <unistd.h>
#include <math.h>

#define BUF_SIZE 10
#define QUI __LINE__,__FILE__


// definizione tipo inmap
// struct insieme di nodi j con archi di destinazione i
typedef struct inmap{
    int *in;
    int size;
    int messi;
} inmap;

// aggiunge nodo j all'inmap di i
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
        
    }

    return i;
}

// crea inmap
inmap create_inmap () {
    inmap i;

    i.size = 1;
    i.messi = 0;
    i.in = malloc(sizeof(int));
    if (i.in == NULL) xtermina("insufficient memory",QUI);

    return i;
}

// distrugge inmap
void destroy_inmap(inmap i) {
    free(i.in);
    return;
}

int compare(const void *a, const void *b) {
    return (*(int *)a-*(int *)b);
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
            xtermina("arco non valido", QUI);

        // controlla se i=j e in caso ignora l'arco 
        if (out != in) {
            xpthread_mutex_lock(mg,QUI);
            // controlla se i è in g.in[j] e in caso contrario lo aggiunge
            if (!bsearch(&out,g->in[in].in,g->in[in].messi,sizeof(int),compare)) {
                g->in[in] = add_node(g->in[in],out);
                //fprintf(stderr,"added arc %d -> %d\n", out+1, in+1);
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

    while(true) {

        if (*(a->numiter) == a->maxiter) break;

        if (*(a->numiter) == 0) {
            fprintf(stderr, "==%d== Iteration %d\n", gettid(), *(a->numiter));

            for(int i=a->start;i<a->end;i++) { // capire come dividere l'array
                a->x[i] = 1/n;
                // fprintf(stderr,"PR nodo %d all'iterazione %d: %.4f\n", i, *(a->numiter), a->x[i]);

                xpthread_mutex_lock(a->mutex,QUI);
                *(a->x_done) += 1;

                if ((*(a->x_done) % g->N) == 0) {
                    *(a->numiter) += 1;
                    fprintf(stderr, "==%d== Array filled.\n", gettid());
                    xpthread_cond_broadcast(a->cv_x,QUI);
                }
                xpthread_mutex_unlock(a->mutex,QUI);
            }

            // attendo che gli altri thread finiscano di riempire l'array
            xpthread_mutex_lock(a->mutex,QUI);
            while ((*(a->x_done) % g->N) > 0) {
                fprintf(stderr, "==%d== I'm waiting for the other threads to finish filling the array.\n", gettid());
                xpthread_cond_wait(a->cv_x,a->mutex,QUI);
            }
            xpthread_mutex_unlock(a->mutex,QUI);
 
        }
        else {
            fprintf(stderr, "==%d== Iteration %d\n", gettid(), *(a->numiter));

            for (int i=a->start;i<a->end;i++) {
                // calcolo y e s
                if (g->out[i] > 0){
                    a->y[i] = a->x[i] / (double)g->out[i];
                    // fprintf(stderr, "==%d== Y(t) nodo %d: %.4f\n", gettid(), i+1, a->y[i]);
                    }
                else {
                    // fprintf(stderr, "X(t) di %d: %.4f\n", i+1, a->x[i]);
                    *(a->s) += a->x[i];
                }
                xpthread_mutex_lock(a->mutex,QUI);
                *(a->y_done) += 1;

                if ((*(a->y_done) % g->N) == 0) {
                    fprintf(stderr, "Yt and St calculated for every node.\n");
                    // fprintf(stderr, "St = %.4f\n", *(a->s));
                    xpthread_cond_broadcast(a->cv_y,QUI);
                }
                xpthread_mutex_unlock(a->mutex,QUI);
            }   
            
            // attendo che gli altri thread calcolo y e s
            xpthread_mutex_lock(a->mutex,QUI);
            while ((*(a->y_done) % g->N) > 0) {
                fprintf(stderr, "==%d== I'm waiting for the other threads to finish calculating Yt and St.\n", gettid());
                xpthread_cond_wait(a->cv_y,a->mutex,QUI);
            }
            xpthread_mutex_unlock(a->mutex,QUI);

            fprintf(stderr, "==%d== I restarted. Calculating X(t+1)...\n", gettid());

            for (int i=a->start;i<a->end;i++) {
                double sumy = 0;
                for (int j=0; j<g->N; j++) {
                    if (bsearch(&j,g->in[i].in,g->in[i].messi,sizeof(int),compare)){
                        sumy += a->y[j];
                    }
                }
                //fprintf(stderr, "Sum of Y(t) for every node in IN(%d): %.4f\n", i+1, sumy);

                a->xnext[i] = a->tp + a->d * sumy + (a->d/n) * (*(a->s));
                // fprintf(stderr,"==%d== PR nodo %d all'iterazione %d: %.4f\n", gettid(), i+1, *(a->numiter), a->xnext[i]);

                *(a->e) += fabs(a->xnext[i] - a->x[i]);
                a->x[i] = a->xnext[i];

                xpthread_mutex_lock(a->mutex,QUI);
                *(a->xnext_done) += 1;

                if ((*(a->xnext_done) % g->N) == 0) {
                    fprintf(stderr, "Errore all'iterazione %d: %.8f\n", *(a->numiter), *(a->e));
                    if (*(a->e) < a->eps) {
                        *(a->stop) = true;
                        xpthread_cond_broadcast(a->cv_xnext,QUI);
                        xpthread_mutex_unlock(a->mutex,QUI);
                        break;
                    }
                    *(a->numiter) += 1;
                    *(a->s) = 0;
                    *(a->e) = 0;
                    fprintf(stderr, "X(t+1) calculated for every node.\n");
                    xpthread_cond_broadcast(a->cv_xnext,QUI);
                }
                xpthread_mutex_unlock(a->mutex,QUI);               
            }

            xpthread_mutex_lock(a->mutex, QUI);

                while ((*(a->xnext_done) % g->N) > 0 && !(*(a->stop))){
                    fprintf(stderr, "==%d== I'm waiting for the other threads to finish calculating X(t+1).\n", gettid());
                    xpthread_cond_wait(a->cv_xnext,a->mutex,QUI);
                }

                if (*(a->stop)) {
                    xpthread_mutex_unlock(a->mutex,QUI);
                    break;
                }

            fprintf(stderr, "==%d== I restarted. Moving to iteration %d...\n", gettid(), *(a->numiter));
            xpthread_mutex_unlock(a->mutex,QUI);
        }
    }

    return (void *) 0;
}


// funzione calcolo pagerank
double *pagerank(grafo *g, double d, double eps, int maxiter, int taux, int *numiter) {
    // thread gestione segnali

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
        // xpthread_mutex_lock(&m,QUI);
        buffer[pindex % BUF_SIZE] = arc;
        // fprintf(stderr, "arc in buffer: %d -> %d\n", buffer[pindex % BUF_SIZE].out,buffer[pindex % BUF_SIZE].in);
        pindex += 1;
        // xpthread_mutex_unlock(&m,QUI);
        xsem_post(&sem_data_items,QUI);
    }
    if(fclose(f)!=0) xtermina("Input file closing error",QUI);

    // inserisce valori di terminazione nel buffer
    for (int i=0; i<taux; i++) {
        arc.out = -1;
        arc.in = -1;
        xsem_wait(&sem_free_slots,QUI);
        buffer[pindex % BUF_SIZE] = arc;
        // fprintf(stderr, "arc in buffer: %d -> %d\n", buffer[pindex % BUF_SIZE].out,buffer[pindex % BUF_SIZE].in);
        pindex += 1;
        xsem_post(&sem_data_items,QUI);

    }

    // join thread lettura dati
    for (int i=0; i<taux; i++) {
        xpthread_join(t[i],NULL,QUI);
    }

    fprintf(stderr, "Finished reading and generating graph.\n");

    // count dead-end nodes and valid arcs
    int deadend = 0;
    int valid = 0;
    for (int i=0; i<g.N; i++) {
        if (g.out[i] == 0) deadend++;
        valid += g.out[i];
    }

    fprintf(stderr, "Calculating Pagerank...\n");
    // chiama pagerank()
    int numiter = 0;
    double *pr = pagerank(&g,d,eps,maxiter,taux,&numiter);
    fprintf(stderr,"Pagerank calculated.\n");

    // calcola somma dei pr
    double prsum = 0;
    for (int i=0; i<g.N; i++) {
        prsum += pr[i];
        fprintf(stderr, "PR nodo %d: %.6f\n", i+1, pr[i]);
    }

    // find top k nodes


    // print risultati
    printf("Number of nodes: %d\n", g.N);
    printf("Number of dead-end nodes: %d\n", deadend);
    printf("Number of valid arcs: %d\n", valid);
    if (numiter < maxiter)
        printf("Converged after %d iterations\n", numiter);
    else
        printf("Did not converge after %d iterations\n", maxiter);
    printf("Sum of ranks: %.4f (should be 1)\n", prsum);
    /*
    printf("Top %d nodes:\n", k);
    for (int i=0, i<k, i++) {
        printf("\t%d %.6f", );
    }*/

    // dealloc
    xsem_destroy(&sem_data_items,QUI);
    xsem_destroy(&sem_free_slots,QUI);
    xpthread_mutex_destroy(&m,QUI);
    xpthread_mutex_destroy(&mg,QUI);
    free(g.out);
    for (int i=0; i<g.N; i++) {
        destroy_inmap(g.in[i]);
    }
    free(g.in);
    free(pr);

    return 0;
}