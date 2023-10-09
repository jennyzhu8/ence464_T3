#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>


/**
 * poisson.c
 * Implementation of a Poisson solver with Dirichlet boundary conditions.
 *
 * This template handles the basic program launch, argument parsing, and memory
 * allocation required to implement the solver *at its most basic level*. You
 * will likely need to allocate more memory, add threading support, account for
 * cache locality, etc...
 *
 * BUILDING:
 * gcc -o poisson poisson.c -lpthread
 *
 * [note: linking pthread isn't strictly needed until you add your
 *        multithreading code]
 *
 * 
 * Run in terminal by entering the folder,
 * $ make (this just compiles it with the default settings in the makefile)
 * A poisson executable will appear (without an extension .c)
 * $ ./poisson to run the executable
 * $ ./poisson -n 101 -i 100 to choose iterations and cube size
 * Change the optimiser flags in the make file
 * $ make clean resolves everything, need to make again after this (use if change header files)
 * 
 * 
 * TODO:
 * 1 - Read through this example, understand what it does and what it gives you
 *     to work with.
 * 2 - Implement the basic algorithm and get a correct output.
 * 3 - Add a timer to track how long your execution takes.
 * 4 - Profile your solution and identify weaknesses.
 * 5 - Improve it!
 * 6 - Remember that this is now *your* code and *you* should modify it however
 *     needed to solve the assignment.
 *
 * See the lab notes for a guide on profiling and an introduction to
 * multithreading (see also threads.c which is reference by the lab notes).
 */


// Global flag
// Set to true when operating in debug mode to enable verbose logging
static bool debug = false;

// Macro to calculate the 1D index by flattening the 3D index
#define to1D(i,j,k,n) (((k) * n * n) + ((j) * n) + (i))

typedef struct
{
    int thread_id;      // Unique id of the worker thread
    int start;          // Start index of the worker thread
    int end;            // End index of the worker thread
    double *curr;
    double *next;
    double *source;
    int n;
    int threads;
    float delta;

} WorkerArgs;

void* worker (void* pargs)
{
    WorkerArgs* args = (WorkerArgs*)pargs;

        // TODO: solve Poisson's equation for the given inputs
        double i_component;
        double j_component;
        double k_component;

    for (int k=args->start; k<args->end; k++)
    {
        for (int j=0; j<args->n; j++)
        {
            for (int i=0; i<args->n; i++)
            {
                if (k==0) {   
                    args->next[to1D(i,j,k,args->n)] = 0.0;
                } else {
                    if (i==0)
                    {   
                        i_component = args->curr[to1D(i+1,j,k,args->n)]*2.0;
                    }
                    else if (i==(args->n-1))
                    {
                        i_component = args->curr[to1D(i-1,j,k,args->n)]*2.0;
                    }
                    else
                    {
                        i_component = args->curr[to1D(i+1,j,k,args->n)] + args->curr[to1D(i-1,j,k,args->n)];
                    }

                    if (j==0)
                    {   
                        j_component = args->curr[to1D(i,j+1,k,args->n)]*2.0;
                    }
                    else if (j==(args->n-1))
                    {
                        j_component = args->curr[to1D(i,j-1,k,args->n)]*2.0;
                    }
                    else
                    {
                        j_component = args->curr[to1D(i,j+1,k,args->n)] + args->curr[to1D(i,j-1,k,args->n)];
                    }
                    
                    if (k==(args->n-1))
                    {
                        k_component = args->curr[to1D(i,j,k-1,args->n)]*2.0;
                    }
                    else
                    {
                        k_component = args->curr[to1D(i,j,k+1,args->n)] + args->curr[to1D(i,j,k-1,args->n)];
                    }

                    args->next[to1D(i,j,k,args->n)] = (i_component + j_component + k_component
                                            - ((args->delta*args->delta)*args->source[to1D(i,j,k,args->n)])) / 6.0;

                }
            }
        }
    }
    return NULL;
}

int main (int argc, char **argv)
{
    // Default settings for solver
    int iterations = 10;
    int n = 5;
    int threads = 2;
    float delta = 1;

    // parse the command line arguments
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp (argv[i], "-h") == 0 || strcmp (argv[i], "--help") == 0)
        {
            printf ("Usage: poisson [-n size] [-i iterations] [-t threads] [--debug]\n");
            return EXIT_SUCCESS;
        }

        if (strcmp (argv[i], "-n") == 0)
        {
            if (i == argc - 1)
            {
                fprintf (stderr, "Error: expected size after -n!\n");
                return EXIT_FAILURE;
            }

            n = atoi (argv[++i]);
        }

        if (strcmp (argv[i], "-i") == 0)
        {
            if (i == argc - 1)
            {
                fprintf (stderr, "Error: expected iterations after -i!\n");
                return EXIT_FAILURE;
            }

            iterations = atoi (argv[++i]);
        }

        if (strcmp (argv[i], "-t") == 0)
        {
            if (i == argc - 1)
            {
                fprintf (stderr, "Error: expected threads after -t!\n");
                return EXIT_FAILURE;
            }

            threads = atoi (argv[++i]);
        }

        if (strcmp (argv[i], "--debug") == 0)
        {
            debug = true;
        }
    }

    // Ensure we have an odd sized cube
    if (n % 2 == 0)
    {
        fprintf (stderr, "Error: n should be an odd number!\n");
        return EXIT_FAILURE;
    }

    // Create a source term with a single point in the centre
    double *source = (double*)calloc (n * n * n, sizeof (double));
    if (source == NULL)
    {
        fprintf (stderr, "Error: failed to allocated source term (n=%i)\n", n);
        return EXIT_FAILURE;
    }

    source[(n * n * n) / 2] = 1;

    // Storage for the thread handles and arguments
    // will exist for the entire lifetime of the program.
    pthread_t threads_list[threads];
    WorkerArgs args[threads];

    double *curr = (double*)calloc (n * n * n, sizeof (double)); 
    double *next = (double*)calloc (n * n * n, sizeof (double)); 
    double *temp; 

    if (debug)
    {
        printf ("Starting solver with:\n"
            "n = %i\n"
            "iterations = %i\n"
            "threads = %i\n"
            "delta = %f\n",
            n, iterations, threads, delta);
    }

    // Ensure we haven't run out of memory
    if (curr == NULL || next == NULL)
    {
        fprintf (stderr, "Error: ran out of memory when trying to allocate %i sized cube\n", n);
        exit (EXIT_FAILURE);
    }
    
    // Iterate over 3 dimensions i, j, and k
    // Need to consider the boundary conditions for each side of the cube
    for (int num=0; num<iterations; num++)
    {   
        // Launch each of the new worker threads
        for (int i = 0; i < threads; i++)
        {
            // Fill in the arguments to the worker
            args[i].thread_id = i;
            args[i].start = (n * i) / threads;
            args[i].end = (n * (i + 1)) / threads;
            args[i].curr = curr;
            args[i].next = next;
            args[i].n = n;
            args[i].threads = threads;
            args[i].delta = delta;
            args[i].source = source;
    
            // Create the worker thread
            if (pthread_create (&threads_list[i], NULL, &worker, &args[i]) != 0)
            {
                fprintf (stderr, "Error creating worker thread!\n");
                return EXIT_FAILURE;
            }
        }

        // Wait for all the threads to finish using join ()
        for (int i = 0; i < threads; i++)
        {
            pthread_join (threads_list[i], NULL);
        }

        // Set next as the new current and to rewrite next in the next iteration.
        temp = curr; // a new pointer to the memory space of curr
        curr = next; // curr now points to the memory space of next
        next = temp; // next now points to memory space of temp which is the old memory space of curr
    }

    // Free one of the buffers and return the correct answer in the other.
    // The caller is now responsible for free'ing the returned pointer.
    free (next);

    if (debug)
    {
        printf ("Finished solving.\n");
    }


    // Print out the middle slice of the cube for validation
    for (int x = 0; x < n; ++x)
    {
        for (int y = 0; y < n; ++y)
        {
            printf ("%0.5f ", curr[((n / 2) * n + y) * n + x]);
        }
        printf ("\n");
    }

    free (source);
    free (curr);

    return EXIT_SUCCESS;
}
