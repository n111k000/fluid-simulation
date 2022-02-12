#include <math.h>
#include <Windows.h>
#include <stdio.h>
#include <iostream>
#include <SDL.h>
#include <vector>
#include <fstream>


#define IX(i,j) ((i)+(N)*(j))



struct FluidSqare {
    int size;
    float dt;
    float diff;
    float visc;

    float* s;
    float* density;

    float* Vx;
    float* Vy;

    float* Vx0;
    float* Vy0;

};
typedef struct FluidSqare FluidSqare;

FluidSqare* FluidSquareCreate(int size, int diffusion, int viscosity, float dt)
{

    FluidSqare* square = (FluidSqare*)malloc(sizeof(*square));
    int N = size;

    square->size = size;
    square->dt = dt;
    square->diff = diffusion;
    square->visc = viscosity;

    square->s = (float*)calloc(N * N, sizeof(float));
    if (square->s == NULL) exit(1);
    square->density = (float*)calloc(N * N, sizeof(float));
    if (square->density == NULL) exit(1);

    square->Vx = (float*)calloc(N * N, sizeof(float));
    if (square->Vx == NULL) exit(1);
    square->Vy = (float*)calloc(N * N, sizeof(float));
    if (square->Vy == NULL) exit(1);

    square->Vx0 = (float*)calloc(N * N, sizeof(float));
    if (square->Vx0 == NULL) exit(1);
    square->Vy0 = (float*)calloc(N * N, sizeof(float));
    if (square->Vy0 == NULL) exit(1);

    return square;
}

void FluidSquareFree(FluidSqare* square)
{
    free(square->s);
    free(square->density);

    free(square->Vx);
    free(square->Vy);

    free(square->Vx0);
    free(square->Vy0);

    free(square);
}

void FluidSquareAddDensity(FluidSqare* square, int x, int y, float amount)
{
    int N = square->size;
    square->density[IX(x, y)] += amount;
}

void FluidSquareAddVelocity(FluidSqare* square, int x, int y, float amountX, float amountY)
{
    int N = square->size;
    int index = IX(x, y);

    square->Vx[index] += amountX;
    square->Vy[index] += amountY;
}

static void set_bnd(int b, float* x, int N)
{

    for (int i = 1; i < N - 1; i++) {
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
    }
    for (int j = 1; j < N - 1; j++) {
        x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
        x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
    }

    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);

    x[IX(0, N - 1)] = 0.5f * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);

    x[IX(N - 1, 0)] = 0.5f * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);

    x[IX(N - 1, N - 1)] = 0.5f * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);


}

static void lin_solve(int b, float* x, float* x0, float a, float c, int iter, int N)
{
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {

        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[IX(i, j)] =
                    (x0[IX(i, j)]
                        + a * (x[IX(i + 1, j)]
                            + x[IX(i - 1, j)]
                            + x[IX(i, j + 1)]
                            + x[IX(i, j - 1)]
                            )) * cRecip;

            }
        }
        set_bnd(b, x, N);
    }
}

static void diffuse(int b, float* x, float* x0, float diff, float dt, int iter, int N)
{
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
}

static void project(float* velocX, float* velocY, float* p, float* div, int iter, int N)
{

    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j)] = -0.5f * (
                velocX[IX(i + 1, j)]
                - velocX[IX(i - 1, j)]
                + velocY[IX(i, j + 1)]
                - velocY[IX(i, j - 1)]
                ) / N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div, N);
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1, 6, iter, N);


    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)]
                - p[IX(i - 1, j)]) * N;
            velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)]
                - p[IX(i, j - 1)]) * N;
        }
    }
    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
}

static void advect(int b, float* d, float* d0, float* velocX, float* velocY, float dt, int N)
{
    float i0, i1, j0, j1;

    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;


    for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
        for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            tmp1 = dtx * velocX[IX(i, j)];
            tmp2 = dty * velocY[IX(i, j)];
            x = ifloat - tmp1;
            y = jfloat - tmp2;

            if (x < 0.5f) x = 0.5f;
            if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
            i0 = floorf(x);
            i1 = i0 + 1.0f;
            if (y < 0.5f) y = 0.5f;
            if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
            j0 = floorf(y);
            j1 = j0 + 1.0f;

            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;

            int i0i = i0;
            int i1i = i1;
            int j0i = j0;
            int j1i = j1;

            d[IX(i, j)] = 
                
                s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)])
                + s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);

        }
    }
    set_bnd(b, d, N);
}

void FluidSquareStep(FluidSqare* square)
{
    int N = square->size;
    float visc = square->visc;
    float diff = square->diff;
    float dt = square->dt;
    float* Vx = square->Vx;
    float* Vy = square->Vy;
    float* Vx0 = square->Vx0;
    float* Vy0 = square->Vy0;
    float* s = square->s;
    float* density = square->density;

    diffuse(1, Vx0, Vx, visc, dt, 4, N);
    diffuse(2, Vy0, Vy, visc, dt, 4, N);

    project(Vx0, Vy0, Vx, Vy, 4, N);

    advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
    advect(2, Vy, Vy0, Vx0, Vy0, dt, N);

    project(Vx, Vy, Vx0, Vy0, 4, N);

    diffuse(0, s, density, diff, dt, 4, N);
    advect(0, density, s, Vx, Vy, dt, N);
}

//pretvaranje arraya u sliku




const int WIDTH = 1000, HEIGHT = 1000;

//sdl rect/(idk)
SDL_Rect newSDL_Rect(int xs, int ys, int widths, int heights) {
    SDL_Rect rectangular;
    rectangular.x = xs;
    rectangular.y = ys;
    rectangular.w = widths;
    rectangular.h = heights;
    return rectangular;
}

int main(int argc, char* argv[]/*, FluidSqare* square*/)
{
    //kreiranje ovoga za simulaciiju
    FluidSqare* square;
    square = FluidSquareCreate(100, 0, 0, .1);
    int N = square->size;

    //sdl defining??

    SDL_Window* window = NULL;
    SDL_Surface* surface = NULL;
    SDL_Renderer* renderer = NULL;
    SDL_Event ev;

    
    SDL_Rect rects[100][100];
    if (SDL_Init(SDL_INIT_VIDEO) < 0) //Init the video driver
    {
        printf("SDL_Error: %s\n", SDL_GetError());
    }
 
    window = SDL_CreateWindow("fluid simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, SDL_WINDOW_SHOWN); //Creates the window
    if (window == NULL)
    {
            printf("SDL_Error: %s\n", SDL_GetError());
    }

    renderer = SDL_CreateRenderer(window, 0, SDL_RENDERER_ACCELERATED); //renderer used to color rects
    if (renderer == NULL)
    {
        printf("SDL_Error %s\n", SDL_GetError());
    }
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);
       
    for (int i = 0; i < 100; i++)
    {
                
        for (int j = 0; j < 100; j++)        
        {
            rects[i][j] = newSDL_Rect(i * 10, j * 10, 10, 10);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderFillRect(renderer, &rects[i][j]);
                
        }

        SDL_UpdateWindowSurface(window);    
        SDL_RenderPresent(renderer);
    
    }

    
    while (true)
    {
        FluidSquareAddDensity(square, 5, 45, 2);

        FluidSquareAddVelocity(square, 5, 45, 5, -5);
   
       
        FluidSquareStep(square);
        
        SDL_PollEvent(&ev);

        if (ev.type == SDL_QUIT)
        {
            break;
        }

        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < 100; j++)
            {
                int dens = square->density[IX(i, j)] * 255;
                rects[i][j] = newSDL_Rect(i * 10, j * 10, 10, 10);
                SDL_SetRenderDrawColor(renderer, dens, dens, 0, 255);
                SDL_RenderFillRect(renderer, &rects[i][j]);
            }
        }

        SDL_UpdateWindowSurface(window);
        SDL_RenderPresent(renderer);
        SDL_Delay(50);

    }

    FluidSquareFree(square);
    SDL_DestroyWindow(window);
    SDL_DestroyRenderer(renderer);
    SDL_Quit();

    return EXIT_SUCCESS;
}
