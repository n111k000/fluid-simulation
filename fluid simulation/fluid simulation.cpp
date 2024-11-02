#include <Windows.h>
#include <iostream>
#include <SDL.h>
#include <thread>


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

FluidSqare* FluidSquareCreate(int size, float diffusion, float viscosity, float dt)
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


//definiranje vlastitih
struct Bounds{
    int* kGRx;
    int* kGRy;
    int brGR;
    bool* GR;
};
typedef struct Bounds Bounds;

Bounds* BoundsCreate(int brGR, int N) {
    Bounds* bounds = (Bounds*)malloc(sizeof(*bounds));

    bounds->brGR = brGR;
    bounds->kGRx = new int[brGR] {};
    bounds->kGRy = new int[brGR] {};
    bounds->GR = (bool*)calloc(N * N, sizeof(bool));

    return bounds;
}

void BoundsDefine(Bounds* bounds, int N) {
    int brGR = bounds->brGR;
    int* kGRx = bounds->kGRx;
    int* kGRy = bounds->kGRy;

    for (int i = 0; i < brGR; i++) {
        bounds->GR[IX(kGRx[i], kGRy[i])] = 1;
    }
}


void BoundsFree(Bounds* bounds) {
    free(bounds->kGRx);
    free(bounds->kGRy);
    free(bounds->GR);

    free(bounds);
}

//postavljanje granica, x-1,y-2,ost-0
static void set_bnd(int b, float* x, int N, Bounds* bounds)
{
    int* kGRx = bounds->kGRx;
    int* kGRy = bounds->kGRy;
    int brGR = bounds->brGR;
    bool* GR = bounds->GR;

    //za stranice koje omeduju
    for (int i = 1; i < N - 1; i++) {
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
    }
    for (int j = 1; j < N - 1; j++) {
        x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
        x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
    }

    //za kuteve koji omeduju
    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);

    x[IX(0, N - 1)] = 0.5f * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);

    x[IX(N - 1, 0)] = 0.5f * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);

    x[IX(N - 1, N - 1)] = 0.5f * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);

    //za vlastite granice

    for (int k = 0; k < brGR; k++) {
        //ravno
        bool u = kGRy[k] != 0;
        bool d = kGRy[k] != N - 1;
        bool r = kGRx[k] != N - 1;
        bool l = kGRx[k] != 0;
        
        
        if (l) {
            x[IX(kGRx[k] - 1, kGRy[k])] = b == 1 ? -abs(x[IX(kGRx[k] - 1, kGRy[k])]) : x[IX(kGRx[k] - 1, kGRy[k])];
        }
        if (r) {
            x[IX(kGRx[k] + 1, kGRy[k])] = b == 1 ? abs(x[IX(kGRx[k] + 1, kGRy[k])]) : x[IX(kGRx[k] + 1, kGRy[k])];
        }
        if (u) {
            x[IX(kGRx[k], kGRy[k] - 1)] = b == 2 ? -abs(x[IX(kGRx[k], kGRy[k] - 1)]) : x[IX(kGRx[k], kGRy[k] - 1)];
        }
        if (d) {
            x[IX(kGRx[k], kGRy[k] + 1)] = b == 2 ? abs(x[IX(kGRx[k], kGRy[k] + 1)]) : x[IX(kGRx[k], kGRy[k] + 1)];
        }
        //kutevi
        

        

        if(l && r && u && d){

            //lijevo gore desno dolje
            bool pl = GR[IX(kGRx[k] - 1, kGRy[k])];
            bool pu = GR[IX(kGRx[k], kGRy[k] - 1)];
            bool pr = GR[IX(kGRx[k] + 1, kGRy[k])];
            bool pd = GR[IX(kGRx[k], kGRy[k] + 1)];

            if (pr == 0 && pd == 0) {
                x[IX(kGRx[k] + 1, kGRy[k] + 1)] = b == 2 ? abs(x[IX(kGRx[k] - 1, kGRy[k] - 1)]) : abs(x[IX(kGRx[k] - 1, kGRy[k])]);
            }
            if (pl == 0 && pd == 0) {
                x[IX(kGRx[k] - 1, kGRy[k] + 1)] = b == 2 ? abs(x[IX(kGRx[k] - 1, kGRy[k] - 1)]) : -abs(x[IX(kGRx[k] - 1, kGRy[k])]);
            }
            if (pu == 0 && pr == 0) {
                x[IX(kGRx[k] + 1, kGRy[k] - 1)] = b == 2 ? -abs(x[IX(kGRx[k] - 1, kGRy[k] - 1)]) : abs(x[IX(kGRx[k] - 1, kGRy[k])]);
            }
            if (pl == 0 && pu == 0) {
                x[IX(kGRx[k] - 1, kGRy[k] - 1)] = b == 2 ? -abs(x[IX(kGRx[k] - 1, kGRy[k] - 1)]) : -abs(x[IX(kGRx[k] - 1, kGRy[k])]);
            }
            
        }
    }
}

//linearno rješavanje
static void lin_solve(int b, float* x, float* x0, float a, float c, int iter, int N, Bounds* bounds)
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
        set_bnd(b, x, N, bounds);
    }
}

//difuzija
static void diffuse(int b, float* x, float* x0, float diff, float dt, int iter, int N, Bounds* bounds)
{
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N, bounds);
}

//projekcija
static void project(float* velocX, float* velocY, float* p, float* div, int iter, int N, Bounds* bounds)
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
    set_bnd(0, div, N, bounds);
    set_bnd(0, p, N, bounds);
    lin_solve(0, p, div, 1, 6, iter, N, bounds);


    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)]
                - p[IX(i - 1, j)]) * N;
            velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)]
                - p[IX(i, j - 1)]) * N;
        }
    }
    set_bnd(1, velocX, N, bounds);
    set_bnd(2, velocY, N, bounds);
}

//advekcija
static void advect(int b, float* d, float* d0, float* velocX, float* velocY, float dt, int N, Bounds* bounds)
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
    set_bnd(b, d, N, bounds);
}

//napravi korak simulcije
void FluidSquareStep(FluidSqare* square, Bounds* bounds)
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

    diffuse(1, Vx0, Vx, visc, dt, 4, N, bounds);
    diffuse(2, Vy0, Vy, visc, dt, 4, N, bounds);

    project(Vx0, Vy0, Vx, Vy, 4, N, bounds);

    advect(1, Vx, Vx0, Vx0, Vy0, dt, N, bounds);
    advect(2, Vy, Vy0, Vx0, Vy0, dt, N, bounds);

    project(Vx, Vy, Vx0, Vy0, 4, N, bounds);

    diffuse(0, s, density, diff, dt, 4, N, bounds);
    advect(0, density, s, Vx, Vy, dt, N, bounds);
}

//pauza

void Pause() {
    SDL_Event ev;
    
    while (true) {
        SDL_PollEvent(&ev);
        if (ev.type == SDL_KEYDOWN && ev.key.keysym.scancode == SDL_SCANCODE_P) {
            break;
        }
    }
}

//prolazak kroz sve tocke nastanka plina i brzina
void GaVAdding(FluidSqare* square, int* SPx, int* SPy, int* SPam, int* VELx, int* VELy, float* VELxam, float* VELyam, int noSP, int noVEL) {

    for (int i = 0; i < noSP; i++) {
        FluidSquareAddDensity(square, SPx[i], SPy[i], SPam[i]);
    }

    for (int i = 0; i < noVEL; i++) {
        FluidSquareAddVelocity(square, VELx[i], VELy[i], VELxam[i], VELyam[i]);
    }

}


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
    square = FluidSquareCreate(100, 0.00000000000001, 0.0000051, 0.1);
    int N = square->size;

    //unosenje podataka
    std::cout << "pri unošenju brojevi koordinata moraju biti od 0 do 100";
    int noSP;
    std::cout << "Broj mjesta nastanka plina: ";
    std::cin >> noSP;

    int* SPx = new int[noSP] {};
    int* SPy = new int[noSP] {};
    int* SPam = new int[noSP] {};

    for (int i = 0; i < noSP; i++) {
        std::cout << "Unesi vrijednost x " << i+1 << ": ";
        std::cin >> SPx[i];
        std::cout << "Unesi vrijednost y " << i + 1 << ": ";
        std::cin >> SPy[i];
        std::cout << "Unesi kolicinu " << i + 1 << ": ";
        std::cin >> SPam[i];
    }

    int noVEL;
    std::cout << "Broj mjesta ubrzanja: ";
    std::cin >> noVEL;

    int* VELx = new int[noVEL] {};
    int* VELy = new int[noVEL] {};
    float* VELxam = new float[noVEL] {};
    float* VELyam = new float[noVEL] {};

    for (int i = 0; i < noVEL; i++) {
        std::cout << "Unesi vrijednost x " << i + 1 << ": ";
        std::cin >> VELx[i];
        std::cout << "Unesi vrijednost y " << i + 1 << ": ";
        std::cin >> VELy[i];
        std::cout << "Unesi kolicinu po x za " << i + 1 << ": ";
        std::cin >> VELxam[i];
        std::cout << "Unesi kolicinu po y za " << i + 1 << ": ";
        std::cin >> VELyam[i];
    }


    //vlastite granice
    Bounds* bounds;
    bounds = BoundsCreate(15, N);


    bounds->kGRx[0] = 48; bounds->kGRy[0] = 48;
    bounds->kGRx[1] = 48; bounds->kGRy[1] = 49;
    bounds->kGRx[2] = 48; bounds->kGRy[2] = 50;
    bounds->kGRx[3] = 48; bounds->kGRy[3] = 51;
    bounds->kGRx[4] = 48; bounds->kGRy[4] = 52;
    bounds->kGRx[10] = 48; bounds->kGRy[10] = 47;
    bounds->kGRx[11] = 48; bounds->kGRy[11] = 46;
    bounds->kGRx[12] = 48; bounds->kGRy[12] = 45;
    bounds->kGRx[13] = 48; bounds->kGRy[13] = 44;
    bounds->kGRx[14] = 48; bounds->kGRy[14] = 43;
    bounds->kGRx[5] = 48; bounds->kGRy[5] = 42;
    bounds->kGRx[6] = 48; bounds->kGRy[6] = 41;
    bounds->kGRx[7] = 48; bounds->kGRy[7] = 40;
    bounds->kGRx[8] = 48; bounds->kGRy[8] = 39;
    bounds->kGRx[9] = 48; bounds->kGRy[9] = 38;

    
    
    

    BoundsDefine(bounds, N);

    //sdl defining??

    SDL_Window* window = NULL;
    SDL_Surface* surface = NULL;
    SDL_Renderer* renderer = NULL;
    SDL_Event ev;

    Uint64 mouseButtons;
    int mouseX, mouseY;
    int mouseXprev = 0, mouseYprev = 0;
    
    
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


    bool QT = false;
    
    while (true)
    {
        
   
        SDL_PollEvent(&ev);
       
        switch (ev.type){
        
        case SDL_QUIT:

            QT = true;

        case SDL_KEYDOWN:

            switch (ev.key.keysym.scancode){

            case SDL_SCANCODE_P:

                Pause();

            }

        }

        if (QT){
            break;
        }

        //napravi da se dodaje upisan broj tocaka za dodavanje plinova

        GaVAdding(square, SPx, SPy, SPam, VELx, VELy, VELxam, VELyam, noSP, noVEL);

               
        
        FluidSquareStep(square, bounds);

        //pretvara iz podataka u sliku
        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < 100; j++)
            {
                int dens = square->density[IX(i, j)] * 255/1.5;
                dens > 255 ? dens = 255 : dens;
                rects[i][j] = newSDL_Rect(i * 10, j * 10, 10, 10);
                SDL_SetRenderDrawColor(renderer, dens, dens, 0, 255);
                SDL_RenderFillRect(renderer, &rects[i][j]);
            }
        }
        for (int k = 0; k < bounds->brGR; k++)
        {
            SDL_SetRenderDrawColor(renderer, 0, 0, 255, 100);
            SDL_RenderFillRect(renderer, &rects[bounds->kGRx[k]][bounds->kGRy[k]]);
        }
        

        SDL_UpdateWindowSurface(window);
        SDL_RenderPresent(renderer);
        SDL_Delay(20);

    }
    delete[] SPx;
    delete[] SPy;
    delete[] SPam;
    delete[] VELx;
    delete[] VELy;
    delete[] VELxam;
    delete[] VELyam;
    FluidSquareFree(square);
    BoundsFree(bounds);
    SDL_DestroyWindow(window);
    SDL_DestroyRenderer(renderer);
    SDL_Quit();
    
    return EXIT_SUCCESS;
}
