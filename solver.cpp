#include "solver.h"
#include "QDebug"
#include "math_api.h"
#include "omp.h"
#include <time.h>
solver::solver()
{

}
void solver::iterate(double *x1,double *y1,int N){

    int const n=100,m=n; //number of lattice nodes in two dimentions
    double f[9][n][m],rho[n][m],feq,x[n],y[m],w[9],cu,rhon;
    Vector2d c[9],u[n][m],sumu;
    double dt=1.0,dy=1.0,dx=dy;
    double u0=0.1,rho0=5.;
    int i,j,k;
    int mstep=40000;//total number of steps
    x[0]=0;
    for(i=1;i<n;i++){
        x[i]=x[i-1]+dx;
    }
    y[0]=0;
    for(i=1;i<m;i++){
        y[i]=y[i-1]+dy;
    }
    double alpha=0.01;
    double Re=u0*m/alpha;
    double sum;
    //double csq=dx*dx/(dt*dt);

    double omega=1.0/(3.*alpha+0.5);
// start wall clock
   double wall_timer = omp_get_wtime();
   clock_t clock_timer = clock();

    w[0]=4./9;
    w[1]=w[2]=w[3]=w[4]=1./9;
    w[5]=w[6]=w[7]=w[8]=1./36;

    c[0]=Vector2d();
    c[1]=Vector2d(1,0);c[2]=Vector2d(0,1);c[3]=Vector2d(-1,0);c[4]=Vector2d(0,-1);
    c[5]=Vector2d(1,1);c[6]=Vector2d(-1,1);c[7]=Vector2d(-1,-1);c[8]=Vector2d(1,-1);




    //initiate solution

    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            rho[i][j]=rho0;
            u[i][j]=Vector2d();
            for(k=0;k<9;k++){
                f[k][i][j]=w[k]*rho[i][j];// zero initial velocity
            }
        }
    }

    //main loop
    for(int t=0;t<mstep;t++){
       // qDebug()<<t;

        for(i=0;i<n;i++){
            for(j=0;j<m;j++){
                //compute rho and u
                sum=0.;
                sumu=Vector2d();

                for(k=0;k<9;k++){
                    sum=sum+f[k][i][j];
                    sumu=sumu+f[k][i][j]*c[k];
                }
                rho[i][j]=sum;
                u[i][j]=sumu*(1/rho[i][j]);
                //collision
                for(k=0;k<9;k++){
                    cu=c[k]*u[i][j];
                    feq=w[k]*rho[i][j]*(1+3*cu+4.5*cu*cu-1.5*(u[i][j]*u[i][j]));
                    f[k][i][j]=omega*feq+(1.-omega)*f[k][i][j];
               }
           }

       }
        //streaming

        for(i=0;i<n;i++){
            for(j=0;j<(m-1);j++){
                f[3][j][i]=f[3][j+1][i];
                f[1][n-j-1][i]=f[1][n-j-2][i];
                f[4][i][j]=f[4][i][j+1];
                f[2][i][m-j-1]=f[2][i][m-j-2];
            }
        }
        for(i=0;i<(n-1);i++){
            for(j=0;j<(m-1);j++){
                f[5][n-i-1][m-j-1]=f[5][n-i-2][m-j-2];
                f[6][i][m-j-1]=f[6][i+1][m-j-2];
                f[7][i][j]=f[7][i+1][j+1];
                f[8][n-i-1][j]=f[8][n-i-2][j+1];
            }
        }
        //boundary conditions

        for(i=0;i<m;i++){
            //x=0 bounce back
            f[1][0][i]=f[3][0][i];
            f[5][0][i]=f[7][0][i];
            f[8][0][i]=f[6][0][i];
            //x=l bounce back
            f[3][n-1][i]=f[1][n-1][i];
            f[7][n-1][i]=f[5][n-1][i];
            f[6][n-1][i]=f[8][n-1][i];
            //y=0 bounce back
            f[5][i][0]=f[7][i][0];
            f[2][i][0]=f[4][i][0];
            f[6][i][0]=f[8][i][0];
            //y=l moving lid
            rhon=f[0][i][m-1]+f[1][i][m-1]+f[3][i][m-1]+2.*(f[2][i][m-1]+f[6][i][m-1]+f[5][i][m-1]);
            f[7][i][m-1]=f[5][i][m-1]-rhon*u0/6.;
            f[4][i][m-1]=f[2][i][m-1];
            f[8][i][m-1]=f[6][i][m-1]+rhon*u0/6.;

        }

    }

   qDebug()<<   " time on wall: " <<  omp_get_wtime() - wall_timer << "\n";

    for (i=0; i<N; i++)//Пробегаем по всем точкам
    {


            x1[i]=x[i];
            y1[i] =u[i][50].y;// rho[i][0];//Формула нашей функции

    }

}
