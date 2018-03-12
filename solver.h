#ifndef SOLVER_H
#define SOLVER_H


class solver
{
private:
    double d;
public:
    solver();
    void iterate(double *x,double *y,int N);

};

#endif // SOLVER_H
