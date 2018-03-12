#pragma once
#ifndef MATH_API_H
#define MATH_API_H

#endif // MATH_API_H
#include <math.h>
#include "QTextStream"
class Vector2d
{
    friend Vector2d operator+(Vector2d a,Vector2d b);
    friend Vector2d operator*(double a,Vector2d b);
    friend Vector2d operator*(Vector2d a,double b);
    friend double operator*(Vector2d a,Vector2d b);
public:
    Vector2d();
    Vector2d(double);
    Vector2d(double a,double b);
    double x,y;
    double mag();
    Vector2d operator % (Vector2d a);
};
QTextStream &operator <<(QTextStream &s, const Vector2d &r);
std::ostream& operator <<(std::ostream& cout_,const Vector2d r);
