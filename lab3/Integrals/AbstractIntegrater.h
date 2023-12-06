#ifndef ABSTRACT_INTEGRATER_H
#define ABSTRACT_INTEGRATER_H


class AbstractIntegrater
{

public:

    virtual double integrate(double (*function)(double ), double a, double b, double h) = 0;

    static double rungeRoombergMethod(double result_h1, double result_h2, int p);

public:

    virtual ~AbstractIntegrater() = default;

};


#endif //ABSTRACT_INTEGRATER_H
