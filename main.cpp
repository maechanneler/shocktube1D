//
//  main.cpp
//  shocktube1D
//
//  Created by Yoshiaki Maehara on 2017/04/23.
//  Copyright © 2017年 Yoshiaki Maehara. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cmath>

//---- definition of fundamental constants
static const int    mdx   = 2001;
static const int    nlast =  625;
static const double gam   =     1.4;

static const double rho1  =     1.0;
static const double p1    =     1.0;
static const double u1    =     0.0;

static const double rho2  =     0.125;
static const double p2    =     0.1;
static const double u2    =     0.0;

//---- definition of grid parameters
static const double xmin =    -0.5;
static const double xmax =     0.5;
static const int    mx   =  1001;
static const double dx   = (xmax-xmin)/(mx-1);
static const double cfl  =     0.4;
static const double dt   = cfl*dx;

//---- field variables
static double x[mdx];
static double rho[mdx], p[mdx], u[mdx];
static double m[mdx], e[mdx], H[mdx];

static double rhop[mdx], up[mdx], Hp[mdx], cp[mdx];
static double Lambda[3], R[3][3], L[3][3], w[3][3];
static double E1[mdx], E2[mdx], E3[mdx];

//---- prototype of functions
void setgrd(void);
void slvflw(void);
void intcnd(void);
void bndcnd(void);

//---- set xy-grid system
void setgrd(void){
    for(int i=0; i<=mx; i++){
        x[i] = xmin+dx*(i-1);
    }
}

//---- initial condition
void intcnd(void){
    //---- uniform flow condition
    for(int i=0; i<=mx; i++){
        if(x[i]<=0.0){
            u[i]   = u1;
            p[i]   = p1;
            rho[i] = rho1;
        }else{
            u[i]   = u2;
            p[i]   = p2;
            rho[i] = rho2;
        }
    }
    
    //---- conservative quantity
    for(int i=0; i<=mx; i++){
        m[i] = rho[i]*u[i];
        e[i] = p[i]/(gam-1.0)+rho[i]*pow(u[i], 2.0)/2.0;
        H[i] = (e[i]+p[i])/rho[i];
    }
}

//---- boundary condition
void bndcnd(void){
    //---- from conservative to state
    for(int i=1; i<=mx-1; i++){
        u[i]   = m[i]/rho[i];
        p[i]   = (gam-1.0)*(e[i]-rho[i]*pow(u[i], 2.0)/2.0);
    }
    
    //---- inflow condition
    u[0]   = u1;
    p[0]   = p1;
    rho[0] = rho1;
    
    //---- outflow condition
    u[mx]   =  2.0*u[mx-1]-u[mx-2];
    p[mx]   =  2.0*p[mx-1]-p[mx-2];
    rho[mx] =  2.0*rho[mx-1]-rho[mx-2];
}

//---- solve 1D-Euler equation
void slvflw(void){
    //---- ROE's average
    for(int i=1; i<=mx-1; i++){
        rhop[i] = sqrt(rho[i]*rho[i+1]);
        up[i]   = (sqrt(rho[i])*u[i]+sqrt(rho[i+1])*u[i+1])/(sqrt(rho[i])+sqrt(rho[i+1]));
        Hp[i]   = (sqrt(rho[i])*H[i]+sqrt(rho[i+1])*H[i])/(sqrt(rho[i])+sqrt(rho[i+1]));
        cp[i]   = sqrt((gam-1.0)*(Hp[i]-up[i]*pow(up[i], 2.0))/2.0);
    }

    //----
    for(int i=1; i<=mx-1; i++){
        Lambda[0] = fabs(up[i]-cp[i]);
        Lambda[1] = fabs(up[i]);
        Lambda[2] = fabs(up[i]+cp[i]);
        
        R[0][0] = 1.0;
        R[0][1] = 1.0;
        R[0][2] = 1.0;
        R[1][0] = up[i]-cp[i];
        R[1][1] = up[i];
        R[1][2] = up[i]+cp[i];
        R[2][0] = Hp[i]-up[i]*cp[i];
        R[2][1] = pow(up[i], 2.0);
        R[2][2] = Hp[i]+up[i]*cp[i];
        
        double b1 = (gam-1.0)/2.0*pow(up[i]/cp[i], 2.0);
        double b2 = (gam-1.0)/pow(cp[i], 2.0);
        
        L[0][0] =  0.5*(b1+up[i]/cp[i]);
        L[0][1] = -0.5*(1.0/cp[i]+b2*up[i]);
        L[0][2] =  0.5*b2;
        L[1][0] =  1.0-b1;
        L[1][1] =  b2*up[i];
        L[1][2] = -b2;
        L[2][0] =  0.5*(b1-up[i]/cp[i]);
        L[2][1] =  0.5*(1.0/cp[i]-b2*up[i]);
        L[2][2] =  0.5*b2;
        
        for(int n1=0; n1<=2; n1++){
            for(int n2=0; n2<=2; n2++){
                w[n1][n2] = 0.0;
                for(int k=0; k<=2; k++){
                    w[n1][n2] = w[n1][n2]+R[n1][k]*Lambda[k]*L[k][n2];
                }
            }
        }
        
        E1[i] = 0.5*(m[i+1]+m[i]
                     -w[0][0]*(rho[i+1]-rho[i])
                     -w[0][1]*(m[i+1]-m[i])
                     -w[0][2]*(e[i+1]-e[i]));
        E2[i] = 0.5*((gam-1.0)*(e[i+1]+e[i])+(3.0-gam)/2.0*(pow(m[i+1], 2.0)/rho[i+1]+pow(m[i], 2.0)/rho[i])
                     -w[1][0]*(rho[i+1]-rho[i])
                     -w[1][1]*(m[i+1]-m[i])
                     -w[1][2]*(e[i+1]-e[i]));
        E3[i] = 0.5*(gam*(e[i+1]*m[i+1]/rho[i+1]+e[i]*m[i]/rho[i])
                     -(gam-1.0)/2.0*(pow(m[i+1]/rho[i+1], 2.0)*m[i+1]+pow(m[i]/rho[i], 2.0)*m[i])
                     -w[2][0]*(rho[i+1]-rho[i])
                     -w[2][1]*(m[i+1]-m[i])
                     -w[2][2]*(e[i+1]-e[i]));
    }

    for(int i=2; i<=mx-2; i++){
        rho[i] = rho[i]-dt/dx*(E1[i]-E1[i-1]);
        m[i]   = m[i]-dt/dx*(E2[i]-E2[i-1]);
        e[i]   = e[i]-dt/dx*(E3[i]-E3[i-1]);
    }
    
    bndcnd();
}

int main(void) {
    setgrd();
    intcnd();

    for(int n=0; n<=nlast; n++){
        slvflw();
        std::cout << n << std::endl;
    }

    for(int i=0; i<=mx; i++){
        std::cout << x[i] << " " <<p[i] << std::endl;
    }

    //---- write results
    std::ofstream ofp("/Users/Yoshiaki/desktop/output/shocktube1D/pressure.txt");
    std::ofstream ofu("/Users/Yoshiaki/desktop/output/shocktube1D/velocity.txt");
    std::ofstream ofr("/Users/Yoshiaki/desktop/output/shocktube1D/density.txt");

    for(int i=0; i<=mx; i++){
        ofp  << x[i] << " " << std::flush;
        ofu  << x[i] << " " << std::flush;
        ofr  << x[i] << " " << std::flush;
    }
    
    ofp  << std::endl;
    ofu  << std::endl;
    ofr  << std::endl;
    
    for(int i=0; i<=mx; i++){
        ofp  << p[i]    << " " << std::flush;
        ofu  << u[i]    << " " << std::flush;
        ofr  << rho[i]  << " " << std::flush;
    }
    
    ofp  << std::endl;
    ofu  << std::endl;
    ofr  << std::endl;
    
    return 0;
}
