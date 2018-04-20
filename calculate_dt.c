/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <stdio.h>
#include<math.h>


double calculate_dt(Re,tau,dt,dx,dy,imax,jmax,U,V){
    temp1=((Re/2)*(pow(2,dx))*(pow(2,dy))/(pow(2,dx)+pow(2,dy));
    temp2=(dx)/abx(U);
    temp3=(dy)/abs(V);
    temp4=fmin(temp1,temp2);
    dt=tau*fmin(temp4,temp3);
    return dt;
}


void boundaryvalues(imax,jmax,U,V){
    
}
