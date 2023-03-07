#!/usr/bin/python3
################ZJU#####
#@auther: Hongfei Ye  
#date: 2016.4.6
#class: Computational aerodynamics
#project name:One-dimensional SOD shock tube 
#teacher:Shi Xing
#platform : python3.6.1 matplotlib 2.0.0 numpy 1.12.1
########################
from numpy import *
import matplotlib.pyplot as  pl
##########initial#########
delta_t=0.0001
delta_x=0.005#you must choose [0.0001 0.005] because of the sick problem
Nt=int(0.20/delta_t)
Nx=int(1/delta_x)
gama=1.4
P=array([0.4]*int(Nx/2)+[0.4]*int(Nx/2))
Den=array([1.0]*int(Nx/2)+[1.0]*int(Nx/2))
u=array([-2.0]*int(Nx/2)+[2.0]*int(Nx/2))

F_p=array([[0.0]*Nx,[0.0]*Nx,[0.0]*Nx])
F_n=array([[0.0]*Nx,[0.0]*Nx,[0.0]*Nx])
Den_u=Den*u
FF_p=array([[0.0]*Nx,[0.0]*Nx,[0.0]*Nx])
FF_n=array([[0.0]*Nx,[0.0]*Nx,[0.0]*Nx])
E=P/(gama-1)+0.5*u**2*Den
#########displacement########
def show(s):
    x=[[i/Nx-0.5] for i in range(Nx)]
    pl.figure(1)
    pl.plot(x,P,label="$Pressure$")
    pl.plot(x,u,label="$velocity$")
    pl.plot(x,Den,label="$Density$")
    pl.legend()
    pl.title(s)
    pl.xlabel("Positon x")
    pl.show()
#########compute#######
def VLsplit(level=1):
    global P,u,C,F_p,F_n,Den,E
    for j in range(Nt):
        C=sqrt(abs(gama*P/Den))
        ma=u/C
        lamta=array([u,u-C,u+C])
        for i in range(Nx):
            if(ma[i]>=1):
                F_p[0][i]=Den[i]*u[i]
                F_p[1][i]=u[i]*(Den[i]*u[i]**2+P[i])
                F_p[2][i]=u[i]*(E[i]+P[i])
                F_n[0][i]=0
                F_n[1][i]=0
                F_n[2][i]=0
            if(ma[i]<=-1):
                F_n[0][i]=Den[i]*u[i]
                F_n[1][i]=u[i]*(Den[i]*u[i]**2+P[i])
                F_n[2][i]=u[i]*(E[i]+P[i])
                F_p[0][i]=0
                F_p[1][i]=0
                F_p[2][i]=0
            else:
                fp=Den[i]*C[i]*((ma[i]+1)/2)**2
                fn=-1*Den[i]*C[i]*((ma[i]-1)/2)**2
                F_p[0][i]=fp
                F_p[1][i]=(fp/gama)*((gama-1)*u[i]+2*C[i])
                F_p[2][i]=(fp/(2*gama**2-2))*((gama-1)*u[i]+2*C[i])**2
                F_n[0][i]=fn
                F_n[1][i]=(fn/gama)*((gama-1)*u[i]-2*C[i])
                F_n[2][i]=(fn/(2*gama**2-2))*((gama-1)*u[i]-2*C[i])**2
        if(level==2):##choose the order of accuracy
            for i in range(Nx-4):
                Den[i+2]   -=(0.5*delta_t/delta_x)*(3*F_p[0][i+2]-4*F_p[0][i+1]+F_p[0][i]-F_n[0][i+4]+4*F_n[0][i+3]-3*F_n[0][i+2])
                Den_u[i+2] -=(0.5*delta_t/delta_x)*(3*F_p[1][i+2]-4*F_p[1][i+1]+F_p[1][i]-F_n[1][i+4]+4*F_n[1][i+3]-3*F_n[1][i+2])
                E[i+2]     -=(0.5*delta_t/delta_x)*(3*F_p[2][i+2]-4*F_p[2][i+1]+F_p[2][i]-F_n[2][i+4]+4*F_n[2][i+3]-3*F_n[2][i+2])
        else:
            for i in range(Nx-2):
                Den[i+1]-=   (delta_t/delta_x)*(F_p[0][i+1]+F_n[0][i+2]-F_p[0][i]-F_n[0][i+1])
                Den_u[i+1]-= (delta_t/delta_x)*(F_p[1][i+1]+F_n[1][i+2]-F_p[1][i]-F_n[1][i+1])
                E[i+1]-=     (delta_t/delta_x)*(F_p[2][i+1]+F_n[2][i+2]-F_p[2][i]-F_n[2][i+1])
        u=Den_u/Den
        P=(gama-1)*(E-0.5*Den*u**2)
    show("Van Leer split")
        
def ROE():
    global P,u,C,F_p,F_n,Den,E
    F1=[Den*u,Den*u**2+P,u*(E+P)]
    for j in range(Nt):
        G=[Den,Den_u,E]
        F=[Den_u,Den*u**2+P,u*(E+P)]
        H=(E+P)/Den
        for i in range(Nx-1):
            aveu=(sqrt(Den[i])*u[i]+sqrt(Den[i+1])*u[i+1])/(sqrt(Den[i])+sqrt(Den[i+1]))
            aveDen=((sqrt(Den[i])+sqrt(Den[i+1]))/2)**2
            aveH=(sqrt(Den[i])*H[i]+sqrt(Den[i+1])*H[i+1])/(sqrt(Den[i])+sqrt(Den[i+1]))
            aveP=((gama-1)/gama)*(aveDen*aveH-0.5*aveDen*aveu**2)
            aveE=aveDen*aveH-aveP
            aveDen_u=aveDen*aveu
            avec=((gama-1)*(aveH-0.5*aveu**2))
            aveA=array([[0,                                       1,                                  0],\
                        [((gama-3)/2)*aveu**2,                    (3-gama)*aveu,                      gama-1],\
                        [(gama-2)*aveu**3/2-aveu*avec/(gama-1),avec/(gama-1)+(3-gama)*aveu**2/2,gama*aveu]])
            aveA_1,aveS=linalg.eig(aveA)
            aveA_1=abs(diag(aveA_1))
            aveS_1=linalg.inv(aveS)
            D=dot(aveS,aveA_1)
            K=0.5*dot(aveS_1,array([x[i+1] for x in G])-array([x[i] for x in G]))
            K=dot(D,K)
            for k in range(3):
                F1[k][i+1]=0.5*(F[k][i]+F[k][i+1])-K[k]
            if(i!=0):
                Den[i]-=   (delta_t/delta_x)*(F1[0][i+1]-F1[0][i])
                Den_u[i]-= (delta_t/delta_x)*(F1[1][i+1]-F1[1][i])
                E[i]-=     (delta_t/delta_x)*(F1[2][i+1]-F1[2][i])
        u=Den_u/Den
        P=(gama-1)*(E-0.5*Den*u**2)
    show("ROE")
 
    
def WENO():
    global P,u,C,F_p,F_n,Den,E,FF_p,FF_n,IS1_p,IS1_n,IS2_p,IS2_n,IS3_p,IS3_n
    for j in range(Nt):
        F=[Den_u,Den*u**2+P,u*(E+P)]
        C=sqrt(abs(gama*P/Den))
        lamta=array([u,u-C,u+C])
        lamta_p=(abs(lamta)+lamta)/2
        lamta_n=lamta-lamta_p
        for i in range(Nx):
            F_p[0][i]=Den[i]/(2*gama)*(2*(gama-1)*lamta_p[0][i]+lamta_p[1][i]+lamta_p[2][i])
            F_p[1][i]=Den[i]/(2*gama)*(2*(gama-1)*lamta_p[0][i]*u[i]+lamta_p[1][i]*(u[i]-C[i])+lamta_p[2][i]*(u[i]+C[i]))
            F_p[2][i]=Den[i]/(2*gama)*((gama-1)*lamta_p[0][i]*u[i]**2+1/2*lamta_p[1][i]*(u[i]-C[i])**2+1/2*lamta_p[2][i]*(u[i]+C[i])**2+(3-gama)/(2*gama-2)*(lamta_p[1][i]+lamta_p[2][i])*(C[i])**2)
            F_n[0][i]=Den[i]/(2*gama)*(2*(gama-1)*lamta_n[0][i]+lamta_n[1][i]+lamta_n[2][i])
            F_n[1][i]=Den[i]/(2*gama)*(2*(gama-1)*lamta_n[0][i]*u[i]+lamta_n[1][i]*(u[i]-C[i])+lamta_n[2][i]*(u[i]+C[i]))
            F_n[2][i]=Den[i]/(2*gama)*((gama-1)*lamta_n[0][i]*u[i]**2+1/2*lamta_n[1][i]*(u[i]-C[i])**2+1/2*lamta_n[2][i]*(u[i]+C[i])**2+(3-gama)/(2*gama-2)*(lamta_n[1][i]+lamta_n[2][i])*(C[i])**2)
        for i in range(Nx-4):
            IS1_p=IS2_p=IS3_p=IS1_n=IS2_n=IS3_n=0
            for k in range(3):
                IS1_p+=1/4*((F_p[k][i]-4*F_p[k][i+1]+3*F_p[k][i+2])**2)+     13/12*((F_p[k][i]  -2*F_p[k][i+1]+F_p[k][i+2])**2)
                IS2_p+=1/4*((F_p[k][i+1]-F_p[k][i+3])**2)+                 13/12*((F_p[k][i+1]-2*F_p[k][i+2]+F_p[k][i+3])**2)
                IS3_p+=1/4*((3*F_p[k][i+2]-4*F_p[k][i+3]+F_p[k][i+4])**2)+   13/12*((F_p[k][i+2]-2*F_p[k][i+3]+F_p[k][i+4])**2)
                IS1_n+=1/4*((F_n[k][i+4]-4*F_n[k][i+3]+3*F_n[k][i+2])**2)+   13/12*((F_n[k][i+4]-2*F_n[k][i+3]+F_n[k][i+2])**2)
                IS2_n+=1/4*((F_n[k][i+3]-F_n[k][i+1])**2)+                 13/12*((F_n[k][i+3]-2*F_n[k][i+2]+F_n[k][i+1])**2)
                IS3_n+=1/4*((3*F_n[k][i+2]-4*F_n[k][i+1]+F_n[k][i])**2)+     13/12*((F_n[k][i+2]-2*F_n[k][i+1]+F_n[k][i])  **2)
            oumiga1_p=(0.1/(IS1_p+1e-6)**2/(0.1/(IS1_p+1e-6)**2+0.6/(IS2_p+1e-6)**2+0.3/(IS3_p+1e-6)**2))#attention,  oumiga is a num instead of an vector
            oumiga2_p=(0.6/(IS2_p+1e-6)**2/(0.1/(IS1_p+1e-6)**2+0.6/(IS2_p+1e-6)**2+0.3/(IS3_p+1e-6)**2))
            oumiga3_p=1-oumiga2_p-oumiga1_p
            oumiga1_n=(0.1/(IS1_n+1e-6)**2/(0.1/(IS1_n+1e-6)**2+0.6/(IS2_n+1e-6)**2+0.3/(IS3_n+1e-6)**2))
            oumiga2_n=(0.6/(IS2_n+1e-6)**2/(0.1/(IS1_n+1e-6)**2+0.6/(IS2_n+1e-6)**2+0.3/(IS3_n+1e-6)**2))
            oumiga3_n=1-oumiga2_n-oumiga1_n
            for k in range(3):
                FF_p[k][i+2]=oumiga1_p*(2*F_p[k][i]/6    -7*F_p[k][i+1]/6   +11*F_p[k][i+2]/6)+oumiga2_p*(-1*F_p[k][i+1]/6 +5*F_p[k][i+2]/6   +2*F_p[k][i+3]/6)+oumiga3_p*(2*F_p[k][i+2]/6  +5*F_p[k][i+3]/6   -1*F_p[k][i+4]/6)
                FF_n[k][i+2]=oumiga1_n*(2*F_n[k][i+4]/6  -7*F_n[k][i+3]/6   +11*F_n[k][i+2]/6)+oumiga2_n*(-1*F_n[k][i+3]/6 +5*F_n[k][i+2]/6   +2*F_n[k][i+1]/6)+oumiga3_n*(2*F_n[k][i+2]/6  +5*F_n[k][i+1]/6     -1*F_n[k][i]/6)
        for i in range(Nx-4):
            if(i>=2):
                Den[i+1]-=   (delta_t/delta_x)*(FF_p[0][i+1]-FF_p[0][i]+FF_n[0][i+2]-FF_n[0][i+1])
                Den_u[i+1]-= (delta_t/delta_x)*(FF_p[1][i+1]-FF_p[1][i]+FF_n[1][i+2]-FF_n[1][i+1])
                E[i+1]-=     (delta_t/delta_x)*(FF_p[2][i+1]-FF_p[2][i]+FF_n[2][i+2]-FF_n[2][i+1])
        u=Den_u/Den
        P=(gama-1)*(E-0.5*Den*u**2)
    show("WENO")
############choose function##########
VLsplit()
#VLsplit(2)
#ROE()
#WENO()
 
