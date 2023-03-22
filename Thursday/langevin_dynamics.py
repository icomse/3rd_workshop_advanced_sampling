# codes for running langevin dynamics and op calculations for 2D potentials
# langevin_dynamics.py
# Modified by PH Minh and Naomi Trampe
# last updated: 03/19/23

import math
import sys
from math import exp
import numpy as np


def potential(pes_type,x,y,well=4):
    # For PES-1
    if pes_type == 1:
        return 0.02*(x**4+y**4) - well*exp(-((x+2)**2 + (y+2)**2)) - well*exp(-((x-2)**2 + (y-2)**2)) + 0.3*(x-y)**2 + 2*well*exp(-8)
    # For PES-2
    elif pes_type == 2:
        return 0.03*(x**4+y**4) - well*exp(-((x+2)**2 + (y+2)**2)) - well*exp(-((x-2)**2 + (y-2)**2)) + 0.4*(x-y)**2 + 4*exp(-(x**2+y**2)) - 2.1245
    # For PES-3
    elif pes_type == 3:
        return 0.02*(x**4+y**4) - 3.73*exp(-((x+2)**2/8 + (y+2)**2/8)) - 3.73*exp(-((x-2)**2/8 + (y-2)**2/8)) + 3*exp(-((x)**2/2 + (y)**2/15)) + 2*exp(-((x)**2/2 + (y)**2/2)) - 0.5085  
    # For PES-4
    elif pes_type == 4:
        # Known constants
        A = [-8, -4, -6.8, 0.6]
        a = [-0.111, -0.111, -0.722, 0.0778]
        b = [0, 0, 1.22, 0.0667]
        c = [-1.11, -1.11, -0.722, 0.0778]
        x0 = [3, 0, -1.5, -3]
        y0 = [-3, -1.5, 1.5, 0]
        # Muller-Brown Potential Equation - Kob 2017 arXiv (https://arxiv.org/pdf/1701.01241.pdf)
        return A[0]*exp(a[0]*((x-x0[0])**2) + b[0]*(x-x0[0])*(y-y0[0]) + c[0]*((y-y0[0])**2)) + A[1]*exp(a[1]*((x-x0[1])**2) + b[1]*(x-x0[1])*(y-y0[1]) + c[1]*((y-y0[1])**2)) + A[2]*exp(a[2]*((x-x0[2])**2) + b[2]*(x-x0[2])*(y-y0[2]) + c[2]*((y-y0[2])**2)) + A[3]*exp(a[3]*((x-x0[3])**2) + b[3]*(x-x0[3])*(y-y0[3]) + c[3]*((y-y0[3])**2)) + 1.6966
    # For invalid PES
    else:
        sys.exit("{} is an invalid choice of PES".format(pes_type))

def kinetic_energy(phasepoint):
    px = phasepoint[2]
    py = phasepoint[3]
    return 0.5*(px*px+py*py) 
 
def force(x,y,px,py,dt,beta,gamma,pes_type,well=4):
    std_dev = math.sqrt(2.0*gamma/(beta*dt))
    # For PES-1
    if pes_type == 1:
        dV_dx = 0.08*x**3 + 2*well*(x-2)*exp(-(x-2)**2-(y-2)**2) + 2*well*(x+2)*exp(-(x+2)**2-(y+2)**2) + 0.6*(x-y)
        dV_dy = -0.6*(x-y) + 2*well*(y-2)*exp(-(x-2)**2-(y-2)**2) + 2*well*(y+2)*exp(-(x+2)**2-(y+2)**2) + 0.08*y**3
    # For PES-2
    elif pes_type == 2:
        dV_dx = 0.12*x**3 + 2*well*(x-2)*exp(-(x-2)**2-(y-2)**2) + 2*well*(x+2)*exp(-(x+2)**2-(y+2)**2) + 0.8*(x-y) - 8*x*exp(-(x**2+y**2))
        dV_dy = 0.12*y**3 + 2*well*(y-2)*exp(-(x-2)**2-(y-2)**2) + 2*well*(y+2)*exp(-(x+2)**2-(y+2)**2) - 0.8*(x-y) - 8*y*exp(-(x**2+y**2))
    # For PES-3
    elif pes_type == 3:
        dV_dx = 0.08*x**3 - 2*x*exp(-x**2/2-y**2/2) - 3*x*exp(-x**2/2-y**2/15) - 0.9325*(-x-2)*exp(-(1/8)*(x+2)**2-(1/8)*(y+2)**2) - 0.9325*(2-x)*exp(-(1/8)*(x-2)**2 - (1/8)*(y-2)**2)
        dV_dy = -2*y*exp(-x**2/2 - y**2/2) - (2/5)*y*exp(-x**2/2 - y**2/15) - 0.9325*(-y-2)*exp(-(1/8)*(x+2)**2 - (1/8)*(y+2)**2) - 0.9325*(2-y)*exp(-(1/8)*(x-2)**2 - (1/8)*(y-2)**2) + 0.08*y**3
    # For PES-4
    elif pes_type == 4:
        A = [-8, -4, -6.8, 0.6]
        a = [-0.111, -0.111, -0.722, 0.0778]
        b = [0, 0, 1.22, 0.0667]
        c = [-1.11, -1.11, -0.722, 0.0778]
        x0 = [3, 0, -1.5, -3]
        y0 = [-3, -1.5, 1.5, 0]
        dV_dx = A[0]*(a[0]*2*(x-x0[0]) + b[0]*(y-y0[0]))*exp(a[0]*((x-x0[0])**2) + b[0]*(x-x0[0])*(y-y0[0]) + c[0]*((y-y0[0])**2)) + A[1]*(a[1]*2*(x-x0[1]) + b[1]*(y-y0[1]))*exp(a[1]*((x-x0[1])**2) + b[1]*(x-x0[1])*(y-y0[1]) + c[1]*((y-y0[1])**2)) + A[2]*(a[2]*2*(x-x0[2]) + b[2]*(y-y0[2]))*exp(a[2]*((x-x0[2])**2) + b[2]*(x-x0[2])*(y-y0[2]) + c[2]*((y-y0[2])**2)) + A[3]*(a[3]*2*(x-x0[3]) + b[3]*(y-y0[3]))*exp(a[3]*((x-x0[3])**2) + b[3]*(x-x0[3])*(y-y0[3]) + c[3]*((y-y0[3])**2)) + A[3]*(a[3]*2*(x-x0[3]) + b[3]*(y-y0[3]))*exp(a[3]*((x-x0[3])**2) + b[3]*(x-x0[3])*(y-y0[3]) + c[3]*((y-y0[3])**2))
        dV_dy = A[0]*(c[0]*2*(y-y0[0]) + b[0]*(x-x0[0]))*exp(a[0]*((x-x0[0])**2) + b[0]*(x-x0[0])*(y-y0[0]) + c[0]*((y-y0[0])**2)) + A[1]*(c[1]*2*(y-y0[1]) + b[1]*(x-x0[1]))*exp(a[1]*((x-x0[1])**2) + b[1]*(x-x0[1])*(y-y0[1]) + c[1]*((y-y0[1])**2)) + A[2]*(c[2]*2*(y-y0[2]) + b[2]*(x-x0[2]))*exp(a[2]*((x-x0[2])**2) + b[2]*(x-x0[2])*(y-y0[2]) + c[2]*((y-y0[2])**2)) + A[3]*(c[3]*2*(y-y0[3]) + b[3]*(x-x0[3]))*exp(a[3]*((x-x0[3])**2) + b[3]*(x-x0[3])*(y-y0[3]) + c[3]*((y-y0[3])**2)) + A[3]*(a[3]*2*(x-x0[3]) + b[3]*(y-y0[3]))*exp(a[3]*((x-x0[3])**2) + b[3]*(x-x0[3])*(y-y0[3]) + c[3]*((y-y0[3])**2))
    else:
        sys.exit("{} is an invalid choice of PES".format(pes_type))
    fx = -dV_dx - gamma*px + np.random.normal(0,std_dev)
    fy = -dV_dy - gamma*py + np.random.normal(0,std_dev)
    return fx,fy
 
def vv_step(phasepoint,dt,beta,gamma,pes_type,well=4):
    x = phasepoint[0]
    y = phasepoint[1]
    px = phasepoint[2]
    py = phasepoint[3]
    fx = phasepoint[4]
    fy = phasepoint[5]
    px = px + (1/2)*dt*fx
    py = py + (1/2)*dt*fy
    x = x + dt*px
    y = y + dt*py
    fx,fy = force(x,y,px,py,dt,beta,gamma,pes_type,well=well)
    px = px + (1/2)*dt*fx
    py = py + (1/2)*dt*fy
    return np.asarray([x,y,px,py,fx,fy])
 
def calc_op(op_type,x,y):
    if op_type == 1:
        return x
    elif op_type == 2:
        return y
    elif op_type == 3:
        return x + y
    else:
        sys.exit("Invalid choice of OP")
    
def calc_op_f(op_xtype,op_ytype,op_xcoef,op_ycoef,x,y):
    if len(op_xcoef) != op_xtype:
        sys.exit("Invalid op_xtype and op_xcoef pairing")
    if len(op_ycoef) != op_ytype:
        sys.exit("Invalid op_ytype and op_ycoef pairing")
    x_vec = [x**i for i in range(1,op_xtype+1)]
    y_vec = [y**i for i in range(1,op_ytype+1)]
    return np.dot(op_xcoef,x_vec) + np.dot(op_ycoef,y_vec)
    
def print_op(op_xtype,op_ytype,op_xcoef,op_ycoef):
    if op_xtype < 4 and op_ytype < 4:
        op_func = "OP = "
        yes = False
        for i in range(len(op_xcoef)):
            if op_xcoef[i] != 0:
                if yes:
                    if op_xcoef[i] == 1:
                        if i == 0:
                            op_func += f"+ x "
                        if i == 1:
                            op_func += f"+ x\u00b2 "
                        if i == 2:
                            op_func += f"+ x\u00b3 "
                    else:
                        if i == 0:
                            op_func += f"+ {op_xcoef[i]}x "
                        if i == 1:
                            op_func += f"+ {op_xcoef[i]}x\u00b2 "
                        if i == 2:
                            op_func += f"+ {op_xcoef[i]}x\u00b3 "
                else:
                    if op_xcoef[i] == 1:
                        if i == 0:
                            op_func += f"x "
                        if i == 1:
                            op_func += f"x\u00b2 "
                        if i == 2:
                            op_func += f"x\u00b3 "
                    else:
                        if i == 0:
                            op_func += f"{op_xcoef[i]}x "
                        if i == 1:
                            op_func += f"{op_xcoef[i]}x\u00b2 "
                        if i == 2:
                            op_func += f"{op_xcoef[i]}x\u00b3 "
                    yes = True
        for i in range(len(op_ycoef)):
            if op_ycoef[i] != 0:
                if yes:
                    if op_ycoef[i] == 1:
                        if i == 0:
                            op_func += f"+ y "
                        if i == 1:
                            op_func += f"+ y\u00b2 "
                        if i == 2:
                            op_func += f"+ y\u00b3 "
                    else:
                        if i == 0:
                            op_func += f"+ {op_ycoef[i]}y "
                        if i == 1:
                            op_func += f"+ {op_ycoef[i]}y\u00b2 "
                        if i == 2:
                            op_func += f"+ {op_ycoef[i]}y\u00b3 "
                else:
                    if op_ycoef[i] == 1:
                        if i == 0:
                            op_func += f"y "
                        if i == 1:
                            op_func += f"y\u00b2 "
                        if i == 2:
                            op_func += f"y\u00b3 "
                    else:
                        if i == 0:
                            op_func += f"{op_ycoef[i]}y "
                        if i == 1:
                            op_func += f"{op_ycoef[i]}y\u00b2 "
                        if i == 2:
                            op_func += f"{op_ycoef[i]}y\u00b3 "
                    yes = True
        print(op_func+"\n")
