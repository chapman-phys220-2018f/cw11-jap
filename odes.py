#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
# Name: Alley Busick
# Student ID: 2293544
# Email: busick@chapman.edu
# Course: PHYS220/MATH220/CPSC220 Fall 2018
# Assignment: CW11
###

import numpy as np
import matplotlib.pyplot as plt
    
def euler(N):
    """ Plots solution with inital conditions 
    x(0) = 1, v(0) = 0 and x'(t) = v(t) v'(t) = -x'(t)
    verifying x(t) = cos(t) and v(t) = -sin(t) using Euler's Method 
    """
    t = np.linspace(0, 5*2*np.pi, N) #domain
    v = np.zeros(N)
    v[0] = 0
    x = np.zeros(N)
    x[0] = 1
    val = (5*2*np.pi)/N
    
    for i in range(1, len(v)):
        v[i] = v[i-1] + val*x[i-1]
        x[i] = x[i-1] - val*v[i-1]
    
    plt.plot(t, v, color="orange", label="v(t)")
    plt.plot(t, np.sin(t), color="red", linestyle="--", label="sin(t)")
    plt.plot(t, x, color="cyan", label="x(t)")
    plt.plot(t, np.cos(t), color="blue", linestyle="--", label="cos(t)")
    plt.title("Euler's Method")
    plt.legend(["v(t)", "sin(t)", "x(t)", "cos(t)"])
    
def heun(N):
    """ Plots solution with inital conditions 
    x(0) = 1, v(0) = 0 and x'(t) = v(t) v'(t) = -x'(t)
    verifying x(t) = cos(t) and v(t) = -sin(t) using Heun's Method 
    """
    t = np.linspace(0, 5*2*np.pi, N) #domain
    v = np.zeros(N)
    v[0] = 0
    x = np.zeros(N)
    x[0] = 1
    val = (5*2*np.pi)/N
    
    for i in range(1, len(t)):
        v1 = x[i-1] - val*v[i-1]
        x1 = v[i-1] + val*x[i-1]
        v[i] = v[i-1] + (val/2)*(x[i-1] + v1)
        x[i] = x[i-1] - (val/2)*(v[i-1] + x1)
    
    plt.plot(t, v, color="orange", label="v(t)")
    plt.plot(t, np.sin(t), color="red", linestyle="--", label="sin(t)")
    plt.plot(t, x, color="cyan", label="x(t)")
    plt.plot(t, np.cos(t), color="blue", linestyle="--", label="cos(t)")
    plt.title("Heun's Method")
    plt.legend(["v(t)", "sin(t)", "x(t)", "cos(t)"])
    
def rk2(N):
    """ Plots solution with inital conditions 
    x(0) = 1, v(0) = 0 and x'(t) = v(t) v'(t) = -x'(t)
    verifying x(t) = cos(t) and v(t) = -sin(t) using 2nd-order Runge Kutta method
    """
    t = np.linspace(0, 5*2*np.pi, N) #domain
    v = np.zeros(N)
    v[0] = 0
    x = np.zeros(N)
    x[0] = 1
    val = (5*2*np.pi)/N
    
    #2nd order
    x1 = 0
    y1 = 0
    x2 = 0
    y2 = 0
    
    for i in range(1, len(v)):
        x1 = val*v[i-1]
        x2 = val*(x[i-1]+x1/2)
        
        y1 = -val*v[i-1]
        y2 = -val*(v[i-1] + y1/2)
        
        v[i] = v[i-1] + x2
        x[i] = x[i-1] + y2
    
    plt.plot(t, v, color="orange", label="v(t)")
    plt.plot(t, np.sin(t), color="red", linestyle="--", label="sin(t)")
    plt.plot(t, x, color="cyan", label="x(t)")
    plt.plot(t, np.cos(t), color="blue", linestyle="--", label="cos(t)")
    plt.title("2nd-order Runge Kutta Method")
    plt.legend(["v(t)", "sin(t)", "x(t)", "cos(t)"])
    
def rk4(N):
    """ Plots solution with inital conditions 
    x(0) = 1, v(0) = 0 and x'(t) = v(t) v'(t) = -x'(t)
    verifying x(t) = cos(t) and v(t) = -sin(t) using 4th-order Runge Kutta method
    """
    t = np.linspace(0, 5*2*np.pi, N) #domain
    v = np.zeros(N)
    v[0] = 0
    x = np.zeros(N)
    x[0] = 1
    val = (5*2*np.pi)/N
    
    #4th order
    x1 = 0
    y1 = 0
    x2 = 0
    y2 = 0
    x3 = 0
    y3 = 0
    x4 = 0
    y4 = 0
    
    for i in range(1, len(v)):
        x1 = val*x[i-1]
        x2 = val*(x[i-1] + x1/2)
        x3 = val*(x[i-1] + x2/2)
        x4 = val*(x[i-1] + x3)
        
        y1 = -val*v[i-1]
        y2 = -val*(v[i-1] + y1/2)
        y3 = -val*(v[i-1] + y2/2)
        y4 = -val*(v[i-1] + y3)
        
        x[i] = x[i-1] + (y1 + (y2*2) + (y3*2) + y4)/6
        v[i] = v[i-1] + (x1 + (x2*2) + (x3*2) + x4)/6
    
    plt.plot(t, v, color="orange", label="v(t)")
    plt.plot(t, np.sin(t), color="red", linestyle="--", label="sin(t)")
    plt.plot(t, x, color="cyan", label="x(t)")
    plt.plot(t, np.cos(t), color="blue", linestyle="--", label="cos(t)")
    plt.title("4th-order Runge Kutta Method")
    plt.legend(["v(t)", "sin(t)", "x(t)", "cos(t)"])
    