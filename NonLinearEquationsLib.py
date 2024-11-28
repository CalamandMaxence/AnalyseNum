#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30

@author: Simone Deparis
"""

import numpy as np

def plotPhi (a,b,phi,label='$\phi$', N = 100) :
    '''
    simple plot of fonction with bisectrice of first quadrant
    useful for preparing grpahics of fixed point iterations
    [a,b] : interval, used for both x and y axis
    phi : funciton to plot
    label : label for the fonction, usually '$\phi$'
    N : number of points for plotting
    '''
    import matplotlib.pyplot as plt

    z = np.linspace(a,b,N)
    plt.plot(z,phi(z),'k-')
    plt.plot([a,b],[a,b],':', linewidth=0.5)
    
    plt.xlabel('x'); plt.ylabel(label);
    # Plot the x,y-axis 
    plt.plot([a,b], [0,0], 'k-',linewidth=0.1)
    plt.plot([0,0], [a,b], 'k-',linewidth=0.1)
    plt.legend([label, 'y=x'])
    plt.title('Graph de la fonction ' + label)


def plotPhiIterations (x) :
    # plot the graphical interpretation of the Fixed Point method
    import matplotlib.pyplot as plt

    plt.plot([x[0],x[0]], [0,x[1]], 'g:')
    
    for k in range(x.size-2) :
        plt.plot([x[k],x[k+1]], [x[k+1],x[k+1]], 'g:')
        plt.plot([x[k+1],x[k+1]], [x[k+1],x[k+2]], 'g:')
        k=x.size-2
        plt.plot([x[k],x[k+1]], [x[k+1],x[k+1]], 'g:')

    # Putting a sign at the initial point
    deltaAxis = plt.gca().axes.get_ylim()[1] - plt.gca().axes.get_ylim()[0]
    if (x[1] < 0) :
        plt.annotate("$x_0$", (x[0], 0.02*deltaAxis) )
    else :
        plt.annotate("$x_0$", (x[0], -0.05*deltaAxis) )


def plotNewtonIterations (a,b,f,x,N=200) :
    # plot the graphical interpretation of the Newton method
    import matplotlib.pyplot as plt

    z = np.linspace(a,b,N)
    plt.plot(z,f(z), 'b-', x,f(x), 'rx')

    # Putting a sign at the initial point
    deltaAxis = plt.gca().axes.get_ylim()[1] - plt.gca().axes.get_ylim()[0]    
    
    # plot the graphical interpretation of the Newton method
    plt.plot([x[0],x[0]], [0,f(x[0])], 'g:')
    for k in range(x.size-1) :
        plt.plot([x[k],x[k+1]], [f(x[k]),0], 'g-')
        plt.plot([x[k+1],x[k+1]], [0,f(x[k+1])], 'g:')

        # Putting a sign at the initial point
        if (f(x[k]) < 0) :
            plt.annotate("$x_"+str(k)+"$", (x[k], 0.02*deltaAxis) )
        else :
            plt.annotate("$x_"+str(k)+"$", (x[k], -0.05*deltaAxis) )


    plt.ylabel('$f$'); plt.xlabel('$x$');

    # Plot the x,y-axis 
    plt.plot([a,b], [0,0], 'k-',linewidth=0.1)
    plt.plot([0,0], [np.min(f(z)),np.max(f(z))], 'k-',linewidth=0.1)

    plt.legend(['$f$','($x_k$,$f(x_k)$)'])

    
def plotBisectionIterations (a,b,f,x,N=200) :
    # plot the graphical interpretation of the Bisection method
    import matplotlib.pyplot as plt

    def putSign(y,f,text) :
        if (f(y) < 0) :
            plt.annotate(text, (y, 0.02*deltaAxis) )
        else :
            plt.annotate(text, (y, -0.05*deltaAxis) )
       
    
    z = np.linspace(a,b,N)
    plt.plot(z,f(z), 'b-')
    
    plt.plot([a,a], [0,f(a)], 'g:')
    plt.plot([b,b], [0,f(b)], 'g:')

    # For putting a sign at the initial point
    deltaAxis = plt.gca().axes.get_ylim()[1] - plt.gca().axes.get_ylim()[0]

    putSign(a,f,"a")
    putSign(b,f,"b")
    
    for k in range(x.size) :
        plt.plot([x[k],x[k]], [0, f(x[k])], 'g:')
        putSign(x[k],f,"$x_"+str(k)+"$")

    # Plot the x,y-axis 
    plt.plot([a,b], [0,0], 'k-',linewidth=0.1)
    plt.plot([0,0], [np.min(f(z)),np.max(f(z))], 'k-',linewidth=0.1)

    plt.ylabel('$f$'); plt.xlabel('$x$');

    return





