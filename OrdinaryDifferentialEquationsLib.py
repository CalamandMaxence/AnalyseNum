import numpy as np

def forwardEuler( fun, interval, y0, N ) :
    #  FORWARDEULER Solve differential equations using the forward Euler method.
    #  [T, U] = FORWARDEULER( FUN, INTERVAL, Y0, N ), with INTERVAL = [T0, TF],
    #  integrates the system of differential equations y'=f(t, y) from time T0 
    #  to time TF, with initial condition Y0, using the forward Euler method on 
    #  an equispaced grid of N intervals. Function FUN(T, Y) must return
    #  a column vector corresponding to f(t, y). Each row in the solution 
    #  array U corresponds to a time returned in the column vector T.
    import numpy as np

    # time step
    h = ( interval[1] - interval[0] ) / N
  
    # time snapshots 
    t = np.linspace( interval[0], interval[1], N+1 )

    # initialize the solution vector
    u = np.zeros(N+1)
    u[0] = y0

    # time loop (n=0,...,n, but array indeces in Matlab start at 1)
    for n in range(N) :
        u[n+1] = u[n] + h * fun( t[n], u[n] )
  
    return t, u


def backwardEuler( fun, interval, y0, N ) :
    #  BACKWARDEULER Solve differential equations using the backward Euler method.
    #  [T, U] = BACKWARDEULER( FUN, INTERVAL, Y0, N ), with INTERVAL = [T0, TF],
    #  integrates the system of differential equations y'=f(t, y) from time T0 
    #  to time TF, with initial condition Y0, using the backward Euler method on 
    #  an equispaced grid of N intervals. Function FUN(T, Y) must return
    #  a column vector corresponding to f(t, y). Each row in the solution 
    #  array U corresponds to a time returned in the column vector T.
    import numpy as np
    from scipy.optimize import fsolve
    
    # time step
    h = ( interval[1] - interval[0] ) / N
  
    # time snapshots 
    t = np.linspace( interval[0], interval[1], N+1 )

    # initialize the solution vector
    u = np.zeros(N+1)
    u[0] = y0

    # time loop (n=0,...,n, but array indeces in Matlab start at 1)
    for n in range(N) :
        # non-linear function
        F = lambda x : u[n] + h * fun(t[n+1],x) - x 
        # solve the non-linear equation using the built-in matlab function "fsolve"
        # to compute u[n+1]  
        u[n+1] = fsolve( F, u[n]);
        # u[n+1] = fsolve( F, u[n]+ h * fun( t[n], u[n] ));

        # NOTE:
        # in the call of fsolve, a more accurate initial guess is obtained 
        # by replacing u[n] with the forward euler method:
        # u[n+1] = fsolve( F, u(n) + h * fun( t(n), u(n) ), options ); 

    return t, u

def Heun( fun, interval, y0, N ) :
    #  Heun Solve differential equations using the Heun method.
    #  [T, U] = Heun( FUN, INTERVAL, Y0, N ), with INTERVAL = [T0, TF],
    #  integrates the system of differential equations y'=f(t, y) from time T0 
    #  to time TF, with initial condition Y0, using the backward Euler method on 
    #  an equispaced grid of N intervals. Function FUN(T, Y) must return
    #  a column vector corresponding to f(t, y). Each row in the solution 
    #  array U corresponds to a time returned in the column vector T.
    import numpy as np
    
    # time step
    h = ( interval[1] - interval[0] ) / N
  
    # time snapshots 
    t = np.linspace( interval[0], interval[1], N+1 )

    # initialize the solution vector
    u = np.zeros(N+1)
    u[0] = y0

    # time loop (n=0,...,n, but array indeces in Matlab start at 1)
    for n in range(N) :
        u[n+1] = u[n] + h/2 * (fun( t[n], u[n] ) + fun( t[n+1], u[n]+h*fun( t[n], u[n] ) ) )

    return t, u


def CrankNicolson( fun, interval, y0, N ) :
    #  CrankNicolson Solve differential equations using the Crank-Nicolson method.
    #  [T, U] = CrankNicolson( FUN, INTERVAL, Y0, N ), with INTERVAL = [T0, TF],
    #  integrates the system of differential equations y'=f(t, y) from time T0 
    #  to time TF, with initial condition Y0, using the backward Euler method on 
    #  an equispaced grid of N intervals. Function FUN(T, Y) must return
    #  a column vector corresponding to f(t, y). Each row in the solution 
    #  array U corresponds to a time returned in the column vector T.
    import numpy as np
    from scipy.optimize import fsolve
    
    # time step
    h = ( interval[1] - interval[0] ) / N
  
    # time snapshots 
    t = np.linspace( interval[0], interval[1], N+1 )

    # initialize the solution vector
    u = np.zeros(N+1)
    u[0] = y0

    # time loop (n=0,...,n, but array indeces in Matlab start at 1)
    for n in range(N) :
        # non-linear function
        F = lambda x : u[n] + h/2 * ( fun(t[n],u[n]) + fun(t[n+1],x) )  - x 
        # solve the non-linear equation using the built-in matlab function "fsolve"
        # to compute u[n+1]  
        u[n+1] = fsolve( F, u[n]);

        # NOTE:
        # in the call of fsolve, a more accurate initial guess is obtained 
        # by replacing u[n] with the forward euler method:
        # u[n+1] = fsolve( F, u(n) + h * fun( t(n), u(n) ), options ); 

    return t, u

def modifiedEuler( fun, interval, y0, N ) :
    #  modifiedEuler Solve differential equations using the modified Euler method.
    #  [T, U] = modifiedEuler( FUN, INTERVAL, Y0, N ), with INTERVAL = [T0, TF],
    #  integrates the system of differential equations y'=f(t, y) from time T0 
    #  to time TF, with initial condition Y0, using the forward Euler method on 
    #  an equispaced grid of N intervals. Function FUN(T, Y) must return
    #  a column vector corresponding to f(t, y). Each row in the solution 
    #  array U corresponds to a time returned in the column vector T.
    import numpy as np

    # time step
    h = ( interval[1] - interval[0] ) / N
  
    # time snapshots 
    t = np.linspace( interval[0], interval[1], N+1 )

    # initialize the solution vector
    u = np.zeros(N+1)
    u[0] = y0

    # time loop (n=0,...,n, but array indeces in Matlab start at 1)
    for n in range(N) :
        u[n+1] = u[n] + h * fun( t[n]+h/2, u[n] + h/2 * fun( t[n], u[n] ) )
  
    return t, u


