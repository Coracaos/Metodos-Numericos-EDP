import numpy as np

def metodoExplicito(c, L, T, h, k, f, a, b):
    
    """
    método numérico para resolver la ecuación del calor
    
    parametros
    ----------
    c: parametro de la ecuacion del calor
    L: longitud de la barra
    T: tiempo final
    h: paso para la partición de la longitud 
    k: paso para la partición del tiempo  
    f: funcion para los valores iniciales de temperatura
    a: funcion para la temperatura en el extremo izquiero
    b: funcion para la temperatura en el extremo derecho
    
    retorna una matriz con las soluciones de la temperatura para cada valor del tiempo
    """
    
        
    r = k*c/h**2 # calculamos el valor de r (parametro de la disusion termica) 
    
    m = round(L/h) + 1 #m puntos a lo largo de la barra
    n = round(T/k) + 1 #n puntos en el tiempo

    x = np.linspace(0, L, m) #particion en la longitud de la barra
    
    sol = np.zeros((n,m)) # inicializamos la matriz de solucion
        
    sol[0] = f(x) # valor inicial de temperatura
    
    sol[:,0] = a(0) #valor de frontera extremo derecho
    sol[:,-1] = b(L) #valor de frontera extremo izquierdo
    
    #aplicamos la formula 
    for j in range(n-1):
        for i in range(1,m-1):
            sol[j+1,i] = r*sol[j,i+1] + (1-2*r)*sol[j,i] + r*sol[j,i-1]
        
    return sol

def tridiag(a,b,c,N):
    A = np.zeros((N,N))
    
    np.fill_diagonal(A[:-1,1:],a)
    np.fill_diagonal(A,b)
    np.fill_diagonal(A[1:,:-1],c)
    
    return A


def metodoExplicitoMatriz(c, L, T, h, k, f, a, b):
    
    r = c*k/h**2
    m = round(L/h) + 1
    n = round(T/k) + 1
    
    x = np.linspace(0, L, m)
    
    sol = np.zeros((n,m))
    sol[0] = f(x)
    sol[:,0] = a(0) 
    sol[:,-1] = b(L)
    
    A = tridiag(r, 1-2*r, r, m-2) 
    
    for i in range(n-1):
        sol[i+1,1:-1] = np.dot(A, sol[i,1:-1])
        sol[i+1,1] += r*sol[i,0]
        sol[i+1,-2] += r*sol[i,-1]
        
    return sol

def metodoImplicito(c, L, T, h, k, f, a, b):
    
    r = c*k/h**2
    m = round(L/h) + 1
    n = round(T/k) + 1
    
    x = np.linspace(0, L, m)
    
    sol = np.zeros((n,m))
    sol[0] = f(x)
    sol[:,0] = a(0) 
    sol[:,-1] = b(L)
    
    A = tridiag(-r, 1+2*r, -r, m-2) 

    A_inv = np.linalg.inv(A)
    
    for i in range(n-1):
        C = sol[i,1:-1].copy() 
        C[0] += r*sol[i,0]
        C[-1] += r*sol[i,-1]
        sol[i+1,1:-1] = np.dot(A_inv, C)
        
    return sol


def metodoCrankNicolson(c, L, T, h, k, f, a, b):
    
    r = c*k/h**2
    m = round(L/h) + 1
    n = round(T/k) + 1
    
    x = np.linspace(0, L, m)
    
    sol = np.zeros((n,m))
    sol[0] = f(x)
    sol[:,0] = a(0) 
    sol[:,-1] = b(L)
    
    A = tridiag(-r/2, 1+r, -r/2, m-2)
    B = tridiag(r/2, 1-r, r/2, m-2)

         
    A_inv = np.linalg.inv(A)


    for i in range(n-1):
        
        C = np.dot(B, sol[i,1:-1])
        C[0]  += r/2*(sol[i+1,0]  + sol[i,0])
        C[-1] += r/2*(sol[i+1,-1] + sol[i,-1])

        sol[i+1,1:-1] = np.dot(A_inv, C)
    
    return sol


def metodoHiperbolica(c, L, T, h, k, f, g, a, b):
     
    r = (c*k**2)/h**2
    n = round(T/k) + 1
    m = round(L/h) + 1

    x = np.linspace(0,L,m)
    
    sol = np.zeros((n,m))

    sol[0] = f(x)
    #primera forma
    #sol[1] = f(x) + k*g(x)

    #segunda forma
    for j in range(1,m-1):
        sol[1,j] = (1-r)*f(x[j]) + (r/2)*f(x[j+1]) + (r/2)*f(x[j-1]) + k*g(x[j])

    sol[:,0] = a(0)
    sol[:,-1] = b(L)

    A = tridiag(r, 2-2*r, r, m-2)

    for i in range(1,n-1):
        sol[i+1,1:-1] = np.dot(A, sol[i,1:-1]) - sol[i-1,1:-1]
        sol[i+1,1] += r*sol[i,0]
        sol[i+1,-2] += r*sol[i,-1]
        
    return sol


def metodoEliptica(a, b, c, d, h, k, f, g_a, g_b, g_c, g_d):
    
    r = (h**2)/k**2

    n = round((b-a)/h) + 1
    m = round((d-c)/k) + 1

    x = np.linspace(a , b, n);
    y = np.linspace(c , d, m);

    p = (n-2)*(m-2)
     
    C = tridiag(r, -2*(1+r), r, n-2)

    A = np.kron(np.eye(m-2), C)
    
    np.fill_diagonal(A[:(p-(n-2)),(n-2):],r)

    np.fill_diagonal(A[(n-2):,:(p-(n-2))],r)

    xx,yy = np.meshgrid(x[1:-1],y[1:-1])
    
    B = np.zeros((m-2,n-2))
    
    for i in range(m-2):
        for j in range(n-2):
            B[i,j] = (h**2) * f(x[j+1], y[-(i+2)]) 

    for i in range(m-2):  
        B[i,0] -= g_a( y[-(i+2)] )
        B[i,-1] -= g_b( y[-(i+2)] )         

    for j in range(n-2):
        B[0,j] -= r*g_d(x[j+1])
        B[-1,j] -= r*g_c(x[j+1]) 
    
    B = np.reshape(B, np.size(B))

    W = np.linalg.solve(A, B)

    #malla de solucion

    sol = np.zeros((n,m))

    sol[0,:]  = [g_c(i) for i in x]
    sol[-1,:] = [g_d(i) for i in x]
    sol[:,0]  = [g_a(i) for i in y]
    sol[:,-1] = [g_b(i) for i in y]

    W = np.reshape(W, (m-2,n-2) )    

    sol[1:-1,1:-1] = W

    return sol


if __name__ == "__main__" :


    W0 = metodoEliptica(0, 4, 0, 4, 1, 1, lambda x, y: 0, lambda x: 80, lambda x: 0.0, lambda x: 20, lambda x: 180)
    print(W0)

    #sol = metodoHiperbolica(1, 1, 0.4, 0.05, 0.1, lambda x: 100*x**2, lambda x: 200*x, lambda x: 100*x**2, lambda x: 100*(1+x)**2)
    #print(sol)

    
    #sol1 = metodoHiperbolica(1, 1, 1, 0.25, 0.25, lambda x: np.sin(np.pi*x), lambda x: 0, lambda x: 0, lambda x:0)
    #sol2 = metodoHiperbolica(1, np.pi, 0.5, np.pi/10.0, 0.05, lambda x: np.sin(x), lambda x: 0, lambda x: 0, lambda x:0)
    #sol3 = metodoHiperbolica(1, np.pi, 0.5, np.pi/20.0, 0.1, lambda x: np.sin(x), lambda x: 0, lambda x: 0, lambda x:0)
    #sol4 = metodoHiperbolica(1, np.pi, 0.5, np.pi/20.0, 0.05, lambda x: np.sin(x), lambda x: 0, lambda x: 0, lambda x:0)
    #sol5 = metodoHiperbolica(1, 1, 0.3, 0.1, 0.1, lambda x: np.sin(2*np.pi*x), lambda x: 2*np.pi*np.sin(2*np.pi*x), lambda x: 0, lambda x:0)
    #print(sol5)
    #sol = metodoExplicito(1, 1, 0.5, 0.1, 0.01, lambda x: np.sin(np.pi*x), lambda x: 0, lambda x: 0 )
    #print(sol[-1])
    #print(metodoExplicitoMatriz(0.1, 1, 0.3, 0.2, 0.1, lambda x: 100, lambda x: 20,
    #lambda x: 40 ) )
           
    












