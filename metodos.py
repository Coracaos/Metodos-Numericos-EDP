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

if __name__ == "__main__" :
	
	print( metodoExplicito(0.1, 1, 0.3, 0.2, 0.1, lambda x: 100, lambda x: 20, lambda x: 40) )		
	print( metodoExplicitoMatriz(0.1, 1, 0.3, 0.2, 0.1, lambda x: 100, lambda x: 20, lambda x: 40) )		

