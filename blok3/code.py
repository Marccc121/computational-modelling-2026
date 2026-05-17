#Uvod u računalno modeliranje - Blok3
#Marko Herceg
#0036556780
#
#
import numpy as N
from matplotlib import pyplot as P
#
#
#3A == Gauss eliminacija s pivotingom i bez njega


import numpy as np

matrixA =np.array([
    [2.0,1.0,-1.0],
    [1.0,-1.0,3.0],
    [-1.0,2.0,1.0]
],dtype=float)

vectorB =np.array([8.0,-1.0,3.0],dtype=float)

#
def pivotGauss( M,v ):
    n =len( v )
    A_copy =M.copy()
    B_copy =v.copy( )
    
    for i in range( 0,n ):
        rowLargest =i
        currMax =abs(A_copy[ i,i ])
        
        for k in range( i+1,n ):
            if abs( A_copy[ k,i ]) >currMax:
                currMax =abs(A_copy[ k,i ])
                rowLargest =k
                
        if rowLargest!=i:
            for col in range( 0,n ):
                tempVal = A_copy[ i,col ]
                A_copy[ i,col ] =A_copy[rowLargest,col]
                A_copy[rowLargest,col] =tempVal
            
            helperB =B_copy[ i ]
            B_copy[ i ] =B_copy[ rowLargest ]
            B_copy[ rowLargest ] =helperB
            
        for j in range( i+1,n ):
            term =A_copy[ j,i ]/A_copy[ i,i ]
            for m in range( i,n ):
                A_copy[ j,m ] =A_copy[ j,m ] -term*A_copy[ i,m ]
            B_copy[ j ] =B_copy[ j ] -term*B_copy[ i ]
#
#
def naiveGauss( M,v ):
    size =len( v )
    copyA =M.copy( )
    copyB =v.copy( )
    
    for i in range( 0,size ):
        for j  in range( i+1,size ):
            factor =copyA[ j,i ]/copyA[ i,i ]
            for k in range( i,size ):
                copyA[ j,k ] =copyA[ j,k ] -factor*copyA[ i,k ]
            copyB[ j] =copyB[ j ] -factor*copyB[ i ]
            
    xSolutions =np.zeros( size )
    for i in range(size-1,-1,-1):
        totSum =0.0
        for j in range(i+1,size):
            totSum =totSum +copyA[ i,j ]*xSolutions[ j ]
        xSolutions[ i ] =(copyB[ i ] -totSum)/copyA[ i,i ]
        
    return xSolutions
#
#      
#final sub
    finalX =np.zeros( n )

    for i in range( n-1,-1, -1 ):
        topSum =0.0
        for  j in range( i+1,n ):
            topSum =topSum +A_copy[ i,j ]*finalX[ j ]
        finalX[ i ] =(B_copy[ i ] -topSum)/A_copy[ i,i ]
    return finalX



builtInAns =np.linalg.solve(matrixA,vectorB)
pivotAns =pivotGauss(matrixA,vectorB)
naiveAns =naiveGauss(matrixA,vectorB)

print("------------3A----------------")
print(f"Pivot:  {pivotAns}")
print(f"Naive: {naiveAns}")
print(f"Numpy check: {builtInAns}")
print("-----------------------------")



#3B == Demonstracija važnosti pivotinga
#
import numpy as np
matrixA2 =np.array([

    [0.0001,1.0001],
    [1.0000,2.0000]

],dtype=float)


vectorB2 =np.array([ 1.0001,2.0000 ],dtype=float)

#This is from 3a
naiveAns2 =naiveGauss(matrixA2,vectorB2)
pivotAns2 =pivotGauss(matrixA2,vectorB2)

print("-------------3B-----------------")
print(f"No pivot: {naiveAns2}")
print(f"Pivot:    {pivotAns2}")



#3C == LU faktorizacija i mreža otpornika

#import numpy as N

#
#L * y =b
def forwardSUB(L, b):
    n =len( b )
    y =np.zeros( n )
    for i in range( 0,n ):
        sumY =0.0
        for j in range( 0,i ):
            sumY =sumY + L[i,j]*y[j]
        y[i] =(b[i] - sumY)/L[i,i]

    return y

#
#U * x =y
def backSUB(U, y):
    n =len( y )
    x =np.zeros( n )
    for i in range(n-1,-1,-1):
        sumX =0.0
        for j  in range(i+1,n):
            sumX =sumX + U[i,j]*x[j]
        x[i] =(y[i] - sumX)/U[i,i]

    return x


resistorA =np.array([
    [10.0,-4.0,-2.0],
    [-4.0,8.0,-2.0],
    [-2.0,-2.0,8.0]

],dtype=float)

resistorB =np.array([ 200.0,0.0,0.0 ],dtype=float)


#
def  myLU(A):
    n =len( A )
    U =np.zeros(( n,n ))
    L = np.zeros(( n,n ))


    for i in range( 0,n ):
        L[i,i] =1.0


    for i in range( 0,n ):
        for j in range( 0,n ):
            if i <=j:
                #calc U
                sumU = 0.0
                for k in range( 0,i ):
                    sumU = sumU + L[i, k] * U[k, j]
                U[i, j] = A[i, j] - sumU
            else:
                #calc L
                sumL =0.0
                for k in range( 0,j ):
                    sumL =sumL + L[i,k]*U[k,j]

                L[i,j] =(A[i,j] - sumL)/U[j,j]

    return L,U


#run LU
matL,matU =myLU( resistorA )
vectorY =forwardSUB( matL,resistorB )
my_currents =backSUB( matU,vectorY )


condNumber =np.linalg.cond(resistorA)
#check
numpyCurrent =np.linalg.solve(resistorA,resistorB)

print("-------------3C---------------")
print(f"Lu currents (I1,I2,I3): {my_currents}")
print(f"Numpy currents:  {numpyCurrent}")
print(f"Matrix condition num:  {condNumber:.4f}")



#
#3D == Gauss-Seidel iterativna metoda

matrixGS =np.array([

    [2.0, 1.0, -1.0],
    [-1.0, 2.0, 1.0],
    [1.0, -1.0, 3.0]

], dtype=float)
vectorGS =np.array([8.0,3.0,-1.0],dtype=float)


#
def iterativeSeidel(mat,rhs,iterationsLimit=30,stopTol= 1e-5):
    errorHistory =[ ]

    totalRows =len( rhs )
    currGuessNum =np.zeros( totalRows )


    for cycle in range(1, iterationsLimit+1):
        prevGuesses =currGuessNum.copy( )

        # Goofy but effective way to bypass the diagonal element
        for  r in range( 0,totalRows ):
            runningSum =0.0
            for c in range( 0,r ):
                runningSum =runningSum + mat[r,c]*currGuessNum[ c ]
            for c in range(r+1,totalRows):
                runningSum =runningSum + mat[r,c] * currGuessNum[ c ]
            currGuessNum[r] =(rhs[ r ] - runningSum)/mat[r,r]


        biggestGap =0.0
        for idx in range(0,totalRows):
            currentDiff =abs(currGuessNum[ idx ] - prevGuesses[ idx ])
            if currentDiff >biggestGap:
                biggestGap =currentDiff

        errorHistory.append(( cycle,biggestGap ))

        if biggestGap <stopTol:
            break

    return currGuessNum, errorHistory
#

#exec
finalSol,trackingLogs =iterativeSeidel( matrixGS,vectorGS )

print("--------------3D-------------")
print(f"Gauss-seidel solutions: {finalSol}")
print(f"Tot cycles needed:  {len(trackingLogs)}")
print("--------------------------------")

P.figure(figsize =( 8,4))
P.semilogy([step[ 0 ] for step in trackingLogs],[step[ 1 ] for step in trackingLogs],"^-",color ="darkblue",
           label ="Seidel process")

P.title("Convergence behavior" )
P.grid(True,linestyle ="--",alpha =0.5)
P.xlabel("Iteration index" )
P.ylabel("Calculated err" )
P.legend( )
P.show( )


