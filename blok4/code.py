#Uvod u računalno modeliranje - Blok4
#Marko Herceg
#0036556780
#
#
#
import numpy as np
from matplotlib import pyplot as P

#Deklaracija
subY =np.array([ 2.1,13.6,61.1,150.7,247.1 ],dtype =float)
fullY =np.array([ 2.1,7.7,13.6,27.2,40.9,61.1,82.3,111.5,150.7,193.6,247.1 ],dtype =float)

subX =np.array([ 0,2,5,8,10 ],dtype =float)
fullX =np.array([ 0,1,2,3,4,5,6,7,8,9,10 ],dtype =float)

#T = HEAT VALUES
#V = VOLT VALUES
V =np.array([ 1.05,1.42,1.98,2.53,3.01,3.62,4.09,4.71,5.12 ],dtype =float)
T =np.array([ 20,30,40,50,60,70,80,90,100 ],dtype=float)
#
#
#4C == Interpolacijski polinomi
#
#
def calcLagrange(dataX,dataY, seekX):
    totalNodes =len( dataX )
    resultY =0.0

    for k in range( 0,totalNodes ):
        runningWeight =dataY[ k ]
        for m in range( 0,totalNodes ):
            if m ==k:
                continue

            numerator =seekX -dataX[ m ]
            denominator =dataX[ k ] -dataX[ m ]
            runningWeight =runningWeight *(numerator/denominator)
        resultY =resultY +runningWeight

    return resultY


def calcNewton(nodesX, nodesY, seekX):
    coefs =[y for y in nodesY]
    numPoints =len( nodesX )


    for step in range(1,numPoints):
        for i in range(numPoints-1,step-1,-1):
            spanX =nodesX[ i ] -nodesX[ i-step ]
            coefs[ i ] =(coefs[ i ] -coefs[i-1])/spanX

    totalVal =coefs[numPoints -1]
    for i in  range(numPoints -2,-1,-1):
        totalVal =coefs[ i ] +(seekX -nodesX[ i ])*totalVal

    return totalVal
#


#
#4B == Usporedba regresijskih modela
#
#
def checkFitMetrics(actual, fitted):
    meanVal = 0.0
    count =len( actual )


    for i in range( 0,count):
        meanVal =meanVal +actual[ i ]
    meanVal =meanVal/count

    ssTot =0.0
    ssRes =0.0
    for i in range( 0,count ):
        ssTot =ssTot +(actual[i] -meanVal) *(actual[i] -meanVal)
        ssRes =ssRes +(actual[i] -fitted[i]) *(actual[i] -fitted[i])

    r2 =1.0 -(ssRes/ssTot)
    rmse =np.sqrt( ssRes/count )
    return r2,rmse

#

coefQuad =np.polyfit( fullX,fullY,2 )
predQuad =np.polyval( coefQuad,fullX )

coefCube =np.polyfit( fullX,fullY,3 )
predCube =np.polyval( coefCube,fullX )

coefLin =np.polyfit( fullX,fullY,1 )
predLin =np.polyval( coefLin,fullX )


weightX2 =0.0
weightLogY =0.0
weightLogXY =0.0
totWeights =0.0
weightX =0.0


for i in range(0,len( fullX )):
    lnY = np.log(fullY[i])
    w =fullY[ i ]

    totWeights =totWeights + w
    weightX =weightX + w * fullX[i]
    weightX2 =weightX2 + w * fullX[i] * fullX[i]
    weightLogY =weightLogY + w * lnY
    weightLogXY =weightLogXY + w * fullX[i] * lnY


detM =totWeights *weightX2 -weightX *weightX
expGrowth =(totWeights *weightLogXY -weightX *weightLogY)/detM
expLogA =(weightX2 *weightLogY -weightX *weightLogXY)/detM

realA =np.exp( expLogA )
predExp =np.zeros(len( fullX ))

for i in  range(0, len( fullX )):
    predExp[ i ] =realA*np.exp(expGrowth*fullX[ i ])

r2Lin,rmseLin =checkFitMetrics( fullY,predLin )
r2Quad,rmseQuad =checkFitMetrics( fullY,predQuad )
r2Cube,rmseCube =checkFitMetrics( fullY,predCube )
r2Exp,rmseExp =checkFitMetrics( fullY,predExp )

print("-----------------4B----------------")
print(f"Linear fit: r2 ={r2Lin:.4f}, RMSE ={rmseLin:.2f}")
print(f"Quadratic fit: r2 ={r2Quad:.4f}, RMSE ={rmseQuad:.2f}")
print(f"Cubic fit:  r2 ={r2Cube:.4f}, RMSE = {rmseCube:.2f}")
print(f"Exponential fit: r2 ={r2Exp:.4f}, RMSE ={rmseExp:.2f}")
print("----------------------------------------")
#

#
#4A == Linearna regresija od nule
#
#
totItems =len( T )
sumT =0.0
sumV =0.0
sumTV =0.0
sumT2 =0.0

for k in range(0, totItems):
    sumT = sumT + T[k]
    sumV = sumV + V[k]
    sumTV = sumTV + T[k]*V[k]
    sumT2 = sumT2 + T[k]*T[k]

coefSlope =(totItems*sumTV -sumT*sumV)/(totItems*sumT2 -sumT*sumT)
avgV =sumV/totItems
avgT =sumT/totItems

coefIntercept =avgV -coefSlope*avgT


linePredictions =np.zeros( totItems )
for i in range( 0,totItems ):
    linePredictions[ i ] = coefIntercept+coefSlope*T[ i ]


topPart =0.0
bottomPart =0.0
for i in range( 0,totItems ):
    topPart =topPart +(V[i] -linePredictions[ i ]) * (V[i] -linePredictions[ i ])
    bottomPart =bottomPart +(V[i] -avgV) * (V[ i ] - avgV)

finalR2 =1.0 -(topPart/bottomPart)

print("----------------4A------------------")
print(f"Slope coefSlope: {coefSlope:.4f}")
print(f"Intercept coefIntercept: {coefIntercept:.4f}")
print(f"R2 Score: {finalR2:.4f}")
print("----------------------------------")
#
#
#
#
#
plotGrid =np.linspace( 0,10,100 )
lagrangeCurve =[calcLagrange( subX,subY,pt ) for pt in plotGrid]
newtonCurve = [calcNewton( subX,subY,pt ) for pt in plotGrid]
#
#graph #1
P.title("4A linear fit check")
P.figure(figsize=( 5,3.5 ))
P.scatter(T,V,color ="red",label ="Experiment data")
P.plot(T,linePredictions,color ="darkred",label ="Calculatted regression")
P.legend( )
#
#graph #2
P.title("4A residuals plot")
residuals_4A =np.zeros( totItems )
for i in range( 0,totItems ):
    residuals_4A[ i ] =V[ i ] -linePredictions[ i ]

P.figure(figsize =( 5,2.5 ))
P.axhline(y=0.0,color ="darkgrey",linestyle ="--")
P.scatter(T,residuals_4A,color ="darkred",marker ="x")
#


#graph 4b
P.title("4B: regression models comparison")
P.figure(figsize=( 6,4 ))
P.scatter(fullX,fullY,zorder =4,label ="Raw data",color ="black")
P.plot(fullX,predLin,label ="Linear model",linestyle =":")
P.plot(fullX,predCube,label ="Cube model",linestyle ="-.")
P.plot(fullX,predQuad,label ="Quad model",linestyle ="--")
P.plot(fullX,predExp,label ="Exponen model",color ="pink")
P.legend( )

#graph 4c
P.title("4C: Interpolation v regression")
P.figure(figsize =( 6,4 ))
P.scatter(subX,subY,color="blue",s=75,zorder =5,label ="5 point subset")
P.plot(fullX, predQuad,label ="4B quad regression reference",linestyle =":",color ="darkgrey")
P.plot(plotGrid,lagrangeCurve,label="Lagrrange method",color ="yellow")
P.plot(plotGrid,newtonCurve,label="Newton method",color ="pink",linestyle ="--")
P.legend( )
#
#
#
P.show( )