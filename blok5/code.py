#Uvod u računalno modeliranje - Blok5
#Marko Herceg
#0036556780
#
#
#

from matplotlib import pyplot as P
import numpy as np
from scipy.integrate import solve_ivp

#5A == Numerička integracija
#
#
xCoords =np.array([ 0,0.5,1.0,1.5,2.0,2.5,3.0 ],dtype =float)
F =np.array([ 0,0.8,2.3,3.7,5.2,6.8,8.9 ],dtype =float)
#
def TRAPEZOID(func, lowerA, upperB, intervalsN):
    stepH =(upperB -lowerA)/intervalsN
    trapSum =0.5 *(func( lowerA ) +func( upperB ))

    for i in range( 1,intervalsN ):
        currX =lowerA +i *stepH
        trapSum =trapSum +func( currX )
    return trapSum *stepH


def SIMPSON(func, lowerA, upperB, intervalsN):
    stepH =(upperB -lowerA)/intervalsN
    simpSum =func( lowerA ) +func( upperB )


    for i in range(1,intervalsN):
        currX =lowerA +i *stepH
        if i %2 ==0:
            simpSum =simpSum +2.0 *func( currX )
        else:
            simpSum =simpSum +4.0 *func( currX )

    return simpSum *(stepH/3.0)

#Calculations + simul.:
testFUNC =lambda x:np.sin( x )
truIntegral =2.0
nValues =[ 2,4,8,16,32,64,128 ]

trapErr =[]
simpsErr =[]

for n in nValues:
    resT =TRAPEZOID( testFUNC,0.0,np.pi,n )
    resS =SIMPSON( testFUNC,0.0,np.pi,n )
    trapErr.append( np.abs( resT -truIntegral ))
    simpsErr.append( np.abs( resS -truIntegral ))

workMade =0.0
for i in range( 0,len( xCoords ) -1 ):
    baseH =xCoords[ i+1 ] -xCoords[ i ]
    avgForce =0.5 *(F[ i ] +F[ i+1 ])
    workMade =workMade +baseH *avgForce

print("-----------------5A----------------")
print(f"Calculated Work from lab data: {workMade:.4f} J")
print("----------------------------------------")
#


#
#5B == Numerička derivacija višeg reda
#
#
t =np.array([ 0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0 ],dtype =float)
V =np.array([ 0,0.98,1.87,2.66,3.29,3.76,4.05,4.19,4.20 ],dtype =float)
#
def diffFORWARD(func, evalPt, stepH):
    return (func( evalPt +stepH ) -func( evalPt ))/stepH


def diffCENTRAL(func, evalPt, stepH):
    return (func( evalPt +stepH ) -func( evalPt -stepH ))/(2.0 *stepH)


def diffHIGH_ACC(func, evalPt, stepH):
    term1 =-3.0 *func( evalPt )
    term2 =4.0 *func( evalPt +stepH )
    term3 =-1.0 *func( evalPt +2.0 *stepH )
    return (term1 +term2 +term3)/(2.0 *stepH)

#Calculations + simul.:
#Derivation check
baseFunc =lambda x:np.exp( x )
pointX =2.0
testH =0.1

valF =diffFORWARD( baseFunc,pointX,testH )
valC =diffCENTRAL( baseFunc,pointX,testH )
valH =diffHIGH_ACC( baseFunc,pointX,testH )

dataPoints =len( t )
A =np.zeros( dataPoints )
deltaH =0.25

for idx in range( 0,dataPoints ):
    if idx ==0:
        A[ idx ] =(V[ idx+1 ] -V[ idx ])/deltaH
    elif idx ==dataPoints -1:
        A[ idx ] =(V[ idx ] -V[ idx-1 ])/deltaH
    else:
        A[ idx ] =(V[ idx+1 ] -V[ idx-1 ])/(2.0 *deltaH)

#


#
#5C == ODE solveri: Euler i RK4 step agregacije
#
#
def runEULER_STEP(func, t, stateY, stepH):
    return stateY +stepH *func( t,stateY )


def runRK4_STEP(func, t, stateY, stepH):
    k1 =func( t,stateY )
    k2 =func( t +0.5 *stepH,stateY +0.5 *stepH *k1 )
    k3 =func( t +0.5 *stepH,stateY +0.5 *stepH *k2 )
    k4 =func( t +stepH,stateY +stepH *k3 )
    return stateY +(stepH/6.0) *(k1 +2.0 *k2 +2.0 *k3 +k4)

testODE =lambda t,y:-2.0 *y
hSteps =[ 1.0,0.5,0.25,0.1,0.05 ]
eulerErr =[ ]
rk4Err =[ ]


#Calculations + simul.:
#Global err
for h in hSteps:

    currY_E =1.0
    t_E =0.0
    while t_E <2.0 -1e-9:
        currY_E =runEULER_STEP( testODE,t_E,currY_E,h )
        t_E =t_E +h


    currY_R =1.0
    t_R =0.0
    while t_R <2.0 -1e-9:

        currY_R =runRK4_STEP( testODE,t_R,currY_R,h )
        t_R =t_R +h

    exactY =np.exp( -4.0 )
    eulerErr.append( np.abs( currY_E -exactY ))
    rk4Err.append( np.abs( currY_R -exactY ))
#


#
#5D ==Model serijskog RLC kruga
#
#
def getRLC_DERIVS(t, stateVector, R_val, L_val, C_val, V_val):
    chargeQ =stateVector[ 0 ]
    currentI =stateVector[ 1 ]

    dqdt =currentI
    didt =(1.0/L_val) *(V_val -R_val *currentI -chargeQ/C_val)
    return np.array([ dqdt,didt ])

#

#imulacija kruga +varijacija otpora
paramSource_C =15.0
timeLimit =10.0
paramInduction_L =3.0
paramCapacitor_C =2.0


scipyTarget =lambda t,y:getRLC_DERIVS( t,y,3.0,paramInduction_L,paramCapacitor_C,paramSource_C )

res ={ }
resistances =[ 1.0,3.0,10.0 ]
fixedH =0.02
timeGrid =np.arange( 0,timeLimit +fixedH,fixedH )

for R in resistances:
    qHistory =[ ]
    iHistory =[ ]
    currentStates =np.array([ 0.0,0.0 ])


    for tInstant in timeGrid:

        qHistory.append( currentStates[ 0 ] )
        iHistory.append( currentStates[ 1] )
        odeSystem =lambda t,y:getRLC_DERIVS( t,y,R,paramInduction_L,paramCapacitor_C,paramSource_C )
        currentStates =runRK4_STEP( odeSystem,tInstant,currentStates,fixedH )

    res[ R ] ={"q":np.array( qHistory ),"i":np.array( iHistory )}

scipyOutput =solve_ivp( scipyTarget,[ 0,timeLimit ],[ 0.0,0.0 ],method="RK45",rtol=1e-6 )
#
#
#
#
#
#graph #1
P.figure(figsize=( 5,4 ))
P.title("5A: Integration methods error slope")

P.loglog(nValues,trapErr,marker ="o",label ="Composite Trapezoid")
P.loglog(nValues,simpsErr,marker ="s",label ="Composite Simpson")
P.legend( )
#
#graph #2
fig,ax1 = P.subplots(figsize =( 5,3.5 ))
P.title("5B:velocity and acceleration plots")

ax1.plot(t,V, color ="darkblue",marker ="o",label ="Velocity v(t)")
ax2 =ax1.twinx( )
ax2.plot(t,A,color ="darkred",marker ="x",linestyle ="--",label ="Calc. acceleration a(t)")
#
#graph #3
P.figure(figsize=( 5,4 ))
P.title("5C:solver  convergence rate analysis")

P.loglog(hSteps,eulerErr,marker ="o",color ="purple",label ="Euler method")
P.loglog(hSteps,rk4Err,marker ="v",color ="darkblue",label ="Rk4 method")
P.legend( )


#graph 4
r3Data =res[ 3.0 ]
vcTrace =r3Data[ "q" ]/paramCapacitor_C
energyTrace =0.5 *paramInduction_L *r3Data[ "i" ]**2 +0.5 *paramCapacitor_C *vcTrace**2

P.figure(figsize=( 6,4 ))
P.title("5D:Rlc internal dynamics(R =3 Ohms)")

P.plot(timeGrid,r3Data[ "i" ],label ="Current[A]")
P.plot(timeGrid,energyTrace,label ="Tot energy[J]", linestyle ="-.")
P.plot(timeGrid,vcTrace,label ="Capacitor volt[V]")
P.legend( )


#graph 5
P.figure(figsize=( 6,4 ))
P.title("5D:current responses acros damping stages")

P.plot(timeGrid,res[ 1.0 ][ "i" ],label ="R=1 (Underdamped)")
P.plot(timeGrid,res[ 3.0 ][ "i" ],label ="R=3 (Critically Damped)")
P.plot(timeGrid,res[ 10.0 ][ "i" ],label ="R=10 (Overdamped)")
P.legend( )


#graph 6
P.figure(figsize=( 6,4 ))
P.title("5E:fixed step vs adaptive solver audit")

P.plot(timeGrid,res[ 3.0 ][ "i" ],label ="Custom Fixed RK4 (h=0.02)",color ="orange")
P.scatter(scipyOutput.t,scipyOutput.y[ 1 ],color ="black",facecolors ="none",label ="Adaptive scipy solve_ivp",zorder =3)
P.legend( )
#
#
#
P.show( )