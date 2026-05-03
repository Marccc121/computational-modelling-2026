#Uvod u računalno modeliranje - Blok2
#Marko Herceg
#0036556780
#
#
import numpy as N
from matplotlib import pyplot as P
#
#
#2A == Implementacija četiri metode za traženje nultočaka

#
def df( x ):
    return (3 * x**2) - (12 * x) + 11
def f( x ):
    return x**3 - (6 * x**2) + (11 * x) - 6
#


#Bisection
def bisectMethd( Xl,Xu,tol,max_it ):
    log =[ ]
    for i in range(1,max_it + 1):
        Xr =( Xl+Xu )/2
        error = abs(f( Xr ))
        log.append((i,Xr,error))
        
        if error<tol or ( Xu-Xl )/2<tol:
            break
        if f( Xl )*f( Xr )<0:
            Xu =Xr
        else:
            Xl =Xr
    return Xr,log

#False position
def falsePos( Xl,Xu,toli,maxI ):
    history =[ ]
    for k  in range( 1,maxI+1 ):

        Xr =Xu -(f( Xu )*( Xl-Xu ))/(f( Xl ) -f( Xu ))
        err =abs(f( Xr ))
        history.append(( k,Xr,err ))
        
        if err <toli:
            break
        if f( Xl )*f( Xr )<0:
            Xu =Xr
        else:
            Xl =Xr
    return Xr,history

#Newton raphson
def newtonRaph( start,t,m ):
    array =[ ]
    Xstari =start
    for j in range( 1,m+1 ):
        Xnovi =Xstari -f( Xstari )/df( Xstari )
        difference =abs(Xnovi-Xstari)
        array.append((j,Xnovi,difference))
        if difference <t:
            break
        Xstari =Xnovi
    return Xnovi,array

#Secant 
def secantAlg(X0,X1,eps, max_k):
    steps =[ ]
    for n in range( 1,max_k+1 ):
        f0 =f( X0 )
        f1 = f( X1 )
        X2 =X1-f1*( X1-X0 )/( f1-f0 )
        err =abs( X2-X1 )
        steps.append( ( n,X2,err ) )
        if err <eps:
            break
        X0 =X1
        X1 =X2
    return X2,steps

#
#Main
#
# 1. Graf funkcije
it_limit =40
eps =1e-5

rez_B,log_B =bisectMethd(0.1,1.9,eps,it_limit)
rez_F,log_F =falsePos(0.5,1.5,eps,it_limit)
rez_N,log_N =newtonRaph(0.5,eps,it_limit)
rez_S,log_S =secantAlg(0.5,0.8,eps,it_limit)

r1, _ =newtonRaph(0.5, eps, it_limit)
r2, _ =newtonRaph(1.7, eps, it_limit)
r3, _ = newtonRaph(3.2, eps, it_limit)

#
print("-------------2A---------------------")
print(f"Bisekcija rjesenje: { rez_B:.6f}")
print(f"False position rjesenje: { rez_F:.6f}")
print(f"Newton rjesenje: { rez_N:.6f}")
print(f"Sekanta rjeseje: { rez_S:.6f}")

print(f"1. nultocka (Metoda 1): {rez_B:.6f}") 
print(f"2. nultocka (Metoda 2): {r2:.6f}")  
print(f"3. nultocka (Metoda 3): {r3:.6f}")
print("--------------------------------------")
#

x_osa =N.linspace( -1,5,100 )
y_osa =f( x_osa )

P.figure(figsize =( 7, 4 ))
P.plot(x_osa,y_osa,label='f(x)')
P.axhline(0, color='red',linestyle='--')

P.plot(r1,0,'ko',label =f'Nultocka (x=1)')
P.plot(r2,0,'ro', label ='Nultocka (x=2)') 
P.plot(r3,0,'go',label ='Nultocka (x=3)') 

P.grid( True )
P.legend()
P.show( )


#convergence 
P.figure( figsize =( 7,4 ))

P.semilogy([p[ 0 ] for p in log_N], [p[ 2 ] for p in log_N],'^-',label ='Newton')
P.semilogy([p[ 0 ] for p in log_S], [p[ 2 ] for p in log_S], 'x-',label ='Sekanta')
P.semilogy([p[ 0 ] for p in log_B], [p[ 2 ] for p in log_B],'o-',label ='Bisekcija',linewidth=5,markersize=10,alpha=0.3)
P.semilogy([p[ 0 ] for p in log_F], [p[ 2 ] for p in log_F],'s-',label ='False pozicija',linewidth=2)

P.xlabel('Broj iteracija')
P.ylabel('err')
P.legend( )
P.title('Usporedba konvergencije')
P.show( )




#2B == Inženjerski primjer: struja u tranzistorskom krugu
#
#
#Početna jednadzba:
#I · e^(0.5·I) - 5 = 0
#f'(I) = e^(0.5I) + I · (0.5 · e^(0.5I))
#f'(I) = e^(0.5I) · (1 + 0.5I)

limit =50
eps_struja =1e-6

def df( I ):
    return N.exp(0.5*I) * (1 +0.5*I)

def f( I ):
    return I * N.exp(0.5*I) -5


sol_1,log_1 =newtonRaph(1.0, eps_struja, limit)
sol_0,log_0 =newtonRaph(0.0, eps_struja, limit)

print("----------------2B------------------")
print(f"Za x0 =1.0   Rješenje: {sol_1:.6f} A (Iteracija: {len(log_1)})")
print(f"Za x0 =0.0   Rješenje: {sol_0:.6f} A (Iteracije: {len(log_0)})")


#7. (za zadatak 2B):
#Krivulja se glatko širi bez ikakvih lokalnih vrhova/udubljenja, 
#stoga je pristup neosjetljiv na odabir x0 (0/1).
#Kako se x0 =1 vizualno približava nultoj točki, 
#primjećuje se da dvije početne vrijednosti brže konvergiraju prema istom rješenju
#,oko 2.37A, što potvrđuje da čvrsta početna aproksimacija smanjuje broj potrebnih ponavljanja.



#2C == Minimizacija jedne varijable
#
#f(x) = x^4 - 3x^3 + 2

def minFunc( x ):
    return x**4-3*x**3+2

# 9. metoda zlatnog reza
def goldSec(funk,lower,upper,tolerancija =0.000001):
    faktor =(1+N.sqrt( 5 ))/2
    ostatak =2-faktor
    
    p1 =lower +ostatak*( upper-lower )
    p2 =upper -ostatak*( upper-lower )
    
    v1 =funk( p1 )
    vrijednost2 =funk( p2 )
    
    historyMetode =[ ]
    for count in range( 1,100):
        sredina =( lower+upper )/2
        difference =abs( upper-lower )
        historyMetode.append((count,sredina,difference))
        
        if  difference <tolerancija:
            break
            
        if v1 <vrijednost2:
            upper =p2
            p2 =p1
            vrijednost2 =v1
            p1 =lower +ostatak*( upper-lower )
            v1 =funk( p1 )
        else:
            lower =p1
            p1 =p2
            v1 =vrijednost2
            p2 =upper-ostatak*( upper-lower )
            vrijednost2 =funk( p2 )
            
    return (lower+upper)/2, historyMetode

def threePparabola( fff,xA,xB,xC,tolic =0.000001,maxit =50):
    lista=[]
    for m in range( 1,maxit+1):
        up =( xB-xA )**2*(fff( xB )-fff( xC ))-(xB-xC)**2*(fff( xB )-fff( xA ))
        down =( xB-xA )*(fff( xB )-fff( xC ))-(xB-xC )*(fff( xB )-fff( xA ))
        
        if down ==0:
            break 
        
        xNext =xB -0.5*( up/down )
        errVal =abs( xNext-xB )
        lista.append((m,xNext,errVal))
        
        if errVal <tolic:
            break
        
        xA =xB
        xB =xNext
        
    return xNext,lista

#run
rezZlato,putanjaZlato =goldSec(minFunc,0,3 )
rezPara,putanjaPara =threePparabola(minFunc,0,1.5,3)

print("-----------------2C-----------------")
print(f"Zlatni rez x: {rezZlato:.6f}")
print(f"Parabola x: { rezPara:.6f}")

#
xGraph =N.linspace( -1,3.5,100 )
yOS_crtez =minFunc( xGraph )

P.figure(figsize =( 7,4 ))
P.plot(xGraph,yOS_crtez,label ='f(x) =x^4 - 3x^3 + 2',color ='blue',linewidth =2.4)

minimum_X =2.25
minimum_Y =minFunc(minimum_X)
P.plot( minimum_X,minimum_Y, 'ro', label=f'Min (x={minimum_X})' ) # Crvena točka

P.title( 'Prikaz funkcije za minimizaciju' )
P.xlabel( 'x os' )
P.ylabel( 'y os' )
P.grid(True,linestyle ='--',alpha =0.5) 
P.legend( )
P.show( )
#



#2D == Inženjerska optimizacija: cilindricni spremnik

def surfaceCyl( r ):
    V =1.0 
    return 2*N.pi*r**2+2*V/r


r_min,logValjak=goldSec( surfaceCyl,0.1,2,1e-6 )

h_valjka=1/( 3.14159*r_min**2 )

print("-----------------2D-----------------")
print(f"r ={ r_min:.4f}")
print(f"h ={ h_valjka:.4f}")


