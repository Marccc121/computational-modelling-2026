#Uvod u računalno modeliranje - Blok1
#Marko Herceg
#0036556780
#
#
import numpy as N
from matplotlib import pyplot as P
#
#
#1A == Modeliranje pada tijela s otporom zraka


g =9.81#m/s2
c_d =0.25#kg/m
m =70#kg
t_min =0; t_max =20#s
v_zero =0#m/s
t = N.linspace(t_min,t_max,200)
v: int

dv_dt =N.sqrt(g *m/c_d) * N.tanh(N.sqrt(g *c_d/m) *t)

x_values_Anal =t
y_values_Anal =dv_dt



dvdt_current =0
P.figure( 1 )
P.figure( 2 )
#
for i in range( 0,3 ):
    dt =( 2**i ) *0.5
    n_points =int(t_max/dt +1)
    t_list =N.zeros(int(t_max/dt +1)) #11/21/41
    v_list =N.zeros(int(t_max/dt +1))
    #
    for j in range(n_points -1):
        dvdt_current =g -( c_d/m )*v_list[ j ]**2

        t_list[ j+1 ] =t_list[ j ] +dt
        v_list[ j+1 ] =(dvdt_current *dt) +v_list[j]

    P.figure( 1 )
    P.plot(t_list,v_list)

    #errors
    dvdt_euler =N.sqrt( g *m/c_d ) * N.tanh(N.sqrt( g *c_d/m ) *t_list)
    dvdt_error =dvdt_euler -dvdt_current
    abs_error =N.abs(dvdt_euler -v_list)
    print('For dt = ',round(dt, 3),' absolute error =',round(abs_error[ -1 ],3))
    rel_error =N.abs(dvdt_euler -v_list) / N.abs(dvdt_euler)
    print('For dt =',round(dt, 3),' relative error =',round(rel_error[ -1 ],3))

    P.figure( 2 )
    P.plot(t_list, abs_error, label =f'dt = {dt}')
    P.title('Apsolutna pogreska kroz vrijeme')
    P.legend()
#

P.figure( 1 )
P.plot(x_values_Anal,y_values_Anal,color ='black')
P.xlabel( 'x-Axis' )
P.ylabel( 'y-Axis' )
P.title('Pad tijela s otporom zraka')
P.legend( )

P.show( )


#1B == Derivacija numerička


h_list =[1,0.5,0.1,0.05,0.01,0.005,0.001,0.0001,0.00001]
central_error =N.zeros( 9 )
forward_error =N.zeros( 9 )


x =1
correct_val =N.cos( x )
h_list_num =9
#
for i in range(0,h_list_num):
    h_current =h_list[ i ]

    central =(N.sin(x + h_current) - N.sin(x -h_current))/(2 *h_current)
    forward =(N.sin(x + h_current) - N.sin(x))/h_current

    central_error[ i ] =N.abs(correct_val -central) #similar to above 1A
    forward_error[ i ] =N.abs(correct_val -forward) #-II-
#

P.figure( 3 )
P.loglog(h_list, forward_error, label ='Forward')
P.loglog(h_list,central_error, label ='Central')
P.title('Err in numerical differentiation')
P.grid(True, which ="both")
P.legend( )

P.show()



#1C == Gubitak značajnih znamenki (roundoff)

x =0
y =[1e-1,1e-2,1e-4,1e-6,1e-8,1e-10]

print(f"{'x' :<10} #{'2 * sin^2(x/2)' :<20} # {'1 - cos(x)' :<20} # {'Razlika' :<20}")
for i in range(0,75):
    print("#", end ='')
print( "#" )
#
for x in y:
    Right = 2 *(N.sin( x/2 ) **2)
    Left =1 - N.cos( x )
    diff =N.abs(Left -Right) #Again same as above 1A,1B

    print(f"{Left :<20.18f} # {Right :<20.18f} # {x:<10} # {diff :<20.18e}")
#