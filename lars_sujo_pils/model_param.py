lou = 1025.0
L = 7.2
B = 2.75

U_ref = 2.057

scale = 0.5*lou

X_scale = 0.5*lou*L*L*L
Y_scale = 0.5*lou*B*B*B
N_scale = 0.5*lou*L*L*L*L*L

Xu_dot  =9.264*0.001*X_scale
Yv_dot  =4.571*0.01*Y_scale
Nr_dot  =1.318*0.001*N_scale
M       =4282.0
I       =9050.0


## Non-dimensional value
nd_Xu      =0.0
nd_Xuu     =-9.842*1e-3
nd_Yv      =-1.81*1e-2
nd_Yvv     =0.0
nd_Yr      =-7.538*1e-3
nd_Yrr     =0.0
nd_Yrrr    =-1.4093*1e-2
nd_Nr      =-4.4614*1e-3
nd_Nrr     =0.0
nd_Nrrr    =1.549*1e-2
nd_Nv      =-4.707*1e-3


# scale = 0.5*lou
# Xu      =0.0*scale
# Xuu     =-9.842*0.001*scale*L*L
# Yv      =-1.81*0.01*scale*B
# Yvv     =0.0*scale*B*B
# Yr      =-7.538*0.001*B
# Yrr     =0.0*scale*B*B
# Yrrr    =-1.4093*0.01*B*B*B
# Nr      =-4.4614*0.001*L
# Nrr     =0.0*scale*L*L
# Nrrr    =1.549*0.01*scale*L*L*L
# Nv      =-4.707*0.001*scale*L


scale = 0.5*lou*U_ref 

Xu      =0.0*scale/(M+Xu_dot)
Xuu     =-98/(M+Xu_dot)
Yv      =-1.81*0.01*scale*B*B/(M+Yv_dot)
Yvv     =0.0*scale*B*B/(M+Yv_dot)
Yr      =-7.538*0.001*B*B*B/(M+Yv_dot)
Yrr     =0.0*scale*B*B*B/(M+Yv_dot)
Yrrr    =-1.4093*0.01*scale*B*B*B/(M+Yv_dot)
Nr      =-4.4614*0.001*scale*L*L*L*L/(I+Nr_dot)
Nrr     =0.0*scale*L*L*L*L/(I+Nr_dot)
Nrrr    =1.549*0.01*scale*L*L*L*L/(I+Nr_dot)
Nv      =-4.707*0.001*scale*L*L*L/(I+Nr_dot)

u_F     =100/(M+Xu_dot)
v_F     =100/(M+Yv_dot)
r_F     =100*3.5/(I+Nr_dot)




print("Drag : ", Xu, Xuu, Yv, Yr, Yrrr, Nr, Nrrr, Nv)

# print("Added Mass : ", Xu_dot, Yv_dot, Nr_dot)

print("Thrust : ",u_F, v_F, r_F)