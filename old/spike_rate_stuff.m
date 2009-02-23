dt=.005;
3*dt
x = log(-log(-(.015-1))/dt)

Nsec = .8
x = .88
p = 1-exp(-exp(x)*dt)
p/dt/Nsec

Nsp = (1-exp(-exp(x)*dt))/(dt*Nsec)
y=1
k=1.1
x = log(-log(-(.2-1))/dt)/k
x = log(-log(-(.8-1))/dt)/k