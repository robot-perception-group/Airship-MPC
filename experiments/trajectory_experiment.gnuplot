psidot=10*pi/180   # rad/s
vmax=4.0           # m/s
vmin=2.0           # m/s
wind=0.3           # m/s
alphastate="disabled"
thetastate="disabled"
phistate="disabled"
rho=90*pi/180
beta=45*pi/180
t0=0
machines=1
stationary=1
spacing_3d=20
spacing_2d=10
limitxy=40
limitz=30
offsetX=0
offsetY=0
offsetZ=0
obstacles=0

vcX=2.31
vcY=0
vcZ=2.12

min(x,y)=x<y?x:y
max(x,y)=x>y?x:y
limits(x,m,n)=max(min(x,n),m)
hardlimits(x,m,n)=(m<limits(x,m,n))<n?limits(x,m,n):1/0
dist(x,y)=sqrt(x*x+y*y)
modulo(d,q)=d-(floor(d/q)*q)



rv(u,v0,w)=v0 + cos(u)*w/psidot

vel(u,v0,w)=sqrt((psidot*rv(u,v0,w))**2 + (sin(u)*w)**2)

#airspeed(u,v0,w)=
#vel(u,rv(u,v0,w),w)+cos(u)*w
#airspeed(u,v0,w)=psidot*rv(u,v0,w)+cos(u)*w
airspeed(u,v0,w)=vel(u,v0,w)+cos(u)*w
#-cos(u)*w

accel(u,v,w)=psidot*airspeed(u,v,w)
phi(u,v,w)=asin(accel(u,v,w)/sqrt(accel(u,v,w)**2+9.81**2))
h(u,v,w)=tan(beta+phi(u,v,w))*rv(u,v,w)

fx(u,v)=v*cos(u)
fy(u,v)=v*sin(u)
fz(u,v)=h(u,v,0)

f2x(u,v,k,j)=(v+k*sin(u+j))*cos(u)
f2y(u,v,k,j)=(v+k*sin(u+j))*sin(u)
f2z(u,v,k,j)=h(u,v+k*sin(u+j),w)

#vdot_speed_integral=-sin(u) * wind  *dt (in m/s)
#vdot_angle_integral=-sin(u) * wind/psidot * du
#v=vinitial + cos(u)*wind/psidot

f3x(u,v,w)=rv(u,v,w)*cos(u)+offsetY
f3y(u,v,w)=rv(u,v,w)*sin(u)-offsetX
f3z(u,v,w)=airspeed(u,v,w)>vmax?1/0:airspeed(u,v,w)<vmin?1/0:h(u,v,w)+offsetZ

f4x(u,v,w)=f3x(u,v,w)
f4y(u,v,w)=f3y(u,v,w)
f4z(u,v,w)=airspeed(0,v,w)>vmax?1/0:airspeed(pi,v,w)<vmin?1/0:h(u,v,w)+offsetZ

f5x(u,v,w)=-f3x(u,v,w)
f5y(u,v,w)=f3y(u,v,w)
f5z(u,v,w)=airspeed(0,v+2*w/psidot,w)>vmax?1/0:airspeed(pi,v-2*w/psidot,w)<vmin?1/0:h(u,v,w)+offsetZ

# no valid solution if
#airspeed(0,v,w)=vmax & airspeed(pi,v,w)<vmin
# maximum wind: w for which
#airspeed(0,v,w)=vmax & airspeed(pi,v,w)=vmin
#psidot*(v+cos(u)*w/psidot)+cos(u)*w=vmax
#psidot*(v+1*w/psidot)+1*w=vmax
#psidot*v=vmax-2w
#psidot*v=vmin+2w
#vmax-2w=vmin+2w
#vmax-vmin=4w
#w=(vmax-vmin)/4


# find the minimum distance
#vmin=vel(u,v0,w)-w  = psidot*rv(u,r0,w) - w = psidot*(r0 - w/psidot) - w = psidot*r0 - 2w
#r0= (vmin + 2w)/psidot
#rmin=r0-w/psidot= (vmin + 2w)/psidot - w/psidot = (vmin+w)/psidot
rmin=(vmin+wind)/psidot
#print sprintf("rmin: %f",rmin)
#print sprintf("r0-min: %f",(vmin+2*wind)/psidot)
rmax=(vmax-wind)/psidot
#print sprintf("rmax: %f",rmax)


# reversing subject/wind direction max wind speed - any point in the trajectory must lay in the reverse trajectory
# this appies especially to extreme points at u=0 and pi
# radius_change_in_one_orbit = 2*w/psidot
#airspeed(0,v+2*w/psidot,w)=vmax & airspeed(pi,v-2*w/psidot,w)=vmin
#psidot*((v+2*w/psidot)+cos(u)*w/psidot)+cos(u)*w=vmax
#psidot*(v+2*w/psidot+w/psidot)+1*w=vmax
#psidot*v=vmax-4w
#psidot*v=vmin+4w
#w=(vmax-vim)/8

#airspeed(0,v+4*w/psidot,w)=vmax & airspeed(pi,v-4*w/psidot,w)=vmin
#psidot*((v+4*w/psidot)+cos(u)*w/psidot)+cos(u)*w=vmax
#psidot*(v+4*w/psidot+w/psidot)+1*w=vmax
#psidot*v=vmax-6w
#psidot*v=vmin+6w
#w=(vmax-vim)/12


persx(u,v)=0.5*sin((v-3)*pi/12)*sin(u)+offsetY
persy(u,v)=0.5*sin((v-3)*pi/12)*cos(u)-offsetX
persz(u,v)=sin(2*(v-3)*pi/12)+offsetZ

cylx(u,v,d)=0.5*d*sin(limits(4*u,0,2*pi))
cyly(u,v,d)=0.5*d*cos(limits(4*u,0,2*pi))
cylz(u,v,d)=limits(4*v,0,40)

qm0(q0,q1,q2,q3,r0,r1,r2,r3)=r0*q0-r1*q1-r2*q2-r3*q3
qm1(q0,q1,q2,q3,r0,r1,r2,r3)=r0*q1+r1*q0-r2*q3+r3*q2
qm2(q0,q1,q2,q3,r0,r1,r2,r3)=r0*q2+r1*q3+r2*q0-r3*q1
qm3(q0,q1,q2,q3,r0,r1,r2,r3)=r0*q3-r1*q2+r2*q1+r3*q0

rot0(v0,v1,v2,q0,q1,q2,q3)=qm1(\
qm0(q0,q1,q2,q3,0,v0,v1,v2),\
qm1(q0,q1,q2,q3,0,v0,v1,v2),\
qm2(q0,q1,q2,q3,0,v0,v1,v2),\
qm3(q0,q1,q2,q3,0,v0,v1,v2),\
q0,-q1,-q2,-q3)
rot1(v0,v1,v2,q0,q1,q2,q3)=qm2(\
qm0(q0,q1,q2,q3,0,v0,v1,v2),\
qm1(q0,q1,q2,q3,0,v0,v1,v2),\
qm2(q0,q1,q2,q3,0,v0,v1,v2),\
qm3(q0,q1,q2,q3,0,v0,v1,v2),\
q0,-q1,-q2,-q3)
rot2(v0,v1,v2,q0,q1,q2,q3)=qm3(\
qm0(q0,q1,q2,q3,0,v0,v1,v2),\
qm1(q0,q1,q2,q3,0,v0,v1,v2),\
qm2(q0,q1,q2,q3,0,v0,v1,v2),\
qm3(q0,q1,q2,q3,0,v0,v1,v2),\
q0,-q1,-q2,-q3)

cX(x,y,z,q0,q1,q2,q3,X,Y,Z)=x + rot0(X,Y,Z,q0,q1,q2,q3)
cY(x,y,z,q0,q1,q2,q3,X,Y,Z)=y + rot1(X,Y,Z,q0,q1,q2,q3)
cZ(x,y,z,q0,q1,q2,q3,X,Y,Z)=z + rot2(X,Y,Z,q0,q1,q2,q3)

r11(q0,q1,q2,q3)=q0*q0 + q1*q1 - q2*q2 - q3*q3
r12(q0,q1,q2,q3)=2*(q1*q2+q0*q3)
qpsi(q0,q1,q2,q3)=atan2(r12(q0,q1,q2,q3),r11(q0,q1,q2,q3))
qpsi3(q0,q1,q2,q3)=modulo(2*pi+qpsi(q0,q1,q2,q3),(2*pi))

t(x)=(x/1000000000)-t0

reference_radius=(vmin+2*wind)/psidot
