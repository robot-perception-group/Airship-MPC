import opengen as og
import casadi.casadi as cs
import matplotlib.pyplot as plt
import numpy as np
import math

# Build parametric optimizer
# ------------------------------------

maxiterations=200
# some parameters:



W={ "TargetExclusion":1.0,
    "VehicleSoftExclusion":8.0,
    "VehicleHardExclusion":4.0,
    "costHVDE":1.0,
    "costVDIST":1.0,
    #"costUWEIGHTS":0.01,
    "costUWEIGHTS":0.0,
    "costXWEIGHTS":1.0,
    "costSOFTCONSTRAINT":0.0,
    "ForceOmega":10*math.pi/180.0,
    "ForceTheta":1.0,
    "ForcePhi":1.0,
    #"Cl":0.7,
    #"Cl":1.0,
    "Cl":0.8,
    #"Cl":10000000,
    }
#note: Tested with minVh 1.0, maxVh 4.0, omega 10°/s - no theta
#      Cl of 0.7 is barely manageable with max AoA of 15°
#      Cl of 0.8 max AoA of 10°
#      Cl of 1.0 -> 8°
#      Cl of 5.0 -> 2°


state_blimp=["Px","Py","Pz","Psi","Vh","Vvp","CamH","CamV","CamDist","Cl"]
parameters_blimp=["phase"]
#parameters_blimp=[]
#update_blimp=["PsiDot","VhDot","VvDot"]
update_blimp=["VhDot","VvpDot"]
#update_blimp=["VhDot","VvpDot","phase"]
ulimits_blimp={"PsiDot":[-25*math.pi/180.,25*math.pi/180.],"VhDot":[-5.,5.],"phase":[0,math.pi],"VvpDot":[-50.,50.]}
plimits_blimp={"phase":[0,math.pi]}
uweights_blimp={"PsiDot":1.0,"VhDot":1.0,"VvpDot":1.0}
xweights_blimp={}
#visualweights_blimp=[0.6,1.0,1.0]
visualweights_blimp=[0.0,1.0,1.0]
#visualweights_blimp=[0.0,1.0,1.0]
#shardlimits_blimp={"Vh":[0.5,5.0],"Vv":[-3.5,3.5],"Px":[-20.,20.],"Py":[-20.,20.],"Pz":[-50.,-2.]}
#ssoftlimits_blimp={"Vh":[0.9,4.6],"Vv":[-3.0,3.0],"Px":[-18.,18.],"Py":[-18.,18.],"Pz":[-48.,-4.]}
shardlimits_blimp={"Vh":[2.0,4.0]}
ssoftlimits_blimp={"Vh":[2.1,3.9]}
#ssoftlimits_blimp={"Vh":[0.7,4.0],"Vv":[-3.0,3.0],"Pz":[-45.,-7.]}

state_wind=["Vx","Vy","Vz"]
state_target=["Px","Py","Pz","Vx","Vy","Vz"]
state_obstacle=["Px","Py","Size"]
blimps=1
obstacles=0

if blimps==2:
    refangle=math.pi/2.0
elif blimps>2:
    refangle=2.0*math.pi/blimps

nu=blimps*len(update_blimp)
nfp=blimps*len(parameters_blimp)
nx=len(state_wind)+obstacles*len(state_obstacle)+blimps*len(state_blimp)+len(state_target)

N=150
Ns=150
dt=0.25

visualweights_blimp[2]=N*visualweights_blimp[2]/Ns
windforce=0.3

# ----------------- CODE STARTS HERE -----------------

# generic state vector object, allows addressing components by name and still do nifty vector math
class State(object):

    def __init__(self, dictionary):
        self._members=dictionary.keys()
        for key in self._members:
            if isinstance(dictionary[key],dict):
                self[key] = State(dictionary[key])
            else:
                self[key] = dictionary[key]

    def __str__(self):
        tmp="{"
        for key in self._members:
            tmp+=str(key)+": "+str(self[key])+", "
        return tmp+"}"

    def __add__(self, other):
        tmp={}
        if (isinstance(other,State)):
            for key in self._members:
                tmp[key]=self[key] + other[key]
        else:
            for key in self._members:
                tmp[key]=self.key + other
        return State(tmp)

    def __iadd__(self, other):
        tmp={}
        if (isinstance(other,State)):
            for key in self._members:
                self[key]+= other[key]
        else:
            for key in self._members:
                self[key]+= other
        return self

    def __mul__(self, other):
        tmp={}
        if (isinstance(other,State)):
            for key in self._members:
                tmp[key]=self[key] * other[key]
        else:
            for key in self._members:
                tmp[key]=self[key] * other
        return State(tmp)

    def __sub__(self, other):
        tmp={}
        if (isinstance(other,State)):
            for key in self._members:
                tmp[key]=self[key] - other[key]
        else:
            for key in self._members:
                tmp[key]=self[key] - other
        return State(tmp)

    def __div__(self, other):
        tmp={}
        if (isinstance(other,State)):
            for key in self._members:
                tmp[key]=self[key] / other[key]
        else:
            for key in self._members:
                tmp[key]=self[key] / other
        return State(tmp)


    def __setitem__(self, key, newvalue):
        if isinstance(newvalue,dict):
            setattr(self,str(key),State(newvalue))
        else:
            setattr(self,str(key),newvalue)

    def __getitem__(self, key):
        return getattr(self,str(key));

def RungeKutta(f,x,u,p,dt):
    k1=f(x,u,p,0)
    k2=f(x+k1*(dt/2),u,p,(dt/2))
    k3=f(x+k2*(dt/2),u,p,(dt/2))
    k4=f(x+k3*dt,u,p,dt)
    return (k1+k2*2+k3*2+k4)*(dt/6)

def initState(z0):
    xi=0
    X={}

    X["wind"]={}
    for key in state_wind:
        X["wind"][key]=z0[xi]
        xi+=1
    X["obstacle"]={}
    for i in range(obstacles):
        X["obstacle"][i]={}
        for key in state_obstacle:
            X["obstacle"][i][key]=z0[xi]
            xi+=1
    X["blimp"]={}
    for i in range(blimps):
        X["blimp"][i]={}
        for key in state_blimp:
            X["blimp"][i][key]=z0[xi]
            xi+=1
    X["target"]={}
    for key in state_target:
        X["target"][key]=z0[xi]
        xi+=1
    return State(X)

def initUpdate(u0):
    ui=0
    U={}
    U["blimp"]={}
    for i in range(blimps):
        U["blimp"][i]={}
        for key in update_blimp:
            U["blimp"][i][key]=u0[ui]
            ui+=1
    return State(U)

def initParameters(p0):
    pi=0
    P={}
    P["blimp"]={}
    for i in range(blimps):
        P["blimp"][i]={}
        for key in parameters_blimp:
            P["blimp"][i][key]=p0[pi]
            pi+=1
    return State(P)

def FDotTarget(x,dt):
    f={}
    f["Px"] = x["Vx"]
    f["Py"] = x["Vy"]
    f["Pz"] = x["Vz"]
    f["Vx"] = 0
    f["Vy"] = 0
    f["Vz"] = 0
    return State(f)

def FDotObstacle(x,dt):
    f={}
    f["Px"] = 0
    f["Py"] = 0
    f["Size"] = 0
    return State(f)

def FDotWind(x,dt):
    f={}
    f["Vx"] = 0
    f["Vy"] = 0
    f["Vz"] = 0
    return State(f)

def Velh(x,u,dt):
    return x["Vh"]+u["VhDot"]*dt

def Alpha(x,u,dt):
    v=Velh(x,u,dt)
    return (u["PsiDot"]*v)/((x["Cl"]*(v**2))-u["VhDot"])

def AlphaDot(x,u,dt):
    return (u["VhDot"]*u["PsiDot"])/(u["VhDot"]-(x["Cl"]*(Velh(x,u,dt)**2)))

def Chi(x,u,dt):
    return x["Psi"]+(dt*u["PsiDot"])-Alpha(x,u,dt)

def Velv(x,u,p,dt):
    return cs.sin(p["phase"]+Chi(x,u,dt))*(x["Vvp"] + u["VvpDot"]*dt)

def Omega(x,u,dt):
    return u["PsiDot"]-AlphaDot(x,u,dt)

def FDotBlimp(x,wind,u,p,dt):
    f={}
    v=Velh(x,u,dt)
    chi=Chi(x,u,dt)
    f["Px"] = v*cs.cos(chi) + wind["Vx"]
    f["Py"] = v*cs.sin(chi) + wind["Vy"]
    f["Pz"] = Velv(x,u,p,dt) + wind["Vz"]
    f["Psi"] = u["PsiDot"]
    f["Vh"] = u["VhDot"]
    f["Vvp"] = u["VvpDot"]
    f["CamH"] = 0
    f["CamV"] = 0
    f["CamDist"] = 0
    f["Cl"] = 0
    return State(f)

def FDot(x,u,p,dt):
    f={}
    f["wind"]=FDotWind(x["wind"],dt)
    obstacles={}
    for o in x["obstacle"]._members:
        obstacles[o]=FDotObstacle(x["obstacle"][o],dt)
    f["obstacle"]=State(obstacles)
    blimps={}
    for b in x["blimp"]._members:
        blimps[b]=FDotBlimp(x["blimp"][b],x["wind"],u["blimp"][b],p["blimp"][b],dt)
    f["blimp"]=State(blimps)
    f["target"]=FDotTarget(x["target"],dt)
    return State(f)


def rotateX(v,a):
    f={}
    f["x"]=v["x"]
    f["y"]=cs.cos(a)*v["y"]-cs.sin(a)*v["z"]
    f["z"]=cs.cos(a)*v["z"]+cs.sin(a)*v["y"]
    return State(f)

def rotateY(v,a):
    f={}
    f["x"]=cs.cos(a)*v["x"]+cs.sin(a)*v["z"]
    f["y"]=v["y"]
    f["z"]=cs.cos(a)*v["z"]-cs.sin(a)*v["x"]
    return State(f)

def rotateZ(v,a):
    f={}
    f["x"]=cs.cos(a)*v["x"]-cs.sin(a)*v["y"]
    f["y"]=cs.cos(a)*v["y"]+cs.sin(a)*v["x"]
    f["z"]=v["z"]
    return State(f)

def rotateIntoCameraFrame(xv,u,p,xt):
    global W
    g=9.81
    phi = cs.arctan2(u["PsiDot"]*Velh(xv,u,0),g)*W["ForcePhi"]
    theta = -cs.arctan2(Velv(xv,u,p,dt),xv["Vh"])*W["ForceTheta"]
    psi = xv["Psi"]
    rho = xv["CamH"]
    beta = xv["CamV"]

    pv={}
    pv["x"]=xv["Px"]
    pv["y"]=xv["Py"]
    pv["z"]=xv["Pz"]
    pv=State(pv)

    pt={}
    pt["x"]=xt["Px"]
    pt["y"]=xt["Py"]
    pt["z"]=xt["Pz"]
    pt=State(pt)

    vectorWorld = pt-pv # vector from vehicle to target

    vectorBody = rotateX(rotateY(rotateZ(vectorWorld,-psi),-theta),-phi)
    return rotateY(rotateZ(vectorBody,-rho),beta) #vectorCamera

def cameraVector(xv,u,p):
    global W
    g=9.81
    phi = cs.arctan2(u["PsiDot"]*Velh(xv,u,0),g)*W["ForcePhi"]
    theta = -cs.arctan2(Velv(xv,u,p,dt),xv["Vh"])*W["ForceTheta"]
    psi = xv["Psi"]
    rho = xv["CamH"]
    beta = xv["CamV"]

    v0=State({"x":1.,"y":0.,"z":0.})

    pv={}
    pv["x"]=xv["Px"]
    pv["y"]=xv["Py"]
    pv["z"]=xv["Pz"]
    pv=State(pv)

    vectorBody = rotateZ(rotateY(v0,-beta),rho)
    return rotateZ(rotateY(rotateX(vectorBody,phi),theta),psi)

# calculates the nearest point (projection) of the camera vector to the target - good to plot errors
def closestPoint(xv,u,p,xt):
    cv = cameraVector(xv,u,p)
    pv={}
    pv["x"]=xv["Px"]
    pv["y"]=xv["Py"]
    pv["z"]=xv["Pz"]
    pv=State(pv)

    pt={}
    pt["x"]=xt["Px"]
    pt["y"]=xt["Py"]
    pt["z"]=xt["Pz"]
    pt=State(pt)


    p0 = pv-pt # vector from vehicle to target

    l = - (cv["x"]*p0["x"]+cv["y"]*p0["y"]+cv["z"]*p0["z"])/(cv["x"]**2+cv["y"]**2+cv["z"]**2)

    return pv + (cv * l)

def costHVDE(x,u,p,t):
        global visualweights_blimp

        tv = rotateIntoCameraFrame(x,u,p,t)
        # this forms an elypsoid with it's axis defined by the main weights - around [camdist,0,0] 
        return ((tv["x"]-x["CamDist"])*visualweights_blimp[0])**2 + (tv["y"]*visualweights_blimp[1])**2 + (tv["z"]*visualweights_blimp[2])**2

def costHVDE2(x,u,p,t):
        global visualweights_blimp

        tv = rotateIntoCameraFrame(x,u,p,t)
        # this forms an elypsoid with it's axis defined by the main weights - around [camdist,0,0] 
        return ((tv["x"]-x["CamDist"])*visualweights_blimp[0])**2 + (tv["y"]*visualweights_blimp[1])**2


u = cs.SX.sym('u', nu*N+nfp)
z0 = cs.SX.sym('z0', nx)

X=initState(z0)

cost = 0
stateconstraint = 0
P = initParameters(u[nu*N:])
for t in range(0, N):
    ui0=t*nu

    U = initUpdate(u[ui0:])

    # force omega
    for i in range(blimps):
        U["blimp"][i]["PsiDot"]=W["ForceOmega"]

    X+=RungeKutta(FDot,X,U,P,dt)

    #cost
    for i in range(blimps):
        # constraints
        #   speed
        for key in shardlimits_blimp:
            stateconstraint += cs.fmax(0,  -(X["blimp"][i][key]-shardlimits_blimp[key][0]))
            stateconstraint += cs.fmax(0,  (X["blimp"][i][key]-shardlimits_blimp[key][1]))

        costSOFTCONSTRAINT = 0
        for key in ssoftlimits_blimp:
            costSOFTCONSTRAINT += cs.fmax(0,  -(X["blimp"][i][key]-ssoftlimits_blimp[key][0]))*1.0/abs(shardlimits_blimp[key][0]-ssoftlimits_blimp[key][0])
            costSOFTCONSTRAINT += cs.fmax(0,  (X["blimp"][i][key]-ssoftlimits_blimp[key][1]))*1.0/abs(shardlimits_blimp[key][1]-ssoftlimits_blimp[key][1])

        #   distance to target
        tvX = X["target"]["Px"] - X["blimp"][i]["Px"]
        tvY = X["target"]["Py"] - X["blimp"][i]["Py"]
        tvXY= cs.sqrt(tvX**2+tvY**2)
        stateconstraint += cs.fmax(0, (W["VehicleHardExclusion"]+W["TargetExclusion"]) - tvXY)**2
        costSOFTCONSTRAINT += (cs.fmax(0, (W["VehicleSoftExclusion"]+W["TargetExclusion"]) - tvXY)**2)*1.0/(abs(W["VehicleHardExclusion"]-W["VehicleSoftExclusion"])**2)

        # don't collide with obstacles
        for j in range(obstacles):
            vvX = X["blimp"][i]["Px"] - X["obstacle"][j]["Px"]
            vvY = X["blimp"][i]["Py"] - X["obstacle"][j]["Py"]
            # constraint - minimum distance between vehicles
            vvXY=cs.sqrt(vvX**2 +  vvY**2)
            stateconstraint += cs.fmax(0, (W["VehicleHardExclusion"]+X["obstacle"][j]["Size"]) - vvXY)**2
            costSOFTCONSTRAINT += (cs.fmax(0, (W["VehicleSoftExclusion"]+X["obstacle"][j]["Size"]) - vvXY)**2)*1.0/(abs(W["VehicleHardExclusion"]-W["VehicleSoftExclusion"])**2)

        costVDIST = 0
        for j in range(blimps):
            if i!=j:
                vvX = X["blimp"][i]["Px"] - X["blimp"][j]["Px"]
                vvY = X["blimp"][i]["Py"] - X["blimp"][j]["Py"]
                vvZ = X["blimp"][i]["Pz"] - X["blimp"][j]["Pz"]
                # constraint - minimum distance between vehicles
                vvXY=cs.sqrt(vvX**2 +  vvY**2)
                stateconstraint += cs.fmax(0, (2.0*W["VehicleHardExclusion"]) - vvXY)**2
                costSOFTCONSTRAINT += (cs.fmax(0, (2.0*W["VehicleSoftExclusion"]) - vvXY)**2)*1.0/(abs(2*(W["VehicleHardExclusion"]-W["VehicleSoftExclusion"]))**2)

                # cost for relative angle to target
                tv2X= X["target"]["Px"] - X["blimp"][j]["Px"]
                tv2Y= X["target"]["Py"] - X["blimp"][j]["Py"]

                angle = cs.acos((tvX*tv2X+tvY*tv2Y)/(tvXY*cs.sqrt(tv2X**2+tv2Y**2)))
                costVDIST +=  cs.sinh(angle-refangle)**2

        costU = 0
        for key in uweights_blimp:
            costU += (uweights_blimp[key] * U["blimp"][i][key])**2
        costX = 0
        for key in xweights_blimp:
            costX += (xweights_blimp[key] * X["blimp"][i][key])**2

        if t<Ns:
            cost += W["costHVDE"]*costHVDE(X["blimp"][i],U["blimp"][i],P["blimp"][i],X["target"])
        else:
            cost += W["costHVDE"]*costHVDE2(X["blimp"][i],U["blimp"][i],P["blimp"][i],X["target"])
        cost += W["costVDIST"]*costVDIST
        cost += W["costUWEIGHTS"]*costU
        cost += W["costXWEIGHTS"]*costX
        cost += W["costSOFTCONSTRAINT"]*costSOFTCONSTRAINT
        #cost += 0.2* ( (cs.sqrt((X["blimp"][i]["Px"])**2+(X["blimp"][i]["Py"])**2) + X["blimp"][i]["Pz"])**2 )
        #cost += (cs.atan2(X["blimp"][i]["Py"],X["blimp"][i]["Px"])-X["blimp"][i]["Psi"])**2
        #cost +=  (-cs.cos(X["blimp"][i]["Psi"])*3+15+X["blimp"][i]["Pz"])**2


# no final cost actually

umin=[]
umax=[]
for i in range(blimps):
    for j in update_blimp:
        umin.append(ulimits_blimp[j][0])
        umax.append(ulimits_blimp[j][1])
umin = umin * N
umax = umax * N
for i in range(blimps):
    for j in parameters_blimp:
        umin.append(plimits_blimp[j][0])
        umax.append(plimits_blimp[j][1])
bounds = og.constraints.Rectangle(umin, umax)
#statebounds = og.constraints.BallInf(None,1e-3)
#statebounds = og.constraints.Zero()

#    .with_aug_lagrangian_constraints(c, statebounds)
problem = og.builder.Problem(u, z0, cost).with_constraints(bounds)\
    .with_penalty_constraints(stateconstraint)
build_config = og.config.BuildConfiguration()\
    .with_build_directory("python_test_build")\
    .with_build_mode("debug")\
    .with_tcp_interface_config()
meta = og.config.OptimizerMeta()\
    .with_optimizer_name("navigation")
solver_config = og.config.SolverConfiguration()\
        .with_tolerance(1e-3)\
        .with_initial_tolerance(1e-2)\
        .with_max_outer_iterations(20)\
        .with_delta_tolerance(1e-2)\
        .with_penalty_weight_update_factor(10.0)\
        .with_initial_penalty(100.0)
builder = og.builder.OpEnOptimizerBuilder(problem, 
                                          meta,
                                          build_config, 
                                          solver_config) \
    .with_verbosity_level(1)
builder.build()


# Use TCP server
# ------------------------------------
mng = og.tcp.OptimizerTcpManager('python_test_build/navigation')
mng.start()



# initial state - this should come from ROS:
x_init=[]

# wind
#x_init += [0.7,0.,0.]
x_init += [windforce,0.,0.]

# no obstacles for now
#x_init += [0.0,-10.0,10.0]
 
# blimps
x_init += [12.0,0.,-15., math.pi/2.0, 3.0, 0., 87*math.pi/180.,45*math.pi/180.,25.0,W["Cl"]]
#x_init += [-20.,-20.,-5.0, 0., 1.0, 0., 90*math.pi/180.,45*math.pi/180.,25.0,10000]
#x_init += [-20.,20.,-5.0, 0., 1.0, 0., 90*math.pi/180.,45*math.pi/180.,25.0,1]
#x_init += [1,20.,-5.0, 0., 1.0, 0., 90*math.pi/180.,45*math.pi/180.,25.0,1]

# target

x_init += [0.,0.,0.,0.,0.,0.]

# init
x_states = [0.0] * (nx*(N+2))
x_states[0:nx] = x_init
initial_guess=[1e-6] * (nu*N) + [1e-6]*nfp
cc= [0.0] * (blimps*(N+2))

plt.ion()
for ITERATION in range(maxiterations):
    x_states[0:nx] = x_init
    mng.ping()
    response = mng.call(x_init, initial_guess=initial_guess)
    #response = mng.call(x_init)
    data=response.get()

    if (not response.is_ok()):
        print(data.message)
        exit(1)
    #print(data.solve_time_ms)
    #print(data.cost)
    #print(data.penalty)
    #print(data.exit_status)

    plt.clf()
    # Plot solution
    # ------------------------------------
    time = np.arange(0, dt*N, dt)
    u_star = data.solution

    P=initParameters(u_star[nu*N:])

    saveVV = blimps*[0.0]
    #plt.subplot(911)
    #for b in range(blimps):
    #    u = u_star[3*b+0:nu*N:nu]
    #    plt.plot(time, u, '-o')
    #plt.ylabel('PsiDot')
    #plt.xlabel('Time')

    #plt.subplot(912)
    #for b in range(blimps):
    #    u = u_star[3*b+1:nu*N:nu]
    #    plt.plot(time, u, '-o')
    #plt.ylabel('VhDot')
    #plt.xlabel('Time')

    #plt.subplot(913)
    #for b in range(blimps):
    #    u = u_star[3*b+2:nu*N:nu]
    #    plt.plot(time, u, '-o')
    #plt.ylabel('VvDot')
    #plt.xlabel('Time')

    # Plot trajectory
    # ------------------------------------

    # simulate for display
    for t in range(0, N):
        #u_t = u_star[t*nu:(t+1)*nu]

        # assign and simulate
        i0=t*nx
        u0=t*nu
        ic=t*blimps

        x=initState(x_states[i0:])
        U=initUpdate(u_star[u0:])

        # force omega
        for i in range(blimps):
            U["blimp"][i]["PsiDot"]=W["ForceOmega"]

        x+=RungeKutta(FDot,x,U,P,dt)
        for i in range(blimps):
            costVDIST = 0
            tvX = x["target"]["Px"] - x["blimp"][i]["Px"]
            tvY = x["target"]["Py"] - x["blimp"][i]["Py"]
            tvXY= cs.sqrt(tvX**2+tvY**2)
            for j in range(blimps):
                if i!=j:
                    # cost for relative angle to target
                    tv2X= x["target"]["Px"] - x["blimp"][j]["Px"]
                    tv2Y= x["target"]["Py"] - x["blimp"][j]["Py"]

                    angle = cs.acos((tvX*tv2X+tvY*tv2Y)/(tvXY*cs.sqrt(tv2X**2+tv2Y**2)))
                    costVDIST +=  cs.sinh(angle-refangle)**2

            costU = 0
            for key in uweights_blimp:
                costU += (uweights_blimp[key] * U["blimp"][i][key])**2
            costX = 0
            for key in xweights_blimp:
                costX += (xweights_blimp[key] * x["blimp"][i][key])**2

            cost = 0
            if t<Ns:
                cost += W["costHVDE"]*costHVDE(x["blimp"][i],U["blimp"][i],P["blimp"][i],x["target"])
            else:
                cost += W["costHVDE"]*costHVDE2(x["blimp"][i],U["blimp"][i],P["blimp"][i],x["target"])
            #print("costHVDE",cost)
            cost += W["costVDIST"]*float(costVDIST)
            #print("costVDIST",W["costVDIST"]*float(costVDIST))
            cost += W["costUWEIGHTS"]*costU
            #print("costUWEIGHTS",W["costUWEIGHTS"]*costU)
            cost += W["costXWEIGHTS"]*costX
            cc[ic]=float(cost)
            #print("=====",cost)
            ic+=1
            #vec = rotateIntoCameraFrame(x["blimp"][i],U["blimp"][i],x["target"])
            #print(vec)
            if (t<10):
                #theta = -cs.arctan2(x["blimp"][i]["Vv"],x["blimp"][i]["Vh"])
                #print(math.sqrt(x["blimp"][i]["Px"]**2+x["blimp"][i]["Py"]**2),x["blimp"][i]["Vh"],theta*180.0/math.pi)
                cv=cameraVector(x["blimp"][i],U["blimp"][i],P["blimp"][i])
                cp=closestPoint(x["blimp"][i],U["blimp"][i],P["blimp"][i],x["target"])
                a=Alpha(x["blimp"][i],U["blimp"][i],0)*180.0/math.pi
                Vv=Velv(x["blimp"][i],U["blimp"][i],P["blimp"][i],0)


                print("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f"%(x["blimp"][i]["Px"],x["blimp"][i]["Py"],x["blimp"][i]["Pz"],cv["x"],cv["y"],cv["z"],cp["x"],cp["y"],cp["z"],x["blimp"][i]["Vh"],Vv,a))



                #print((Chi(x["blimp"][i],U["blimp"][i],0)-x["blimp"][i]["Psi"])*180.0/math.pi)
                Xdot=FDotBlimp(x["blimp"][i],x["wind"],U["blimp"][i],P["blimp"][i],0)
                dist=math.sqrt(x["blimp"][i]["Px"]**2+x["blimp"][i]["Py"]**2)
                pdot=math.sqrt(Xdot["Px"]**2+Xdot["Py"]**2)
                phi=math.atan2(x["blimp"][i]["Px"],x["blimp"][i]["Py"])
                r0_exp=(shardlimits_blimp["Vh"][0]+2*windforce)/W["ForceOmega"]
                rdot=math.sin(phi)*windforce
                rExp=r0_exp+math.cos(phi)*windforce/W["ForceOmega"]
                vExp=math.sqrt((dist*U["blimp"][i]["PsiDot"])**2+rdot**2)
                delta=abs(vExp-pdot)
#                print("Absolute Velocity:",pdot)
#                print("Expected Velocity:",vExp)
#                print("Error:",delta)
#                print("Distance:",dist,rExp)
                #delta=x["blimp"][i]["Vh"]-vh
                #-math.sqrt(x["wind"]["Vx"]**2+x["wind"]["Vy"]**2))
                #print(delta,U["blimp"][i]["PsiDot"]*180.0/math.pi)
        #        print(rs,pdot,x["blimp"][i]["Vh"])
            elif (t==10):
                saveVV[i] = Velv(x["blimp"][i],U["blimp"][i],P["blimp"][i],0)

        # store
        ix=i0+nx
        for key in state_wind:
            x_states[ix]=x["wind"][key]
            ix+=1
        for i in range(obstacles):
            for key in state_obstacle:
                x_states[ix]=x["obstacle"][i][key]
                ix+=1
        for i in range(blimps):
            for key in state_blimp:
                x_states[ix]=x["blimp"][i][key]
                ix+=1
        for key in state_target:
            x_states[ix]=x["target"][key]
            ix+=1


    xx = x_states[0:nx*N:nx]
    xy = x_states[1:nx*N:nx]

    #print(x_states)
    #print(xx)
    plt.subplot(311)
    #target
    plt.ylabel('Position')
    plt.xlabel('Time')

    plt.gca().set_aspect('equal',adjustable='box')
    plt.gca().set_xlim(-50,50)
    plt.gca().set_ylim(-50,50)
    #blimps
    for i in range(blimps):
        xx = x_states[len(state_wind)+obstacles*len(state_obstacle)+i*len(state_blimp):nx*N:nx]
        xy = x_states[len(state_wind)+obstacles*len(state_obstacle)+i*len(state_blimp)+1:nx*N:nx]
        plt.plot(xx, xy, '-o')
    xx = x_states[len(state_wind)+obstacles*len(state_obstacle)+blimps*len(state_blimp):nx*N:nx]
    xy = x_states[len(state_wind)+obstacles*len(state_obstacle)+blimps*len(state_blimp)+1:nx*N:nx]
    plt.plot(xx, xy, '-o')


    plt.subplot(312)
    plt.ylabel('Altitude')
    plt.gca().set_ylim(0,-50)
    plt.xlabel('Time')
    for i in range(blimps):
        xx = x_states[len(state_wind)+obstacles*len(state_obstacle)+i*len(state_blimp)+2:nx*N:nx]
        plt.plot(time, xx,'-o')
    xx = x_states[len(state_wind)+obstacles*len(state_obstacle)+blimps*len(state_blimp)+2:nx*N:nx]
    plt.plot(time, xx,'-o')

    #plt.subplot(916)
    #plt.ylabel('Chi')
    #plt.xlabel('Time')
    #for i in range(blimps):
    #    xx = x_states[len(state_wind)+obstacles*len(state_obstacle)+i*len(state_blimp)+3:nx*N:nx]
    #    plt.plot(time, xx,'-o')

    #plt.subplot(917)
    #plt.ylabel('Vh')
    #plt.xlabel('Time')
    #for i in range(blimps):
    #    xx = x_states[len(state_wind)+obstacles*len(state_obstacle)+i*len(state_blimp)+4:nx*N:nx]
    #    plt.plot(time, xx,'-o')

    #plt.subplot(918)
    #plt.ylabel('Vv')
    #plt.xlabel('Time')
    #for i in range(blimps):
    #    xx = x_states[len(state_wind)+obstacles*len(state_obstacle)+i*len(state_blimp)+5:nx*N:nx]
    #    plt.plot(time, xx,'-o')

    plt.subplot(313)
    plt.ylabel('COST')
    plt.xlabel('Time')
    plt.gca().set_ylim(0,5)
    for i in range(blimps):
        xx = cc[i:N*blimps:blimps]
        plt.plot(time, xx,'-o')



    x_init=x_states[10*nx:11*nx]
    # simulate loss of periodic information:
    for i in range(blimps):
        #x_init[len(state_wind)+obstacles*len(state_obstacle)+i*len(state_blimp)+5]=saveVV[i] # put all vertical in the linear part
        #x_init[len(state_wind)+obstacles*len(state_obstacle)+i*len(state_blimp)+6]=0.0 # and none in the periodic (because we can't ask the vehicle for it
        #print(saveVV[i],"sucks")
        pass

    #initial_guess = u_star[10*nu:nu*N] + ([1e-6]*(nu*10))
    plt.pause(0.01)

mng.kill()
