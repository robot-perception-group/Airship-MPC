#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#

import opengen as og
import casadi.casadi as cs
import matplotlib.pyplot as plt
import numpy as np
import math

import blimp_nmpc_wrapper_node.nodes.formation_config as FormationConfig

W=FormationConfig.W
W_blimp=FormationConfig.W_blimp
blimps=FormationConfig.blimps
obstacles=FormationConfig.obstacles
state_blimp_cost=FormationConfig.state_blimp_const
N=FormationConfig.N
dt=FormationConfig.dt
# Build parametric optimizer
# ------------------------------------


# some parameters:

state_blimp_dynamic=["Px","Py","Pz","Psi","PsiDot","Vh","Vvl","Vvp"]
state_blimp = state_blimp_dynamic
parameters_blimp=["phase"]
update_blimp=["PsiDot","VhDot","VvlDot","VvpDot"]
ulimits_blimp={"PsiDot":[-W["PsiDotMax"],W["PsiDotMax"]],"VhDot":[W["VhDotMin"],W["VhDotMax"]],"VvlDot":[W["VvlDotMin"],W["VvlDotMax"]],"VvpDot":[W["VvpDotMin"],W["VvpDotMax"]]}
plimits_blimp={"phase":[0,math.pi]}
uweights_blimp={"PsiDot":1.0,"VhDot":1.0,"VvlDot":1.0,"VvpDot":1.0}
xweights_blimp={"Vvl":1.0}
pweights_blimp={"phase":1.0}
shardlimits_blimp={"Vh":["VhMinHard","VhMaxHard"],"Vv":["VvMinHard","VhMaxHard"],"Pz":["PzMinHard","PzMaxHard"]}
ssoftlimits_blimp={"Vh":["VhMinSoft","VhMaxSoft"],"Vv":["VvMinSoft","VhMaxSoft"],"Pz":["PzMinSoft","PzMaxSoft"]}

state_wind=["Vx","Vy","Vz"]
state_target=["Px","Py","Pz","Vx","Vy","Vz"]
state_obstacle=["Px","Py","Size"]

if blimps==2:
    refangle=math.pi/2.0
elif blimps>2:
    refangle=2.0*math.pi/blimps

nu=blimps*len(update_blimp)
nfp=blimps*len(parameters_blimp)
nx=len(state_wind)+obstacles*len(state_obstacle)+blimps*len(state_blimp)+len(state_target)

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
        # scaling to limits works here, so we can alter ulimits at runtime
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

def FDotConst(x,dt):
    f={}
    for key in x._members:
        f[key] = 0
    return State(f)

def Velh(x,u,dt):
    return x["Vh"]+u["VhDot"]*dt

def Alpha(x,u,dt,Wb):
    v=Velh(x,u,dt)
    return (u["PsiDot"]*v)/((Wb["Cl"]*(v**2))-u["VhDot"])

def AlphaDot(x,u,dt,Wb):
    return (u["VhDot"]*u["PsiDot"])/(u["VhDot"]-(Wb["Cl"]*(Velh(x,u,dt)**2)))

def Chi(x,u,dt,Wb):
    return x["Psi"]+(dt*u["PsiDot"])-Alpha(x,u,dt,Wb)

def Velv(x,u,p,dt,Wb):
    return (x["Vvl"] + u["VvlDot"]*dt) + Wb["PFactor"]*cs.sin(p["phase"]+Chi(x,u,dt,Wb))*(x["Vvp"] + u["VvpDot"]*dt)

def Omega(x,u,dt,Wb):
    return u["PsiDot"]-AlphaDot(x,u,dt,Wb)

def FDotBlimp(x,wind,u,p,dt,Wb):
    f={}
    v=Velh(x,u,dt)
    chi=Chi(x,u,dt,Wb)
    f["Px"] = v*cs.cos(chi) + wind["Vx"]
    f["Py"] = v*cs.sin(chi) + wind["Vy"]
    f["Pz"] = Velv(x,u,p,dt,Wb) + wind["Vz"]
    f["Psi"] = u["PsiDot"]
    f["PsiDot"] = u["PsiDot"]-x["PsiDot"]
    f["Vh"] = u["VhDot"]
    f["Vvl"] = u["VvlDot"]
    f["Vvp"] = u["VvpDot"]
    return State(f)

def FDot(x,u,p,dt):
    f={}
    f["wind"]=FDotConst(x["wind"],dt)
    obstacles={}
    for o in x["obstacle"]._members:
        obstacles[o]=FDotConst(x["obstacle"][o],dt)
    f["obstacle"]=State(obstacles)
    blimps={}
    for b in x["blimp"]._members:
        blimps[b]=FDotBlimp(x["blimp"][b],x["wind"],u["blimp"][b],p["blimp"][b],dt,W_blimp[b])
    f["blimp"]=State(blimps)
    f["target"]=FDotTarget(x["target"],dt)
    return State(f)

def bound(v,minimum,maximum):
    return cs.fmin(cs.fmax(v,minimum),maximum)

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

def rotateIntoCameraFrame(xv,u,p,xt,Wb):
    global W_blimp
    g=9.81
    phi = cs.arctan2(u["PsiDot"]*Velh(xv,u,0),g)*Wb["PhiFactor"]
    theta = bound(-cs.arctan2(Velv(xv,u,p,dt,Wb),xv["Vh"])*Wb["ThetaFactor"],Wb["ThetaMin"],Wb["ThetaMax"])
    psi = xv["Psi"]
    rho = Wb["CamH"]
    beta = Wb["CamV"]

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

def cameraVector(xv,u,p,Wb):
    g=9.81
    phi = cs.arctan2(u["PsiDot"]*Velh(xv,u,0),g)*Wb["PhiFactor"]
    theta = bound(-cs.arctan2(Velv(xv,u,p,dt,Wb),xv["Vh"])*Wb["ThetaFactor"],Wb["ThetaMin"],Wb["ThetaMax"])
    psi = xv["Psi"]
    rho = Wb["CamH"]
    beta = Wb["CamV"]

    v0=State({"x":1.,"y":0.,"z":0.})

    pv={}
    pv["x"]=xv["Px"]
    pv["y"]=xv["Py"]
    pv["z"]=xv["Pz"]
    pv=State(pv)

    vectorBody = rotateZ(rotateY(v0,-beta),rho)
    return rotateZ(rotateY(rotateX(vectorBody,phi),theta),psi)

# calculates the nearest point (projection) of the camera vector to the target - good to plot errors
def closestPoint(xv,u,p,xt,Wb):
    cv = cameraVector(xv,u,p,Wb)
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

def costHVDE(x,u,p,t,Wb):
        #global visualweights_blimp
        global W

        tv = rotateIntoCameraFrame(x,u,p,t,Wb)
        # this forms an elypsoid with it's axis defined by the main weights - around [camdist,0,0] 
        return ((tv["x"]-Wb["CamDist"])*W["WeightDepthFactor"])**2 + (tv["y"])**2 + (tv["z"])**2

def costSoftConstraintMin(value,softMinimum,hardMinimum,weight):
    return (cs.fmax( 0, softMinimum-value )*weight/cs.fabs(hardMinimum-softMinimum))**2

u = cs.SX.sym('u', nu*N+nfp)
z0 = cs.SX.sym('z0', nx)

X=initState(z0)

cost = 0
stateconstraint = 0
P = initParameters(u[nu*N:])
for t in range(0, N):
    ui0=t*nu

    U = initUpdate(u[ui0:])

    X+=RungeKutta(FDot,X,U,P,dt)

    #cost
    for i in range(blimps):
        # constraints
        #   speed
        for key in shardlimits_blimp:
            if key!="Vv":
                l=X["blimp"][i][key]
            else:
                l=Velv(X["blimp"][i],U["blimp"][i],P["blimp"][i],0,W_blimp[i])
            stateconstraint += cs.fmax(0,  -(l-W_blimp[i][shardlimits_blimp[key][0]]))
            stateconstraint += cs.fmax(0,  (l-W_blimp[i][shardlimits_blimp[key][1]]))

        costSOFTCONSTRAINT = 0
        for key in ssoftlimits_blimp:
            if key!="Vv":
                l=X["blimp"][i][key]
            else:
                l=Velv(X["blimp"][i],U["blimp"][i],P["blimp"][i],0,W_blimp[i])
            costSOFTCONSTRAINT += costSoftConstraintMin(l, W_blimp[i][ssoftlimits_blimp[key][0]], W_blimp[i][shardlimits_blimp[key][0]], 1.0)
            costSOFTCONSTRAINT += costSoftConstraintMin(-l, -W_blimp[i][ssoftlimits_blimp[key][1]], -W_blimp[i][shardlimits_blimp[key][1]], 1.0)

        # limit change in PsiDot (hack!!!)
        stateconstraint += cs.fmax(0, cs.fabs(X["blimp"][i]["PsiDot"]-U["blimp"][i]["PsiDot"])-(W_blimp[i]["PsiDotDotMaxHard"]*dt))
        costSOFTCONSTRAINT += costSoftConstraintMin( -cs.fabs(X["blimp"][i]["PsiDot"]-U["blimp"][i]["PsiDot"]), -W_blimp[i]["PsiDotDotMaxSoft"]*dt , -W_blimp[i]["PsiDotDotMaxHard"]*dt,1.0)

        #   distance to target
        tvX = X["target"]["Px"] - X["blimp"][i]["Px"]
        tvY = X["target"]["Py"] - X["blimp"][i]["Py"]
        tvXY= cs.sqrt(tvX**2+tvY**2)
        #stateconstraint += cs.fmax(0, (W["VehicleHardExclusion"]+W["TargetExclusion"]) - tvXY)**2
        costSOFTCONSTRAINT += costSoftConstraintMin( tvXY, W["VehicleSoftExclusion"]+W["TargetExclusion"], W["VehicleHardExclusion"]+W["TargetExclusion"], W["AntiCollisionWeightFactor"])
        
        costSOFTCONSTRAINT += costSoftConstraintMin( -tvXY, -W["MaxDistanceSoft"], -W["MaxDistanceHard"], W["AntiCollisionWeightFactor"])

        # don't collide with obstacles
        for j in range(obstacles):
            vvX = X["blimp"][i]["Px"] - X["obstacle"][j]["Px"]
            vvY = X["blimp"][i]["Py"] - X["obstacle"][j]["Py"]
            # constraint - minimum distance between vehicles
            vvXY=cs.sqrt(vvX**2 +  vvY**2)
            stateconstraint += cs.fmax(0, (W["VehicleHardExclusion"]+X["obstacle"][j]["Size"]) - vvXY)**2
            costSOFTCONSTRAINT += costSoftConstraintMin( vvXY, W["VehicleSoftExclusion"]+X["obstacle"][j]["Size"], W["VehicleHardExclusion"]+X["obstacle"][j]["Size"], W["AntiCollisionWeightFactor"])

        costVDIST = 0
        for j in range(blimps):
            if i!=j:
                vvX = X["blimp"][i]["Px"] - X["blimp"][j]["Px"]
                vvY = X["blimp"][i]["Py"] - X["blimp"][j]["Py"]
                vvZ = X["blimp"][i]["Pz"] - X["blimp"][j]["Pz"]
                # constraint - minimum distance between vehicles
                vvXY=cs.sqrt(vvX**2 +  vvY**2)
                stateconstraint += cs.fmax(0, (2.0*W["VehicleHardExclusion"]) - vvXY)**2
                costSOFTCONSTRAINT += costSoftConstraintMin( vvXY, 2.0*W["VehicleSoftExclusion"], 2.0*W["VehicleHardExclusion"], W["AntiCollisionWeightFactor"])

                # cost for relative angle to target
                tv2X= X["target"]["Px"] - X["blimp"][j]["Px"]
                tv2Y= X["target"]["Py"] - X["blimp"][j]["Py"]

                angle = cs.acos((tvX*tv2X+tvY*tv2Y)/(tvXY*cs.sqrt(tv2X**2+tv2Y**2)))
                if blimps<=3:
                    costVDIST +=  (refangle-angle)**2
                else:
                    costVDIST += cs.fmax(0,(refangle-angle))**2

        costU = 0
        for key in uweights_blimp:
            costU += (uweights_blimp[key] * U["blimp"][i][key])**2
        costX = 0
        for key in xweights_blimp:
            costX += (xweights_blimp[key] * X["blimp"][i][key])**2

        cost += W["costHVDE"]*costHVDE(X["blimp"][i],U["blimp"][i],P["blimp"][i],X["target"],W_blimp[i])
        cost += W["costVDIST"]*costVDIST
        cost += W["costUWEIGHTS"]*costU
        cost += W["costXWEIGHTS"]*costX
        cost += W["costSOFTCONSTRAINT"]*costSOFTCONSTRAINT

costP = 0
for i in range(blimps):
    for key in pweights_blimp:
        costP += (pweights_blimp[key] * P["blimp"][i][key])**2
cost += W["costPWEIGHTS"]*costP

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
print ("defining problem:")
problem = og.builder.Problem(u, z0, cost).with_constraints(bounds)\
    .with_penalty_constraints(stateconstraint)
ros_config = og.config.RosConfiguration()\
    .with_package_name("nmpc_blimp_formation_planner")\
    .with_node_name("nmpc_blimp_planner")\
    .with_rate(W["UpdateRate"])\
    .with_queue_sizes(3)\
    .with_description("nmpc optimizing for blimp formation")
build_config = og.config.BuildConfiguration()\
    .with_build_directory("codegen")\
    .with_build_mode("release")\
    .with_build_c_bindings()  \
    .with_ros(ros_config)
meta = og.config.OptimizerMeta()\
    .with_optimizer_name("nmpc_blimp_formation_optimizer")\
    .with_version('0.1.1')\
    .with_licence('GPLv3')
solver_config = og.config.SolverConfiguration()\
        .with_tolerance(1e-3)\
        .with_initial_tolerance(1e-2)\
        .with_max_outer_iterations(20)\
        .with_delta_tolerance(1e-2)\
        .with_penalty_weight_update_factor(10.0)\
        .with_initial_penalty(100.0)\
        .with_max_duration_micros(1.0/W["UpdateRate"]*1000000.0)
builder = og.builder.OpEnOptimizerBuilder(problem, 
                                          meta,
                                          build_config, 
                                          solver_config) \
    .with_verbosity_level(1)
print ("building the builder:")
builder.build()




