#!/usr/bin/env python

import math
import rospy
from tf import transformations
from geometry_msgs.msg import *
from sensor_msgs.msg import *
from std_msgs.msg import *
from visualization_msgs.msg import *
from uav_msgs.msg import *
from nmpc_blimp_formation_planner.msg import *
import sys
import copy
import os
import numpy as np

sys.path.append(os.path.dirname(__file__))
import formation_config as FormationConfig

W=FormationConfig.W
W_blimp=FormationConfig.W_blimp
blimps=FormationConfig.blimps
obstacles=FormationConfig.obstacles
olocations=FormationConfig.olocations
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

def initState(z0,selfID):
    xi=0
    X={}
    X["wind"]={}
    for key in state_wind:
        X["wind"][key]=z0[xi]
        xi+=1
    xi+=obstacles*len(state_obstacle)
    X["blimp"]={}
    for i in range(blimps):
        if i==selfID-1:
            X["blimp"][i]={}
            for key in state_blimp:
                X["blimp"][i][key]=z0[xi]
                xi+=1
        else:
            xi+=len(state_blimp)
    xi+=len(state_target)
    return State(X)

def initUpdate(u0,selfID):
    ui=0
    U={}
    U["blimp"]={}
    for i in range(blimps):
        U["blimp"][i]={}
        for key in update_blimp:
            U["blimp"][i][key]=u0[ui]
            ui+=1
    return State(U)

def initParameters(p0,selfID):
    pi=0
    P={}
    P["blimp"]={}
    for i in range(blimps):
        if i==selfID-1:
            P["blimp"][i]={}
            for key in parameters_blimp:
                P["blimp"][i][key]=p0[pi]
                pi+=1
        else:
            pi+=len(parameters_blimp)
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
    return (x["Vvl"] + u["VvlDot"]*dt) + Wb["PFactor"]*math.sin(p["phase"]+Chi(x,u,dt,Wb))*(x["Vvp"] + u["VvpDot"]*dt)

def Omega(x,u,dt,Wb):
    return u["PsiDot"]-AlphaDot(x,u,dt,Wb)

def FDotBlimp(x,wind,u,p,dt,Wb):
    f={}
    v=Velh(x,u,dt)
    chi=Chi(x,u,dt,Wb)
    f["Px"] = v*math.cos(chi) + wind["Vx"]
    f["Py"] = v*math.sin(chi) + wind["Vy"]
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
    #obstacles={}
    #for o in x["obstacle"]._members:
    #    obstacles[o]=FDotConst(x["obstacle"][o],dt)
    #f["obstacle"]=State(obstacles)
    blimps={}
    for b in x["blimp"]._members:
        blimps[b]=FDotBlimp(x["blimp"][b],x["wind"],u["blimp"][b],p["blimp"][b],dt,W_blimp[b])
    f["blimp"]=State(blimps)
    #f["target"]=FDotTarget(x["target"],dt)
    return State(f)

def bound(v,minimum,maximum):
    return min(max(v,minimum),maximum)

def rotateX(v,a):
    f={}
    f["x"]=v["x"]
    f["y"]=math.cos(a)*v["y"]-math.sin(a)*v["z"]
    f["z"]=math.cos(a)*v["z"]+math.sin(a)*v["y"]
    return State(f)

def rotateY(v,a):
    f={}
    f["x"]=math.cos(a)*v["x"]+math.sin(a)*v["z"]
    f["y"]=v["y"]
    f["z"]=math.cos(a)*v["z"]-math.sin(a)*v["x"]
    return State(f)

def rotateZ(v,a):
    f={}
    f["x"]=math.cos(a)*v["x"]-math.sin(a)*v["y"]
    f["y"]=math.cos(a)*v["y"]+math.sin(a)*v["x"]
    f["z"]=v["z"]
    return State(f)

def rotateIntoCameraFrame(xv,u,p,xt,Wb):
    global W_blimp
    g=9.81
    phi = math.atan2(u["PsiDot"]*Velh(xv,u,0),g)*Wb["PhiFactor"]
    theta = bound(-math.atan2(Velv(xv,u,p,dt,Wb),xv["Vh"])*Wb["ThetaFactor"],Wb["ThetaMin"],Wb["ThetaMax"])
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
    phi = math.atan2(u["PsiDot"]*Velh(xv,u,0),g)*Wb["PhiFactor"]
    theta = bound(-math.atan2(Velv(xv,u,p,dt,Wb),xv["Vh"])*Wb["ThetaFactor"],Wb["ThetaMin"],Wb["ThetaMax"])
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
    return (max( 0, softMinimum-value )*weight/math.fabs(hardMinimum-softMinimum))**2


def dopredict(z0,u,selfID):
    trajectory=[]
    X=initState(z0,selfID)*1.0

    cost = 0
    stateconstraint = 0
    P = initParameters(u[nu*N:],selfID)
    for t in range(0, N):
        ui0=t*nu

        U = initUpdate(u[ui0:],selfID)

        X+=RungeKutta(FDot,X,U,P,dt)
        trajectory.append(X*1.0)
    return trajectory

COMMANDTOPIC="command"
POSETOPIC="throttledUAVPose"
TARGETPOSETOPIC="target_tracker/pose"
TARGETTWISTTOPIC="target_tracker/twist"
MPCRESULTTOPIC="nmpc_blimp_planner/result"
MPCSTATETOPIC="nmpc_blimp_planner/parameters"
VISTOPIC="commandvis"
OBSTOPIC="obstaclevis"
WINDTOPIC="winddebug"
NAMESPACE_PREFIX="/machine_"
MACHINES=FormationConfig.blimps #  MUST match solver definition
OBSTACLECOUNT=0 # MUST match solver definition
SELF=1
NAME='blimp_nmpc_wrapper_node'
TIMESTEP=FormationConfig.dt # 0.25 seconds per set - MUST mach solver definition
#LOOKAHEAD=int(2.0/TIMESTEP) #2 seconds into the future for command
MPC=None
PI=math.pi
RATE=W["UpdateRate"]
OLDSTATE=[]
commandpub=None
statepub=None
vispub=None
windpub=None

def param(s,default):
    if rospy.has_param(s):
        return rospy.get_param(s)
    return default


command_blimp=["PsiDot"]

class MPC(object):
    blimps=[]
    obstacles=[]
    target=[]
    wind=[]
    windlowpassalpha=0.1
    psidotlowpassalpha=0.2
    minAirspeedEpsilon=0.2
    alphacutoff=25*math.pi/180.0  # above this AoA do not perform wind calculation

    oldstate=None
    oldcommand=[]
    newstate=None

    def __init__(self,selfid,machines,numobstacles):
        global W,NAME,W_blimp
        self.selfid = selfid
        self.machines = machines
        self.numobstacles = numobstacles

        for a in range(machines):
            x=[10.0,10.0*a,W_blimp[a]["PzMaxSoft"],0.0,0.0,W_blimp[a]["VhMinSoft"],0.0,0.0,0.0,0.0,0.0]
            y=[0.0]*len(command_blimp)
            self.blimps.append(x)
            self.oldcommand.append(y)
        for a in range(self.numobstacles):
            self.obstacles.append(olocations[a])
        self.target=[0.0]*len(state_target)
        self.wind=[0.0]*len(state_wind)
        #self.windlowpassalpha=self.windlowpassalpha/self.machines

    def set_target_pose(self,x,y,z):
        self.target[0]=x
        self.target[1]=y
        self.target[2]=z

    def set_target_twist(self,x,y,z):
        self.target[3]=x
        self.target[4]=y
        self.target[5]=z

    def set_blimp_pose(self,rID,x,y,z,psi,psidot,vh,vv,phi,theta):
        global CAMH,CAMV,CAMD,NAME,CL,PHIFACTOR,THETAFACTOR,PFACTOR
        self.blimps[rID-1][0]=x
        self.blimps[rID-1][1]=y
        self.blimps[rID-1][2]=z
        self.blimps[rID-1][3]=psi
        self.blimps[rID-1][4]=bound(self.oldcommand[rID-1][0],-W_blimp[rID-1]["PsiDotMax"],W_blimp[rID-1]["PsiDotMax"])
        self.blimps[rID-1][5]=vh
        self.blimps[rID-1][6]=vv
        self.blimps[rID-1][7]=0.0
        self.blimps[rID-1][8]=phi
        self.blimps[rID-1][9]=theta
        self.blimps[rID-1][10]=psidot
        if (rID==self.selfid and self.oldstate is None):
            self.optimize()


    def set_wind(self,x,y,z):
        self.wind[0]=self.windlowpassalpha*x+(1.0-self.windlowpassalpha)*self.wind[0]
        self.wind[1]=self.windlowpassalpha*y+(1.0-self.windlowpassalpha)*self.wind[1]
        self.wind[2]=self.windlowpassalpha*z+(1.0-self.windlowpassalpha)*self.wind[2]

    def quat_mult(self,q,r):
        return [r[0]*q[0]-r[1]*q[1]-r[2]*q[2]-r[3]*q[3],
                r[0]*q[1]+r[1]*q[0]-r[2]*q[3]+r[3]*q[2],
                r[0]*q[2]+r[1]*q[3]+r[2]*q[0]-r[3]*q[1],
                r[0]*q[3]-r[1]*q[2]+r[2]*q[1]+r[3]*q[0]]

    def rotate_by_q(self,v,q):
        r = [0.]+v
        q_conj = [q[0],-1*q[1],-1*q[2],-1*q[3]]
        return self.quat_mult(self.quat_mult(q,r),q_conj)[1:]

    def process_blimppose(self,rID,pose):
        global W,W_blimp
        (phi,theta,psi)=transformations.euler_from_quaternion([pose.orientation.y,-pose.orientation.x,pose.orientation.z,pose.orientation.w])
        centerpoint =  self.rotate_by_q([W_blimp[rID-1]["VehicleCenterX"],W_blimp[rID-1]["VehicleCenterY"],W_blimp[rID-1]["VehicleCenterZ"]], [pose.orientation.w,pose.orientation.x,pose.orientation.y,pose.orientation.z])
        psidot = pose.angVelocity.z * math.pi/180.
        psidot = psidot*self.psidotlowpassalpha + self.blimps[rID-1][10]*(1.0-self.psidotlowpassalpha) # this is required, both for smoothing and since alpha does build up gradually, not instantly with changes in psidot
        motion_direction=math.atan2(pose.velocity.y,pose.velocity.x)
        real_alpha=psi-motion_direction
        real_airspeed=math.sqrt(pose.velocity.x**2 + pose.velocity.y**2 + pose.velocity.z**2 )
        real_Cl=psidot/(real_alpha*real_airspeed)
        airspeed=pose.POI.x
        airspeedH=airspeed # first order approximation, this is likely too low due to angle of attack
        if airspeed>=self.minAirspeedEpsilon and airspeed>abs(pose.velocity.z):
            for a in range(3):  # iteratively estimate horizontal and vertical angle of attack
                alpha=psidot/(W_blimp[rID-1]["Cl"]*airspeed)
                valpha = theta-math.asin(pose.velocity.z/airspeed)
                airspeed=pose.POI.x/(math.cos(alpha)*math.cos(valpha))
            if abs(alpha)<self.alphacutoff and abs(valpha)<self.alphacutoff:
                # estimate wind only if estimated angle of attack is sane.
                airspeedH=airspeed/math.cos(theta-valpha)
                airvecX=math.cos(psi-alpha)
                airvecY=math.sin(psi-alpha)
                windX = pose.velocity.x - airvecX*airspeedH
                windY = pose.velocity.y - airvecY*airspeedH
                if rID==self.selfid and airspeed>=self.minAirspeedEpsilon:
                    self.set_wind(windX,windY,0) # each vehicle does its own wind estimation
        vH = bound(airspeedH,W_blimp[rID-1]["VhMinHard"],W_blimp[rID-1]["VhMaxHard"]) # by forcing state within bounds, solver will converge even if vehicle out of bounds
        vV = bound(pose.velocity.z,W_blimp[rID-1]["VvMinHard"],W_blimp[rID-1]["VvMaxHard"])
        pZ = bound(pose.position.z+centerpoint[2],W_blimp[rID-1]["PzMinHard"],W_blimp[rID-1]["PzMaxHard"])
        self.set_blimp_pose(rID,pose.position.x+centerpoint[0],pose.position.y+centerpoint[1],pZ,psi,psidot,vH,vV,phi,theta)

    def optimize(self):
        global statepub,N
        self.newstate=[]
        self.newstate+=self.wind
        for o in self.obstacles:
            self.newstate+=o
        for b in self.blimps:
            self.newstate+=b[0:8]
        self.newstate+=self.target
        self.oldstate=self.newstate
        params = OptimizationParameters()
        params.parameter=self.oldstate
        params.initial_guess=[1e-6]*(nu*N+nfp)
        params.initial_penalty=100.0
        statepub.publish(params)

    def predict(self,updates,status):
        #global LOOKAHEAD,NAME,vispub,windpub,N
        global NAME,vispub,windpub,N
        if self.oldstate is None:
            return None

        x=initState(self.oldstate,self.selfid)*1.0
        p=initParameters(self.oldstate,self.selfid)
        u=initUpdate(updates,self.selfid)
        target={}
        ti=0
        for key in state_target:
            target[key]=self.target[ti]
            ti+=1
        target=State(target)
        state_trajectory=dopredict(self.oldstate,updates,self.selfid)
        for i in range(self.machines):
            if not math.isnan(u["blimp"][i]["PsiDot"]):
                self.oldcommand[i][0] = u["blimp"][i]["PsiDot"]
        result=(
                x["blimp"][self.selfid-1]["Vh"] + u["blimp"][self.selfid-1]["VhDot"] * W_blimp[self.selfid-1]["LowLevelOversteerVh"],
                (u["blimp"][self.selfid-1]["PsiDot"] - self.blimps[self.selfid-1][10]) * W_blimp[self.selfid-1]["LowLevelOversteerPsiDot"],
                (x["blimp"][self.selfid-1]["Vvl"] + u["blimp"][self.selfid-1]["VvlDot"] * W_blimp[self.selfid-1]["LowLevelOversteerVv"])
                    + math.sin( p["blimp"][self.selfid-1]["phase"] + W_blimp[self.selfid-1]["PFactor"] * Chi(x["blimp"][self.selfid-1],u["blimp"][self.selfid-1],0,W_blimp[self.selfid-1]))
                        * (x["blimp"][self.selfid-1]["Vvp"] + u["blimp"][self.selfid-1]["VvpDot"] * W_blimp[self.selfid-1]["LowLevelOversteerVv"]),
                self.target[0],
                self.target[1],
                self.target[2],
                closestPoint(x["blimp"][self.selfid-1],u["blimp"][self.selfid-1],p["blimp"][self.selfid-1],target,W_blimp[self.selfid-1])
                )
        trajectory=Marker()
        trajectory.header.stamp=rospy.Time.now()
        trajectory.header.frame_id="world"
        trajectory.type=Marker.LINE_STRIP
        trajectory.color.a=1.0
        trajectory.color.r=1.0
        trajectory.color.g=1.0
        trajectory.color.b=1.0
        trajectory.scale.x=0.1
        trajectory.scale.y=0
        trajectory.scale.z=0
        trajectory.points=[]
        trajectory.colors=[]
        def addpoint(x,y,z,i,speed):
            a=Point()
            a.x=x
            a.y=y
            a.z=z
            trajectory.points.append(a)
            a=ColorRGBA()
            a.a=1.0
            if (status==0):
                a.r=(1.0*i)/N
                a.g=1.0 - (1.0*i)/N
                a.b=speed/W_blimp[self.selfid-1]["VhMaxHard"]
            else:
                a.r=1.0
                a.g=0.0
                a.b=0.0
            trajectory.colors.append(a)
        addpoint(x["blimp"][self.selfid-1]["Px"],x["blimp"][self.selfid-1]["Py"],x["blimp"][self.selfid-1]["Pz"],0,x["blimp"][self.selfid-1]["Vh"])
        for i in range(N):
            if (i%(max(N/10,1))==0):
                addpoint(state_trajectory[i]["blimp"][self.selfid-1]["Px"],state_trajectory[i]["blimp"][self.selfid-1]["Py"],state_trajectory[i]["blimp"][self.selfid-1]["Pz"],0,state_trajectory[i]["blimp"][self.selfid-1]["Vh"])
        vispub.publish(trajectory)
        ma=MarkerArray()
        oid=0
        for o in self.obstacles:
            mo =  Marker()
            mo.ns="nmpsobstacles"
            mo.id=oid
            oid+=1
            mo.header.stamp=rospy.Time.now()
            mo.header.frame_id="world"
            mo.type=Marker.CYLINDER
            mo.color.a=1.0
            mo.color.r=1.0
            mo.color.g=1.0
            mo.color.b=1.0
            mo.scale.x=o[2]*2.0
            mo.scale.y=o[2]*2.0
            mo.scale.z=50.0
            mo.pose.position.x=o[0]
            mo.pose.position.y=o[1]
            mo.pose.position.z=-25
            mo.pose.orientation.w=1.0
            ma.markers.append(mo)
        obspub.publish(ma)

        wind=PointStamped()
        wind.header.frame_id="world"
        wind.header.stamp=rospy.Time.now()
        wind.point.x=self.wind[0]
        wind.point.y=self.wind[1]
        wind.point.z=self.wind[2]
        windpub.publish(wind)
        self.oldstate=None
        return result
     
def PoseCallback(msg, rID):
    global MPC
    MPC.process_blimppose(rID,msg)

def TargetPoseCallback(msg):
    global MPC
    MPC.set_target_pose(msg.pose.pose.position.x,msg.pose.pose.position.y,msg.pose.pose.position.z)

def TargetTwistCallback(msg):
    global MPC
    MPC.set_target_twist(msg.twist.twist.linear.x,msg.twist.twist.linear.y,msg.twist.twist.linear.z)

def OptimizationResultCallback(msg):
    global MPC,commandpub,W_blimp

    print("Received result: ",msg.status)
    newstate=MPC.predict(msg.solution,msg.status)
    if newstate is None or msg.status>=2:
        print("Ignoring")
        return

    command=uav_pose()
    command.header.stamp=rospy.Time.now()
    command.header.frame_id="world"
    command.flightmode=4 #ROSBRIDGEMESSAGE_FLIGHTCONTROL_MODE_ROSCONTROL = 4
    command.position.x=newstate[3]
    command.position.y=newstate[4]
    command.position.z=newstate[5]
    command.velocity.x=newstate[0]
    command.velocity.y=newstate[1]*180.0/math.pi
    command.velocity.z=newstate[2]
    command.POI.x=newstate[3]
    command.POI.y=newstate[4]
    command.POI.z=newstate[5]
    pv=newstate[6]
    command.covariance[0]=pv["x"]
    command.covariance[1]=pv["y"]
    command.covariance[2]=pv["z"]
    commandpub.publish(command)

if __name__ == '__main__':
    rospy.init_node(NAME)
    NAME=rospy.get_name()
    SELF=param(NAME+"/robotID",SELF)
    MACHINES=param(NAME+"/numRobots",MACHINES)
    NAMESPACE_PREFIX=param(NAME+"/PREFIX",NAMESPACE_PREFIX)
    # topic names set via params, since direct override isn't possible due to multi subscription
    POSETOPIC=param(NAME+"/POSETOPIC",POSETOPIC)
    MPC=MPC(SELF,MACHINES,obstacles)

    subscribers=[]
    for robot in range(1,MACHINES+1):
        subscribers.append( rospy.Subscriber(NAMESPACE_PREFIX+str(robot)+"/"+POSETOPIC,uav_pose,PoseCallback,robot,queue_size=3))
    subscribers.append( rospy.Subscriber(TARGETPOSETOPIC,PoseWithCovarianceStamped,TargetPoseCallback,queue_size=3))
    subscribers.append( rospy.Subscriber(TARGETTWISTTOPIC,TwistWithCovarianceStamped,TargetTwistCallback,queue_size=3))
    subscribers.append( rospy.Subscriber(MPCRESULTTOPIC,OptimizationResult,OptimizationResultCallback,queue_size=3))

    commandpub = rospy.Publisher(COMMANDTOPIC,uav_pose,queue_size=3)
    statepub = rospy.Publisher(MPCSTATETOPIC,OptimizationParameters,queue_size=3)
    vispub = rospy.Publisher(VISTOPIC,Marker,queue_size=3)
    obspub = rospy.Publisher(OBSTOPIC,MarkerArray,queue_size=3)
    windpub = rospy.Publisher(WINDTOPIC,PointStamped,queue_size=3)
    rospy.spin()
    
