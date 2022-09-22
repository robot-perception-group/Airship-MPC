import math

# W stores the majority of our data
W={ "TargetExclusion":1.0,                  # how large the subject is - don't fly over there
    "VehicleSoftExclusion":10.0,            # soft exclusion zone (soft constrain, cost starts raising if closer
    "VehicleHardExclusion":6.0,             # hard exclusion zone, state constraine, no vehicle may be closer than that
    "costHVDE":1.0,                         # relative cost calculated by distance between subject and point where the sensor points at in 3 dimensions
    "costVDIST":100.0,                      # relative cost of angular separation between vehicles error
    "costUWEIGHTS":0.01,                    # relative cost of control updates, near 0
    "costPWEIGHTS":0.01,                    # relative cost of free parameters, near 0
    "costXWEIGHTS":1.0,                     # relative cost of state penalties
    "costSOFTCONSTRAINT":100.0,             # relative cost of soft-constraints
    "CamH":82.*math.pi/180.,                # default sensor rho angle
    "CamV":30.*math.pi/180.,                # default sensor beta angle
    "CamDist":15.0,                         # default optimum viewing distance for sensor
    "ThetaMax":10.*math.pi/180.0,           # maximum theta the vehicle will achieve. steeper vertical trajectories will nto be reflected in increased theta
    "ThetaMin":-10.*math.pi/180.0,          # minimum theta
    "PsiDotDotMaxHard":3.0*math.pi/180,     # hard constraint for rotation change rate
    "PsiDotDotMaxSoft":1.0*math.pi/180,     # soft constraint for rotation change rate
    "PsiDotMax":18.0*math.pi/180,           # maximum rotation rate (hard constraint)
    "VhDotMin":-0.5,                        # maximum break rate of the vehicle (default)(m/s^2)
    "VhDotMax":1.0,                         # maximum acceleration rate of the vehicle (default)
    "VvlDotMin":-1.0,                       # linear vertical acceleration rate minimum (default)(m/s^2)
    "VvlDotMax":1.0,                        
    "VvpDotMin":-1.0,                       # periodic vertical acceleration rate minimum - see paper
    "VvpDotMax":1.0,
    "VhMinSoft":1.0,                        # soft constraint minimum velocity
    "VhMaxSoft":3.5,                        # soft constraint maximum velocity
    "VhMinHard":0.5,                        # hard constraint
    "VhMaxHard":4.0,
    "VvMinSoft":-1.0,                       # same for vertical velocity
    "VvMaxSoft":1.0,
    "VvMinHard":-1.5,
    "VvMaxHard":1.5,
    "PzMinSoft":-45.0,                      # maximum altitude allowed soft constraint
    "PzMaxSoft":-7.0,                       # minimum altitude allowed
    "PzMinHard":-50.0,                      # hard constraints
    "PzMaxHard":-2.0,
    "ThetaFactor":1.0,                      # theta as a function of vertical motion angle, 1.0 means vehicle points where it moves (limited by ThetaMax)
    "PhiFactor":1.0,                        # phi as a function of lateral acceleration versus gravity. 0.0 for ground vehicles
    "PFactor":1.0,                          # periodic factor in vertical speed. set to 0 to disable periodic component. linear component can not be disabled
    "WeightDepthFactor":0.6,                # 1.0 : distance is weighed as much as horizontal and vertical sensor displacement. 0: means distance is unconstrained, sensor can point to positive and negative infinity (including behind the sensor)
    "AntiCollisionWeightFactor":5.0,        # cost weights for anti-collision soft constraints get boosted by that much, multiplied with costSOFTCONSTRAINT above
    "MaxDistanceSoft": 40.0,                # soft constraint for max distance from subject - prevent flying too far away
    "MaxDistanceHard": 60.0,                # hard constraint (currently not really a hard constraint, only set to scale the soft constraint gradient)
    "UpdateRate": 4.0,                      # update rate of solver (frequenzy)
    "LowLevelOversteerVh": 4.0,             # oversteer low level controller by that much to get prompt response for sudden control changes
    "LowLevelOversteerVv": 4.0,
    "LowLevelOversteerPsiDot": 1.7,
    "VehicleCenterX": 2.31,                 # sensor location relative to vehicle position
    "VehicleCenterY": 0.0,                  # sensor location relative to vehicle position
    "VehicleCenterZ": 2.12,                 # sensor location relative to vehicle position
    #"Cl":10000000,
    "Cl":0.24,                               # combined lift+drag coefficient (default) see below
    }
#note: Tested with minVh 1.0, maxVh 4.0, omega 10deg/s - no theta
#      Cl of 0.24 measured in simulation
#      Cl of 0.7 is barely manageable with max AoA of 15deg
#      Cl of 0.8 max AoA of 10deg
#      Cl of 1.0 -> 8deg
#      Cl of 5.0 -> 2deg

state_blimp_const=["CamH","CamV","CamDist","Cl","PhiFactor","ThetaFactor","PFactor","PsiDotMax","LowLevelOversteerVh","LowLevelOversteerVv","LowLevelOversteerPsiDot", "VhDotMin","VhDotMax","VvlDotMin","VvlDotMax","VvpDotMin","VvpDotMax","PsiDotDotMaxSoft","VhMinSoft","VhMaxSoft","VvMinSoft","VvMaxSoft","PzMinSoft","PzMaxSoft","PsiDotDotMaxHard","VhMinHard","VhMaxHard","VvMinHard","VvMaxHard","PzMinHard","PzMaxHard","ThetaMax","ThetaMin","VehicleCenterX","VehicleCenterY","VehicleCenterZ"]

blimps=3
W_blimp=[]
for i in range(blimps):
    w={}
    for k in state_blimp_const:
        w[k]=W[k]
    W_blimp.append(w)
obstacles=0
olocations=[ [-7.,-14.,5.],[3.,13.,5.],[6.,-12.,5.],[6.0,40.0,5.],[-22.0,22.0,5.0] ]

N=10
dt=1.25

