#TestTrackR
#Port of "Track Test" from Fortran90 to R
#"Test Track" originally written by D. Nikezic and K.N. Yu
#See original program here: http://cpc.cs.qub.ac.uk/summaries/ADWT
#Reverse engineered by T.J. Wiegman
#Dept of Physics/Engineering at Point Loma Nazarene University, CA, USA

#Input data ----
alpha <- 1      #Particle energy in MeV
etchTime <- 8.5 #Etching time in hours
Vb <- 1.1       #Bulk etching rate in um/h
thetas <- 90    #Incident angle in degrees
equation <- 1   #Which model to use (1, 2, or 3)

#Calculated inputs ----
removed <- Vb*etchTime
if (thetas > 88) {thetas <- 88}
theta <- thetas*pi/180
if (equation == 1) {
  VtVb <- function(y) {
    #Using default constants: a1=11.45, a2=4, b1=0.339, b2=0.044, b3=0.58
    return(1 + (11.45*exp(-0.339*y) + 4*exp(-0.044*y))*(1 - exp(-0.58*y)))
  }
}

#"Detnik8a" Performing the actual track calculations
#Import SRIM2003 data ----
#energy in MeV, range in um (in CR39)
srim <- read.csv(text = {"Energy, Range
0.010, 0.1601
0.011, 0.1745
0.012, 0.1886
0.013, 0.2023
0.014, 0.2157
0.015, 0.2287
0.016, 0.2415
0.017, 0.2540
0.018, 0.2662
0.020, 0.2899
0.022, 0.3183
0.025, 0.3455
0.027, 0.3716
0.030, 0.3968
0.032, 0.4211
0.035, 0.4447
0.037, 0.4676
0.040, 0.4898
0.045, 0.5326
0.050, 0.5733
0.055, 0.6124
0.060, 0.6499
0.065, 0.6861
0.070, 0.7211
0.080, 0.7879
0.090, 0.8511
0.100, 0.9114
0.110, 0.9690
0.120, 1.02
0.130, 1.08
0.140, 1.13
0.150, 1.18
0.160, 1.23
0.170, 1.28
0.180, 1.32
0.200, 1.41
0.225, 1.52
0.250, 1.62
0.275, 1.72
0.300, 1.81
0.325, 1.91
0.350, 2.00
0.375, 2.09
0.400, 2.18
0.450, 2.35
0.500, 2.52
0.550, 2.69
0.600, 2.85
0.650, 3.02
0.700, 3.19
0.800, 3.53
0.900, 3.88
1.00, 4.24
1.10, 4.60
1.20, 4.98
1.30, 5.38
1.40, 5.78
1.50, 6.20
1.60, 6.63
1.70, 7.07
1.80, 7.52
2.00, 8.47
2.25, 9.73
2.50, 11.07
2.75, 12.49
3.00, 13.99
3.25, 15.57
3.50, 17.22
3.75, 18.95
4.00, 20.75
4.50, 24.52
5.00, 28.52
5.50, 32.85
6.00, 37.48
6.50, 42.42
7.00, 47.67
8.00, 59.02
9.00, 71.52
10.00, 85.15"}, header = TRUE)

#Check bounds of things ----
#Finding maximum of Vt/Vb
ytest <- seq(from = 0, to = 10, by = 0.001)
vmax <- max(VtVb(ytest))

#Check if theta is less than the critical angle
thetaLimit <- asin(1/vmax)
if (theta < thetaLimit || alpha < 0.1) {
  simpleError(message = "No track formation")
}

#90 Degree angles won't work, so approximate with ~88 degrees
if (theta > 0.99*pi/2) theta <- theta*0.98 #force an angle less than 90 degrees

#Get index of minimum srim energy that is greater than alpha
minE <- which.max(1/(srim$Energy - alpha))
if (srim$Energy[minE] < alpha) minE <- minE + 1

#Calculate range ----
#Linear interpolation of srim data
ran <- srim$Range[minE-1] + (((alpha - srim$Energy[minE-1])
                              *(srim$Range[minE] - srim$Range[minE-1]))
                              /(srim$Energy[minE] - srim$Energy[minE-1]))

#"dint3" Numerical integration of 1/VtVb(R-x)
dInt3 <- function(initial, final) {
  # #Calculate values
  # stepSize <- 0.01
  # x <- seq(from = initial, to = final, by = stepSize)
  # y <- 1/(VtVb(ran - x))
  # 
  # #Reimann Sum
  # p <- x[-1]
  # for (i in 2:length(x)) p[i-1] <- stepSize*mean(y[i-1], y[i])
  # return(sum(p))
  
  #R's built-in integrator is much, much faster than above
  return(integrate(f = function(x) {1/VtVb(x)},
                   lower = initial,
                   upper = final)$value)
}

#Vertical distance from detector surface to end of track
remLayer <- dInt3(initial = 0, final =  ran)

#If origin below surface ----
if(VtVb(ran)*sin(theta) < 1) {
  warning("Track does not start at detector surface")
  
  #Find the origin depth of track
  origin <- 0
  while(origin < ran) {
    if (VtVb(ran - origin)*sin(theta) >= 1) break
    origin <- origin + 0.01
    if (origin >= ran) simpleError(message = "No track formation")
  }
  hc <- origin*sin(theta)
  if (hc >= removed) simpleError(message = "No track formation")
  
  #Recalculating range based on new "surface"
  ran <- removed - hc
  etchTime <- etchTime - (hc/Vb)
  remLayer <- dInt3(initial = 0, final = ran)
}

#Etch time required to reach end of track
timeC <- remLayer/Vb

#Calculate depth ----
depth <- 0.01
while (depth < ran) {
  if (depth >= ran) {
    depth <- ran
    break
  }
  
  Db <- dInt3(initial = 0, final = depth)
  
  if ((Db/Vb) >= etchTime) {
    depth <- depth - 0.01
    break
  }
  
  depth <- depth + 0.01
}

print(depth)
