from telnetlib import X3PAD
from vpython import *
import numpy as np

# constants
g = 9.8                     # Gravitational constant
k = 10                      # spring constant
m = 0.1                     # mass of the bobs

tht1 = 20* np.pi/180        # angle of first mass
tht2 = 20* np.pi/180         # angle of second mass
tht3 = 20* np.pi/180
tht4 = 20* np.pi/180
tht5 = 20* np.pi/180

l = 2                    # lenght of pendula
d = 2                    # separation between masses

T  = 20                  # simulation time
dt = 0.001               # accuracy

# Coordinates of bobs
x1, y1 =     l* np.sin(tht1), -l* np.cos(tht1)
x2, y2 = d + l* np.sin(tht2), -l* np.cos(tht2)
x3, y3 = 2*d + l* np.sin(tht3), -l* np.cos(tht3)
x4, y4 = 3*d + l* np.sin(tht4), -l* np.cos(tht4)
x5, y5 = 4*d + l* np.sin(tht5), -l* np.cos(tht5)

# Canvas settings
scene = canvas(title = 'Coupled Oscillator',
               width = 600, height = 300,
               center= vector(0.5,-1,0), background = color.white)

Hook = box( pos = vector(2*d, 0, 0), size = vector(4*d + 0.5, 0.5, 0.5),
             color = color.yellow)

ball1 = sphere( pos = vector(x1, y1, 0), radius = 0.2, color = color.blue)
ball2 = sphere( pos = vector(x2, y2, 0), radius = 0.2, color = color.blue)
ball3 = sphere( pos = vector(x3, y3, 0), radius = 0.2, color = color.blue)
ball4 = sphere( pos = vector(x4, y4, 0), radius = 0.2, color = color.blue)
ball5 = sphere( pos = vector(x5, y5, 0), radius = 0.2, color = color.blue)

rod1 = cylinder(pos = vector(0, 0, 0), axis = vector(l* np.sin(tht1), -l* np.cos(tht1), 0), radius = 0.05)
rod2 = cylinder(pos = vector(d, 0, 0), axis = vector(l* np.sin(tht2), -l* np.cos(tht2), 0), radius = 0.05)
rod3 = cylinder(pos = vector(2*d, 0, 0), axis = vector(l* np.sin(tht3), -l* np.cos(tht3), 0), radius = 0.05)
rod4 = cylinder(pos = vector(3*d, 0, 0), axis = vector(l* np.sin(tht4), -l* np.cos(tht4), 0), radius = 0.05)
rod5 = cylinder(pos = vector(4*d, 0, 0), axis = vector(l* np.sin(tht5), -l* np.cos(tht5), 0), radius = 0.05)

Spring1 = cylinder(pos  = vector(l* np.sin(tht1)/2,-l* np.cos(tht1)/2, 0),
                  axis = vector(d + l* np.sin(tht2)/2 - l* np.sin(tht1)/2,-l* np.cos(tht2)/2 + l* np.cos(tht1)/2, 0),
                  radius = 0.05,
                  color  = color.orange)

Spring2 = cylinder(pos  = vector(l* np.sin(tht1)/2,-l* np.cos(tht1)/2, 0),
                  axis = vector(-d + l* np.sin(tht3)/2 - l* np.sin(tht2)/2,-l* np.cos(tht3)/2 + l* np.cos(tht2)/2, 0),
                  radius = 0.05,
                  color  = color.orange)

Spring3 = cylinder(pos  = vector(l* np.sin(tht1)/2,-l* np.cos(tht1)/2, 0),
                  axis = vector(2*d + l* np.sin(tht4)/2 - l* np.sin(tht3)/2,-l* np.cos(tht4)/2 + l* np.cos(tht3)/2, 0),
                  radius = 0.05,
                  color  = color.orange)

Spring4 = cylinder(pos  = vector(l* np.sin(tht1)/2,-l* np.cos(tht1)/2, 0),
                  axis = vector(3*d + l* np.sin(tht5)/2 - l* np.sin(tht4)/2,-l* np.cos(tht5)/2 + l* np.cos(tht4)/2, 0),
                  radius = 0.05,
                  color  = color.orange)

# Graph
gd = graph( width = 600, height = 300,
             title = '<b>Phase Space</b>',
             xtitle = '<i>Theta</i>', ytitle = '<i>Omega</i>',
             foreground = color.black, background = color.white)

phase_1 = gcurve( color = color.red,  label = 'm1' )
phase_2 = gcurve( color = color.blue, label = 'm2' )
phase_3 = gcurve( color = color.green, label = 'm3')
phase_4 = gcurve( color = color.yellow, label = 'm4')
phase_5 = gcurve( color = color.black, label = 'm5')

# Solver
omega_1 = 0
omega_2 = 0
omega_3 = 0
omega_4 = 0
omega_5 = 0

t = 0

while (t < T):

    rate(500)

    # angular velocity update
    omega_1 = omega_1 + ( -g* np.sin(tht1)/l + k* np.sin(2*(tht2-tht1))/ (8* m) )* dt
    omega_2 = omega_2 + ( -g* np.sin(tht2)/l - k* np.sin(2*(tht2-tht1))/ (8* m) + k* np.sin(2*(tht3-tht2))/ (8* m))* dt
    omega_3 = omega_3 + ( -g* np.sin(tht3)/l - k* np.sin(2*(tht3-tht2))/ (8* m) + k* np.sin(2*(tht4-tht3))/ (8* m))* dt
    omega_4 = omega_4 + ( -g* np.sin(tht4)/l - k* np.sin(2*(tht4-tht3))/ (8* m) + k* np.sin(2*(tht5-tht4))/ (8* m))* dt
    omega_5 = omega_5 + ( -g* np.sin(tht5)/l - k* np.sin(2*(tht5-tht4))/ (8* m) )* dt

    # angle update
    tht1 = tht1 + omega_1* dt
    tht2 = tht2 + omega_2* dt
    tht3 = tht3 + omega_3* dt
    tht4 = tht4 + omega_4* dt
    tht5 = tht5 + omega_5* dt    

    # time update
    t = t + dt

    # 1st mass
    ball1.pos.x =   l* np.sin(tht1)
    ball1.pos.y = - l* np.cos(tht1)

    rod1.axis.x =   l* np.sin(tht1)
    rod1.axis.y = - l* np.cos(tht1)

    # 2nd mass
    ball2.pos.x =  l* np.sin(tht2) + d
    ball2.pos.y =- l* np.cos(tht2)

    rod2.axis.x =   l* np.sin(tht2)
    rod2.axis.y = - l* np.cos(tht2)

    # 3rd mass
    ball3.pos.x =  l* np.sin(tht3) + 2*d
    ball3.pos.y =- l* np.cos(tht3)

    rod3.axis.x =   l* np.sin(tht3)
    rod3.axis.y = - l* np.cos(tht3)

    # 4th mass
    ball4.pos.x =  l* np.sin(tht4) + 3*d
    ball4.pos.y =- l* np.cos(tht4)

    rod4.axis.x =   l* np.sin(tht4)
    rod4.axis.y = - l* np.cos(tht4)

    # 5th mass
    ball5.pos.x =  l* np.sin(tht5) + 4*d
    ball5.pos.y =- l* np.cos(tht5)

    rod5.axis.x =   l* np.sin(tht5)
    rod5.axis.y = - l* np.cos(tht5)


    # Spring1
    Spring1.pos.x  =   l* np.sin(tht1)/2
    Spring1.pos.y  = - l* np.cos(tht1)/2
    Spring1.axis.x = d + l* np.sin(tht2)/2 - l* np.sin(tht1)/2
    Spring1.axis.y =   - l* np.cos(tht2)/2 + l* np.cos(tht1)/2

    # Spring2
    Spring2.pos.x  =   d + l* np.sin(tht2)/2
    Spring2.pos.y  = - l* np.cos(tht2)/2
    Spring2.axis.x = d + l* np.sin(tht3)/2 - l* np.sin(tht2)/2
    Spring2.axis.y =   - l* np.cos(tht3)/2 + l* np.cos(tht2)/2

    # Spring3
    Spring3.pos.x  =   2*d + l* np.sin(tht3)/2
    Spring3.pos.y  = - l* np.cos(tht4)/2
    Spring3.axis.x = d + l* np.sin(tht4)/2 - l* np.sin(tht3)/2
    Spring3.axis.y =   - l* np.cos(tht4)/2 + l* np.cos(tht3)/2

    # Spring4
    Spring4.pos.x  =   3*d + l* np.sin(tht4)/2
    Spring4.pos.y  = - l* np.cos(tht5)/2
    Spring4.axis.x = d + l* np.sin(tht5)/2 - l* np.sin(tht4)/2
    Spring4.axis.y =   - l* np.cos(tht5)/2 + l* np.cos(tht4)/2

    # # Phase space plot
    phase_1.plot( pos=(tht1, omega_1) )
    phase_2.plot( pos=(tht2, omega_2) )
    phase_3.plot( pos=(tht3, omega_3) )
    phase_4.plot( pos=(tht4, omega_4) )
    phase_5.plot( pos=(tht5, omega_5) )