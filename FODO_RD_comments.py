#import all the things from warp
from warp import *

#To make or not to make movies of 3D plots right now
l_movieplot = False
l_movieplot3d = False

#Set descriptions for the plots
top.pline2   = "Example 3D beam in a FODO lattice"
top.pline1   = "Semi-Gaussian cigar beam. 32x32x128"
top.runmaker = "David P. Grote"

#Do the plotting setup to make a cgm file at the end
setup()

#Create the beam species. This is an instance of the Species class
#This is a beam of Potassium ions
beam = Species(type=Potassium,charge_state=+1,name="Beam species")

#Set up some of the beam parameters
#All parameters are in mks units
beam.a0    = 8.760439903086566*mm #initial beam width in x direction
beam.b0    = 15.599886448447793*mm #initial beam width in y direction
beam.emit  = 6.247186343204832e-05 #transverse beam emittance
beam.ap0   = 0. #initial beam envelope vx/vz
beam.bp0   = 0. #initial beam envelope vy/vz
beam.ibeam = 2.*mA #current in the beam
beam.vbeam = 0. #beam speed
beam.ekin  = 80.*kV #beam kinetic energy

#This command, which is short for derivequantity, will derive some quantities
#based on the beam parameters given above
#For example, it will calculate beam.vbeam from beam.ekin
derivqty()

#Set the longitudinal thermal velocity spread based on the emittance
beam.vthz = 0.5*beam.vbeam*beam.emit/sqrt(beam.a0*beam.b0)

#Set up the FODO lattice!
hlp     = 36.*cm   #half of a  lattice period length
piperad = 3.445*cm #radius of the beam pipe
quadlen = 11.*cm   #length of the main quadrupoles

#Set the quadrupole magnetic field gradient
#This will determine the phase advance of the lattice
dbdx = 0.93230106124518164/quadlen

#Set up a focusing and defocusing quadrupole for one lattice period
addnewquad(zs= 0.5*hlp - quadlen/2.,
           ze= 0.5*hlp + quadlen/2.,
           db=+dbdx)
addnewquad(zs=1.5*hlp - quadlen/2.,
           ze=1.5*hlp + quadlen/2.,
           db=-dbdx)

#This parameter defines the start of the periodicity of the lattice
#it is relative to the quadrupole position
top.zlatstrt  = 0.

#This is the length of one lattice period
top.zlatperi  = 2.0*hlp

# ------------------------------------------------------------------------
#This section sets up and runs the envelope solver

#length of the tune per lattice
top.tunelen = 2.*hlp

#Begginning and end of the envelope calculation
#this calculation covers the length of the lattice
env.zl = -2.5*hlp #z direction lower limit
env.zu = -env.zl  #z direction upper limit
env.dzenv = top.tunelen/100. #calculate the step size for the envelope solver

#Call the envelope solver
package("env")
generate()
step()

#Make a plot of the solution
penv()
fma()

# ------------------------------------------------------------------------
#Set up some parameters that describe the 3D simulation


steps_p_perd = 50 #size of the tiem step

#specify that the time step should take the above number of time steps
#for each lattice period
top.dt = (top.tunelen/steps_p_perd)/beam.vbeam

#number of grid cells for each spatial dimension
w3d.nx = 32
w3d.ny = 32
w3d.nz = 128

#Specify the upper and lower limits of the grid in all three dimensions
#in this case the x and y upper and lower limits are set to be the
#the length of the beam pipe radius
w3d.xmmin = -piperad
w3d.xmmax =  piperad
w3d.ymmin = -piperad
w3d.ymmax =  piperad
w3d.zmmin = -hlp*2
w3d.zmmax = +hlp*2

#Set the boundary conditions for the edges of the grid
#this can be set to dirichlet, periodic, or neumann
w3d.bound0 = dirichlet
w3d.boundnz = dirichlet
w3d.boundxy = dirichlet

#Boundary conditions for the particles on the outer edges of the grid
#Possible values are absorb, periodic, or reflect
top.pbound0 = absorb
top.pboundnz = absorb
top.pboundxy = absorb

#Specify the length of the beam pulse
#in this case, it is 80% of the length of the grid
beam.zimin = w3d.zmmin*0.8
beam.zimax = w3d.zmmax*0.8

#Numer of particles to simulate
top.npmax = 200000

#Specify the distribution of the Beam
#the options are semigauss, KV, or WB
w3d.distrbtn = "semigaus"

#Specigy the distribution of the beam longitudinal velocity
w3d.distr_l = "gaussian"

#This turns on the cigar option, if selected this make the beams sort of taper
#off at the ends of the beam, which will make a cigar shape, hence the name
#it will also adjust the beam envelope to stay mathced
#See the cigar plots made by this script
w3d.cigarld = true

#Specify the fraction of the beam that is in the middle and not part of the
#cigar tapering
beam.straight = 0.5

#Set up the filed solver and specity the species of the FFT solver
top.fstype = 0

#This is an option for imposing symmetry. If this parameter were true, the
#fields would only be calculated in the transverse quadrant and mirrored in
#the rest of the quadrants, but in this case we don't want that
w3d.l4symtry = false

# --- Setup various diagnostics and plots.
# --- By default, Warp calculates all 1st and 2nd order moments of the particles
# --- as a function of z position.


#zwindows specifies the locations where histories of the momnets are saved by
#warp. In this case, the momnents at the center of the window are saved.
top.zwindows[:,1] = [-0.35, -0.3]
top.zwindows[:,2] = [-0.25, 0.25]
top.zwindows[:,3] = [0.3, 0.35]

#Turn on the saving of the time histories for some values that are not
#saved by default
top.lhxrmsz = true #save X rms values
top.lhyrmsz = true #save Y rms values
top.lhepsnxz = true #save normalized X emittance values
top.lhepsnyz = true #save normalized Y emittance values
top.lhcurrz = true #save beam current values

#The history of the values specified above are saved every nhist time steps
top.nhist = 1

#zzplalways = [zstart, zend, zperiod, extra_z_values, ...]
#specify some plots to generate and the frequency with which they are generated
#This results in the same kind of plot being made multiple times,
#at different time steps/positions in the lattice
top.zzplalways[0:4] = [0.,100000.,2*hlp,0.]

#ipzxy are space plot of the spatial locations of the particles in the x axis vs z axis
#and the y axis vs z axis plotted on the same page

#ipzvz is a plot of the particles velocity in the z direction vs their
#position in the z direction

top.ipzxy[-2] = always #call ppzxy as specified by zzplalways
top.ipzvz[-2] = always #call ipzvs as specified by zzplalways


# @callfromafterstep is a python decorator that specifies that this function
#gets called after every time step
@callfromafterstep
def runtimeplots(nsteps=steps_p_perd):
    "Make user defined plots, every steps_p_perd steps"
    if top.it%nsteps != 0:
        return
    # --- Create overlaid plots in subframes of the plot window.
    plsys(9)
    pfzx(cellarray=1, contours=0, centering='cell')
    pzxedges(color=red, titles=False)
    plsys(10)
    pfzy(cellarray=1, contours=0, centering='cell')
    pzyedges(color=red, titles=False)
    fma()

    # --- Make plots of the transverse distribution in two zwindows.
    plsys(1)
    ppxy(centering='cell')
    limits(-0.02, +0.02, -0.02, +0.02)
    fma()
    plsys(1)
    ppxxp(centering='cell')
    limits(-0.02, +0.02, -0.04, +0.04)
    fma()
    plsys(1)
    ppxy(centering='cell')
    limits(-0.02, +0.02, -0.02, +0.02)
    fma()
    plsys(1)
    ppxxp(centering='cell')
    limits(-0.02, +0.02, -0.04, +0.04)
    fma()

# --- Switch to the w3d package, which runs the 3D PIC model.
# --- The generate command does the initialization, including creating
# --- the particles, doing the initial Poisson solve, and calculating
# --- initial diagnostics and moments.
package("w3d")
generate()

# --- Directly call the user defined function, producing plots of the initial conditions.
runtimeplots()

# --- Run for 50 time steps.
# --- Note that after each time step, the routine runtimeplots will be automatically called.
step(50)

# --- Make various post processing diagnostic plots.
ptitles('Beam X envelope history, in beam frame', 'Lattice periods', 'Beam frame location (m)', 'Envelope is 2*Xrms')
ppgeneric(gridt=2.*top.hxrmsz[:,:top.jhist,0], xmin=0., xmax=top.zbeam/(2.*hlp), ymin=w3d.zmmin, ymax=w3d.zmmax)
fma()

ptitles('Beam X normalized emittance history, in beam frame', 'Lattice periods', 'Beam frame location (m)')
ppgeneric(gridt=top.hepsnxz[:,:top.jhist,0], xmin=0., xmax=top.zbeam/(2.*hlp), ymin=w3d.zmmin, ymax=w3d.zmmax)
fma()

hpepsnx()
hpepsny(titles=0)
fma()
