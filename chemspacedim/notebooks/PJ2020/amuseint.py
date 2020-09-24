"""
Class to facilitate integration of chemically tagged cluster, with some 
helper functions. Modified from AMUSE script found on
https://galpy.readthedocs.io/en/latest/potential.html.

Price-Jones, UofT, 2019 

"""
import copy
from tqdm import tqdm
import pandas as pd
from amuse.lab import nbody_system,new_plummer_sphere,BHTree,units
from amuse.couple import bridge
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.datamodel import Particles

from astropy import units as u
from astropy.coordinates import SkyCoord
from galpy.orbit import Orbit
from galpy import potential
from galpy.util import bovy_conversion,bovy_coords
import numpy as np

def get_potential(bar=False,dw_arm=False,trans_arm=False,trans_dwarm=False,
                  Nwaves=0, tstart=None,ro=8.,vo=220.,nomwpot=False,
                  transwrap=True,steady_arm=True):
    #Example, for SS paper the transient arm potential is called via:
    #pot=get_potential(bar=True,dw_arm=False,trans_arm=True,trans_dwarm=False,Nwaves=0, tstart=-5000.,ro=8.,vo=220.,nomwpot=False,transwrap=True)
    #tstart assumes Myr!
    """
    Create new galpy potential.

    bar:            Boolean control to include bar potential
    dw_arm:         Boolean control to include density wave spiral
    trans_arm:      Boolean control to include transient spiral
    trans_dwarm:    Boolean control to include transient density wave spiral
    Nwaves:         Number of spiral waves
    tstart:         Start time of simulation [Myr]
    ro:             Radial location of the Sun [kpc]
    vo:             Rotational speed at the Sun's position [km/s]
    nomwpot:        Boolean control for whether to include base MW potential
    transwrap:

    Returns potential
    """

    if nomwpot:
        pot=[]
    else:
        pot=copy.deepcopy(potential.MWPotential2014)

    if bar:
        # Initiate bar potential
        #tform= -10. # bar periods
        tform=-1.0e9 #bar periods
        tsteady= 5. # bar periods, slowly grown
        omega= 1.3 # Pattern speed of the bar

        angle=25./180.*np.pi # Bar Angle
        length=5. # Length of the bar
        dp= potential.DehnenBarPotential(omegab=omega,rb=length/8.,Af=(1./75.),tform=tform,tsteady=tsteady,barphi=angle)

        pot.append(dp)

    if dw_arm:
        # Add density wave spiral arm potential (http://galpy.readthedocs.io/en/latest/reference/potentialspiralarms.html)
        density_wave_spiral_arm=potential.DehnenSmoothWrapperPotential(pot=potential.SpiralArmsPotential(N=4,amp=1.,phi_ref=np.pi/4.,alpha=np.deg2rad(12.),omega=0.79),tform=dp.tform())
        pot.append(density_wave_spiral_arm)

    elif trans_dwarm:
        #WIP - no finished - ignore
        #Motivated by De Simone et al. 2004
        #Authors vary number of spiral arms, number of waves, density of spiral waves and pitch angle manually 
        sigma_s=1. #(Galpy units)
        sigma_c=0.25 #(Galpy units)

        tstart=(tstart/1000.)/bovy_conversion.time_in_Gyr(ro=ro,vo=vo)
        tpot=[]

        for i in range(0,Nwaves):
            theta=np.random.random()*2.0*np.pi
            to=np.random.random()*(tstart-4.0*sigma_s)+2.0*sigma_s

            xc=np.random.normal(loc=1.,scale=sigma_c)
            omega=(potential.vcirc(potential.MWPotential2014,xc,ro=ro,vo=vo)/vo)/xc

            print(i,to,tstart,to/tstart,theta,omega)

            sp= potential.SpiralArmsPotential(N=2,amp=0.75,phi_ref=theta,alpha=np.deg2rad(25.),omega=omega)
            csp= potential.GaussianAmplitudeWrapperPotential(pot=sp,to=to,sigma=sigma_s)
            tpot.append(csp)
            
        if transwrap:
            pot.append(TransientWrapperPotential(amp=1.,pot=tpot,_init=True))
        else:
            pot=pot+tpot

    elif trans_arm:
        # Add a series of three transient winding spirals
        #to is peak (--> Initialize 3 per 500 Myr)
        #lifetime is gaussian as grows and decays
        #Loop over spirals to
        to=-12.9885852455
        lifetimeS=1.3
        betap=-0.1
        omega= 1.3 # Pattern speed of the bar

        sp= potential.SpiralArmsPotential(N=2,amp=0.75,phi_ref=25./180.*np.pi,alpha=np.deg2rad(25.))

        if tstart==None:
            
            csp= potential.GaussianAmplitudeWrapperPotential(\
                   pot=potential.CorotatingRotationWrapperPotential(pot=sp,vpo=1.,to=to,beta=betap,pa=omega*to),
                   to=to,sigma=lifetimeS)
            pot.append(csp)
            csp3= potential.GaussianAmplitudeWrapperPotential(\
                   pot=potential.CorotatingRotationWrapperPotential(pot=sp,vpo=1.,to=to/2.,beta=betap,pa=omega*to),
                   to=to/2.,sigma=lifetimeS)
            pot.append(csp3)
            csp5= potential.GaussianAmplitudeWrapperPotential(\
                   pot=potential.CorotatingRotationWrapperPotential(pot=sp,vpo=1.,to=0.,beta=betap,pa=omega*to),
                   to=0.,sigma=lifetimeS)
            pot.append(csp5)
        else:
            #Assume tstart is in Myr and setup so there are three transient winding spirals per 500 Myr
            tstart=(tstart/1000.)/bovy_conversion.time_in_Gyr(ro=ro,vo=vo)
            t=0.
            tpot=[]
            while t > tstart:
                csp= potential.GaussianAmplitudeWrapperPotential(\
                       pot=potential.CorotatingRotationWrapperPotential(pot=sp,vpo=1.,to=t,beta=betap,pa=omega*to),
                       to=t,sigma=lifetimeS)
                tpot.append(csp)
                t+=to/2.
            
            if transwrap:
                pot.append(TransientWrapperPotential(amp=1.,pot=tpot,_init=True))
            else:
                pot=pot+tpot
    elif steady_arm:
        to=-12.9885852455
        lifetimeS=1.3
        betap=-0.1
        omega= 1.3 # Pattern speed of the bar

        sp= potential.SpiralArmsPotential(N=2,amp=0.75,phi_ref=25./180.*np.pi,alpha=np.deg2rad(25.))
        pot.append(potential.CorotatingRotationWrapperPotential(pot=sp,vpo=1.,to=to,beta=betap,pa=omega*to))


    return pot

class TransientWrapperPotential(potential.Potential):
    def __init__(self,amp=1.,pot=None,ro=None,vo=None,_init=None,**kwargs):
        """
        NAME:

           __init__

        PURPOSE:

           initialize a WrapperPotential, a super-class for wrapper potentials

        INPUT:

           amp - amplitude to be applied to the potential (default: 1.)

           pot - Potential instance or list thereof; the amplitude of this will be grown by this wrapper

        OUTPUT:

           (none)

        HISTORY:

           2017-06-26 - Started - Bovy (UofT)

        """
        if not _init: return None # Don't run __init__ at the end of setup
        potential.Potential.__init__(self,amp=amp,ro=ro,vo=vo)
        self._pot= pot
        self.isNonAxi= potential._isNonAxi(self._pot)
        self.ro=ro
        self.vo=vo
        
        self.times=np.array([])
        self.tend=np.array([])
        for i in range(0,len(self._pot)):
            self.times=np.append(self.times,self._pot[i]._to-10.*np.sqrt(self._pot[i]._sigma2))
            self.tend=np.append(self.tend,self._pot[i]._to+10.*np.sqrt(self._pot[i]._sigma2))

        #print(self.times)
        #print(self.tend)

        self.isNonAxi = True

    def __getattr__(self,attribute):
        if attribute == '_evaluate' \
                or attribute == '_Rforce' or attribute == '_zforce' \
                or attribute == '_phiforce' \
                or attribute == '_R2deriv' or attribute == '_z2deriv' \
                or attribute == '_Rzderiv' or attribute == '_phi2deriv' \
                or attribute == '_Rphideriv' or attribute == '_dens' \
                or attribute == 'isNonAxi':
            return lambda R,Z,phi=0.,t=0.: \
                self._wrap(attribute,R,Z,phi=phi,t=t)
        else:
            return super(TransientWrapperPotential,self).__getattr__(attribute)

    def _wrap_pot_func(self,attribute):
        if attribute == '_evaluate':
            return potential.evaluatePotentials
        elif attribute == '_dens':
            return potential.evaluateDensities
        elif attribute == '_Rforce':
            return potential.evaluateRforces
        elif attribute == '_zforce':
            return potential.evaluatezforces
        elif attribute == '_phiforce':
            return potential.evaluatephiforces
        elif attribute == '_R2deriv':
            return potential.evaluateR2derivs
        elif attribute == '_z2deriv':
            return potential.evaluatez2derivs
        elif attribute == '_Rzderiv':
            return potential.evaluateRzderivs
        elif attribute == '_phi2deriv':
            return lambda p,R,Z,phi=0.,t=0.: \
                potential.evaluatePotentials(p,R,Z,phi=phi,t=t,dphi=2)
        elif attribute == '_Rphideriv':
            return lambda p,R,Z,phi=0.,t=0.: \
                potential.evaluatePotentials(p,R,Z,phi=phi,t=t,dR=1,dphi=1)
        else: #pragma: no cover
            raise AttributeError("Attribute %s not found in for this TransientWrapperPotential" % attribute)
    
    def _evaluate(self,R,z,phi=0.,t=0.):
        mindx=self.match_time(t)
        if len(mindx)>0:
            return potential.evaluatePotentials(self._pot[mindx[0]:mindx[-1]+1],R,z,phi=phi,t=t)
        else:
            return 0.

    def _Rforce(self,R,z,phi=0.,t=0.):
        mindx=self.match_time(t)
        if len(mindx)>0:
            return potential.evaluateRforces(self._pot[mindx[0]:mindx[-1]+1],R,z,phi=phi,t=t)
        else:
            return 0.

    def _zforce(self,R,z,phi=0.,t=0.):
        mindx=self.match_time(t)
        if len(mindx)>0:
            return potential.evaluatezforces(self._pot[mindx[0]:mindx[-1]+1],R,z,phi=phi,t=t)
        else:
            return 0.

    def _phiforce(self,R,z,phi=0.,t=0.):
        mindx=self.match_time(t)
        if len(mindx)>0:
            return potential.evaluatephiforces(self._pot[mindx[0]:mindx[-1]+1],R,z,phi=phi,t=t)
        else:
            return 0.

    def match_time(self,t):
        
        mindx=(self.times<=t) * (self.tend>=t)
        indx=np.linspace(0,len(self.times)-1,len(self.times),dtype=int)

        #print((self.times<=t))
        #print((self.tend>=t))
        #print('MINDX: ',np.ndarray.tolist(indx[mindx]))
        return np.ndarray.tolist(indx[mindx])

class clusterIntegrate(object):
    def __init__(self,ra,dec,distance,pm_ra_cosdec,pm_dec,rv,ages,
                 cluster_id,pot=potential.MWPotential2014,
                 potname='MW2014'):
        """
        CLUSTER PROPERTIES
        The following quantities are arrays that must share the same length,
        they give the properties of each cluster member in sequence.

        The following quantities are assumed to have astropy units already
        assigned to them.
        RA:             right ascension  
        dec:            declination 
        distance:       distances 
        pm_ra_cosdec:   proper motions in right ascension (times cos(dec))
        pm_dec:         proper motions in declination
        rv:             radial velocities
        ages:           ages 

        SIMULATION PARAMETERS
        pot:            base potential in which to integrate
        """
        # Store cluster identifier
        self.cid = cluster_id
        # Store potential identifier
        self.pid = potname
        # Create coordinate object with stellar positions
        self.ra = ra
        self.dec = dec
        self.dis = distance
        self.pm_ra_cosdec = pm_ra_cosdec
        self.pm_dec = pm_dec
        self.rv = rv
        self.ages = ages
        self.current_coords = SkyCoord(ra=self.ra, dec=self.dec, 
                                       distance=self.dis, 
                                       pm_ra_cosdec=self.pm_ra_cosdec, 
                                       pm_dec=self.pm_dec,
                                       radial_velocity=self.rv)
        self.nstars = len(self.ra)
        # Set global base potential
        self.pot = pot

        return None

    def reverse_integrate_star(self,starind,age):
        """
        Reverse integrate a star's position in the base potential.

        starind:    integer index of star to integrate
        age:        amount of time in the past to integrate
        """
        self.starind=starind
        # Create stellar orbit
        self.o = Orbit(self.current_coords[self.starind])
        # Create timesteps
        self.ts = -np.linspace(0,age.value,int(1e4))*u.Gyr
        if isinstance(self.pot,tuple):
            print('EVALUATING POTENTIAL')
            potargs = self.pot[1]
            potargs['tstart'] = -age.to(u.Myr).value
            pot = self.pot[0](**potargs)
        elif isinstance(self.pot,list):
            pot = self.pot
        # Export potential to AMUSE format
        self.galaxy_code = potential.to_amuse(pot, tgalpy = -age.value | units.Gyr, reverse=False)
        print('POTENTIAL',pot)
        # Integrate
        self.o.integrate(self.ts,pot)
        return None

    def Nbody_initial_conditions(self,N=1000,Mcluster=1000.,Rcluster=10.):
        """
        Set up the initial conditions for the cluster in the Nbody simuation.

        N:          number of particles
        Mcluster:   mass of the cluster (in solar masses)
        Rcluster:   radius of the cluster (in pc)

        """

        # Number of particles in the cluster
        self.N = N
        # Mass of the clusters
        self.Mcluster= Mcluster| units.MSun
        # Tidal radius of the cluster
        self.Rcluster = Rcluster | units.parsec
        # Initial positions
        self.Rinit=[self.o.x(self.ts[-1].value).value,
                    self.o.y(self.ts[-1].value).value,
                    self.o.z(self.ts[-1].value).value] | units.kpc
        # Initial velocities
        self.Vinit=[self.o.vx(self.ts[-1].value).value,
                    self.o.vy(self.ts[-1].value).value,
                    self.o.vz(self.ts[-1].value).value] | units.km/units.s
        return None

    def setup_cluster(self):
        """
        Initialize cluster in phase space, assuming Plummer distribution.
        """

        #Create star cluster with origin at 0,0,0 and no net velocity
        self.converter=nbody_system.nbody_to_si(self.Mcluster,self.Rcluster)
        self.stars=new_plummer_sphere(self.N,self.converter)

        #Place star cluster in Galactocentric position
        self.stars.x += self.Rinit[0]
        self.stars.y += self.Rinit[1]
        self.stars.z += self.Rinit[2]

        self.stars.vx += self.Vinit[0]
        self.stars.vy += self.Vinit[1]
        self.stars.vz += self.Vinit[2]

        return None

    def evolve_cluster_in_galaxy(self,dt=1.0, dtout=5.0, tend=1000.0,saveout=True,plotout=False):
        """
        Cycle through the timesteps and update cluster stars according to internal interactions and
        the influence of the external potential.
        
        dt          Timestep in Myr
        dtout       Output timestep (is this used)
        tend        Maximum time to simulate in Myr
        saveout     Boolean controller for saving output
        plotout     Boolean controller for generating plots during simulation (deprecated)
        """
        # Save initial cluster properties
        savename = 'cid{0}_pot{1}_N{2}_Mc{3}_Rc{4}.hdf5'.format(self.cid,self.pid,self.N,self.Mcluster.value_in(units.MSun),self.Rcluster.value_in(units.parsec))
        if saveout:
            # Read out positions
            xpos = self.stars.x.value_in(units.kpc)*u.kpc
            ypos = self.stars.y.value_in(units.kpc)*u.kpc
            zpos = self.stars.z.value_in(units.kpc)*u.kpc
            # Read out velocities
            vx = self.stars.vx.value_in(units.km/units.s)*(u.km/u.s)
            vy = self.stars.vy.value_in(units.km/units.s)*(u.km/u.s)
            vz = self.stars.vz.value_in(units.km/units.s)*(u.km/u.s)
            # Galactocentric position and velocity
            Rpos,phi,z = bovy_coords.rect_to_cyl(xpos,ypos,zpos) #R not being expressed with units
            Rpos = Rpos*u.kpc
            vR,vT,vz = bovy_coords.rect_to_cyl_vec(vx,vy,vz,xpos,ypos,zpos)
            # Calculate actions
            o = Orbit([Rpos,vR,vT,z,vz,phi])
            jr=o.jr(pot=potential.MWPotential2014,type='staeckel',delta=0.45,c=True)#,ro=r0,vo=v0)
            jphi=o.jp(pot=potential.MWPotential2014,type='staeckel',delta=0.45,c=True)#,ro=r0,vo=v0)
            jz=o.jz(pot=potential.MWPotential2014,type='staeckel',delta=0.45,c=True)#,ro=r0,vo=v0)
            # Create dataframe row, with columns for each star
            propkeys = ['xp','yp','zp','vx','vy','vz','Rp','ph','vR',
                        'vT','jr','jp','jz']
            proplist = [xpos.value.astype(float),ypos.value.astype(float),
                        zpos.value.astype(float),vx.value.astype(float),
                        vy.value.astype(float),vz.value.astype(float),
                        Rpos.value.astype(float),phi.value.astype(float),
                        vR.value.astype(float),vT.value.astype(float),
                        jr.value.astype(float),jphi.value.astype(float),
                        jz.value.astype(float)]
            # Create dataframe for initial timestep
            df = pd.DataFrame(data=np.array(proplist).T,index=np.arange(self.N),columns=propkeys)
            savekey = 'star{0}_t{2:04}_tt{1:04}'.format(self.starind,
                                                  int(tend),0)
            df.to_hdf(savename,savekey,mode='a',format='table',data_columns=True)
           
        #Setup star cluster simulation
        total = tend/dt
        pbar = tqdm(total+1)
        #Simulation end time
        tend = tend | units.Myr #np.round(np.max(ages).to(u.Myr).value) | units.Myr
        #Frequency of data output
        dtout=dtout | units.Myr
        time= dt | tend.unit
        #Frequency of star cluster gravity calculation
        dt = dt | units.Myr

        #Setup cluster
        cluster_code=BHTree(self.converter,number_of_workers=1)     #Change number of workers depending on CPUS available
        cluster_code.parameters.epsilon_squared = (3. | units.parsec)**2
        cluster_code.parameters.opening_angle=0.6
        cluster_code.parameters.timestep=dt
        cluster_code.particles.add_particles(self.stars)

        #Setup channels between stars particle dataset and the cluster code
        channel_from_stars_to_cluster_code=self.stars.new_channel_to(cluster_code.particles, attributes=["mass", "x", "y", "z", "vx", "vy", "vz"])    
        channel_from_cluster_code_to_stars=cluster_code.particles.new_channel_to(self.stars, attributes=["mass", "x", "y", "z", "vx", "vy", "vz"])

        #Setup gravity bridge
        gravity=bridge.Bridge(use_threading=False)
        #stars in cluster_code depend on gravity from galaxy_code
        gravity.add_system(cluster_code, (self.galaxy_code,))
        #galaxy_code still needs to be added to system so it evolves with time
        gravity.add_system(self.galaxy_code,)
        #Set how often to update external potential
        gravity.timestep = cluster_code.parameters.timestep/2.

        while time<tend:
            # Update tqdm bar
            pbar.update(1)
            # Save cluster properties for this timestep
            if saveout or plotout:
                channel_from_cluster_code_to_stars.copy()
                # Read out positions
                xpos = self.stars.x.value_in(units.kpc)*u.kpc
                ypos = self.stars.y.value_in(units.kpc)*u.kpc
                zpos = self.stars.z.value_in(units.kpc)*u.kpc
                # Read out velocities
                vx = self.stars.vx.value_in(units.km/units.s)*(u.km/u.s)
                vy = self.stars.vy.value_in(units.km/units.s)*(u.km/u.s)
                vz = self.stars.vz.value_in(units.km/units.s)*(u.km/u.s)
                # Galactocentric position and velocity
                Rpos,phi,z = bovy_coords.rect_to_cyl(xpos,ypos,zpos) #R not being expressed with units
                Rpos = Rpos*u.kpc
                vR,vT,vz = bovy_coords.rect_to_cyl_vec(vx,vy,vz,xpos,ypos,zpos)
                # Calculate actions
                o = Orbit([Rpos,vR,vT,z,vz,phi])
                jr=o.jr(pot=potential.MWPotential2014,type='staeckel',delta=0.45,c=True)#,ro=r0,vo=v0)
                jphi=o.jp(pot=potential.MWPotential2014,type='staeckel',delta=0.45,c=True)#,ro=r0,vo=v0)
                jz=o.jz(pot=potential.MWPotential2014,type='staeckel',delta=0.45,c=True)#,ro=r0,vo=v0)
                # Create dataframe row, with columns for each star
                propkeys = ['xp','yp','zp','vx','vy','vz','Rp','ph','vR',
                            'vT','jr','jp','jz']
                proplist = [xpos.value.astype(float),ypos.value.astype(float),
                            zpos.value.astype(float),vx.value.astype(float),
                            vy.value.astype(float),vz.value.astype(float),
                            Rpos.value.astype(float),phi.value.astype(float),
                            vR.value.astype(float),vT.value.astype(float),
                            jr.value.astype(float),jphi.value.astype(float),
                            jz.value.astype(float)]
                # Create dataframe for this timestep
                df = pd.DataFrame(data=np.array(proplist).T,index=np.arange(self.N),columns=propkeys)
                savekey = 'star{0}_t{2:04}_tt{1:04}'.format(self.starind,
                                                      int(tend.value_in(units.Myr)),
                                                      int(time.value_in(units.Myr)))
                df.to_hdf(savename,savekey,mode='a',format='table',data_columns=True)
                
            # Continue time evolution
            gravity.evolve_model(time+dt)

            #You need to copy stars from cluster_code to output or analyse:

            #channel_from_cluster_code_to_stars.copy()
            #Output/Analyse
                
            #If you edited the stars particle set, lets say to remove stars from the array because they have 
            #been kicked far from the cluster, you need to copy the array back to cluster_code:
                
            #channel_from_stars_to_cluster_code.copy()
            
            time = gravity.model_time

        #Copy back to stars for final dataset
        channel_from_cluster_code_to_stars.copy()

        # Read out positions
        xpos = self.stars.x.value_in(units.kpc)*u.kpc
        ypos = self.stars.y.value_in(units.kpc)*u.kpc
        zpos = self.stars.z.value_in(units.kpc)*u.kpc
        # Read out velocities
        vx = self.stars.vx.value_in(units.km/units.s)*(u.km/u.s)
        vy = self.stars.vy.value_in(units.km/units.s)*(u.km/u.s)
        vz = self.stars.vz.value_in(units.km/units.s)*(u.km/u.s)
        # Galactocentric position and velocity
        Rpos,phi,z = bovy_coords.rect_to_cyl(xpos,ypos,zpos) #R not being expressed with units
        Rpos = Rpos*u.kpc
        vR,vT,vz = bovy_coords.rect_to_cyl_vec(vx,vy,vz,xpos,ypos,zpos)
        # Calculate actions
        o = Orbit([Rpos,vR,vT,z,vz,phi])
        jr=o.jr(pot=potential.MWPotential2014,type='staeckel',delta=0.45,c=True)#,ro=r0,vo=v0)
        jphi=o.jp(pot=potential.MWPotential2014,type='staeckel',delta=0.45,c=True)#,ro=r0,vo=v0)
        jz=o.jz(pot=potential.MWPotential2014,type='staeckel',delta=0.45,c=True)#,ro=r0,vo=v0)
        # Create dataframe row, with columns for each star
        propkeys = ['xp','yp','zp','vx','vy','vz','Rp','ph','vR',
                    'vT','jr','jp','jz']
        proplist = [xpos.value.astype(float),ypos.value.astype(float),
                    zpos.value.astype(float),vx.value.astype(float),
                    vy.value.astype(float),vz.value.astype(float),
                    Rpos.value.astype(float),phi.value.astype(float),
                    vR.value.astype(float),vT.value.astype(float),
                    jr.value.astype(float),jphi.value.astype(float),
                    jz.value.astype(float)]
        # Create dataframe for this timestep
        df = pd.DataFrame(data=np.array(proplist).T,index=np.arange(self.N),columns=propkeys)
        savekey = 'star{0}_t{2:04}_tt{1:04}'.format(self.starind,
                                              int(tend.value_in(units.Myr)),
                                              int(time.value_in(units.Myr)))
        df.to_hdf(savename,savekey,mode='a',format='table',data_columns=True)


        gravity.stop()

        return None

# Units for each key
punit = {'xp':u.kpc,'yp':u.kpc,'zp':u.kpc,
         'vx':u.km/u.s,'vy':u.km/u.s,'vz':u.km/u.s,
         'Rp':u.kpc,'ph':u.rad,
         'vR':u.km/u.s,'vT':u.km/u.s,
         'jr':u.kpc*u.km/u.s,'jp':u.kpc*u.km/u.s,'jz':u.kpc*u.km/u.s}

# Plot labels for each key
plabs = {'xp':'x [kpc]','yp':'y [kpc]','zp':'z [kpc]',
         'vx':'vx [km/s]','vy':'vy [km/s]','vz':'vz [km/s]',
         'Rp':'R [kpc]','ph':'phi [rad]',
         'vR':'vR [km/s]','vT':'vT [km/s]',
         'jr':'jR [kpc km/s]','jp':'jphi [kpc km/s]','jz':'jz [kpc km/s]'}

# Plot bounds for each key (needs update)
bounds = {'xp':(-20,20),'yp':(-20,20),'zp':(-5,5),
          'vx':(-300,300),'vy':(-300,300),'vz':(-20,20),
          'Rp':(0,20),'ph':(-1,-1),
          'jr':(-1,-1),'jp':(-1,-1),'jz':(-1,-1)}


if __name__ == '__main__':
    print('write tests, write clean exit with ^C')

# 
