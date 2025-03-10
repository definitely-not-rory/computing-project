#IMPORTS
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.rcParams.update({'font.size': 14})

import matplotlib.font_manager as font_manager
plt.rcParams["axes.formatter.use_mathtext"]=True
plt.rcParams['font.family']='serif'
cmfont = font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
plt.rcParams['font.serif']=cmfont.get_name()
plt.rcParams['mathtext.fontset']='cm'
plt.rcParams['axes.unicode_minus']=False

#CONSTANTS (I'm lazy and don't want to google them everytime)
grav_const=np.float64(6.67430*(10**-11))
au = np.float64(1.496*10**11)
day_in_seconds=np.float64(24*60*60)
year_in_seconds=day_in_seconds*365.242374
sun_jupiter_distance=5.2*au

#VECTOR PROCESSING FUNCTIONS
def scalar_separation(pos1,pos2): #function to get scalar separation of two point/position vectors
    total_separation=0 #setting accumulating separation stat to 0
    if len(pos1)!=len(pos2): #checking vectors are of same length
        print('position vectors are of different lengths')
        return None
    else:
        current_dimension=0 #variable indicating which dimension is being processesd
        while current_dimension<len(pos1): #iterating through dimensions
            dimension_diff=pos2[current_dimension]-pos1[current_dimension]
            total_separation=total_separation+dimension_diff**2 #summing squared separated distances
            current_dimension+=1 #increment to next dimension
        total_separation=np.sqrt(total_separation) #square rooting and returning true separation
        return total_separation

def unit_vector(vector): #function to generate unit vector of a given input vector
    sum=0
    for i in vector:
        sum+=i**2
    sum=np.sqrt(sum)
    return vector/sum

def magnitude(vector):
    sum=0
    for i in vector:
        sum+=i**2
    return np.sqrt(sum)

#PARTICLE INTERACTIONS
def gravity(mass1,mass2,pos1,pos2,softening): #basic particle-particle gravity
    separation=pos2-pos1 #determine separation vector of 2 inputted particles
    scalar_force=(grav_const*mass1*mass2)/(((scalar_separation(pos1,pos2))**2+softening**2)**(3/2)) #calculate scalar force from F = G*m_1*m_2/(r_1-r_2)^2
    vector_force=scalar_force*separation #calculate vector force by generating unit vector of separation
    #vector_force=(grav_const*mass1*mass2)*unit_vector(separation)/(scalar_separation(pos1,pos2))**2
    return vector_force

def potential_energy(mass1,mass2,pos1,pos2,softening):
    separation=scalar_separation(pos1,pos2)
    gpe=-1*(mass1*mass2*grav_const)/(np.sqrt(separation**2+softening**2))
    return gpe

def centre_of_mass(particles):
    com_contributions = np.float64(0.0)
    total_mass = np.float64(0.0)
    for particle in particles:
        com_contributions = com_contributions +particle.mass * particle.pos
        total_mass = total_mass +particle.mass

    return com_contributions / total_mass

def momentum(particles):
    momentum=np.array([0.0,0.0,0.0])
    for particle in particles:
        momentum+=particle.mass*particle.velocity
    return momentum


#PARTICLE CLASS
class Particle:
    def __init__(self, mass, pos, velocity,force, potential, dt):
        self.mass=mass #assigning variables for each particle
        self.pos=pos
        self.velocity=velocity
        self.force=force
        self.dt=dt
        self.potential=potential
    
    def calculate_forces(self, particles, softening): #calculate all gravitational forces on a particle
        resultant_force=np.array([0.0,0.0,0.0])
        for particle in particles:
            if particle != self:
                resultant_force =resultant_force + gravity(self.mass,particle.mass,self.pos,particle.pos,softening)
        self.force=resultant_force
        return resultant_force

    def calculate_potential_energy(self,particles,softening): #calculate current gravitational potential on a particle
        total_potential=0.0
        for particle in particles:
            if particle != self:
                total_potential= total_potential+potential_energy(self.mass, particle.mass,self.pos,particle.pos,softening)
        self.potential=total_potential
        return total_potential
        
    def iterate(self, particles,force, softening): #add function to complete one full iteration step (originally had v and pos iterate in different steps (too hard, got rid))
        self.potential=self.calculate_potential_energy(particles,softening)
        self.velocity=self.velocity+(force/self.mass)*self.dt #iterate velocity first/do n-1/2 step
        self.pos = self.pos +self.velocity*self.dt #iterate position second/do n step
    
    #get functions for position, velocity and forces to make plotting easier
    def get_pos(self):
        return self.pos

    def get_velocity(self):
        return self.velocity
    
    def get_force(self):
        return self.force

    def get_potential(self):
        return self.potential

#MILESTONE PROJECT
def milestone(display_vals, plot_orbits, show_animation, earth_sun_separation, jupiter_sun_separation, energy_plots, difference_from_circular,com_plot):
    
    testing_dt=10*day_in_seconds #testing values in line with milestone brief for quick referencing
    testing_dt_in_days=testing_dt/day_in_seconds
    total_sim_time=365.242374*day_in_seconds*100
    testing_iterations=int(total_sim_time/testing_dt)

    #Calculations for velocity vectors at -1/2dt
    e_angle=(2*np.pi/365.242374)*(testing_dt_in_days)/2

    e_omega=2*np.pi/(365.242374*day_in_seconds)

    e_vel=au*e_omega

    e_vel_x=e_vel*np.cos(np.pi/2-e_angle)
    e_vel_y=e_vel*np.sin(np.pi/2-e_angle)

    j_angle=(2*np.pi/(11.86*365.242374))*(testing_dt_in_days)/2

    j_omega=2*np.pi/(11.86*365.242374*day_in_seconds)

    j_vel=5.2*au*j_omega

    j_vel_x=j_vel*np.cos(np.pi/2-j_angle)
    j_vel_y=j_vel*np.sin(np.pi/2-j_angle)


    earth_pos=np.array([au,0.0,0.0])
    sun_pos=np.array([0.0,0.0,0.0])
    jupiter_pos=np.array([5.2*au,0.0,0.0])

    earth_vel=np.array([e_vel_x,e_vel_y,0.0])
    sun_vel=np.array([0.0,0.0,0.0])
    jupiter_vel=np.array([j_vel_x,j_vel_y,0.0])

    earth=Particle(np.float64(5.972*10**24),earth_pos,earth_vel,np.array([0.0,0.0,0.0]),0,testing_dt) #setting earth, sun and jupiter particle classes to correct variables
    sun=Particle(np.float64(1.989*10**30),sun_pos,sun_vel,np.array([0.0,0.0,0.0]),0,testing_dt)
    jupiter=Particle(np.float64(1.898*10**27),jupiter_pos,jupiter_vel,np.array([0.0,0.0,0.0]),0,testing_dt)

    particles=[sun,jupiter,earth] #defining list of active particles
    initial_momentum=np.array([0.0,0.0,0.0])
    for particle in particles:
        if particle!=sun:
            initial_momentum=initial_momentum+particle.mass*particle.velocity
    
    sun.velocity=-1*initial_momentum/sun.mass

    initial_com=centre_of_mass(particles)

    for particle in particles:
        particle.pos=particle.pos-initial_com

    '''new_com=centre_of_mass(particles)
    
    
    print(new_com)

    print(sun.pos)
    print(earth.pos)
    print(jupiter.pos)

    for particle in particles:
        particle.pos=particle.pos-new_com'''

    com_vel=0
    total_mass=0
    for particle in particles:
        com_vel+=particle.mass*particle.velocity
        total_mass+=particle.mass
    
    com_vel=com_vel/total_mass

    for particle in particles:
        particle.velocity-=com_vel
    
    earths=np.zeros((testing_iterations+1,3)) #setting storage arrays to correct size to store positions of each particle
    suns=np.zeros((testing_iterations+1,3))
    jupiters=np.zeros((testing_iterations+1,3))

    earth_vels=np.zeros((testing_iterations+1,3))
    sun_vels=np.zeros((testing_iterations+1,3))
    jupiter_vels=np.zeros((testing_iterations+1,3))

    earth_potentials=np.zeros((testing_iterations+1))
    sun_potentials=np.zeros((testing_iterations+1))
    jupiter_potentials=np.zeros((testing_iterations+1))

    coms=[]
    momenta=[]

    time=0 #setting t=0 to start simulation
    for i in range(testing_iterations+1): #iteration loop
        counter=0 #counting variable to track where in list 'particles iteration cycle is, resets to 0 at the end of one iteration of each particle
        coms.append(centre_of_mass(particles))
        momenta.append(momentum(particles))
        sun_forces=sun.calculate_forces(particles,0)
        jupiter_forces=jupiter.calculate_forces(particles,0)
        earth_forces=earth.calculate_forces(particles,0)
        for particle in particles:
            if counter==2: #checking which particle is currently selected using counting variable
                particle.iterate(particles,earth_forces,0) #do iteration with selected softening
                current_pos=particle.get_pos() #stores current position as local variable
                current_velocity=particle.get_velocity() #stores current velocity as local variable
                current_potential=particle.get_potential()
                name='Earth' #setting name if displaying all values of position is needed
                earths[int(time/testing_dt_in_days)]=current_pos #storing current position in array
                earth_vels[int(time/testing_dt_in_days)]=current_velocity #storing current velocity in array
                earth_potentials[int(time/testing_dt_in_days)]=current_potential   
            elif counter==0:
                particle.iterate(particles,sun_forces,0) #do iteration with selected softening
                current_pos=particle.get_pos() #stores current position as local variable
                current_velocity=particle.get_velocity() #stores current velocity as local variable
                current_potential=particle.get_potential()
                name='Sun'
                suns[int(time/testing_dt_in_days)]=current_pos 
                sun_vels[int(time/testing_dt_in_days)]=current_velocity
                sun_potentials[int(time/testing_dt_in_days)]=current_potential
            else:
                name='Jupiter'
                particle.iterate(particles,jupiter_forces,0) #do iteration with selected softening
                current_pos=particle.get_pos() #stores current position as local variable
                current_velocity=particle.get_velocity() #stores current velocity as local variable
                current_potential=particle.get_potential()
                jupiters[int(time/testing_dt_in_days)]=current_pos
                jupiter_vels[int(time/testing_dt_in_days)]=current_velocity 
                jupiter_potentials[int(time/testing_dt_in_days)]=current_potential

            if display_vals==True: #if requested to display all position values in function input, print each dt cycles position for each particle
                print('~~~ Day '+str(time)+' ~~~')
                print(name+': '+str(current_pos))
            counter+=1 #increment particle counter by 1 to move to next in list
            
        time+=testing_dt_in_days #advance by timestep (note: time is purely a cosmetic variable, testing_dt is used in iterative simulations)

    if plot_orbits==True: #check if plot requested, if true, produce position plot
        fig, ax=plt.subplots()
        ax.set_aspect('equal')
        plt.title('')
        plt.plot(earths[:,0],earths[:,1])
        plt.plot(suns[:,0],suns[:,1],c='y')
        plt.plot(jupiters[:,0],jupiters[:,1])
        x_values=[]
        y_values=[]
        for com in coms:
            x_values.append(com[0])
            y_values.append(com[1])
        plt.plot(x_values,y_values,c='r')
        plt.show()
        
    if show_animation==True:

        fig, ax = plt.subplots() #generate square x-y grid to plot animation on
        ax.set_aspect('equal')
        fig.set_figheight(15)
        max_x = max([pos[0] for pos in jupiters]) + 2.1*10**11 #define axis limits
        max_y = max([pos[1] for pos in jupiters]) + 2.1*10**11
        min_x = min([pos[0] for pos in jupiters]) - 2.1*10**11
        min_y = min([pos[1] for pos in jupiters]) - 2.1*10**11
        
        plt.xlim(min_x, max_x) #set axis limits
        plt.ylim(min_y, max_y)
        ax.set_facecolor('black') #set black background

        
        x_earths = earths[:,0] #set initial earth position
        y_earths = earths[:,1]
        line_earths,=ax.plot(x_earths,y_earths) #generate orbit line plot
        scat_earths = ax.scatter(1, 0) #plot psuedo-initial earth position
        earth_label=ax.annotate('Earth',xy=(x_earths[0]+2*10**10,y_earths[0]+2*10**10),c='w') #Assign label to earth

        def animate_earth(i):
            line_earths.set_data(x_earths[:i],y_earths[:i]) #update orbit line data
            scat_earths.set_offsets((x_earths[i], y_earths[i])) #for a given point 'i' in earths array, set earth plot to that value
            earth_label.xy = (x_earths[i]+2*10**10,y_earths[i]+10**10) #update label positioning for earth
            earth_label.set_position((x_earths[i]+2*10**10,y_earths[i]+2*10**10))
            return (scat_earths,),(line_earths,) #return as unpacked array so can be plotted
            

        ani_earth = animation.FuncAnimation(fig, animate_earth, repeat=True, frames=len(x_earths) - 1, interval=10) #animate earth movement
        
        x_jupiters = jupiters[:,0]
        y_jupiters = jupiters[:,1]
        line_jupiters,=ax.plot(x_jupiters,y_jupiters)
        scat_jupiters = ax.scatter(1, 0) #all same code as earth but for jupiter
        jupiter_label=ax.annotate('Jupiter',xy=(x_jupiters[0]+2*10**10,y_jupiters[0]+2*10**10),c='w') 

        def animate_jupiter(i):
            line_jupiters.set_data(x_jupiters[:i],y_jupiters[:i])
            scat_jupiters.set_offsets((x_jupiters[i], y_jupiters[i]))
            jupiter_label.xy = (x_jupiters[i]+2*10**10,y_jupiters[i]+2*10**10) 
            jupiter_label.set_position((x_jupiters[i]+2*10**10,y_jupiters[i]+2*10**10))
            return (scat_jupiters,),(line_jupiters,)

        ani_jupiter = animation.FuncAnimation(fig, animate_jupiter, repeat=True, frames=len(x_jupiters) - 1, interval=10)

        scat_suns = ax.scatter(1, 0, c='yellow') #all same code as earth but for sun
        x_suns = suns[:,0]
        y_suns = suns[:,1]
        sun_label=ax.annotate('Sun',xy=(x_suns[0]+2*10**10,y_suns[0]+2*10**10),c='white') 

        def animate_sun(i):
            scat_suns.set_offsets((x_suns[i], y_suns[i]))
            sun_label.xy = (x_suns[i]+2*10**10,y_suns[i]+2*10**10) 
            sun_label.set_position((x_suns[i]+2*10**10,y_suns[i]+2*10**10))
            return scat_suns,

        ani_sun = animation.FuncAnimation(fig, animate_sun, repeat=True, frames=len(x_suns) - 1, interval=1)
        
        plt.show()

    if earth_sun_separation==True:
        separation=np.zeros(len(earths)) #creates storage list for separation values
        for i in range(len(earths)):
            separation[i]=scalar_separation(earths[i],suns[i]) #calculate separation using previous function
        x_values=np.linspace(0,len(earths)*10,len(earths)) #generate x-axis in terms of days
        separation_au=separation/(1.496*10**11) #turn distance into AU
        fig = plt.figure()
        plt.plot(x_values,separation_au)
        plt.xlabel('Time (days)')
        plt.ylabel('Earth-Sun Separation (AU)')
        plt.show()

    if jupiter_sun_separation==True:
        separation=np.zeros(len(jupiters)) #creates storage list for separation values
        for i in range(len(jupiters)):
            separation[i]=scalar_separation(jupiters[i],suns[i]) #calculate separation using previous function
        x_values=np.linspace(0,len(jupiters)*10,len(jupiters)) #generate x-axis in terms of days
        separation_au=separation/(4.31*1.496*10**11) #turn distance into AU
        fig = plt.figure()
        plt.plot(x_values,separation_au)
        plt.xlabel('Time (days)')
        plt.ylabel('Jupiter-Sun Separation (AU)')
        plt.show()    
    
    if energy_plots==True:
        ke_earth=np.zeros(len(earths)) #generates storage array for kinetic energy
        ke_sun=np.zeros(len(earths)) #generates storage array for kinetic energy
        ke_jupiter=np.zeros(len(earths)) #generates storage array for kinetic energy
        
        count=0 #counting variable
        for velocity in earth_vels:
            ke_earth[count]=ke_earth[count]+1/2*earth.mass*(velocity[0]**2+velocity[1]**2+velocity[2]**2) #generate kinetic energy from velocity array
            count+=1
        
        count=0
        for velocity in sun_vels:
            ke_sun[count]=ke_sun[count]+1/2*sun.mass*(velocity[0]**2+velocity[1]**2+velocity[2]**2)
            count+=1

        count=0
        for velocity in jupiter_vels:
            ke_jupiter[count]=ke_jupiter[count]+1/2*jupiter.mass*(velocity[0]**2+velocity[1]**2+velocity[2]**2)
            count+=1

        x_values=np.linspace(0,len(earths)*10,len(earths)) #generate x-axis in terms of days

        ke=ke_earth+ke_jupiter+ke_sun

        pe=earth_potentials+jupiter_potentials+sun_potentials

        total_energy=ke_earth+ke_jupiter+ke_sun+earth_potentials+jupiter_potentials+sun_potentials

        change_in_energy=(total_energy[0]-total_energy[-1])/total_energy[0]

        print('Change in Energy: '+str(change_in_energy*100)+'%')

        fig, axs = plt.subplots(7)

        axs[0].plot(x_values,ke_earth)
        axs[0].set_xlabel('Time (days)')
        axs[0].set_ylabel('Kinetic Energy (J)')
        axs[1].plot(x_values,earth_potentials)
        axs[1].set_xlabel('Time (days)')
        axs[1].set_ylabel('Potential Energy (J)')

        axs[2].plot(x_values,ke_sun,c='y')
        axs[2].set_xlabel('Time (days)')
        axs[2].set_ylabel('Kinetic Energy (J)')
        axs[3].plot(x_values,sun_potentials,c='y')
        axs[3].set_xlabel('Time (days)')
        axs[3].set_ylabel('Potential Energy (J)')

        axs[4].plot(x_values,ke_jupiter,c='orange')
        axs[4].set_xlabel('Time (days)')
        axs[4].set_ylabel('Kinetic Energy (J)')
        axs[5].plot(x_values,jupiter_potentials,c='orange')
        axs[5].set_xlabel('Time (days)')
        axs[5].set_ylabel('Potential Energy (J)')

        axs[6].plot(x_values,total_energy,c='r')

        fig,axs=plt.subplots(3)
        axs[0].plot(x_values,ke)
        axs[0].set_xlabel('Time (days)')
        axs[0].set_ylabel('Kinetic Energy (J)')
        axs[1].plot(x_values,pe)
        axs[1].set_xlabel('Time (days)')
        axs[1].set_ylabel('Potential Energy (J)')
        axs[2].plot(x_values,total_energy)
        axs[2].set_xlabel('Time (days)')
        axs[2].set_ylabel('Total Energy (J)')
        '''

        fig.set_figheight(30)
        fig.tight_layout(pad=2)'''
        plt.show()
    
    if difference_from_circular==True:
        separation=np.zeros(len(earths)) #creates storage list for separation values
        for i in range(len(earths)):
            separation[i]=scalar_separation(earths[i],suns[i]) #calculate separation using previous function
        x_values=np.linspace(0,len(earths)*10,len(earths)) #generate x-axis in terms of days
        separation_au=separation/(1.496*10**11) #turn distance into AU
        radii_au=np.ones(len(earths)) #generate circular radius values for 1AU
        difference=separation_au-radii_au #calculate difference from circular orbit
        total_variation = abs(max(difference)-min(difference))
        percentage_variation = total_variation *100
        print('Circ Change: '+str(percentage_variation)+'%')
        fig=plt.figure()
        plt.plot(x_values,difference)
        plt.xlabel('Time (days)')
        plt.ylabel('Variation from Circular Earth Orbit (AU)')
        plt.show()

    if com_plot==True:
        x_values=[]
        y_values=[]
        for com in coms:
            x_values.append(com[0])
            y_values.append(com[1])
        plt.plot(x_values,y_values)
        plt.show()
        

milestone(False,True,False,False,False,True,True,False) 

'''

Extensions:
Solar System
Galaxies
Planetary Rings
Globular Cluster
Open Cluster



'''