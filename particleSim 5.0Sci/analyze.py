import pandas as pd
import matplotlib.pyplot as plt

##get particle number, k constant, and height of spring assume that (nP, k, and h are the same for rk4 and Euler)
with open("particleInfo5.0.csv","r") as f:
    nP_line=f.readline()
    k_line=f.readline()
    h_line=f.readline()

    #get number of particles, k, and h
    nP=int(nP_line.split(":")[1].strip())
    k=double(k_line.split(":")[1].strip())
    h=double(h_line.split(":")[1].strip())

df=pd.read_csv("particleInfo5.0.csv",skiprows=3)

df2=pd.read_csv("EulerParticle.csv",skiprows=3)

mass=0.5
g=9.8

print(f"Loaded Simulation with {nP} particles. Spring Constant: {k}")

time_steps=[]
total_energies=[]

##RK4 validation

##for every timestep go through every row
for t,group in df.groupby('time(s)'):
    frame_total_energy=0.0

    for _,row in group.iterrows():
        ##convert y, vx, and vy to numbers w NaN error handling
        y=pd.to_numeric(row['y'], errors='coerce')
        vx=pd.to_numeric(row['vx'],errors='coerce')
        vy=pd.to_numeric(row['vy'],errors='coerce')

        ##if any of them are not numbers then continue
        if pd.isna(y) or pd.isna(vx) or pd.isna(vy):
            continue

        ##calculate kinetic, potential, and spring energy to add to total energy for that particle at a timestep
        ek=0.5*mass*((vx*vx)+(vy*vy))

        ep=mass*g*y

        if(y<=h):
            es=0.5*k*((y-h)*(y-h))
        else:
            es=0.0

        frame_total_energy+=ek+ep+es
  
    #add the time to the timestep and the total energy arrays    
    time_steps.append(t);
    total_energies.append(frame_total_energy)

##create a table of timesteps and total energies
energy_df=pd.DataFrame({
    'time':time_steps,
    'total_energy':total_energies
    
    })


initial_energy=double(energy_df['total_energy'].iloc[0])

##calculate energy drift % and put it in a new row of the data frame
energy_df['drift_percent']=((energy_df['total_energy']-initial_energy)/initial_energy)*100

##clear timesteps and total_energies list
time_steps = []
total_energies = []

##Euler validation

##for every timestep go through every row
for t,group in df2.groupby('time(s)'):
    frame_total_energy=0.0

    for _,row in group.iterrows():
        ##convert y, vx, and vy to numbers w NaN error handling
        y=pd.to_numeric(row['y'], errors='coerce')
        vx=pd.to_numeric(row['vx'],errors='coerce')
        vy=pd.to_numeric(row['vy'],errors='coerce')

        ##if any of them are not numbers then continue
        if pd.isna(y) or pd.isna(vx) or pd.isna(vy):
            continue

        ##calculate kinetic, potential, and spring energy to add to total energy for that particle at a timestep
        ek=0.5*mass*((vx*vx)+(vy*vy))

        ep=mass*g*y

        if(y<=h):
            es=0.5*k*((y-h)*(y-h))
        else:
            es=0.0

        frame_total_energy+=ek+ep+es
  
    #add the time to the timestep and the total energy arrays    
    time_steps.append(t);
    total_energies.append(frame_total_energy)

##create a table of timesteps and total energies
energy_df2=pd.DataFrame({
    'time':time_steps,
    'total_energy':total_energies
    
    })


initial_energy=double(energy_df2['total_energy'].iloc[0])

##calculate energy drift % and put it in a new row of the data frame
energy_df2['drift_percent']=((energy_df2['total_energy']-initial_energy)/initial_energy)*100




plt.figure(figsize=(9,5))

##plot the time and drift percent for Rk4 and Euler
plt.plot(energy_df['time'],energy_df['drift_percent'],color='green',linewidth=1.5,label='RK4 Energy Drift')
plt.plot(energy_df2['time'],energy_df2['drift_percent'],color='blue',linewidth=1.5,label='Euler Energy Drift')


plt.title("Particle System Energy Conservation Validation Over Time")
plt.xlabel('Time(s)')
plt.ylabel('Energy Drift(%)')

##make the ideal energy conservation flat because energy should always be conserved in an ideal world
plt.axhline(0,color='red',linestyle='--',alpha=0.6,label='Ideal Energy Conservation(0% drift)')

plt.legend()
plt.grid(True)

plt.show()
