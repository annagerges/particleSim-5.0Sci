import pandas as pd
import matplotlib.pyplot as plt

#collect number of particles
with open("particleInfo5.0.csv", "r") as f:
    nP_line=f.readline()

    nP=int(nP_line.split(":")[1].strip())

df=pd.read_csv("particleInfo5.0.csv", skiprows=3)

print("Loading Particle Trajectories")

target_particles=[]

num_chosen=0

n_part=0

##lets user pick as many particles to track as they want
while True:
    n_part=int(input(f"How many particle trajectories do you want to graph (1-{nP}): "))
    if(n_part<=0 or n_part>nP):
        print("Invalid Amount")
    else:
        break

#let user pick the particles they want to track 
while(num_chosen!=n_part):
    part_id=int(input(f"Choose a particle to track(1-{nP}):"))
    if part_id<=0 or part_id>nP:
        print("Choose an ID from 1-",nP)
    elif part_id in target_particles:
        print(f"You already chose that one earlier. Choose a different number.")
    else:
        num_chosen+=1;
        target_particles.append(part_id)


#create slots in the dictionary that stores coordinates for the target particles
particle_paths={pid:{'x':[],'y':[]} for pid in target_particles}

##go through every row
for _, row in df.iterrows():
    #extract particle number
    raw_pid=row['Particle num']

    try:
        ##check if particle num is valid and convert it to an integer
        p_id=int(float(raw_pid))
    except (ValueError, TypeError):
        p_id=None

    ##if the particle number is a number and is one of the target particles convert the number to an integer, the coordinates to number (with NaN error handling), and append it to the dictionary if the coordinates are numbers
    if p_id is not None and p_id in target_particles:
        px=pd.to_numeric(row['x'], errors='coerce')
        py=pd.to_numeric(row['y'], errors='coerce')

        if not (pd.isna(px) or pd.isna(py)):
            particle_paths[p_id]['x'].append(px)
            particle_paths[p_id]['y'].append(py)

###Debugging###
##print("Columns Pandas found:", df.columns.tolist())
##print("Data points collected:", {pid: len(coords['x']) for pid, coords in particle_paths.items()})


#create the graph
plt.figure(figsize=(8,6))

# plot only particles that actually have collected data points
for particle, coords in particle_paths.items():
    if len(coords['x']) > 0:
        plt.plot(coords['x'], coords['y'], label=f'Particle {particle}', linewidth=1.5, marker='o', markersize=3)
    else:
        print(f"Warning: Particle {particle} was selected, but no valid coordinate data was found in the file.")

plt.title("Sample Particle Trajectories")
plt.xlabel('X Position(m)')
plt.ylabel('Y Position(m)')
plt.legend()
plt.grid(True)
plt.show()

