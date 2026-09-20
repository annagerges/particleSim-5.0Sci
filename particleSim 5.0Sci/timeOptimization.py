import pandas as pd
import matplotlib.pyplot as plt

#Read the CSV file
df = pd.read_csv('particleInfo4.0.csv', skiprows=3)

#a list of times
unique_times = df['time(s)'].unique()

#create a data table with steps rounded as an integer and simulation times
steps_df = pd.DataFrame({
    'total_steps': (unique_times / 0.001).round().astype(int),
    'simulation_time': unique_times
})

##do the sme thing for the other file
df2 = pd.read_csv('particleInfo.csv', skiprows=3)

unique_times = df2['time(s)'].unique()

steps_df2 = pd.DataFrame({
    'total_steps': (unique_times / 0.001).round().astype(int),
    'simulation_time': unique_times
})

# 4. Plot Steps vs. Time
plt.figure(figsize=(8, 5))

##lines
plt.plot(steps_df['total_steps'], steps_df['simulation_time'], label='Version 4.0 (Spatial Hash)', color='blue', linewidth=4)
plt.plot(steps_df2['total_steps'], steps_df2['simulation_time'], label='Version 3.0 (Vector Grid)', color='green')

plt.ylabel('Simulation Steps')
plt.xlabel('Simulation Time (s)')
plt.title('Simulation Steps vs. Time')
plt.grid(True)
plt.legend()
plt.show()