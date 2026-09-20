// Executes particle simulation where a user defined amount of particles (100-1000) move around indefinitely. 
// Uses vectors and hashmap based spatial partitioning to store particles and optimize collision detection. Euclidean approximation is used to simulate movement over a fixed time step.

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <limits>
#include <fstream>
#include <unordered_map>
#include "Particles.h"

using namespace std;

//timestep
const double dt = 0.001f;

int main() {

    int nP, row, col, stepCount = 0, cellKey;
    double accumulator, totalTime = 0;

    string buffer;
    buffer.reserve(128000);

    Spring s;


    random_device myEngine;

    //prompts user to enter valid amount of particles
    cout << "How many particles do you want(100-1000): ";
    cin >> nP;

    while (nP < 100 || nP>1000) {
        cout << "\nInvalid Input. Enter a valid amount(100-1000): ";
        cin >> nP;
    }

    cout << "Enter dampening coefficent(0-0.2): ";
    double damp;
    cin >> damp;

    while (damp < 0 || damp>0.2) {
        cout << "\nInvalid Input. Enter a valid input(0-0.2): ";
        cin >> damp;
    }

    s.setDamp(damp);

    double h;
    cout << "Set height of spring(m): ";
    cin >> h;

    while (h < 0) {
        cout << "\nInvalid Input. Enter a positive number: ";
        cin >> h;
    }

    int nBox = max(5, (int)sqrt(nP / 5));
    int width = 800 / nBox;

    //creates a vector of particles with the valid size the user specified to keep track of every particle
    vector<Particles>particles(nP);

    //sets up random number generator for particle position and velocity
    uniform_real_distribution<double>randPos(h + 1, 799);
    uniform_real_distribution<double>randVelo(1, 30);

    unordered_map <int, vector<Particles*>> hash;

    for (int index = 0; index < nP; index++) {

        //Randomly assigns x and y to be from 1 to 799 because having the user decide would be tedious
        particles[index].setX(randPos(myEngine));
        particles[index].setY(randPos(myEngine));

        //Randomly assigns vx and vy to be from 1 to 30
        particles[index].setVx(randVelo(myEngine));
        particles[index].setVy(randVelo(myEngine));

        //finds cell coordinates based on particle position
        row = particles[index].getY() / width;
        col = particles[index].getX() / width;

        if (row < 0) {
            row = 0;
        }
        else if (row >= nBox) {
            row = nBox - 1;
        }
        if (col < 0) {
            col = 0;
        }
        else if (col >= nBox) {
            col = nBox - 1;
        }

        particles[index].setRow(row);
        particles[index].setCol(col);

        //cell (any paritcle w same row and col will have the same cell key)
        cellKey = (row * nBox) + col;

        hash[cellKey].push_back(&particles[index]);

    }

    fstream file("particleInfo5.0.csv", ios::out);

    s.setK(((nP * 9.8 * particles[0].getMass()) / 0.2) * 4);

    //writing k and num of particles into the file for it to be analyzed using python but not seen
    file << "# nP: " << nP << "\n";
    file << "# k: " << s.getK() << "\n";
    file << "# h: " << s.getHeight() << "\n";

    //file rows
    file << "Particle num,x,y,vx,vy,ay,cellRow,cellCol,time(s)" << "\n";

    //sets accumulator to 0
    accumulator = 0;

    //start frame timer
    auto previousTime = chrono::high_resolution_clock::now();


    //infinite loop
    while (true) {
        //end of frame timer
        auto currentTime = chrono::high_resolution_clock::now();

        //duration of frame in seconds
        chrono::duration<double> frameDur = currentTime - previousTime;

        //assign the start of the next frame to be the time the previous one ended
        previousTime = currentTime;

        //increment accumulator by frame duration
        accumulator += frameDur.count();

        //so that the CPU doesn't have to run around and infinitely do stepCount to check for the 1000 particles.
        if (accumulator > 0.1f) {
            accumulator = 0.1f;
        }

        while (accumulator >= dt) {
            totalTime += dt;

            //update position, wall collision checks, and clear and update the hashmap 
            updatePos(particles, s);

            //debugging purposes
            //cout << particles[0].getY() << endl;

            wallCollis(particles);

            // Empty the vectors but keep their memory capacity in RAM
            for (auto& pair : hash) {
                pair.second.clear();
            }

            fix(particles, hash, width, nBox);
            particleCollis(hash, nBox);

            stepCount++;

            if (stepCount % 10 == 0) {
                if (buffer.size() > 100000) {
                    file << buffer;
                    buffer.clear();
                }
            }

            //decrement accumulator by the timestep
            accumulator -= dt;


        }
    }

    file.close();

    return 0;
}
