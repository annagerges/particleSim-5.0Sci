//cpp file for the program's functions

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include<string>
#include <unordered_map>
#include "Particles.h"

using namespace std;

const double dt = 0.001f;

//updates position using rk4 approximation
void updatePos(vector<Particles>& particles, Spring& s) {
	for (auto& p : particles) {
		// Determine if particle is on spring once at start of step
		bool onSpring = (p.getY() <= s.getHeight());

		// RK4 stages
		double k1y = p.getVy();

		//if it's on the spring change it's accelaration to k*dy/m-(dampening coeffiecient*velo/mass)-g and if not than keep accelaration to -9.8
		double k1v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - p.getY()) - (s.getDamp() * k1y) / p.getMass() - 9.8f : -9.8f;

		double k2y = p.getVy() + 0.5f * k1v * dt;
		double k2v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - (p.getY() + 0.5f * k1y * dt)) - (s.getDamp() * k2y) / p.getMass() - 9.8f : -9.8f;

		double k3y = p.getVy() + 0.5f * k2v * dt;
		double k3v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - (p.getY() + 0.5f * k2y * dt)) - (s.getDamp() * k3y) / p.getMass() - 9.8f : -9.8f;

		double k4y = p.getVy() + k3v * dt;
		double k4v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - (p.getY() + k3y * dt)) - (s.getDamp() * k4y) / p.getMass() - 9.8f : -9.8f;

		// Update
		p.setY(p.getY() + (dt / 6.0f) * (k1y + 2 * k2y + 2 * k3y + k4y));
		p.setVy(p.getVy() + (dt / 6.0f) * (k1v + 2 * k2v + 2 * k3v + k4v));
		p.setX(p.getX() + p.getVx() * dt);

	}
}


//checks if any of the particles collided with the wall
void  wallCollis(vector<Particles>& part) {
	for (int index = 0; index < part.size(); index++) {
		//if the particle is beyond or colliding with the left edge: reverse its x velocity and change its x position to 1
		if (part[index].getX() <= 0) {
			part[index].setVx(abs(part[index].getVx()));
			part[index].setX(1);

		}

		//if the particle is beyond or colliding with the right edge: reverse its x velocity and change its x position to 799
		else if (part[index].getX() >= 800) {
			part[index].setVx(abs(part[index].getVx()) * -1);
			part[index].setX(799);
		}

		//if the particle is beyond or colliding with the upper edge: reverse its y velocity and change its y position to 799
		if (part[index].getY() >= 800) {
			part[index].setVy(abs(part[index].getVy()) * -1);
			part[index].setY(799);
		}
		// if the particle is beyond or colliding with the floor: reverse its y velocity and clamp it
		else if (part[index].getY() <= 1) {
			part[index].setVy(abs(part[index].getVy()) * -1);
			part[index].setY(1);
		}
	}
}

//clears and places the updated particles into their correct keys
void fix(vector<Particles>& part, unordered_map<int, vector<Particles*>>& grid, int width, int nBox) {
	int row, col, cellKey;


	//assign every particle a grid cell based on its position
	for (int index = 0; index < part.size(); index++) {
		row = part[index].getY() / width;
		col = part[index].getX() / width;

		//safeguards to make sure that every particle gets assigned a valid cell 
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

		part[index].setRow(row);
		part[index].setCol(col);

		//cell (any paritcle w same row and col will have the same cell key)
		cellKey = (row * nBox) + col;

		grid[cellKey].push_back(&part[index]);
	}
}

//checks if particles collided with each other by passing the hashmap and looping through its particles
void particleCollis(unordered_map<int, vector<Particles*>>& grid, int nBox) {
	double d, dx, dy;

	for (const auto& pair : grid) {
		auto& cellKey = pair.first;

		int row = cellKey / nBox;
		int col = cellKey % nBox;

		const auto& cellParticles = pair.second;
		for (size_t start = 0; start < cellParticles.size(); start++) {
			for (size_t index = start + 1; index < cellParticles.size(); index++) {
				//compute dy and dx to figure out if they collided on x or y axis
				dy = cellParticles[index]->getY() - cellParticles[start]->getY();
				dx = cellParticles[index]->getX() - cellParticles[start]->getX();

				//absolute distance
				double absDx = abs(dx);
				double absDy = abs(dy);

				//if distance^2 is less than the particle radius^2: then they collided. Using squared to budget CPU resources and be accurate at the same time
				if (absDx * absDx + absDy * absDy < 100) {
					// Calculate relative velocity
					double dvx = cellParticles[index]->getVx() - cellParticles[start]->getVx();
					double dvy = cellParticles[index]->getVy() - cellParticles[start]->getVy();

					//only swap velo if they are moving towards eachother. Determines if they are pointing to eachother and acts accordingly
					if (dx * dvx + dy * dvy < 0) {
						if (absDx <= absDy) {
							double tempVy = cellParticles[start]->getVy();
							cellParticles[start]->setVy(cellParticles[index]->getVy());
							cellParticles[index]->setVy(tempVy);
						}
						else {
							double tempVx = cellParticles[start]->getVx();
							cellParticles[start]->setVx(cellParticles[index]->getVx());
							cellParticles[index]->setVx(tempVx);
						}
					}
				}
			}
		}

		//check cells that are foreward and the left and right lower diagonals instead of the traditional up, right, left, down so that the particles don't swap velocities twice
		int neighbors[4][2] = {
			{0,1},
			{ 1,1},
			{ 1,0 },
			{1,-1}
		};

		//loops through the vector of neighbor cell positions
		for (auto neighborJump : neighbors) {
			//row and col of neighbor cell
			int nRow = row + neighborJump[0];
			int nCol = col + neighborJump[1];

			if (nRow >= 0 && nRow < nBox && nCol >= 0 && nCol < nBox) {
				int neighborkey = nRow * nBox + nCol;

				//if the cell isn't the end iterator of the hash map see if the neigbor particles collided with any of the original cell neigbors
				if (grid.find(neighborkey) != grid.end()) {
					for (auto* ogCellPart : cellParticles) {

						//loop through the neighbors particles
						for (auto* neighborParts : grid[neighborkey]) {

							dy = neighborParts->getY() - ogCellPart->getY();
							dx = neighborParts->getX() - ogCellPart->getX();

							//absolute distance
							double absDx = abs(dx);
							double absDy = abs(dy);

							//if distance^2 is less than the particle radius^2: then they collided. Using squared to budget CPU resources and be accurate at the same time
							if (absDx * absDx + absDy * absDy < 100) {
								// Calculate relative velocity
								double dvx = neighborParts->getVx() - ogCellPart->getVx();
								double dvy = neighborParts->getVy() - ogCellPart->getVy();

								//only swap velo if they are moving towards eachother. Determines if they are pointing to eachother and acts accordingly
								if (dx * dvx + dy * dvy < 0) {
									if (absDx <= absDy) {
										double tempVy = neighborParts->getVy();
										neighborParts->setVy(ogCellPart->getVy());
										ogCellPart->setVy(tempVy);
									}
									else {
										double tempVx = neighborParts->getVx();
										neighborParts->setVx(ogCellPart->getVx());
										ogCellPart->setVx(tempVx);
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

//saves particle state into csv file
void csvDump(vector<Particles>& part, string& buffer, double time) {
	for (int index = 0; index < part.size(); index++) {
		buffer += to_string(index + 1) + "," + to_string(part[index].getX()) + "," + to_string(part[index].getY()) + "," + to_string(part[index].getVx())
			+ "," + to_string(part[index].getVy()) + "," + to_string(part[index].getA()) + "," + to_string(part[index].getRow()) + "," + to_string(part[index].getCol()) + ", " + to_string(time) + "\n";
	}
}

//checks accelaration at different points
//checks accelaration at different points
double checkAccel(double h, Spring& s, double m) {
	if (h <= s.getHeight()) {
		return (s.getK() / m) * (s.getHeight() - h) - 9.8f;
	}
	else {
		return -9.8f;
	}
}
