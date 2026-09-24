#ifndef PARTICLES_H
#define PARTICLES_H

#include<vector>
#include <fstream>
#include <unordered_map>
#include<string>

//Particle class to keep track of every particle's position, velocity, acceleration, and energy
class Particles {
private:
	double vx = 0.0f;
	double vy;

	double x;
	double y;

	double a = -9.8f;

	double mass = 0.5f;

	int cellRow;
	int cellCol;


public:

	//setter and getter vx
	void setVx(double num) {
		vx = num;
	}

	double getVx() const {
		return vx;
	}

	//setter and getter vy
	void setVy(double num) {
		vy = num;
	}

	double getVy() const {
		return vy;
	}

	//setter and getter x and y
	void setX(double num) {
		x = num;
	}

	double getX() const {
		return x;
	}

	void setY(double num) {
		y = num;
	}

	double getY() const {
		return y;
	}

	//setter and getter acceleration
	void setA(double num) {
		a = num;
	}

	double getA() const {
		return a;
	}

	//setter and getter mass (for when user can dynamically customize mass)
	void setMass(double num) {
		mass = num;
	}

	double getMass() const {
		return mass;
	}

	//setter and getter for grid coordinates

	void setRow(int r) {
		cellRow = r;
	}

	int getRow() {
		return cellRow;
	}

	void setCol(int c) {
		cellCol = c;
	}

	int getCol() {
		return cellCol;
	}

};


class Spring {
private:
	double k;
	double h = 6;
	double dampen;

public:
	double getK() const {
		return k;
	}

	void setK(double num) {
		k = num;
	}

	double getHeight() const {
		return h;
	}

	//setter and getter for damping coefficient
	void setDamp(double d)
	{
		dampen = d;
	}

	double getDamp()
	{
		return dampen;
	}

};

//functions
void updatePos(std::vector<Particles>&, Spring&);
void wallCollis(std::vector<Particles>&);
void fix(std::vector<Particles>&, std::unordered_map<int, std::vector<Particles*>>&, int, int);
void particleCollis(std::unordered_map<int, std::vector<Particles*>>&, int);
void csvDump(std::vector<Particles>&, std::string&, double);

#endif


