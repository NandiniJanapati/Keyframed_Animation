#pragma  once
#ifndef HELICOPTER_H
#define HELICOPTER_H

#include "Shape.h"
#include "MatrixStack.h"
#include "Program.h"

#include <string>
#include <vector>
#include <memory>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/quaternion.hpp>


class Helicopter {
public:
	
	
	std::vector<std::shared_ptr<Shape>> allshapes;

	~Helicopter();
	Helicopter(const std::string& dir);

	void drawHeli(std::shared_ptr<Program>& prog, std::shared_ptr<MatrixStack>& MV, std::shared_ptr<MatrixStack>& P, float time);
	

private:
	
};

#endif