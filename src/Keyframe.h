#pragma  once
#ifndef KEYFRAME_H
#define KEYFRAME_H

#include <vector>
#include <string>
#include <memory>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/quaternion.hpp>


class Keyframe {
public:
	glm::vec3 position;
	glm::quat rotation;

	Keyframe(glm::vec3 pos, glm::quat rot);
	~Keyframe();
	Keyframe();

	void setPos(glm::vec3 pos);

	void setRot(glm::quat rot);

private:

};

#endif
