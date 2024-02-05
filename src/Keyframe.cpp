#include "Keyframe.h"

using namespace std;


Keyframe::Keyframe(glm::vec3 pos, glm::quat rot) {
	position = pos;
	rotation = rot;
}

Keyframe::~Keyframe() {

}

Keyframe::Keyframe() {

}


void Keyframe::setPos(glm::vec3 pos) {
	position = pos;
}

void Keyframe::setRot(glm::quat rot) {
	rotation = rot;
}