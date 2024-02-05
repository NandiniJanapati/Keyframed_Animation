#include "Helicopter.h"

using namespace std;


Helicopter::~Helicopter() {

}

Helicopter::Helicopter(const std::string& dir) {
	for (int i = 0; i < 4; i++) {
		allshapes.push_back(make_shared<Shape>());
	}

	allshapes[0]->loadMesh(dir + "helicopter_body1.obj");
	allshapes[1]->loadMesh(dir + "helicopter_body2.obj");
	allshapes[2]->loadMesh(dir + "helicopter_prop1.obj");
	allshapes[3]->loadMesh(dir + "helicopter_prop2.obj");

	for (int i = 0; i < 4; i++) {
		allshapes[i]->init();
	}
}

void Helicopter::drawHeli(std::shared_ptr<Program>& prog, std::shared_ptr<MatrixStack>& MV, std::shared_ptr<MatrixStack>& P, float t) {

	glUniform3f(prog->getUniform("kd"), 1.0f, 0.0f, 0.0f);
	MV->pushMatrix();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->popMatrix();
	allshapes[0]->draw(prog);//draw heli body/carriage

	glUniform3f(prog->getUniform("kd"), 1.0f, 1.0f, 0.0f);
	MV->pushMatrix();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->popMatrix();
	allshapes[1]->draw(prog);//draw heli tail


	glUniform3f(prog->getUniform("kd"), 0.5f, 0.5f, 0.5f);
	glm::vec3 topPropAxis(0.0f, 0.4819f, 0.0f);
	MV->pushMatrix();
	MV->rotate(t, glm::vec3(0, 1, 0));
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->popMatrix();
	allshapes[2]->draw(prog);//draw top propeller


	glm::vec3 tailPropAxis(0.6228f, 0.1179f, 0.1365f);
	MV->pushMatrix();
	MV->translate(tailPropAxis);
	MV->rotate(t, glm::vec3(0, 0, 1)); 
	MV->translate(-1.0f * tailPropAxis);
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->popMatrix();
	allshapes[3]->draw(prog);//draw tail propeller

}


