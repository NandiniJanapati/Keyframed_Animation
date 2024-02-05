#include <iostream>
#include <vector>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/quaternion.hpp>

#include "Camera.h"
#include "GLSL.h"
#include "Program.h"
#include "MatrixStack.h"
#include "Shape.h"
#include "Helicopter.h"
#include "Keyframe.h"

using namespace std;

GLFWwindow *window; // Main application window
string RESOURCE_DIR = ""; // Where the resources are loaded from

int keyPresses[256] = {0}; // only for English keyboards!

shared_ptr<Program> prog;
shared_ptr<Camera> camera;
//shared_ptr<Shape> shape;
//---------
//vector<shared_ptr<Shape>> allshapes;
shared_ptr<Helicopter> helicopter;
// Control points

vector<glm::vec3> keyframePos;
vector<glm::quat> keyframeRots;
vector<shared_ptr<Keyframe>> keyframes;
vector<glm::vec3> cps;
vector<glm::quat> cquats;
bool displayKeyframes;

vector<pair<float, float> > usTable;
float tmax = 5;
float smax;

enum TimeControl {
	NO_ARC_LENGTH = 0,
	ARC_LENGTH,
	EASE_IN_OUT,
	MY_FUNCTION
};
TimeControl timectrl = NO_ARC_LENGTH;

static void error_callback(int error, const char *description)
{
	cerr << description << endl;
}

static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	}
	if (action == GLFW_PRESS) {
		switch (key) {
			case GLFW_KEY_K: //107
				//keyPresses[(unsigned)'k'];
				break;
			case GLFW_KEY_S: //115
				timectrl = (TimeControl) ((keyPresses[(unsigned)'s'] + 1) % 4);
				break;
			case GLFW_KEY_SPACE:
				if (!((keyPresses[(unsigned)' '] + 1) % 2)) {
					camera->rotations = glm::vec2(0.0, 0.0);
					camera->translations = glm::vec3(0, 0, -5);
				}
				else {
					camera->rotations = glm::vec2(0.0, 0.0);
					camera->translations = glm::vec3(0, 0, -5);
				}
		}
	}

}

static void char_callback(GLFWwindow *window, unsigned int key)
{
	keyPresses[key]++;
}

static void cursor_position_callback(GLFWwindow* window, double xmouse, double ymouse)
{
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if(state == GLFW_PRESS && !((keyPresses[(unsigned)' ']) % 2)) { //if the camera is not set to helicopter's
		camera->mouseMoved(xmouse, ymouse);
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	// Get the current mouse position.
	double xmouse, ymouse;
	glfwGetCursorPos(window, &xmouse, &ymouse);
	// Get current window size.
	int width, height;
	glfwGetWindowSize(window, &width, &height);
	if(action == GLFW_PRESS) {
		bool shift = mods & GLFW_MOD_SHIFT;
		bool ctrl  = mods & GLFW_MOD_CONTROL;
		bool alt   = mods & GLFW_MOD_ALT;
		camera->mouseClicked(xmouse, ymouse, shift, ctrl, alt);
	}
}

//uses cquats
void addQuat(glm::quat q1) {
	if (cquats.size() > 0) {
		glm::quat q0 = cquats.back();
		if (glm::dot(q0, q1) < 0) {
			q1 = -1.0f * q1;
		}
	}
	cquats.push_back(q1);
}

//uses usTable, smax, cps
void buildTable()
{
	usTable.clear();
	// INSERT CODE HERE
	int ncps = (int)cps.size();
	float u = 0.0f;
	float s = 0.0f;

	usTable.push_back(make_pair(u, s));

	glm::mat4 G; //you have to switch the g matrix depending on what ctrl points your looking at.
	G[0] = glm::vec4(cps[0], 0.0f);
	G[1] = glm::vec4(cps[0 + 1], 0.0f);
	G[2] = glm::vec4(cps[0 + 2], 0.0f);
	G[3] = glm::vec4(cps[0 + 3], 0.0f);

	glm::mat4 B;
	B[0] = glm::vec4(0, 2, 0, 0);
	B[1] = glm::vec4(-1, 0, 1, 0);
	B[2] = glm::vec4(2, -5, 4, -1);
	B[3] = glm::vec4(-1, 3, -3, 1);
	B = 0.5f * B;
	glm::vec4 uVec(1.0f, 0, 0 * 0, 0 * 0 * 0);

	glm::vec4 Pu0 = G * (B * uVec); //point 0 in delta s = || P(u1) - P(u0) ||

	while (u < (ncps - 3) && (abs(ncps - 3 - u) > 0.01)) {
		u += 0.1f; //increment u
		int leftmostCrtlPt = floor(u);
		float uhat = u - leftmostCrtlPt; //un-concatenized u value

		//making the g matrix
		G[0] = glm::vec4(cps[leftmostCrtlPt], 0.0f);
		G[1] = glm::vec4(cps[leftmostCrtlPt + 1], 0.0f);
		G[2] = glm::vec4(cps[leftmostCrtlPt + 2], 0.0f);
		G[3] = glm::vec4(cps[leftmostCrtlPt + 3], 0.0f);

		uVec = glm::vec4(1.0f, uhat, uhat * uhat, uhat * uhat * uhat);
		// Compute position at u
		glm::vec4 Pu1 = G * (B * uVec);
		s += glm::length(Pu1 - Pu0);
		usTable.push_back(make_pair(u, s));
		Pu0 = Pu1;
	}
	smax = usTable.back().second;

}

static void init()
{
	GLSL::checkVersion();

	// Set background color
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	// Enable z-buffer test
	glEnable(GL_DEPTH_TEST);

	keyPresses[(unsigned)'c'] = 1;

	prog = make_shared<Program>();
	prog->setShaderNames(RESOURCE_DIR + "phong_vert.glsl", RESOURCE_DIR + "phong_frag.glsl");
	prog->setVerbose(true);
	prog->init();
	prog->addUniform("P");
	prog->addUniform("MV");
	prog->addUniform("lightPos");
	prog->addUniform("ka");
	prog->addUniform("kd");
	prog->addUniform("ks");
	prog->addUniform("s");
	prog->addAttribute("aPos");
	prog->addAttribute("aNor");
	prog->setVerbose(false);

	camera = make_shared<Camera>();

	//for (int i = 0; i < 4; i++) {
	//	allshapes.push_back(make_shared<Shape>());
	//}
	////shape = make_shared<Shape>();
	//allshapes[0]->loadMesh(RESOURCE_DIR + "helicopter_body1.obj");
	//// ----
	//allshapes[1]->loadMesh(RESOURCE_DIR + "helicopter_body2.obj");
	//allshapes[2]->loadMesh(RESOURCE_DIR + "helicopter_prop1.obj");
	//allshapes[3]->loadMesh(RESOURCE_DIR + "helicopter_prop2.obj");
	//// ----
	//for (int i = 0; i < 4; i++) {
	//	allshapes[i]->init();
	//}
	//shape->init();

	helicopter = make_shared<Helicopter>(RESOURCE_DIR);


	keyframePos.push_back(glm::vec3(0.0f, 0.0f, 0.0f) * 0.5f);
	keyframePos.push_back(glm::vec3(-6.0f, 2.0f, 0.0f) * 0.5f);
	keyframePos.push_back(glm::vec3(-4.0f, 6.0f, 2.0f) * 0.5f);
	keyframePos.push_back(glm::vec3(4.0f, 6.0f, -2.0f) * 0.5f);
	keyframePos.push_back(glm::vec3(6.0f, 2.0f, 0.0f) * 0.5f);

	keyframeRots.push_back(glm::angleAxis(2.0f, glm::normalize(glm::vec3(1.0f, 1.0f, 1.0f))));
	keyframeRots.push_back(glm::angleAxis(3.0f, glm::normalize(glm::vec3(0.0f, 0.0f, 1.0f))));
	keyframeRots.push_back(glm::angleAxis(6.0f, glm::normalize(glm::vec3(0.0f, 0.0f, 1.0f))));
	keyframeRots.push_back(glm::angleAxis(2.5f, glm::normalize(glm::vec3(1.0f, 0.0f, 0.0f))));
	keyframeRots.push_back(glm::angleAxis(2.7f, glm::normalize(glm::vec3(1.0f, -1.0f, 0.0f))));

	for (int i = 0; i < 5; i++) {
		keyframes.push_back(make_shared<Keyframe>());
		keyframes[i]->setPos(keyframePos[i]);
		if (i > 0) {
			glm::quat q0 = keyframes[i-1]->rotation;
			glm::quat q1 = keyframeRots[i];
			if (glm::dot(q0, q1) < 0) {
				q1 = -1.0f * q1;
			}
			keyframeRots[i] = q1;
		}
		keyframes[i]->setRot(keyframeRots[i]);
	}


	if (keyframes.size() >= 4) {
		int ncps = keyframes.size() + 3;
		for (int i = 0; i < ncps; i++) {
			if ((i == 0) || (i > (ncps - 3))) {
				cps.push_back(keyframes[0]->position);
				addQuat(keyframes[0]->rotation);
			}
			else {
				cps.push_back(keyframes[i - 1]->position);
				addQuat(keyframes[i - 1]->rotation);
			}
		}
	}

	buildTable();

	/*
	cps.push_back(keyframes[0]->position); //cps 0 = keyframe 0
	cps.push_back(keyframes[0]->position); //cps 1 = keyframe 0
	cps.push_back(keyframes[1]->position); //cps 2 = keyframe 1
	cps.push_back(keyframes[2]->position); //cps 3 = keyframe 2
	cps.push_back(keyframes[3]->position); //cps 4 = keyframe 3
	cps.push_back(keyframes[4]->position); //cps 5 = keyframe 4
	cps.push_back(keyframes[0]->position); //cps 6 = keyframe 0
	cps.push_back(keyframes[0]->position); //cps 7 = keyframe 0


	addQuat(keyframes[0]->rotation); //cps 0 = keyframe 0
	addQuat(keyframes[0]->rotation); //cps 1 = keyframe 0
	addQuat(keyframes[1]->rotation); //cps 2 = keyframe 1
	addQuat(keyframes[2]->rotation); //cps 3 = keyframe 2
	addQuat(keyframes[3]->rotation); //cps 4 = keyframe 3
	addQuat(keyframes[4]->rotation); //cps 5 = keyframe 4
	addQuat(keyframes[0]->rotation); //cps 6 = keyframe 0
	addQuat(keyframes[0]->rotation); //cps 7 = keyframe 0
	*/

	/*
	cquats.push_back(keyframes[0]->rotation); //cps 0 = keyframe 0
	cquats.push_back(keyframes[0]->rotation); //cps 1 = keyframe 0
	cquats.push_back(keyframes[1]->rotation); //cps 2 = keyframe 1
	cquats.push_back(keyframes[2]->rotation); //cps 3 = keyframe 2
	cquats.push_back(keyframes[3]->rotation); //cps 4 = keyframe 3
	cquats.push_back(keyframes[4]->rotation); //cps 5 = keyframe 4
	cquats.push_back(keyframes[0]->rotation); //cps 6 = keyframe 0
	cquats.push_back(keyframes[0]->rotation); //cps 7 = keyframe 0
	*/

	// Initialize time.
	glfwSetTime(0.0);

	// If there were any OpenGL errors, this will print something.
	// You can intersperse this line in your code to find the exact location
	// of your OpenGL error.
	GLSL::checkError(GET_FILE_LINE);
}

//uses usTable
float s2u(float s)
{
	// INSERT CODE HERE
	if (s < usTable[0].second) {
		return usTable[0].first;
	}

	for (int i = 0; i < usTable.size(); i++) {
		if (s < usTable[i].second) {
			float s0 = usTable[i - 1].second;
			float s1 = usTable[i].second;
			float alpha = (s - s0) / (s1 - s0);

			float u0 = usTable[i - 1].first;
			float u1 = usTable[i].first;
			float u = ((1 - alpha) * u0) + (alpha * u1);
			return u;
		}
	}
	return usTable[usTable.size() - 1].first;
}

//uses keyframes, cps
void drawSpline() { 

	glLineWidth(2.0f);
	glColor3f(0.0f, 0.8f, 0.8f);
	glm::mat4 B; //B matrix
	glm::mat4 G; //G matrix
	B[0] = glm::vec4(0, 2, 0, 0);
	B[1] = glm::vec4(-1, 0, 1, 0);
	B[2] = glm::vec4(2, -5, 4, -1);
	B[3] = glm::vec4(-1, 3, -3, 1);
	B = 0.5f * B;

	float u = 0.0f;
	glBegin(GL_LINE_STRIP);
	int nKFs = keyframes.size();
	while ((u < nKFs + 0.05) && (abs(nKFs - u + 0.1) > 0.01)) {
		int k = floor(u);
		float uhat = u-k;
		for (int i = 0; i < 4; i++) {
			G[i] = glm::vec4(cps[i+k], 0.0f);
		}

		glm::vec4 uVec(1.0f, uhat, uhat * uhat, uhat * uhat * uhat);
		glm::vec4 p = G * (B * uVec);
		glVertex3f(p.x, p.y, p.z);// draw a point at this location
		u += 0.1f; //increment u
	}
	glEnd();


}

//uses tmax, timectrl, smax, cquats, cps
glm::mat4 interpolatePos(float t) 
{  
	//glm::mat4 interpolatePos(float u) 
	float alpha = std::fmod(0.5f * t, 1.0f);
	float umax = keyframes.size();

	float u = std::fmod(t, umax); //if timectrl == 0
	if (timectrl == ARC_LENGTH) { //if timectrl == 1
		float tNorm = std::fmod(t, tmax) / tmax;
		float sNorm = tNorm;
		float s = smax * sNorm;
		u = s2u(s);
	}
	else if (timectrl == EASE_IN_OUT) { //ease in out
		float tNorm = std::fmod(t, tmax) / tmax;
		float sNorm = (- 2 * tNorm * tNorm * tNorm) + (3 * tNorm * tNorm);
		float s = smax * sNorm;
		u = s2u(s);
	}
	else if (timectrl == MY_FUNCTION) {
		float tNorm = std::fmod(t, tmax) / tmax;
		glm::mat4 A;
		glm::vec4 b(0, 0.5, 0, 1);
		// Fill A and b
		A[0] = glm::vec4(0, 0.3 * 0.3 * 0.3, 0.7 * 0.7 * 0.7, 1);
		A[1] = glm::vec4(0, 0.3 * 0.3, 0.7 * 0.7, 1);
		A[2] = glm::vec4(0, 0.3, 0.7, 1);
		A[3] = glm::vec4(1, 1, 1, 1);
		// Solve for x
		glm::vec4 x = glm::inverse(A) * b;
		float sNorm = (x.x * tNorm * tNorm * tNorm) + (x.y * tNorm * tNorm) + (x.z * tNorm) + x.w;//(10.71428 * tNorm * tNorm * tNorm) + (-14.88095 * tNorm * tNorm) + (5.1666 * tNorm) + 0;
		float s = smax * sNorm;
		u = s2u(s);
	}

	glm::mat4 B; //B matrix
	glm::mat4x3 G; //G matrix coll by row for some rsn??
	glm::mat4 Gq; //G matrix for quats

	B[0] = glm::vec4(0, 2, 0, 0);
	B[1] = glm::vec4(-1, 0, 1, 0);
	B[2] = glm::vec4(2, -5, 4, -1);
	B[3] = glm::vec4(-1, 3, -3, 1);
	B = 0.5f * B;

	int k = floor(u);
	float uhat = u - k;
	for (int i = 0; i < 4; i++) {
		G[i] = glm::vec3(cps[i + k]);
		Gq[i] = glm::vec4(cquats[i + k].x, cquats[i + k].y, cquats[i + k].z, cquats[i + k].w);
	}

	glm::vec4 uVec(1.0f, uhat, uhat * uhat, uhat * uhat * uhat);
	glm::vec3 p = G * (B * uVec);

	glm::vec4 qVec = Gq * (B * uVec);
	glm::quat q(qVec[3], qVec[0], qVec[1], qVec[2]); // Constructor argument order: (w, x, y, z)
	glm::mat4 E = glm::mat4_cast(glm::normalize(q)); // Creates a rotation matrix

	E[3] = glm::vec4(p, 1.0f); // Puts the position into the last column

	return E;

}

void render()
{
	// Update time.
	double t = glfwGetTime();

	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	glViewport(0, 0, width, height);

	// Use the window size for camera.
	glfwGetWindowSize(window, &width, &height);
	camera->setAspect((float)width/(float)height);

	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(keyPresses[(unsigned)'c'] % 2) {
		glEnable(GL_CULL_FACE);
	} else {
		glDisable(GL_CULL_FACE);
	}
	if(keyPresses[(unsigned)'z'] % 2) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	auto P = make_shared<MatrixStack>();
	auto MV = make_shared<MatrixStack>();

	// Apply camera transforms
	P->pushMatrix();
	camera->applyProjectionMatrix(P);
	MV->pushMatrix();

	camera->applyViewMatrix(MV);

	//drawing the interpolated position over time
	glm::mat4 transform = interpolatePos(t);
	if (((keyPresses[(unsigned)' ']) % 2)) {
		glm::vec3 heliPos = transform[3];
		glm::mat4 view(1);
		view = glm::translate(view, heliPos);

		glm::mat3 heliRot = transform;
		glm::mat4 heliRot4x4 = heliRot;
		view = view * heliRot4x4;

		//view = glm::translate(view, glm::vec3(5,0,-5));
		view = glm::rotate(view, 1.5708f, glm::vec3(0, 1, 0));
		view = glm::inverse(view);

		//MV->translate(transform[3]);				--> (A) Translate to hel pos
		//MV->multMatrix((glm::mat3)transform);		--> (B)	Rotate to heli rotation
		MV->multMatrix(view);

	}

	prog->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	glUniform3f(prog->getUniform("kd"), 1.0f, 0.0f, 0.0f);
	//shape->draw(prog);

	//draw the keyframes
	if(keyPresses[(unsigned)'k'] % 2) {
		int nKFs = keyframes.size();
		for (int i = 0; i < nKFs; i++) {

			glm::mat4 R = glm::mat4_cast(keyframes[i]->rotation);

			MV->pushMatrix();
			MV->translate(keyframes[i]->position);
			MV->multMatrix(R);

			helicopter->drawHeli(prog, MV, P, 0);

			MV->popMatrix();

		}
	}

	//drawing the interpolated position over time
	/*
	float alpha = std::fmod(0.5f * t, 1.0f);
	float umax = keyframes.size();

	float u = std::fmod(t, umax); //if timectrl == 0
	if (timectrl == ARC_LENGTH) { //if timectrl == 1
		float tNorm = std::fmod(t, tmax) / tmax;
		float sNorm = tNorm;
		float s = smax * sNorm;
		u = s2u(s);
	}
	else if (timectrl == EASE_IN_OUT) { //ease in out
		float tNorm = std::fmod(t, tmax) / tmax;
		float sNorm = (- 2 * tNorm * tNorm * tNorm) + (3 * tNorm * tNorm);
		float s = smax * sNorm;
		u = s2u(s);
	}
	else if (timectrl == MY_FUNCTION) {
		float tNorm = std::fmod(t, tmax) / tmax;
		glm::mat4 A;
		glm::vec4 b(0, 0.5, 0, 1);
		// Fill A and b
		A[0] = glm::vec4(0, 0.3 * 0.3 * 0.3, 0.7 * 0.7 * 0.7, 1);
		A[1] = glm::vec4(0, 0.3 * 0.3, 0.7 * 0.7, 1);
		A[2] = glm::vec4(0, 0.3, 0.7, 1);
		A[3] = glm::vec4(1, 1, 1, 1);
		// Solve for x
		glm::vec4 x = glm::inverse(A) * b;
		float sNorm = (x.x * tNorm * tNorm * tNorm) + (x.y * tNorm * tNorm) + (x.z * tNorm) + x.w;//(10.71428 * tNorm * tNorm * tNorm) + (-14.88095 * tNorm * tNorm) + (5.1666 * tNorm) + 0;
		float s = smax * sNorm;
		u = s2u(s);
	}
	*/

	//glm::mat4 transform = interpolatePos(u);
	//glm::mat4 transform = interpolatePos(t); //moved to before the camera set up so I can use the transform on the camera too
	glm::vec3 positionToDraw = transform[3];
	glm::mat3 rotationToDraw = transform;
	MV->pushMatrix();
	MV->translate(positionToDraw);
	MV->multMatrix(rotationToDraw);
	helicopter->drawHeli(prog, MV, P, 5*t);
	MV->popMatrix();

	prog->unbind();

	// Draw the frame and the grid with OpenGL 1.x (no GLSL)

	// Setup the projection matrix
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadMatrixf(glm::value_ptr(P->topMatrix()));

	// Setup the modelview matrix
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixf(glm::value_ptr(MV->topMatrix()));

	// Draw frame
	glLineWidth(2);
	glBegin(GL_LINES);
	glColor3f(1, 0, 0);
	glVertex3f(0, 0, 0);
	glVertex3f(1, 0, 0);
	glColor3f(0, 1, 0);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 1, 0);
	glColor3f(0, 0, 1);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, 1);
	glEnd();
	glLineWidth(1);

	//draw keyframe locations
	glPointSize(5.0f);
	glColor3f(0.0f, 0.0f, 0.0f);
	glBegin(GL_POINTS);
	for (int i = 0; i < keyframes.size(); ++i) {
		glVertex3f(keyframes[i]->position.x, keyframes[i]->position.y, keyframes[i]->position.z);
	}
	glEnd();

	if (timectrl > 0) { //draw equally spaced points
		float ds = 0.5;
		glColor3f(1.0f, 0.0f, 0.0f);
		glPointSize(10.0f);
		glBegin(GL_POINTS);
		glm::mat4 B;
		B[0] = glm::vec4(0, 2, 0, 0);
		B[1] = glm::vec4(-1, 0, 1, 0);
		B[2] = glm::vec4(2, -5, 4, -1);
		B[3] = glm::vec4(-1, 3, -3, 1);
		B = 0.5f * B;
		for (float s = 0.0f; s < smax; s += ds) {
			// Convert from s to (concatenated) u
			float uu = s2u(s);
			// Convert from concatenated u to the usual u between 0 and 1.
			float kfloat;
			float uhat = std::modf(uu, &kfloat);
			// k is the index of the starting control point
			int k = (int)std::floor(kfloat);
			// Compute spline point at u
			glm::mat4 Gk;
			for (int i = 0; i < 4; ++i) {
				Gk[i] = glm::vec4(cps[k + i], 0.0f);
			}
			glm::vec4 uVec(1.0f, uhat, uhat * uhat, uhat * uhat * uhat);
			glm::vec3 P(Gk * (B * uVec));
			glVertex3fv(&P[0]);
		}
		glEnd();
	}

	//draw spline
	drawSpline();

	// Draw grid
	float gridSizeHalf = 20.0f;
	int gridNx = 40;
	int gridNz = 40;
	glLineWidth(1);
	glColor3f(0.8f, 0.8f, 0.8f);
	glBegin(GL_LINES);
	for(int i = 0; i < gridNx+1; ++i) {
		float alpha = i / (float)gridNx;
		float x = (1.0f - alpha) * (-gridSizeHalf) + alpha * gridSizeHalf;
		glVertex3f(x, 0, -gridSizeHalf);
		glVertex3f(x, 0,  gridSizeHalf);
	}
	for(int i = 0; i < gridNz+1; ++i) {
		float alpha = i / (float)gridNz;
		float z = (1.0f - alpha) * (-gridSizeHalf) + alpha * gridSizeHalf;
		glVertex3f(-gridSizeHalf, 0, z);
		glVertex3f( gridSizeHalf, 0, z);
	}
	glEnd();

	// Pop modelview matrix
	glPopMatrix();

	// Pop projection matrix
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	// Pop stacks
	MV->popMatrix();
	P->popMatrix();

	GLSL::checkError(GET_FILE_LINE);
}

int main(int argc, char **argv)
{
	if(argc < 2) {
		cout << "Please specify the resource directory." << endl;
		return 0;
	}
	RESOURCE_DIR = argv[1] + string("/");

	// Set error callback.
	glfwSetErrorCallback(error_callback);
	// Initialize the library.
	if(!glfwInit()) {
		return -1;
	}
	// Create a windowed mode window and its OpenGL context.
	window = glfwCreateWindow(640, 480, "YOUR NAME", NULL, NULL);
	if(!window) {
		glfwTerminate();
		return -1;
	}
	// Make the window's context current.
	glfwMakeContextCurrent(window);
	// Initialize GLEW.
	glewExperimental = true;
	if(glewInit() != GLEW_OK) {
		cerr << "Failed to initialize GLEW" << endl;
		return -1;
	}
	glGetError(); // A bug in glewInit() causes an error that we can safely ignore.
	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	// Set vsync.
	glfwSwapInterval(1);
	// Set keyboard callback.
	glfwSetKeyCallback(window, key_callback);
	// Set char callback.
	glfwSetCharCallback(window, char_callback);
	// Set cursor position callback.
	glfwSetCursorPosCallback(window, cursor_position_callback);
	// Set mouse button callback.
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	// Initialize scene.
	init();
	// Loop until the user closes the window.
	while(!glfwWindowShouldClose(window)) {
		if(!glfwGetWindowAttrib(window, GLFW_ICONIFIED)) {
			// Render scene.
			render();
			// Swap front and back buffers.
			glfwSwapBuffers(window);
		}
		// Poll for and process events.
		glfwPollEvents();
	}
	// Quit program.
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
