#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include <math.h>
#include <algorithm>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::vec2;
using glm::ivec2;
using glm::mat3;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 200; //500
const int SCREEN_HEIGHT = 200; //500
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
vec3 cameraPos(0, 0, -3.001);
float fl = SCREEN_WIDTH;
mat3 R; //The rotation matrix
float yaw = 0; // Yaw angle controlling camera rotation around y-axis
vec3 currentColor;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
struct Pixel
{
	int x;
	int y;
	float zinv;
	//6 Illumination:
	//vec3 illumination;//to compute the illumination on each pixel and each vertex
	//6.2
	vec3 pos3d;
};
//6 ILLUMINATION
struct Vertex
{
	vec3 position;
	//vec3 normal;
	//vec3 reflectance;
};

vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 14.1f*vec3(1, 1, 1); //ERROR IN THE LAB DOC-> 1.1f*vec3(1, 1, 1); IT'S OUT OF RANGE
vec3 indirectLightPowerPerArea = 0.5f*vec3(1, 1, 1);
//6.2
vec3 currentNormal;
vec3 currentReflectance;

//4 Filled Triangles
//vector<ivec2> leftPixels(ROWS);
//vector<ivec2> rightPixels(ROWS);

// ----------------------------------------------------------------------------
// FUNCTIONS
void Update();
void Draw();
void VertexShader(const Vertex& v, Pixel& p);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
//void DrawLineSDL(SDL_Surface* surface, Pixel a, Pixel b, vec3 color);
//void DrawPolygonEdges(const vector<vec3>& vertices);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void DrawPolygon(const vector<Vertex>& vertices);
//6 ILLUMINATION
void PixelShader(const Pixel& p);

int main(int argc, char* argv[])
{
	LoadTestModel(triangles);
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.
						/**
						//Interpolation example
						vector<Pixel> result(4);
						Pixel a(4, 2);
						Pixel b(1, 8);
						Interpolate(a, b, result);

						//ComputePolygonsRow() example
						vector<Pixel> vertexPixels(3);
						vertexPixels[0] = Pixel(10, 5);
						vertexPixels[1] = Pixel(5, 10);
						vertexPixels[2] = Pixel(15, 15);
						vector<Pixel> leftPixels;
						vector<Pixel> rightPixels;
						ComputePolygonRows(vertexPixels, leftPixels, rightPixels);**/

						/**
						for (int row = 0; row<leftPixels.size(); ++row)
						{
						cout << "Start: ("
						<< leftPixels[row].x << ","
						<< leftPixels[row].y << "). "
						<< "End: ("
						<< rightPixels[row].x << ","
						<< rightPixels[row].y << "). " << endl;
						}**/

	while (NoQuitMessageSDL())
	{
		Update();
		Draw();
	}

	SDL_SaveBMP(screen, "screenshot.bmp");
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	//cout << "Render time: " << dt << " ms." << endl;

	Uint8* keystate = SDL_GetKeyState(0);

	if (keystate[SDLK_UP])
		// Move camera forward
		cameraPos.z += 0.1;

	if (keystate[SDLK_DOWN])
		// Move camera backward
		cameraPos.z -= 0.1;

	if (keystate[SDLK_RIGHT])
		// Move camera to the right
		yaw -= 0.1;

	if (keystate[SDLK_LEFT])
		// Move camera to the left
		yaw += 0.1;

	if (keystate[SDLK_RSHIFT])
		;

	if (keystate[SDLK_RCTRL])
		;

	//LIGHT MOVEMENT (ILLUMINATION)
	if (keystate[SDLK_w])
	{
		//Move light forward
		//lightPos.z += 0.1;
	}
	if (keystate[SDLK_s])
	{
		//Move light backward
		//lightPos.z -= 0.1;
	}
	if (keystate[SDLK_a])
	{
		//Move light left
		//lightPos.x -= 0.1;
	}
	if (keystate[SDLK_d])
	{
		//Move light right
		//lightPos.x += 0.1;
	}
	if (keystate[SDLK_q])
	{
		//Move light up
		//lightPos.y += 0.1;
	}
	if (keystate[SDLK_e])
	{
		//Move light down
		//lightPos.y -= 0.1;
	}
	R = mat3(cos(yaw), 0, sin(yaw),
		0, 1, 0,
		-sin(yaw), 0, cos(yaw));

}

void Draw()
{
	SDL_FillRect(screen, 0, 0);

	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	//Clearing the depth buffer
	for (int y = 0; y<SCREEN_HEIGHT; ++y)
		for (int x = 0; x<SCREEN_WIDTH; ++x)
			depthBuffer[y][x] = 0;

	for (int i = 0; i<triangles.size(); ++i)
	{
		//6.2
		currentColor = triangles[i].color;
		vector<Vertex> vertices(3);
		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;
		currentReflectance = vec3(1, 1, 1);
		currentNormal = triangles[i].normal;

		DrawPolygon(vertices);

		/**6.1
		currentColor = triangles[i].color;
		vector<Vertex> vertices(3);
		vertices[0].position = triangles[i].v0;
		vertices[0].normal = triangles[i].normal;
		vertices[0].reflectance = vec3(1, 1, 1);

		vertices[1].position = triangles[i].v1;
		vertices[1].normal = triangles[i].normal;
		vertices[1].reflectance = vec3(1, 1, 1);

		vertices[2].position = triangles[i].v2;
		vertices[2].normal = triangles[i].normal;
		vertices[2].reflectance = vec3(1, 1, 1);

		DrawPolygon(vertices);**/

		/**5
		currentColor = triangles[i].color;
		vector<Vertex> vertices(3);
		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;
		DrawPolygon(vertices);**/


		/** 3/4
		vector<vec3> vertices(3);

		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;
		DrawPolygonEdges(vertices);
		**/

		// Add drawing: 
		/**2 Drawing Vertices
		for (int v = 0; v<3; ++v)
		{
		ivec2 projPos;
		VertexShader(vertices[v], projPos);
		vec3 color(1, 1, 1);
		PutPixelSDL(screen, projPos.x, projPos.y, color);
		}**/
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}
//It's called for each Vertex and computes a Pixel.
void VertexShader(const Vertex& v, Pixel& p) {
	vec3 res;
	res = ((v.position - cameraPos) * R);
	p.x = (fl * (res.x / res.z)) + SCREEN_WIDTH / 2;
	p.y = (fl * (res.y / res.z)) + SCREEN_HEIGHT / 2;
	p.zinv = 1 / res.z;//managing the inverse of z

					   /**6.1
					   // Compute pixel illumination
					   float r = glm::distance(lightPos, v.position);
					   vec3 rhat = glm::normalize(lightPos - v.position);
					   vec3 normal = v.normal;//normal
					   float prod = glm::dot(rhat, normal);
					   if (prod < 0)
					   prod = 0;
					   vec3 D = (lightPower * prod) / (4 * 3.1415926f*pow(r, 2));//equ 10
					   p.illumination = v.reflectance *(D + indirectLightPowerPerArea); //equ 11**/

					   //6.2
	p.pos3d = v.position;//store the 3D position of the Vertex to the corresponding variable in Pixel
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result)
{
	int N = result.size();
	float x = (float)(b.x - a.x) / float(glm::max(N - 1, 1));
	float y = (float)(b.y - a.y) / float(glm::max(N - 1, 1));
	float z = (float)(b.zinv - a.zinv) / float(glm::max(N - 1, 1));
	//6.1 light
	//vec3 light = vec3(b.illumination - a.illumination) / float(glm::max(N - 1, 1));
	vec3 position = vec3(b.pos3d - a.pos3d) / float(glm::max(N - 1, 1));
	Pixel current(a);

	for (int i = 0; i<N; ++i)
	{
		result[i] = current;
		current.x = x * i + a.x;
		current.y = y * i + a.y;
		current.zinv = z * i + a.zinv;
		current.pos3d += position; //light
	}
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels)
{
	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int max = glm::max(glm::max(vertexPixels[0].y, vertexPixels[1].y), vertexPixels[2].y);
	int min = glm::min(glm::min(vertexPixels[0].y, vertexPixels[1].y), vertexPixels[2].y);
	int ROWS = max - min + 1;

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.
	leftPixels.resize(ROWS);
	rightPixels.resize(ROWS);

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.

	/** Filling arrays with values representing the polygon by initializing the start arrays
	with really big values and the end array with really small values ***/
	for (int i = 0; i<ROWS; ++i)
	{
		leftPixels[i].x = +numeric_limits<int>::max();
		leftPixels[i].y = i + min;
		rightPixels[i].x = -numeric_limits<int>::max();
		rightPixels[i].y = i + min;
	}
	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	Pixel v; //actual vertex
	Pixel nextV; //next vertex
	int pixels;
	for (size_t i = 0; i < vertexPixels.size(); i++)
	{
		int j = (i + 1) % vertexPixels.size(); // to obtain the next vertex
		v = vertexPixels[i]; // obtaining the position of actual vertex
		nextV = vertexPixels[j];
		Pixel d;
		d.x = glm::abs(vertexPixels[i].x - vertexPixels[j].x);
		d.y = glm::abs(vertexPixels[i].y - vertexPixels[j].y);
		pixels = glm::max(d.y, d.y) + 1; // we use "y" because we are computing the rows
		vector<Pixel> line(pixels);

		//Interpolation between the two positions
		Interpolate(v, nextV, line);

		// Updating for each line the pixels
		for (size_t j = 0; j < line.size(); j++)
		{
			int tmp = line[j].y - min;
			if (line[j].x > rightPixels[tmp].x) {
				rightPixels[tmp].x = line[j].x;
				rightPixels[tmp].zinv = line[j].zinv;
				//6 Illumination
				//rightPixels[tmp].illumination = line[j].illumination;
				//6.2
				rightPixels[tmp].pos3d = line[j].pos3d;
			}
			if (line[j].x <= leftPixels[tmp].x) {
				leftPixels[tmp].x = line[j].x;
				leftPixels[tmp].zinv = line[j].zinv;
				//6 Illumination
				//leftPixels[tmp].illumination = line[j].illumination;
				//6.2
				leftPixels[tmp].pos3d = line[j].pos3d;
			}
		}
	}

}

void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels)
{
	for (size_t i = 0; i < leftPixels.size(); i++) //for all rows...
	{
		Pixel d;

		d.x = glm::abs(leftPixels[i].x - rightPixels[i].x);
		d.y = glm::abs(leftPixels[i].y - rightPixels[i].y);

		int pixels = glm::max(d.x, d.y) + 1;

		vector<Pixel> line(pixels);

		Interpolate(leftPixels[i], rightPixels[i], line);

		for (size_t i = 0; i < line.size(); i++)
		{
			// We check for the limitations of the depthBuffer (recommended by a friend)
			if ((line[i].zinv > depthBuffer[line[i].x][line[i].y]) && (line[i].x >= 0 && line[i].y >= 0) && (line[i].x < SCREEN_WIDTH && line[i].y < SCREEN_HEIGHT))
			{
				//PutPixelSDL(screen, line[i].x, line[i].y, currentColor);
				PixelShader(line[i]);
				depthBuffer[line[i].x][line[i].y] = line[i].zinv;
			}
		}
	}
}
//Copied as appeared in the documentation only changing ivec2 with Pixel
void DrawPolygon(const vector<Vertex>& vertices)
{
	int V = vertices.size();

	vector<Pixel> vertexPixels(V);
	for (int i = 0; i<V; ++i) {
		VertexShader(vertices[i], vertexPixels[i]);
	}
	vector<Pixel> leftPixels(3);
	vector<Pixel> rightPixels(3);
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawPolygonRows(leftPixels, rightPixels);
}
//Determines the final color of the image at that position
void PixelShader(const Pixel& p) {
	int x = p.x;
	int y = p.y;

	//6.1 
	//if (p.zinv > depthBuffer[y][x])
	//{
	//depthBuffer[y][x] = p.zinv;
	//PutPixelSDL(screen, x, y, currentColor); => without accessing the illumination property
	//PutPixelSDL(screen, x, y, p.illumination); => If we want a colorless illuminated scene!
	//PutPixelSDL(screen, x, y, p.illumination * currentColor);
	//}

	//6.2
	float r = glm::distance(lightPos, p.pos3d);
	vec3 rhat = glm::normalize(lightPos - p.pos3d);
	vec3 normal = currentNormal;
	float prod = glm::dot(rhat, normal);

	if (prod < 0)
	{
		prod = 0;
	}
	vec3 D = (lightPower * prod) / (4 * 3.1415926f*pow(r, 2));//equ 10
	vec3 illumination = D * currentReflectance + indirectLightPowerPerArea;//equ 11

	PutPixelSDL(screen, x, y, illumination * currentColor);
}

/**
void DrawLineSDL(SDL_Surface* surface, Pixel a, Pixel b, vec3 color) {
Pixel d.x = glm::abs(a.x - b.x);
Pixel d.y = glm::abs(a.y - b.y);
int pixels = glm::max(d.x, d.y) + 1;
vector<Pixel> line(pixels);
Interpolate(a, b, line);
//We loop to calculate all the pixels of a line
for (size_t i = 0; i < line.size(); i++)
{
PutPixelSDL(screen, line[i].x, line[i].y, color);
}
}

void DrawPolygonEdges(const vector<vec3>& vertices)
{
int V = vertices.size();
// Transform each vertex from 3D world position to 2D image position:
vector<Pixel> projectedVertices(V);
for (int i = 0; i<V; ++i)
{
VertexShader(vertices[i], projectedVertices[i]);
}
// Loop over all vertices and draw the edge from it to the next vertex:
for (int i = 0; i<V; ++i)
{
int j = (i + 1) % V; // The next vertex
vec3 color(1, 1, 1);
DrawLineSDL(screen, projectedVertices[i], projectedVertices[j],
color);
}
}**/