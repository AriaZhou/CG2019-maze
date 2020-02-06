#include "OpenGLWidget.h"
#include <iostream>
#include "MazeWidget.h"
#include <gl\gl.h>
#include <gl\GLU.h>
#include<QDebug>

OpenGLWidget::OpenGLWidget(QWidget *parent) : QGLWidget(parent)
{
	
	top_z = 1.5f;
	but_z = -1;

	QDir dir("Pic");
	if(dir.exists())
		pic_path = "Pic/";
	else
		pic_path = "../x64/Release/Pic/";
}
void OpenGLWidget::initializeGL()
{
	glClearColor(0,0,0,1);
	glEnable(GL_TEXTURE_2D);
	loadTexture2D(pic_path + "grass.png",grass_ID);
	loadTexture2D(pic_path + "sky.png",sky_ID);
	//loadTexture2D(pic_path + "wall.png", wall_ID);
}
void OpenGLWidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(MazeWidget::maze!=NULL)
	{
		//View 1
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glViewport(0, 0, MazeWidget::w / 2, (((MazeWidget::w / 2) * MazeWidget::maze->max_yp) / MazeWidget::maze->max_xp));
		//glViewport(0 , 0 , MazeWidget::w/2 , MazeWidget::h);
		glOrtho (-0.1, MazeWidget::maze->max_xp +0.1, -0.1 , MazeWidget::maze->max_yp +0.1, 0 , 10);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		Mini_Map();
		
		//View 2
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glViewport(MazeWidget::w / 2, 0, MazeWidget::w / 2, MazeWidget::h);
		//glViewport(MazeWidget::w/2,0, MazeWidget::w/2, (((MazeWidget::w / 2) * MazeWidget::maze->max_yp) / MazeWidget::maze->max_xp));
		/*gluPerspective ©w¸q³zµø
		//µø³¥¤j¤p, nearplane, farplane, distance
		//Note: You shouldn't use this function to get view matrix, otherwise you will get 0.
		*/
		//gluPerspective(MazeWidget::maze->viewer_fov, 1 , 0.01 , 200);

		/* gluLookAt
		//­ì¥»¬Û¾÷¦ì¸m
		//¬Ýªº¤è¦V
		//­þÃä¬O¤W­±
		//Note: You shouldn't use this function to get view matrix, otherwise you will get 0.
		*//*
		float viewerPosX = MazeWidget::maze->viewer_posn[Maze::X];
		float viewerPosY = MazeWidget::maze->viewer_posn[Maze::Y];
		float viewerPosZ = MazeWidget::maze->viewer_posn[Maze::Z];*/

		/*gluLookAt(viewerPosX, viewerPosZ, viewerPosY,
			viewerPosX + cos(degree_change(MazeWidget::maze->viewer_dir)), viewerPosZ, viewerPosY + sin(degree_change(MazeWidget::maze->viewer_dir)),
			0.0, -1.0, 0.0);*/
		glOrtho(-1, 1, -1, 1, -10, 10);
		//glOrtho(-0.1, MazeWidget::maze->max_xp + 0.1, -0.1, MazeWidget::maze->max_yp + 0.1, 0, 10);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		Map_3D();
	}
}
void OpenGLWidget::resizeGL(int w,int h)
{
}

//Draw Left Part
void OpenGLWidget::Mini_Map()	
{					
	glBegin(GL_LINES);

		float viewerPosX = MazeWidget::maze->viewer_posn[Maze::X];
		float viewerPosY = MazeWidget::maze->viewer_posn[Maze::Y];
		float viewerPosZ = MazeWidget::maze->viewer_posn[Maze::Z];

		for(int i = 0 ; i < (int)MazeWidget::maze->num_edges; i++)
		{
			float edgeStartX = MazeWidget::maze->edges[i]->endpoints[Edge::START]->posn[Vertex::X];
			float edgeStartY = MazeWidget::maze->edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
			float edgeEndX = MazeWidget::maze->edges[i]->endpoints[Edge::END]->posn[Vertex::X];
			float edgeEndY = MazeWidget::maze->edges[i]->endpoints[Edge::END]->posn[Vertex::Y];

			glColor3f(MazeWidget::maze->edges[i]->color[0] , MazeWidget::maze->edges[i]->color[1], MazeWidget::maze->edges[i]->color[2]);
			if(MazeWidget::maze->edges[i]->opaque)
			{
				glVertex2f(edgeStartX, edgeStartY);
				glVertex2f(edgeEndX, edgeEndY);
			}
		}

		//draw frustum
		float len = 10;
		glColor3f(1, 1, 1);
		glVertex2f(viewerPosX, viewerPosY);
		glVertex2f(viewerPosX + len * cos(degree_change(MazeWidget::maze->viewer_dir - MazeWidget::maze->viewer_fov / 2)),
			viewerPosY + len * sin(degree_change(MazeWidget::maze->viewer_dir - MazeWidget::maze->viewer_fov / 2)));

		glVertex2f(viewerPosX, viewerPosY);
		glVertex2f(viewerPosX + len * cos(degree_change(MazeWidget::maze->viewer_dir + MazeWidget::maze->viewer_fov / 2)),
			viewerPosY + len * sin(degree_change(MazeWidget::maze->viewer_dir + MazeWidget::maze->viewer_fov / 2)));

		/*
		glVertex2f(viewerPosX, viewerPosY);
		glVertex2f(viewerPosX + (MazeWidget::maze->max_xp) * len * cos(degree_change(MazeWidget::maze->viewer_dir - MazeWidget::maze->viewer_fov/2)) ,
			viewerPosY + (MazeWidget::maze->max_yp) * len * sin(degree_change(MazeWidget::maze->viewer_dir - MazeWidget::maze->viewer_fov/2)));

		glVertex2f(viewerPosX, viewerPosY);
		glVertex2f(viewerPosX + (MazeWidget::maze->max_xp) * len * cos(degree_change(MazeWidget::maze->viewer_dir + MazeWidget::maze->viewer_fov/2)) ,
			viewerPosY + (MazeWidget::maze->max_yp) * len *  sin(degree_change(MazeWidget::maze->viewer_dir + MazeWidget::maze->viewer_fov/2)));*/
	glEnd();
}

void OpenGLWidget::findVisibleEdgesForCell(Cell *currentCell, float upAngle, float downAngle, int enterEdge) {
	float viewerPosX = MazeWidget::maze->viewer_posn[Maze::X];
	float viewerPosY = MazeWidget::maze->viewer_posn[Maze::Y];
	float viewerPosZ = MazeWidget::maze->viewer_posn[Maze::Z];
	
	//shijiaoyi
	Edge** newEdges = new Edge*[4];
	for (int i = 0; i < 4; i++) {
		double edgeStartX = currentCell->edges[i]->endpoints[Edge::START]->posn[Vertex::X];
		double edgeStartY = currentCell->edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
		double edgeEndX = currentCell->edges[i]->endpoints[Edge::END]->posn[Vertex::X];
		double edgeEndY = currentCell->edges[i]->endpoints[Edge::END]->posn[Vertex::Y];

		double k1 = (double)(edgeEndY - edgeStartY) / (double)(edgeEndX - edgeStartX);
		double b1 = (double)edgeStartY - k1 * (double)edgeStartX;
		double k2 = (double)tan(upAngle);
		double b2 = (double)viewerPosY - k2 * (double)viewerPosX;
		double k3 = (double)tan(upAngle+3.14159/2.);
		double b3 = (double)viewerPosY - k3 * (double)viewerPosX;

		double sp = ((edgeStartX - viewerPosX)*sin(upAngle)
			+ (edgeStartY - viewerPosY)*(-cos(upAngle)));
		double ep = ((edgeEndX - viewerPosX)*sin(upAngle)
			+ (edgeEndY - viewerPosY)*(-cos(upAngle)));

		if (sp <= 0 && ep <= 0) {
			newEdges[i] = NULL;
		}
		else if (sp > 0 && ep < 0) {
			double ptix, ptiy;
			if (edgeStartX - edgeEndX == 0) {
				ptix = (double)edgeStartX;
				ptiy = k2 * ptix + b2;
			}
			else if (edgeEndY - edgeStartY == 0) {
				ptiy = (double)edgeStartY;
				ptix = (ptiy - b2) / k2;
			}
			else {
				ptix = (b2 - b1) / (k1 - k2);
				ptiy = k2 * ptix + b2;
			}
			newEdges[i] = new Edge(currentCell->edges[i]->index, currentCell->edges[i]->endpoints[Edge::START], new Vertex(Edge::END, ptix, ptiy),
				currentCell->edges[i]->color[0], currentCell->edges[i]->color[1], currentCell->edges[i]->color[2]);
		}
		else if (sp < 0 && ep > 0) {
			double ptix, ptiy;
			if (edgeStartX - edgeEndX == 0) {
				ptix = edgeStartX;
				ptiy = k2 * ptix + b2;
			}
			else if (edgeEndY - edgeStartY == 0) {
				ptiy = edgeStartY;
				ptix = (ptiy - b2) / k2;
			}
			else {
				ptix = (b2 - b1) / (k1 - k2);
				ptiy = k2 * ptix + b2;
			}
			newEdges[i] = new Edge(currentCell->edges[i]->index, new Vertex(Edge::START, ptix, ptiy), currentCell->edges[i]->endpoints[Edge::END],
				currentCell->edges[i]->color[0], currentCell->edges[i]->color[1], currentCell->edges[i]->color[2]);
		}
		else {
			newEdges[i] = currentCell->edges[i];
		}

		if (newEdges[i] == NULL)
			continue;

	}

	//shijiaoer
	for (int i = 0; i < 4; i++) {
		if (newEdges[i] == NULL)
			continue;

		double edgeStartX = newEdges[i]->endpoints[Edge::START]->posn[Vertex::X];
		double edgeStartY = newEdges[i]->endpoints[Edge::START]->posn[Vertex::Y];
		double edgeEndX = newEdges[i]->endpoints[Edge::END]->posn[Vertex::X];
		double edgeEndY = newEdges[i]->endpoints[Edge::END]->posn[Vertex::Y];


		double k1 = (double)(edgeEndY - edgeStartY) / (double)(edgeEndX - edgeStartX);
		double b1 = (double)edgeStartY - k1 * (double)edgeStartX;
		double k2 = (double)tan(downAngle);
		double b2 = (double)viewerPosY - k2 * (double)viewerPosX;
		double k3 = (double)tan(downAngle + 3.14159 / 2.);
		double b3 = (double)viewerPosY - k2 * (double)viewerPosX;

		double sp = ((edgeStartX - viewerPosX)*(sin(downAngle))
			+ (edgeStartY - viewerPosY)*(-cos(downAngle)));
		double ep = ((edgeEndX - viewerPosX)*(sin(downAngle))
			+ (edgeEndY - viewerPosY)*(-cos(downAngle)));

		//qDebug() << i << " " << sp << "*****" << ep;
		if (sp >= 0 && ep >= 0) {
			newEdges[i] = NULL;
		}
		else if (sp < 0 && ep > 0) {
			double ptix, ptiy;
			if (edgeStartX - edgeEndX == 0) {
				ptix = (double)edgeStartX;
				ptiy = k2 * ptix + b2;
			}
			else if (edgeEndY - edgeStartY == 0) {
				ptiy = (double)edgeStartY;
				ptix = (ptiy - b2) / k2;
			}
			else {
				ptix = (b2 - b1) / (k1 - k2);
				ptiy = k2 * ptix + b2;
			}
			newEdges[i] = new Edge(currentCell->edges[i]->index, newEdges[i]->endpoints[Edge::START], new Vertex(Edge::END, ptix, ptiy),
				currentCell->edges[i]->color[0], currentCell->edges[i]->color[1], currentCell->edges[i]->color[2]);

		}
		else if (sp > 0 && ep < 0) {
			double ptix, ptiy;
			if (edgeStartX - edgeEndX == 0) {
				ptix = edgeStartX;
				ptiy = k2 * ptix + b2;
			}
			else if (edgeEndY - edgeStartY == 0) {
				ptiy = edgeStartY;
				ptix = (ptiy - b2) / k2;
			}
			else {
				ptix = (b2 - b1) / (k1 - k2);
				ptiy = k2 * ptix + b2;
			}
			newEdges[i] = new Edge(currentCell->edges[i]->index, new Vertex(Edge::START, ptix, ptiy), newEdges[i]->endpoints[Edge::END],
				currentCell->edges[i]->color[0], currentCell->edges[i]->color[1], currentCell->edges[i]->color[2]);
		}
		else {
			newEdges[i] = newEdges[i];
		}

		if (newEdges[i] == NULL)
			continue;
	}

	/*
	for (int i = 0; i < 4; i++)
	{
		if (newEdges[i] == NULL)
			continue;
		
		float edgeStartX = newEdges[i]->endpoints[Edge::START]->posn[Vertex::X];
		float edgeStartY = newEdges[i]->endpoints[Edge::START]->posn[Vertex::Y];
		float edgeEndX = newEdges[i]->endpoints[Edge::END]->posn[Vertex::X];
		float edgeEndY = newEdges[i]->endpoints[Edge::END]->posn[Vertex::Y];

		glColor3f(newEdges[i]->color[0], newEdges[i]->color[1], newEdges[i]->color[2]);
		if (currentCell->edges[i]->opaque)
		{
			glVertex2f(edgeStartX, edgeStartY);
			glVertex2f(edgeEndX, edgeEndY);
		}
		currentCell->counter = 1;
	}
	*/
	transformation(newEdges);
	
	for (int i = 0; i < 4; i++)
	{
		if (newEdges[i] != NULL && currentCell->edges[i]->Neighbor(currentCell)) {
			float edgeStartX = newEdges[i]->endpoints[Edge::START]->posn[Vertex::X];
			float edgeStartY = newEdges[i]->endpoints[Edge::START]->posn[Vertex::Y];
			float edgeEndX = newEdges[i]->endpoints[Edge::END]->posn[Vertex::X];
			float edgeEndY = newEdges[i]->endpoints[Edge::END]->posn[Vertex::Y];

			if (!currentCell->edges[i]->opaque && currentCell->edges[i]->index != enterEdge) {
				float downAngleT = atan2(edgeStartY - viewerPosY, edgeStartX - viewerPosX);
				float upAngleT = atan2(edgeEndY - viewerPosY, edgeEndX - viewerPosX);
				
				if (currentCell->edges[i]->Point_Side(viewerPosX, viewerPosY) == Edge::LEFT) {
					upAngle = upAngleT;
					downAngle = downAngleT;
				}
				else {
					downAngle = upAngleT;
					upAngle = downAngleT;
				}
				
				findVisibleEdgesForCell(currentCell->edges[i]->Neighbor(currentCell), upAngle, downAngle, currentCell->edges[i]->index);
			}
		}
	}
}

void OpenGLWidget::transformation(Edge** visibleEdge) {

	float viewerPosX = MazeWidget::maze->viewer_posn[Maze::X];
	float viewerPosY = MazeWidget::maze->viewer_posn[Maze::Y];
	float viewerPosZ = MazeWidget::maze->viewer_posn[Maze::Z];
	float fov = MazeWidget::maze->viewer_fov;
	float direction = MazeWidget::maze->viewer_dir;
	float d = 1/tan(degree_change(MazeWidget::maze->viewer_fov / 2));     

	float rotateMatrix[16] = {cos(degree_change(direction)), 0, -sin(degree_change(direction)), 0,
							  0, 1, 0, 0,
							  -sin(degree_change(direction)),0, -cos(degree_change(direction)), 0,
							  0, 0, 0, 1 };

	float transformMatrix[16] = { 1, 0, 0, -viewerPosY,
								  0, 1, 0, 0,
								  0, 0, 1, -viewerPosX,
								  0, 0, 0, 1 };

	float projectionMatrix[16] = { 1, 0, 0, 0,
								   0, 1, 0, 0,
								   0, 0, 1, 0,
								   0, 0, 1/d, 0 };

	for (int i = 0; i < 4; i++) {
		if (visibleEdge[i] == NULL || !MazeWidget::maze->edges[visibleEdge[i]->index]->opaque)
			continue;

		glBegin(GL_QUADS);
		glColor3f(visibleEdge[i]->color[0], visibleEdge[i]->color[1], visibleEdge[i]->color[2]);
		float coorUpStart[4] = { visibleEdge[i]->endpoints[Edge::START]->posn[Vertex::Y], 1, visibleEdge[i]->endpoints[Edge::START]->posn[Vertex::X], 1 };
		float coorDownStart[4] = { visibleEdge[i]->endpoints[Edge::START]->posn[Vertex::Y], -1, visibleEdge[i]->endpoints[Edge::START]->posn[Vertex::X], 1};
		float coorUpEnd[4] = { visibleEdge[i]->endpoints[Edge::END]->posn[Vertex::Y], 1, visibleEdge[i]->endpoints[Edge::END]->posn[Vertex::X], 1 };
		float coorDownEnd[4] = { visibleEdge[i]->endpoints[Edge::END]->posn[Vertex::Y], -1, visibleEdge[i]->endpoints[Edge::END]->posn[Vertex::X], 1 };

		float* coorList[4] = { coorDownStart, coorUpStart, coorUpEnd, coorDownEnd};
		for (int a = 0; a < 4; a++) {
			
			/*
			| ux uy uz 0 | | 1 0 0 -ex |
			| vx vy vz 0 | | 0 1 0 -ey |
			| wx wy wz 0 | | 0 0 1 -ez |
			|  0  0  0 1 | | 0 0 0   1 |
			*/
			float RTMatrix[16] = { 0 };
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					for (int l = 0; l < 4; l++) {
						RTMatrix[j*4+k] += rotateMatrix[j * 4 + l] * transformMatrix[l * 4 + k];
					}
				}
			}

			/*
			| ux uy uz -u·ex | | x |
			| vx vy vz -v·ey | | y |
			| wx wy wz -w·ez | | z |
			|  0  0  0     1 | | 1 |
			*/
			float viewCoor[4] = { 0 };
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					viewCoor[j] += RTMatrix[j * 4 + k] * coorList[a][k];
				}
			}

			float screenCoor[4] = { 0 };
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					screenCoor[j] += projectionMatrix[j * 4 + k] * viewCoor[k];
				}
			}

			/*float dis = (MazeWidget::maze->edges[visibleEdge[i]->index]->endpoints[Edge::START]->posn[Vertex::X]
				- visibleEdge[i]->endpoints[Edge::START]->posn[Vertex::X]) / MazeWidget::maze->edges[visibleEdge[i]->index]->endpoints[Edge::START]->posn[Vertex::X];
			switch (a) {
			case 0:glTexCoord2f(1 - dis, 0.0f); break;
			case 1:glTexCoord2f(1 - dis, 1.0f); break;
			case 2:glTexCoord2f(1.0f, 1.0f); break;
			case 3:glTexCoord2f(1.0f, 0.0f); break;
			}*/
			
			glVertex2f((screenCoor[0] / screenCoor[3]), (screenCoor[1] / screenCoor[3]));
		}
		glEnd();
	}
}

//**********************************************************************
//
// * Draws the first-person view of the maze.
//   THIS IS THE FUINCTION YOU SHOULD MODIFY.
//
//Note: You must not use any openGL build-in function to set model matrix, view matrix and projection matrix.
//		ex: gluPerspective, gluLookAt, glTraslatef, glRotatef... etc.
//		Otherwise, You will get 0 !
//======================================================================
void OpenGLWidget::Map_3D()
{
	// µe¥kÃä°Ï¶ôªº©Ò¦³ªF¦è
	float viewerPosX = MazeWidget::maze->viewer_posn[Maze::X];
	float viewerPosY = MazeWidget::maze->viewer_posn[Maze::Y];
	float viewerPosZ = MazeWidget::maze->viewer_posn[Maze::Z];
	float upAngle = MazeWidget::maze->viewer_dir + MazeWidget::maze->viewer_fov / 2;
	float downAngle = MazeWidget::maze->viewer_dir - MazeWidget::maze->viewer_fov / 2;
	Edge **visibleEdges = new Edge*[MazeWidget::maze->num_edges];

	Cell *currentCell = MazeWidget::maze->cells[0];
	while (!(currentCell->Point_In_Cell(viewerPosX, viewerPosY, viewerPosZ, currentCell))) {
		if (currentCell == 0) {
			// The viewer is outside the top or bottom of the maze.
			throw new MazeException("Maze: View not in maze\n");
		}
	}

	//glEnable(GL_TEXTURE_2D);
	//glBindTexture(GL_TEXTURE_2D, grass_ID);
	//glBegin(GL_QUADS);
	//glTexCoord2f(0.0f, 0.0f); glVertex2f(-1, -1);
	//glTexCoord2f(0.0f, 1.0f); glVertex2f(-1, 0);
	//glTexCoord2f(0.25f, 1.0f); glVertex2f(1, 0);
	//glTexCoord2f(0.25f, 0.0f); glVertex2f(1, -1);
	//glEnd();

	//glBindTexture(GL_TEXTURE_2D, sky_ID);
	//glBegin(GL_QUADS);
	//glTexCoord2f(0.0f, 0.0f); glVertex2f(-1, 0);
	//glTexCoord2f(0.0f, 1.0f); glVertex2f(-1, 1);
	//glTexCoord2f(0.25f, 1.0f); glVertex2f(1, 1);
	//glTexCoord2f(0.25f, 0.0f); glVertex2f(1, 0);
	//glEnd();
	//glDisable(GL_TEXTURE_2D);
	//
	//glBindTexture(GL_TEXTURE_2D, wall_ID);
	

	//find visible edges cell by cell
	findVisibleEdgesForCell(currentCell, degree_change(upAngle), degree_change(downAngle), -1);
	for (int i = 0; i < MazeWidget::maze->num_cells; i++)
		MazeWidget::maze->cells[i]->counter = 0;

	/*­Y¦³¿³½ìªº¸Ü, ¥i¥H¬°¦aªO©Î°g®c¤W¶K¹Ï, ¦¹¶µ¥Ø¤£¼vÅTµû¤À*/
	//glBindTexture(GL_TEXTURE_2D, sky_ID); 
	
	// µe¶K¹Ï & ºâ UV
	glDisable(GL_TEXTURE_2D);
}
void OpenGLWidget::loadTexture2D(QString str,GLuint &textureID)
{
	glGenTextures(1, &textureID);
	glBindTexture(GL_TEXTURE_2D, textureID);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	
	QImage img(str);
	QImage opengl_grass = QGLWidget::convertToGLFormat(img);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, opengl_grass.width(), opengl_grass.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, opengl_grass.bits());
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glDisable(GL_TEXTURE_2D);
}
float OpenGLWidget::degree_change(float num)
{
	return num /180.0f * 3.14159f;
}