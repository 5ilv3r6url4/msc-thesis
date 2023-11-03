#include "utils/gui/util_obj_writer.h"
#include "utils/math/util_math_defs.h"

#include <iostream>
#include <QVector>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>

using namespace std;

void Util_Obj_Writer::writePolylines(QVector< QVector<point_t> > lines)
{
	// Get user input for file name
	QString filename = QFileDialog::getSaveFileName(0, "Save OBJ file", "./", "Wavefront OBJ (*.obj)");

	// Try to open the file
	QFile file(filename);
	if (!file.open(QIODevice::WriteOnly)) return;

	// The file is now open; create a stream to it
	QTextStream ts(&file);

	// Keeps track of starting vertex index for each spline
	int numVertices = 0;

	// Print object name
	ts << "o lineset" << endl;

	// Print each polyline
	for (int i=0; i<lines.size(); i++)
	{
		// Print vertices
		QVector<point_t> points = lines[i];
		for (int j=0; j<points.size(); j++)
		{
			point_t p = points[j];
			ts << "v " << p.x << " " << p.y << " " << p.z << endl;
		}

		// Print lines connecting vertices
		ts << "g polyline" << i << endl;
		for (int j=numVertices; j<numVertices+points.size()-1; j++)
		{
			ts << "l " << j+1 << " " << j+2 << endl;
		}

		// Note the number of vertices used so the next line's vertices
		// can be properly indexed
		numVertices += points.size();
	}

	file.close();
}

