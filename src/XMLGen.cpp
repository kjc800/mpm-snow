#pragma once

#include "snowSimulator.h"
#include "particle.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void SnowSimulator::outputXML(size_t frame_count)
{
	const string path = "../../../renders/frames/" + to_string(frame_count) + ".xml";
	ofstream frame(path);
	frame << "<?xml version='1.0' encoding='utf-8'?>\n\n";

	frame << "<scene version=\"0.5.0\">\n";

	// CREATE A PATHTRACER
	frame << "\t<integrator type=\"path\">\n"  
			 "\t\t<integer name=\"maxDepth\" value=\"10\"/>\n" 
			 "\t</integrator>\n\n";

	// WRITE LIGHT SOURCE INFORMATION (ENVIRONMENT MAP)
	frame << "\t<emitter type=\"envmap\" id=\"envmapLight\">\n"
		     "\t\t<string name =\"filename\" value=\"textures/envmap.exr\"/>\n"
			 "\t\t<transform name=\"toWorld\">\n"
			 "\t\t\t<rotate y=\"1\" angle=\"75\"/>\n"
			 "\t\t</transform>\n"
			 "\t\t<float name=\"scale\" value =\"2.75\"/>\n"
			 "\t</emitter>\n\n";

	// CREATE A PINHOLE CAMERA
	frame << "\t<sensor type=\"perspective\">\n"
			 "\t\t<transform name=\"toWorld\">\n"
			 "\t\t\t<lookat origin=\"0, 1, -3\" target=\"0, 0, 0\" up=\"0, 0, 1\"/>\n"
			 "\t\t</transform>\n"
			 "\t\t<sampler type=\"independent\">\n"
			 "\t\t\t<integer name=\"sampleCount\" value=\"64\"/>\n"
			 "\t\t</sampler>\n"
			 "\t\t<film type=\"hdrfilm\">\n"
			 "\t\t\t<integer name=\"width\" value=\"500\"/>\n"
			 "\t\t\t<integer name=\"height\" value=\"500\"/>\n"
			 "\t\t</film>\n"
			 "\t\t<float name=\"fov\" value=\"50\"/>\n"
			 "\t</sensor>\n\n";

	// DRAW THE FLOOR
	frame << "\t<shape type=\"rectangle\">\n"
			 "\t\t<transform name=\"toWorld\">\n"
			 "\t\t\t<scale value=\"1\"/>\n"
			 "\t\t\t<rotate angle=\"270\" x=\"1\"/>\n"
			 "\t\t\t<translate y=\"0.15\"/>\n"
			 "\t\t</transform>\n"
			 "\t\t<bsdf type=\"diffuse\">\n"
			 "\t\t\t<rgb name=\"reflectance\" value=\"0.30, 0.30, 0.30\"/>\n"
			 "\t\t</bsdf>\n"
			 "\t</shape>\n\n";
	
	// FLOOR DIFFUSE MATERIAL BSDF
	frame << "\t<bsdf type=\"diffuse\" id=\"diffmat\">\n"
			 "\t\t<rgb name=\"reflectance\" value=\"0.4 0.4 0.4\"/>\n"
			 "\t</bsdf>\n\n";

	// FLOOR OBJECT
	frame << "\t<shape type=\"serialized\" id=\"floor\">\n"
			 "\t\t<string name=\"filename\" value=\"objects/floor.serialized\"/>\n"
			 "\t\t<integer name=\"shapeIndex\" value=\"0\"/>\n"
			 "\t\t<transform name=\"toWorld\">\n"
			 "\t\t\t<rotate x=\"1\" angle=\"90\"/>\n"
			 "\t\t\t<rotate z=\"1\" angle=\"180\"/>\n"
			 "\t\t\t<translate x=\"0.5\" y=\"0\" z=\"0.5\"/>\n"
			 "\t\t\t<scale x=\"45\" y=\"10\" z=\"45\"/>\n"
			 "\t\t</transform>\n"
			 "\t\t<ref name=\"bsdf\" id=\"diffmat\"/>\n"
			 "\t</shape>\n\n";

	drawSnowParticles(frame);

	frame << "</scene>";

	frame.close();
}

void SnowSimulator::drawSnowParticles(ofstream& frame)
{
	for (Particle* p : particles)
	{
		string x = to_string(p->position[0]);
		string y = to_string(p->position[2]);
		string z = to_string(p->position[1]);

		frame << "\t<shape type=\"sphere\">\n"
				 "\t\t<transform name=\"toWorld\">\n"
				 "\t\t\t<scale value=\"0.0125\"/>\n"
				 "\t\t\t<rotate angle=\"45\" y=\"1\"/>\n"
				 "\t\t\t<translate x=\"" + x + "\" y=\"" + y + "\" z=\"" + z + "\"/>\n"
				 "\t\t</transform>\n"
				 "\t\t<bsdf type=\"diffuse\">\n"
				 "\t\t\t<rgb name = \"reflectance\" value=\"0.2, 0.25, 0.7\"/>\n"
				 "\t\t</bsdf>\n"
				 "\t</shape>\n\n";
	}
}

