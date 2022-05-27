#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"
#include "grid.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
    int num_height_points, float thickness) {
    this->width = width;
    this->height = height;
    this->num_width_points = num_width_points;
    this->num_height_points = num_height_points;
    this->thickness = thickness;

    buildGrid();
    buildClothMesh();
}

Cloth::~Cloth() {
    point_masses.clear();
    springs.clear();

    if (clothMesh) {
        delete clothMesh;
    }
}

bool isPinned(const vector<vector<int>>& pinned, int x, int y)
{
    for (auto coords : pinned)
    {
        if (coords[0] == x && coords[1] == y)
        {
            return true;
        }
    }
    return false;
}

void Cloth::buildGrid() {
    // TODO (Part 1): Build a grid of masses and springs.

    // Point Mass Creation
    double deltaX = width / (double)(num_width_points - 1);
    double deltaY = height / (double)(num_height_points - 1);

    double widthTrack = 0.0;
    double heightTrack = 0.0;

    for (int j = 0; j < num_height_points; j++)
    {
        for (int i = 0; i < num_width_points; i++)
        {
            Vector3D position;

            if (orientation == HORIZONTAL)
            {
                position = Vector3D(widthTrack, 1.0, heightTrack);
            }
            else
            {
                double randZ = (double((rand() % 2) - 1)) / 1000.0;
                position = Vector3D(widthTrack, heightTrack, randZ);
            }

            if (isPinned(pinned, i, j))
            {
                PointMass mass = PointMass(position, true);
                point_masses.emplace_back(mass);
            }
            else
            {
                PointMass mass = PointMass(position, false);
                point_masses.emplace_back(mass);
            }

            widthTrack += deltaX;
        }

        heightTrack += deltaY;
        widthTrack = 0;
    }

    // Spring Generation
    for (int j = 0; j < num_height_points; j++)
    {
        for (int i = 0; i < num_width_points; i++)
        {
            if (j == 0 && i == 0)
            {
                continue;
            }
            else if (j == 0)
            {
                int m = j * num_width_points + i;

                //structural
                int l = m - 1;
                Spring struct_spring = Spring(&point_masses[m], &point_masses[l], STRUCTURAL);
                springs.emplace_back(struct_spring);

                //bending
                int l2 = m - 2;
                if (i != 1 && i != 0)
                {
                    Spring bend_spring = Spring(&point_masses[m], &point_masses[l2], BENDING);
                    springs.emplace_back(bend_spring);
                }
            }
            else
            {
                int m = j * num_width_points + i;

                //structural
                int l = m - 1;
                int a = m - num_width_points;

                Spring abv_struct_spring = Spring(&point_masses[m], &point_masses[a], STRUCTURAL);
                springs.emplace_back(abv_struct_spring);

                if (i != 0)
                {
                    Spring l_struct_spring = Spring(&point_masses[m], &point_masses[l], STRUCTURAL);
                    springs.emplace_back(l_struct_spring);
                }

                //shearing
                int ul = a - 1;
                int ur = a + 1;

                if (i != 0)
                {
                    Spring ul_shear_spring = Spring(&point_masses[m], &point_masses[ul], SHEARING);
                    springs.emplace_back(ul_shear_spring);
                }

                if (i != num_width_points - 1)
                {
                    Spring ur_shear_spring = Spring(&point_masses[m], &point_masses[ur], SHEARING);
                    springs.emplace_back(ur_shear_spring);
                }

                //bending
                int l2 = m - 2;
                int a2 = a - num_width_points;

                if (j != 1)
                {
                    Spring a2_bend_spring = Spring(&point_masses[m], &point_masses[a2], BENDING);
                    springs.emplace_back(a2_bend_spring);
                }

                if (i != 0 && i != 1)
                {
                    Spring l2_bend_spring = Spring(&point_masses[m], &point_masses[l2], BENDING);
                    springs.emplace_back(l2_bend_spring);
                }
            }
        }
    }




}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters* cp,
    vector<Vector3D> external_accelerations,
    vector<CollisionObject*>* collision_objects) {
    double mass = width * height * cp->density / num_width_points / num_height_points;
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    // TODO (Part 2): Compute total force acting on each point mass.

    Vector3D force;
    Vector3D total_acceleration;

    for (Vector3D acceleration : external_accelerations)
    {
        force += mass * acceleration;
        total_acceleration += acceleration;
    }

    // Apply universal accelerations
    for (int i = 0; i < point_masses.size(); i++)
    {
        point_masses[i].forces = force;
    }

    // Apply Spring corrections
    for (int i = 0; i < springs.size(); i++)
    {
        double spring_force_magnitude = cp->ks * ((springs[i].pm_b->position - springs[i].pm_a->position).norm() - springs[i].rest_length);

        if (springs[i].spring_type == STRUCTURAL)
        {
            if (!cp->enable_structural_constraints)
            {
                continue;
            }

            Vector3D ab_force = spring_force_magnitude * (springs[i].pm_b->position - springs[i].pm_a->position).unit();
            Vector3D ba_force = spring_force_magnitude * (springs[i].pm_a->position - springs[i].pm_b->position).unit();
            springs[i].pm_a->forces += ab_force;
            springs[i].pm_b->forces += ba_force;
        }
        else if (springs[i].spring_type == SHEARING)
        {
            if (!cp->enable_shearing_constraints)
            {
                continue;
            }

            Vector3D ab_force = spring_force_magnitude * (springs[i].pm_b->position - springs[i].pm_a->position).unit();
            Vector3D ba_force = spring_force_magnitude * (springs[i].pm_a->position - springs[i].pm_b->position).unit();
            springs[i].pm_a->forces += ab_force;
            springs[i].pm_b->forces += ba_force;
        }
        else
        {
            if (!cp->enable_bending_constraints)
            {
                continue;
            }

            Vector3D ab_force = 0.2 * spring_force_magnitude * (springs[i].pm_b->position - springs[i].pm_a->position).unit();
            Vector3D ba_force = 0.2 * spring_force_magnitude * (springs[i].pm_a->position - springs[i].pm_b->position).unit();
            springs[i].pm_a->forces += ab_force;
            springs[i].pm_b->forces += ba_force;
        }
    }


    // TODO (Part 2): Use Verlet integration to compute new point mass positions
    for (int i = 0; i < point_masses.size(); i++)
    {
        if (point_masses[i].pinned)
        {
            continue;
        }

        Vector3D new_position = point_masses[i].position;
        new_position += (1 - (cp->damping / 100.0)) * (point_masses[i].position - point_masses[i].last_position);
        new_position += (point_masses[i].forces / mass * delta_t * delta_t);

        point_masses[i].last_position = point_masses[i].position;
        point_masses[i].position = new_position;
    }


    //// TODO (Part 4): Handle self-collisions.
    //build_spatial_map();
    //for (int i = 0; i < point_masses.size(); i++)
    //{
    //    self_collide(point_masses[i], simulation_steps);
    //}


    // TODO (Part 3): Handle collisions with other primitives.

    for (CollisionObject* obj : *collision_objects)
    {
        for (int i = 0; i < point_masses.size(); i++)
        {
            obj->collide(point_masses[i]);
        }
    }


    // TODO (Part 2): Constrain the changes to be such that the spring does not change
    // in length more than 10% per timestep [Provot 1995].

    for (int i = 0; i < springs.size(); i++)
    {
        if ((springs[i].pm_a->position - springs[i].pm_b->position).norm() < 1.10 * springs[i].rest_length)
        {
            continue;
        }

        double delta_len = (springs[i].pm_a->position - springs[i].pm_b->position).norm() - (1.10 * springs[i].rest_length);

        Vector3D b_a_unit = (springs[i].pm_a->position - springs[i].pm_b->position).unit();
        Vector3D a_b_unit = (springs[i].pm_b->position - springs[i].pm_a->position).unit();

        if (springs[i].pm_a->pinned)
        {
            springs[i].pm_b->position += delta_len * b_a_unit;
        }
        else if (springs[i].pm_b->pinned)
        {
            springs[i].pm_a->position += delta_len * a_b_unit;
        }
        else
        {
            springs[i].pm_b->position += (delta_len / 2.0) * b_a_unit;
            springs[i].pm_a->position += (delta_len / 2.0) * a_b_unit;
        }
    }

}

void Cloth::build_spatial_map() {
    for (const auto& entry : map) {
        delete(entry.second);
    }
    map.clear();

    // TODO (Part 4): Build a spatial map out of all of the point masses.
    for (int i = 0; i < point_masses.size(); i++)
    {
        float hash_key = hash_position(point_masses[i].position);
        vector<PointMass*>* point_mass_ptr = map[hash_key];

        if (point_mass_ptr != nullptr)
        {
            map[hash_key]->push_back(&point_masses[i]);
        }
        else
        {
            vector<PointMass*>* mass_vec = new vector<PointMass*>;
            mass_vec->push_back(&point_masses[i]);
            map[hash_key] = mass_vec;
        }
    }
}

void Cloth::self_collide(PointMass& pm, double simulation_steps) {
    // TODO (Part 4): Handle self-collision for a given point mass.
    float hash_key = hash_position(pm.position);
    vector<PointMass*>* point_masses_vec_ptr = map[hash_key];

    if (point_masses_vec_ptr == nullptr)
    {
        return;
    }
    vector<PointMass*> point_masses_vec = *point_masses_vec_ptr;

    Vector3D total_correction = Vector3D();
    int num_corrections = 0;

    for (auto pt_mass : point_masses_vec)
    {
        if (pt_mass == &pm)
        {
            continue;
        }

        double dist = (pm.position - pt_mass->position).norm();

        // Not clipping against each other so ignore these
        if (dist > 2 * thickness)
        {
            continue;
        }

        Vector3D correction = 2 * thickness * (pm.position - pt_mass->position).unit();
        total_correction += correction;
        num_corrections++;
    }

    if (num_corrections == 0)
    {
        return;
    }

    total_correction = total_correction / ((double)num_corrections);
    total_correction = total_correction / simulation_steps;
    pm.position += total_correction;
}

float CantorPairing(float x, float y)
{
    return 0.5 * (x + y) * (x + y + 1) + y;
}

float Cloth::hash_position(Vector3D pos) {
    // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.

    double t = max((3.0 * width / num_width_points), (3.0 * height / num_height_points));

    int box_x = floor(pos[0] / (3.0 * width / num_width_points));
    int box_y = floor(pos[1] / (3.0 * height / num_height_points));
    int box_z = floor(pos[2] / t);

    //float hash_key = CantorPairing(CantorPairing(box_x, box_y), box_z);
    float hash_key = box_x * 31 * 31 + box_y * 31 + box_z;
    return hash_key;
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
    PointMass* pm = &point_masses[0];
    for (int i = 0; i < point_masses.size(); i++) {
        pm->position = pm->start_position;
        pm->last_position = pm->start_position;
        pm++;
    }
}

void Cloth::buildClothMesh() {
    if (point_masses.size() == 0) return;

    ClothMesh* clothMesh = new ClothMesh();
    vector<Triangle*> triangles;

    // Create vector of triangles
    for (int y = 0; y < num_height_points - 1; y++) {
        for (int x = 0; x < num_width_points - 1; x++) {
            PointMass* pm = &point_masses[y * num_width_points + x];
            // Get neighboring point masses:
            /*                      *
             * pm_A -------- pm_B   *
             *             /        *
             *  |         /   |     *
             *  |        /    |     *
             *  |       /     |     *
             *  |      /      |     *
             *  |     /       |     *
             *  |    /        |     *
             *      /               *
             * pm_C -------- pm_D   *
             *                      *
             */

            float u_min = x;
            u_min /= num_width_points - 1;
            float u_max = x + 1;
            u_max /= num_width_points - 1;
            float v_min = y;
            v_min /= num_height_points - 1;
            float v_max = y + 1;
            v_max /= num_height_points - 1;

            PointMass* pm_A = pm;
            PointMass* pm_B = pm + 1;
            PointMass* pm_C = pm + num_width_points;
            PointMass* pm_D = pm + num_width_points + 1;

            Vector3D uv_A = Vector3D(u_min, v_min, 0);
            Vector3D uv_B = Vector3D(u_max, v_min, 0);
            Vector3D uv_C = Vector3D(u_min, v_max, 0);
            Vector3D uv_D = Vector3D(u_max, v_max, 0);


            // Both triangles defined by vertices in counter-clockwise orientation
            triangles.push_back(new Triangle(pm_A, pm_C, pm_B,
                uv_A, uv_C, uv_B));
            triangles.push_back(new Triangle(pm_B, pm_C, pm_D,
                uv_B, uv_C, uv_D));
        }
    }

    // For each triangle in row-order, create 3 edges and 3 internal halfedges
    for (int i = 0; i < triangles.size(); i++) {
        Triangle* t = triangles[i];

        // Allocate new halfedges on heap
        Halfedge* h1 = new Halfedge();
        Halfedge* h2 = new Halfedge();
        Halfedge* h3 = new Halfedge();

        // Allocate new edges on heap
        Edge* e1 = new Edge();
        Edge* e2 = new Edge();
        Edge* e3 = new Edge();

        // Assign a halfedge pointer to the triangle
        t->halfedge = h1;

        // Assign halfedge pointers to point masses
        t->pm1->halfedge = h1;
        t->pm2->halfedge = h2;
        t->pm3->halfedge = h3;

        // Update all halfedge pointers
        h1->edge = e1;
        h1->next = h2;
        h1->pm = t->pm1;
        h1->triangle = t;

        h2->edge = e2;
        h2->next = h3;
        h2->pm = t->pm2;
        h2->triangle = t;

        h3->edge = e3;
        h3->next = h1;
        h3->pm = t->pm3;
        h3->triangle = t;
    }

    // Go back through the cloth mesh and link triangles together using halfedge
    // twin pointers

    // Convenient variables for math
    int num_height_tris = (num_height_points - 1) * 2;
    int num_width_tris = (num_width_points - 1) * 2;

    bool topLeft = true;
    for (int i = 0; i < triangles.size(); i++) {
        Triangle* t = triangles[i];

        if (topLeft) {
            // Get left triangle, if it exists
            if (i % num_width_tris != 0) { // Not a left-most triangle
                Triangle* temp = triangles[i - 1];
                t->pm1->halfedge->twin = temp->pm3->halfedge;
            }
            else {
                t->pm1->halfedge->twin = nullptr;
            }

            // Get triangle above, if it exists
            if (i >= num_width_tris) { // Not a top-most triangle
                Triangle* temp = triangles[i - num_width_tris + 1];
                t->pm3->halfedge->twin = temp->pm2->halfedge;
            }
            else {
                t->pm3->halfedge->twin = nullptr;
            }

            // Get triangle to bottom right; guaranteed to exist
            Triangle* temp = triangles[i + 1];
            t->pm2->halfedge->twin = temp->pm1->halfedge;
        }
        else {
            // Get right triangle, if it exists
            if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
                Triangle* temp = triangles[i + 1];
                t->pm3->halfedge->twin = temp->pm1->halfedge;
            }
            else {
                t->pm3->halfedge->twin = nullptr;
            }

            // Get triangle below, if it exists
            if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
                Triangle* temp = triangles[i + num_width_tris - 1];
                t->pm2->halfedge->twin = temp->pm3->halfedge;
            }
            else {
                t->pm2->halfedge->twin = nullptr;
            }

            // Get triangle to top left; guaranteed to exist
            Triangle* temp = triangles[i - 1];
            t->pm1->halfedge->twin = temp->pm2->halfedge;
        }

        topLeft = !topLeft;
    }

    clothMesh->triangles = triangles;
    this->clothMesh = clothMesh;
}
