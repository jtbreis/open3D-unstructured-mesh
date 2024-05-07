#include "open3d/t/geometry/TriangleMesh.h"
#include "open3d/io/TriangleMeshIO.h"
#include "open3d/t/geometry/PointCloud.h"
#include <iterator>
#include <open3d/geometry/PointCloud.h>
#include <open3d/geometry/TriangleMesh.h>
#include<string>
#include<iostream>
#include<cmath>

open3d::geometry::TriangleMesh create_plane(double radius, int nelements, double offset) {
    double spacing = 2*M_PI*radius / nelements;
    int number_of_points = 2*radius / spacing;
   
    std::vector<Eigen::Vector3d> points;
    std::vector<long unsigned int> removeidx;
    for (int i = 0; i <= number_of_points; i++) {
        for (int j = 0; j <= number_of_points; j++) {
            double x = -radius + i * spacing;
            double y = -radius + j * spacing;
            Eigen::Vector3d pos = {x,y,offset};
            points.push_back(pos);
            if (pos.squaredNorm()-offset > radius-spacing/2) {
                removeidx.push_back(points.size()-1);
            }
        }
    }
    // Define faces for the plane
    std::vector<Eigen::Vector3i> faces;
    for (int i = 0; i < number_of_points; i++) {
        for (int j = 0; j < number_of_points; j++) {
            int index0 = j + (number_of_points+1) * i;
            int index1 = index0 + 1;
            int index2 = index0 + number_of_points + 1;
            int index3 = index2 + 1;
            faces.push_back(Eigen::Vector3i(index0, index1, index2));
            faces.push_back(Eigen::Vector3i(index2, index3, index1));
        }
    }
    open3d::geometry::TriangleMesh mesh;
    mesh.vertices_ = points;
    mesh.triangles_ = faces;
    mesh.RemoveVerticesByIndex(removeidx);
    return mesh;
}

std::vector<long unsigned int> get_outside_vertices(open3d::geometry::TriangleMesh& mesh) {
    std::vector<int> vertex_counts(mesh.vertices_.size(), 0);
    for (const auto& triangle : mesh.triangles_) {
        for (int i = 0; i < 3; ++i) {
            vertex_counts[triangle[i]]++;
        }
    }

    std::vector<long unsigned int> outside_vertices;
    for (size_t i = 0; i < vertex_counts.size(); ++i) {
        if (vertex_counts[i]<6) {
            outside_vertices.push_back(i);
        }
    }
    return outside_vertices;
}

std::vector<long unsigned int> get_closest_neighbors(std::vector<Eigen::Vector3d> inside_points, Eigen::Vector3d point) {
    std::vector<long unsigned int> closest_neighbors(2);
    std::vector<double> closest_distances = {5, 5};
    int idx = 0;
    for (Eigen::Vector3d& curr_point : inside_points) {
        Eigen::Vector3d distance = curr_point - point; 
        double length = distance.squaredNorm();
        if (length > closest_distances[1]) {
            ++idx;
            continue;
        }
        if (length < closest_distances[0]) {
            closest_neighbors[1] = closest_neighbors[0];
            closest_distances[1] = closest_distances[0];
            closest_neighbors[0] = idx;
            closest_distances[0] = length;
        } else {
            closest_neighbors[1] = idx;
            closest_distances[1] = length;
        }
        ++idx;
    }
    return closest_neighbors;
}

template <typename T>
std::vector<T> extractElements(const std::vector<T>& largerVector, const std::vector<size_t>& indices) {
    std::vector<T> extractedElements;
    for (size_t index : indices) {
        // Check if the index is within the bounds of the larger vector
        if (index < largerVector.size()) {
            // Add the element at the specified index to the extracted elements
            extractedElements.push_back(largerVector[index]);
        }
    }
    return extractedElements;
}

std::vector<Eigen::Vector3d> create_circle_positions(double radius, double offset, int number_of_elements) {
    double arc_spacing = 2*M_PI / (number_of_elements);
    std::vector<Eigen::Vector3d> positions;
    for (double angle = 0; angle < 2*M_PI; angle+=arc_spacing) {
        double x = std::cos(angle) * radius;
        double y = std::sin(angle) * radius;
        double z = offset;
        positions.push_back(Eigen::Vector3d(x, y, z));
    }
    return positions;
}

void make_circle(open3d::geometry::TriangleMesh& mesh, std::vector<long unsigned int> outside_vertices, 
                                           std::vector<Eigen::Vector3d> outside_positions, int number_of_elements) {
    std::vector<Eigen::Vector3d> circle_positions = create_circle_positions(1.0, 1.0, number_of_elements);
    std::vector<Eigen::Vector3i> faces;
    std::vector<int> connections(circle_positions.size());
    int curr_idx;
    int idx_offset = mesh.vertices_.size();
    long unsigned int i = 0;
    for (Eigen::Vector3d& pos : circle_positions) {
        std::vector<long unsigned int> neighbor_idx = get_closest_neighbors(outside_positions, pos);
        std::vector<long unsigned int> outside_idx = extractElements(outside_vertices, neighbor_idx);
        curr_idx = i + idx_offset;
        long unsigned int index_1 = (i == number_of_elements-1) ? idx_offset : curr_idx+1;
        if (i == 0) {
            connections[i]=outside_idx[0];
            faces.push_back(Eigen::Vector3i(curr_idx, outside_idx[0], outside_idx[1]));
            faces.push_back(Eigen::Vector3i(curr_idx, index_1, outside_idx[0]));
            connections[index_1-idx_offset]=outside_idx[0];
            i++;
            continue;
        }
        if (connections[i]==outside_idx[0]) {
            faces.push_back(Eigen::Vector3i(curr_idx, index_1, outside_idx[0]));
            connections[i]=outside_idx[0];
            connections[index_1-idx_offset]=outside_idx[0];
        } else {
            long unsigned int idx = (i == number_of_elements) ? connections[0] : connections[i];
            faces.push_back(Eigen::Vector3i(curr_idx, outside_idx[0], idx));
            faces.push_back(Eigen::Vector3i(curr_idx, index_1, outside_idx[0]));
            connections[i]=outside_idx[0];
            connections[index_1-idx_offset]=outside_idx[0];
        }
        i++;
    }

    mesh.vertices_.insert(mesh.vertices_.end(), circle_positions.begin(), circle_positions.end());
    mesh.triangles_.insert(mesh.triangles_.end(), faces.begin(), faces.end());
}

void repositionVertices(open3d::geometry::TriangleMesh& mesh, std::vector<long unsigned int> outside_vertices, double radius, double offset, int nelements) {
    double spacing = 2*M_PI*radius / nelements;
    int number_of_outside_vertices = outside_vertices.size();
    double arc_spacing = 2 * M_PI / number_of_outside_vertices;
    double angle = 0;
    long unsigned int last_repositioned = 0;
    std::vector<Eigen::Vector3d> circle_positions = create_circle_positions(radius-spacing/2, offset, 4*nelements);
    circle_positions.erase(circle_positions.end()-10, circle_positions.end()-4);
    for (long unsigned int& idx : outside_vertices) {
        std::vector<long unsigned int> closest_neighbors = get_closest_neighbors(circle_positions, mesh.vertices_[idx]);
        mesh.vertices_[idx] = circle_positions[closest_neighbors[0]];
    }
}

open3d::geometry::TriangleMesh create_sides(double radius, int number_of_elements, double offset) {
    open3d::geometry::TriangleMesh plane = create_plane(radius, number_of_elements, offset);

    std::vector<long unsigned int> outside_vertices = get_outside_vertices(plane);
    repositionVertices(plane, outside_vertices, radius, offset, number_of_elements);
    std::vector<Eigen::Vector3d> pos_outside_vertices = extractElements(plane.vertices_, outside_vertices);

    make_circle(plane, outside_vertices, pos_outside_vertices, number_of_elements);
    return plane;
}


int main() {
    //auto hull = open3d::geometry::TriangleMesh::CreateCylinder(1.0,2.0,20,1);
    //std::vector<long unsigned int> removeidx= {0,1};
    //hull->RemoveVerticesByIndex(removeidx);
    //open3d::geometry::PointCloud hull_points;
    //hull_points.points_ = hull->vertices_;
    //hull_points.normals_ = hull->vertex_normals_;
    //open3d::geometry::TriangleMesh bottom_points = create_plane(1.0, 0.1, 1.0);
    //open3d::geometry::PointCloud cylinder;
    //cylinder += hull_points;
    //cylinder += top_points;
    //cylinder += bottom_points;

//    for (int i = 8; i <= 26; i++) {
    int i = 42;
        std::string filename = "test" + std::to_string(i) + ".obj";
        open3d::geometry::TriangleMesh mesh = create_sides(1.0, i, 1.0);
        open3d::io::WriteTriangleMeshToOBJ(filename, mesh, false, false, true, false, false, false);
//    }

    return 0;
}
