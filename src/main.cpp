#include "open3d/t/geometry/TriangleMesh.h"
#include "open3d/io/TriangleMeshIO.h"
#include "open3d/t/geometry/PointCloud.h"
#include <iterator>
#include <limits>
#include <open3d/geometry/PointCloud.h>
#include <open3d/geometry/TriangleMesh.h>
#include<string>
#include<iostream>
#include<cmath>

void create_plane(open3d::geometry::TriangleMesh& mesh, const double& radius, const int& nelements, const double& offset) {
    const double spacing = 2*M_PI*radius / (nelements*1.2);
    const int number_of_points = static_cast<int>(2*radius / spacing) + 2;
    const double cutoff = radius + spacing;
   
    std::vector<Eigen::Vector3d> points;
    std::vector<long unsigned int> removeidx;
    for (int i = 0; i <= number_of_points; ++i) {
        for (int j = 0; j <= number_of_points; ++j) {
            double x = -radius + i * spacing - spacing;
            double y = -radius + j * spacing - spacing;
            Eigen::Vector3d pos = {x,y,offset};
            points.push_back(pos);
            if (pos.squaredNorm()-offset > cutoff) {
                removeidx.push_back(points.size()-1);
            }
        }
    }
    // Define faces for the plane
    std::vector<Eigen::Vector3i> faces;
    for (int i = 0; i < number_of_points; ++i) {
        for (int j = 0; j < number_of_points; ++j) {
            int index0 = j + (number_of_points+1) * i;
            int index1 = index0 + 1;
            int index2 = index0 + number_of_points + 1;
            int index3 = index2 + 1;
            if ((i <= number_of_points/2 && j <= number_of_points/2) || (i > number_of_points/2 && j > number_of_points/2)) {
                faces.push_back(Eigen::Vector3i(index0, index2, index3));
                faces.push_back(Eigen::Vector3i(index0, index1, index3));
            } else {
                faces.push_back(Eigen::Vector3i(index0, index1, index2));
                faces.push_back(Eigen::Vector3i(index1, index2, index3));
            }
        }
    }
    mesh.vertices_ = points;
    mesh.triangles_ = faces;
    mesh.RemoveVerticesByIndex(removeidx);
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

std::vector<long unsigned int> get_closest_neighbors(std::vector<Eigen::Vector3d>& inside_points, const Eigen::Vector3d& point) {
    std::vector<long unsigned int> closest_neighbors(3,0);
    std::vector<double> closest_distances(3,std::numeric_limits<double>::max());;
    int idx = 0;
    for (Eigen::Vector3d& curr_point : inside_points) {
        Eigen::Vector3d distance = curr_point - point; 
        double length = distance.squaredNorm();
        if (length > closest_distances[2]) {
            ++idx;
            continue;
        }
        if (length < closest_distances[0]) {
            closest_neighbors[1] = closest_neighbors[0];
            closest_distances[1] = closest_distances[0];
            closest_neighbors[0] = idx;
            closest_distances[0] = length;
        } else if (length < closest_distances[1]) {
            closest_neighbors[2] = closest_neighbors[1];
            closest_distances[2] = closest_distances[1];
            closest_neighbors[1] = idx;
            closest_distances[1] = length;
        } else {
            closest_neighbors[2] = idx;
            closest_distances[2] = length;
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

std::vector<Eigen::Vector3d> create_circle_positions(const double& radius, const double& offset, const int& number_of_elements) {
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

std::vector<long unsigned int> find_triangle(open3d::geometry::TriangleMesh& mesh, std::vector<long unsigned int>& vertices) {
    unsigned long int idx = 0;
    std::vector<long unsigned int> triangle_idx; 
    for (auto& triangle : mesh.triangles_) {
        if ((vertices[0] == triangle[0] && vertices[1] == triangle[1]) ||
            (vertices[0] == triangle[1] && vertices[1] == triangle[0]) ||
            (vertices[0] == triangle[1] && vertices[1] == triangle[2]) ||
            (vertices[0] == triangle[2] && vertices[1] == triangle[1]) ||
            (vertices[0] == triangle[2] && vertices[1] == triangle[0]) ||
            (vertices[0] == triangle[0] && vertices[1] == triangle[2])) {
            triangle_idx.push_back(idx);
            return triangle_idx;
        }
        idx++;
    }
    return triangle_idx;
}

long unsigned int findMissingValue(const Eigen::Vector3i& values, std::vector<long unsigned int>& points) {
    for (long unsigned int value : values.array()) {
        if (value != points[0] && value != points[1]) {
            return value;
        }
    }
    std::cout << "Value not found\n";
    return -1; // Value not found
}

void addMissingConnections(open3d::geometry::TriangleMesh& mesh, std::vector<Eigen::Vector3d>& circle_positions, std::vector<long unsigned int>& connections) {
    for (long unsigned int idx = 0; idx < connections.size(); ++idx) {
        if (connections[idx] != 0) {
            continue;
        }
        std::vector<long unsigned int> closest_neighbors = get_closest_neighbors(mesh.vertices_, circle_positions[idx]);
        std::vector<long unsigned int> triangle_idx = find_triangle(mesh, closest_neighbors);
        Eigen::Vector3i triangle = mesh.triangles_[triangle_idx[0]];
        long unsigned int tip_index = findMissingValue(triangle, closest_neighbors);

        if (!triangle_idx.empty()) {
            mesh.RemoveTrianglesByIndex(triangle_idx);
        }
        mesh.vertices_.push_back(circle_positions[idx]);
        mesh.triangles_.push_back(Eigen::Vector3i(closest_neighbors[0], mesh.vertices_.size()-1, tip_index));
        mesh.triangles_.push_back(Eigen::Vector3i(closest_neighbors[1], mesh.vertices_.size()-1, tip_index));

        triangle_idx.clear();
    }
}

void repositionVertices(open3d::geometry::TriangleMesh& mesh, std::vector<long unsigned int>& outside_vertices, const double& radius, const double& offset, const int& nelements) {
    std::vector<Eigen::Vector3d> circle_positions = create_circle_positions(radius, offset, nelements);
    std::vector<long unsigned int> connections(circle_positions.size(), 0);
    for (long unsigned int& idx : outside_vertices) {
        std::vector<long unsigned int> closest_neighbors = get_closest_neighbors(circle_positions, mesh.vertices_[idx]);
        mesh.vertices_[idx] = circle_positions[closest_neighbors[0]];
        connections[closest_neighbors[0]]++;
    }
    addMissingConnections(mesh, circle_positions, connections);
}

void create_sides(open3d::geometry::TriangleMesh& mesh, const double& radius, const int& number_of_elements, const double& offset) {
    create_plane(mesh, radius, number_of_elements, offset);
    std::vector<long unsigned int> outside_vertices = get_outside_vertices(mesh);
    repositionVertices(mesh, outside_vertices, radius, offset, number_of_elements);
}


int main() {
    for (int i = 30; i <= 50; i++) {
        open3d::geometry::TriangleMesh* mesh = new open3d::geometry::TriangleMesh();
        std::string filename = "test" + std::to_string(i) + ".obj";
        create_sides(*mesh, 1.0, i, 1.0);
        open3d::io::WriteTriangleMeshToOBJ(filename,*mesh, false, false, true, false, false, false);
    }
    return 0;
}
