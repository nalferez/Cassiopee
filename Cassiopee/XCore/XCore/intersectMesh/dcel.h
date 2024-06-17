#pragma once

#include <vector>
#include <cstddef>

#include "queue.h"

struct Vertex;
struct Hedge;
struct Face;
struct Segment;
struct Smesh;
struct Cycle;

struct Dcel {
    std::vector<Vertex *> V;
    std::vector<Hedge *> H;
    std::vector<Face *> F;
    std::vector<Cycle *> C;

    Queue Q; // Filters out duplicate vertices

    Face *f_unbounded[2];

    static E_Int RED;
    static E_Int BLACK;
    static E_Int NO_IDEA;

    Dcel(const Smesh &M0, const Smesh &M1);
    ~Dcel();
    
    void init_vertices(const Smesh &M0, const Smesh &M1);

    void init_hedges_and_faces(const Smesh &M, E_Int color);

    static E_Int check_hedges(const std::vector<Hedge *> &H);

    static E_Int check_faces(const std::vector<Hedge *> &H,
        const std::vector<Face *> &F);

    void find_intersections();

    static void resolve(Vertex *p, const std::vector<Segment *> &L,
        const std::vector<Segment *> &C, const std::vector<Segment *> &U,
        std::vector<Hedge *> &H);
    
    void make_cycles();

    void set_face_labels(std::vector<Face *> &F);

    Hedge *get_hedge_of_color(Face *f, E_Int color);

    std::vector<Face *> make_cycle_faces(const std::vector<Cycle *> &C);

    void update_hedge_faces(const std::vector<Face *> &F);

    void set_cycles_inout();

    std::vector<E_Int> extract_indices_of_type(E_Int inout);
    
    std::vector<Face *> extract_faces_of_indices(
        const std::vector<E_Int> &indices);

    void write_ngon(const char *fname, const std::vector<Face *> &faces) const;

    void write_degen_faces(const char *fname);
    
    void write_outer_faces(const char *fname);
    
    void write_inner_faces(const char *fname);

    static std::vector<Vertex *> get_face_vertices(Face *f);
};