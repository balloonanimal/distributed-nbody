#include "bh_tree.h"
#include <stdlib.h>

void build_tree(Simulation *sim) {
    // alloc tree
    sim->tree = calloc(1, sizeof(BHTreeNode));

    for (int i = 1; i < sim->parray.len; i++) {
        add_pt_to_node(sim->tree, &sim->parray.particles[i]);
    }
}

void add_pt_to_node(BHTreeNode* node, Particle* pt) {
    // Node is terminal
    if (node->pcount == 0) {
        node->pt = pt;
        node->COM_mass = pt->mass;
        node->COM_position = pt->position;
        node->pcount = 1;
    }
    // Node is terminal but non-empty
    else if (node->pcount == 1) {
        // relocate node->pt
        int oct = get_octant(node, node->pt);
        if (node->children[oct] == NULL) {
            create_child(node, oct);
        }
        add_pt_to_node(node->children[oct], node->pt);
        node->pt = NULL;
        // relocate pt
        oct = get_octant(node, pt);
        if (node->children[oct] == NULL) {
            create_child(node, oct);
        }
        add_pt_to_node(node->children[oct], pt);
        // recalc COM
        node->COM_mass = (node->COM_mass + pt->mass) / 2;
        Vec3d COM_position = {
            (node->position.x + pt->position.x) / 2,
            (node->position.y + pt->position.y) / 2,
            (node->position.z + pt->position.z) / 2};
        node->COM_position = COM_position;
        node->pcount = 2;
    }
    // Node is not terminal
    else {
        int oct = get_octant(node, pt);
        if (node->children[oct] == NULL) {
            create_child(node, oct);
        }
        add_pt_to_node(node->children[oct], pt);
        // recalc COM
        node->COM_mass = (node->pcount * node->COM_mass
                          + pt->mass) / (node->pcount + 1);
        Vec3d COM_position = {
            (node->pcount * node->position.x + pt->position.x) / (node->pcount + 1),
            (node->pcount * node->position.y + pt->position.y) / (node->pcount + 1),
            (node->pcount * node->position.z + pt->position.z) / (node->pcount + 1)};
        node->COM_position = COM_position;
        node->pcount += 1;
    }
}

void create_child(BHTreeNode* node, int oct) {
    BHTreeNode* child = calloc(1, sizeof(BHTreeNode));
    child->width = node->width / 2;
    // position calculations
    int xdir, ydir, zdir;
    if (oct == 0 || oct == 1 || oct == 2 || oct == 3) {
        xdir = 1;
    } else {
        xdir = -1;
    }
    if (oct == 0 || oct == 1 || oct == 4 || oct == 5) {
        ydir = 1;
    } else {
        ydir = -1;
    }
    if (oct == 0 || oct == 2 || oct == 4 || oct == 6) {
        zdir = 1;
    } else {
        zdir = -1;
    }
    node->position.x = node->position.x + xdir * child->width;
    node->position.y = node->position.y + ydir * child->width;
    node->position.z = node->position.z + zdir * child->width;
    node->children[oct] = child;
}

int get_octant(BHTreeNode* node, Particle* pt) {
    int oct = 0;
    if (pt->position.x < node->position.x) {
        oct += 4;
    }
    if (pt->position.y < node->position.y) {
        oct += 2;
    }
    if (pt->position.z < node->position.z) {
        oct += 1;
    }
    return oct;
}

void free_tree(Simulation* sim) {
    if (sim->tree == NULL) {
        return;
    }
    free_node(sim->tree);
    sim->tree = NULL;
}

void free_node(BHTreeNode* node) {
    for (int oct = 0; oct < 8; oct++) {
        if (node->children[oct] != NULL) {
            free_node(node->children[oct]);
        }
    }
    free(node);
}
