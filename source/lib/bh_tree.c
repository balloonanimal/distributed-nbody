#include "bh_tree.h"
#include <stdlib.h>

#define MIN_WIDTH 0.01

void build_tree(Simulation *sim) {
  // alloc tree
  sim->tree = calloc(1, sizeof(BHTreeNode));
  sim->tree->width = sim->width;

  for (int i = 0; i < sim->parray.len; i++) {
    add_pt_to_node(sim->tree, &sim->parray.particles[i]);
  }
}

void add_pt_to_node(BHTreeNode *node, Particle *pt) {
  // Node is terminal
  if (node->pcount == 0) {
    node->pt = pt;
    node->COM_mass = pt->mass;
    node->COM_x = pt->pos_x;
    node->COM_y = pt->pos_y;
    node->COM_z = pt->pos_z;
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
    int oct_ = get_octant(node, pt);
    // NOTE: offset the octant of a new particle if the tree partition is
    // getting very small. This is to stop the tree from overflowing if two
    // particles with practically identical positions are placed in the same
    // tree.
    // This can cause physical innacuracies, but since it only happens for
    // extrememly close together points it shouldn't be an issue
    if (oct == oct_ && node->width < MIN_WIDTH) {
      oct_ = (oct_ + 1) % 8;
    }
    if (node->children[oct_] == NULL) {
      create_child(node, oct_);
    }
    add_pt_to_node(node->children[oct_], pt);
    // recalc COM
    node->COM_mass = node->COM_mass + pt->mass;
    node->COM_x = (node->COM_x + pt->pos_x) / 2;
    node->COM_y = (node->COM_y + pt->pos_y) / 2;
    node->COM_z = (node->COM_z + pt->pos_z) / 2;
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
    node->COM_mass = node->COM_mass + pt->mass;
    node->COM_x =
        ((node->pcount * node->COM_x + pt->pos_x) / (node->pcount + 1));
    node->COM_y =
        ((node->pcount * node->COM_y + pt->pos_y) / (node->pcount + 1));
    node->COM_z =
        ((node->pcount * node->COM_z + pt->pos_z) / (node->pcount + 1));
    node->pcount += 1;
  }
}

void create_child(BHTreeNode *node, int oct) {
  BHTreeNode *child = calloc(1, sizeof(BHTreeNode));
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
  child->x = node->x + xdir * child->width;
  child->y = node->y + ydir * child->width;
  child->z = node->z + zdir * child->width;
  node->children[oct] = child;
}

int get_octant(BHTreeNode *node, Particle *pt) {
  int oct = 0;
  if (pt->pos_x < node->x) {
    oct += 4;
  }
  if (pt->pos_y < node->y) {
    oct += 2;
  }
  if (pt->pos_z < node->z) {
    oct += 1;
  }
  return oct;
}

void free_tree(Simulation *sim) {
  if (sim->tree == NULL) {
    return;
  }
  free_node(sim->tree);
  sim->tree = NULL;
}

void free_node(BHTreeNode *node) {
  for (int oct = 0; oct < 8; oct++) {
    if (node->children[oct] != NULL) {
      free_node(node->children[oct]);
    }
  }
  free(node);
}
