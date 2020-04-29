#pragma once
#include "types.h"

void build_tree(Simulation* sim);

void add_pt_to_node(BHTreeNode* node, Particle* pt);
void create_child(BHTreeNode* node, int oct);
int get_octant(BHTreeNode* node, Particle* pt);

void free_tree(Simulation* sim);
void free_node(BHTreeNode* node);
