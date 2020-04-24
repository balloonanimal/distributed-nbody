#pragma once
#include "types.h"

void build_tree(Simulation* sim);

void add_pt_to_node(BHTree* node, Particle* pt);
void create_child(BHTreeNode* node, int oct);
int get_octant(BHTreeNode* node, Particle* pt);
