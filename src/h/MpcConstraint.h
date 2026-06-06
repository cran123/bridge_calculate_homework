#pragma once

#include <vector>

struct CMpcTerm
{
	unsigned int node;
	unsigned int dof;
	double coefficient;
};

struct CMpcConstraint
{
	std::vector<CMpcTerm> terms;
	std::vector<unsigned int> equations;
};
