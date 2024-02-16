/*
   This example is free and distributed under
   Mozilla Public License Version 2.0.
*/

// clear && astyle so-rastrigin.cpp && g++ so-rastrigin.cpp -o so-rastrigin.exe && so-rastrigin.exe

#include <cmath>
#include <string>
#include <vector>
#include <climits>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "openGA.hpp"

const int SPACE_SIZE = 100;
const int POPULATION_SIZE = 100;
const int MAXIMUM_GENERATIONS = 1000;

struct MySolution {
    std::vector < double > x;

    std::string to_string() const {
        std::ostringstream out;

        out << "{";
        for (unsigned long i = 0; i < x.size(); i++) {
            out << (i ? "," : "") << std::setprecision(10) << x[i];
        }
        out << "}";

        return out.str();
    }
};

struct MyMiddleCost {
    // This is where the results of simulation is stored but not yet finalized.
    double cost;
};

typedef EA::Genetic < MySolution, MyMiddleCost > GA_Type;
typedef EA::GenerationType < MySolution, MyMiddleCost > Generation_Type;

void init_genes(MySolution & p, const std:: function < double(void) > & rnd01) {
    // Straight twin.
    for (int i = 0; i < SPACE_SIZE; i++) {
        p.x.push_back(2.0 * (rnd01() - 0.5));
    }
    // Inverted twin.
    for (int i = 0; i < SPACE_SIZE; i++) {
        p.x.push_back( -p.x[i] );
    }
}

bool eval_solution(const MySolution & p, MyMiddleCost & c) {
    double cost1, cost2;

    // DNA-insired pairs.
    cost1 = cost2 = 10.0 * p.x.size() / 2;
    for (int i = 0; i < p.x.size()/2; i++) {
        cost1 += p.x[i] * p.x[i] - 10.0 * std::cos(M_2_PI * p.x[i]);
    }
    for (int j = p.x.size()/2; j < p.x.size(); j++) {
        cost2 += p.x[j] * p.x[j] - 10.0 * std::cos(M_2_PI * p.x[j]);
    }

    c.cost = std::min(cost1, cost2);

    return true;
}

MySolution mutate(
    const MySolution & base, const std:: function < double(void) > & rnd01, double shrink_scale) {
    MySolution next = base;

    //Toggle bit and random postiton.
    int index = rand() % (next.x.size() / 2);
    int parts = sizeof(base.x[index]) / sizeof(unsigned short);
    unsigned short *mutated = (unsigned short *)(&next.x[ index ]);
    mutated += rand() % parts;
    unsigned short mask = 1;
    mask <<= rand() % (8 * sizeof(mask));
    (*mutated) ^= mask;

    //Handle the state of the inverse twin.
    next.x[ next.x.size() / 2 + index ] = -next.x[ index ];

    return next;
}

MySolution crossover(
    const MySolution & A, const MySolution & B, const std:: function < double(void) > & rnd01) {
    MySolution next;

    int size = std::min(A.x.size() / 2, B.x.size() / 2);
    for (int i = 0; i < size; i++) {
        next.x.push_back((rand()%2) ? A.x[i] : B.x[i]);
    }
    for (int i = 0; i < size; i++) {
        next.x.push_back( -next.x[i] );
    }

    return next;
}

double calculate_SO_total_fitness(const GA_Type::thisChromosomeType & X) {
    // Finalize the cost.
    return X.middle_costs.cost;
}

std::ofstream output_file;

void SO_report_generation(
    int generation_number,
    const EA::GenerationType < MySolution, MyMiddleCost > & last_generation, const MySolution & best_genes) {
    std::cout <<
              "Generation [" << generation_number << "], " <<
              "Best=" << last_generation.best_total_cost << ", " <<
              "Average=" << last_generation.average_cost << ", " <<
              "Best genes=(" << best_genes.to_string() << ")" << ", " <<
              "Exe_time=" << last_generation.exe_time <<
              std::endl;

    output_file <<
                generation_number << "\t" <<
                last_generation.average_cost << "\t" <<
                last_generation.best_total_cost << "\t";

    for(int i=0; i<best_genes.x.size(); i++) {
        output_file << best_genes.x[i] << "\t";
    }

    output_file << std::endl;
}

int main() {
    srand(time(NULL));

    output_file.open("result-so-rastrigin.txt");

    output_file <<
                "step" << "\t" <<
                "cost_avg" << "\t" <<
                "cost_best" << "\t";

    for(int i=0; i<2*SPACE_SIZE; i++) {
        output_file << "best(" << (i+1) << ")" << "\t";
    }

    output_file << std::endl;

    EA::Chronometer timer;
    timer.tic();

    GA_Type ga_obj;

    ga_obj.problem_mode = EA::GA_MODE::SOGA;
    ga_obj.multi_threading = false;
    ga_obj.dynamic_threading = false;
    ga_obj.idle_delay_us = 0; // Switch between threads quickly.
    ga_obj.verbose = false;
    ga_obj.population = POPULATION_SIZE;
    ga_obj.generation_max = MAXIMUM_GENERATIONS;
    ga_obj.calculate_SO_total_fitness = calculate_SO_total_fitness;
    ga_obj.init_genes = init_genes;
    ga_obj.eval_solution = eval_solution;
    ga_obj.mutate = mutate;
    ga_obj.crossover = crossover;
    ga_obj.SO_report_generation = SO_report_generation;
    ga_obj.best_stall_max = 20;
    ga_obj.average_stall_max = 20;
    ga_obj.tol_stall_best = 1e-6;
    ga_obj.tol_stall_average = 1e-6;
    ga_obj.elite_count = 10;
    ga_obj.crossover_fraction = 0.7;
    ga_obj.mutation_rate = 0.1;

    ga_obj.solve();

    std::cout << "The problem is optimized in " << timer.toc() << " seconds." << std::endl;

    output_file.close();

    return EXIT_SUCCESS;
}
