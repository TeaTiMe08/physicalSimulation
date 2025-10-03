# physicalSimulation — MPI
This submodule implements a parallel N-body / orbital simulation using **MPI** (Message Passing Interface) in C.  
It allows distributed computation of gravitational interactions (or similar physics) across multiple processes.

#### Background
This is a project i had done for my Computer Science degree in a course about Parallel computing.

<img width="30%" alt="exampleImageRendered" src="https://github.com/user-attachments/assets/994e6bed-9c48-40f3-a2b9-9243964d0049" />


## Overview

The MPI module is designed to parallelize the computationally expensive parts of physics simulation — e.g. computing forces in an N-body system or evolving orbits — by distributing work across multiple processes.  

Each process computes partial results, communicates boundary or shared data, and then aggregates to maintain global state. This approach enables scaling to more bodies or higher resolution than a sequential version.

---

## Features

- Parallel force / acceleration computation using MPI  
- Domain decomposition or workload partitioning across ranks  
- Aggregation and synchronization via MPI collective calls  
- Modular design: existing sequential code (if any) can be adapted  
- Basic timing / performance measurement hooks (if included)  

---

## Prerequisites

- A C compiler supporting MPI (e.g. `mpicc`, or `gcc` with MPI wrappers)  
- MPI implementation installed (e.g. Open MPI, MPICH)  
- Make / build tools (make, cmake, etc. depending on setup)  
- (Optional) Ability to run across multiple nodes or cores  

---

## Building & Running

Here’s a generic workflow. Adjust according to your build system (Makefile, CMake, etc.):

```bash
cd mpi
make             # or your build command
```
To run with, say, 4 processes:
```
mpirun -np 4 ./your_simulation_binary [arguments...]
```

## Configuration / Parameters

Your simulation executable likely accepts arguments such as:

| Parameter         | Meaning                                      | Example             |
|-------------------|----------------------------------------------|---------------------|
| number_of_bodies  | Total count of simulated particles           | `1000`              |
| time_steps        | Number of discrete steps                      | `10000`             |
| output_interval   | How often to record / dump state              | every `100` steps   |
| initial_positions | (Optional) file or mode for initial conditions| random, file path   |
| seed              | Random seed for reproducibility               | `42`                |

## Example Usage
```
mpirun -np 8 ./nbody_mpi --bodies 5000 --steps 20000 --out-interval 500
```
This runs the simulation with 5000 bodies, 20k steps, and outputs state every 500 steps using 8 MPI ranks.

🍵
