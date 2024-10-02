# Chasqui: Allosteric Pathways Between Protein Residues

Finds an allosteric pathways between 2 protein residues based on an interaction map (for instance, a Cieplak
interaction map), and the [Gonum](https://www.gonum.org) implementation of the Dijkstra's algorithm. Connections
between residues are weighted by the energy of the corresponding interaction, and those weights are used to choose
the "best" pathway between two residues. 

Additional pathways can be found by imposing penalties (by default, random penalties) to some of the
interactions in order to force the algorithm to take different routes.

## Installation

Chasqui is a pure Go program and supports Modules. Thus, after cloning the repository, and assuming the [Go](https://go.dev/) toolchain is installed, issuing the command:

```
go mod tidy
```

Followed by:

```
go build
```

Should produce a functioning binary. 

## Usage

```
chasqui [flags] geomtry.pdb ResidueID1 ResidueID2 Chain1 Chain2 
```

Several flags are available. Use the following for details:

```
chasqui -h
```


LICENSE

Copyright (c) 2024 Raul Mera <rmeraa{at}academicosDOTutaDOTcl>


This program, including its documentation,
is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2.0 as
published by the Free Software Foundation.

This program and its documentation is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General
Public License along with this program.  If not, see
<http://www.gnu.org/licenses/>.

