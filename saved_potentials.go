/*
saved_potentials.go, part of Chasqui



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

*/

package main

import (
	"fmt"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/scu"
)

func CieplakIsSymmetric(j int, R []*chem.Bond) bool {
	for i, v := range R {
		if i == j {
			continue
		}

		if R[j].At1.ID == v.At2.ID && R[j].At2.ID == v.At1.ID {
			return true
		}
	}
	return false
}

func CieplakIsThere(r *chem.Bond, R []*chem.Bond) bool {
	for _, v := range R {
		if r.At1.ID == v.At2.ID && r.At2.ID == v.At1.ID {
			return true
		}
	}
	return false
}

func IsContiguous(MolID1, MolID2, limit int, chain1, chain2 string) bool {
	dif := MolID1 - MolID2
	if dif < 0 {
		dif *= -1
	}
	if chain1 == chain2 && dif < limit {
		return true
	}
	return false
}

func CieplakPotRead(mol chem.Atomer, mapname string, exclude int, e float64, lastID int) []*chem.Bond {
	if e <= 0 {
		e = 12.0 * chem.KJ2Kcal
	}
	inp, err := scu.NewMustReadFile(mapname)
	scu.QErr(err)
	defer inp.Close()
	set := make([]*chem.Bond, 0, 150)
	reading := false
	for i := inp.Next(); i != "EOF"; i = inp.Next() {
		if strings.Contains(i, "    I1  AA  C I(PDB)    I2  AA  C I(PDB)    DISTANCE       CMs    rCSU    aSurf    rSurf    nSurf") {
			inp.Next() //skip a line. The file should definitely not end here, so, if this is EOF, we'll let it panic on the next round
			reading = true
			continue
		}
		if reading && i == "\n" {
			break
		}
		if !reading {
			continue
		}
		line := strings.Fields(i)
		l := len(line)
		ov, err := strconv.Atoi(line[l-8])
		scu.QErr(err)
		rscu, err := strconv.Atoi(line[l-5])
		scu.QErr(err)
		if ov != 1 && rscu != 1 {
			continue //the original rule by Poma et al.
		}
		a1, err := strconv.Atoi(line[5])
		scu.QErr(err)
		a2, err := strconv.Atoi(line[9])
		scu.QErr(err)
		chain1 := line[l-15]
		chain2 := line[l-11]
		R, err := strconv.ParseFloat(line[l-9], 64)
		scu.QErr(err)
		if IsContiguous(a1, a2, exclude, chain1, chain2) {
			continue
		}
		index1, err := chem.MolIDNameChain2Index(mol, a1, "CA", chain1)
		scu.QErr(err)
		index2, err := chem.MolIDNameChain2Index(mol, a2, "CA", chain2)

		r := &chem.Bond{Index: 0, Dist: R, Energy: e, At1: mol.Atom(index1), At2: mol.Atom(index2)}
		if ov != 1 {
			r.Order = 1.0 //I'm abusing this field to signal an rCSU-only term.
		}
		set = append(set, r)
	}
	set2 := make([]*chem.Bond, 0, len(set))
	repeatedTest := make([][2]int, 0, 30)
	for i, v := range set {
		if v.Order > 0 && !CieplakIsSymmetric(i, set) {
			fmt.Println("Non-symmetric potential excluded", v)
			continue
		}
		if isinArray([2]int{v.At1.ID, v.At2.ID}, repeatedTest) {
			continue
		}
		repeatedTest = append(repeatedTest, [2]int{v.At1.ID, v.At2.ID})
		v.Order = 0 //I was abusing this field, now I reset it to the value it should actually have, all these are short range
		set2 = append(set2, v)
	}

	set3 := make([]*chem.Bond, 0, len(set2))
	for _, v := range set2 {
		if CieplakIsThere(v, set3) {
			continue
		}
		set3 = append(set3, v)
	}
	for i, v := range set3 {
		v.At1.Bonds = append(v.At1.Bonds, v)
		v.At2.Bonds = append(v.At2.Bonds, v)
		v.Index = i + lastID
	}
	return set3

}

func isinArray(test [2]int, set [][2]int) bool {
	t := test
	for _, v := range set {
		if (t[0] == v[0] && t[1] == v[1]) || (t[0] == v[1] && t[1] == v[0]) {
			return true
		}

	}
	return false

}
