/*

Chasqui, an allosteric-path finder.

This program makes extensive use of the goChem Computational Chemistry library.
If you use this program, we kindly ask you support it by to citing the library as:

R. Mera-Adasme, G. Savasci and J. Pesonen, "goChem, a library for computational chemistry", http://www.gochem.org.


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
	"flag"
	"fmt"
	"math"
	"math/rand"
	"os"
	"sort"
	"strings"

	chem "github.com/rmera/gochem"
	path "gonum.org/v1/gonum/graph/path"

	"github.com/rmera/gochem/chemgraph"
	"github.com/rmera/scu"
)

func Qerr(err error) {
	if err != nil {
		panic(err.Error())
	}
}

/*
type Bond struct {
	Index  int
	At1    *Atom
	At2    *Atom
	Dist   float64
	Energy float64 //One might prefer to use energy insteaf of order.
	//NOTE
	//I'm not sure if I leave just the order and let the user "Abuse" that field for energy
	//or anything else they want, or if I keep this extra field.
	Order float64 //Order 0 means undetermined
}

*/

func parseline(mol chem.Atomer, line string, dcutoff, ecutoff float64) *chem.Bond {
	rep := strings.Replace
	f := strings.Fields(line)
	chain1 := f[1]
	molid1 := scu.MustAtoi(strings.Replace(f[3], ":", "", -1))
	chain2 := f[5]
	//cutoff := 2.5 ////////////
	molid2 := scu.MustAtoi(strings.Replace(f[7], ":", "", -1))
	distance := scu.MustParseFloat(strings.Replace(f[len(f)-1], ",", "", -1)) //the replace is not really needed
	energy := scu.MustParseFloat(rep(rep(f[len(f)-2], ",", "", -1), ":", "", -1))
	//	fmt.Printf("Distance: %f, Energy: %f\n", distance, energy) //////////////////////
	if distance > dcutoff || math.Abs(energy) < ecutoff {
		return nil
	}
	index1, err := chem.MolIDNameChain2Index(mol, molid1, "CA", chain1)
	Qerr(err)
	index2, err := chem.MolIDNameChain2Index(mol, molid2, "CA", chain2)
	Qerr(err)
	b := &chem.Bond{Index: 0, Dist: distance, Energy: energy, At1: mol.Atom(index1), At2: mol.Atom(index2)}
	//Go-like cutoff
	gocutoff := 12 * 4.184
	if energy > gocutoff {
		b.Energy = gocutoff
	}
	return b

}

// Adds bonds from a file, plus, backbone bonds with a given energy
func QMBonds(mol *chem.Molecule, bondsf string, bbenergy, dcutoff, ecutoff float64) ([]*chem.Bond, error) {
	//var t1, t2 *v3.Matrix
	bonds := make([]*chem.Bond, 0, 10)
	bondslist, err := scu.NewMustReadFile(bondsf)
	Qerr(err)
	cont := 0
	for i := bondslist.Next(); i != "EOF"; i = bondslist.Next() {
		b := parseline(mol, i, dcutoff, ecutoff)
		if b == nil {
			continue
		}
		b.Index = cont
		at1 := mol.Atoms[b.At1.Index()]
		at2 := mol.Atoms[b.At2.Index()]
		at1.Bonds = append(at1.Bonds, b)
		at2.Bonds = append(at2.Bonds, b)
		bonds = append(bonds, b) //just to easily keep track of them.
		cont++
	}
	return bonds, nil
} //Ill not add all backbone bonds

func BBBonds(mol *chem.Molecule, bbenergy float64, lastbondid int) ([]*chem.Bond, error) {
	back := make([]*chem.Bond, 10)
	cont := lastbondid
	var lastMolID, lastCA int = 1000, -1
	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		if at.Name == "CA" {
			if at.MolID > lastMolID {
				back[len(back)-1].At2 = at
				at.Bonds = append(at.Bonds, back[len(back)-1])
			} else if lastCA >= 0 {
				//This means we finished reading a chain, and started reading another
				//there is a "dangling" bond that goes nowhere, in the previous chain, we got to remove it now.
				back = back[:len(back)-1]
				at2 := mol.Atom(lastCA)
				at2.Bonds = at2.Bonds[:len(at2.Bonds)-1]
			}
			//For the last CA of a chain, this bond will not have a "partner" so it needs to be removed.
			//Above we take care of the case of a chain that comes before another chain.
			//below we take care of the last chain (i.e. one that doesn't come before another)
			b := &chem.Bond{Index: cont, Dist: 1, Energy: bbenergy, At1: at, At2: nil}
			at.Bonds = append(at.Bonds, b)
			back = append(back, b)
			//update the counters
			lastCA = i
			lastMolID = at.MolID
			cont++
		}
	}
	//We fix the last dangling bond
	back = back[:len(back)-1]
	at2 := mol.Atom(lastCA)
	at2.Bonds = at2.Bonds[:len(at2.Bonds)-1]
	return back, nil
}

func checkmustbe(path string, mustbe []string) bool {
	if mustbe == nil {
		return true
	}
	for _, v := range mustbe {
		if !strings.Contains(path, v) {
			return false
		}
	}
	return true
}

func CaminitoAgreste(top *chemgraph.Topology, index1, index2 int, O *opts) ([]string, []float64) {
	if O == nil {
		panic("nil option structure. Definitely shouldn't happen")
	}
	weights := make([]float64, 0, O.maxpath)
	ret := make([]string, 0, O.maxpath)
	ori := top.Atoms.Atoms[index1]
	oldpath := ""
	var prevpens map[int]float64
	for i := 0; i < O.maxpath; {
		Paths := path.DijkstraFrom(ori, top)
		Path, _ := Paths.To(int64(index2))
		if len(Path) == 0 {
			break
		}
		pathstr := make([]string, 0, len(Path))
		pens := O.penfunc(len(Path)-1, O.worst, O.least)
		for j, v := range Path {
			index := int(v.ID())
			at := top.Atoms.Atoms[index]
			s := fmt.Sprintf("%s%d%s", at.MolName, at.MolID, at.Chain)
			pathstr = append(pathstr, s)
			if j%2 != 0 {
				a2 := top.Atoms.Atoms[int(Path[j-1].ID())]
				pr := false
				toundo := setpenalties(at, a2, top.Bonds, pens[j-1], pr)
				if O.undo {
					undopenalties(prevpens, top.Bonds)
				}
				prevpens = toundo
			}
		}
		newpath := strings.Join(pathstr, "-->")
		if newpath != oldpath {
			if checkmustbe(newpath, O.mustbe) {
				ret = append(ret, newpath)
				weights = append(weights, Paths.WeightTo(int64(index2)))
			}
			oldpath = newpath
			i++
		}
	}
	return ret, weights
}

func undopenalties(pen map[int]float64, B chemgraph.Bonds) {
	for k, v := range pen {
		B[k].Energy /= v
	}
}

type caminos struct {
	s []string
	w []float64
}

func (a caminos) Len() int { return len(a.s) }
func (a caminos) Less(i, j int) bool {
	//i.e. the number of elements
	return strings.Count(a.s[i], "-->") < strings.Count(a.s[j], "-->")
}
func (a caminos) Swap(i, j int) {
	a.s[i], a.s[j] = a.s[j], a.s[i]
	a.w[i], a.w[j] = a.w[j], a.w[i]
}

// applies function f to every bond in B that connects a1 and a2, in either order.
func setpenalties(a1, a2 *chemgraph.Atom, B chemgraph.Bonds, penalty float64, pr bool) map[int]float64 {
	ret := make(map[int]float64)
	for i, v := range B {
		if (v.At1.ID() == a1.ID() && v.At2.ID() == a2.ID()) || (v.At1.ID() == a2.ID() && v.At2.ID() == a1.ID()) {
			v.Energy = v.Energy * penalty
			//	fmt.Println("Changed the weeeight!!!") ///////
			ret[i] = v.Energy
			if pr {
				fmt.Println(v.Weight())
			}
		}
	}
	return ret
}

type opts struct {
	penfunc      func(int, float64, float64) []float64
	undo         bool
	maxpath      int
	worst, least float64
	mustbe       []string
}

func main() {
	pen := flag.String("pen", "random", "Penalty function to use. Options: wedge, random")
	undo := flag.Bool("undo", false, "Undo the penalties associated to a path after the next path is found")
	worst := flag.Float64("worst", 0.8, "Worst penalty, between 0 and 1, the lower, the worse")
	least := flag.Float64("least", 0.95, "Least penalty, between 0 and 1, the lower, the worse")
	maxpath := flag.Int("maxpath", 10, "Maximum number of paths to obtain (less could be shown)")
	exclude := flag.Int("exclude", 2, "Exclude interactions between residues separated by less than this number of residues. Only for Cieplak potentials")
	bbenergy := flag.Float64("bbenergy", 1.0, "Energy to assign to backbone bonds")
	crosschains := flag.Bool("crosschains", false, "Don't discourage paths from crossing chains")
	dcutoff := flag.Float64("dcutoff", 2.8, "Distance cutoff for interactions")
	ecutoff := flag.Float64("ecutoff", 0.6, "Energy cutoff for interaction (the default is ~RT)")
	sortpaths := flag.Bool("sort", false, "Sort paths by length")
	cieplake := flag.Float64("cieplake", 2.87, "Energy for Cieplak's Go contacts")
	cieplakfile := flag.String("cieplakfile", "", "File with Cieplak's Go contacts")
	qmfile := flag.String("qmfile", "", "File with QM contacts")
	mustbestr := flag.String("mustbe", "", "Residues that must be in the shown paths, RESNAMERESIDCHAIN (ex. HIS192A), separated by commas")
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Chasqui: Allosteric pathways between 2 residues.\n Usage:\n  %s [flags] geomtry.pdb ResidueID1 ResidueID2 Chain1 Chain2 \n\nFlags:\n", os.Args[0])
		flag.PrintDefaults()
	}
	flag.Parse()
	args := flag.Args()
	//we first prepare the system
	mol, err := chem.PDBFileRead(args[0])
	Qerr(err)
	mol.FillIndexes()
	cbonds := make([]*chem.Bond, 0)
	if *qmfile != "" {
		cbonds, err = QMBonds(mol, *qmfile, *bbenergy, *dcutoff, *ecutoff)
		Qerr(err)
	}
	if *cieplakfile != "" {
		lastid := 0
		if len(cbonds) != 0 {
			lastid = cbonds[len(cbonds)-1].Index
		}
		sbonds := CieplakPotRead(mol, *cieplakfile, *exclude, *cieplake, lastid)
		cbonds = append(cbonds, sbonds...)
	}
	lastid := 0
	if len(cbonds) != 0 {
		lastid = cbonds[len(cbonds)-1].Index
	}
	_, err = BBBonds(mol, *bbenergy, lastid)
	scu.QErr(err)
	top := chemgraph.TopologyFromChem(mol, nil, func(B *chemgraph.Bond) float64 { return 1 / math.Abs(B.Energy) }) //by default the weight is 1/Abs(Energy)
	top.SymmetrizePath()
	atoi := scu.MustAtoi
	index1, err := chem.MolIDNameChain2Index(mol, atoi(args[1]), "CA", args[3])
	index2, err := chem.MolIDNameChain2Index(mol, atoi(args[2]), "CA", args[4])
	bonds := make([]*chemgraph.Bond, 0)
	if !*crosschains {
		for _, v := range top.Bonds {
			if v.At1.Chain == v.At2.Chain {
				bonds = append(bonds, v)
			}
		}
		top.Bonds = bonds

	}
	var mustbe []string
	if *mustbestr != "" {
		mustbe = strings.Split(*mustbestr, ",")
	}
	opts := &opts{penfunc: penaltymap[*pen], undo: *undo, worst: *worst, least: *least, maxpath: *maxpath, mustbe: mustbe}
	Paths, weights := CaminitoAgreste(top, index1, index2, opts)
	c := caminos{s: Paths, w: weights}
	if *sortpaths {
		sort.Sort(c)
	}
	for i, v := range c.s {
		fmt.Printf("%s  abs(w):%3.2f\n", v, c.w[i])
	}

}

//penalty functions

var penaltymap = map[string]func(int, float64, float64) []float64{
	"random": randompen,
	"wedge":  wedgepen,
}

func randompen(listlen int, worst, least float64) []float64 {
	//	fmt.Println("begin random penalty") ///////////////
	pen := make([]float64, listlen)
	for i := 0; i < listlen; i++ {
		fac := rand.Float64()
		pen[i] = (1-fac)*least - fac*worst
	}
	//	fmt.Println("end random penalty", pen) ///////////////
	return pen
}

// create penalties that get worse the closest you are to the center of the list ina wedgie shape "^"
func wedgepen(listlen int, worst, least float64) []float64 {
	disp := 0
	pen := make([]float64, listlen)
	if listlen%2 != 0 {
		pen[listlen/2] = worst
		disp = 1
	}
	halfway := listlen/2 - 1
	fn := func(x int) float64 {
		return least + worst*float64(x)/float64(halfway+disp)
	}
	for i := 0; i < halfway; i++ {
		pen[i] = fn(i)
		pen[listlen-i-1] = fn(i)
	}
	return pen

}
