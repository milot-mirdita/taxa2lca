package main

import (
	"bytes"
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
	"time"

	"github.com/milot-mirdita/mmseqs-web-backend/dbreader"
)

var (
	nodesflag  string
	namesflag  string
	taxdbflag  string
	outflag    string
	taxidxflag string
	taxlevel   string
	helpflag   bool
)

func init() {
	flag.StringVar(&nodesflag, "nodes", "nodes.dmp", "nodes.dmp file of taxonomy")
	flag.StringVar(&namesflag, "names", "names.dmp", "names.dmp file of taxonomy")
	flag.StringVar(&taxdbflag, "database", "taxons_db", "Input database name")
	flag.StringVar(&taxidxflag, "index", "taxons_db.index", "Input index name")
	flag.StringVar(&outflag, "output", "taxa.tsv", "Output tsv file")
	flag.StringVar(&taxlevel, "levels", "", "Desired LCA taxonomical levels [optional]")
	flag.BoolVar(&helpflag, "help", false, "Print USAGE and exits")
	flag.Parse()

	// The blast file is the first unparsed argument
	if helpflag {
		fmt.Printf("taxa2lca\n")
		flag.Usage()
		os.Exit(0)
	}
}

func main() {
	levs := bytes.Split([]byte(taxlevel), []byte{':'})
	taxDB, err := NewTaxonomy(nodesflag, namesflag)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Taxonomy error: %s\n", err)
		os.Exit(1)
	}

	t1 := time.Now()

	reader := dbreader.Reader{}
	reader.Make(taxdbflag, taxidxflag)

	defer reader.Delete()

	f, err := os.Create(outflag)
	if err != nil {
		log.Panic("Could not create output file!")
	}
	defer f.Close()

	var i int64
	for i = 0; i < reader.Size(); i++ {
		data := reader.Data(i)
		if data == "" {
			continue
		}

		split := strings.Split(data, "\n")
		taxa := make([]int, len(split))

		query := strconv.FormatUint(uint64(reader.Key(i)), 10)
		for i, e := range strings.Split(data, "\n") {
			taxon, _ := strconv.ParseInt(e, 10, 64)
			taxa[i] = int(taxon)
		}

		var atLevs [][]byte
		var allLevs []byte
		lcaNode, err := taxDB.LCA(taxa...)
		if err != nil {
			atLevs = make([][]byte, 1)
			atLevs[0] = Unknown
		} else {
			if len(levs[0]) != 0 {
				atLevs = taxDB.AtLevels(lcaNode, levs...)
			}
		}
		allLevs = bytes.Join(atLevs, []byte{';'})
		f.WriteString(fmt.Sprintf("%s\t%d\t%s\t%s\t%s\n", query, lcaNode.Taxid, lcaNode.Name, lcaNode.Taxon, allLevs))
	}

	t2 := time.Now()
	dur := t2.Sub(t1)
	secs := dur.Seconds()
	log.Printf("Done in %.3f seconds\n", secs)
}
