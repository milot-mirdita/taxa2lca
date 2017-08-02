package main

import (
	"bytes"
	"flag"
	"fmt"
	"github.com/milot-mirdita/mmseqs-web-backend/dbreader"
	"log"
	"os"
	"runtime"
	"strconv"
	"strings"
	"time"
)

var (
	procsflag  int
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
	flag.StringVar(&taxdbflag, "taxid database", "taxons_db", "Input database name")
	flag.StringVar(&taxidxflag, "taxid index", "taxons_db.index", "Input index name")
	flag.StringVar(&outflag, "Output files", "taxa.tsv", "Output tsv files")
	flag.IntVar(&procsflag, "nprocs", 4, "Number of cpus for multithreading [optional]")
	flag.StringVar(&taxlevel, "levels", "", "Desired LCA taxonomical levels [optional]")
	flag.BoolVar(&helpflag, "help", false, "Print USAGE and exits")
	flag.Parse()

	// The blast file is the first unparsed argument
	if helpflag {
		fmt.Printf("taxa2lca\n")
		flag.Usage()
		os.Exit(0)
	}

	runtime.GOMAXPROCS(procsflag)
}

type Job struct {
	Rank  int64
	Start int64
	Size  int64
}

func lca(input <-chan Job, taxDB *Taxonomy, reader *dbreader.Reader, levs [][]byte, outbase string) {

	for {
		select {
		case job, ok := <-input:
			if ok {
				f, err := os.Create(outbase + "." + strconv.FormatInt(job.Rank, 10))
				if err != nil {
					log.Panic("Could not create output file!")
				}
				defer f.Close()

				log.Printf("Computing split from %s to %s", job.Start, job.Start+job.Size)
				for i := job.Start; i < job.Start+job.Size; i++ {
					data := reader.Data(i)
					split := strings.Split(data, "\n")
					taxa := make([]int, len(split))
					var query string
					for i, e := range strings.Split(data, "\n") {
						values := strings.Split(e, "\t")
						query = values[0]
						taxon, _ := strconv.ParseInt(values[1], 10, 64)
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
					f.WriteString(fmt.Sprintf("%s\t%s\t%s\t%s\n", query, lcaNode.Name, lcaNode.Taxon, allLevs))
				}
			} else {
				return
			}
		default:
		}
	}
}

func decomposeDomain(domain_size, world_rank, world_size int64) Job {
	if world_size > domain_size {
		return Job{world_rank, 0, 0}
	}

	subdomain_start := domain_size / world_size * world_rank
	subdomain_size := domain_size / world_size

	if world_rank == world_size-1 {
		subdomain_size += domain_size % world_size
	}

	return Job{world_rank, subdomain_start, subdomain_size}

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

	jobs := make(chan Job, 1)

	for i := 0; i < procsflag; i++ {
		jobs <- decomposeDomain(reader.Size(), int64(i), int64(procsflag))
	}

	go lca(jobs, taxDB, &reader, levs, outflag)

	reader.Delete()

	t2 := time.Now()
	dur := t2.Sub(t1)
	secs := dur.Seconds()
	log.Printf("Done in %.3f seconds\n", secs)
}
