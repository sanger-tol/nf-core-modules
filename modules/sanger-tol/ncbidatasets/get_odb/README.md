# get_odb.py

Originally written by Matthieu Muffato (mm49)

Modified by Damon-Lee Pointon (dp24)

## Folder Structure

```
.
├── resources
│   └── usr
│       └── assets
│       │   ├── mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt
│       │   └── mapping_taxids-busco_dataset_name.eukaryota_odb12.2025-01-15.txt
│       └── bin
│           └── get_odb.py
├── tests
│   ├── main.nf.test
│   ├── main.nf.snapshot
│   └── nextflow.config
│
├── environment.yml
├── README.md
├── main.nf
└── meta.yml
```

## Sources

ODB mapping files are provided by the [BUSCO](https://busco.ezlab.org/) project at the `busco-data.ezlab.org` data repo:

- [odb10 file](https://busco-data.ezlab.org/v5/data/placement_files/mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt) - Created 2019-12-16
- [odb12 file](https://busco-data.ezlab.org/v5/data/placement_files/mapping_taxids-busco_dataset_name.eukaryota_odb12.2025-01-15.txt) - Created 2025-01-15
