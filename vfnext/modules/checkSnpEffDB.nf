process checkSnpEffDB{
    publishDir "${params.outDir}", mode: "copy"
    label "singlethread"
    
    input:
        val(genome_code)

    output:
        path("snpEffDB_entry_found.log")
$/
#!/usr/bin/env python

with open("${params.snpEffDBCatalog}", 'r') as srch_fl:
    dct_lst = []
    for line in srch_fl:
        d_line = line.replace(' ', '').split('\t')
        dct = {
            "gnm": d_line[0],
            "organism": d_line[1],
            "status": d_line[2],
            "bundle": d_line[3],
            "download-link":d_line[4]
        }

        if "${genome_code}" in dct["gnm"]:
            dct_lst.append(dct)
    
    out_fl = open("snpEffDB_entry_found.log", 'w')
    header = "gnm,organism,status,bundle,download-link\n"
    out_fl.write(header)
    
    n_founds = len(dct_lst)
    if n_founds == 0:
        print("ERROR: ${genome_code} is not available at SnpEFF Database")
        exit(1)
    if n_founds == 1:
        data = dct_lst[0]
        line = f"{data['gnm']},{data['organism']},{data['status']},{data['bundle']},{data['download-link']}\n"
        out_fl.write(line)
    if n_founds > 1:
        print("ERROR: more than one entry found for ${genome_code}")
        print("       check ${params.outDir}/snpEffDB_entry_found.log")
        for data in dct_lst:
            line = f"{data['gnm']},{data['organism']},{data['status']},{data['bundle']},{data['download-link']}\n"
            out_fl.write(line)
        exit(1)
/$
}
