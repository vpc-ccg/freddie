#!/bin/awk -f

BEGIN{
    s=2;
    p=3;
    g=4;
    OFS="\t";
}

FNR == NR {
    if ($p != "None") {
        pg[$p]=$g;
    }
    if ($g != "None") {
        gp[$g]=$p;
    }
}

FNR != NR {
    read_g = substr($1, 1, match($1,"_")-1)
    if (!(read_g in gp)) {
        gp[read_g] = "OFF_CHROM";
    }
    print $1,$2,pg[$2],gp[read_g], $2==gp[read_g];
}

END {
}
