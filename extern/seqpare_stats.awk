#!/bin/awk -f

func print_record (predicted,ground,scores,s_sum,tool) {
    printf tool"\t";
    printf "%d\t",length(ground);
    printf "%d\t",length(predicted);
    printf "%d\t",length(scores);
    printf "%.2f%%\t",100.0*length(scores)/length(ground);
    printf "%.2f%%\t",100.0*length(scores)/length(predicted);
    printf "%1.5f\t",s_sum/length(scores);
    printf "%1.5f\t",scores[int(length(scores)/2)];
    printf "%d\t",length(ground)-length(scores);
    printf "%.2f%%\t",100.0*(length(ground)-length(scores))/length(ground);
    printf "%d\t",length(predicted)-length(scores);
    printf "%.2f%%\n",100.0*(length(predicted)-length(scores))/length(predicted);
}

BEGIN{
    tool="";
    s=2;
    p=3;
    g=4;
    printf "tool\t"
    printf "gr_tot\t"
    printf "pr_tot\t"
    printf "pairs\t"
    printf "gr_pair_prct\t"
    printf "pr_pair_prct\t"
    printf "sp_avg\t"
    printf "sp_med\t"
    printf "gr_unpair\t"
    printf "gr_unpair_prct\t"
    printf "pr_unpair\t"
    printf "pr_unpair_prct\n"
}

FNR == 1 {
    if (tool!="") {
        print_record(predicted,ground,scores,s_sum,tool);
    }
    delete predicted;
    delete ground;
    delete scores;
    s_sum = 0;
    tool = split(FILENAME,a,"/"); f=a[length(a)]; tool=substr(f, 1, match(f, "pairings")-2);
}

$p != "None" {
    predicted[FNR];
}
$g != "None" {
    ground[FNR];
}
$p != "None" && $g != "None" {
    scores[FNR]=$s;
    s_sum += $s;
}

END {
    print_record(predicted,ground,scores,s_sum,tool);
}
