import sys

def usage():
    print( "Usage: cat the tsv from stdin, give output tsv path as  first parameter the min score with the second parameter.!!",file = sys.stderr)
    return -1

def main():
    if len(sys.argv) < 3:
        return usage()
    out_tsv = sys.argv[1]
    min_score = int(sys.argv[2])

    with open(out_tsv, 'w') as thand:
        for line in sys.stdin:
            fields = line.rstrip().split("\t")
            score = int(fields[4])
            polyAbegin = int(fields[2])
            print(">",end="")
            print(fields[0])
            if score >= min_score:
                print(fields[5][:polyAbegin])
                print("{}\t1".format(fields[0]),file=thand)
            else:
                print(fields[5])
                print("{}\t0".format(fields[0]),file=thand)

    return 0

if __name__ == "__main__":
    exit(main())
