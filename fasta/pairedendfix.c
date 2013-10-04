#include <stdio.h>
#include <stdlib.h>
#include "fasta.h"
#include "strmap.h"

int main(int argc, char **argv){

    /* aux */
    void showHelp(void);
    typedef enum { false, true } bool;
    int x=0, nf=0, f=0, cut=0;
    bool help=false;

    /* Files */
    char in1[255],in2[255], out1[255],out2[255], outf[255];
    FILE *_out1, *_out2, *_outf;

    /* fasta */
    FASTAFILE *ffp;
    char *seq;
    char *name, *nametmp;
    int L, i;

    /* hash */
    StrMap *sm;
    char buf[2];
    int result;

    memset(in1, 0, sizeof(in1));
    memset(in2, 0, sizeof(in2));
    memset(out1, 0, sizeof(out1));
    memset(out2, 0, sizeof(out2));
    memset(outf, 0, sizeof(outf));

    if(argc==0 || argc==1){
        help = true;
        showHelp();
        return;
    }else{
        for(i = 0; i < argc; i++){
            /*printf("arg %d: %s\n", i, argv[i]);*/
            if(!strcmp("-1", argv[i])){
                i++;
                strcpy (in1,argv[i]);
                continue;
            }

            if(!strcmp("-2", argv[i])){
                i++;
                strcpy (in2,argv[i]);
                continue;
            }

            if(!strcmp("-o1", argv[i])){
                i++;
                strcpy (out1,argv[i]);
                continue;
            }

            if(!strcmp("-o2", argv[i])){
                i++;
                strcpy (out2,argv[i]);
                continue;
            }

            if(!strcmp("-of", argv[i])){
                i++;
                strcpy (outf,argv[i]);
                continue;
            }

            if(!strcmp("-c", argv[i])){
                i++;
                cut = atoi(argv[i]);
                continue;
            }

            if(!strcmp("--help", argv[i]) || !strcmp("-h", argv[i])){
                help = true;
                showHelp();
                return;
            }

        }

    }

	if ((fopen(in1,"r")) == NULL) {
		printf("Error - invalid file: %s\n",in1);
		return;
	}

	if ((fopen(in2,"r")) == NULL) {
		printf("Error - invalid file: %s\n",in2);
		return;
	}

	if ((_out1 = fopen(out1,"w")) == NULL) {
		printf("Error - invalid file: %s\n",out1);
		return;
	}

	if ((_out2 = fopen(out2,"w")) == NULL) {
		printf("Error - invalid file: %s\n",out2);
		return;
	}

	if ((_outf = fopen(outf,"w")) == NULL) {
		printf("Error - invalid file: %s\n",outf);
		return;
	}


    printf("Input:\n   %s\n   %s\nOutput:\n   %s\n   %s\n   %s\n", in1, in2, out1, out2, outf);

    /* Count reads */
    ffp = OpenFASTA(in1);
    while (ReadFASTA(ffp, &seq, &name, &L)){
        x++;
        free(name);
        free(seq);
    }
    CloseFASTA(ffp);

    /* inicia hash */
    sm = sm_new(x);

    /* Carrega dados na hash */
    ffp = OpenFASTA(in1);
    while (ReadFASTA(ffp, &seq, &name, &L)){
        name=strndup(name, strlen(name)-cut);
        sm_put(sm, name, "1");
        free(name);
        free(seq);
    }
    CloseFASTA(ffp);

    /* Verifica com arquivo 2 */
    ffp = OpenFASTA(in2);
    while (ReadFASTA(ffp, &seq, &name, &L)){
        nametmp = strndup(name, strlen(name)-cut);
        result = sm_get(sm, nametmp, buf, sizeof(buf));
        /* Caso não encontre, coloca read em arquivo de fragmentos */
        if(result==0){
            nf++;
            fprintf(_outf,">%s\n%s", name, seq);
        }else{ /* Caso encontre, coloca read em arquivo correto e marca como 0*/
            f++;
            fprintf(_out2,">%s\n%s", name, seq);
            sm_put(sm, nametmp, "0");
        }
        free(seq);
        free(name);
        free(nametmp);
    }
    CloseFASTA(ffp);

    /* Verifica com arquivo 1 */
    ffp = OpenFASTA(in1);
    while (ReadFASTA(ffp, &seq, &name, &L)){
        nametmp = strndup(name, strlen(name)-cut);
        result = sm_get(sm, nametmp, buf, sizeof(buf));
        /* Caso esteja marcado com 0, coloca read em arquivo correto */
        if(buf[0]=='0'){
            fprintf(_out1,">%s\n%s", name, seq);
        }else{ /* Caso esteja marcado com 1, não encontrou no passo anterior e coloca read em arquivo de fragmentos */
            nf++;
            fprintf(_outf,">%s\n%s", name, seq);
        }
        free(seq);
        free(name);
        free(nametmp);
    }
    CloseFASTA(ffp);

    sm_delete(sm);
    fclose(_out1);
    fclose(_out2);
    fclose(_outf);

    printf("Paired sequences: %d (%d each file)\nFragment sequences: %d\n", f*2, f, nf);

}

void showHelp(){
    printf("pairedendfix -1 input_file1 -2 input_file2 -o1 output_file1 -o2 output_file2 -of output_fragment -c cut_header_len\n");
    printf("\n -c INT \t Length to cut the header to compare");
    printf("\n -h/--help \t Show this help\n");
}
