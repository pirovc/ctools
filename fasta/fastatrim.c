#include <stdio.h>
#include <stdlib.h>
#include "fasta.h"

int main(int argc, char **argv){

    /* aux */
    void showHelp(void);
    typedef enum { false, true } bool;
    int trim=0, i=0;
    char f[1];
    bool find=false, help=false, brk=false;

    /* fasta */
    FASTAFILE *ffp;
    int   L;
    char *name, *seq;

    /* Files */
    char in[255];

    if(argc==0 || argc==1){
        help = true;
        showHelp();
        return;
    }else{
        for(i = 0; i < argc; i++){
            if(!strcmp("-i", argv[i])){
                i++;
                strcpy (in,argv[i]);
                continue;
            }

            if(!strcmp("-t", argv[i])){
                i++;
                trim = atoi(argv[i]);
                continue;
            }

            if(!strcmp("-f", argv[i])){
                i++;
                strcpy (f,argv[i]);
                find=true;
                continue;
            }

            if(!strcmp("--help", argv[i]) || !strcmp("-h", argv[i])){
                help = true;
                showHelp();
                return;
            }
        }
    }


	if ((fopen(in,"r")) == NULL) {
		printf("Error - invalid file: %s\n",in);
		return;
	}

    ffp = OpenFASTA(in);
    while (ReadFASTA(ffp, &seq, &name, &L)){
        if(trim==0) {trim=strlen(seq);}

        for(i=0;i<trim;i++){
            if(find==true){
                if(seq[i]==f[0]){
                    brk=true;
                    break;
                }
            }
        }

        if(brk==false){
            printf(">%s\n", name);
            for(i=0;i<trim;i++){
                printf("%c", seq[i]);
            }
            printf("\n");
        }

        brk=false;
        free(name);
        free(seq);
    }
    CloseFASTA(ffp);

}

void showHelp(){
    printf("fastatrim -i FILE input_file [-t INT new_length_of_sequence] [-f CHAR remove_sequences_with_specific_char] > out.fasta\n");
}
