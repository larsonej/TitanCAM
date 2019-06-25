/*  @(#) byteorder.c  McKie  Sep-1996
 *  This program converts an aerosol model's binary history 
 *  file written in standard unix fortran sequential
 *  record organization between little endian (e.g. Intel)
 *  and big endian (e.g. Sun) byte ordering.
 */


/*  Include system stuff  */

#include <stdio.h>


/*  Define symbolic constants  */

#define PROGNAME	"byteorder"
#define USAGE		"Usage: byteorder old_history_file new_history_file"

#define MIN_BUF_GET	1024    /* min # bytes to get for i/o buffer from system */

#define NBYTE_INT	4	/* # bytes in an fortran INTEGER */
#define NBYTE_REAL	4	/* # bytes in a fortran REAL */
#define NBYTE_DOUBLE	8       /* # bytes in a fortran DOUBLE PRECISION */
#define NBYTE_LOGICAL	4       /* # bytes in a fortran LOGICAL */


/*  Declare 4 byte integer for old & new record byte count bit pattern  */

int n_rec_old;
int n_rec_new;


/*  Declare ptr to i/o buffer to hold each fortran record  */

unsigned char *buf;


/*  Declare & initialize current size of i/o buffer  */

unsigned int nbuf = 0; 


/*  Declare file handles for the input & output streams  */

FILE *in_file;
FILE *ou_file;


/*  Declare ptr to file name for input & output history files  */

char *in_name;
char *ou_name;


/*  Declare & initialize record counter  */

int k_record = 0;


/*  Declare unions to overlap byte storage with ints or doubles  */

union { char b[NBYTE_INT];    int i;    } itime;
union { char b[NBYTE_DOUBLE]; double d; } time;




main(argc, argv)
int argc;
char *argv[];
{
 int irec;
 int i;


 /*  Check for expected # cmd line args (cmd name counts as 1 arg)  */

 if( argc != 3 )
 {
  printf("%s--Error: Unexpected number of cmd line arguments\n", PROGNAME);
  printf("%s\n", USAGE);
  exit(1);
 }


 /*  Point to names of input & output history files  */

 in_name = argv[1];
 ou_name = argv[2];


 /*  Attempt to open old history file for input  */

 if( ( in_file = fopen( in_name, "r" ) ) == (FILE *)(NULL) )
 {
  printf("%s--Error: Cant open input history file %s\n", PROGNAME, in_name);
  exit(1);
 }


 /*  Attempt to open new history file for output  */

 if( ( ou_file = fopen( ou_name, "w" ) ) == (FILE *)(NULL) )
 {
  printf("%s--Error: Cant open output history file %s\n", PROGNAME, ou_name);
  exit(1);
 }


 /*  Debug print the # bytes in each record.  Normally commented  */

 /*

 while( getrec() != EOF )
 {
  printf("Debug. Record # %d:   %7d bytes\n", k_record, n_rec_old); 
 }
 fseek(in_file, 0, SEEK_SET);
 k_record = 0;

 */


 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 


 /*  Input, convert, output header records  */

 printf("%s--Begin converting header records ...\n", PROGNAME);


 /*  Program name & version  */

 if( getrec() == EOF ) { badeof("Header: prog name & version"); } 
 outbuf();


 /*  simtitle  */

 if( getrec() == EOF ) { badeof("Header: simtitle"); } 
 outbuf();


 /*  ibtime, ietime, nhist  */

 if( getrec() == EOF ) { badeof("Header: ibtime, ietime, nhist"); } 
 revbuf( NBYTE_INT );
 outbuf();


 /*  MAXNX, MAXNY, ...  */

 if( getrec() == EOF ) { badeof("Header: MAXNX, MAXNY, ..."); } 
 revbuf( NBYTE_INT );
 outbuf();


 /*  nx, ny, ...  */

 if( getrec() == EOF ) { badeof("Header: nx, ny, ..."); } 
 revbuf( NBYTE_INT );
 outbuf();


 /*  itype(), igelem(), ...  */

 for(irec=0; irec<5; irec++)
 {
  if( getrec() == EOF ) { badeof("Header: itype(), igelem(), ..."); } 
  revbuf( NBYTE_INT );
  outbuf();
 }


 /*  groupname(), elemname(), ...  */

 for(irec=0; irec<4; irec++)
 {
  if( getrec() == EOF ) { badeof("Header: names"); } 
  outbuf();
 }


 /*  r, dr, rmass, dm, alt_mid3, xc, u ; each in its own record  */

 for(irec=0; irec<7; irec++)
 {
  if( getrec() == EOF ) { badeof("Header: r, dr, rmasss, dm, alt_mid3, xc, u"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();
 }


 /*  Report end of header records processing  */

 printf("%s--End converting header records at record # %d\n", PROGNAME, k_record);


 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 


 /*  Input, convert, output timepoint records until eof  */

 while( getrec() != EOF )
 {

  /*  itime, time  */

  for(i=0; i<NBYTE_INT; i++) { itime.b[i] = buf[i]; } 
  for(i=0; i<NBYTE_DOUBLE; i++) { time.b[i] = buf[NBYTE_INT+i]; } 
  printf("%s--Converting timepoint: rec %5d  itime=%7d  time=%f\n",
          PROGNAME, k_record, itime.i, time.d);

  revbytes( &buf[0], NBYTE_INT, 1 );
  revbytes( &buf[NBYTE_INT], NBYTE_DOUBLE, 1 );
  outbuf();


  /*  pc3  */

  if( getrec() == EOF ) { badeof("pc3"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();


  /*  gc3  */

  if( getrec() == EOF ) { badeof("gc3"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();


  /*  p3  */

  if( getrec() == EOF ) { badeof("p3"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();


  /*  t3  */

  if( getrec() == EOF ) { badeof("t3"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();


  /*  pt3  */

  if( getrec() == EOF ) { badeof("pt3"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();


  /*  pti3  */

  if( getrec() == EOF ) { badeof("pti3"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();


  /*  supsatl3  */

  if( getrec() == EOF ) { badeof("supsatl3"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();


  /*  supsati3  */

  if( getrec() == EOF ) { badeof("supsati3"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();


  /*  w3  */

  if( getrec() == EOF ) { badeof("w3"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();


  /*  dkv3  */

  if( getrec() == EOF ) { badeof("dkv3"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();


  /*  rnuclg  */

  if( getrec() == EOF ) { badeof("rnuclg"); } 
  revbuf( NBYTE_DOUBLE );
  outbuf();


  /*  Go attempt to input next timepoint records  */

 }


 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 


 /*  Close input & output history files  */

 fclose( in_file );
 fclose( ou_file );


 /*  Terminate normally  */

 exit(0); 

}


/*===================================================================================*/

/*  @(#) getrec.c  McKie  Sep-1996
 *  This routine inputs the next fortran record, checks header and trailer
 *  byte counts, expands the i/o buffer if necessary, counts records
 *
 *  Returns:  # data bytes in record
 */

int getrec()
{
 int trailer; 
 int n_buf_get;


 /*  Increment record counter to indicate # of this attempted input  */

 k_record++;


 /*  Input record size  */

 if( fread(&n_rec_old, 4, 1, in_file) != 1 )
 {
  printf("%s--EOF at record # %d\n", PROGNAME, k_record);
  return EOF;
 }


 /*  Make i/o buffer bigger if current record is bigger than old buffer  */

 if( n_rec_old > nbuf )
 {
  if( nbuf > 0 )
  {
   free( buf );
  }
  n_buf_get = ( n_rec_old > MIN_BUF_GET ) ? n_rec_old : MIN_BUF_GET;
  if( ( buf = (char *)malloc( n_buf_get ) ) == NULL )
  {
   printf("%s--Error: Cant get memory for i/o buffer [%d bytes]\n",
                   PROGNAME, n_buf_get);
   printf("%s--At record # %d\n", PROGNAME, k_record);
   exit(1);
  }
  nbuf = n_buf_get;
  /* printf("Debug. getrec. New buffer size: %d\n", nbuf); */
 } 


 /*  Input the data bytes in the record  */

 if( fread(buf, 1, n_rec_old, in_file) != n_rec_old )
 {
  printf("%s--Error: Cant input data record [%d bytes]\n",
                  PROGNAME, n_rec_old);
  printf("%s--At record # %d\n", PROGNAME, k_record);
  exit(1);
 }


 /*  Input 4 byte integer trailer record byte count */

 if( fread(&trailer, 4, 1, in_file) != 1 )
 {
  printf("%s--Error: Cant input record trailer byte count\n", PROGNAME);
  printf("%s--At record # %d\n", PROGNAME, k_record);
  exit(1);
  if( trailer != n_rec_old )
  {
   printf("%s--Error: Trailer byte count [%d] not same as header count [%d]\n",
                   PROGNAME, trailer, n_rec_old);
   printf("%s--At record # %d\n", PROGNAME, k_record);
   exit(1);
  }
 }


 /*  Return to caller with byte count of record   */

 return n_rec_old;

}


/*===================================================================================*/

/*  @(#) badeof.c  McKie  Sep-1996
 *  This routine handles premature eof condition.
 *
 *  Input:   msg = message string
 *
 *  Returns:  { no return }
 */

int badeof( msg )
char *msg;
{
 printf("%s--Unexpected eof at record %d\n", PROGNAME, k_record);
 printf("%s--%s\n", PROGNAME, msg);
 exit(1);
}


/*===================================================================================*/

/*  @(#) outbuf.c  McKie  Sep-1996
 *  This routine outputs the data bytes in the current i/o buffer
 *  as a fortran sequential binary record to the output file.
 *
 *  Returns:  { nothing meaningful }
 */

int outbuf()
{

 /*  Convert header byte count & output to file  */

 n_rec_new = n_rec_old;
 revbytes(&n_rec_new, 4, 1);
 fwrite(&n_rec_new, 4, 1, ou_file);


 /*  Output data bytes to file  */

 fwrite(buf, 1, n_rec_old, ou_file);


 /*  Append trailer byte count to file  */

 fwrite(&n_rec_new, 4, 1, ou_file);


 /*  Return to caller with buffer output as a fortran binary record  */

 return 0;

}


/*===================================================================================*/

/*  @(#) revbytes.c  McKie  Sep-1996
 *  This routine converts the data bytes in an array, high order bytes
 *  reversed with low order bytes in each word.
 *
 *  Input:
 *               a[i] = i-th byte in array of words
 *          word_size = # bytes per word in the buffer
 *            n_words = # words in array
 *
 *  Output:
 *               a[i] = i-th new byte value in array of words
 *
 *  Returns:  { nothing meaningful }
 */

int revbytes(a, word_size, n_words)
char a[];
int word_size;
int n_words;
{
 int nw;
 int ib;
 int ie;
 int iw;
 int word_size_m1;
 int word_size_half;
 int i, j;
 char t;

 /*  Return with no action if no reversal is necessary  */

 if( word_size < 2 ) { return 0; }


 /*  Compute loop constants  */

 word_size_m1 = word_size - 1; 
 word_size_half = word_size / 2;


 /*  Visit each word in the array, reversing bytes in each word  */

 for(i=0; i<n_words; i++)
 {
  ib = i * word_size;
  ie = ib + word_size_m1;
  for(j=0; j<word_size_half; j++)
  {
   t = a[ib];
   a[ib] = a[ie];
   a[ie] = t;
   ib++;
   ie--;
  } 
 }


 /*  Return to caller with array bytes converted  */

 return 0;

}


/*===================================================================================*/

/*  @(#) revbuf.c  McKie  Sep-1996
 *  This routine converts the data bytes in the current buffer, high order bytes
 *  reversed with low order bytes in each word.
 *
 *  Input:
 *          word_size = # bytes per word in the buffer
 *
 *  Returns:  { nothing meaningful }
 */

int revbuf( word_size )
int word_size;
{
 int n_words;

 /*  Compute number of words  */

 n_words = n_rec_old / word_size;


 /*  Reverse the bytes in each word for each word in the buffer  */

 revbytes(buf, word_size, n_words);


 /*  Return to caller with array bytes converted  */

 return 0;

}
