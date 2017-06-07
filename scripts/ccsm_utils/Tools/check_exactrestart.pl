#!/usr/bin/env perl
$ENV{SHELL}="/bin/csh";  # Some of the comands assume csh syntax
#### compare the two log files with check_exactrestart.pl  $file1 $file2 

if($#ARGV != 1){
   print "Usage: check_exactrestart.pl  file1 file2\n";
   print "ERROR\n";
   exit;
}
## compare the two log files
   $f1=shift(@ARGV);
   $f2=shift(@ARGV);
   if($f1 eq $f2){
      print "Error: $f1 and $f2 are the same file.\n";
      print "FAIL \n";
      exit(0);
   }
   if($f2=~"\.gz"||$f1=~"\.gz") {
      print "Error: Files still gzipped\n";
      print "ERROR\n";
      exit(0);
   }
   compare($f1,$f2);
   print "log files match!\n";
   print "PASS \n";
   exit(0);

###end of the main program

#sub compare
##compare the last 500 lines of the log files.
sub compare{
  local($f1,$f2)=@_;
  local(@tmp1,@tmp2,@tt1,@tt2);
  if(!(-f $f1)){
     print" Error: file not exist: $f1\n";
     print "FAIL \n";
     exit(0);
  }
  if(!(-f $f2)){
     print " Error: file not exist: $f2\n";
     print "FAIL \n";
     exit(0);
  }
  @tmp1=`grep tStamp_write $f1|tail -5`;
  if($#tmp1<1){
     @tmp1=`grep -a tStamp_write $f1|tail -5`;
  }
  if($#tmp1<1){
     print "Error: output failed in $f1\n";
     print "FAIL \n";
     exit(0);
  }
  @tmp2=`grep tStamp_write $f2|tail -5`;
  if($#tmp2<1){
     @tmp2=`grep -a tStamp_write $f2|tail -5`;
  }
  if($#tmp2<1){
     print "Error: output failed in $f2\n";
     print "FAIL \n";
     exit(0);
  }
  @tt1=split("wall",$tmp1[$#tmp1]);
  @tt2=split("wall",$tmp2[$#tmp2]);
  if(!($tt1[0] eq $tt2[0])){ 
     print "The last date does not match. \n  $tmp1[$#tmp1]  $tmp2[$#tmp2]\n";
     print "FAIL \n";
     exit(0);
  }
  @tmp1=`grep comm_diag $f1|tail -500`;
  if($#tmp1<1){
     @tmp1=`grep -a comm_diag $f1|tail -500`;
  }
#  if($#tmp1<1){
#     print "Error: comm_diag not in output file in $f1 \n";
#     print "UNDEF\n";
#     exit(0);
#  }
  @tmp2=`grep comm_diag $f2|tail -500`;
  if($#tmp2<1){
     @tmp2=`grep -a comm_diag $f2|tail -500`;
  }
#  if($#tmp2<1){
#     print "Error: comm_diag not in output file in $f2 \n";
#     print "UNDEF\n";
#     exit(0);
#  }
  if($#tmp2!=$#tmp1){
     print "Error: $f1 and $f2 are different\n";
     print "FAIL \n";
     exit(0);
  }
  for($i=0;$i<=$#tmp1;$i++){
     if($tmp1[$i] eq $tmp2[$i]) {next;}
     else{
       print "Error: $f1 and $f2 are different.\n";
       print ">$tmp1[$i]\n";
       print "<$tmp2[$i]\n";
       print "FAIL \n";
       exit(0);
     }
  }
}##end of compare
