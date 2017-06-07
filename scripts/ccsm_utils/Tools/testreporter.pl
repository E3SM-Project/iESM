#!/usr/bin/env perl
use Getopt::Long;
use Data::Dumper;
#-------------------------------------------------------------------------------
# testreporter.pl
# Perl script that watches the CESM tests as they progress, and emails the test reports
# to the specified email address. 

#-------------------------------------------------------------------------------
# Constants. Filenames we look in for test results, and the time we sleep between
# reporting test results.  
#-------------------------------------------------------------------------------
my $teststatusfilename = "TestStatus";
my $iopstatusfilename = "TestStatus.IOP";
my $casestatusfilename = "CaseStatus";
my $sleeptime = 120;
# Options and global variables.
#-------------------------------------------------------------------------------
# root of the test suite currently running. 
my $testroot = undef;
# The tag name you are testing. 
my $tagname = undef;
# the testid parameter specified for ./create_test_suite.  This script uses this parameter to 
# find the tests you are running.  
my $testid = undef;
# the email address to send test reports to. 
my $email = undef;
# the test type: one of prealpha, prebeta, or prerelease. 
my $testtype = undef;
my $debug = 0;

#-------------------------------------------------------------------------------
# Main
# Get the options first. 
# Then, get the test directories, get the test suite info, get the test status for all the tests, 
# and send the results. 
#-------------------------------------------------------------------------------
&opts;
my @testdirs;
my %suiteinfo;
while(1)
{
	@testdirs = &getTestDirs($testroot, $testid);
	%suiteinfo = &getTestSuiteInfo(\@testdirs);
	my %teststatus;
	%teststatus = getTestStatus(\@testdirs, $tagname, $testid);
	&Debug( eval { Dumper \%teststatus} );
	&Debug( eval { Dumper \%suiteinfo } );
	&sendresults(\%teststatus, \%suiteinfo);
	sleep($sleeptime);
}
#-------------------------------------------------------------------------------
# End Main
#-------------------------------------------------------------------------------
  

#-------------------------------------------------------------------------------
# Get the options. 
#-------------------------------------------------------------------------------
sub opts
{
	my $opt_help;
	GetOptions("testroot=s" => \$testroot,
			   "tagname=s" => \$tagname,
 			   "testid=s" => \$testid,
			   "email=s" => \$email,
			   "debug|d" => \$debug,	
			   "testtype=s" => \$testtype,
			   "help" => \$opt_help,
			   );
	# Show usage if the required options aren't specified. 
	&help if (defined $opt_help);
	&usage if ( (! defined $testroot) ||(! defined $tagname) || (! defined $testid) || (!defined $email) || (!defined $testtype) );
}

#-------------------------------------------------------------------------------
# Show the usage, and exit.  
#-------------------------------------------------------------------------------
sub usage
{
print <<'END';
Usage: 
./testreporter --testroot /glade/scratch/$user/tests/cesm1_1_alphaXX --tagname cesm1_1_alpha15c --testid testid --email emailaddr --testtype prealpha|prebeta|prerelease
END
exit(1);
}

sub help
{
	
print <<'END';
This is the CESM test reporter script, intended to be used to simplify the reporting of big test sets.  
Usage is as follows:
./testreporter --testroot /glade/scratch/$user/tests/cesm1_1_alphaXX 
--tagname cesm1_1_alpha15c --testid testid --email emailaddr --testtype prealpha
It gathers all the test results found in the --testroot, gets the relevant test 
status fields, and emails the results to the specified email address.
There is a companion script that should be run under 
/web/web-data/cseg/testing/tags/cesm1_1/$tag_under_test called testcatcher, which will 
read the test results from the specified email account, and automatically put 
the test results up on the cesg test page.  

Options:
--testroot This is the root of the test suite you are running. 
--tagname  The name of the tag you are testing. 
           For example, if you are testing a sandbox that will eventilally 
           become cesm1_1_alpha16c, then put cesm1_1_alpha16c. If you are 
           running a suite that will become a beta tag, put the last alpha
		   tag used.  
--testid   The testid you specified to create_test_suite.  
--testtype The type of test you are running: prealpha, prebeta, or prerelease. 
--email	   the email address you want the test results sent to.  
           
END
}



#-------------------------------------------------------------------------------
# Show debugging information if desired. 
#-------------------------------------------------------------------------------
sub Debug
{
	if($debug)
	{
			my ($msg) = @_;
			chomp $msg;
			print "Debug: $msg\n";
	}
}

#-------------------------------------------------------------------------------
# Using the testroot, find the test directories that end with the specified testid.
# If no matching directories are found, then exit. 
#-------------------------------------------------------------------------------
sub getTestDirs
{
	my ($testd, $tid) = @_;

	# Abort if the testroot does not exist.  
	if ( ! -d $testd)
	{
		print STDERR "The testroot does not exist! Aborting...\n";
		exit(1);
	}
	# open the testroot, find the test directories.  
	opendir(my $DIR, $testd) or die "can't open $testd, error was $!";
	my @testdirs = grep { $_ =~ /($tid)$/ } readdir($DIR);
	closedir $DIR;
	&Debug("in gettestdirs: test directories: \n");
	&Debug( eval { Dumper \@testdirs} );
	if(@testdirs)
	{
		return @testdirs;
	}
	else
	{
		print STDERR "It appears that the test root exists, but there aren't any test directories\n";
		print STDERR "in the testroot. Aborting...\n"; 
		exit(1);
	}
	return @testdirs;
}

#-------------------------------------------------------------------------------
# Get the suite info from config_definition, and send it back.  
#-------------------------------------------------------------------------------
sub getTestSuiteInfo
{
	my $testlist = shift;
	my %caseinfo; 
	my $firsttest = (@$testlist)[0];
	my $abspath = $testroot . "/" . $firsttest;
	&Debug("abs path:  $abspath\n");
	my @dirs = ( $abspath, $abspath . "/Tools");
	unshift @INC, @dirs;
	require XML::Lite;
	require ConfigCase;
	my $caseenv = ConfigCase->new("$abspath/Tools/config_definition.xml", "$abspath/env_case.xml");
	my $runenv = ConfigCase->new("$abspath/Tools/config_definition.xml", "$abspath/env_run.xml");
	my $buildenv = ConfigCase->new("$abspath/Tools/config_definition.xml", "$abspath/env_build.xml");
	$caseinfo{'ccsm_repotag'} = $runenv->get('CCSM_REPOTAG');
	$caseinfo{'mach'} = $caseenv->get('MACH');
	$caseinfo{'ccsmuser'} = $caseenv->get('CCSMUSER');
	$caseinfo{'compiler'} = $buildenv->get('COMPILER');
	
	&Debug("caseinfo: " . eval { Dumper \%caseinfo} );
	return %caseinfo;
}

#-------------------------------------------------------------------------------
# Get the test status. open the $testroot , look for all the test directories, 
# then get the test status. 
#-------------------------------------------------------------------------------
sub getTestStatus
{
	#my ($testd, $tag, $tid) = @_;
	my ($testdirs, $tag)  = @_;

	my %teststatushash;

	my $time = localtime;
	print "$time\n";
	
	# Iterate through each of the test directories, and get the requisite test information. 
	foreach my $testcase(@$testdirs)
	{
		# Get the test status 
	    my $statusfile = $testroot  . "/"  . $testcase .  "/" . $teststatusfilename; 
        if( ! -e $statusfile)
        {
        	warn("$statusfile does not exist, skipping to next test.");
            $teststatushash{$testcase}{'status'} = "TFAIL";
            $tetstatushash{$testcase}{'comment'} = "TestStatus file could not be found!";
            next;
        }
		&Debug( "Status file: $statusfile\n");
		open (my $teststatusfile, "<", $statusfile) or die "cannot open TestStatus file for $testcase, $!";
		my $teststatus = <$teststatusfile>;
		chomp $teststatus;
		$teststatus = (split(/\s+/, $teststatus))[0];
		&Debug("Testcase:   $testcase\n");
		&Debug( "Teststatus: $teststatus\n"); 
		$teststatushash{$testcase}{'status'} = $teststatus;

		# Now go through the TestStats getting the memleak, compare, baseline tag, throughput, and comments if any. 
		my @statuslines = <$teststatusfile>;
		my @memleaklines = grep { /memleak/ } @statuslines;
		my $memleakstatus = (split(/\s+/, $memleaklines[0]))[0];
		$teststatushash{$testcase}{'memleak'} = $memleakstatus;

		my @comparelines = grep { /compare_hist/} @statuslines;
		my ($comparestatus,$comparetest)  = split(/\s+/, $comparelines[0]);
		$teststatushash{$testcase}{'compare'} = $comparestatus;
		my $comparetag = (split(/\./, $comparetest))[-1];		
		$teststatushash{$testcase}{'baselinetag'} = $comparetag;

		my @memcomplines = grep { /memcomp/} @statuslines;
		my $memcompstatus = (split(/\s+/, $memcomplines[0]))[0];
		$teststatushash{$testcase}{'memcomp'} = $memcompstatus;

		my @tputcomplines = grep { /tputcomp/ } @statuslines;
		my $tputcompstatus = (split(/\s+/, $tputcomplines[0]))[0];
		$teststatushash{$testcase}{'tputcomp'} = $tputcompstatus;

		my @commentlines = grep { /COMMENT/ } @statuslines;
		my $comment = (split(/\s+/, $commentlines[0], 2) )[1];		
		chomp $comment;
		$teststatushash{$testcase}{'comment'} = $comment;
		
		close $teststatusfile;
			
		# Check the CaseStatus, and print out the last line...
		my $casestatusfile = $testroot . "/"  . $testcase . "/" . $casestatusfilename;
		open (my $casestatusfile, "<", $casestatusfile) or die "cannot open CaseStatusfile for $testcase, $!";

		my $lastline;
		while(<$casestatusfile>)
		{
			$lastline = $_ if eof;
		}
		close $casestatusfile;
		chomp $lastline;
		&Debug ("last line of CaseStatus: $lastline\n");
		$teststatushash{$testcase}{'casestatus'} = $lastline;

		# If the test is an IOP test, set a flag in the test status hash indicating it as such.
		if($testcase =~ /IOP\./)
		{
			$teststatushash{$testcase}{'isioptest'} = "true";
		}

		# Get the IOP test status if the file exists.   
		# If so, create separate iop* entries in the teststatus hash for this test.  
		my $iopstatusfile = $testroot . "/" . $testcase . "/" . $iopstatusfilename;
		if( -e $iopstatusfile)
		{
			open (my $iopfh, "<", $iopstatusfile) or die " cannot open IOP status file for $testcase, $!";
			my $iopstatus = <$iopfh>;
			chomp $iopstatus;
			$iopstatus = (split(/\s+/, $iopstatus))[0];
			$teststatushash{$testcase}{'iopstatus'} = $iopstatus;
			@statuslines = <$iopstatusfile>;
 			
			@memleaklines = grep { /memleak/ } @statuslines;
			$memleakstatus = (split(/\s+/, $memleaklines[0]))[0];
			$teststatushash{$testcase}{'iopmemleak'} = $memleakstatus;

			@comparelines = grep { /compare_hist/} @statuslines;
			($comparestatus,$comparetest)  = split(/\s+/, $comparelines[0]);
			$teststatushash{$testcase}{'compare'} = $comparestatus;
			$comparetag = (split(/\./, $comparetest))[-1];
			$teststatushash{$testcase}{'iopbaselinetag'} = $comparetag;

			@memcomplines = grep { /memcomp/} @statuslines;
			$memcompstatus = (split(/\s+/, $memcomplines[0]))[0];
			$teststatushash{$testcase}{'iopmemcomp'} = $memcompstatus;

			@tputcomplines = grep { /tputcomp/ } @statuslines;
			$tputcompstatus = (split(/\s+/, $tputcomplines[0]))[0];
			$teststatushash{$testcase}{'ioptputcomp'} = $tputcompstatus;

			@commentlines = grep { /COMMENT/ } @statuslines;
			$comment = (split(/\s+/, $commentlines[0], 2) )[1];
			chomp $comment;
			$teststatushash{$testcase}{'iopcomment'} = $comment;
		
			$teststatushash{$testcase}{'iopcasestatus'} = $teststatushash{$testcase}{'casestatus'};
			
			close $iopfh;
		}
			

	}
	return %teststatushash;
}

#-------------------------------------------------------------------------------
# Email the test report to the email address specified in the arguments. 
#-------------------------------------------------------------------------------
sub sendresults
{
	my ($testresults, $suiteinfo) = @_;
	my $resultsstr = "";
	&Debug ("in sendresults, email; $email\n");

	$resultsstr .= "repotag:$tagname\n"; 
	$resultsstr .= "testroot:$testroot\n";
	$resultsstr .= "mach:$$suiteinfo{'mach'}\n";
	$resultsstr .= "compiler:$$suiteinfo{'compiler'}\n";
	$resultsstr .= "testtype:$testtype\n";
	foreach my $test(sort keys %$testresults)
	{
		foreach my $detail(sort keys %{$$testresults{$test}})
		{
			$resultsstr .= "$test:$detail:$$testresults{$test}{$detail}\n";
		}
	}
	&Debug("resultsstr: \n $resultsstr\n");

	my $subject = "CESM TEST REPORT FOR $tagname $$suiteinfo{'mach'} $$suiteinfo{'compiler'}";
	open (my $mh, "|-", "mail -s '$subject' $email") or die "$!";
	print $mh $resultsstr ;
	close $mh; 
}
