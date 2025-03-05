die unless $#ARGV=2;
$min = $ARGV[0];
$max = $ARGV[1];
#my int ($linenum :64bit, $filterout :64bit, $totalreads :64bit, $filterednum :64bit, $outputnum :64bit);
$linenum=0;
$filterout=0;
$totalreads=0;
$filterednum=0;
$outputnum=0;
while(<STDIN>){
  $line = $_; chomp($line);
  $record[$linenum]=$line;
  if ($linenum==1){
    $filterout = (length($line)<$min) || (length($line)>$max);
    if ($filterout) {
      $filterednum=$filterednum+1;
    }else{
      $outputnum=$outputnum+1;
    }
    $totalreads=$totalreads+1;
  }
  if ($linenum==3 && (!$filterout)){
    for ($cx=0;$cx<4;$cx++){
      $lineinrecord = $record[$cx];
      print($lineinrecord."\n");
    }
  }
  $linenum=($linenum+1)%4;
}
printf STDERR "$totalreads processed.\n$filterednum reads filtered out.\n$outputnum reads accepted.\n";
