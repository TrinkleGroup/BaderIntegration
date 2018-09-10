(NR==1){x=$1; y=$2;}
{if ($1!=x) {
  print x,y;
  x=$1; y=$2;
 } else {
  if ($2>y) y=$2;
 }
}
END{print x,y}
