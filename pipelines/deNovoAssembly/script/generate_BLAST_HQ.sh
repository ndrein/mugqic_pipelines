awk ' {if ($11 < 0.00001 && $3 > 90) { print $0}} ' $1;

