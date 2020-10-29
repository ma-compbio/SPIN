#!/bin/bash

bin=$1
state=$2
output=$3

#merge states and write in bed format
paste $bin $state|awk '{
	if($1==lchrom && $5==lname && $2 == lend)
		{
			lend = $3
		}else{
			if(lchrom) {
				print lchrom"\t"lstart"\t"lend"\t"lname;
			}; 
			lchrom=$1; lstart=$2; lend=$3; lname=$5
		}
}
END{
	print lchrom"\t"lstart"\t"lend"\t"lname
}' > $output















