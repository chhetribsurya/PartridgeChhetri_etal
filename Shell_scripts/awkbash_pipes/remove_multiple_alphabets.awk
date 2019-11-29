#!/usr/bin/awk -f

BEGIN {
	header_printed = 0
	in_ignore_section = 0
}

/^MEME version/ {
	if ( header_printed == 1 )
		in_ignore_section = 1
}

/^ A/ {
	header_printed = 1
	if ( in_ignore_section == 1 ) {
		in_ignore_section = 0
		next
	}
}

{
	if ( in_ignore_section == 0 )
		print
}
