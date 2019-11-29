#!/usr/bin/awk -f

BEGIN { in_ignore_section = 0 }

/^log-odds matrix/ { in_ignore_section = 1 }/^$/ { in_ignore_section = 0 }

# general case
{ 
	if (in_ignore_section == 0)
		print
}
