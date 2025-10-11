#!/usr/bin/env python3
"""
Fix gene sequence headers to match metadata accessions.
Extracts the accession number from headers like 'NC_012920.1:14747-15887'
"""

import sys
import re

for line in sys.stdin:
    if line.startswith('>'):
        # Extract accession from headers like '>NC_012920.1:14747-15887 MT-CYB ...'
        match = re.match(r'>([^:]+)(:[0-9-]+)?\s', line)
        if match:
            accession = match.group(1)
            # Replace the header with just the accession
            line = f'>{accession}\n'
    sys.stdout.write(line)