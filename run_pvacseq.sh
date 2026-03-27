#!/bin/bash
python3 -c "
import sys
import pvactools.tools.pvacseq.run as run
run.main(sys.argv[1:])
" "$@"
