# import
## batteries
import os
import sys
import pytest

# test/data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')

def test_help_docs(script_runner):
    ret = script_runner.run('phylophlan', '-h')
    assert ret.success

