
import os
import os.path as op
import tempfile

# We import the `getoutput()` function provided by IPython.
# It allows us to do system calls from Python.
from IPython.utils.process import getoutput

def exec_cpp(code):
    """Compile, execute C++ code, and return the standard
    output."""
    # We create a temporary directory. This directory will
    # be deleted at the end of the 'with' context.
    # All created files will be in this directory.
    with tempfile.TemporaryDirectory() as tmpdir:

        # We define the source and executable filenames.
        source_path = op.join(tmpdir, 'temp.cpp')
        program_path = op.join(tmpdir, 'temp')

        # We write the code to the C++ file.
        with open(source_path, 'w') as f:
            f.write(code)

        # We compile the C++ code into an executable.
        os.system("g++ {0:s} -o {1:s}".format(
            source_path, program_path))

        # We execute the program and return the output.
        return getoutput(program_path)
"""C++ wrapper kernel."""
from ipykernel.kernelbase import Kernel

class CppKernel(Kernel):

    # Kernel information.
    implementation = 'C++'
    implementation_version = '1.0'
    language = 'c++'
    language_version = '1.0'
    language_info = {'name': 'c++',
                     'mimetype': 'text/plain'}
    banner = "C++ kernel"

    def do_execute(self, code, silent,
                   store_history=True,
                   user_expressions=None,
                   allow_stdin=False):
        """This function is called when a code cell is
        executed."""
        if not silent:
            # We run the C++ code and get the output.
            output = exec_cpp(code)

            # We send back the result to the frontend.
            stream_content = {'name': 'stdout',
                              'text': output}
            self.send_response(self.iopub_socket,
                               'stream', stream_content)
        return {'status': 'ok',
                # The base class increments the execution
                # count
                'execution_count': self.execution_count,
                'payload': [],
                'user_expressions': {},
               }

if __name__ == '__main__':
    from ipykernel.kernelapp import IPKernelApp
    IPKernelApp.launch_instance(kernel_class=CppKernel)