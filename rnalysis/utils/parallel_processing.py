import subprocess


def start_ipcluster(n_engines: int = 'default'):
    """
    Start an ipyparallel ipcluster in order to perform parallelized computation.

    :type n_engines: int or 'default'
    :param n_engines: if 'default', will initiate the default amount of engines. \
    Otherwise, will initiate n_engines engines.
    """
    assert (isinstance(n_engines,
                       int) and n_engines > 0) or n_engines == 'default', f"Invalid number of engines {n_engines}"
    if n_engines == 'default':
        return subprocess.Popen("ipcluster start", stderr=subprocess.PIPE, shell=True)
    else:
        return subprocess.Popen(["ipcluster", "start", "-n={:d}".format(n_engines)], stderr=subprocess.PIPE, shell=True)


def stop_ipcluster():
    """
    Stop a previously started ipyparallel ipcluster.

    """
    subprocess.Popen("ipcluster stop", stderr=subprocess.PIPE, shell=True)
