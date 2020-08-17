import subprocess
import time


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


def start_parallel_session(n_engines: int = 'default'):
    """
    Stop previous ipyparallel ipcluster and start a new one in order to perform parallelized computation.

    :type n_engines: int or 'default'
    :param n_engines: if 'default', will initiate the default amount of engines. \
    Otherwise, will initiate n_engines engines.

    :Examples:
    >>> from rnalysis import general
    >>> general.start_parallel_session()
    Starting parallel session...
    Parallel session started successfully
    """
    print("Starting parallel session...")
    start_time = time.time()
    try:
        stop_ipcluster()
    except FileNotFoundError:
        pass
    time.sleep(0.5)
    stream = start_ipcluster(n_engines)
    while True:
        line = stream.stderr.readline().decode('utf8')
        if 'Engines appear to have started successfully' in line:
            break
        elif 'Cluster is already running' in line:
            stop_ipcluster()
            time.sleep(0.5)
            stream = start_ipcluster(n_engines)
        elif line != '':
            print(line.replace('\n', ''))

        if time.time() - start_time > 180:
            raise ConnectionError("Failed to start parallel session for over 3 minutes.")
    print('Parallel session started successfully')
