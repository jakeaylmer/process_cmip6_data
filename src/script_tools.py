"""Tools for scripts, mainly parsing command-line arguments
using argparse.
"""

from argparse import ArgumentParser

from . import metadata as md


def default_cmd_argument_parser(**kwargs):
    """Default command line arguments -- returns the argparse
    ArgumentParser instance so that ad-hoc arguments can be
    added. Keyword arguments are passed to ArgumentParser().
    """
    
    prsr = ArgumentParser(**kwargs)
    
    prsr.add_argument("-m", "--model", type=str,
                      default=md.default_model_id,
                      choices=md.defined_models,
                      help="Model ID (case sensitive)")
    
    prsr.add_argument("-x", "--experiment", type=str,
                      default=md.default_experiment_id,
                      choices=md.defined_experiments,
                      help="Experiment ID (case sensitive)")
    
    return prsr



def default_cmd_args(**kwargs):
    """Default command line arguments -- returns the parsed
    arguments directly. Keyword arguments are passed to
    default_cmd_argument_parser().
    """
    return default_cmd_argument_parser(**kwargs).parse_args()



def qc_argument_parser(diagnostic_choices=["none"], **kwargs):
    """Command line arguments parser with additional arguments
    for quality control (qc) plots. Parameter diagnostic_choices
    is a list of valid string arguments for the flag "-d" or
    "--diagnostic" (default is ["none"]). Keyword arguments are
    passed to default_cmd_argument_parser(). Returns the
    ArgumentParser instance so that ad-hoc arguments can be
    added.
    """
    
    prsr = default_cmd_argument_parser(**kwargs)
    
    prsr.add_argument("-d", "--diagnostic", type=str,
                      choices=diagnostic_choices,
                      default=diagnostic_choices[0],
                      help="Diagnostic")
    
    prsr.add_argument("--seldiag", type=int, default=0,
                      help="Diagnostic alias select if multiple"
                           + " available for --model (and "
                           + "--diagnostic, if applicable)")
    
    prsr.add_argument("--x2", type=str, default="none",
                      choices=["none"] + md.defined_experiments,
                      help="Second experiment (only valid with "
                           + "--experiment historical)")
    
    prsr.add_argument("-e", "--ensmemberid", type=int,
                      default=-1,
                      help="Realisation ID for member (or 0 "
                           + "for whichever is first, -1 for "
                           + "ensemble mean, -2 to show all)")
    
    prsr.add_argument("-l", "--latitude", type=float,
                      default=65.0,
                      help="Degrees north/south latitude to "
                           + "evaluate at for profiles or "
                           + "area integrals/means")
    
    prsr.add_argument("--movavg", type=int, default=21,
                      help="Moving average filter width (years)")
    
    prsr.add_argument("--yravgs", type=int, default=21,
                      help="Averaging periods for piControl "
                           + "(years)")
    
    prsr.add_argument("--yravg", type=int, nargs=2,
                      default=(1, 21),
                      help="Averaging period (year start, "
                           + "year end)")
    
    prsr.add_argument("--yravg1", type=int, nargs=2,
                      default=(1980, 2000),
                      help="Average period 1 (year start, "
                           + "year end)")
    
    prsr.add_argument("--yravg2", type=int, nargs=2,
                      default=(2001, 2021),
                      help="Average period 2 (year start, "
                           + "year end)")
    
    prsr.add_argument("-s", "--savefigs", action="store_true",
                      help="Save figures rather than show")
    
    return prsr



def _check_qc_experiment_args(cmd):
    """Check that two provided experiments are compatible for
    concatenation, raising a ValueError if not. Also adjusts the
    default year averaging arguments as appropriate.
    """
    
    x1 = cmd.experiment
    x2 = cmd.x2
    
    if x2 != "none":
        if x1 == "piControl":
            raise ValueError("Cannot combine experiment \'"
                + x2 + "\' with \'" + x1 + "\'")
        elif x1 == "historical" and "ssp" not in x2:
            raise ValueError("Cannot combine experiment \'"
                + x2 + "\' with \'" + x1 + "\'")
    
    if x1 != "piControl" and cmd.yravg == (1, 21):
        cmd.yravg = (1850, 1870)



def parse_qc_cmd_args(diagnostic_choices=["none"], **kwargs):
    """Command line arguments parser with additional arguments
    for quality control (qc) plots. Keyword arguments are passed
    to qc_argument_parser(). Returns the parsed args.
    """
    
    prsr = qc_argument_parser(diagnostic_choices, **kwargs)
    cmd = prsr.parse_args()
    _check_qc_experiment_args(cmd)
    
    return cmd
