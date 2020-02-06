import sys

from classes.Printer import Printer


def main():
    p = Printer()
    printer = {
        "times_and_obj": p.print_times_and_obj,
        "polygons": p.print_polygons,
        "iterations": p.print_iterations,
        "comparison": p.print_comparison,
        "calibration": p.print_calibration
    }
    printer[sys.argv[1]]()


if __name__ == "__main__":
    main()
