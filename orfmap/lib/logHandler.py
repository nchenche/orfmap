import logging
import logging.config


def get_logger(name, outpath=None):

    if outpath:
        config = get_config(outpath=outpath)
        logging.config.dictConfig(config)

    return logging.getLogger(name)


def get_config(outpath='./'):
    config = {
      "version": 1,
      "disable_existing_loggers": False,
      "formatters": {
        "detailed": {
          "format": "%(asctime)s %(name)s %(levelname)s:\t%(message)s",
          "datefmt": "%H:%M:%S"
        },
        "simple": {
          "format": "%(message)s"
        }
      },

      "handlers": {
        "console": {
          "class": "logging.StreamHandler",
          "level": "DEBUG",
          "formatter": "detailed"
        },

        "info_handler": {
          "class": "logging.FileHandler",
          "level": "DEBUG",
          "formatter": "detailed",
          "filename": outpath+"info.log",
          "mode": "w+",
          "encoding": "utf8",
          "delay": False
        },

        "summary_handler": {
          "class": "logging.FileHandler",
          "level": "INFO",
          "formatter": "simple",
          "filename": outpath+"summary.log",
          "mode": "w+",
          "encoding": "utf8",
          "delay": False
        }
      },

      "loggers": {
        "orfmap.lib": {
          "level": "DEBUG",
          "handlers": ["console", "info_handler"],
          "propagate": False
        }
      },

      "root": {
        "level": "NOTSET",
        "handlers": ["console", "info_handler"],
        "propagate": False
      }
    }

    return config
