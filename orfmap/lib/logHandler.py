import logging
import logging.config


class Logger:
  def __init__(self, name, outpath=None):
    if outpath:
        config = get_config(outpath=outpath)
        logging.config.dictConfig(config)
    self.logger = logging.getLogger(name)

  def deco(self, log):
    return '-' * len(log)

  def info(self, log, decoration=False):
    if decoration:
      deco = self.deco(log)
      self.logger.info(deco)
      self.logger.info(log)
      self.logger.info(deco)
    else:
      self.logger.info(log)


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
