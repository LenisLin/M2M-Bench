from __future__ import annotations

import logging


def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """
    Create or return a configured logger.
    - No duplicate handlers on repeated calls.
    - Formatter: [TIME] [LEVEL] name: message
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        handler = logging.StreamHandler()
        formatter = logging.Formatter("[%(asctime)s] [%(levelname)s] %(name)s: %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger
