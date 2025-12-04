import logging
import logging.handlers
import multiprocessing
from multiprocessing import Queue

# Module-level logger
_logger = None
_queue = None
_listener = None


def get_logger():
    """Get the CNS logger instance."""
    global _logger
    if _logger is None:
        _logger = logging.getLogger("cns")
        if not _logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter("%(levelname)s: %(message)s")
            handler.setFormatter(formatter)
            _logger.addHandler(handler)
            _logger.setLevel(logging.WARNING)  # Default to WARNING
    return _logger


def set_verbose(verbose):
    """Set the logging level based on verbosity flag.
    
    Parameters
    ----------
    verbose : bool
        If True, sets logging level to INFO. Otherwise, sets to WARNING.
    """
    logger = get_logger()
    logger.setLevel(logging.INFO if verbose else logging.WARNING)


def setup_mp_logging():
    """Setup logging for multiprocessing.
    
    Returns a queue that worker processes can use to send log records
    to the main process.
    
    Returns
    -------
    multiprocessing.Queue
        Queue for sending log records from workers to main process.
    """
    global _queue, _listener
    if _queue is None:
        _queue = multiprocessing.Manager().Queue(-1)
        handler = logging.StreamHandler()
        formatter = logging.Formatter("%(levelname)s: %(message)s")
        handler.setFormatter(formatter)
        _listener = logging.handlers.QueueListener(_queue, handler, respect_handler_level=True)
        _listener.start()
    return _queue


def stop_mp_logging():
    """Stop the multiprocessing logging listener."""
    global _listener, _queue
    if _listener is not None:
        _listener.stop()
        _listener = None
        _queue = None


def configure_worker_logging(queue, verbose):
    """Configure logging in a worker process.
    
    Parameters
    ----------
    queue : multiprocessing.Queue
        Queue for sending log records to main process.
    verbose : bool
        If True, sets logging level to INFO.
    """
    global _logger
    _logger = logging.getLogger("cns")
    # Remove any existing handlers
    _logger.handlers.clear()
    if queue is not None:
        handler = logging.handlers.QueueHandler(queue)
        _logger.addHandler(handler)
    else:
        handler = logging.StreamHandler()
        formatter = logging.Formatter("%(levelname)s: %(message)s")
        handler.setFormatter(formatter)
        _logger.addHandler(handler)
    _logger.setLevel(logging.INFO if verbose else logging.WARNING)


def log_info(text):
    """Log an info message.
    
    The message will only be displayed if verbose mode is enabled via set_verbose(True).
    
    Parameters
    ----------
    text : str
        The message to log.
    """
    get_logger().info(text)


def log_warn(text):
    """Log a warning message.
    
    Parameters
    ----------
    text : str
        The warning message to log.
    """
    get_logger().warning(text)


def log_error(text):
    """Log an error message.
    
    Parameters
    ----------
    text : str
        The error message to log.
    """
    get_logger().error(text)


def suppress_errors(suppress=True):
    """Suppress error messages from being displayed.
    
    Parameters
    ----------
    suppress : bool
        If True, sets logging level to CRITICAL (suppressing ERROR and below).
        If False, restores to WARNING level.
    """
    logger = get_logger()
    logger.setLevel(logging.CRITICAL if suppress else logging.WARNING)