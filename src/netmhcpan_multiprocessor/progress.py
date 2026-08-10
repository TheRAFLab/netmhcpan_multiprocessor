"""
Progress reporting, as a protocol rather than a display.

A long prediction run wants a progress bar in a terminal and complete silence
inside a pipeline. Keeping the reporting to a three-method protocol means the
library core never imports rich, and a caller can plug in a logger, a queue, or
a web socket without this package knowing.
"""

from typing import Protocol, runtime_checkable


@runtime_checkable
class ProgressReporter(Protocol):
    """Receives progress as netMHCpan's predictions stream back."""

    def add_run(self, label: str, total: int) -> object:
        """
        Register a run, and return a handle naming it in later calls.

        Args:
            label (str): a human-readable name for the run.
            total (int): how many peptides the run will predict.

        Returns:
            object: a handle, passed back to advance() and fail().
        """
        ...

    def advance(self, handle: object, amount: int = 1) -> None:
        """
        Record that a run has predicted some more peptides.

        Args:
            handle (object): the handle from add_run.
            amount (int): how many peptides have been predicted since the last call.

        Returns:
            None
        """
        ...

    def fail(self, handle: object, message: str) -> None:
        """
        Record that a run has failed and will not advance again.

        Args:
            handle (object): the handle from add_run.
            message (str): why it failed.

        Returns:
            None
        """
        ...


class GuardedProgressReporter:
    """
    Wraps a reporter so that a display failure cannot fail a prediction.

    Progress is cosmetic and a reporter is supplied by the caller: a logger, a
    queue, a web socket. None of them should be able to lose the results of a
    forty minute run by raising, so every call through here is swallowed. The
    first failure is kept in `error` for a caller that wants to know.
    """

    def __init__(self, reporter: ProgressReporter):
        self.reporter = reporter
        self.error: str | None = None

    def add_run(self, label: str, total: int) -> object:
        """
        Register a run, falling back to the label if the reporter raises.

        Args:
            label (str): a human-readable name for the run.
            total (int): how many peptides the run will predict.

        Returns:
            object: the reporter's handle, or the label if it raised.
        """
        try:
            return self.reporter.add_run(label, total)
        except Exception as reporter_error:
            self._record(reporter_error)
            return label

    def advance(self, handle: object, amount: int = 1) -> None:
        """
        Advance a run, swallowing anything the reporter raises.

        Args:
            handle (object): the handle from add_run.
            amount (int): how many peptides have been predicted since the last call.

        Returns:
            None
        """
        try:
            self.reporter.advance(handle, amount)
        except Exception as reporter_error:
            self._record(reporter_error)

    def fail(self, handle: object, message: str) -> None:
        """
        Report a failed run, swallowing anything the reporter raises.

        Args:
            handle (object): the handle from add_run.
            message (str): why it failed.

        Returns:
            None
        """
        try:
            self.reporter.fail(handle, message)
        except Exception as reporter_error:
            self._record(reporter_error)

    def _record(self, reporter_error: Exception) -> None:
        if self.error is None:
            self.error = f"{type(reporter_error).__name__}: {reporter_error}"


class NullProgressReporter:
    """The default: reports nothing, so importing this package is silent."""

    def add_run(self, label: str, total: int) -> object:
        """
        Register a run.

        Args:
            label (str): ignored.
            total (int): ignored.

        Returns:
            object: the label, so a handle is still distinguishable in a debugger.
        """
        return label

    def advance(self, handle: object, amount: int = 1) -> None:
        """
        Do nothing.

        Args:
            handle (object): ignored.
            amount (int): ignored.

        Returns:
            None
        """

    def fail(self, handle: object, message: str) -> None:
        """
        Do nothing.

        Args:
            handle (object): ignored.
            message (str): ignored.

        Returns:
            None
        """
