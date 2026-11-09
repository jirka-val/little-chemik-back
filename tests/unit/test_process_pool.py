"""
Pool výpočetních procesů (app/core/process_pool.py) se skutečnými procesy:
výsledek, výjimky (i ty, které nejde přenést celé) a zprávy konzole.
"""

import logging
import os

import pytest

from app.core.exceptions import AppBaseException, BadRequestError
from app.core.process_pool import HeavyProcessPool

pytestmark = pytest.mark.unit


def square_with_pid(x):
    logging.getLogger("app.console").info(f"computing {x}")
    return x * x, os.getpid()


def fail_with_code():
    raise BadRequestError("bad input", code="some_code")


class _Unpicklable(AppBaseException):
    status_code = 409
    code = "missing_dof"

    def __init__(self, handle):           # jako ForgeMissingDOFError: konstruktor s objekty
        super().__init__("needs a decision", payload={"chain": "A"})
        self.handle = handle


def fail_unpicklable():
    raise _Unpicklable(lambda: None)       # lambda nejde přes pickle


@pytest.fixture(scope="module")
def pool():
    p = HeavyProcessPool(workers=1)
    yield p
    p.shutdown()


async def test_result_comes_from_another_process(pool):
    value, pid = await pool.run(square_with_pid, 7)
    assert value == 49 and pid != os.getpid()


async def test_console_messages_are_forwarded(pool):
    # vlastní handler místo caplog - check_all.sh spouští pytest s -p no:logging
    messages = []
    handler = logging.Handler(level=logging.INFO)
    handler.emit = lambda record: messages.append(record.getMessage())
    console = logging.getLogger("app.console")
    console.addHandler(handler)
    try:
        await pool.run(square_with_pid, 3)
    finally:
        console.removeHandler(handler)
    assert "computing 3" in messages


async def test_app_error_keeps_its_code(pool):
    with pytest.raises(BadRequestError) as exc:
        await pool.run(fail_with_code)
    assert exc.value.code == "some_code" and exc.value.message == "bad input"


async def test_unpicklable_app_error_keeps_status_code_and_payload(pool):
    with pytest.raises(AppBaseException) as exc:
        await pool.run(fail_unpicklable)
    assert (exc.value.status_code, exc.value.code, exc.value.payload) == (409, "missing_dof", {"chain": "A"})


async def test_disabled_pool_runs_in_a_thread():
    value, pid = await HeavyProcessPool(workers=0).run(square_with_pid, 2)
    assert value == 4 and pid == os.getpid()
