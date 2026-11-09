"""
Zápis sdílené FF cache (forcefield_service._write_if_changed): jen při změně
a atomicky - souběžní uživatelé ji čtou (viz tests/load/).
"""

import os
import threading

import pytest

from app.services.forcefield_service import _write_if_changed

pytestmark = pytest.mark.unit


def test_unchanged_content_is_not_rewritten(tmp_path):
    path = tmp_path / "ff.rtp"
    _write_if_changed(path, b"[ atoms ]\n")
    os.utime(path, (1_000_000, 1_000_000))
    _write_if_changed(path, b"[ atoms ]\n")
    assert path.stat().st_mtime == 1_000_000


def test_changed_content_replaces_the_file_without_leftovers(tmp_path):
    path = tmp_path / "ff.rtp"
    _write_if_changed(path, b"old")
    _write_if_changed(path, b"new content")
    assert path.read_bytes() == b"new content"
    assert [p.name for p in tmp_path.iterdir()] == ["ff.rtp"]


def test_readers_never_see_a_partial_file(tmp_path):
    path = tmp_path / "ff.itp"
    versions = [bytes([65 + i]) * 200_000 for i in range(4)]
    _write_if_changed(path, versions[0])
    stop = threading.Event()
    seen = set()

    def reader():
        while not stop.is_set():
            try:
                seen.add(len(path.read_bytes()))
            except PermissionError:  # Windows: soubor se zrovna nahrazuje
                pass
            # Windows nenahradí soubor, který má někdo otevřený - čtenář bez
            # pauzy by ho držel otevřený pořád (v provozu čte jen občas).
            stop.wait(0.002)

    threads = [threading.Thread(target=reader, daemon=True) for _ in range(3)]
    for t in threads:
        t.start()
    try:
        for i in range(40):
            _write_if_changed(path, versions[i % 4])
    finally:
        stop.set()
        for t in threads:
            t.join()
    assert seen == {200_000}
