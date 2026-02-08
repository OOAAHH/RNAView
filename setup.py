from __future__ import annotations

from setuptools import setup


def _cmdclass():
    try:
        from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
    except Exception:  # pragma: no cover
        return {}

    class bdist_wheel(_bdist_wheel):
        def finalize_options(self) -> None:
            super().finalize_options()
            # Bundle `rnaview-hotcore` as package data => wheel must be platform-specific.
            self.root_is_pure = False

        def get_tag(self):  # type: ignore[override]
            python, abi, plat = super().get_tag()
            return "py3", "none", plat

    return {"bdist_wheel": bdist_wheel}


setup(cmdclass=_cmdclass())

