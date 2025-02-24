import os
import typing
import inspect
import functools


PathStr = typing.Union[os.PathLike, str]


def is_type(instance, type_info):
    if type_info == typing.Any:
        return True
    if hasattr(type_info, "__origin__"):
        if type_info.__origin__ in {typing.Union, typing.Optional}:
            return any(is_type(instance, arg) for arg in type_info.__args__)
    return isinstance(instance, type_info)


def enforce_types(func):
    types_info = typing.get_type_hints(func)
    signature = inspect.signature(func)

    @functools.wraps(func)
    def inner(*args, **kwargs):
        sig = signature.bind(*args, **kwargs)
        sig.apply_defaults()
        for name, value in sig.arguments.items():
            expected_type = types_info.get(name, typing.Any)
            if not is_type(value, expected_type):
                raise TypeError(
                    f"{name} is not of the correct type: {type(value).__name__}. "
                    f"Expected {expected_type.__name__}."
                )
        return func(*args, **kwargs)

    return inner
