def get_lines(source):
    """
    Reads lines from a given source, which can be a file path or a file-like object.

    Args:
        source (str or file-like object): The input source. It can be a string representing
            a file path or a file-like object (binary or text).

    Yields:
        str: Lines from the source as strings. Generator

    Raises:
        TypeError: If the input source is neither a string nor a file-like object.
    """

    # Accepts file path or file-like object, yields lines as str
    if isinstance(source, str):
        with open(source, 'r') as f:
            for line in f:
                yield line

    
    elif isinstance(source, (bytes, bytearray, memoryview)):
        # Raw binary data or memoryview
        if isinstance(source, memoryview):
            source = source.tobytes()
        text = source.decode('utf-8')
        for line in text.splitlines(keepends=True):
            yield line

    elif hasattr(source, 'read'):
        # file-like object, could be binary or text
        # If binary, decode to str
        first = source.read(0)
        if hasattr(source, 'encoding'):
            # TextIO
            source.seek(0)
            for line in source:
                yield line
        else:
            # BinaryIO
            source.seek(0)
            for line in source:
                yield line.decode('utf-8')
    else:
        raise TypeError("Unsupported input type. Must be file path or file-like object.")
