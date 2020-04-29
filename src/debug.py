import debug_globals as dbgg
def printname(func, module='', color='', print_this=True):
    # ---------------------------------------------
    # this wrapper is for bind func info (like __name__) to wrapped local function, instead of the wrapper itself
    # @functools.wraps(func)
    # ---------------------------------------------
    def wrapper(*args, **kwargs):
        if dbgg.print_this and print_this:
            if len(color) > 0:
                  print(color, end='')  
            # ---------------------------------------------
            # --- print starting char, and put it to '' ---
            if dbgg.previous_mod == module and dbgg.previous_func == func:
                dbgg.starting_char = ''
                print('.', end = '')
            else:
                if dbgg.starting_char == '':
                    print('\n', end = '')
                print("@", module, "@ §", func.__name__, "§")
                dbgg.starting_char = '\n'
            dbgg.previous_mod = module
            dbgg.previous_func = func
            if len(color) > 0:
                    print(dbgg.RESET, end='')
        dbgg.print_this = True
        return func(*args, **kwargs)
    return wrapper

def printempty(func):
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper