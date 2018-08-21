#!/usr/bin/env python

import string
import random
import os

def password_generator(size=10, chars=string.ascii_letters + string.digits + '!@#$%^&*.()'):
    return ''.join(random.choice(chars) for i in range(size))


if __name__ == '__main__':
    # generate a strong password automatically if you are lazy
    pw = password_generator()

    # copy this password to your clipboard
    os.system("echo '%s' | pbcopy" % pw)

    print(pw)
