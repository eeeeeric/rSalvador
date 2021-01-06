# Dec 28, 2020


def prune (x, big=5000):

   if (x<big):

      return x

   else:

      return big

def is_number(s):

   try:

      float(s)
      return True

   except ValueError:
      return False


 


