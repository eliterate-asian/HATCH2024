
def hello():
    print("Hello, there!")

def howareyou():
    print("How are you doing?")

def goodbye():
    print("See you later!")
    exit(0)

#### SWITCHER ####
    
def Choice(i, switcher):
    func = switcher.get(i, 'Invalid')
    if func != 'Invalid':
        func()

##################
        
menu = {
    "1": hello,
    "2": howareyou,
    "3": goodbye
}

while True:
    print("\n\tThis is a very cool menu yes yes\n\t1 - Hi\n\t2 - Hmm?\n\t3 - Bye")
    tabby = input("\n\t\t>>> ")
    Choice(tabby, menu)