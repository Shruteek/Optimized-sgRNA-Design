class MetaGenome:
    instanceVariable = "This is an instance variable."
    def __init__(self):
        self.instanceVariable = "This is a redefined instance variable."
        methodVariable = "This is a method variable."
        print("This is the instantiation method.")
        print(self.instanceVariable)
        print(methodVariable)

    def method(self):
        print("This is a method.")
