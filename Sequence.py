class Genome:
    instanceVariable = "This is an instance variable."
    def __init__(self):
        self.instanceVariable = "This is a redefined instance variable."
        methodVariable = "This is a method variable."
        print(self.instanceVariable)
        print(methodVariable)

    def findOffTargets(self):
        return ["ACTG", "TACC"]
