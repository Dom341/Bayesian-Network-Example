__author__ = 'joshkang'

"""Node class. Holds necessary data types for Neapolitan's propagation algorithm.
Uniquely identified by single char, nodeChar. The event is described by the string, event.
Probability is stored as a list and the probability of the node's event occurring
is described as the syntax D = node, D0 = probability of D not occurring, D1 = probability of D occurring
assuming events are binary probability variables
"""
class Node(object):

    """initializes node: sets event name, its node char, and its given probability P.
    Sets all other values by default to -1 and parent as None
    """
    def __init__(self, event, nodeChar, prob):
        self._event = event
        self._nodeChar = nodeChar
        self._children = []

        #for probabilities of child node, the 4 probabilities are stored
        #in a list of 2 lists where P(c0|d0) = prob[0][0] and P(c1|d0) = prob[1][0]
        self._prob = prob
        self._belief = [-1.0, -1.0]
        #for pi, lambda, pi message, and lambda message, syntax is like pi(d0,d1)
        self._pValue = [-1.0, -1.0]
        self._lValue = [-1.0, -1.0]
        self._pMsg = [-1.0, -1.0]
        self._lMsg = [-1.0, -1.0]
        self._parent = None

    @property
    def event(self):
        return self._event

    @property
    def nodeChar(self):
        return self._nodeChar

    @property
    def children(self):
        return self._children

    @property
    def prob(self):
        return self._prob

    @property
    def belief(self):
        return self._belief
    @belief.setter
    def belief(self, belief):
        self._belief = belief

    @property
    def pValue(self):
        return self._pValue
    @pValue.setter
    def pValue(self, pValue):
        self._pValue = pValue

    @property
    def lValue(self):
        return self._lValue
    @lValue.setter
    def lValue(self, lValue):
        self._lValue = lValue

    @property
    def pMsg(self):
        return self._pMsg
    @pMsg.setter
    def pMsg(self, pMsg):
        self._pMsg = pMsg

    @property
    def lMsg(self):
        return self._lMsg
    @lMsg.setter
    def lMsg(self, lMsg):
        self._lMsg = lMsg

    @property
    def parent(self):
        return self._parent
    @parent.setter
    def parent(self, parent):
        self._parent = parent


    def add_child(self, child):
        self._children.append(child)
        child.parent = self


    def __str__(self):

        beliefStr = "P'(%s0) = %.3f \nP'(%s1) = %.3f" %(self.nodeChar, self.belief[0], self.nodeChar, self.belief[1])
        probStr = None

        #is root node, thus has 2 probabilities
        if self.parent == None:
            probStr = "P(%s0) = %.3f \nP(%s1) = %.3f" %(self.nodeChar, self.prob[0], self.nodeChar, self.prob[1])

        #is not root node, thus has 4 probabilities because each child has only one parent
        else:
            probStr = "P(%s0|%s0) = %.3f \nP(%s0|%s1) = %.3f" \
                      "\nP(%s1|%s0) = %.3f \nP(%s1|%s1) = %.3f" \
                      %(self.nodeChar, self.parent.nodeChar, self.prob[0][0], self.nodeChar, self.parent.nodeChar, self.prob[0][1],
                        self.nodeChar, self.parent.nodeChar, self.prob[1][0], self.nodeChar, self.parent.nodeChar, self.prob[1][1])

        return "Node: %s \nEvent: %s \n%s \n%s" \
               "\nPi value: (%.3f, %.3f) \nLambda value: (%.3f, %.3f)  \n" \
               "Pi message: (%.3f, %.3f) \nLambda message: (%.3f, %.3f)" \
               %(self.nodeChar, self.event, probStr, beliefStr, self.pValue[0], self.pValue[1],
                 self.lValue[0], self.lValue[1], self.pMsg[0], self.pMsg[1], self.lMsg[0], self.lMsg[1])

"""Checks tree for given node to instantiate, handles 'C' and 'c' as the same as all nodes should be unique (one letter) """
def searchTree(node, nodeChar, result):

    if node.nodeChar == nodeChar or node.nodeChar.lower() == nodeChar:
        updateI(node, result)
        return

    else:
        for child in node.children:
            return searchTree(child, nodeChar, result)

        print("Node not found")


"""Prints the state of all nodes in the tree"""
def printTree(node):
    print(node)
    print("\n")#

    for child in node.children:
        printTree(child)

"""Sets all lambda messages and lambda values of tree to 1"""
def initSetL(node):
    node.lValue = [1.0, 1.0]
    node.lMsg = [1.0, 1.0]

    for child in node.children:
        initSetL(child)

"""Initialization of tree according to Neapolitan's probability propagation algorithm"""
def initialize(node):

    if node.parent == None:

        initSetL(node)

        for i in range(len(node.prob)):
            node.belief[i] = node.prob[i]
            node.pValue[i] = node.prob[i]


        for child in node.children:
            updateP(child, node.belief)

    else:
        print("Can only initialize from root of a tree")

"""Starts propagation of pi and lambda messages from instantiated node
If result = 0, the event did not occur. If result = 1, the event did occur
"""
def updateI(node, result):
    if result == 0:
        node.belief[0] = 1
        node.belief[1] = 0
    elif result == 1:
        node.belief[0] = 0
        node.belief[1] = 1

    calcLambda(node)

    if node.parent != None:
        updateL(node, node.parent)

    for child in node.children:
        updateP(child, node.belief)



"""Recursive formula - propagates pi and lambda messages throughout tree"""
def updateL(child, node):

    calcLambdaMsg(child)

    #checks if node is already instantiated
    if node.belief[0] != 0 and node.belief[1] != 0:
        calcLambda(node)
        calcNewP(node)

        if node.parent != None:
            updateL(node, node.parent)

        #call updateP for children
        for currChild in node.children:
            if currChild != child:
                updateP(currChild, node.belief)


"""Recursive formula - propagates pi messages throughout tree"""
def updateP(node, pBelief):

    calcPiMsg(node, pBelief)

    #checks if node is already instantiated
    if node.belief[0] != 0 and node.belief[1] != 0:
        calcPi(node)
        calcNewP(node)

        for child in node.children:
            updateP(child, node.belief)

"""Operative Formula 1. Calculates the lambda message from Child B to Parent A.
Posterior evidence vector- stores in Child B
"""
def calcLambdaMsg(node):
    node.lMsg[0] = node.prob[0][0] * node.lValue[0] + node.prob[1][0] * node.lValue[1]

    node.lMsg[1] = node.prob[0][1] * node.lValue[0] + node.prob[1][1] * node.lValue[1]

"""Operative Formula 2. Calculates the pi message from Parent A to Child B.
Prior evidence vector- stores in Child B
"""
def calcPiMsg(node, pBelief):
    if pBelief[0] == 1:
        node.pMsg[0] = 1
        node.pMsg[1] = 0
    elif pBelief[1] == 1:
        node.pMsg[0] = 0
        node.pMsg[1] = 1
    else:
        node.pMsg[0] = pBelief[0] / node.lValue[0]
        node.pMsg[1] = pBelief[1] / node.lValue[1]

"""Operative Formula 3. Calculates the lambda value of a node"""
def calcLambda(node):
    if node.belief[0] == 1:
        node.lValue[0] = 1
        node.lValue[1] = 0
    elif node.belief[1] == 1:
        node.lValue[0] = 0
        node.lValue[1] = 1
    else:
        for child in node.children:
            node.lValue[0] *= child.lMsg[0]
            node.lValue[1] *= child.lMsg[1]


"""Operative Formula 4. Calculates the pi value of a node"""
def calcPi(node):
    node.pValue[0] = node.prob[0][0] * node.pMsg[0] + node.prob[0][1] * node.pMsg[1]

    node.pValue[1] = node.prob[1][0] * node.pMsg[0] + node.prob[1][1] * node.pMsg[1]

"""Operative Formula 5. Calculates the conditional probability of a node, P'
Calculates the belief of a node, how likely the probability of it not occurring/occurring
given information from parent and children
"""
def calcNewP(node):
    temp = node.lValue[0] * node.pValue[0]
    temp1 = node.lValue[1] * node.pValue[1]

    #normalization
    if (temp + temp1) != 1:
        total = temp + temp1
        alpha = 1 / total
        temp *= alpha
        temp1 *= alpha

    node.belief[0] = temp
    node.belief[1] = temp1

"""Examples to run Neapolitan's algorithm for probability propagation in trees
Drug study example from class
"""
def example1():
    #setting up tree
    drug = Node("Being in drug study", 'D', [0.9, 0.1])
    cured = Node("Patient cured", 'C', [[0.5, 0.25],[0.5, 0.75]])
    drug.add_child(cured)

    print("___EXAMPLE 1___")
    print("**Before initialization** \n")
    printTree(drug)

    initialize(drug)
    print("**After initialization** \n")
    printTree(drug)

    #Instantiates node C as occurring, aka doctor encounters cured patient
    #sets P': C1 = 1 C0 = 0

    searchTree(drug, 'C', 1)
    print("**After instantiation** \n")
    printTree(drug)



"""Cheating spouse example from Neapolitan's"""
def example2():
    #setting up tree
    spouseA = Node("Spouse cheating", 'A', [0.9, 0.1])
    spouseB = Node("Spouse dining with another", 'B', [[0.8, 0.3], [0.2, 0.7]])
    spouseC = Node("Spouse reported seen dining with another", 'C', [[0.999, 0.6], [0.001, 0.4]])
    spouseD = Node("Strange man/lady calls on the phone", 'D', [[0.6, 0.2], [0.4, 0.8]])

    spouseA.add_child(spouseB)
    spouseA.add_child(spouseD)
    spouseB.add_child(spouseC)

    print("___EXAMPLE 2___")
    print("**Before initialization** \n")
    printTree(spouseA)

    initialize(spouseA)
    print("**After initialization** \n")
    printTree(spouseA)

    #Instantiates node C as occurring, aka spouse is reported seen dining with another
    #sets P': C1 = 1 C0 = 0

    searchTree(spouseA, 'C', 1)##
    print("**After instantiation** \n")
    printTree(spouseA)


def __main__():
    example1()
    example2()


if __name__ == '__main__':__main__()
