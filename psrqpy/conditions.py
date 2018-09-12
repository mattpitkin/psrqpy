"""
Parse logical conditions.

Adapted from https://gist.github.com/leehsueh/1290686
"""

from __future__ import print_function, division

import re


NUM_TT, VAR_TT, GT_TT, GTE_TT, LT_TT, LTE_TT, EQ_TT, NEQ_TT, LP_TT, RP_TT, AND_TT, OR_TT, NOT_TT = range(13)


class TreeNode(object):
    tokenType = None
    value = None
    left = None
    right = None
    this = None

    def __init__(self, tokenType):
        self.tokenType = tokenType


class Tokenizer(object):
    expression = None
    tokens = None
    tokenTypes = None
    i = 0

    def __init__(self, exp):
        self.expression = exp

    def next(self):
        self.i += 1
        return self.tokens[self.i-1]

    @property
    def idx(self):
        return self.i

    def peek(self):
        return self.tokens[self.i]
	
    def hasNext(self):
        return self.i < len(self.tokens)

    def nextTokenType(self):
        return self.tokenTypes[self.i]
	
    def nextTokenTypeIsOperator(self):
        t = self.tokenTypes[self.i]
        return (t == GT_TT or t == GTE_TT or t == LT_TT or t == LTE_TT or
                t == EQ_TT or t == NEQ_TT)

    def tokenize(self):
        reg = re.compile(r'(\bAND\b|\band\b|\&\&|\bOR\b|\bor\b|\|\||!=|==|<=|>=|<|>|\(|\)|\bNOT\b|\bnot\b|!|~)')
        self.tokens = reg.split(self.expression)
        self.tokens = [t.strip() for t in self.tokens if t.strip() != '']

        self.tokenTypes = []
        for t in self.tokens:
            if t.upper() == 'AND' or t == '&&':
                self.tokenTypes.append(AND_TT)
            elif t.upper() == 'OR' or t == '||':
                self.tokenTypes.append(OR_TT)
            elif t.upper() == 'NOT' or t == '!' or t == '~':
                self.tokenTypes.append(NOT_TT)
            elif t == '(':
                self.tokenTypes.append(LP_TT)
            elif t == ')':
                self.tokenTypes.append(RP_TT)
            elif t == '<':
                self.tokenTypes.append(LT_TT)
            elif t == '<=':
                self.tokenTypes.append(LTE_TT)
            elif t == '>':
                self.tokenTypes.append(GT_TT)
            elif t == '>=':
                self.tokenTypes.append(GTE_TT)
            elif t == '==':
                self.tokenTypes.append(EQ_TT)
            elif t == '!=':
                self.tokenTypes.append(NEQ_TT)
            else:
                try:
                    number = float(t)
                    self.tokenTypes.append(NUM_TT)
                except ValueError:
                    if re.search('^[a-zA-Z0-9_]+$', t):
                        self.tokenTypes.append(VAR_TT)
                    else:
                        self.tokenTypes.append(None)


class ConditionParser(object):
    tokenizer = None
    root = None

    def __init__(self, exp):
        self.tokenizer = Tokenizer(exp)
        self.tokenizer.tokenize()
        self.parse()

    def parse(self):
        self.root = self.parseExpression()

    def parseExpression(self):
        andTerm1 = self.parseAndTerm()
        while self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == OR_TT:
            self.tokenizer.next()
            andTermX = self.parseAndTerm()
            andTerm = TreeNode(OR_TT)
            andTerm.left = andTerm1
            andTerm.right = andTermX
            andTerm1 = andTerm
        return andTerm1

    def parseAndTerm(self):
        condition1 = self.parseCondition()
        while self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == AND_TT:
            self.tokenizer.next()
            conditionX = self.parseCondition()
            condition = TreeNode(AND_TT)
            condition.left = condition1
            condition.right = conditionX
            condition1 = condition
        return condition1

    def parseCondition(self):
        if self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == LP_TT:
            self.tokenizer.next()
            expression = self.parseExpression()
            if self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == RP_TT:
                self.tokenizer.next()
                return expression
            else:
                raise Exception("Closing ) expected, but got " + self.tokenizer.next())

        if self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == NOT_TT:
            condition = TreeNode(self.tokenizer.nextTokenType())
            self.tokenizer.next()
            condition.this = self.parseCondition()
            return condition

        terminal1 = self.parseTerminal()
        if self.tokenizer.hasNext() and self.tokenizer.nextTokenTypeIsOperator():
            condition = TreeNode(self.tokenizer.nextTokenType())
            self.tokenizer.next()
            terminal2 = self.parseTerminal()
            condition.left = terminal1
            condition.right = terminal2
            return condition
        else:
            raise Exception('Operator expected, but got ' + self.tokenizer.next())
	
    def parseTerminal(self):
        if self.tokenizer.hasNext():
            tokenType = self.tokenizer.nextTokenType()
            if tokenType == NUM_TT:
                n = TreeNode(tokenType)
                n.value = float(self.tokenizer.next())
                return n
            elif tokenType == VAR_TT:
                n = TreeNode(tokenType)
                n.value = self.tokenizer.next()
                return n
            else:
                raise Exception('NUM, STR, or VAR expected, but got ' + self.tokenizer.next())
        else:
            raise Exception('NUM, STR, or VAR expected, but got ' + self.tokenizer.next())
	
    def evaluate(self, table):
        return self.evaluateRecursive(self.root, table)
	
    def evaluateRecursive(self, treeNode, table):
        if treeNode.tokenType == NUM_TT:
            return treeNode.value
        if treeNode.tokenType == VAR_TT:
            if treeNode.value not in table.colnames:
                raise KeyError("Parameter '{}' not in table".format(treeNode.value))
            return table[treeNode.value]

        if treeNode.this is not None:
            this = self.evaluateRecursive(treeNode.this, table)

            if treeNode.tokenType == NOT_TT:
                return ~this

        left = self.evaluateRecursive(treeNode.left, table)
        right = self.evaluateRecursive(treeNode.right, table)
        if treeNode.tokenType == GT_TT:
            return left > right
        elif treeNode.tokenType == GTE_TT:
            return left >= right
        elif treeNode.tokenType == LT_TT:
            return left < right
        elif treeNode.tokenType == LTE_TT:
            return left <= right
        elif treeNode.tokenType == EQ_TT:
            return left == right
        elif treeNode.tokenType == NEQ_TT:
            return left != right
        elif treeNode.tokenType == AND_TT:
            return left & right
        elif treeNode.tokenType == OR_TT:
            return left | right
        else:
            raise TypeError('Unexpected type ' + str(treeNode.tokenType))
