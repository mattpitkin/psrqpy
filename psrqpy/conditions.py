"""
Parse logical conditions.

Adapted from https://gist.github.com/leehsueh/1290686
"""

from __future__ import print_function, division

import re

class TokenType(object):
    NUM, STR, VAR, GT, GTE, LT, LTE, EQ, NEQ, LP, RP, AND, OR, NOT = range(14)


class TreeNode(object):
    tokenType = None
    value = None
    left = None
    right = None

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
        return (t == TokenType.GT or t == TokenType.GTE
                or t == TokenType.LT or t == TokenType.LTE
                or t == TokenType.EQ or t == TokenType.NEQ)

    def tokenize(self):
        reg = re.compile(r'(\bAND\b|\band\b|\&\&|\bOR\b|\bor\b|\|\||!=|==|<=|>=|<|>|\(|\)|\bNOT\b|\bnot\b|!|~)')
        self.tokens = reg.split(self.expression)
        self.tokens = [t.strip() for t in self.tokens if t.strip() != '']

        self.tokenTypes = []
        for t in self.tokens:
            if t.upper() == 'AND' or t == '&&':
                self.tokenTypes.append(TokenType.AND)
            elif t.upper() == 'OR' or t == '||':
                self.tokenTypes.append(TokenType.OR)
            elif t.upper() == 'NOT' or t == '!' or t == '~':
                self.tokenTypes.append(TokenType.NOT)
            elif t == '(':
                self.tokenTypes.append(TokenType.LP)
            elif t == ')':
                self.tokenTypes.append(TokenType.RP)
            elif t == '<':
                self.tokenTypes.append(TokenType.LT)
            elif t == '<=':
                self.tokenTypes.append(TokenType.LTE)
            elif t == '>':
                self.tokenTypes.append(TokenType.GT)
            elif t == '>=':
                self.tokenTypes.append(TokenType.GTE)
            elif t == '==':
                self.tokenTypes.append(TokenType.EQ)
            elif t == '!=':
                self.tokenTypes.append(TokenType.NEQ)
            else:
                # number of string or variable
                if t[0] == '"' and t[-1] == '"':
                    self.tokenTypes.append(TokenType.STR)
                else:
                    try:
                        number = float(t)
                        self.tokenTypes.append(TokenType.NUM)
                    except ValueError:
                        if re.search('^[a-zA-Z_]+$', t):
                            self.tokenTypes.append(TokenType.VAR)
                        else:
                            self.tokenTypes.append(None)


class BooleanParser(object):
    tokenizer = None
    root = None

    def __init__(self, exp):
        self.tokenizer = Tokenizer(exp)
        self.tokenizer.tokenize()
        self.parse()

    def parse(self):
        self.root = self.parseExpression()

    def parseExpression(self):
        if self.tokenizer.idx == 0 and self.tokenizer.nextTokenType() == TokenType.NOT:
            orTerm1 = None
        else:
            orTerm1 = self.parseOrTerm()

        while self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == TokenType.NOT:
            self.tokenizer.next()
            orTermX = self.parseOrTerm()
            orTerm = TreeNode(TokenType.NOT)
            orTerm.left = orTerm1
            orTerm.right = orTermX
            orTerm1 = orTerm

        if orTerm1 is None:
            raise Exception("Error passing expression")

        return orTerm1

    def parseOrTerm(self):
        andTerm1 = self.parseAndTerm()
        while self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == TokenType.OR:
            self.tokenizer.next()
            andTermX = self.parseAndTerm()
            andTerm = TreeNode(TokenType.OR)
            andTerm.left = andTerm1
            andTerm.right = andTermX
            andTerm1 = andTerm
        return andTerm1

    def parseAndTerm(self):
        condition1 = self.parseCondition()
        while self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == TokenType.AND:
            self.tokenizer.next()
            conditionX = self.parseCondition()
            condition = TreeNode(TokenType.AND)
            condition.left = condition1
            condition.right = conditionX
            condition1 = condition
        return condition1

    def parseCondition(self):
        if self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == TokenType.LP:
            self.tokenizer.next()
            expression = self.parseExpression()
            if self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == TokenType.RP:
                self.tokenizer.next()
                return expression
            else:
                raise Exception("Closing ) expected, but got " + self.tokenizer.next())

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
            if tokenType == TokenType.NUM:
                n = TreeNode(tokenType)
                n.value = float(self.tokenizer.next())
                return n
            elif tokenType == TokenType.STR or tokenType == TokenType.VAR:
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
        if treeNode.tokenType == TokenType.NUM or treeNode.tokenType == TokenType.STR:
            return treeNode.value
        if treeNode.tokenType == TokenType.VAR:
            return table[treeNode.value]

        right = self.evaluateRecursive(treeNode.right, table)
        if treeNode.tokenType == TokenType.NOT and treeNode.left is None:
            return ~right

        left = self.evaluateRecursive(treeNode.left, table)
        if treeNode.tokenType == TokenType.GT:
            return left > right
        elif treeNode.tokenType == TokenType.GTE:
            return left >= right
        elif treeNode.tokenType == TokenType.LT:
            return left < right
        elif treeNode.tokenType == TokenType.LTE:
            return left <= right
        elif treeNode.tokenType == TokenType.EQ:
            return left == right
        elif treeNode.tokenType == TokenType.NEQ:
            return left != right
        elif treeNode.tokenType == TokenType.AND:
            return left & right
        elif treeNode.tokenType == TokenType.OR:
            return left | right
        elif treeNode.tokenType == TokenType.NOT:
            return ~right
        else:
            raise Exception('Unexpected type ' + str(treeNode.tokenType))
