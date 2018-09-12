"""
Parser for logical conditions.

Adapted from https://gist.github.com/leehsueh/1290686
"""

from __future__ import print_function, division

import re
import numpy as np
from .config import *



# string of logical expressions for use in regex parser
LOGEXPRS = (r'(\bAND\b'     # logical AND
            r'|\band\b'     # logical AND
            r'|\&\&'        # logical AND
            r'|\bOR\b'      # logical OR
            r'|\bor\b'      # logical OR
            r'|\|\|'        # logical OR
            r'|!='          # not equal to
            r'|=='          # equal to
            r'|<='          # less than or equal to
            r'|>='          # greater than or equal to
            r'|<'           # less than
            r'|>'           # greater than
            r'|\('          # left opening bracket
            r'|\)'          # right closing bracket
            r'|\bNOT\b'     # logical NOT
            r'|\bnot\b'     # logical NOT
            r'|!'           # logical NOT
            r'|~'           # logical NOT
            r'|\bASSOC\b'   # pulsar association
            r'|\bassoc\b'   # pulsar association
            r'|\bTYPE\b'    # pulsar type
            r'|\btype\b)')  # pulsar type


class TokenType(object):
    """
    Object containing an enumerated list of token types parsed from the
    expression string.
    """

    NUM = 0       # a numeric value
    VAR = 1       # a named variable
    GT = 2        # greater than '>'
    GTE = 3       # greater than or equal to '>='
    LT = 4        # less than '<'
    LTE = 5       # less than or equal to '<='
    EQ = 6        # equal to '=='
    NEQ = 7       # not equal to '!='
    LP = 8        # left bracket '('
    RP = 9        # right bracket ')'
    AND = 10      # logical AND
    OR = 11       # logical OR
    NOT = 12      # logical NOT
    ASSOC = 13    # pulsar association
    TYPE = 14     # pulsar type


class TreeNode(object):
    tokenType = None
    value = None
    left = None
    right = None
    this = None
    isOp = False
    isGroup = False
    isMatch = False

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
        return (t == TokenType.GT or t == TokenType.GTE or t == TokenType.LT or
                t == TokenType.LTE or t == TokenType.EQ or t == TokenType.NEQ)

    def nextTokenTypeIsGroup(self):
        t = self.tokenTypes[self.i]
        return (t == TokenType.NOT or self.nextTokenTypeIsMatch())

    def nextTokenTypeIsMatch(self):
        t = self.tokenTypes[self.i]
        return (t == TokenType.ASSOC or t == TokenType.TYPE)

    def tokenize(self):
        reg = re.compile(LOGEXPRS)
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
            elif t.upper() == 'ASSOC':
                self.tokenTypes.append(TokenType.ASSOC)
            elif t.upper() == 'TYPE':
                self.tokenTypes.append(TokenType.TYPE)
            else:
                try:
                    number = float(t)
                    self.tokenTypes.append(TokenType.NUM)
                except ValueError:
                    if re.search('^[a-zA-Z0-9_]+$', t):
                        self.tokenTypes.append(TokenType.VAR)
                    else:
                        self.tokenTypes.append(None)


class ConditionParser(object):
    tokenizer = None
    root = None

    def __init__(self, expression):
        self.tokenizer = Tokenizer(expression)
        self.tokenizer.tokenize()
        self.parse()

    def parse(self):
        self.root = self.parseExpression()

    def parseExpression(self):
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
        tokenType = self.tokenizer.nextTokenType()

        if self.tokenizer.hasNext() and tokenType == TokenType.LP:
            self.tokenizer.next()
            expression = self.parseExpression()
            if self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == TokenType.RP:
                self.tokenizer.next()
                return expression
            else:
                raise Exception("Closing ')' expected, but got '{}'".format(self.tokenizer.next()))

        if self.tokenizer.hasNext() and self.tokenizer.nextTokenTypeIsGroup():
            matchType = self.tokenizer.nextTokenTypeIsMatch()
            condition = TreeNode(tokenType)
            self.tokenizer.next()
            condition.this = self.parseCondition()
            condition.isGroup = True

            if matchType:
                condition.isMatch = True
                # a "matching" expression, e.g., an association
                if condition.this.tokenType != TokenType.VAR:
                    raise Exception("Problem passing 'matchType' expression")

            return condition

        terminal1 = self.parseTerminal()
        if not self.tokenizer.hasNext() or self.tokenizer.nextTokenType() == TokenType.RP:
            return terminal1
        elif self.tokenizer.hasNext() and self.tokenizer.nextTokenTypeIsOperator():
            condition = TreeNode(self.tokenizer.nextTokenType())
            self.tokenizer.next()
            terminal2 = self.parseTerminal()
            condition.left = terminal1
            condition.right = terminal2
            condition.isOp = True
            return condition
        else:
            raise Exception("Operator expected, but got '{}'".format(self.tokenizer.next()))

    def parseTerminal(self):
        if self.tokenizer.hasNext():
            tokenType = self.tokenizer.nextTokenType()
            if tokenType == TokenType.NUM:
                n = TreeNode(tokenType)
                n.value = float(self.tokenizer.next())
                return n
            elif tokenType == TokenType.VAR:
                n = TreeNode(tokenType)
                n.value = self.tokenizer.next()
                return n
            else:
                raise Exception('NUM, STR, or VAR expected, but got ' + self.tokenizer.next())
        else:
            raise Exception('NUM, STR, or VAR expected, but got ' + self.tokenizer.next())
	
    def evaluate(self, table, exactMatch=False):
        return self.evaluateRecursive(self.root, table, exactMatch=exactMatch)
	
    def evaluateRecursive(self, treeNode, table, exactMatch):
        if treeNode.tokenType == TokenType.NUM:
            return treeNode.value
        if treeNode.tokenType == TokenType.VAR:
            if treeNode.value not in table.colnames:
                raise KeyError("Parameter '{}' not in table".format(treeNode.value))
            return table[treeNode.value]

        if treeNode.this is not None:
            if treeNode.isGroup:
                if treeNode.tokenType == TokenType.NOT:
                    this = self.evaluateRecursive(treeNode.this, table, exactMatch)
                    return ~this.astype(np.bool)  # make sure its a boolean
                elif treeNode.isMatch:
                    return self.evaluateMatch(treeNode, table, exactMatch)
                else:
                    raise TypeError("Unexpected type '{}'".format(str(treeNode.tokenType)))

        left = self.evaluateRecursive(treeNode.left, table, exactMatch)
        right = self.evaluateRecursive(treeNode.right, table, exactMatch)
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
        else:
            raise TypeError("Unexpected type '{}'".format(str(treeNode.tokenType)))

    def evaluateMatch(self, treeNode, table, exactMatch):
        if treeNode.tokenType == TokenType.ASSOC:
            if 'ASSOC' in table.colnames:
                assocs = table['ASSOC']

                if exactMatch:
                    # only return exact matches
                    return assocs == treeNode.this.value
                else:
                    # return True if containing the value
                    return np.array([treeNode.this.value in a for a in assocs])
            else:
                raise KeyError("No 'ASSOC' in table")
        elif treeNode.tokenType == TokenType.TYPE:
            if 'TYPE' in table.colnames:
                types = table['TYPE']

                if treeNode.this.value not in PSR_TYPES:
                    raise KeyError("Type '{}' is not a valid type".format(treeNode.this.value))

                if exactMatch:
                    # only return exact matches
                    return types == treeNode.this.value
                else:
                    # return True if containing the value
                    return np.array([treeNode.this.value in t for t in types])
            else:
                raise KeyError("No 'TYPE' in table")
        else:
            raise TypeError("Unexpected type '{}'".format(str(treeNode.tokenType)))
