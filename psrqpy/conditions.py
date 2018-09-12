"""
Parser for logical conditions.

Adapted from https://gist.github.com/leehsueh/1290686
"""

from __future__ import print_function, division

import re
import numpy as np
from .config import *

NUM_TT, VAR_TT, GT_TT, GTE_TT, LT_TT, LTE_TT, EQ_TT, NEQ_TT, LP_TT, RP_TT, AND_TT, OR_TT, NOT_TT, ASSOC_TT, TYPE_TT = range(15)

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
        return (t == GT_TT or t == GTE_TT or t == LT_TT or t == LTE_TT or
                t == EQ_TT or t == NEQ_TT)

    def nextTokenTypeIsGroup(self):
        t = self.tokenTypes[self.i]
        return (t == NOT_TT or self.nextTokenTypeIsMatch())

    def nextTokenTypeIsMatch(self):
        t = self.tokenTypes[self.i]
        return (t == ASSOC_TT or t == TYPE_TT)

    def tokenize(self):
        reg = re.compile(LOGEXPRS)
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
            elif t.upper() == 'ASSOC':
                self.tokenTypes.append(ASSOC_TT)
            elif t.upper() == 'TYPE':
                self.tokenTypes.append(TYPE_TT)
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

    def __init__(self, expression):
        self.tokenizer = Tokenizer(expression)
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
        tokenType = self.tokenizer.nextTokenType()

        if self.tokenizer.hasNext() and tokenType == LP_TT:
            self.tokenizer.next()
            expression = self.parseExpression()
            if self.tokenizer.hasNext() and self.tokenizer.nextTokenType() == RP_TT:
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
                if condition.this.tokenType != VAR_TT:
                    raise Exception("Problem passing 'matchType' expression")

            return condition

        terminal1 = self.parseTerminal()
        if not self.tokenizer.hasNext() or self.tokenizer.nextTokenType() == RP_TT:
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
	
    def evaluate(self, table, exactMatch=False):
        return self.evaluateRecursive(self.root, table, exactMatch=exactMatch)
	
    def evaluateRecursive(self, treeNode, table, exactMatch):
        if treeNode.tokenType == NUM_TT:
            return treeNode.value
        if treeNode.tokenType == VAR_TT:
            if treeNode.value not in table.colnames:
                raise KeyError("Parameter '{}' not in table".format(treeNode.value))
            return table[treeNode.value]

        if treeNode.this is not None:
            if treeNode.isGroup:
                if treeNode.tokenType == NOT_TT:
                    this = self.evaluateRecursive(treeNode.this, table, exactMatch)
                    return ~this.astype(np.bool)  # make sure its a boolean
                elif treeNode.isMatch:
                    return self.evaluateMatch(treeNode, table, exactMatch)
                else:
                    raise TypeError("Unexpected type '{}'".format(str(treeNode.tokenType)))

        left = self.evaluateRecursive(treeNode.left, table, exactMatch)
        right = self.evaluateRecursive(treeNode.right, table, exactMatch)
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
            raise TypeError("Unexpected type '{}'".format(str(treeNode.tokenType)))

    def evaluateMatch(self, treeNode, table, exactMatch):
        if treeNode.tokenType == ASSOC_TT:
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
        elif treeNode.tokenType == TYPE_TT:
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
