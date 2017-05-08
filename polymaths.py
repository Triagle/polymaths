
class OrderError(Exception):
    """ An exception raised if any expressions that an operator is
    operating on have unequal orders. """

    def __init__(self, order_a, order_b):
        self.order_a = order_a
        self.order_b = order_b

    def __str__(self):
        return 'Order {} does not match order {}'.format(self.order_a,
                                                         self.order_b)


class Term():
    """ Represents a single term in a polynomial. """
    def __init__(self, order, coefficient):
        """ Build term from order and coefficient. """
        self.order = order
        self.coefficient = coefficient
    def add(self, other, more=[]):
        """ Add two (or more) terms together. """
        # Guard against unequal orders
        well_ordered = any(self.order != other.order for other in more)
        if self.order != other.order or not well_ordered:
            raise OrderError(self.order, other.order)
        rhs = sum(term.coefficient for term in more + [other])
        coefficient = self.coefficient + rhs
        return Term(self.order, coefficient)

    def negate(self):
        """ Return the negated variant of self. """
        return Term(self, self.coefficient * -1)

    def sub(self, other, more=[]):
        """ Subtract two (or more) terms. """
        return self.add(other.negate(), more=more)

    def divide(self, other):
        """ Divide two terms. """
        self.order -= other.order
        self.coefficient /= other.coefficient

    def multiply(self, other):
        """ Multiply two terms. """
        self.order += other.order
        self.coefficient *= other.coefficient

    def empty():
        """ Return an empty term. """
        return Term(0, 0)

    def is_zero(self):
        """ Returns true if self is a zero term. """
        return self == Term.empty()


class Polynomial():
    """ Base polynomial to factorize. """
    def __init__(self, terms):
        """ Build a polynomial on supplied terms. """
        self.terms = terms

    def order(self):
        return max(term.order for term in self.terms)

    def map_to_order(self, order):
        """ Return a polynomial mapped to order. Works only if desired
        order is higher than the current polynomials order. """
        if self.order() == order:
            return self
        elif self.order() < order:
            # Pad zeros to the end of the order
            pad = (order - self.order()) * [Term.empty()]
            Polynomial(order, pad + self.terms)
        else:
            return None

    def _bucket_terms(self):
        """ Return a dictionary (or bucket) where each key associates
        an order with some terms. """
        terms_dict = {}
        for term in self.terms:
            term_order = term.order
            if term_order not in terms_dict:
                terms_dict[term_order] = [term]
            else:
                terms_dict[term_order].append(term)
        return terms_dict

    def normalized_application(self, fn, other):
        """ Apply fn to the normalized representation of self and other. """

        to_order = max(self.order(), other.order())
        lhs = self.map_to_order(to_order)
        rhs = other.map_to_order(to_order)
        return [fn(a, b) for a, b in zip(lhs.terms,
                                         rhs.terms)]

    def subtract(self, other):
        """ Return the result of subtracting two polynomials """
        subtract_terms = self.normalized_application(Term.sub,
                                                     other)
        return Polynomial(subtract_terms)

    def add(self, other):
        """ Return the result of adding two polynomials """
        add_terms = self.normalized_application(Term.add,
                                                other)
        return Polynomial(add_terms)

    def collect_terms(self):
        """ Return a polynomial where all like terms are collected. """
        if len(self.terms) <= 1:
            return self
        else:
            term_map = self._bucket_terms()
            term_results = []
            for terms in term_map.values():
                if len(terms) == 1:
                    term_results.append(terms[0])
                else:
                    term_results.append(terms[0].add(terms[1:]))
                    term_results.sort(key=lambda term: term.order, reverse=True)
            return Polynomial(term_results)

    def multiply(self, other):
        """ multiply two polynomials together. """
        dot_product = [Term.multiply(i, j) for (i, j) in
                       zip(self.terms, other.terms)]
        polynomial = Polynomial(dot_product)
        return polynomial.collect_terms()


    def index_of_order(self, order):
        """ Return the index of the variable at order in the
        coefficient array. """
        if order > self.order():
            return None
        else:
            return self.order() - order

    def zero_polynomial(self):
        """ Returns True if polynomial is zero. """
        return all(term.is_zero() for term in self.terms)

    def term_at(self, order):
        """ Return term at order. """
        index = self.index_of_order(order)
        return self.terms[index]
