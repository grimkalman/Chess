# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 10:35:01 2020

@author: Grim Kalman
"""

"""
Open pre-calculated masks to avoid unnecessary calculations and facilitate
table lookup instead of on the spot calculations for speed. The attack masks
was calculated using the o-2s trick and the idea to store the masks in
dictionaries within a dictionary was taken from:
"Avoiding rotated bitboards with direct lookup" (S. Tannous, 2007).
"""

with open("file_mask_table.txt", "r") as file:
    fmt = eval(file.read())
    file.close()
with open("rank_mask_table.txt", "r") as file:
    rmt = eval(file.read())
    file.close()
with open("diagonal_mask_table.txt", "r") as file:
    dmt = eval(file.read())
    file.close()
with open("anti_diagonal_mask_table.txt", "r") as file:
    admt = eval(file.read())
    file.close()
with open("bitboard_to_pos.txt", "r") as file:
    bitboard_to_pos = eval(file.read())
    file.close()
with open("file_attacks.txt", "r") as file:
    fa = eval(file.read())
    file.close()
with open("rank_attacks.txt", "r") as file:
    ra = eval(file.read())
    file.close()
with open("diagonal_attacks.txt", "r") as file:
    da = eval(file.read())
    file.close()
with open("anti_diagonal_attacks.txt", "r") as file:
    ada = eval(file.read())
    file.close()
with open("king_attacks.txt", "r") as file:
    king_attacks = eval(file.read())
    file.close()
with open("knight_attacks.txt", "r") as file:
    knight_attacks = eval(file.read())
    file.close()
with open("white_pawn_attacks.txt", "r") as file:
    white_pawn_attacks = eval(file.read())
    file.close()
with open("black_pawn_attacks.txt", "r") as file:
    black_pawn_attacks = eval(file.read())
    file.close()
with open("piece_value_table.txt", "r") as file:
    pst = eval(file.read())
    file.close()

pieces = ["P", "R", "N", "B", "Q", "K", "p", "r", "n", "b", "q", "k"]

rows = ["1", "2", "3", "4", "5", "6", "7", "8"]

cols = ["a", "b", "c", "d", "e", "f", "g", "h"]


class Game_State:
    """Class that represents the state of the chessboard."""

    def __init__(self):
        self.white_to_move = True
        self.bitboards = [0xFF00,               # 0 : white pawns
                          0x81,                 # 1 : white rooks
                          0x42,                 # 2 : white knights
                          0x24,                 # 3 : white bishops
                          0x10,                 # 4 : white queen
                          0x8,                  # 5 : white king
                          0xFF000000000000,     # 6 : black pawns
                          0x8100000000000000,   # 7 : black rooks
                          0x4200000000000000,   # 8 : black knights
                          0x2400000000000000,   # 9 : black bishops
                          0x1000000000000000,  # 10 : black queen
                          0x800000000000000,   # 11 : black king
                          0x000000000000FFFF,  # 12 : white pieces
                          0xFFFF000000000000,  # 13 : black pieces
                          0xFFFF00000000FFFF,  # 14 : all pieces
                          0,                   # 15 : en passant
                          0x8900000000000089   # 16 : castle rights
                          ]
        self.score = 0
        self.tp = {}
        self.history = []

    def load_FEN(self, string):
        """Load a position from a FEN-string."""
        self.bitboards = [0 for bitboards in self.bitboards]
        fen = string.split(" ")
        pos = fen[0].split("/")

        for row, rank in enumerate(pos):
            char = list(rank)
            col = 0
            for item in char:
                if item.isnumeric():
                    col += int(item)
                else:
                    self.bitboards[pieces.index(item)] += (2**(8*(7 - row) +
                                                               (7 - col)))
                    col += 1

        self.bitboards[12] = (self.bitboards[0] | self.bitboards[1] |
                              self.bitboards[2] | self.bitboards[3] |
                              self.bitboards[4] | self.bitboards[5])
        self.bitboards[13] = (self.bitboards[6] | self.bitboards[7] |
                              self.bitboards[8] | self.bitboards[9] |
                              self.bitboards[10] | self.bitboards[11])
        self.bitboards[14] = self.bitboards[12] | self.bitboards[13]
        self.bitboards[16] |= 0x9 if "K" in list(fen[2]) else 0
        self.bitboards[16] |= 0x88 if "Q" in list(fen[2]) else 0
        self.bitboards[16] |= 0x900000000000000 if "k" in list(fen[2]) else 0
        self.bitboards[16] |= 0x8800000000000000 if "q" in list(fen[2]) else 0

        if fen[1] == "w":
            self.white_to_move = True
            if fen[3] == '-':
                self.bitboards[15] = 0
            else:
                self.bitboards[15] = (2**(8*(rows.index(list(fen[3])[1]) - 1) +
                                          (7 - cols.index(list(fen[3])[0]))))
        else:
            self.white_to_move = False
            if fen[3] == '-':
                self.bitboards[15] = 0
            else:
                self.bitboards[15] = (2**(8*(rows.index(list(fen[3])[1]) + 1) +
                                          (7 - cols.index(list(fen[3])[0]))))

    def print_board(self):
        """Print a ASCII representation of the board in the console."""
        chessBoard = ["." for i in range(64)]
        for i in range(12):
            board = "{:064b}".format(self.bitboards[i])
            for y in range(8):
                for x in range(8):
                    if board[8 * y + x] == "1":
                        chessBoard[8 * y + x] = pieces[i]
        for i in range(8):
            print(str(8 - i), " ".join(chessBoard[8*i:8*(i + 1)]))
        print("  a b c d e f g h")

    def generate_moves(self):
        """Generate all pseudo-legal moves."""
        if self.white_to_move:
            return (king_moves(self, self.bitboards[14], self.bitboards[16],
                               self.bitboards[12], self.bitboards[5]) +
                    queen_moves(self.bitboards[14], self.bitboards[12],
                                self.bitboards[4]) +
                    rook_moves(self.bitboards[14], self.bitboards[12],
                               self.bitboards[1]) +
                    bishop_moves(self.bitboards[14], self.bitboards[12],
                                 self.bitboards[3]) +
                    knight_moves(self.bitboards[12], self.bitboards[2]) +
                    white_pawn_moves(self.bitboards[15], self.bitboards[13],
                                     self.bitboards[14], self.bitboards[0]))
        else:
            return (king_moves(self, self.bitboards[14], self.bitboards[16],
                               self.bitboards[13], self.bitboards[11]) +
                    queen_moves(self.bitboards[14], self.bitboards[13],
                                self.bitboards[10]) +
                    rook_moves(self.bitboards[14], self.bitboards[13],
                               self.bitboards[7]) +
                    bishop_moves(self.bitboards[14], self.bitboards[13],
                                 self.bitboards[9]) +
                    knight_moves(self.bitboards[13], self.bitboards[8]) +
                    black_pawn_moves(self.bitboards[15], self.bitboards[12],
                                     self.bitboards[14], self.bitboards[6]))

    def move(self, move):
        """Update the game state given a move.

        :param move: string representing a move, eg. 'e2e4'.
        """
        self.tp.clear()
        from_col, from_row, to_col, to_row = list(move)
        from_sq = 2**((7 - cols.index(from_col)) + 8*rows.index(from_row))
        to_sq = 2**((7 - cols.index(to_col)) + 8*rows.index(to_row))

        if self.white_to_move:
            for i in range(6):
                if from_sq & self.bitboards[i] != 0:
                    piece = i
                    break
        else:
            for i in range(6, 12):
                if from_sq & self.bitboards[i] != 0:
                    piece = i - 6
                    break

        self.make_move((from_sq, to_sq, piece))
        self.print_board()
        self.make_move(self.alpha_beta(6, -float('Inf'), float('Inf'))[1])
        print()
        self.print_board()

    def make_move(self, move):
        """Update the game state given a move.

        :param move: (from, to, piece) tuple. Eg. a pawn move from e2 to e4
                      would be written as: (2048, 134217728, 0).
        """
        self.history.append((self.white_to_move, self.bitboards[:],
                             self.score))
        move_from, move_to, piece = move

        if self.white_to_move:
            self.bitboards[piece] += move_to - move_from
            self.bitboards[12] += move_to - move_from
            self.score += (pst.get(piece).get(move_to) -
                           pst.get(piece).get(move_from))
            if (self.bitboards[13] & move_to):
                for i in range(6, 12):
                    if (self.bitboards[i] & move_to):
                        self.bitboards[i] -= move_to
                        self.bitboards[13] -= move_to
                        self.bitboards[15] = 0
                        self.score += pst.get(i).get(move_to)
                        break
            if piece == 0:
                if (move_from << 16) == move_to:
                    self.bitboards[15] = move_to
                elif move_to == self.bitboards[15]:
                    self.bitboards[6] -= move_to >> 8
                    self.bitboards[13] -= move_to >> 8
                    self.bitboards[15] = 0
                elif move_to & 0xff00000000000000 != 0:
                    self.bitboards[0] -= move_to
                    self.bitboards[4] += move_to
                else:
                    self.bitboards[15] = 0
            elif piece == 5:
                if move_from >> 2 == move_to:
                    self.bitboards[1] += (move_from >> 1) - (move_to >> 1)
                    self.bitboards[12] += (move_from >> 1) - (move_to >> 1)
                elif move_from << 2 == move_to:
                    self.bitboards[1] += (move_from << 1) - (move_to << 2)
                    self.bitboards[12] += (move_from << 1) - (move_to << 2)
                else:
                    self.bitboards[15] = 0
            else:
                self.bitboards[15] = 0
            self.white_to_move = False
        else:
            self.bitboards[piece + 6] += move_to - move_from
            self.bitboards[13] += move_to - move_from
            self.score -= (pst.get(piece + 6).get(move_to) -
                           pst.get(piece + 6).get(move_from))
            if (self.bitboards[12] & move_to):
                for i in range(6):
                    if (self.bitboards[i] & move_to):
                        self.bitboards[i] -= move_to
                        self.bitboards[12] -= move_to
                        self.bitboards[15] = 0
                        self.score -= pst.get(i).get(move_to)
                        break
            if piece == 0:
                if (move_from >> 16) == move_to:
                    self.bitboards[15] = move_to
                elif move_to == self.bitboards[15]:
                    self.bitboards[0] -= move_to << 8
                    self.bitboards[12] -= move_to << 8
                    self.bitboards[15] = 0
                elif move_to & 0xff != 0:
                    self.bitboards[6] -= move_to
                    self.bitboards[10] += move_to
                else:
                    self.bitboards[15] = 0
            elif piece == 5:
                if move_from >> 2 == move_to:
                    self.bitboards[7] += (move_from >> 1) - (move_to >> 1)
                    self.bitboards[13] += (move_from >> 1) - (move_to >> 1)
                elif move_from << 2 == move_to:
                    self.bitboards[7] += (move_from << 1) - (move_to << 2)
                    self.bitboards[13] += (move_from >> 1) - (move_to >> 1)
                else:
                    self.bitboards[15] = 0
            else:
                self.bitboards[15] = 0
            self.white_to_move = True

        self.bitboards[14] = self.bitboards[12] | self.bitboards[13]
        self.bitboards[16] = ((move_from & self.bitboards[16]) ^
                              self.bitboards[16])
        self.bitboards[16] = ((move_to & self.bitboards[16]) ^
                              self.bitboards[16])

    def undo_move(self):
        """Undo the previous move."""
        self.white_to_move, self.bitboards, self.score = self.history.pop()

    def is_check(self):
        """Determine if the king is in check.

        :return: bitboard containing all attacking pieces or 1 if the king is
                 taken.
        """
        if self.white_to_move and self.bitboards[5] != 0:
            (col, row), (diag, adiag) = bitboard_to_pos.get(self.bitboards[11])
            return ((self.bitboards[4] | self.bitboards[3]) & (da.get(self.bitboards[11]).get(self.bitboards[14] & dmt[diag]) | ada.get(self.bitboards[11]).get(self.bitboards[14] & admt[adiag]))) | ((self.bitboards[4] | self.bitboards[1]) & (ra.get(self.bitboards[11]).get(self.bitboards[14] & rmt[row]) | fa.get(self.bitboards[11]).get(self.bitboards[14] & fmt[col]))) | (self.bitboards[2] & knight_attacks.get(self.bitboards[11])) | (self.bitboards[0] & black_pawn_attacks.get(self.bitboards[11])) | (self.bitboards[5] & king_attacks.get(self.bitboards[11]))
        elif not self.white_to_move and self.bitboards[11] != 0:
            (col, row), (diag, adiag) = bitboard_to_pos.get(self.bitboards[5])
            return ((self.bitboards[10] | self.bitboards[9]) & (da.get(self.bitboards[5]).get(self.bitboards[14] & dmt[diag]) | ada.get(self.bitboards[5]).get(self.bitboards[14] & admt[adiag]))) | ((self.bitboards[10] | self.bitboards[7]) & (ra.get(self.bitboards[5]).get(self.bitboards[14] & rmt[row]) | fa.get(self.bitboards[5]).get(self.bitboards[14] & fmt[col]))) | (self.bitboards[8] & knight_attacks.get(self.bitboards[5])) | (self.bitboards[6] & white_pawn_attacks.get(self.bitboards[5])) | (self.bitboards[11] & king_attacks.get(self.bitboards[5]))
        else:
            return 1

    def perft(self, depth):
        """Traverse the game tree, mainly for debugging purposes."""
        if depth > 0:
            nodes = 0
            for move in self.generate_moves():
                self.make_move(move)
                if not self.is_check():
                    nodes += self.perft(depth - 1)
                self.undo_move()
            return nodes
        else:
            return 1

    def alpha_beta(self, depth, alpha, beta):
        """Preform a alpha-beta search to the desired depth."""
        if depth == 0:
            return self.score, None
        else:
            best_move = None
            if self.white_to_move:
                move_list = self.sort(self.generate_moves(), -1)
                for move in move_list:
                    self.make_move(move[1])
                    if not self.is_check():
                        score = self.alpha_beta(depth - 1, alpha, beta)[0]
                        if score > alpha:  # white maximizes her score
                            alpha = score
                            best_move = move[1]
                            self.undo_move()
                            if alpha >= beta:  # alpha-beta cutoff
                                break
                        else:
                            self.undo_move()
                    else:
                        self.undo_move()
                return (alpha, best_move)
            else:
                move_list = self.sort(self.generate_moves(), 1)
                for move in move_list:
                    self.make_move(move[1])
                    if not self.is_check():
                        score = self.alpha_beta(depth - 1, alpha, beta)[0]
                        if score < beta:  # black minimizes his score
                            beta = score
                            best_move = move[1]
                            self.undo_move()
                            if alpha >= beta:  # alpha-beta cutoff
                                break
                        else:
                            self.undo_move()
                    else:
                        self.undo_move()
                return (beta, best_move)

    def sort(self, move_list, turn):
        sorted_list = []
        for move in move_list:
            self.make_move(move)
            sorted_list.append((self.score, move))
            self.undo_move()
        return sorted(sorted_list, key=lambda item: turn*item[0])

    def generate_hash(self):
        """Generate a hash-key from the position."""
        return hash(tuple(self.bitboards)) + int(self.white_to_move)


def is_attacked(gs, s):
    """
    Determine if the square is attacked by enemy pieces.

    :param gs: a Game_State class object.
    :param s: bitboard with the square to be checked set to 1.
    :return: bitboard containing all attacking pieces.
    """
    (col, row), (diag, adiag) = bitboard_to_pos.get(s)

    if gs.white_to_move:
        return (((gs.bitboards[10] | gs.bitboards[9]) &
                (da.get(s).get(gs.bitboards[14] & dmt[diag]) |
                 ada.get(s).get(gs.bitboards[14] & admt[adiag]))) |
                ((gs.bitboards[10] | gs.bitboards[7]) &
                 (ra.get(s).get(gs.bitboards[14] & rmt[row]) |
                  fa.get(s).get(gs.bitboards[14] & fmt[col]))) |
                (gs.bitboards[8] & knight_attacks.get(s)) |
                (gs.bitboards[6] & white_pawn_attacks.get(s)) |
                (gs.bitboards[11] & king_attacks.get(s)))
    else:
        return (((gs.bitboards[4] | gs.bitboards[3]) &
                 (da.get(s).get(gs.bitboards[14] & dmt[diag]) |
                  ada.get(s).get(gs.bitboards[14] & admt[adiag]))) |
                ((gs.bitboards[4] | gs.bitboards[1]) &
                 (ra.get(s).get(gs.bitboards[14] & rmt[row]) |
                  fa.get(s).get(gs.bitboards[14] & fmt[col]))) |
                (gs.bitboards[2] & knight_attacks.get(s)) |
                (gs.bitboards[0] & black_pawn_attacks.get(s)) |
                (gs.bitboards[5] & king_attacks.get(s)))


def white_pawn_moves(ep, bp, o, s):
    """
    Generate a bitboard of available pawn moves for white.

    :param ep: bitboard containing the available en passant captures.
    :param bp: bitboard containing all black pieces.
    :param o: bitboard containing all pieces.
    :param s: bitboard of all white pawns.
    :return: list of pseudo-legal moves.
    """
    move_list = []
    clippedL = s & 0x7F7F7F7F7F7F7F7F
    clippedR = s & 0xFEFEFEFEFEFEFEFE
    push = (s << 8) & (o ^ 0xFFFFFFFFFFFFFFFF)
    double_push = ((push & 0x0000000000FF0000) << 8) & (o ^ 0xFFFFFFFFFFFFFFFF)
    captures = (clippedL << 9 & bp | clippedR << 7 & bp |
                (clippedL << 1 & ep) << 8 | (clippedR >> 1 & ep) << 8)

    move_to = push & -push
    while move_to > 0:
        move_list.append((move_to >> 8, move_to, 0))
        push = push & (push - 1)
        move_to = push & -push

    move_to = double_push & -double_push
    while move_to > 0:
        move_list.append((move_to >> 16, move_to, 0))
        double_push = double_push & (double_push - 1)
        move_to = double_push & -double_push

    move_to = captures & -captures
    while move_to > 0:
        source = s & black_pawn_attacks.get(move_to)
        move_from = source & -source
        while move_from > 0:
            move_list.append((move_from, move_to, 0))
            source = source & (source - 1)
            move_from = source & -source
        captures = captures & (captures - 1)
        move_to = captures & -captures

    return move_list


def black_pawn_moves(ep, wp, o, s):
    """
    Generate a bitboard of available pawn moves for black.

    :param ep: bitboard containing the available en passant captures.
    :param wp: bitboard containing all white pieces.
    :param o: bitboard containing all pieces.
    :param s: bitboard of all white pawns.
    :return: list of pseudo-legal moves.
    """
    move_list = []
    clippedL = s & 0x7F7F7F7F7F7F7F7F
    clippedR = s & 0xFEFEFEFEFEFEFEFE
    push = (s >> 8) & (o ^ 0xFFFFFFFFFFFFFFFF)
    double_push = ((push & 0xFF0000000000) >> 8) & (o ^ 0xFFFFFFFFFFFFFFFF)
    captures = (clippedR >> 9 & wp | clippedL >> 7 & wp |
                (clippedL << 1 & ep) >> 8 | (clippedR >> 1 & ep) >> 8)

    move_to = push & -push
    while move_to > 0:
        move_list.append((move_to << 8, move_to, 0))
        push = push & (push - 1)
        move_to = push & -push

    move_to = double_push & -double_push
    while move_to > 0:
        move_list.append((move_to << 16, move_to, 0))
        double_push = double_push & (double_push - 1)
        move_to = double_push & -double_push

    move_to = captures & -captures
    while move_to > 0:
        source = s & white_pawn_attacks.get(move_to)
        move_from = source & -source
        while move_from > 0:
            move_list.append((move_from, move_to, 0))
            source = source & (source - 1)
            move_from = source & -source
        captures = captures & (captures - 1)
        move_to = captures & -captures

    return move_list


def king_moves(gs, o, c, b, s):
    """
    Generate a bitboard of available king moves.

    :param gs: a Game_State class object.
    :param o: bitboard containing all pieces.
    :param c: biboard containing castle relevant pieces that has not moved.
    :param b: bitboard of all pieces of own color.
    :param s: bitboard of the king
    :return: list of pseudo-legal moves.
    """
    move_list = []
    moves = king_attacks.get(s) & (b ^ 0xFFFFFFFFFFFFFFFF)

    move_to = moves & -moves
    while move_to > 0:
        move_list.append((s, move_to, 5))
        moves = moves & (moves - 1)
        move_to = moves & -moves

    if s | c == c and s >> 3 & c != 0 and s >> 1 & o == 0 and s >> 2 & o == 0:
        if not is_attacked(gs, s) and not is_attacked(gs, s >> 1):
            move_list.append((s, s >> 2, 5))  # O-O

    if (s | c == c and s << 4 & c != 0 and s << 1 & o == 0 and
        s << 2 & o == 0 and s << 3 & o == 0):
        if not is_attacked(gs, s) and not is_attacked(gs, s << 1):
            move_list.append((s, s << 2, 5))  # O-O-O

    return move_list


def knight_moves(b, s):
    """Generate a bitboard of available knight moves.

    :param b: bitboard of all pieces of own color.
    :param s: bitboard of the king
    :return: list of pseudo-legal moves.
    """
    move_list = []

    move_from = s & -s
    while move_from > 0:
        moves = knight_attacks.get(move_from) & (b ^ 0xFFFFFFFFFFFFFFFF)

        move_to = moves & -moves
        while move_to > 0:
            move_list.append((move_from, move_to, 2))
            moves = moves & (moves - 1)
            move_to = moves & -moves
        s = s & (s - 1)
        move_from = s & -s

    return move_list


def rook_moves(o, b, s):
    """Generate a bitboard of available rook moves.

    :param o: bitboard of all pieces
    :param b: bitboard of all pieces of own color.
    :param s: bitboard of the king
    :return: list of pseudo-legal moves.
    """
    move_list = []

    move_from = s & -s
    while move_from > 0:
        (col, row), (_) = bitboard_to_pos.get(move_from)
        moves = ((ra.get(move_from).get(o & rmt[row]) |
                 fa.get(move_from).get(o & fmt[col])) &
                 (b ^ 0xFFFFFFFFFFFFFFFF))

        move_to = moves & -moves
        while move_to > 0:
            move_list.append((move_from, move_to, 1))
            moves = moves & (moves - 1)
            move_to = moves & -moves
        s = s & (s - 1)
        move_from = s & -s

    return move_list


def bishop_moves(o, b, s):
    """Generate a bitboard of available bishop moves.

    :param o: bitboard of all pieces
    :param b: bitboard of all pieces of own color.
    :param s: bitboard of the king
    :return: list of pseudo-legal moves.
    """
    move_list = []

    move_from = s & -s
    while move_from > 0:
        (_), (diag, adiag) = bitboard_to_pos.get(move_from)
        moves = ((da.get(move_from).get(o & dmt[diag]) |
                  ada.get(move_from).get(o & admt[adiag])) &
                 (b ^ 0xFFFFFFFFFFFFFFFF))
        move_to = moves & -moves

        while move_to > 0:
            move_list.append((move_from, move_to, 3))
            moves = moves & (moves - 1)
            move_to = moves & -moves
        s = s & (s - 1)
        move_from = s & -s

    return move_list


def queen_moves(o, b, s):
    """Generate a bitboard of available queen moves.

    :param o: bitboard of all pieces
    :param b: bitboard of all pieces of own color.
    :param s: bitboard of the king
    :return: list of pseudo-legal moves.
    """
    move_list = []
    move_from = s & -s

    while move_from > 0:
        (col, row), (diag, adiag) = bitboard_to_pos.get(move_from)
        moves = ((ra.get(move_from).get(o & rmt[row]) |
                 fa.get(move_from).get(o & fmt[col]) |
                 da.get(move_from).get(o & dmt[diag]) |
                 ada.get(move_from).get(o & admt[adiag])) &
                 (b ^ 0xFFFFFFFFFFFFFFFF))
        move_to = moves & -moves

        while move_to > 0:
            move_list.append((s, move_to, 4))
            moves = moves & (moves - 1)
            move_to = moves & -moves
        s = s & (s - 1)
        move_from = s & -s

    return move_list


def print_bitboard(bitboard):
    """Print a bitboard, mainly for development purposes."""
    board = "{:064b}".format(bitboard)
    for y in range(8):
        for x in range(8):
            print(board[8 * y + x] + " ", end="")
        print("")


if __name__ == "__main__":
    game = Game_State()
    game.print_board()
    while abs(game.score) < 15000:
        print('Input move:')
        move = input()
        game.move(move)
    print("Game over!")

    a = Game_State()
