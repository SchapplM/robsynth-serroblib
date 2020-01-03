% Calculate inertial parameters regressor of potential energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:31:21
% EndTime: 2020-01-03 11:31:21
% DurationCPUTime: 0.21s
% Computational Cost: add. (132->71), mult. (169->95), div. (0->0), fcn. (168->10), ass. (0->41)
t72 = sin(pkin(9));
t83 = t72 * pkin(3) + qJ(2);
t100 = g(1) * pkin(5);
t73 = sin(pkin(8));
t99 = g(1) * t73;
t77 = sin(qJ(1));
t69 = t77 * pkin(1);
t98 = g(2) * t69;
t74 = cos(pkin(9));
t64 = t74 * pkin(3) + pkin(2);
t71 = pkin(9) + qJ(4);
t67 = qJ(5) + t71;
t62 = sin(t67);
t96 = t77 * t62;
t63 = cos(t67);
t95 = t77 * t63;
t65 = sin(t71);
t94 = t77 * t65;
t66 = cos(t71);
t93 = t77 * t66;
t92 = t77 * t72;
t91 = t77 * t74;
t78 = cos(qJ(1));
t90 = t78 * t62;
t89 = t78 * t63;
t88 = t78 * t65;
t87 = t78 * t66;
t86 = t78 * t72;
t85 = t78 * t74;
t76 = -pkin(6) - qJ(3);
t84 = pkin(4) * t65 + t83;
t82 = -g(2) * t77 + g(3) * t78;
t75 = cos(pkin(8));
t81 = pkin(2) * t75 + qJ(3) * t73;
t59 = pkin(4) * t66 + t64;
t70 = -pkin(7) + t76;
t80 = t59 * t75 - t70 * t73;
t79 = t64 * t75 - t73 * t76;
t61 = g(2) * t78 + g(3) * t77;
t58 = g(1) * t75 + t82 * t73;
t1 = [0, 0, 0, 0, 0, 0, t82, -t61, -g(1), -t100, 0, 0, 0, 0, 0, 0, t82 * t75 - t99, -t58, t61, -t100 - g(2) * (-t78 * qJ(2) + t69) - g(3) * (-t78 * pkin(1) - t77 * qJ(2)), 0, 0, 0, 0, 0, 0, -t74 * t99 - g(2) * (t75 * t91 - t86) - g(3) * (-t75 * t85 - t92), t72 * t99 - g(2) * (-t75 * t92 - t85) - g(3) * (t75 * t86 - t91), t58, -g(1) * (t73 * pkin(2) - t75 * qJ(3) + pkin(5)) - t98 + (-g(2) * t81 + g(3) * qJ(2)) * t77 + (g(2) * qJ(2) - g(3) * (-pkin(1) - t81)) * t78, 0, 0, 0, 0, 0, 0, -t66 * t99 - g(2) * (t75 * t93 - t88) - g(3) * (-t75 * t87 - t94), t65 * t99 - g(2) * (-t75 * t94 - t87) - g(3) * (t75 * t88 - t93), t58, -g(1) * (t73 * t64 + t75 * t76 + pkin(5)) - t98 + (-g(2) * t79 + g(3) * t83) * t77 + (g(2) * t83 - g(3) * (-pkin(1) - t79)) * t78, 0, 0, 0, 0, 0, 0, -t63 * t99 - g(2) * (t75 * t95 - t90) - g(3) * (-t75 * t89 - t96), t62 * t99 - g(2) * (-t75 * t96 - t89) - g(3) * (t75 * t90 - t95), t58, -g(1) * (t73 * t59 + t75 * t70 + pkin(5)) - t98 + (-g(2) * t80 + g(3) * t84) * t77 + (g(2) * t84 - g(3) * (-pkin(1) - t80)) * t78;];
U_reg = t1;
