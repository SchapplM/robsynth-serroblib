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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:44:53
% EndTime: 2019-12-05 17:44:54
% DurationCPUTime: 0.18s
% Computational Cost: add. (132->73), mult. (169->98), div. (0->0), fcn. (168->10), ass. (0->42)
t101 = g(1) * pkin(5);
t75 = sin(pkin(8));
t100 = g(1) * t75;
t79 = sin(qJ(1));
t99 = g(2) * t79;
t74 = sin(pkin(9));
t98 = t74 * pkin(3);
t76 = cos(pkin(9));
t64 = t76 * pkin(3) + pkin(2);
t80 = cos(qJ(1));
t97 = t75 * t80;
t77 = cos(pkin(8));
t96 = t77 * t80;
t73 = pkin(9) + qJ(4);
t67 = qJ(5) + t73;
t62 = sin(t67);
t95 = t79 * t62;
t63 = cos(t67);
t94 = t79 * t63;
t65 = sin(t73);
t93 = t79 * t65;
t66 = cos(t73);
t92 = t79 * t66;
t91 = t79 * t74;
t90 = t79 * t76;
t89 = t80 * t62;
t88 = t80 * t63;
t87 = t80 * t65;
t86 = t80 * t66;
t85 = t80 * t74;
t84 = t80 * t76;
t78 = -pkin(6) - qJ(3);
t83 = t80 * pkin(1) + t79 * qJ(2);
t82 = -g(3) * t80 + t99;
t81 = pkin(2) * t77 + qJ(3) * t75;
t72 = -pkin(7) + t78;
t70 = t80 * qJ(2);
t61 = g(2) * t80 + g(3) * t79;
t60 = pkin(4) * t65 + t98;
t59 = pkin(4) * t66 + t64;
t58 = g(1) * t77 + t82 * t75;
t1 = [0, 0, 0, 0, 0, 0, t82, t61, -g(1), -t101, 0, 0, 0, 0, 0, 0, t82 * t77 - t100, -t58, -t61, -t101 - g(2) * (-t79 * pkin(1) + t70) - g(3) * t83, 0, 0, 0, 0, 0, 0, -t76 * t100 - g(2) * (-t77 * t90 + t85) - g(3) * (t77 * t84 + t91), t74 * t100 - g(2) * (t77 * t91 + t84) - g(3) * (-t77 * t85 + t90), t58, -g(1) * (t75 * pkin(2) - t77 * qJ(3) + pkin(5)) - g(2) * t70 - g(3) * (t81 * t80 + t83) - (-pkin(1) - t81) * t99, 0, 0, 0, 0, 0, 0, -t66 * t100 - g(2) * (-t77 * t92 + t87) - g(3) * (t77 * t86 + t93), t65 * t100 - g(2) * (t77 * t93 + t86) - g(3) * (-t77 * t87 + t92), t58, -g(1) * (t75 * t64 + t77 * t78 + pkin(5)) - g(2) * (pkin(3) * t85 + t70) - g(3) * (t64 * t96 - t78 * t97 + t83) + (-g(2) * (-t64 * t77 + t75 * t78 - pkin(1)) - g(3) * t98) * t79, 0, 0, 0, 0, 0, 0, -t63 * t100 - g(2) * (-t77 * t94 + t89) - g(3) * (t77 * t88 + t95), t62 * t100 - g(2) * (t77 * t95 + t88) - g(3) * (-t77 * t89 + t94), t58, -g(1) * (t75 * t59 + t77 * t72 + pkin(5)) - g(2) * (t80 * t60 + t70) - g(3) * (t59 * t96 - t72 * t97 + t83) + (-g(2) * (-t59 * t77 + t72 * t75 - pkin(1)) - g(3) * t60) * t79;];
U_reg = t1;
