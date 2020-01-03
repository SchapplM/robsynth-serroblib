% Calculate inertial parameters regressor of potential energy for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:51
% EndTime: 2019-12-31 22:29:51
% DurationCPUTime: 0.19s
% Computational Cost: add. (132->75), mult. (169->97), div. (0->0), fcn. (168->10), ass. (0->38)
t100 = g(3) * pkin(5);
t83 = -pkin(8) - pkin(7);
t78 = sin(qJ(2));
t99 = g(3) * t78;
t77 = sin(qJ(3));
t98 = t77 * pkin(3);
t80 = cos(qJ(3));
t66 = t80 * pkin(3) + pkin(2);
t75 = -pkin(9) + t83;
t97 = t75 * t78;
t96 = t78 * t83;
t79 = sin(qJ(1));
t95 = t79 * t77;
t81 = cos(qJ(2));
t94 = t79 * t81;
t76 = qJ(3) + qJ(4);
t69 = qJ(5) + t76;
t64 = sin(t69);
t82 = cos(qJ(1));
t93 = t82 * t64;
t65 = cos(t69);
t92 = t82 * t65;
t67 = sin(t76);
t91 = t82 * t67;
t68 = cos(t76);
t90 = t82 * t68;
t89 = t82 * t77;
t88 = t82 * t80;
t87 = t82 * pkin(1) + t79 * pkin(6);
t71 = t79 * pkin(1);
t86 = -t82 * pkin(6) + t71;
t85 = pkin(2) * t81 + pkin(7) * t78;
t84 = g(1) * t82 + g(2) * t79;
t63 = g(1) * t79 - g(2) * t82;
t62 = pkin(4) * t67 + t98;
t61 = pkin(4) * t68 + t66;
t60 = -g(3) * t81 + t78 * t84;
t1 = [0, 0, 0, 0, 0, 0, -t84, t63, -g(3), -t100, 0, 0, 0, 0, 0, 0, -t81 * t84 - t99, t60, -t63, -g(1) * t87 - g(2) * t86 - t100, 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t88 + t95) - g(2) * (t80 * t94 - t89) - t80 * t99, -g(1) * (t79 * t80 - t81 * t89) - g(2) * (-t77 * t94 - t88) + t77 * t99, -t60, -g(1) * (t82 * t85 + t87) - g(2) * (t79 * t85 + t86) - g(3) * (t78 * pkin(2) - t81 * pkin(7) + pkin(5)), 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t67 + t81 * t90) - g(2) * (t68 * t94 - t91) - t68 * t99, -g(1) * (t79 * t68 - t81 * t91) - g(2) * (-t67 * t94 - t90) + t67 * t99, -t60, -g(1) * (pkin(3) * t95 + t87) - g(2) * (t66 * t94 - t79 * t96 + t71) - g(3) * (t78 * t66 + t81 * t83 + pkin(5)) + (-g(1) * (t66 * t81 - t96) - g(2) * (-pkin(6) - t98)) * t82, 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t64 + t81 * t92) - g(2) * (t65 * t94 - t93) - t65 * t99, -g(1) * (t79 * t65 - t81 * t93) - g(2) * (-t64 * t94 - t92) + t64 * t99, -t60, -g(1) * (t79 * t62 + t87) - g(2) * (t61 * t94 - t79 * t97 + t71) - g(3) * (t78 * t61 + t81 * t75 + pkin(5)) + (-g(1) * (t61 * t81 - t97) - g(2) * (-pkin(6) - t62)) * t82;];
U_reg = t1;
