% Calculate inertial parameters regressor of potential energy for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:29:54
% EndTime: 2019-03-08 20:29:54
% DurationCPUTime: 0.30s
% Computational Cost: add. (376->102), mult. (841->157), div. (0->0), fcn. (1059->14), ass. (0->58)
t100 = sin(pkin(6));
t132 = pkin(7) * t100;
t101 = cos(pkin(11));
t102 = cos(pkin(6));
t105 = sin(qJ(2));
t125 = t102 * t105;
t82 = pkin(2) * t125 + (-pkin(7) - qJ(3)) * t100;
t108 = cos(qJ(2));
t92 = pkin(2) * t108 + pkin(1);
t99 = sin(pkin(11));
t131 = t101 * t82 + t99 * t92;
t130 = t102 * pkin(7) + qJ(1);
t129 = cos(pkin(12));
t104 = sin(qJ(4));
t128 = t100 * t104;
t127 = t100 * t105;
t107 = cos(qJ(4));
t126 = t100 * t107;
t124 = t102 * t108;
t98 = sin(pkin(12));
t83 = -t105 * t129 - t108 * t98;
t81 = t83 * t102;
t120 = t108 * t129;
t84 = -t105 * t98 + t120;
t70 = -t101 * t81 + t84 * t99;
t123 = t70 * pkin(3) + t131;
t103 = sin(qJ(5));
t122 = pkin(5) * t103 + pkin(8);
t121 = t101 * t92 - t99 * t82;
t119 = pkin(2) * t127 + t102 * qJ(3) + t130;
t72 = t101 * t84 + t81 * t99;
t118 = t72 * pkin(3) + t121;
t80 = t83 * t100;
t117 = -t80 * pkin(3) + t119;
t116 = g(1) * t99 - g(2) * t101;
t110 = t102 * t84;
t69 = t101 * t110 + t99 * t83;
t115 = -t69 * pkin(8) + t123;
t71 = t101 * t83 - t99 * t110;
t114 = -t71 * pkin(8) + t118;
t63 = t101 * t126 + t104 * t70;
t65 = t104 * t72 - t99 * t126;
t73 = -t102 * t107 - t104 * t80;
t113 = g(1) * t65 + g(2) * t63 + g(3) * t73;
t79 = -t100 * t120 + t98 * t127;
t112 = g(1) * t71 + g(2) * t69 - g(3) * t79;
t111 = t79 * pkin(8) + t117;
t109 = -pkin(10) - pkin(9);
t106 = cos(qJ(5));
t97 = qJ(5) + qJ(6);
t95 = cos(t97);
t94 = sin(t97);
t91 = pkin(5) * t106 + pkin(4);
t78 = -g(3) * t102 - t116 * t100;
t74 = t102 * t104 - t107 * t80;
t66 = t107 * t72 + t99 * t128;
t64 = -t101 * t128 + t107 * t70;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t101 - g(2) * t99, t116, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t108 - t99 * t125) - g(2) * (t101 * t125 + t108 * t99) - g(3) * t127, -g(1) * (-t101 * t105 - t99 * t124) - g(2) * (t101 * t124 - t105 * t99) - g(3) * t100 * t108, t78, -g(1) * (pkin(1) * t101 + t99 * t132) - g(2) * (pkin(1) * t99 - t101 * t132) - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * t72 - g(2) * t70 + g(3) * t80, -t112, t78, -g(1) * t121 - g(2) * t131 - g(3) * t119, 0, 0, 0, 0, 0, 0, -g(1) * t66 - g(2) * t64 - g(3) * t74, t113, t112, -g(1) * t114 - g(2) * t115 - g(3) * t111, 0, 0, 0, 0, 0, 0, -g(1) * (-t103 * t71 + t106 * t66) - g(2) * (-t103 * t69 + t106 * t64) - g(3) * (t103 * t79 + t106 * t74) -g(1) * (-t103 * t66 - t106 * t71) - g(2) * (-t103 * t64 - t106 * t69) - g(3) * (-t103 * t74 + t106 * t79) -t113, -g(1) * (pkin(4) * t66 + pkin(9) * t65 + t114) - g(2) * (pkin(4) * t64 + pkin(9) * t63 + t115) - g(3) * (pkin(4) * t74 + pkin(9) * t73 + t111) 0, 0, 0, 0, 0, 0, -g(1) * (t66 * t95 - t71 * t94) - g(2) * (t64 * t95 - t69 * t94) - g(3) * (t74 * t95 + t79 * t94) -g(1) * (-t66 * t94 - t71 * t95) - g(2) * (-t64 * t94 - t69 * t95) - g(3) * (-t74 * t94 + t79 * t95) -t113, -g(1) * (-t109 * t65 - t122 * t71 + t66 * t91 + t118) - g(2) * (-t109 * t63 - t122 * t69 + t64 * t91 + t123) - g(3) * (-t73 * t109 + t122 * t79 + t74 * t91 + t117);];
U_reg  = t1;
