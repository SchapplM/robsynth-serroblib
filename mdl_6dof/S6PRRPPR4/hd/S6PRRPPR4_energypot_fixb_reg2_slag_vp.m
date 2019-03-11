% Calculate inertial parameters regressor of potential energy for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:16:25
% EndTime: 2019-03-08 21:16:25
% DurationCPUTime: 0.23s
% Computational Cost: add. (337->88), mult. (795->126), div. (0->0), fcn. (997->12), ass. (0->57)
t138 = cos(qJ(3));
t103 = sin(pkin(6));
t137 = pkin(7) * t103;
t136 = -pkin(9) + qJ(4);
t102 = sin(pkin(10));
t105 = cos(pkin(10));
t135 = pkin(1) * t105 + t102 * t137;
t107 = sin(qJ(3));
t126 = t103 * t138;
t110 = cos(qJ(2));
t108 = sin(qJ(2));
t131 = cos(pkin(6));
t125 = t108 * t131;
t88 = t102 * t110 + t105 * t125;
t78 = t105 * t126 + t107 * t88;
t134 = t78 * qJ(4);
t90 = -t102 * t125 + t105 * t110;
t80 = -t102 * t126 + t107 * t90;
t133 = t80 * qJ(4);
t129 = t103 * t108;
t91 = t107 * t129 - t131 * t138;
t132 = t91 * qJ(4);
t130 = t103 * t107;
t128 = t103 * t110;
t127 = pkin(7) * t131 + qJ(1);
t124 = t110 * t131;
t123 = pkin(1) * t102 - t105 * t137;
t122 = g(1) * t102 - g(2) * t105;
t89 = t102 * t124 + t105 * t108;
t121 = pkin(2) * t90 + t89 * pkin(8) + t135;
t81 = t102 * t130 + t138 * t90;
t120 = pkin(3) * t81 + t121;
t119 = pkin(2) * t129 - pkin(8) * t128 + t127;
t101 = sin(pkin(11));
t104 = cos(pkin(11));
t79 = -t105 * t130 + t138 * t88;
t87 = t102 * t108 - t105 * t124;
t69 = t101 * t79 - t104 * t87;
t71 = t101 * t81 - t104 * t89;
t92 = t107 * t131 + t108 * t126;
t76 = t101 * t92 + t104 * t128;
t118 = g(1) * t71 + g(2) * t69 + g(3) * t76;
t66 = g(1) * t80 + g(2) * t78 + g(3) * t91;
t117 = pkin(3) * t92 + t119;
t116 = pkin(2) * t88 + t87 * pkin(8) + t123;
t115 = -g(1) * t89 - g(2) * t87 + g(3) * t128;
t114 = pkin(3) * t79 + t116;
t72 = t101 * t89 + t104 * t81;
t113 = pkin(4) * t72 + t71 * qJ(5) + t120;
t77 = -t101 * t128 + t104 * t92;
t112 = pkin(4) * t77 + t76 * qJ(5) + t117;
t70 = t101 * t87 + t104 * t79;
t111 = pkin(4) * t70 + t69 * qJ(5) + t114;
t109 = cos(qJ(6));
t106 = sin(qJ(6));
t64 = -g(1) * t72 - g(2) * t70 - g(3) * t77;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t102, t122, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t90 - g(2) * t88 - g(3) * t129, -t115, -g(3) * t131 - t103 * t122, -g(1) * t135 - g(2) * t123 - g(3) * t127, 0, 0, 0, 0, 0, 0, -g(1) * t81 - g(2) * t79 - g(3) * t92, t66, t115, -g(1) * t121 - g(2) * t116 - g(3) * t119, 0, 0, 0, 0, 0, 0, t64, t118, -t66, -g(1) * (t120 + t133) - g(2) * (t114 + t134) - g(3) * (t117 + t132) 0, 0, 0, 0, 0, 0, t64, -t66, -t118, -g(1) * (t113 + t133) - g(2) * (t111 + t134) - g(3) * (t112 + t132) 0, 0, 0, 0, 0, 0, -g(1) * (t106 * t71 + t109 * t72) - g(2) * (t106 * t69 + t109 * t70) - g(3) * (t106 * t76 + t109 * t77) -g(1) * (-t106 * t72 + t109 * t71) - g(2) * (-t106 * t70 + t109 * t69) - g(3) * (-t106 * t77 + t109 * t76) t66, -g(1) * (t72 * pkin(5) + t136 * t80 + t113) - g(2) * (t70 * pkin(5) + t136 * t78 + t111) - g(3) * (t77 * pkin(5) + t136 * t91 + t112);];
U_reg  = t1;
