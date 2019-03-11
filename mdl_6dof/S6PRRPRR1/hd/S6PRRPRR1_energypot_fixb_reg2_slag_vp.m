% Calculate inertial parameters regressor of potential energy for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:54:05
% EndTime: 2019-03-08 21:54:06
% DurationCPUTime: 0.23s
% Computational Cost: add. (333->107), mult. (459->158), div. (0->0), fcn. (534->14), ass. (0->54)
t104 = sin(pkin(11));
t134 = g(1) * t104;
t106 = cos(pkin(11));
t133 = g(2) * t106;
t110 = sin(qJ(3));
t132 = t110 * pkin(3);
t108 = -qJ(4) - pkin(8);
t113 = cos(qJ(3));
t94 = t113 * pkin(3) + pkin(2);
t105 = sin(pkin(6));
t130 = t104 * t105;
t131 = t106 * pkin(1) + pkin(7) * t130;
t129 = t105 * t106;
t128 = t105 * t110;
t111 = sin(qJ(2));
t127 = t105 * t111;
t126 = t105 * t113;
t114 = cos(qJ(2));
t125 = t105 * t114;
t107 = cos(pkin(6));
t124 = t107 * t110;
t123 = t107 * t111;
t122 = t107 * t114;
t121 = t107 * pkin(7) + qJ(1);
t103 = qJ(3) + pkin(12);
t98 = t104 * pkin(1);
t120 = -pkin(7) * t129 + t98;
t102 = -pkin(9) + t108;
t82 = t104 * t122 + t106 * t111;
t83 = -t104 * t123 + t106 * t114;
t96 = cos(t103);
t86 = pkin(4) * t96 + t94;
t95 = sin(t103);
t87 = pkin(4) * t95 + t132;
t119 = -t82 * t102 + t87 * t130 + t83 * t86 + t131;
t118 = t102 * t125 + t107 * t87 + t86 * t127 + t121;
t117 = -t133 + t134;
t81 = t104 * t114 + t106 * t123;
t97 = qJ(5) + t103;
t92 = sin(t97);
t93 = cos(t97);
t68 = t93 * t129 + t81 * t92;
t70 = -t93 * t130 + t83 * t92;
t74 = -t107 * t93 + t92 * t127;
t116 = g(1) * t70 + g(2) * t68 + g(3) * t74;
t80 = t104 * t111 - t106 * t122;
t115 = t81 * t86 - t80 * t102 + t98 + (-pkin(7) - t87) * t129;
t67 = -g(1) * t82 - g(2) * t80 + g(3) * t125;
t112 = cos(qJ(6));
t109 = sin(qJ(6));
t75 = t107 * t92 + t93 * t127;
t71 = t92 * t130 + t83 * t93;
t69 = -t92 * t129 + t81 * t93;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t106 - g(2) * t104, t117, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t127, -t67, -g(3) * t107 - t117 * t105, -g(1) * t131 - g(2) * t120 - g(3) * t121, 0, 0, 0, 0, 0, 0, -g(1) * (t104 * t128 + t83 * t113) - g(2) * (-t106 * t128 + t81 * t113) - g(3) * (t111 * t126 + t124) -g(1) * (t104 * t126 - t83 * t110) - g(2) * (-t106 * t126 - t81 * t110) - g(3) * (t107 * t113 - t110 * t127) t67, -g(1) * (t83 * pkin(2) + t82 * pkin(8) + t131) - g(2) * (t81 * pkin(2) + t80 * pkin(8) + t120) - g(3) * ((pkin(2) * t111 - pkin(8) * t114) * t105 + t121) 0, 0, 0, 0, 0, 0, -g(1) * (t95 * t130 + t83 * t96) - g(2) * (-t95 * t129 + t81 * t96) - g(3) * (t107 * t95 + t96 * t127) -g(1) * (t96 * t130 - t83 * t95) - g(2) * (-t96 * t129 - t81 * t95) - g(3) * (t107 * t96 - t95 * t127) t67, -g(1) * (-t82 * t108 + t83 * t94 + t131) - g(2) * (-t80 * t108 + t81 * t94 + t98) - g(3) * (pkin(3) * t124 + t121) + (-t132 * t134 - g(3) * (t108 * t114 + t111 * t94) - (-pkin(7) - t132) * t133) * t105, 0, 0, 0, 0, 0, 0, -g(1) * t71 - g(2) * t69 - g(3) * t75, t116, t67, -g(1) * t119 - g(2) * t115 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t109 + t71 * t112) - g(2) * (t80 * t109 + t69 * t112) - g(3) * (-t109 * t125 + t75 * t112) -g(1) * (-t71 * t109 + t82 * t112) - g(2) * (-t69 * t109 + t80 * t112) - g(3) * (-t75 * t109 - t112 * t125) -t116, -g(1) * (t71 * pkin(5) + t70 * pkin(10) + t119) - g(2) * (t69 * pkin(5) + t68 * pkin(10) + t115) - g(3) * (t75 * pkin(5) + t74 * pkin(10) + t118);];
U_reg  = t1;
