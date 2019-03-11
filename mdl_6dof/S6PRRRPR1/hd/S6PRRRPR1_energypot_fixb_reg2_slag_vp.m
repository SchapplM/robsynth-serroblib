% Calculate inertial parameters regressor of potential energy for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:03:21
% EndTime: 2019-03-08 23:03:21
% DurationCPUTime: 0.25s
% Computational Cost: add. (333->107), mult. (459->158), div. (0->0), fcn. (534->14), ass. (0->54)
t113 = -pkin(9) - pkin(8);
t103 = sin(pkin(11));
t133 = g(1) * t103;
t105 = cos(pkin(11));
t132 = g(2) * t105;
t108 = sin(qJ(3));
t131 = t108 * pkin(3);
t111 = cos(qJ(3));
t93 = t111 * pkin(3) + pkin(2);
t104 = sin(pkin(6));
t128 = t103 * t104;
t130 = t105 * pkin(1) + pkin(7) * t128;
t106 = cos(pkin(6));
t129 = t106 * pkin(7) + qJ(1);
t127 = t104 * t105;
t126 = t104 * t108;
t109 = sin(qJ(2));
t125 = t104 * t109;
t124 = t104 * t111;
t112 = cos(qJ(2));
t123 = t104 * t112;
t122 = t106 * t108;
t121 = t106 * t109;
t120 = t106 * t112;
t102 = qJ(3) + qJ(4);
t97 = t103 * pkin(1);
t119 = -pkin(7) * t127 + t97;
t101 = -qJ(5) + t113;
t81 = t103 * t120 + t105 * t109;
t82 = -t103 * t121 + t105 * t112;
t96 = cos(t102);
t85 = pkin(4) * t96 + t93;
t95 = sin(t102);
t86 = pkin(4) * t95 + t131;
t118 = -t81 * t101 + t86 * t128 + t82 * t85 + t130;
t117 = t101 * t123 + t106 * t86 + t85 * t125 + t129;
t116 = -t132 + t133;
t80 = t103 * t112 + t105 * t121;
t94 = pkin(12) + t102;
t91 = sin(t94);
t92 = cos(t94);
t67 = t92 * t127 + t80 * t91;
t69 = -t92 * t128 + t82 * t91;
t73 = -t106 * t92 + t91 * t125;
t115 = g(1) * t69 + g(2) * t67 + g(3) * t73;
t79 = t103 * t109 - t105 * t120;
t114 = t80 * t85 - t79 * t101 + t97 + (-pkin(7) - t86) * t127;
t66 = -g(1) * t81 - g(2) * t79 + g(3) * t123;
t110 = cos(qJ(6));
t107 = sin(qJ(6));
t74 = t106 * t91 + t92 * t125;
t70 = t91 * t128 + t82 * t92;
t68 = -t91 * t127 + t80 * t92;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t103, t116, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t80 - g(3) * t125, -t66, -g(3) * t106 - t116 * t104, -g(1) * t130 - g(2) * t119 - g(3) * t129, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t126 + t82 * t111) - g(2) * (-t105 * t126 + t80 * t111) - g(3) * (t109 * t124 + t122) -g(1) * (t103 * t124 - t82 * t108) - g(2) * (-t105 * t124 - t80 * t108) - g(3) * (t106 * t111 - t108 * t125) t66, -g(1) * (t82 * pkin(2) + t81 * pkin(8) + t130) - g(2) * (t80 * pkin(2) + t79 * pkin(8) + t119) - g(3) * ((pkin(2) * t109 - pkin(8) * t112) * t104 + t129) 0, 0, 0, 0, 0, 0, -g(1) * (t95 * t128 + t82 * t96) - g(2) * (-t95 * t127 + t80 * t96) - g(3) * (t106 * t95 + t96 * t125) -g(1) * (t96 * t128 - t82 * t95) - g(2) * (-t96 * t127 - t80 * t95) - g(3) * (t106 * t96 - t95 * t125) t66, -g(1) * (-t81 * t113 + t82 * t93 + t130) - g(2) * (-t79 * t113 + t80 * t93 + t97) - g(3) * (pkin(3) * t122 + t129) + (-t131 * t133 - g(3) * (t109 * t93 + t112 * t113) - (-pkin(7) - t131) * t132) * t104, 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68 - g(3) * t74, t115, t66, -g(1) * t118 - g(2) * t114 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t107 + t70 * t110) - g(2) * (t79 * t107 + t68 * t110) - g(3) * (-t107 * t123 + t74 * t110) -g(1) * (-t70 * t107 + t81 * t110) - g(2) * (-t68 * t107 + t79 * t110) - g(3) * (-t74 * t107 - t110 * t123) -t115, -g(1) * (t70 * pkin(5) + t69 * pkin(10) + t118) - g(2) * (t68 * pkin(5) + t67 * pkin(10) + t114) - g(3) * (t74 * pkin(5) + t73 * pkin(10) + t117);];
U_reg  = t1;
