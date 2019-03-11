% Calculate inertial parameters regressor of potential energy for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:27
% EndTime: 2019-03-08 18:47:27
% DurationCPUTime: 0.38s
% Computational Cost: add. (539->111), mult. (1376->172), div. (0->0), fcn. (1753->16), ass. (0->65)
t112 = sin(pkin(7));
t117 = cos(pkin(7));
t118 = cos(pkin(6));
t113 = sin(pkin(6));
t115 = cos(pkin(12));
t149 = t113 * t115;
t133 = t112 * t149 - t117 * t118;
t111 = sin(pkin(11));
t147 = t113 * t117;
t110 = sin(pkin(12));
t116 = cos(pkin(11));
t150 = t111 * t118;
t93 = -t110 * t116 - t115 * t150;
t135 = t111 * t147 - t112 * t93;
t156 = cos(qJ(3));
t155 = cos(qJ(4));
t151 = t111 * t113;
t153 = pkin(1) * t116 + qJ(2) * t151;
t152 = t110 * t113;
t148 = t113 * t116;
t146 = t116 * t118;
t144 = qJ(2) * t118 + qJ(1);
t109 = sin(pkin(13));
t141 = pkin(5) * t109 + pkin(9);
t140 = t112 * t156;
t139 = t117 * t156;
t138 = t113 * t140;
t137 = g(1) * t111 - g(2) * t116;
t136 = pkin(1) * t111 - qJ(2) * t148;
t91 = -t110 * t111 + t115 * t146;
t134 = t112 * t91 + t116 * t147;
t120 = sin(qJ(4));
t121 = sin(qJ(3));
t92 = t110 * t146 + t111 * t115;
t76 = t92 * t156 + (-t112 * t148 + t117 * t91) * t121;
t69 = t120 * t76 + t134 * t155;
t94 = -t110 * t150 + t115 * t116;
t78 = t94 * t156 + (t112 * t151 + t117 * t93) * t121;
t71 = t120 * t78 - t135 * t155;
t85 = t118 * t112 * t121 + (t115 * t117 * t121 + t110 * t156) * t113;
t79 = t120 * t85 + t133 * t155;
t132 = g(1) * t71 + g(2) * t69 + g(3) * t79;
t75 = t116 * t138 + t121 * t92 - t139 * t91;
t77 = -t111 * t138 + t121 * t94 - t139 * t93;
t84 = -t118 * t140 + t121 * t152 - t139 * t149;
t131 = g(1) * t77 + g(2) * t75 + g(3) * t84;
t130 = t94 * pkin(2) + pkin(8) * t135 + t153;
t129 = pkin(3) * t78 + t130;
t128 = pkin(2) * t152 - pkin(8) * t133 + t144;
t127 = pkin(3) * t85 + t128;
t126 = pkin(9) * t77 + t129;
t125 = pkin(9) * t84 + t127;
t124 = pkin(2) * t92 - pkin(8) * t134 + t136;
t123 = pkin(3) * t76 + t124;
t122 = t75 * pkin(9) + t123;
t119 = -pkin(10) - qJ(5);
t114 = cos(pkin(13));
t108 = pkin(13) + qJ(6);
t104 = cos(t108);
t103 = sin(t108);
t102 = pkin(5) * t114 + pkin(4);
t80 = -t120 * t133 + t155 * t85;
t72 = t120 * t135 + t155 * t78;
t70 = -t120 * t134 + t155 * t76;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t116 - g(2) * t111, t137, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t94 - g(2) * t92 - g(3) * t152, -g(1) * t93 - g(2) * t91 - g(3) * t149, -g(3) * t118 - t113 * t137, -g(1) * t153 - g(2) * t136 - g(3) * t144, 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t76 - g(3) * t85, t131, -g(1) * t135 + g(2) * t134 + g(3) * t133, -g(1) * t130 - g(2) * t124 - g(3) * t128, 0, 0, 0, 0, 0, 0, -g(1) * t72 - g(2) * t70 - g(3) * t80, t132, -t131, -g(1) * t126 - g(2) * t122 - g(3) * t125, 0, 0, 0, 0, 0, 0, -g(1) * (t109 * t77 + t114 * t72) - g(2) * (t109 * t75 + t114 * t70) - g(3) * (t109 * t84 + t114 * t80) -g(1) * (-t109 * t72 + t114 * t77) - g(2) * (-t109 * t70 + t114 * t75) - g(3) * (-t109 * t80 + t114 * t84) -t132, -g(1) * (pkin(4) * t72 + qJ(5) * t71 + t126) - g(2) * (t70 * pkin(4) + t69 * qJ(5) + t122) - g(3) * (pkin(4) * t80 + qJ(5) * t79 + t125) 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t77 + t104 * t72) - g(2) * (t103 * t75 + t104 * t70) - g(3) * (t103 * t84 + t104 * t80) -g(1) * (-t103 * t72 + t104 * t77) - g(2) * (-t103 * t70 + t104 * t75) - g(3) * (-t103 * t80 + t104 * t84) -t132, -g(1) * (t72 * t102 - t71 * t119 + t141 * t77 + t129) - g(2) * (t70 * t102 - t69 * t119 + t141 * t75 + t123) - g(3) * (t80 * t102 - t79 * t119 + t141 * t84 + t127);];
U_reg  = t1;
