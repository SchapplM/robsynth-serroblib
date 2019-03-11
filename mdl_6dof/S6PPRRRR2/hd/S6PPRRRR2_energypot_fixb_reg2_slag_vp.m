% Calculate inertial parameters regressor of potential energy for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:05:23
% EndTime: 2019-03-08 19:05:23
% DurationCPUTime: 0.37s
% Computational Cost: add. (539->111), mult. (1376->172), div. (0->0), fcn. (1753->16), ass. (0->65)
t112 = sin(pkin(7));
t116 = cos(pkin(7));
t117 = cos(pkin(6));
t113 = sin(pkin(6));
t114 = cos(pkin(13));
t151 = t113 * t114;
t134 = t112 * t151 - t116 * t117;
t111 = sin(pkin(12));
t149 = t113 * t116;
t110 = sin(pkin(13));
t115 = cos(pkin(12));
t152 = t111 * t117;
t94 = -t110 * t115 - t114 * t152;
t136 = t111 * t149 - t112 * t94;
t157 = cos(qJ(3));
t156 = cos(qJ(4));
t154 = t110 * t113;
t153 = t111 * t113;
t150 = t113 * t115;
t148 = t115 * t117;
t146 = qJ(2) * t117 + qJ(1);
t145 = pkin(1) * t115 + qJ(2) * t153;
t118 = sin(qJ(5));
t142 = pkin(5) * t118 + pkin(9);
t141 = t112 * t157;
t140 = t116 * t157;
t139 = t113 * t141;
t138 = g(1) * t111 - g(2) * t115;
t137 = pkin(1) * t111 - qJ(2) * t150;
t92 = -t110 * t111 + t114 * t148;
t135 = t112 * t92 + t115 * t149;
t119 = sin(qJ(4));
t120 = sin(qJ(3));
t93 = t110 * t148 + t111 * t114;
t77 = t93 * t157 + (-t112 * t150 + t116 * t92) * t120;
t70 = t119 * t77 + t135 * t156;
t95 = -t110 * t152 + t114 * t115;
t79 = t95 * t157 + (t112 * t153 + t116 * t94) * t120;
t72 = t119 * t79 - t136 * t156;
t86 = t117 * t112 * t120 + (t114 * t116 * t120 + t110 * t157) * t113;
t80 = t119 * t86 + t134 * t156;
t133 = g(1) * t72 + g(2) * t70 + g(3) * t80;
t76 = t115 * t139 + t120 * t93 - t140 * t92;
t78 = -t111 * t139 + t120 * t95 - t140 * t94;
t85 = -t117 * t141 + t120 * t154 - t140 * t151;
t132 = g(1) * t78 + g(2) * t76 + g(3) * t85;
t131 = t95 * pkin(2) + pkin(8) * t136 + t145;
t130 = pkin(3) * t79 + t131;
t129 = pkin(2) * t154 - pkin(8) * t134 + t146;
t128 = pkin(3) * t86 + t129;
t127 = pkin(9) * t78 + t130;
t126 = pkin(9) * t85 + t128;
t125 = pkin(2) * t93 - pkin(8) * t135 + t137;
t124 = pkin(3) * t77 + t125;
t123 = t76 * pkin(9) + t124;
t122 = -pkin(11) - pkin(10);
t121 = cos(qJ(5));
t109 = qJ(5) + qJ(6);
t106 = cos(t109);
t105 = sin(t109);
t103 = pkin(5) * t121 + pkin(4);
t81 = -t119 * t134 + t156 * t86;
t73 = t119 * t136 + t156 * t79;
t71 = -t119 * t135 + t156 * t77;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t115 - g(2) * t111, t138, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t154, -g(1) * t94 - g(2) * t92 - g(3) * t151, -g(3) * t117 - t113 * t138, -g(1) * t145 - g(2) * t137 - g(3) * t146, 0, 0, 0, 0, 0, 0, -g(1) * t79 - g(2) * t77 - g(3) * t86, t132, -g(1) * t136 + g(2) * t135 + g(3) * t134, -g(1) * t131 - g(2) * t125 - g(3) * t129, 0, 0, 0, 0, 0, 0, -g(1) * t73 - g(2) * t71 - g(3) * t81, t133, -t132, -g(1) * t127 - g(2) * t123 - g(3) * t126, 0, 0, 0, 0, 0, 0, -g(1) * (t118 * t78 + t121 * t73) - g(2) * (t118 * t76 + t121 * t71) - g(3) * (t118 * t85 + t121 * t81) -g(1) * (-t118 * t73 + t121 * t78) - g(2) * (-t118 * t71 + t121 * t76) - g(3) * (-t118 * t81 + t121 * t85) -t133, -g(1) * (pkin(4) * t73 + pkin(10) * t72 + t127) - g(2) * (t71 * pkin(4) + t70 * pkin(10) + t123) - g(3) * (pkin(4) * t81 + pkin(10) * t80 + t126) 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t78 + t106 * t73) - g(2) * (t105 * t76 + t106 * t71) - g(3) * (t105 * t85 + t106 * t81) -g(1) * (-t105 * t73 + t106 * t78) - g(2) * (-t105 * t71 + t106 * t76) - g(3) * (-t105 * t81 + t106 * t85) -t133, -g(1) * (t73 * t103 - t72 * t122 + t142 * t78 + t130) - g(2) * (t71 * t103 - t70 * t122 + t142 * t76 + t124) - g(3) * (t81 * t103 - t80 * t122 + t142 * t85 + t128);];
U_reg  = t1;
