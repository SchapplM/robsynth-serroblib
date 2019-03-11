% Calculate inertial parameters regressor of potential energy for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:08
% EndTime: 2019-03-08 18:58:08
% DurationCPUTime: 0.33s
% Computational Cost: add. (576->95), mult. (1520->149), div. (0->0), fcn. (1949->14), ass. (0->67)
t116 = sin(pkin(7));
t120 = cos(pkin(7));
t121 = cos(pkin(6));
t117 = sin(pkin(6));
t118 = cos(pkin(12));
t154 = t117 * t118;
t138 = t116 * t154 - t121 * t120;
t114 = sin(pkin(12));
t119 = cos(pkin(11));
t115 = sin(pkin(11));
t155 = t115 * t121;
t102 = -t119 * t114 - t118 * t155;
t152 = t117 * t120;
t139 = -t102 * t116 + t115 * t152;
t160 = cos(qJ(3));
t159 = cos(qJ(4));
t157 = t114 * t117;
t156 = t115 * t117;
t153 = t117 * t119;
t151 = t119 * t121;
t149 = t121 * qJ(2) + qJ(1);
t148 = t119 * pkin(1) + qJ(2) * t156;
t145 = t116 * t160;
t144 = t120 * t160;
t143 = t117 * t145;
t142 = g(1) * t115 - g(2) * t119;
t141 = t115 * pkin(1) - qJ(2) * t153;
t100 = -t115 * t114 + t118 * t151;
t140 = t100 * t116 + t119 * t152;
t122 = sin(qJ(5));
t125 = cos(qJ(5));
t123 = sin(qJ(4));
t101 = t114 * t151 + t115 * t118;
t124 = sin(qJ(3));
t84 = t101 * t160 + (t100 * t120 - t116 * t153) * t124;
t73 = -t140 * t123 + t84 * t159;
t83 = -t100 * t144 + t101 * t124 + t119 * t143;
t66 = t73 * t122 - t83 * t125;
t103 = -t114 * t155 + t119 * t118;
t86 = t103 * t160 + (t102 * t120 + t116 * t156) * t124;
t75 = t139 * t123 + t86 * t159;
t85 = -t102 * t144 + t103 * t124 - t115 * t143;
t68 = t75 * t122 - t85 * t125;
t94 = t121 * t116 * t124 + (t118 * t120 * t124 + t160 * t114) * t117;
t88 = -t138 * t123 + t94 * t159;
t93 = -t121 * t145 + t124 * t157 - t144 * t154;
t76 = t88 * t122 - t93 * t125;
t137 = g(1) * t68 + g(2) * t66 + g(3) * t76;
t72 = t84 * t123 + t140 * t159;
t74 = t86 * t123 - t139 * t159;
t87 = t94 * t123 + t138 * t159;
t136 = g(1) * t74 + g(2) * t72 + g(3) * t87;
t135 = g(1) * t85 + g(2) * t83 + g(3) * t93;
t134 = t103 * pkin(2) + t139 * pkin(8) + t148;
t133 = pkin(2) * t157 - t138 * pkin(8) + t149;
t132 = t86 * pkin(3) + t85 * pkin(9) + t134;
t131 = t94 * pkin(3) + t93 * pkin(9) + t133;
t130 = t75 * pkin(4) + t74 * pkin(10) + t132;
t129 = t101 * pkin(2) - t140 * pkin(8) + t141;
t128 = t88 * pkin(4) + t87 * pkin(10) + t131;
t127 = t84 * pkin(3) + t83 * pkin(9) + t129;
t126 = t73 * pkin(4) + t72 * pkin(10) + t127;
t77 = t93 * t122 + t88 * t125;
t69 = t85 * t122 + t75 * t125;
t67 = t83 * t122 + t73 * t125;
t64 = -g(1) * t69 - g(2) * t67 - g(3) * t77;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t119 - g(2) * t115, t142, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t103 - g(2) * t101 - g(3) * t157, -g(1) * t102 - g(2) * t100 - g(3) * t154, -g(3) * t121 - t142 * t117, -g(1) * t148 - g(2) * t141 - g(3) * t149, 0, 0, 0, 0, 0, 0, -g(1) * t86 - g(2) * t84 - g(3) * t94, t135, -g(1) * t139 + g(2) * t140 + g(3) * t138, -g(1) * t134 - g(2) * t129 - g(3) * t133, 0, 0, 0, 0, 0, 0, -g(1) * t75 - g(2) * t73 - g(3) * t88, t136, -t135, -g(1) * t132 - g(2) * t127 - g(3) * t131, 0, 0, 0, 0, 0, 0, t64, t137, -t136, -g(1) * t130 - g(2) * t126 - g(3) * t128, 0, 0, 0, 0, 0, 0, t64, -t136, -t137, -g(1) * (t69 * pkin(5) + t68 * qJ(6) + t130) - g(2) * (t67 * pkin(5) + t66 * qJ(6) + t126) - g(3) * (t77 * pkin(5) + t76 * qJ(6) + t128);];
U_reg  = t1;
