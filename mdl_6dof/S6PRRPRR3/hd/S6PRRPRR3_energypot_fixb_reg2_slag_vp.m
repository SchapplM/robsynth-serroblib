% Calculate inertial parameters regressor of potential energy for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:07:42
% EndTime: 2019-03-08 22:07:42
% DurationCPUTime: 0.35s
% Computational Cost: add. (574->109), mult. (1452->172), div. (0->0), fcn. (1852->16), ass. (0->71)
t128 = sin(pkin(12));
t132 = cos(pkin(12));
t138 = sin(qJ(2));
t134 = cos(pkin(6));
t142 = cos(qJ(2));
t160 = t134 * t142;
t115 = -t128 * t138 + t132 * t160;
t129 = sin(pkin(7));
t130 = sin(pkin(6));
t133 = cos(pkin(7));
t164 = t130 * t133;
t100 = -t115 * t129 - t132 * t164;
t117 = -t128 * t160 - t132 * t138;
t101 = -t117 * t129 + t128 * t164;
t162 = t130 * t142;
t114 = -t129 * t162 + t134 * t133;
t92 = -g(1) * t101 - g(2) * t100 - g(3) * t114;
t137 = sin(qJ(3));
t171 = pkin(3) * t137;
t167 = pkin(9) + qJ(4);
t166 = t128 * t130;
t165 = t130 * t132;
t163 = t130 * t138;
t161 = t134 * t138;
t159 = t132 * pkin(1) + pkin(8) * t166;
t158 = t134 * pkin(8) + qJ(1);
t124 = t128 * pkin(1);
t157 = -pkin(8) * t165 + t124;
t156 = g(1) * t128 - g(2) * t132;
t112 = t129 * t171 + t167 * t133;
t113 = -t167 * t129 + t133 * t171;
t118 = -t128 * t161 + t132 * t142;
t141 = cos(qJ(3));
t123 = t141 * pkin(3) + pkin(2);
t155 = t112 * t166 + t117 * t113 + t118 * t123 + t159;
t154 = t134 * t112 + t113 * t162 + t123 * t163 + t158;
t127 = sin(pkin(13));
t131 = cos(pkin(13));
t153 = t141 * t127 + t137 * t131;
t120 = -t137 * t127 + t141 * t131;
t136 = sin(qJ(5));
t140 = cos(qJ(5));
t109 = t153 * t129;
t111 = t153 * t133;
t116 = t128 * t142 + t132 * t161;
t87 = -t109 * t165 + t115 * t111 + t116 * t120;
t80 = -t100 * t140 + t87 * t136;
t89 = t109 * t166 + t117 * t111 + t118 * t120;
t82 = -t101 * t140 + t89 * t136;
t95 = t134 * t109 + (t111 * t142 + t120 * t138) * t130;
t90 = -t114 * t140 + t95 * t136;
t152 = g(1) * t82 + g(2) * t80 + g(3) * t90;
t110 = t120 * t133;
t150 = t120 * t129;
t149 = t130 * t150;
t86 = t115 * t110 - t116 * t153 - t132 * t149;
t88 = t117 * t110 - t118 * t153 + t128 * t149;
t94 = (t142 * t110 - t138 * t153) * t130 + t134 * t150;
t151 = g(1) * t88 + g(2) * t86 + g(3) * t94;
t148 = t116 * t123 + t124 + t115 * t113 + (-pkin(8) - t112) * t165;
t147 = g(1) * t118 + g(2) * t116 + g(3) * t163;
t146 = t89 * pkin(4) - t88 * pkin(10) + t155;
t145 = t95 * pkin(4) - t94 * pkin(10) + t154;
t144 = t87 * pkin(4) - t86 * pkin(10) + t148;
t143 = -g(1) * (t117 * t133 + t129 * t166) - g(2) * (t115 * t133 - t129 * t165) - g(3) * (t129 * t134 + t133 * t162);
t139 = cos(qJ(6));
t135 = sin(qJ(6));
t91 = t114 * t136 + t95 * t140;
t83 = t101 * t136 + t89 * t140;
t81 = t100 * t136 + t87 * t140;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t132 - g(2) * t128, t156, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t147, -g(1) * t117 - g(2) * t115 - g(3) * t162, -g(3) * t134 - t156 * t130, -g(1) * t159 - g(2) * t157 - g(3) * t158, 0, 0, 0, 0, 0, 0, t143 * t137 - t147 * t141, t147 * t137 + t143 * t141, t92, -g(1) * (t118 * pkin(2) + t159) - g(2) * (t116 * pkin(2) + t157) - g(3) * (pkin(2) * t163 + t158) + t92 * pkin(9), 0, 0, 0, 0, 0, 0, -g(1) * t89 - g(2) * t87 - g(3) * t95, -t151, t92, -g(1) * t155 - g(2) * t148 - g(3) * t154, 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t91, t152, t151, -g(1) * t146 - g(2) * t144 - g(3) * t145, 0, 0, 0, 0, 0, 0, -g(1) * (-t88 * t135 + t83 * t139) - g(2) * (-t86 * t135 + t81 * t139) - g(3) * (-t94 * t135 + t91 * t139) -g(1) * (-t83 * t135 - t88 * t139) - g(2) * (-t81 * t135 - t86 * t139) - g(3) * (-t91 * t135 - t94 * t139) -t152, -g(1) * (t83 * pkin(5) + t82 * pkin(11) + t146) - g(2) * (t81 * pkin(5) + t80 * pkin(11) + t144) - g(3) * (t91 * pkin(5) + t90 * pkin(11) + t145);];
U_reg  = t1;
