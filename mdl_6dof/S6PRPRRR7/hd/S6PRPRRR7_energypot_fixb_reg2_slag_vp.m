% Calculate inertial parameters regressor of potential energy for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_energypot_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:55:52
% EndTime: 2019-03-08 20:55:53
% DurationCPUTime: 0.46s
% Computational Cost: add. (927->112), mult. (2553->182), div. (0->0), fcn. (3324->18), ass. (0->75)
t151 = sin(pkin(7));
t156 = cos(pkin(7));
t157 = cos(pkin(6));
t152 = sin(pkin(6));
t164 = cos(qJ(2));
t190 = t152 * t164;
t136 = -t151 * t190 + t156 * t157;
t149 = sin(pkin(13));
t154 = cos(pkin(13));
t161 = sin(qJ(2));
t186 = t157 * t164;
t139 = -t149 * t186 - t154 * t161;
t192 = t152 * t156;
t131 = -t139 * t151 + t149 * t192;
t148 = sin(pkin(14));
t153 = cos(pkin(14));
t189 = t156 * t164;
t194 = t151 * t157;
t128 = t153 * t194 + (-t148 * t161 + t153 * t189) * t152;
t150 = sin(pkin(8));
t155 = cos(pkin(8));
t117 = -t128 * t150 + t136 * t155;
t187 = t157 * t161;
t140 = -t149 * t187 + t154 * t164;
t195 = t149 * t152;
t176 = t139 * t156 + t151 * t195;
t120 = -t140 * t148 + t153 * t176;
t111 = -t120 * t150 + t131 * t155;
t138 = t149 * t164 + t154 * t187;
t137 = -t149 * t161 + t154 * t186;
t193 = t152 * t154;
t177 = t137 * t156 - t151 * t193;
t118 = -t138 * t148 + t153 * t177;
t130 = -t137 * t151 - t154 * t192;
t110 = -t118 * t150 + t130 * t155;
t203 = cos(qJ(4));
t191 = t152 * t161;
t185 = pkin(1) * t154 + pkin(9) * t195;
t184 = pkin(9) * t157 + qJ(1);
t181 = t150 * t203;
t180 = t155 * t203;
t179 = pkin(1) * t149 - pkin(9) * t193;
t178 = g(1) * t149 - g(2) * t154;
t119 = t138 * t153 + t148 * t177;
t160 = sin(qJ(4));
t102 = t119 * t203 + (t118 * t155 + t130 * t150) * t160;
t159 = sin(qJ(5));
t163 = cos(qJ(5));
t93 = t102 * t159 - t110 * t163;
t121 = t140 * t153 + t148 * t176;
t104 = t121 * t203 + (t120 * t155 + t131 * t150) * t160;
t95 = t104 * t159 - t111 * t163;
t129 = t153 * t191 + (t152 * t189 + t194) * t148;
t109 = t129 * t203 + (t128 * t155 + t136 * t150) * t160;
t97 = t109 * t159 - t117 * t163;
t175 = g(1) * t95 + g(2) * t93 + g(3) * t97;
t101 = -t118 * t180 + t119 * t160 - t130 * t181;
t103 = -t120 * t180 + t121 * t160 - t131 * t181;
t108 = -t128 * t180 + t129 * t160 - t136 * t181;
t174 = g(1) * t103 + g(2) * t101 + g(3) * t108;
t173 = t140 * pkin(2) + qJ(3) * t131 + t185;
t172 = pkin(2) * t191 + qJ(3) * t136 + t184;
t171 = t121 * pkin(3) + pkin(10) * t111 + t173;
t170 = pkin(2) * t138 + qJ(3) * t130 + t179;
t169 = t129 * pkin(3) + pkin(10) * t117 + t172;
t168 = pkin(4) * t104 + pkin(11) * t103 + t171;
t167 = pkin(4) * t109 + pkin(11) * t108 + t169;
t166 = t119 * pkin(3) + pkin(10) * t110 + t170;
t165 = pkin(4) * t102 + t101 * pkin(11) + t166;
t162 = cos(qJ(6));
t158 = sin(qJ(6));
t98 = t109 * t163 + t117 * t159;
t96 = t104 * t163 + t111 * t159;
t94 = t102 * t163 + t110 * t159;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t154 - g(2) * t149, t178, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t140 - g(2) * t138 - g(3) * t191, -g(1) * t139 - g(2) * t137 - g(3) * t190, -g(3) * t157 - t152 * t178, -g(1) * t185 - g(2) * t179 - g(3) * t184, 0, 0, 0, 0, 0, 0, -g(1) * t121 - g(2) * t119 - g(3) * t129, -g(1) * t120 - g(2) * t118 - g(3) * t128, -g(1) * t131 - g(2) * t130 - g(3) * t136, -g(1) * t173 - g(2) * t170 - g(3) * t172, 0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t102 - g(3) * t109, t174, -g(1) * t111 - g(2) * t110 - g(3) * t117, -g(1) * t171 - g(2) * t166 - g(3) * t169, 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t98, t175, -t174, -g(1) * t168 - g(2) * t165 - g(3) * t167, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t158 + t162 * t96) - g(2) * (t101 * t158 + t162 * t94) - g(3) * (t108 * t158 + t162 * t98) -g(1) * (t103 * t162 - t158 * t96) - g(2) * (t101 * t162 - t158 * t94) - g(3) * (t108 * t162 - t158 * t98) -t175, -g(1) * (pkin(5) * t96 + pkin(12) * t95 + t168) - g(2) * (t94 * pkin(5) + t93 * pkin(12) + t165) - g(3) * (pkin(5) * t98 + pkin(12) * t97 + t167);];
U_reg  = t1;
