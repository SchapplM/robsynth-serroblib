% Calculate inertial parameters regressor of potential energy for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_energypot_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:25:02
% EndTime: 2019-03-09 01:25:02
% DurationCPUTime: 0.47s
% Computational Cost: add. (927->112), mult. (2553->183), div. (0->0), fcn. (3324->18), ass. (0->75)
t149 = sin(pkin(7));
t153 = cos(pkin(7));
t154 = cos(pkin(6));
t150 = sin(pkin(6));
t163 = cos(qJ(2));
t189 = t150 * t163;
t135 = -t149 * t189 + t154 * t153;
t147 = sin(pkin(14));
t151 = cos(pkin(14));
t159 = sin(qJ(2));
t185 = t154 * t163;
t138 = -t147 * t185 - t151 * t159;
t191 = t150 * t153;
t130 = -t138 * t149 + t147 * t191;
t158 = sin(qJ(3));
t162 = cos(qJ(3));
t188 = t153 * t163;
t193 = t149 * t154;
t127 = t162 * t193 + (-t158 * t159 + t162 * t188) * t150;
t148 = sin(pkin(8));
t152 = cos(pkin(8));
t116 = -t127 * t148 + t135 * t152;
t186 = t154 * t159;
t139 = -t147 * t186 + t151 * t163;
t194 = t147 * t150;
t175 = t138 * t153 + t149 * t194;
t119 = -t139 * t158 + t175 * t162;
t110 = -t119 * t148 + t130 * t152;
t137 = t147 * t163 + t151 * t186;
t136 = -t147 * t159 + t151 * t185;
t192 = t150 * t151;
t176 = t136 * t153 - t149 * t192;
t117 = -t137 * t158 + t176 * t162;
t129 = -t136 * t149 - t151 * t191;
t109 = -t117 * t148 + t129 * t152;
t202 = cos(qJ(4));
t190 = t150 * t159;
t184 = t151 * pkin(1) + pkin(9) * t194;
t183 = t154 * pkin(9) + qJ(1);
t180 = t148 * t202;
t179 = t152 * t202;
t178 = t147 * pkin(1) - pkin(9) * t192;
t177 = g(1) * t147 - g(2) * t151;
t118 = t137 * t162 + t176 * t158;
t157 = sin(qJ(4));
t101 = t118 * t202 + (t117 * t152 + t129 * t148) * t157;
t156 = sin(qJ(5));
t161 = cos(qJ(5));
t92 = t101 * t156 - t109 * t161;
t120 = t139 * t162 + t175 * t158;
t103 = t120 * t202 + (t119 * t152 + t130 * t148) * t157;
t94 = t103 * t156 - t110 * t161;
t128 = t158 * t193 + (t158 * t188 + t159 * t162) * t150;
t108 = t128 * t202 + (t127 * t152 + t135 * t148) * t157;
t98 = t108 * t156 - t116 * t161;
t174 = g(1) * t94 + g(2) * t92 + g(3) * t98;
t100 = -t117 * t179 + t118 * t157 - t129 * t180;
t102 = -t119 * t179 + t120 * t157 - t130 * t180;
t107 = -t127 * t179 + t128 * t157 - t135 * t180;
t173 = g(1) * t102 + g(2) * t100 + g(3) * t107;
t172 = t139 * pkin(2) + t130 * pkin(10) + t184;
t171 = pkin(2) * t190 + t135 * pkin(10) + t183;
t170 = t120 * pkin(3) + t110 * pkin(11) + t172;
t169 = t137 * pkin(2) + t129 * pkin(10) + t178;
t168 = t128 * pkin(3) + t116 * pkin(11) + t171;
t167 = t103 * pkin(4) + t102 * pkin(12) + t170;
t166 = t108 * pkin(4) + t107 * pkin(12) + t168;
t165 = t118 * pkin(3) + t109 * pkin(11) + t169;
t164 = t101 * pkin(4) + t100 * pkin(12) + t165;
t160 = cos(qJ(6));
t155 = sin(qJ(6));
t99 = t108 * t161 + t116 * t156;
t95 = t103 * t161 + t110 * t156;
t93 = t101 * t161 + t109 * t156;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t151 - g(2) * t147, t177, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t139 - g(2) * t137 - g(3) * t190, -g(1) * t138 - g(2) * t136 - g(3) * t189, -g(3) * t154 - t177 * t150, -g(1) * t184 - g(2) * t178 - g(3) * t183, 0, 0, 0, 0, 0, 0, -g(1) * t120 - g(2) * t118 - g(3) * t128, -g(1) * t119 - g(2) * t117 - g(3) * t127, -g(1) * t130 - g(2) * t129 - g(3) * t135, -g(1) * t172 - g(2) * t169 - g(3) * t171, 0, 0, 0, 0, 0, 0, -g(1) * t103 - g(2) * t101 - g(3) * t108, t173, -g(1) * t110 - g(2) * t109 - g(3) * t116, -g(1) * t170 - g(2) * t165 - g(3) * t168, 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t99, t174, -t173, -g(1) * t167 - g(2) * t164 - g(3) * t166, 0, 0, 0, 0, 0, 0, -g(1) * (t102 * t155 + t95 * t160) - g(2) * (t100 * t155 + t93 * t160) - g(3) * (t107 * t155 + t99 * t160) -g(1) * (t102 * t160 - t95 * t155) - g(2) * (t100 * t160 - t93 * t155) - g(3) * (t107 * t160 - t99 * t155) -t174, -g(1) * (t95 * pkin(5) + t94 * pkin(13) + t167) - g(2) * (t93 * pkin(5) + t92 * pkin(13) + t164) - g(3) * (t99 * pkin(5) + t98 * pkin(13) + t166);];
U_reg  = t1;
