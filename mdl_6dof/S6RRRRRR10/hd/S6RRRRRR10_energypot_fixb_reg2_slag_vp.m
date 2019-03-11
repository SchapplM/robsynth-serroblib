% Calculate inertial parameters regressor of potential energy for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energypot_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 06:16:47
% EndTime: 2019-03-10 06:16:47
% DurationCPUTime: 0.46s
% Computational Cost: add. (927->112), mult. (2553->180), div. (0->0), fcn. (3324->18), ass. (0->76)
t166 = cos(pkin(6));
t171 = sin(qJ(2));
t177 = cos(qJ(1));
t199 = t177 * t171;
t172 = sin(qJ(1));
t176 = cos(qJ(2));
t200 = t172 * t176;
t152 = -t166 * t200 - t199;
t162 = sin(pkin(7));
t165 = cos(pkin(7));
t163 = sin(pkin(6));
t206 = t163 * t172;
t144 = -t152 * t162 + t165 * t206;
t205 = t163 * t176;
t149 = -t162 * t205 + t166 * t165;
t170 = sin(qJ(3));
t175 = cos(qJ(3));
t203 = t165 * t176;
t208 = t162 * t166;
t141 = t175 * t208 + (-t170 * t171 + t175 * t203) * t163;
t161 = sin(pkin(8));
t164 = cos(pkin(8));
t128 = -t141 * t161 + t149 * t164;
t198 = t177 * t176;
t201 = t172 * t171;
t153 = -t166 * t201 + t198;
t189 = t152 * t165 + t162 * t206;
t133 = -t153 * t170 + t189 * t175;
t124 = -t133 * t161 + t144 * t164;
t151 = t166 * t199 + t200;
t150 = t166 * t198 - t201;
t204 = t163 * t177;
t190 = t150 * t165 - t162 * t204;
t131 = -t151 * t170 + t190 * t175;
t143 = -t150 * t162 - t165 * t204;
t123 = -t131 * t161 + t143 * t164;
t217 = cos(qJ(4));
t216 = t166 * pkin(10) + pkin(9);
t207 = t163 * t171;
t197 = t177 * pkin(1) + pkin(10) * t206;
t194 = t161 * t217;
t193 = t164 * t217;
t192 = t172 * pkin(1) - pkin(10) * t204;
t191 = g(1) * t172 - g(2) * t177;
t132 = t151 * t175 + t190 * t170;
t169 = sin(qJ(4));
t115 = t132 * t217 + (t131 * t164 + t143 * t161) * t169;
t168 = sin(qJ(5));
t174 = cos(qJ(5));
t106 = t115 * t168 - t123 * t174;
t134 = t153 * t175 + t189 * t170;
t117 = t134 * t217 + (t133 * t164 + t144 * t161) * t169;
t108 = t117 * t168 - t124 * t174;
t142 = t170 * t208 + (t170 * t203 + t171 * t175) * t163;
t120 = t142 * t217 + (t141 * t164 + t149 * t161) * t169;
t110 = t120 * t168 - t128 * t174;
t188 = g(1) * t108 + g(2) * t106 + g(3) * t110;
t114 = -t131 * t193 + t132 * t169 - t143 * t194;
t116 = -t133 * t193 + t134 * t169 - t144 * t194;
t119 = -t141 * t193 + t142 * t169 - t149 * t194;
t187 = g(1) * t116 + g(2) * t114 + g(3) * t119;
t186 = t153 * pkin(2) + t144 * pkin(11) + t197;
t185 = pkin(2) * t207 + t149 * pkin(11) + t216;
t184 = t134 * pkin(3) + t124 * pkin(12) + t186;
t183 = t151 * pkin(2) + t143 * pkin(11) + t192;
t182 = t142 * pkin(3) + t128 * pkin(12) + t185;
t181 = t117 * pkin(4) + t116 * pkin(13) + t184;
t180 = t120 * pkin(4) + t119 * pkin(13) + t182;
t179 = t132 * pkin(3) + t123 * pkin(12) + t183;
t178 = t115 * pkin(4) + t114 * pkin(13) + t179;
t173 = cos(qJ(6));
t167 = sin(qJ(6));
t111 = t120 * t174 + t128 * t168;
t109 = t117 * t174 + t124 * t168;
t107 = t115 * t174 + t123 * t168;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t177 - g(2) * t172, t191, -g(3), -g(3) * pkin(9), 0, 0, 0, 0, 0, 0, -g(1) * t153 - g(2) * t151 - g(3) * t207, -g(1) * t152 - g(2) * t150 - g(3) * t205, -g(3) * t166 - t191 * t163, -g(1) * t197 - g(2) * t192 - g(3) * t216, 0, 0, 0, 0, 0, 0, -g(1) * t134 - g(2) * t132 - g(3) * t142, -g(1) * t133 - g(2) * t131 - g(3) * t141, -g(1) * t144 - g(2) * t143 - g(3) * t149, -g(1) * t186 - g(2) * t183 - g(3) * t185, 0, 0, 0, 0, 0, 0, -g(1) * t117 - g(2) * t115 - g(3) * t120, t187, -g(1) * t124 - g(2) * t123 - g(3) * t128, -g(1) * t184 - g(2) * t179 - g(3) * t182, 0, 0, 0, 0, 0, 0, -g(1) * t109 - g(2) * t107 - g(3) * t111, t188, -t187, -g(1) * t181 - g(2) * t178 - g(3) * t180, 0, 0, 0, 0, 0, 0, -g(1) * (t109 * t173 + t116 * t167) - g(2) * (t107 * t173 + t114 * t167) - g(3) * (t111 * t173 + t119 * t167) -g(1) * (-t109 * t167 + t116 * t173) - g(2) * (-t107 * t167 + t114 * t173) - g(3) * (-t111 * t167 + t119 * t173) -t188, -g(1) * (t109 * pkin(5) + t108 * pkin(14) + t181) - g(2) * (t107 * pkin(5) + t106 * pkin(14) + t178) - g(3) * (t111 * pkin(5) + t110 * pkin(14) + t180);];
U_reg  = t1;
