% Calculate inertial parameters regressor of potential energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:26
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR14_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energypot_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:20:39
% EndTime: 2019-01-03 10:20:40
% DurationCPUTime: 0.48s
% Computational Cost: add. (927->112), mult. (2553->179), div. (0->0), fcn. (3324->18), ass. (0->76)
t168 = cos(pkin(6));
t172 = sin(qJ(2));
t177 = cos(qJ(1));
t199 = t177 * t172;
t173 = sin(qJ(1));
t176 = cos(qJ(2));
t200 = t173 * t176;
t152 = -t168 * t200 - t199;
t163 = sin(pkin(7));
t167 = cos(pkin(7));
t164 = sin(pkin(6));
t206 = t164 * t173;
t144 = -t152 * t163 + t167 * t206;
t205 = t164 * t176;
t149 = -t163 * t205 + t168 * t167;
t161 = sin(pkin(14));
t165 = cos(pkin(14));
t203 = t167 * t176;
t208 = t163 * t168;
t141 = t165 * t208 + (-t161 * t172 + t165 * t203) * t164;
t162 = sin(pkin(8));
t166 = cos(pkin(8));
t128 = -t141 * t162 + t149 * t166;
t198 = t177 * t176;
t201 = t173 * t172;
t153 = -t168 * t201 + t198;
t189 = t152 * t167 + t163 * t206;
t133 = -t153 * t161 + t189 * t165;
t124 = -t133 * t162 + t144 * t166;
t151 = t168 * t199 + t200;
t150 = t168 * t198 - t201;
t204 = t164 * t177;
t190 = t150 * t167 - t163 * t204;
t131 = -t151 * t161 + t190 * t165;
t143 = -t150 * t163 - t167 * t204;
t123 = -t131 * t162 + t143 * t166;
t217 = cos(qJ(4));
t216 = t168 * pkin(10) + pkin(9);
t207 = t164 * t172;
t197 = t177 * pkin(1) + pkin(10) * t206;
t194 = t162 * t217;
t193 = t166 * t217;
t192 = t173 * pkin(1) - pkin(10) * t204;
t191 = g(1) * t173 - g(2) * t177;
t132 = t151 * t165 + t190 * t161;
t171 = sin(qJ(4));
t115 = t132 * t217 + (t131 * t166 + t143 * t162) * t171;
t170 = sin(qJ(5));
t175 = cos(qJ(5));
t106 = t115 * t170 - t123 * t175;
t134 = t153 * t165 + t189 * t161;
t117 = t134 * t217 + (t133 * t166 + t144 * t162) * t171;
t108 = t117 * t170 - t124 * t175;
t142 = t165 * t207 + (t164 * t203 + t208) * t161;
t120 = t142 * t217 + (t141 * t166 + t149 * t162) * t171;
t110 = t120 * t170 - t128 * t175;
t188 = g(1) * t108 + g(2) * t106 + g(3) * t110;
t114 = -t131 * t193 + t132 * t171 - t143 * t194;
t116 = -t133 * t193 + t134 * t171 - t144 * t194;
t119 = -t141 * t193 + t142 * t171 - t149 * t194;
t187 = g(1) * t116 + g(2) * t114 + g(3) * t119;
t186 = t153 * pkin(2) + t144 * qJ(3) + t197;
t185 = pkin(2) * t207 + t149 * qJ(3) + t216;
t184 = t134 * pkin(3) + t124 * pkin(11) + t186;
t183 = t151 * pkin(2) + t143 * qJ(3) + t192;
t182 = t142 * pkin(3) + t128 * pkin(11) + t185;
t181 = t117 * pkin(4) + t116 * pkin(12) + t184;
t180 = t120 * pkin(4) + t119 * pkin(12) + t182;
t179 = t132 * pkin(3) + t123 * pkin(11) + t183;
t178 = t115 * pkin(4) + t114 * pkin(12) + t179;
t174 = cos(qJ(6));
t169 = sin(qJ(6));
t111 = t120 * t175 + t128 * t170;
t109 = t117 * t175 + t124 * t170;
t107 = t115 * t175 + t123 * t170;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t177 - g(2) * t173, t191, -g(3), -g(3) * pkin(9), 0, 0, 0, 0, 0, 0, -g(1) * t153 - g(2) * t151 - g(3) * t207, -g(1) * t152 - g(2) * t150 - g(3) * t205, -g(3) * t168 - t191 * t164, -g(1) * t197 - g(2) * t192 - g(3) * t216, 0, 0, 0, 0, 0, 0, -g(1) * t134 - g(2) * t132 - g(3) * t142, -g(1) * t133 - g(2) * t131 - g(3) * t141, -g(1) * t144 - g(2) * t143 - g(3) * t149, -g(1) * t186 - g(2) * t183 - g(3) * t185, 0, 0, 0, 0, 0, 0, -g(1) * t117 - g(2) * t115 - g(3) * t120, t187, -g(1) * t124 - g(2) * t123 - g(3) * t128, -g(1) * t184 - g(2) * t179 - g(3) * t182, 0, 0, 0, 0, 0, 0, -g(1) * t109 - g(2) * t107 - g(3) * t111, t188, -t187, -g(1) * t181 - g(2) * t178 - g(3) * t180, 0, 0, 0, 0, 0, 0, -g(1) * (t109 * t174 + t116 * t169) - g(2) * (t107 * t174 + t114 * t169) - g(3) * (t111 * t174 + t119 * t169) -g(1) * (-t109 * t169 + t116 * t174) - g(2) * (-t107 * t169 + t114 * t174) - g(3) * (-t111 * t169 + t119 * t174) -t188, -g(1) * (t109 * pkin(5) + t108 * pkin(13) + t181) - g(2) * (t107 * pkin(5) + t106 * pkin(13) + t178) - g(3) * (t111 * pkin(5) + t110 * pkin(13) + t180);];
U_reg  = t1;
