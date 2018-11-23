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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
% StartTime: 2018-11-23 11:16:22
% EndTime: 2018-11-23 11:16:22
% DurationCPUTime: 0.59s
% Computational Cost: add. (3285->132), mult. (3339->187), div. (0->0), fcn. (3324->30), ass. (0->94)
t182 = pkin(6) - qJ(2);
t171 = cos(t182) / 0.2e1;
t181 = pkin(6) + qJ(2);
t175 = cos(t181);
t163 = t171 + t175 / 0.2e1;
t193 = sin(qJ(2));
t194 = sin(qJ(1));
t200 = cos(qJ(1));
t152 = -t194 * t163 - t200 * t193;
t184 = sin(pkin(7));
t187 = cos(pkin(7));
t185 = sin(pkin(6));
t226 = t185 * t194;
t144 = -t152 * t184 + t187 * t226;
t169 = sin(t181) / 0.2e1;
t173 = sin(t182);
t158 = t169 + t173 / 0.2e1;
t188 = cos(pkin(6));
t149 = -t158 * t184 + t188 * t187;
t179 = pkin(7) + qJ(3);
t168 = sin(t179) / 0.2e1;
t180 = pkin(7) - qJ(3);
t172 = sin(t180);
t156 = t168 + t172 / 0.2e1;
t170 = cos(t180) / 0.2e1;
t174 = cos(t179);
t161 = t170 + t174 / 0.2e1;
t164 = t171 - t175 / 0.2e1;
t192 = sin(qJ(3));
t137 = t188 * t156 + t158 * t161 - t164 * t192;
t183 = sin(pkin(8));
t186 = cos(pkin(8));
t128 = -t137 * t183 + t149 * t186;
t159 = t169 - t173 / 0.2e1;
t199 = cos(qJ(2));
t153 = -t194 * t159 + t200 * t199;
t133 = t152 * t161 - t153 * t192 + t156 * t226;
t124 = -t133 * t183 + t144 * t186;
t150 = t200 * t163 - t194 * t193;
t151 = t200 * t159 + t194 * t199;
t225 = t185 * t200;
t131 = t150 * t161 - t151 * t192 - t156 * t225;
t143 = -t150 * t184 - t187 * t225;
t123 = -t131 * t183 + t143 * t186;
t235 = t188 * pkin(10) + pkin(9);
t223 = t200 * pkin(1) + pkin(10) * t226;
t222 = pkin(8) - qJ(4);
t221 = pkin(8) + qJ(4);
t219 = cos(t221);
t218 = sin(t222);
t217 = t194 * pkin(1) - pkin(10) * t225;
t216 = cos(t222) / 0.2e1;
t215 = sin(t221) / 0.2e1;
t214 = g(1) * t194 - g(2) * t200;
t157 = t168 - t172 / 0.2e1;
t162 = t170 - t174 / 0.2e1;
t198 = cos(qJ(3));
t132 = t150 * t157 + t151 * t198 - t162 * t225;
t155 = t215 - t218 / 0.2e1;
t160 = t216 - t219 / 0.2e1;
t197 = cos(qJ(4));
t115 = t131 * t155 + t132 * t197 + t143 * t160;
t190 = sin(qJ(5));
t196 = cos(qJ(5));
t106 = t115 * t190 - t123 * t196;
t134 = t152 * t157 + t153 * t198 + t162 * t226;
t117 = t133 * t155 + t134 * t197 + t144 * t160;
t108 = t117 * t190 - t124 * t196;
t138 = t158 * t157 + t188 * t162 + t164 * t198;
t120 = t137 * t155 + t138 * t197 + t149 * t160;
t110 = t120 * t190 - t128 * t196;
t213 = g(1) * t108 + g(2) * t106 + g(3) * t110;
t191 = sin(qJ(4));
t208 = t215 + t218 / 0.2e1;
t209 = t216 + t219 / 0.2e1;
t114 = -t131 * t209 + t132 * t191 - t143 * t208;
t116 = -t133 * t209 + t134 * t191 - t144 * t208;
t119 = -t137 * t209 + t138 * t191 - t149 * t208;
t212 = g(1) * t116 + g(2) * t114 + g(3) * t119;
t211 = t164 * pkin(2) + t149 * pkin(11) + t235;
t210 = t153 * pkin(2) + t144 * pkin(11) + t223;
t207 = t138 * pkin(3) + t128 * pkin(12) + t211;
t206 = t134 * pkin(3) + t124 * pkin(12) + t210;
t205 = t151 * pkin(2) + t143 * pkin(11) + t217;
t204 = t120 * pkin(4) + t119 * pkin(13) + t207;
t203 = t117 * pkin(4) + t116 * pkin(13) + t206;
t202 = t132 * pkin(3) + t123 * pkin(12) + t205;
t201 = t115 * pkin(4) + t114 * pkin(13) + t202;
t195 = cos(qJ(6));
t189 = sin(qJ(6));
t111 = t120 * t196 + t128 * t190;
t109 = t117 * t196 + t124 * t190;
t107 = t115 * t196 + t123 * t190;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t200 - g(2) * t194, t214, -g(3), -g(3) * pkin(9), 0, 0, 0, 0, 0, 0, -g(1) * t153 - g(2) * t151 - g(3) * t164, -g(1) * t152 - g(2) * t150 - g(3) * t158, -g(3) * t188 - t214 * t185, -g(1) * t223 - g(2) * t217 - g(3) * t235, 0, 0, 0, 0, 0, 0, -g(1) * t134 - g(2) * t132 - g(3) * t138, -g(1) * t133 - g(2) * t131 - g(3) * t137, -g(1) * t144 - g(2) * t143 - g(3) * t149, -g(1) * t210 - g(2) * t205 - g(3) * t211, 0, 0, 0, 0, 0, 0, -g(1) * t117 - g(2) * t115 - g(3) * t120, t212, -g(1) * t124 - g(2) * t123 - g(3) * t128, -g(1) * t206 - g(2) * t202 - g(3) * t207, 0, 0, 0, 0, 0, 0, -g(1) * t109 - g(2) * t107 - g(3) * t111, t213, -t212, -g(1) * t203 - g(2) * t201 - g(3) * t204, 0, 0, 0, 0, 0, 0, -g(1) * (t109 * t195 + t116 * t189) - g(2) * (t107 * t195 + t114 * t189) - g(3) * (t111 * t195 + t119 * t189) -g(1) * (-t109 * t189 + t116 * t195) - g(2) * (-t107 * t189 + t114 * t195) - g(3) * (-t111 * t189 + t119 * t195) -t213, -g(1) * (t109 * pkin(5) + t108 * pkin(14) + t203) - g(2) * (t107 * pkin(5) + t106 * pkin(14) + t201) - g(3) * (t111 * pkin(5) + t110 * pkin(14) + t204);];
U_reg  = t1;
