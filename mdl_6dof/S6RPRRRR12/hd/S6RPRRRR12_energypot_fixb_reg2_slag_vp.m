% Calculate inertial parameters regressor of potential energy for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_energypot_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:00:40
% EndTime: 2019-03-09 08:00:40
% DurationCPUTime: 0.46s
% Computational Cost: add. (927->112), mult. (2553->180), div. (0->0), fcn. (3324->18), ass. (0->76)
t172 = cos(pkin(6));
t165 = sin(pkin(14));
t181 = cos(qJ(1));
t203 = t181 * t165;
t169 = cos(pkin(14));
t177 = sin(qJ(1));
t204 = t177 * t169;
t156 = -t172 * t204 - t203;
t167 = sin(pkin(7));
t171 = cos(pkin(7));
t168 = sin(pkin(6));
t209 = t168 * t177;
t148 = -t156 * t167 + t171 * t209;
t210 = t168 * t169;
t153 = -t167 * t210 + t172 * t171;
t176 = sin(qJ(3));
t180 = cos(qJ(3));
t207 = t169 * t171;
t211 = t167 * t172;
t145 = t180 * t211 + (-t165 * t176 + t180 * t207) * t168;
t166 = sin(pkin(8));
t170 = cos(pkin(8));
t132 = -t145 * t166 + t153 * t170;
t202 = t181 * t169;
t205 = t177 * t165;
t157 = -t172 * t205 + t202;
t193 = t156 * t171 + t167 * t209;
t137 = -t157 * t176 + t193 * t180;
t128 = -t137 * t166 + t148 * t170;
t155 = t172 * t203 + t204;
t154 = t172 * t202 - t205;
t208 = t168 * t181;
t194 = t154 * t171 - t167 * t208;
t135 = -t155 * t176 + t194 * t180;
t147 = -t154 * t167 - t171 * t208;
t127 = -t135 * t166 + t147 * t170;
t221 = cos(qJ(4));
t220 = t172 * qJ(2) + pkin(9);
t212 = t165 * t168;
t201 = t181 * pkin(1) + qJ(2) * t209;
t198 = t166 * t221;
t197 = t170 * t221;
t196 = g(1) * t177 - g(2) * t181;
t195 = t177 * pkin(1) - qJ(2) * t208;
t136 = t155 * t180 + t194 * t176;
t175 = sin(qJ(4));
t119 = t136 * t221 + (t135 * t170 + t147 * t166) * t175;
t174 = sin(qJ(5));
t179 = cos(qJ(5));
t110 = t119 * t174 - t127 * t179;
t138 = t157 * t180 + t193 * t176;
t121 = t138 * t221 + (t137 * t170 + t148 * t166) * t175;
t112 = t121 * t174 - t128 * t179;
t146 = t176 * t211 + (t165 * t180 + t176 * t207) * t168;
t124 = t146 * t221 + (t145 * t170 + t153 * t166) * t175;
t114 = t124 * t174 - t132 * t179;
t192 = g(1) * t112 + g(2) * t110 + g(3) * t114;
t118 = -t135 * t197 + t136 * t175 - t147 * t198;
t120 = -t137 * t197 + t138 * t175 - t148 * t198;
t123 = -t145 * t197 + t146 * t175 - t153 * t198;
t191 = g(1) * t120 + g(2) * t118 + g(3) * t123;
t190 = t157 * pkin(2) + t148 * pkin(10) + t201;
t189 = pkin(2) * t212 + t153 * pkin(10) + t220;
t188 = t138 * pkin(3) + t128 * pkin(11) + t190;
t187 = t146 * pkin(3) + t132 * pkin(11) + t189;
t186 = t155 * pkin(2) + t147 * pkin(10) + t195;
t185 = t121 * pkin(4) + t120 * pkin(12) + t188;
t184 = t124 * pkin(4) + t123 * pkin(12) + t187;
t183 = t136 * pkin(3) + t127 * pkin(11) + t186;
t182 = t119 * pkin(4) + t118 * pkin(12) + t183;
t178 = cos(qJ(6));
t173 = sin(qJ(6));
t115 = t124 * t179 + t132 * t174;
t113 = t121 * t179 + t128 * t174;
t111 = t119 * t179 + t127 * t174;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t181 - g(2) * t177, t196, -g(3), -g(3) * pkin(9), 0, 0, 0, 0, 0, 0, -g(1) * t157 - g(2) * t155 - g(3) * t212, -g(1) * t156 - g(2) * t154 - g(3) * t210, -g(3) * t172 - t196 * t168, -g(1) * t201 - g(2) * t195 - g(3) * t220, 0, 0, 0, 0, 0, 0, -g(1) * t138 - g(2) * t136 - g(3) * t146, -g(1) * t137 - g(2) * t135 - g(3) * t145, -g(1) * t148 - g(2) * t147 - g(3) * t153, -g(1) * t190 - g(2) * t186 - g(3) * t189, 0, 0, 0, 0, 0, 0, -g(1) * t121 - g(2) * t119 - g(3) * t124, t191, -g(1) * t128 - g(2) * t127 - g(3) * t132, -g(1) * t188 - g(2) * t183 - g(3) * t187, 0, 0, 0, 0, 0, 0, -g(1) * t113 - g(2) * t111 - g(3) * t115, t192, -t191, -g(1) * t185 - g(2) * t182 - g(3) * t184, 0, 0, 0, 0, 0, 0, -g(1) * (t113 * t178 + t120 * t173) - g(2) * (t111 * t178 + t118 * t173) - g(3) * (t115 * t178 + t123 * t173) -g(1) * (-t113 * t173 + t120 * t178) - g(2) * (-t111 * t173 + t118 * t178) - g(3) * (-t115 * t173 + t123 * t178) -t192, -g(1) * (t113 * pkin(5) + t112 * pkin(13) + t185) - g(2) * (t111 * pkin(5) + t110 * pkin(13) + t182) - g(3) * (t115 * pkin(5) + t114 * pkin(13) + t184);];
U_reg  = t1;
