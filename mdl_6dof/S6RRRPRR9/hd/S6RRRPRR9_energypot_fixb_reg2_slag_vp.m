% Calculate inertial parameters regressor of potential energy for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:10:15
% EndTime: 2019-03-09 19:10:15
% DurationCPUTime: 0.35s
% Computational Cost: add. (574->109), mult. (1452->169), div. (0->0), fcn. (1852->16), ass. (0->72)
t143 = cos(pkin(6));
t152 = cos(qJ(2));
t153 = cos(qJ(1));
t170 = t152 * t153;
t147 = sin(qJ(2));
t148 = sin(qJ(1));
t172 = t148 * t147;
t126 = t143 * t170 - t172;
t139 = sin(pkin(7));
t142 = cos(pkin(7));
t140 = sin(pkin(6));
t175 = t140 * t153;
t111 = -t126 * t139 - t142 * t175;
t171 = t148 * t152;
t174 = t147 * t153;
t128 = -t143 * t171 - t174;
t173 = t148 * t140;
t112 = -t128 * t139 + t142 * t173;
t176 = t140 * t152;
t125 = -t139 * t176 + t142 * t143;
t103 = -g(1) * t112 - g(2) * t111 - g(3) * t125;
t146 = sin(qJ(3));
t183 = pkin(3) * t146;
t179 = pkin(10) + qJ(4);
t178 = t143 * pkin(9) + pkin(8);
t177 = t140 * t147;
t169 = t153 * pkin(1) + pkin(9) * t173;
t136 = t148 * pkin(1);
t168 = -pkin(9) * t175 + t136;
t167 = g(1) * t148 - g(2) * t153;
t123 = t139 * t183 + t179 * t142;
t124 = -t179 * t139 + t142 * t183;
t151 = cos(qJ(3));
t134 = pkin(3) * t151 + pkin(2);
t166 = t143 * t123 + t124 * t176 + t134 * t177 + t178;
t129 = -t143 * t172 + t170;
t165 = t123 * t173 + t128 * t124 + t129 * t134 + t169;
t138 = sin(pkin(13));
t141 = cos(pkin(13));
t164 = t138 * t151 + t141 * t146;
t131 = -t138 * t146 + t141 * t151;
t120 = t164 * t139;
t122 = t164 * t142;
t127 = t143 * t174 + t171;
t100 = -t120 * t175 + t126 * t122 + t127 * t131;
t145 = sin(qJ(5));
t150 = cos(qJ(5));
t91 = t100 * t145 - t111 * t150;
t102 = t120 * t173 + t122 * t128 + t129 * t131;
t93 = t102 * t145 - t112 * t150;
t106 = t120 * t143 + (t122 * t152 + t131 * t147) * t140;
t97 = t106 * t145 - t125 * t150;
t163 = g(1) * t93 + g(2) * t91 + g(3) * t97;
t162 = t131 * t139;
t121 = t131 * t142;
t160 = t140 * t162;
t101 = t128 * t121 - t129 * t164 + t148 * t160;
t105 = (t152 * t121 - t147 * t164) * t140 + t143 * t162;
t99 = t121 * t126 - t127 * t164 - t153 * t160;
t161 = g(1) * t101 + g(2) * t99 + g(3) * t105;
t159 = t126 * t124 + t127 * t134 + t136 + (-pkin(9) - t123) * t175;
t158 = g(1) * t129 + g(2) * t127 + g(3) * t177;
t157 = t106 * pkin(4) - pkin(11) * t105 + t166;
t156 = t102 * pkin(4) - pkin(11) * t101 + t165;
t155 = t100 * pkin(4) - t99 * pkin(11) + t159;
t154 = -g(1) * (t128 * t142 + t139 * t173) - g(2) * (t126 * t142 - t139 * t175) - g(3) * (t139 * t143 + t142 * t176);
t149 = cos(qJ(6));
t144 = sin(qJ(6));
t98 = t106 * t150 + t125 * t145;
t94 = t102 * t150 + t112 * t145;
t92 = t100 * t150 + t111 * t145;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t153 - g(2) * t148, t167, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -t158, -g(1) * t128 - g(2) * t126 - g(3) * t176, -g(3) * t143 - t167 * t140, -g(1) * t169 - g(2) * t168 - g(3) * t178, 0, 0, 0, 0, 0, 0, t154 * t146 - t158 * t151, t158 * t146 + t154 * t151, t103, -g(1) * (pkin(2) * t129 + t169) - g(2) * (t127 * pkin(2) + t168) - g(3) * (pkin(2) * t177 + t178) + t103 * pkin(10), 0, 0, 0, 0, 0, 0, -g(1) * t102 - g(2) * t100 - g(3) * t106, -t161, t103, -g(1) * t165 - g(2) * t159 - g(3) * t166, 0, 0, 0, 0, 0, 0, -g(1) * t94 - g(2) * t92 - g(3) * t98, t163, t161, -g(1) * t156 - g(2) * t155 - g(3) * t157, 0, 0, 0, 0, 0, 0, -g(1) * (-t101 * t144 + t149 * t94) - g(2) * (-t144 * t99 + t149 * t92) - g(3) * (-t105 * t144 + t149 * t98) -g(1) * (-t101 * t149 - t144 * t94) - g(2) * (-t144 * t92 - t149 * t99) - g(3) * (-t105 * t149 - t144 * t98) -t163, -g(1) * (pkin(5) * t94 + pkin(12) * t93 + t156) - g(2) * (t92 * pkin(5) + t91 * pkin(12) + t155) - g(3) * (pkin(5) * t98 + pkin(12) * t97 + t157);];
U_reg  = t1;
