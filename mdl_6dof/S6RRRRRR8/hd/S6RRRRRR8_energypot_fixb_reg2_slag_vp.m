% Calculate inertial parameters regressor of potential energy for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:06:43
% EndTime: 2019-03-10 05:06:43
% DurationCPUTime: 0.37s
% Computational Cost: add. (502->110), mult. (1179->168), div. (0->0), fcn. (1486->16), ass. (0->65)
t144 = cos(pkin(6));
t148 = sin(qJ(2));
t153 = cos(qJ(1));
t172 = t153 * t148;
t149 = sin(qJ(1));
t152 = cos(qJ(2));
t173 = t149 * t152;
t125 = -t144 * t173 - t172;
t141 = sin(pkin(7));
t143 = cos(pkin(7));
t142 = sin(pkin(6));
t178 = t142 * t149;
t115 = -t125 * t141 + t143 * t178;
t177 = t142 * t152;
t122 = -t141 * t177 + t144 * t143;
t185 = cos(qJ(3));
t184 = t144 * pkin(9) + pkin(8);
t171 = t153 * t152;
t174 = t149 * t148;
t123 = t144 * t171 - t174;
t176 = t142 * t153;
t114 = -t123 * t141 - t143 * t176;
t146 = sin(qJ(4));
t183 = t114 * t146;
t182 = t115 * t146;
t181 = t122 * t146;
t179 = t142 * t148;
t170 = t153 * pkin(1) + pkin(9) * t178;
t167 = t141 * t185;
t166 = t143 * t185;
t165 = t142 * t167;
t164 = t149 * pkin(1) - pkin(9) * t176;
t163 = g(1) * t149 - g(2) * t153;
t124 = t144 * t172 + t173;
t147 = sin(qJ(3));
t103 = t124 * t185 + (t123 * t143 - t141 * t176) * t147;
t140 = qJ(4) + qJ(5);
t135 = sin(t140);
t136 = cos(t140);
t92 = t103 * t135 - t114 * t136;
t126 = -t144 * t174 + t171;
t105 = t126 * t185 + (t125 * t143 + t141 * t178) * t147;
t94 = t105 * t135 - t115 * t136;
t113 = t144 * t141 * t147 + (t143 * t147 * t152 + t185 * t148) * t142;
t98 = t113 * t135 - t122 * t136;
t162 = g(1) * t94 + g(2) * t92 + g(3) * t98;
t102 = -t123 * t166 + t124 * t147 + t153 * t165;
t104 = -t125 * t166 + t126 * t147 - t149 * t165;
t112 = -t144 * t167 + t147 * t179 - t166 * t177;
t161 = g(1) * t104 + g(2) * t102 + g(3) * t112;
t160 = t126 * pkin(2) + t115 * pkin(10) + t170;
t159 = pkin(2) * t179 + t122 * pkin(10) + t184;
t151 = cos(qJ(4));
t134 = t151 * pkin(4) + pkin(3);
t154 = -pkin(12) - pkin(11);
t158 = pkin(4) * t182 - t104 * t154 + t105 * t134 + t160;
t157 = pkin(4) * t181 - t112 * t154 + t113 * t134 + t159;
t156 = t124 * pkin(2) + t114 * pkin(10) + t164;
t155 = pkin(4) * t183 - t102 * t154 + t103 * t134 + t156;
t150 = cos(qJ(6));
t145 = sin(qJ(6));
t99 = t113 * t136 + t122 * t135;
t95 = t105 * t136 + t115 * t135;
t93 = t103 * t136 + t114 * t135;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t153 - g(2) * t149, t163, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t126 - g(2) * t124 - g(3) * t179, -g(1) * t125 - g(2) * t123 - g(3) * t177, -g(3) * t144 - t163 * t142, -g(1) * t170 - g(2) * t164 - g(3) * t184, 0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t103 - g(3) * t113, t161, -g(1) * t115 - g(2) * t114 - g(3) * t122, -g(1) * t160 - g(2) * t156 - g(3) * t159, 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t151 + t182) - g(2) * (t103 * t151 + t183) - g(3) * (t113 * t151 + t181) -g(1) * (-t105 * t146 + t115 * t151) - g(2) * (-t103 * t146 + t114 * t151) - g(3) * (-t113 * t146 + t122 * t151) -t161, -g(1) * (t105 * pkin(3) + t104 * pkin(11) + t160) - g(2) * (t103 * pkin(3) + t102 * pkin(11) + t156) - g(3) * (t113 * pkin(3) + t112 * pkin(11) + t159) 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t99, t162, -t161, -g(1) * t158 - g(2) * t155 - g(3) * t157, 0, 0, 0, 0, 0, 0, -g(1) * (t104 * t145 + t95 * t150) - g(2) * (t102 * t145 + t93 * t150) - g(3) * (t112 * t145 + t99 * t150) -g(1) * (t104 * t150 - t95 * t145) - g(2) * (t102 * t150 - t93 * t145) - g(3) * (t112 * t150 - t99 * t145) -t162, -g(1) * (t95 * pkin(5) + t94 * pkin(13) + t158) - g(2) * (t93 * pkin(5) + t92 * pkin(13) + t155) - g(3) * (t99 * pkin(5) + t98 * pkin(13) + t157);];
U_reg  = t1;
