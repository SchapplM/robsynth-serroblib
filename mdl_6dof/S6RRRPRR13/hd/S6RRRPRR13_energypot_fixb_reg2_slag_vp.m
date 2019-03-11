% Calculate inertial parameters regressor of potential energy for
% S6RRRPRR13
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
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR13_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:05:07
% EndTime: 2019-03-09 20:05:08
% DurationCPUTime: 0.37s
% Computational Cost: add. (502->110), mult. (1179->168), div. (0->0), fcn. (1486->16), ass. (0->65)
t147 = cos(pkin(6));
t151 = sin(qJ(2));
t155 = cos(qJ(1));
t173 = t155 * t151;
t152 = sin(qJ(1));
t154 = cos(qJ(2));
t174 = t152 * t154;
t126 = -t147 * t174 - t173;
t143 = sin(pkin(7));
t146 = cos(pkin(7));
t144 = sin(pkin(6));
t179 = t144 * t152;
t116 = -t126 * t143 + t146 * t179;
t178 = t144 * t154;
t123 = -t143 * t178 + t147 * t146;
t186 = cos(qJ(3));
t185 = t147 * pkin(9) + pkin(8);
t172 = t155 * t154;
t175 = t152 * t151;
t124 = t147 * t172 - t175;
t177 = t144 * t155;
t115 = -t124 * t143 - t146 * t177;
t142 = sin(pkin(13));
t184 = t115 * t142;
t183 = t116 * t142;
t182 = t123 * t142;
t180 = t144 * t151;
t171 = t155 * pkin(1) + pkin(9) * t179;
t168 = t143 * t186;
t167 = t146 * t186;
t166 = t144 * t168;
t165 = t152 * pkin(1) - pkin(9) * t177;
t164 = g(1) * t152 - g(2) * t155;
t125 = t147 * t173 + t174;
t150 = sin(qJ(3));
t104 = t125 * t186 + (t124 * t146 - t143 * t177) * t150;
t141 = pkin(13) + qJ(5);
t136 = sin(t141);
t137 = cos(t141);
t93 = t104 * t136 - t115 * t137;
t127 = -t147 * t175 + t172;
t106 = t127 * t186 + (t126 * t146 + t143 * t179) * t150;
t95 = t106 * t136 - t116 * t137;
t114 = t147 * t143 * t150 + (t146 * t150 * t154 + t186 * t151) * t144;
t99 = t114 * t136 - t123 * t137;
t163 = g(1) * t95 + g(2) * t93 + g(3) * t99;
t103 = -t124 * t167 + t125 * t150 + t155 * t166;
t105 = -t126 * t167 + t127 * t150 - t152 * t166;
t113 = -t147 * t168 + t150 * t180 - t167 * t178;
t162 = g(1) * t105 + g(2) * t103 + g(3) * t113;
t161 = t127 * pkin(2) + t116 * pkin(10) + t171;
t160 = pkin(2) * t180 + t123 * pkin(10) + t185;
t145 = cos(pkin(13));
t135 = t145 * pkin(4) + pkin(3);
t148 = -pkin(11) - qJ(4);
t159 = pkin(4) * t183 - t105 * t148 + t106 * t135 + t161;
t158 = pkin(4) * t182 - t113 * t148 + t114 * t135 + t160;
t157 = t125 * pkin(2) + t115 * pkin(10) + t165;
t156 = pkin(4) * t184 - t103 * t148 + t104 * t135 + t157;
t153 = cos(qJ(6));
t149 = sin(qJ(6));
t100 = t114 * t137 + t123 * t136;
t96 = t106 * t137 + t116 * t136;
t94 = t104 * t137 + t115 * t136;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t155 - g(2) * t152, t164, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t127 - g(2) * t125 - g(3) * t180, -g(1) * t126 - g(2) * t124 - g(3) * t178, -g(3) * t147 - t164 * t144, -g(1) * t171 - g(2) * t165 - g(3) * t185, 0, 0, 0, 0, 0, 0, -g(1) * t106 - g(2) * t104 - g(3) * t114, t162, -g(1) * t116 - g(2) * t115 - g(3) * t123, -g(1) * t161 - g(2) * t157 - g(3) * t160, 0, 0, 0, 0, 0, 0, -g(1) * (t106 * t145 + t183) - g(2) * (t104 * t145 + t184) - g(3) * (t114 * t145 + t182) -g(1) * (-t106 * t142 + t116 * t145) - g(2) * (-t104 * t142 + t115 * t145) - g(3) * (-t114 * t142 + t123 * t145) -t162, -g(1) * (t106 * pkin(3) + t105 * qJ(4) + t161) - g(2) * (t104 * pkin(3) + t103 * qJ(4) + t157) - g(3) * (t114 * pkin(3) + t113 * qJ(4) + t160) 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t100, t163, -t162, -g(1) * t159 - g(2) * t156 - g(3) * t158, 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t149 + t96 * t153) - g(2) * (t103 * t149 + t94 * t153) - g(3) * (t100 * t153 + t113 * t149) -g(1) * (t105 * t153 - t96 * t149) - g(2) * (t103 * t153 - t94 * t149) - g(3) * (-t100 * t149 + t113 * t153) -t163, -g(1) * (t96 * pkin(5) + t95 * pkin(12) + t159) - g(2) * (t94 * pkin(5) + t93 * pkin(12) + t156) - g(3) * (t100 * pkin(5) + t99 * pkin(12) + t158);];
U_reg  = t1;
