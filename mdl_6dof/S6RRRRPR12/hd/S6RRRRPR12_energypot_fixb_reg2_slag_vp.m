% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:44:39
% EndTime: 2019-03-09 23:44:39
% DurationCPUTime: 0.37s
% Computational Cost: add. (502->110), mult. (1179->168), div. (0->0), fcn. (1486->16), ass. (0->65)
t145 = cos(pkin(6));
t150 = sin(qJ(2));
t155 = cos(qJ(1));
t173 = t155 * t150;
t151 = sin(qJ(1));
t154 = cos(qJ(2));
t174 = t151 * t154;
t126 = -t145 * t174 - t173;
t142 = sin(pkin(7));
t144 = cos(pkin(7));
t143 = sin(pkin(6));
t179 = t143 * t151;
t116 = -t126 * t142 + t144 * t179;
t178 = t143 * t154;
t123 = -t142 * t178 + t145 * t144;
t186 = cos(qJ(3));
t185 = t145 * pkin(9) + pkin(8);
t172 = t155 * t154;
t175 = t151 * t150;
t124 = t145 * t172 - t175;
t177 = t143 * t155;
t115 = -t124 * t142 - t144 * t177;
t148 = sin(qJ(4));
t184 = t115 * t148;
t183 = t116 * t148;
t182 = t123 * t148;
t180 = t143 * t150;
t171 = t155 * pkin(1) + pkin(9) * t179;
t168 = t142 * t186;
t167 = t144 * t186;
t166 = t143 * t168;
t165 = t151 * pkin(1) - pkin(9) * t177;
t164 = g(1) * t151 - g(2) * t155;
t125 = t145 * t173 + t174;
t149 = sin(qJ(3));
t104 = t125 * t186 + (t124 * t144 - t142 * t177) * t149;
t141 = qJ(4) + pkin(13);
t136 = sin(t141);
t137 = cos(t141);
t93 = t104 * t136 - t115 * t137;
t127 = -t145 * t175 + t172;
t106 = t127 * t186 + (t126 * t144 + t142 * t179) * t149;
t95 = t106 * t136 - t116 * t137;
t114 = t145 * t142 * t149 + (t144 * t149 * t154 + t186 * t150) * t143;
t99 = t114 * t136 - t123 * t137;
t163 = g(1) * t95 + g(2) * t93 + g(3) * t99;
t103 = -t124 * t167 + t125 * t149 + t155 * t166;
t105 = -t126 * t167 + t127 * t149 - t151 * t166;
t113 = -t145 * t168 + t149 * t180 - t167 * t178;
t162 = g(1) * t105 + g(2) * t103 + g(3) * t113;
t161 = t127 * pkin(2) + t116 * pkin(10) + t171;
t160 = pkin(2) * t180 + t123 * pkin(10) + t185;
t153 = cos(qJ(4));
t135 = t153 * pkin(4) + pkin(3);
t146 = -qJ(5) - pkin(11);
t159 = pkin(4) * t183 - t105 * t146 + t106 * t135 + t161;
t158 = pkin(4) * t182 - t113 * t146 + t114 * t135 + t160;
t157 = t125 * pkin(2) + t115 * pkin(10) + t165;
t156 = pkin(4) * t184 - t103 * t146 + t104 * t135 + t157;
t152 = cos(qJ(6));
t147 = sin(qJ(6));
t100 = t114 * t137 + t123 * t136;
t96 = t106 * t137 + t116 * t136;
t94 = t104 * t137 + t115 * t136;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t155 - g(2) * t151, t164, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t127 - g(2) * t125 - g(3) * t180, -g(1) * t126 - g(2) * t124 - g(3) * t178, -g(3) * t145 - t164 * t143, -g(1) * t171 - g(2) * t165 - g(3) * t185, 0, 0, 0, 0, 0, 0, -g(1) * t106 - g(2) * t104 - g(3) * t114, t162, -g(1) * t116 - g(2) * t115 - g(3) * t123, -g(1) * t161 - g(2) * t157 - g(3) * t160, 0, 0, 0, 0, 0, 0, -g(1) * (t106 * t153 + t183) - g(2) * (t104 * t153 + t184) - g(3) * (t114 * t153 + t182) -g(1) * (-t106 * t148 + t116 * t153) - g(2) * (-t104 * t148 + t115 * t153) - g(3) * (-t114 * t148 + t123 * t153) -t162, -g(1) * (t106 * pkin(3) + t105 * pkin(11) + t161) - g(2) * (t104 * pkin(3) + t103 * pkin(11) + t157) - g(3) * (t114 * pkin(3) + t113 * pkin(11) + t160) 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t100, t163, -t162, -g(1) * t159 - g(2) * t156 - g(3) * t158, 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t147 + t96 * t152) - g(2) * (t103 * t147 + t94 * t152) - g(3) * (t100 * t152 + t113 * t147) -g(1) * (t105 * t152 - t96 * t147) - g(2) * (t103 * t152 - t94 * t147) - g(3) * (-t100 * t147 + t113 * t152) -t163, -g(1) * (t96 * pkin(5) + t95 * pkin(12) + t159) - g(2) * (t94 * pkin(5) + t93 * pkin(12) + t156) - g(3) * (t100 * pkin(5) + t99 * pkin(12) + t158);];
U_reg  = t1;
