% Calculate inertial parameters regressor of potential energy for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:34:08
% EndTime: 2019-03-09 07:34:08
% DurationCPUTime: 0.38s
% Computational Cost: add. (502->110), mult. (1179->168), div. (0->0), fcn. (1486->16), ass. (0->65)
t146 = cos(pkin(6));
t141 = sin(pkin(13));
t153 = cos(qJ(1));
t172 = t153 * t141;
t144 = cos(pkin(13));
t150 = sin(qJ(1));
t173 = t150 * t144;
t125 = -t146 * t173 - t172;
t142 = sin(pkin(7));
t145 = cos(pkin(7));
t143 = sin(pkin(6));
t177 = t143 * t150;
t115 = -t125 * t142 + t145 * t177;
t178 = t143 * t144;
t122 = -t142 * t178 + t146 * t145;
t185 = cos(qJ(3));
t184 = t146 * qJ(2) + pkin(8);
t171 = t153 * t144;
t174 = t150 * t141;
t123 = t146 * t171 - t174;
t176 = t143 * t153;
t114 = -t123 * t142 - t145 * t176;
t148 = sin(qJ(4));
t183 = t114 * t148;
t182 = t115 * t148;
t181 = t122 * t148;
t179 = t141 * t143;
t170 = t153 * pkin(1) + qJ(2) * t177;
t167 = t142 * t185;
t166 = t145 * t185;
t165 = t143 * t167;
t164 = g(1) * t150 - g(2) * t153;
t163 = t150 * pkin(1) - qJ(2) * t176;
t124 = t146 * t172 + t173;
t149 = sin(qJ(3));
t103 = t124 * t185 + (t123 * t145 - t142 * t176) * t149;
t140 = qJ(4) + qJ(5);
t136 = sin(t140);
t137 = cos(t140);
t92 = t103 * t136 - t114 * t137;
t126 = -t146 * t174 + t171;
t105 = t126 * t185 + (t125 * t145 + t142 * t177) * t149;
t94 = t105 * t136 - t115 * t137;
t113 = t146 * t142 * t149 + (t144 * t145 * t149 + t185 * t141) * t143;
t98 = t113 * t136 - t122 * t137;
t162 = g(1) * t94 + g(2) * t92 + g(3) * t98;
t102 = -t123 * t166 + t124 * t149 + t153 * t165;
t104 = -t125 * t166 + t126 * t149 - t150 * t165;
t112 = -t146 * t167 + t149 * t179 - t166 * t178;
t161 = g(1) * t104 + g(2) * t102 + g(3) * t112;
t160 = t126 * pkin(2) + t115 * pkin(9) + t170;
t159 = pkin(2) * t179 + t122 * pkin(9) + t184;
t152 = cos(qJ(4));
t134 = t152 * pkin(4) + pkin(3);
t154 = -pkin(11) - pkin(10);
t158 = pkin(4) * t182 - t104 * t154 + t105 * t134 + t160;
t157 = pkin(4) * t181 - t112 * t154 + t113 * t134 + t159;
t156 = t124 * pkin(2) + t114 * pkin(9) + t163;
t155 = pkin(4) * t183 - t102 * t154 + t103 * t134 + t156;
t151 = cos(qJ(6));
t147 = sin(qJ(6));
t99 = t113 * t137 + t122 * t136;
t95 = t105 * t137 + t115 * t136;
t93 = t103 * t137 + t114 * t136;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t153 - g(2) * t150, t164, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t126 - g(2) * t124 - g(3) * t179, -g(1) * t125 - g(2) * t123 - g(3) * t178, -g(3) * t146 - t164 * t143, -g(1) * t170 - g(2) * t163 - g(3) * t184, 0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t103 - g(3) * t113, t161, -g(1) * t115 - g(2) * t114 - g(3) * t122, -g(1) * t160 - g(2) * t156 - g(3) * t159, 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t152 + t182) - g(2) * (t103 * t152 + t183) - g(3) * (t113 * t152 + t181) -g(1) * (-t105 * t148 + t115 * t152) - g(2) * (-t103 * t148 + t114 * t152) - g(3) * (-t113 * t148 + t122 * t152) -t161, -g(1) * (t105 * pkin(3) + t104 * pkin(10) + t160) - g(2) * (t103 * pkin(3) + t102 * pkin(10) + t156) - g(3) * (t113 * pkin(3) + t112 * pkin(10) + t159) 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t99, t162, -t161, -g(1) * t158 - g(2) * t155 - g(3) * t157, 0, 0, 0, 0, 0, 0, -g(1) * (t104 * t147 + t95 * t151) - g(2) * (t102 * t147 + t93 * t151) - g(3) * (t112 * t147 + t99 * t151) -g(1) * (t104 * t151 - t95 * t147) - g(2) * (t102 * t151 - t93 * t147) - g(3) * (t112 * t151 - t99 * t147) -t162, -g(1) * (t95 * pkin(5) + t94 * pkin(12) + t158) - g(2) * (t93 * pkin(5) + t92 * pkin(12) + t155) - g(3) * (t99 * pkin(5) + t98 * pkin(12) + t157);];
U_reg  = t1;
