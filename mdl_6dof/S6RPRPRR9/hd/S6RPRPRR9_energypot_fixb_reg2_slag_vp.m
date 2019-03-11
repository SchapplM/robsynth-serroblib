% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:05:04
% EndTime: 2019-03-09 04:05:05
% DurationCPUTime: 0.35s
% Computational Cost: add. (574->109), mult. (1452->170), div. (0->0), fcn. (1852->16), ass. (0->71)
t145 = cos(pkin(12));
t141 = sin(pkin(12));
t151 = sin(qJ(1));
t173 = t151 * t141;
t147 = cos(pkin(6));
t155 = cos(qJ(1));
t174 = t147 * t155;
t128 = t145 * t174 - t173;
t142 = sin(pkin(7));
t146 = cos(pkin(7));
t143 = sin(pkin(6));
t175 = t143 * t155;
t113 = -t128 * t142 - t146 * t175;
t172 = t151 * t145;
t130 = -t141 * t155 - t147 * t172;
t176 = t143 * t151;
t114 = -t130 * t142 + t146 * t176;
t177 = t143 * t145;
t127 = -t142 * t177 + t146 * t147;
t105 = -g(1) * t114 - g(2) * t113 - g(3) * t127;
t150 = sin(qJ(3));
t184 = pkin(3) * t150;
t180 = pkin(9) + qJ(4);
t179 = t147 * qJ(2) + pkin(8);
t178 = t141 * t143;
t171 = t155 * pkin(1) + qJ(2) * t176;
t170 = g(1) * t151 - g(2) * t155;
t138 = t151 * pkin(1);
t169 = -qJ(2) * t175 + t138;
t125 = t142 * t184 + t180 * t146;
t126 = -t180 * t142 + t146 * t184;
t154 = cos(qJ(3));
t136 = pkin(3) * t154 + pkin(2);
t168 = t147 * t125 + t126 * t177 + t136 * t178 + t179;
t131 = t145 * t155 - t147 * t173;
t167 = t125 * t176 + t130 * t126 + t131 * t136 + t171;
t140 = sin(pkin(13));
t144 = cos(pkin(13));
t166 = t140 * t154 + t144 * t150;
t133 = -t140 * t150 + t144 * t154;
t122 = t166 * t142;
t124 = t166 * t146;
t129 = t141 * t174 + t172;
t102 = -t122 * t175 + t128 * t124 + t129 * t133;
t149 = sin(qJ(5));
t153 = cos(qJ(5));
t93 = t102 * t149 - t113 * t153;
t104 = t122 * t176 + t124 * t130 + t131 * t133;
t95 = t104 * t149 - t114 * t153;
t108 = t122 * t147 + (t124 * t145 + t133 * t141) * t143;
t99 = t108 * t149 - t127 * t153;
t165 = g(1) * t95 + g(2) * t93 + g(3) * t99;
t164 = t133 * t142;
t123 = t133 * t146;
t162 = t143 * t164;
t101 = t123 * t128 - t129 * t166 - t155 * t162;
t103 = t130 * t123 - t131 * t166 + t151 * t162;
t107 = (t123 * t145 - t141 * t166) * t143 + t147 * t164;
t163 = g(1) * t103 + g(2) * t101 + g(3) * t107;
t161 = g(1) * t131 + g(2) * t129 + g(3) * t178;
t160 = t108 * pkin(4) - pkin(10) * t107 + t168;
t159 = t128 * t126 + t129 * t136 + t138 + (-qJ(2) - t125) * t175;
t158 = t104 * pkin(4) - pkin(10) * t103 + t167;
t157 = t102 * pkin(4) - t101 * pkin(10) + t159;
t156 = -g(1) * (t130 * t146 + t142 * t176) - g(2) * (t128 * t146 - t142 * t175) - g(3) * (t142 * t147 + t146 * t177);
t152 = cos(qJ(6));
t148 = sin(qJ(6));
t100 = t108 * t153 + t127 * t149;
t96 = t104 * t153 + t114 * t149;
t94 = t102 * t153 + t113 * t149;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t155 - g(2) * t151, t170, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -t161, -g(1) * t130 - g(2) * t128 - g(3) * t177, -g(3) * t147 - t170 * t143, -g(1) * t171 - g(2) * t169 - g(3) * t179, 0, 0, 0, 0, 0, 0, t156 * t150 - t161 * t154, t161 * t150 + t156 * t154, t105, -g(1) * (pkin(2) * t131 + t171) - g(2) * (t129 * pkin(2) + t169) - g(3) * (pkin(2) * t178 + t179) + t105 * pkin(9), 0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t102 - g(3) * t108, -t163, t105, -g(1) * t167 - g(2) * t159 - g(3) * t168, 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t100, t165, t163, -g(1) * t158 - g(2) * t157 - g(3) * t160, 0, 0, 0, 0, 0, 0, -g(1) * (-t103 * t148 + t152 * t96) - g(2) * (-t101 * t148 + t152 * t94) - g(3) * (t100 * t152 - t107 * t148) -g(1) * (-t103 * t152 - t148 * t96) - g(2) * (-t101 * t152 - t148 * t94) - g(3) * (-t100 * t148 - t107 * t152) -t165, -g(1) * (pkin(5) * t96 + pkin(11) * t95 + t158) - g(2) * (t94 * pkin(5) + t93 * pkin(11) + t157) - g(3) * (pkin(5) * t100 + pkin(11) * t99 + t160);];
U_reg  = t1;
