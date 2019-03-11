% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR13_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:25:19
% EndTime: 2019-03-09 04:25:20
% DurationCPUTime: 0.31s
% Computational Cost: add. (453->95), mult. (1177->145), div. (0->0), fcn. (1483->14), ass. (0->64)
t144 = cos(pkin(6));
t139 = sin(pkin(12));
t151 = cos(qJ(1));
t173 = t151 * t139;
t142 = cos(pkin(12));
t148 = sin(qJ(1));
t174 = t148 * t142;
t126 = -t144 * t174 - t173;
t140 = sin(pkin(7));
t143 = cos(pkin(7));
t141 = sin(pkin(6));
t179 = t141 * t148;
t117 = -t126 * t140 + t143 * t179;
t180 = t141 * t142;
t123 = -t140 * t180 + t144 * t143;
t185 = cos(qJ(3));
t184 = t144 * qJ(2) + pkin(8);
t182 = t139 * t141;
t147 = sin(qJ(3));
t181 = t140 * t147;
t178 = t141 * t151;
t177 = t143 * t147;
t175 = t148 * t139;
t172 = t151 * t142;
t171 = t151 * pkin(1) + qJ(2) * t179;
t168 = t140 * t185;
t167 = t143 * t185;
t166 = t141 * t168;
t165 = g(1) * t148 - g(2) * t151;
t164 = t148 * pkin(1) - qJ(2) * t178;
t124 = t144 * t172 - t175;
t116 = -t124 * t140 - t143 * t178;
t112 = -t144 * t168 + t147 * t182 - t167 * t180;
t146 = sin(qJ(5));
t150 = cos(qJ(5));
t104 = -t112 * t150 + t123 * t146;
t125 = t144 * t173 + t174;
t106 = -t124 * t167 + t125 * t147 + t151 * t166;
t95 = -t106 * t150 + t116 * t146;
t127 = -t144 * t175 + t172;
t108 = -t126 * t167 + t127 * t147 - t148 * t166;
t97 = -t108 * t150 + t117 * t146;
t163 = g(1) * t97 + g(2) * t95 + g(3) * t104;
t162 = g(1) * t108 + g(2) * t106 + g(3) * t112;
t107 = t124 * t177 + t125 * t185 - t178 * t181;
t109 = t127 * t185 + (t126 * t143 + t140 * t179) * t147;
t113 = t144 * t181 + (t185 * t139 + t142 * t177) * t141;
t161 = g(1) * t109 + g(2) * t107 + g(3) * t113;
t160 = t127 * pkin(2) + t117 * pkin(9) + t171;
t159 = pkin(2) * t182 + t123 * pkin(9) + t184;
t158 = t109 * pkin(3) + t108 * qJ(4) + t160;
t157 = t113 * pkin(3) + t112 * qJ(4) + t159;
t156 = t125 * pkin(2) + t116 * pkin(9) + t164;
t155 = t117 * pkin(4) + t109 * pkin(10) + t158;
t154 = t123 * pkin(4) + t113 * pkin(10) + t157;
t153 = t107 * pkin(3) + t106 * qJ(4) + t156;
t152 = t116 * pkin(4) + t107 * pkin(10) + t153;
t149 = cos(qJ(6));
t145 = sin(qJ(6));
t105 = t112 * t146 + t123 * t150;
t99 = -g(1) * t117 - g(2) * t116 - g(3) * t123;
t98 = t108 * t146 + t117 * t150;
t96 = t106 * t146 + t116 * t150;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t151 - g(2) * t148, t165, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t127 - g(2) * t125 - g(3) * t182, -g(1) * t126 - g(2) * t124 - g(3) * t180, -g(3) * t144 - t165 * t141, -g(1) * t171 - g(2) * t164 - g(3) * t184, 0, 0, 0, 0, 0, 0, -t161, t162, t99, -g(1) * t160 - g(2) * t156 - g(3) * t159, 0, 0, 0, 0, 0, 0, t99, t161, -t162, -g(1) * t158 - g(2) * t153 - g(3) * t157, 0, 0, 0, 0, 0, 0, -g(1) * t98 - g(2) * t96 - g(3) * t105, t163, -t161, -g(1) * t155 - g(2) * t152 - g(3) * t154, 0, 0, 0, 0, 0, 0, -g(1) * (t109 * t145 + t98 * t149) - g(2) * (t107 * t145 + t96 * t149) - g(3) * (t105 * t149 + t113 * t145) -g(1) * (t109 * t149 - t98 * t145) - g(2) * (t107 * t149 - t96 * t145) - g(3) * (-t105 * t145 + t113 * t149) -t163, -g(1) * (t98 * pkin(5) + t97 * pkin(11) + t155) - g(2) * (t96 * pkin(5) + t95 * pkin(11) + t152) - g(3) * (t105 * pkin(5) + t104 * pkin(11) + t154);];
U_reg  = t1;
