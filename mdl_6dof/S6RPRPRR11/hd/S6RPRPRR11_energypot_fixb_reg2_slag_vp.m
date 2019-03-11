% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR11
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
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:15:11
% EndTime: 2019-03-09 04:15:11
% DurationCPUTime: 0.36s
% Computational Cost: add. (502->110), mult. (1179->168), div. (0->0), fcn. (1486->16), ass. (0->65)
t152 = cos(pkin(6));
t146 = sin(pkin(12));
t158 = cos(qJ(1));
t176 = t158 * t146;
t150 = cos(pkin(12));
t156 = sin(qJ(1));
t177 = t156 * t150;
t129 = -t152 * t177 - t176;
t147 = sin(pkin(7));
t151 = cos(pkin(7));
t148 = sin(pkin(6));
t181 = t148 * t156;
t119 = -t129 * t147 + t151 * t181;
t182 = t148 * t150;
t126 = -t147 * t182 + t152 * t151;
t189 = cos(qJ(3));
t188 = t152 * qJ(2) + pkin(8);
t175 = t158 * t150;
t178 = t156 * t146;
t127 = t152 * t175 - t178;
t180 = t148 * t158;
t118 = -t127 * t147 - t151 * t180;
t145 = sin(pkin(13));
t187 = t118 * t145;
t186 = t119 * t145;
t185 = t126 * t145;
t183 = t146 * t148;
t174 = t158 * pkin(1) + qJ(2) * t181;
t171 = t147 * t189;
t170 = t151 * t189;
t169 = t148 * t171;
t168 = g(1) * t156 - g(2) * t158;
t167 = t156 * pkin(1) - qJ(2) * t180;
t155 = sin(qJ(3));
t117 = t152 * t147 * t155 + (t150 * t151 * t155 + t189 * t146) * t148;
t144 = pkin(13) + qJ(5);
t139 = sin(t144);
t140 = cos(t144);
t102 = t117 * t139 - t126 * t140;
t128 = t152 * t176 + t177;
t107 = t128 * t189 + (t127 * t151 - t147 * t180) * t155;
t96 = t107 * t139 - t118 * t140;
t130 = -t152 * t178 + t175;
t109 = t130 * t189 + (t129 * t151 + t147 * t181) * t155;
t98 = t109 * t139 - t119 * t140;
t166 = g(1) * t98 + g(2) * t96 + g(3) * t102;
t106 = -t127 * t170 + t128 * t155 + t158 * t169;
t108 = -t129 * t170 + t130 * t155 - t156 * t169;
t116 = -t152 * t171 + t155 * t183 - t170 * t182;
t165 = g(1) * t108 + g(2) * t106 + g(3) * t116;
t164 = t130 * pkin(2) + t119 * pkin(9) + t174;
t163 = pkin(2) * t183 + t126 * pkin(9) + t188;
t149 = cos(pkin(13));
t138 = t149 * pkin(4) + pkin(3);
t153 = -pkin(10) - qJ(4);
t162 = pkin(4) * t186 - t108 * t153 + t109 * t138 + t164;
t161 = pkin(4) * t185 - t116 * t153 + t117 * t138 + t163;
t160 = t128 * pkin(2) + t118 * pkin(9) + t167;
t159 = pkin(4) * t187 - t106 * t153 + t107 * t138 + t160;
t157 = cos(qJ(6));
t154 = sin(qJ(6));
t103 = t117 * t140 + t126 * t139;
t99 = t109 * t140 + t119 * t139;
t97 = t107 * t140 + t118 * t139;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t158 - g(2) * t156, t168, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t130 - g(2) * t128 - g(3) * t183, -g(1) * t129 - g(2) * t127 - g(3) * t182, -g(3) * t152 - t168 * t148, -g(1) * t174 - g(2) * t167 - g(3) * t188, 0, 0, 0, 0, 0, 0, -g(1) * t109 - g(2) * t107 - g(3) * t117, t165, -g(1) * t119 - g(2) * t118 - g(3) * t126, -g(1) * t164 - g(2) * t160 - g(3) * t163, 0, 0, 0, 0, 0, 0, -g(1) * (t109 * t149 + t186) - g(2) * (t107 * t149 + t187) - g(3) * (t117 * t149 + t185) -g(1) * (-t109 * t145 + t119 * t149) - g(2) * (-t107 * t145 + t118 * t149) - g(3) * (-t117 * t145 + t126 * t149) -t165, -g(1) * (t109 * pkin(3) + t108 * qJ(4) + t164) - g(2) * (t107 * pkin(3) + t106 * qJ(4) + t160) - g(3) * (t117 * pkin(3) + t116 * qJ(4) + t163) 0, 0, 0, 0, 0, 0, -g(1) * t99 - g(2) * t97 - g(3) * t103, t166, -t165, -g(1) * t162 - g(2) * t159 - g(3) * t161, 0, 0, 0, 0, 0, 0, -g(1) * (t108 * t154 + t99 * t157) - g(2) * (t106 * t154 + t97 * t157) - g(3) * (t103 * t157 + t116 * t154) -g(1) * (t108 * t157 - t99 * t154) - g(2) * (t106 * t157 - t97 * t154) - g(3) * (-t103 * t154 + t116 * t157) -t166, -g(1) * (t99 * pkin(5) + t98 * pkin(11) + t162) - g(2) * (t97 * pkin(5) + t96 * pkin(11) + t159) - g(3) * (t103 * pkin(5) + t102 * pkin(11) + t161);];
U_reg  = t1;
