% Calculate inertial parameters regressor of potential energy for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:31:32
% EndTime: 2019-03-09 05:31:32
% DurationCPUTime: 0.36s
% Computational Cost: add. (502->110), mult. (1179->168), div. (0->0), fcn. (1486->16), ass. (0->65)
t150 = cos(pkin(6));
t145 = sin(pkin(12));
t158 = cos(qJ(1));
t176 = t158 * t145;
t148 = cos(pkin(12));
t155 = sin(qJ(1));
t177 = t155 * t148;
t129 = -t150 * t177 - t176;
t146 = sin(pkin(7));
t149 = cos(pkin(7));
t147 = sin(pkin(6));
t181 = t147 * t155;
t119 = -t129 * t146 + t149 * t181;
t182 = t147 * t148;
t126 = -t146 * t182 + t150 * t149;
t189 = cos(qJ(3));
t188 = t150 * qJ(2) + pkin(8);
t175 = t158 * t148;
t178 = t155 * t145;
t127 = t150 * t175 - t178;
t180 = t147 * t158;
t118 = -t127 * t146 - t149 * t180;
t153 = sin(qJ(4));
t187 = t118 * t153;
t186 = t119 * t153;
t185 = t126 * t153;
t183 = t145 * t147;
t174 = t158 * pkin(1) + qJ(2) * t181;
t171 = t146 * t189;
t170 = t149 * t189;
t169 = t147 * t171;
t168 = g(1) * t155 - g(2) * t158;
t167 = t155 * pkin(1) - qJ(2) * t180;
t154 = sin(qJ(3));
t117 = t150 * t146 * t154 + (t148 * t149 * t154 + t189 * t145) * t147;
t144 = qJ(4) + pkin(13);
t139 = sin(t144);
t140 = cos(t144);
t102 = t117 * t139 - t126 * t140;
t128 = t150 * t176 + t177;
t107 = t128 * t189 + (t127 * t149 - t146 * t180) * t154;
t96 = t107 * t139 - t118 * t140;
t130 = -t150 * t178 + t175;
t109 = t130 * t189 + (t129 * t149 + t146 * t181) * t154;
t98 = t109 * t139 - t119 * t140;
t166 = g(1) * t98 + g(2) * t96 + g(3) * t102;
t106 = -t127 * t170 + t128 * t154 + t158 * t169;
t108 = -t129 * t170 + t130 * t154 - t155 * t169;
t116 = -t150 * t171 + t154 * t183 - t170 * t182;
t165 = g(1) * t108 + g(2) * t106 + g(3) * t116;
t164 = t130 * pkin(2) + t119 * pkin(9) + t174;
t163 = pkin(2) * t183 + t126 * pkin(9) + t188;
t157 = cos(qJ(4));
t138 = t157 * pkin(4) + pkin(3);
t151 = -qJ(5) - pkin(10);
t162 = pkin(4) * t186 - t108 * t151 + t109 * t138 + t164;
t161 = pkin(4) * t185 - t116 * t151 + t117 * t138 + t163;
t160 = t128 * pkin(2) + t118 * pkin(9) + t167;
t159 = pkin(4) * t187 - t106 * t151 + t107 * t138 + t160;
t156 = cos(qJ(6));
t152 = sin(qJ(6));
t103 = t117 * t140 + t126 * t139;
t99 = t109 * t140 + t119 * t139;
t97 = t107 * t140 + t118 * t139;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t158 - g(2) * t155, t168, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t130 - g(2) * t128 - g(3) * t183, -g(1) * t129 - g(2) * t127 - g(3) * t182, -g(3) * t150 - t168 * t147, -g(1) * t174 - g(2) * t167 - g(3) * t188, 0, 0, 0, 0, 0, 0, -g(1) * t109 - g(2) * t107 - g(3) * t117, t165, -g(1) * t119 - g(2) * t118 - g(3) * t126, -g(1) * t164 - g(2) * t160 - g(3) * t163, 0, 0, 0, 0, 0, 0, -g(1) * (t109 * t157 + t186) - g(2) * (t107 * t157 + t187) - g(3) * (t117 * t157 + t185) -g(1) * (-t109 * t153 + t119 * t157) - g(2) * (-t107 * t153 + t118 * t157) - g(3) * (-t117 * t153 + t126 * t157) -t165, -g(1) * (t109 * pkin(3) + t108 * pkin(10) + t164) - g(2) * (t107 * pkin(3) + t106 * pkin(10) + t160) - g(3) * (t117 * pkin(3) + t116 * pkin(10) + t163) 0, 0, 0, 0, 0, 0, -g(1) * t99 - g(2) * t97 - g(3) * t103, t166, -t165, -g(1) * t162 - g(2) * t159 - g(3) * t161, 0, 0, 0, 0, 0, 0, -g(1) * (t108 * t152 + t99 * t156) - g(2) * (t106 * t152 + t97 * t156) - g(3) * (t103 * t156 + t116 * t152) -g(1) * (t108 * t156 - t99 * t152) - g(2) * (t106 * t156 - t97 * t152) - g(3) * (-t103 * t152 + t116 * t156) -t166, -g(1) * (t99 * pkin(5) + t98 * pkin(11) + t162) - g(2) * (t97 * pkin(5) + t96 * pkin(11) + t159) - g(3) * (t103 * pkin(5) + t102 * pkin(11) + t161);];
U_reg  = t1;
