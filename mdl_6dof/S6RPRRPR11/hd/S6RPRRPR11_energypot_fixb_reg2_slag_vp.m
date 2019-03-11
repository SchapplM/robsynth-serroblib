% Calculate inertial parameters regressor of potential energy for
% S6RPRRPR11
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:44:30
% EndTime: 2019-03-09 05:44:30
% DurationCPUTime: 0.40s
% Computational Cost: add. (539->111), mult. (1376->169), div. (0->0), fcn. (1753->16), ass. (0->66)
t140 = cos(pkin(6));
t134 = sin(pkin(12));
t145 = cos(qJ(1));
t170 = t145 * t134;
t138 = cos(pkin(12));
t144 = sin(qJ(1));
t171 = t144 * t138;
t117 = -t140 * t171 - t170;
t135 = sin(pkin(7));
t139 = cos(pkin(7));
t136 = sin(pkin(6));
t175 = t136 * t144;
t158 = -t117 * t135 + t139 * t175;
t176 = t136 * t138;
t157 = t135 * t176 - t140 * t139;
t181 = cos(qJ(3));
t180 = cos(qJ(4));
t179 = t140 * qJ(2) + pkin(8);
t177 = t134 * t136;
t174 = t136 * t145;
t172 = t144 * t134;
t169 = t145 * t138;
t168 = t145 * pkin(1) + qJ(2) * t175;
t133 = sin(pkin(13));
t165 = pkin(5) * t133 + pkin(10);
t164 = t135 * t181;
t163 = t139 * t181;
t162 = t136 * t164;
t161 = g(1) * t144 - g(2) * t145;
t160 = t144 * pkin(1) - qJ(2) * t174;
t115 = t140 * t169 - t172;
t159 = t115 * t135 + t139 * t174;
t116 = t140 * t170 + t171;
t143 = sin(qJ(3));
t102 = t116 * t181 + (t115 * t139 - t135 * t174) * t143;
t142 = sin(qJ(4));
t93 = t102 * t142 + t159 * t180;
t118 = -t140 * t172 + t169;
t104 = t118 * t181 + (t117 * t139 + t135 * t175) * t143;
t95 = t104 * t142 - t158 * t180;
t109 = t140 * t135 * t143 + (t138 * t139 * t143 + t181 * t134) * t136;
t99 = t109 * t142 + t157 * t180;
t156 = g(1) * t95 + g(2) * t93 + g(3) * t99;
t101 = -t115 * t163 + t116 * t143 + t145 * t162;
t103 = -t117 * t163 + t118 * t143 - t144 * t162;
t108 = -t140 * t164 + t143 * t177 - t163 * t176;
t155 = g(1) * t103 + g(2) * t101 + g(3) * t108;
t154 = t118 * pkin(2) + t158 * pkin(9) + t168;
t153 = pkin(2) * t177 - t157 * pkin(9) + t179;
t152 = t104 * pkin(3) + t154;
t151 = t109 * pkin(3) + t153;
t150 = t103 * pkin(10) + t152;
t149 = t108 * pkin(10) + t151;
t148 = t116 * pkin(2) - t159 * pkin(9) + t160;
t147 = t102 * pkin(3) + t148;
t146 = t101 * pkin(10) + t147;
t141 = -pkin(11) - qJ(5);
t137 = cos(pkin(13));
t132 = pkin(13) + qJ(6);
t128 = cos(t132);
t127 = sin(t132);
t126 = t137 * pkin(5) + pkin(4);
t100 = t109 * t180 - t157 * t142;
t96 = t104 * t180 + t158 * t142;
t94 = t102 * t180 - t159 * t142;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t145 - g(2) * t144, t161, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t118 - g(2) * t116 - g(3) * t177, -g(1) * t117 - g(2) * t115 - g(3) * t176, -g(3) * t140 - t161 * t136, -g(1) * t168 - g(2) * t160 - g(3) * t179, 0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t102 - g(3) * t109, t155, -g(1) * t158 + g(2) * t159 + g(3) * t157, -g(1) * t154 - g(2) * t148 - g(3) * t153, 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t100, t156, -t155, -g(1) * t150 - g(2) * t146 - g(3) * t149, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t133 + t96 * t137) - g(2) * (t101 * t133 + t94 * t137) - g(3) * (t100 * t137 + t108 * t133) -g(1) * (t103 * t137 - t96 * t133) - g(2) * (t101 * t137 - t94 * t133) - g(3) * (-t100 * t133 + t108 * t137) -t156, -g(1) * (t96 * pkin(4) + t95 * qJ(5) + t150) - g(2) * (t94 * pkin(4) + t93 * qJ(5) + t146) - g(3) * (t100 * pkin(4) + t99 * qJ(5) + t149) 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t127 + t96 * t128) - g(2) * (t101 * t127 + t94 * t128) - g(3) * (t100 * t128 + t108 * t127) -g(1) * (t103 * t128 - t96 * t127) - g(2) * (t101 * t128 - t94 * t127) - g(3) * (-t100 * t127 + t108 * t128) -t156, -g(1) * (t165 * t103 + t96 * t126 - t95 * t141 + t152) - g(2) * (t165 * t101 + t94 * t126 - t93 * t141 + t147) - g(3) * (t100 * t126 + t165 * t108 - t99 * t141 + t151);];
U_reg  = t1;
