% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR14
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR14_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:25:28
% EndTime: 2019-03-10 00:25:29
% DurationCPUTime: 0.37s
% Computational Cost: add. (539->111), mult. (1376->169), div. (0->0), fcn. (1753->16), ass. (0->66)
t138 = cos(pkin(6));
t142 = sin(qJ(2));
t145 = cos(qJ(1));
t170 = t145 * t142;
t143 = sin(qJ(1));
t144 = cos(qJ(2));
t171 = t143 * t144;
t117 = -t138 * t171 - t170;
t134 = sin(pkin(7));
t137 = cos(pkin(7));
t135 = sin(pkin(6));
t176 = t135 * t143;
t158 = -t117 * t134 + t137 * t176;
t175 = t135 * t144;
t157 = t134 * t175 - t137 * t138;
t181 = cos(qJ(3));
t180 = cos(qJ(4));
t179 = pkin(9) * t138 + pkin(8);
t177 = t135 * t142;
t174 = t135 * t145;
t172 = t143 * t142;
t169 = t145 * t144;
t168 = pkin(1) * t145 + pkin(9) * t176;
t133 = sin(pkin(13));
t165 = pkin(5) * t133 + pkin(11);
t164 = t134 * t181;
t163 = t137 * t181;
t162 = t135 * t164;
t161 = pkin(1) * t143 - pkin(9) * t174;
t160 = g(1) * t143 - g(2) * t145;
t115 = t138 * t169 - t172;
t159 = t115 * t134 + t137 * t174;
t116 = t138 * t170 + t171;
t141 = sin(qJ(3));
t102 = t116 * t181 + (t115 * t137 - t134 * t174) * t141;
t140 = sin(qJ(4));
t93 = t102 * t140 + t159 * t180;
t118 = -t138 * t172 + t169;
t104 = t118 * t181 + (t117 * t137 + t134 * t176) * t141;
t95 = t104 * t140 - t158 * t180;
t109 = t138 * t134 * t141 + (t137 * t141 * t144 + t142 * t181) * t135;
t99 = t109 * t140 + t157 * t180;
t156 = g(1) * t95 + g(2) * t93 + g(3) * t99;
t101 = -t115 * t163 + t116 * t141 + t145 * t162;
t103 = -t117 * t163 + t118 * t141 - t143 * t162;
t108 = -t138 * t164 + t141 * t177 - t163 * t175;
t155 = g(1) * t103 + g(2) * t101 + g(3) * t108;
t154 = t118 * pkin(2) + pkin(10) * t158 + t168;
t153 = pkin(2) * t177 - pkin(10) * t157 + t179;
t152 = pkin(3) * t104 + t154;
t151 = pkin(3) * t109 + t153;
t150 = pkin(11) * t103 + t152;
t149 = pkin(11) * t108 + t151;
t148 = pkin(2) * t116 - pkin(10) * t159 + t161;
t147 = pkin(3) * t102 + t148;
t146 = t101 * pkin(11) + t147;
t139 = -pkin(12) - qJ(5);
t136 = cos(pkin(13));
t132 = pkin(13) + qJ(6);
t128 = cos(t132);
t127 = sin(t132);
t126 = pkin(5) * t136 + pkin(4);
t100 = t109 * t180 - t140 * t157;
t96 = t104 * t180 + t140 * t158;
t94 = t102 * t180 - t140 * t159;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t145 - g(2) * t143, t160, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t118 - g(2) * t116 - g(3) * t177, -g(1) * t117 - g(2) * t115 - g(3) * t175, -g(3) * t138 - t135 * t160, -g(1) * t168 - g(2) * t161 - g(3) * t179, 0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t102 - g(3) * t109, t155, -g(1) * t158 + g(2) * t159 + g(3) * t157, -g(1) * t154 - g(2) * t148 - g(3) * t153, 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t100, t156, -t155, -g(1) * t150 - g(2) * t146 - g(3) * t149, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t133 + t136 * t96) - g(2) * (t101 * t133 + t136 * t94) - g(3) * (t100 * t136 + t108 * t133) -g(1) * (t103 * t136 - t133 * t96) - g(2) * (t101 * t136 - t133 * t94) - g(3) * (-t100 * t133 + t108 * t136) -t156, -g(1) * (pkin(4) * t96 + qJ(5) * t95 + t150) - g(2) * (t94 * pkin(4) + t93 * qJ(5) + t146) - g(3) * (pkin(4) * t100 + qJ(5) * t99 + t149) 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t127 + t128 * t96) - g(2) * (t101 * t127 + t128 * t94) - g(3) * (t100 * t128 + t108 * t127) -g(1) * (t103 * t128 - t127 * t96) - g(2) * (t101 * t128 - t127 * t94) - g(3) * (-t100 * t127 + t108 * t128) -t156, -g(1) * (t103 * t165 + t96 * t126 - t95 * t139 + t152) - g(2) * (t101 * t165 + t94 * t126 - t93 * t139 + t147) - g(3) * (t100 * t126 + t108 * t165 - t99 * t139 + t151);];
U_reg  = t1;
