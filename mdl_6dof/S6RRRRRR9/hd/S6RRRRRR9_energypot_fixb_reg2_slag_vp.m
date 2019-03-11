% Calculate inertial parameters regressor of potential energy for
% S6RRRRRR9
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:39:12
% EndTime: 2019-03-10 05:39:12
% DurationCPUTime: 0.37s
% Computational Cost: add. (539->111), mult. (1376->169), div. (0->0), fcn. (1753->16), ass. (0->66)
t136 = cos(pkin(6));
t140 = sin(qJ(2));
t144 = cos(qJ(1));
t170 = t144 * t140;
t141 = sin(qJ(1));
t143 = cos(qJ(2));
t171 = t141 * t143;
t117 = -t136 * t171 - t170;
t133 = sin(pkin(7));
t135 = cos(pkin(7));
t134 = sin(pkin(6));
t176 = t134 * t141;
t158 = -t117 * t133 + t135 * t176;
t175 = t134 * t143;
t157 = t133 * t175 - t136 * t135;
t181 = cos(qJ(3));
t180 = cos(qJ(4));
t179 = t136 * pkin(9) + pkin(8);
t177 = t134 * t140;
t174 = t134 * t144;
t172 = t141 * t140;
t169 = t144 * t143;
t168 = t144 * pkin(1) + pkin(9) * t176;
t137 = sin(qJ(5));
t165 = pkin(5) * t137 + pkin(11);
t164 = t133 * t181;
t163 = t135 * t181;
t162 = t134 * t164;
t161 = t141 * pkin(1) - pkin(9) * t174;
t160 = g(1) * t141 - g(2) * t144;
t115 = t136 * t169 - t172;
t159 = t115 * t133 + t135 * t174;
t116 = t136 * t170 + t171;
t139 = sin(qJ(3));
t102 = t116 * t181 + (t115 * t135 - t133 * t174) * t139;
t138 = sin(qJ(4));
t93 = t102 * t138 + t159 * t180;
t118 = -t136 * t172 + t169;
t104 = t118 * t181 + (t117 * t135 + t133 * t176) * t139;
t95 = t104 * t138 - t158 * t180;
t109 = t136 * t133 * t139 + (t135 * t139 * t143 + t140 * t181) * t134;
t99 = t109 * t138 + t157 * t180;
t156 = g(1) * t95 + g(2) * t93 + g(3) * t99;
t101 = -t115 * t163 + t116 * t139 + t144 * t162;
t103 = -t117 * t163 + t118 * t139 - t141 * t162;
t108 = -t136 * t164 + t139 * t177 - t163 * t175;
t155 = g(1) * t103 + g(2) * t101 + g(3) * t108;
t154 = t118 * pkin(2) + t158 * pkin(10) + t168;
t153 = pkin(2) * t177 - t157 * pkin(10) + t179;
t152 = t104 * pkin(3) + t154;
t151 = t109 * pkin(3) + t153;
t150 = t103 * pkin(11) + t152;
t149 = t108 * pkin(11) + t151;
t148 = t116 * pkin(2) - pkin(10) * t159 + t161;
t147 = t102 * pkin(3) + t148;
t146 = t101 * pkin(11) + t147;
t145 = -pkin(13) - pkin(12);
t142 = cos(qJ(5));
t132 = qJ(5) + qJ(6);
t128 = cos(t132);
t127 = sin(t132);
t126 = t142 * pkin(5) + pkin(4);
t100 = t109 * t180 - t138 * t157;
t96 = t104 * t180 + t138 * t158;
t94 = t102 * t180 - t138 * t159;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t144 - g(2) * t141, t160, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t118 - g(2) * t116 - g(3) * t177, -g(1) * t117 - g(2) * t115 - g(3) * t175, -g(3) * t136 - t134 * t160, -g(1) * t168 - g(2) * t161 - g(3) * t179, 0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t102 - g(3) * t109, t155, -g(1) * t158 + g(2) * t159 + g(3) * t157, -g(1) * t154 - g(2) * t148 - g(3) * t153, 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t100, t156, -t155, -g(1) * t150 - g(2) * t146 - g(3) * t149, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t137 + t96 * t142) - g(2) * (t101 * t137 + t94 * t142) - g(3) * (t100 * t142 + t108 * t137) -g(1) * (t103 * t142 - t96 * t137) - g(2) * (t101 * t142 - t94 * t137) - g(3) * (-t100 * t137 + t108 * t142) -t156, -g(1) * (t96 * pkin(4) + t95 * pkin(12) + t150) - g(2) * (t94 * pkin(4) + t93 * pkin(12) + t146) - g(3) * (t100 * pkin(4) + t99 * pkin(12) + t149) 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t127 + t96 * t128) - g(2) * (t101 * t127 + t94 * t128) - g(3) * (t100 * t128 + t108 * t127) -g(1) * (t103 * t128 - t96 * t127) - g(2) * (t101 * t128 - t94 * t127) - g(3) * (-t100 * t127 + t108 * t128) -t156, -g(1) * (t103 * t165 + t96 * t126 - t95 * t145 + t152) - g(2) * (t101 * t165 + t94 * t126 - t93 * t145 + t147) - g(3) * (t100 * t126 + t108 * t165 - t99 * t145 + t151);];
U_reg  = t1;
