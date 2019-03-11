% Calculate inertial parameters regressor of potential energy for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:41:11
% EndTime: 2019-03-09 06:41:11
% DurationCPUTime: 0.34s
% Computational Cost: add. (527->100), mult. (1376->151), div. (0->0), fcn. (1753->14), ass. (0->65)
t134 = cos(pkin(6));
t129 = sin(pkin(12));
t141 = cos(qJ(1));
t166 = t141 * t129;
t132 = cos(pkin(12));
t139 = sin(qJ(1));
t167 = t139 * t132;
t116 = -t134 * t167 - t166;
t130 = sin(pkin(7));
t133 = cos(pkin(7));
t131 = sin(pkin(6));
t171 = t131 * t139;
t154 = -t116 * t130 + t133 * t171;
t172 = t131 * t132;
t153 = t130 * t172 - t134 * t133;
t177 = cos(qJ(3));
t176 = cos(qJ(4));
t175 = t134 * qJ(2) + pkin(8);
t173 = t129 * t131;
t170 = t131 * t141;
t168 = t139 * t129;
t165 = t141 * t132;
t164 = t141 * pkin(1) + qJ(2) * t171;
t136 = sin(qJ(5));
t161 = pkin(5) * t136 + pkin(10);
t160 = t130 * t177;
t159 = t133 * t177;
t158 = t131 * t160;
t157 = g(1) * t139 - g(2) * t141;
t156 = t139 * pkin(1) - qJ(2) * t170;
t114 = t134 * t165 - t168;
t155 = t114 * t130 + t133 * t170;
t115 = t134 * t166 + t167;
t138 = sin(qJ(3));
t101 = t115 * t177 + (t114 * t133 - t130 * t170) * t138;
t137 = sin(qJ(4));
t92 = t101 * t137 + t155 * t176;
t117 = -t134 * t168 + t165;
t103 = t117 * t177 + (t116 * t133 + t130 * t171) * t138;
t94 = t103 * t137 - t154 * t176;
t108 = t134 * t130 * t138 + (t132 * t133 * t138 + t177 * t129) * t131;
t98 = t108 * t137 + t153 * t176;
t152 = g(1) * t94 + g(2) * t92 + g(3) * t98;
t100 = -t114 * t159 + t115 * t138 + t141 * t158;
t102 = -t116 * t159 + t117 * t138 - t139 * t158;
t107 = -t134 * t160 + t138 * t173 - t159 * t172;
t151 = g(1) * t102 + g(2) * t100 + g(3) * t107;
t150 = t117 * pkin(2) + t154 * pkin(9) + t164;
t149 = pkin(2) * t173 - t153 * pkin(9) + t175;
t148 = t103 * pkin(3) + t150;
t147 = t108 * pkin(3) + t149;
t146 = t102 * pkin(10) + t148;
t145 = t107 * pkin(10) + t147;
t144 = t115 * pkin(2) - t155 * pkin(9) + t156;
t143 = t101 * pkin(3) + t144;
t142 = t100 * pkin(10) + t143;
t140 = cos(qJ(5));
t135 = -qJ(6) - pkin(11);
t125 = t140 * pkin(5) + pkin(4);
t99 = t108 * t176 - t153 * t137;
t95 = t103 * t176 + t154 * t137;
t93 = t101 * t176 - t155 * t137;
t90 = -g(1) * (t102 * t136 + t95 * t140) - g(2) * (t100 * t136 + t93 * t140) - g(3) * (t107 * t136 + t99 * t140);
t89 = -g(1) * (t102 * t140 - t95 * t136) - g(2) * (t100 * t140 - t93 * t136) - g(3) * (t107 * t140 - t99 * t136);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t141 - g(2) * t139, t157, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t117 - g(2) * t115 - g(3) * t173, -g(1) * t116 - g(2) * t114 - g(3) * t172, -g(3) * t134 - t157 * t131, -g(1) * t164 - g(2) * t156 - g(3) * t175, 0, 0, 0, 0, 0, 0, -g(1) * t103 - g(2) * t101 - g(3) * t108, t151, -g(1) * t154 + g(2) * t155 + g(3) * t153, -g(1) * t150 - g(2) * t144 - g(3) * t149, 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t99, t152, -t151, -g(1) * t146 - g(2) * t142 - g(3) * t145, 0, 0, 0, 0, 0, 0, t90, t89, -t152, -g(1) * (t95 * pkin(4) + t94 * pkin(11) + t146) - g(2) * (t93 * pkin(4) + t92 * pkin(11) + t142) - g(3) * (t99 * pkin(4) + t98 * pkin(11) + t145) 0, 0, 0, 0, 0, 0, t90, t89, -t152, -g(1) * (t161 * t102 + t95 * t125 - t94 * t135 + t148) - g(2) * (t161 * t100 + t93 * t125 - t92 * t135 + t143) - g(3) * (t161 * t107 + t99 * t125 - t98 * t135 + t147);];
U_reg  = t1;
