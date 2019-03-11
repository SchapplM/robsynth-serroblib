% Calculate inertial parameters regressor of potential energy for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:53:36
% EndTime: 2019-03-10 02:53:36
% DurationCPUTime: 0.34s
% Computational Cost: add. (527->100), mult. (1376->151), div. (0->0), fcn. (1753->14), ass. (0->65)
t132 = cos(pkin(6));
t138 = sin(qJ(1));
t140 = cos(qJ(2));
t166 = t138 * t140;
t137 = sin(qJ(2));
t141 = cos(qJ(1));
t168 = t137 * t141;
t116 = -t132 * t166 - t168;
t129 = sin(pkin(7));
t131 = cos(pkin(7));
t130 = sin(pkin(6));
t172 = t130 * t138;
t154 = -t116 * t129 + t131 * t172;
t171 = t130 * t140;
t153 = t129 * t171 - t131 * t132;
t177 = cos(qJ(3));
t176 = cos(qJ(4));
t175 = t132 * pkin(9) + pkin(8);
t173 = t130 * t137;
t170 = t130 * t141;
t167 = t138 * t137;
t165 = t140 * t141;
t164 = t141 * pkin(1) + pkin(9) * t172;
t134 = sin(qJ(5));
t161 = pkin(5) * t134 + pkin(11);
t160 = t129 * t177;
t159 = t131 * t177;
t158 = t130 * t160;
t157 = t138 * pkin(1) - pkin(9) * t170;
t156 = g(1) * t138 - g(2) * t141;
t114 = t132 * t165 - t167;
t155 = t114 * t129 + t131 * t170;
t115 = t132 * t168 + t166;
t136 = sin(qJ(3));
t101 = t115 * t177 + (t114 * t131 - t129 * t170) * t136;
t135 = sin(qJ(4));
t92 = t101 * t135 + t155 * t176;
t117 = -t132 * t167 + t165;
t103 = t117 * t177 + (t116 * t131 + t129 * t172) * t136;
t94 = t103 * t135 - t154 * t176;
t108 = t132 * t129 * t136 + (t131 * t136 * t140 + t177 * t137) * t130;
t98 = t108 * t135 + t153 * t176;
t152 = g(1) * t94 + g(2) * t92 + g(3) * t98;
t100 = -t114 * t159 + t115 * t136 + t141 * t158;
t102 = -t116 * t159 + t117 * t136 - t138 * t158;
t107 = -t132 * t160 + t136 * t173 - t159 * t171;
t151 = g(1) * t102 + g(2) * t100 + g(3) * t107;
t150 = t117 * pkin(2) + t154 * pkin(10) + t164;
t149 = pkin(2) * t173 - t153 * pkin(10) + t175;
t148 = t103 * pkin(3) + t150;
t147 = t108 * pkin(3) + t149;
t146 = pkin(11) * t102 + t148;
t145 = pkin(11) * t107 + t147;
t144 = t115 * pkin(2) - t155 * pkin(10) + t157;
t143 = t101 * pkin(3) + t144;
t142 = t100 * pkin(11) + t143;
t139 = cos(qJ(5));
t133 = -qJ(6) - pkin(12);
t125 = pkin(5) * t139 + pkin(4);
t99 = t108 * t176 - t153 * t135;
t95 = t103 * t176 + t154 * t135;
t93 = t101 * t176 - t155 * t135;
t90 = -g(1) * (t102 * t134 + t139 * t95) - g(2) * (t100 * t134 + t139 * t93) - g(3) * (t107 * t134 + t139 * t99);
t89 = -g(1) * (t102 * t139 - t134 * t95) - g(2) * (t100 * t139 - t134 * t93) - g(3) * (t107 * t139 - t99 * t134);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t141 - g(2) * t138, t156, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t117 - g(2) * t115 - g(3) * t173, -g(1) * t116 - g(2) * t114 - g(3) * t171, -g(3) * t132 - t156 * t130, -g(1) * t164 - g(2) * t157 - g(3) * t175, 0, 0, 0, 0, 0, 0, -g(1) * t103 - g(2) * t101 - g(3) * t108, t151, -g(1) * t154 + g(2) * t155 + g(3) * t153, -g(1) * t150 - g(2) * t144 - g(3) * t149, 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t99, t152, -t151, -g(1) * t146 - g(2) * t142 - g(3) * t145, 0, 0, 0, 0, 0, 0, t90, t89, -t152, -g(1) * (pkin(4) * t95 + pkin(12) * t94 + t146) - g(2) * (t93 * pkin(4) + t92 * pkin(12) + t142) - g(3) * (pkin(4) * t99 + pkin(12) * t98 + t145) 0, 0, 0, 0, 0, 0, t90, t89, -t152, -g(1) * (t161 * t102 + t125 * t95 - t133 * t94 + t148) - g(2) * (t161 * t100 + t93 * t125 - t92 * t133 + t143) - g(3) * (t161 * t107 + t125 * t99 - t98 * t133 + t147);];
U_reg  = t1;
