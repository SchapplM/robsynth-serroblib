% Calculate inertial parameters regressor of potential energy for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:59:59
% EndTime: 2019-03-09 00:59:59
% DurationCPUTime: 0.37s
% Computational Cost: add. (502->110), mult. (1179->171), div. (0->0), fcn. (1486->16), ass. (0->64)
t130 = sin(pkin(7));
t133 = cos(pkin(7));
t134 = cos(pkin(6));
t131 = sin(pkin(6));
t141 = cos(qJ(2));
t163 = t131 * t141;
t110 = -t130 * t163 + t134 * t133;
t129 = sin(pkin(13));
t132 = cos(pkin(13));
t138 = sin(qJ(2));
t160 = t134 * t141;
t113 = -t129 * t160 - t132 * t138;
t165 = t131 * t133;
t103 = -t113 * t130 + t129 * t165;
t172 = cos(qJ(3));
t111 = -t129 * t138 + t132 * t160;
t102 = -t111 * t130 - t132 * t165;
t136 = sin(qJ(4));
t171 = t102 * t136;
t170 = t103 * t136;
t169 = t110 * t136;
t167 = t129 * t131;
t166 = t131 * t132;
t164 = t131 * t138;
t161 = t134 * t138;
t159 = t132 * pkin(1) + pkin(8) * t167;
t158 = t134 * pkin(8) + qJ(1);
t155 = t130 * t172;
t154 = t133 * t172;
t153 = t131 * t155;
t152 = t129 * pkin(1) - pkin(8) * t166;
t151 = g(1) * t129 - g(2) * t132;
t128 = qJ(4) + qJ(5);
t123 = sin(t128);
t124 = cos(t128);
t112 = t129 * t141 + t132 * t161;
t137 = sin(qJ(3));
t91 = t112 * t172 + (t111 * t133 - t130 * t166) * t137;
t80 = -t102 * t124 + t91 * t123;
t114 = -t129 * t161 + t132 * t141;
t93 = t114 * t172 + (t113 * t133 + t130 * t167) * t137;
t82 = -t103 * t124 + t93 * t123;
t101 = t134 * t130 * t137 + (t133 * t137 * t141 + t172 * t138) * t131;
t86 = t101 * t123 - t110 * t124;
t150 = g(1) * t82 + g(2) * t80 + g(3) * t86;
t100 = -t134 * t155 + t137 * t164 - t154 * t163;
t90 = -t111 * t154 + t112 * t137 + t132 * t153;
t92 = -t113 * t154 + t114 * t137 - t129 * t153;
t149 = g(1) * t92 + g(2) * t90 + g(3) * t100;
t148 = t114 * pkin(2) + t103 * pkin(9) + t159;
t147 = pkin(2) * t164 + t110 * pkin(9) + t158;
t140 = cos(qJ(4));
t122 = t140 * pkin(4) + pkin(3);
t142 = -pkin(11) - pkin(10);
t146 = pkin(4) * t170 + t93 * t122 - t92 * t142 + t148;
t145 = pkin(4) * t169 - t100 * t142 + t101 * t122 + t147;
t144 = t112 * pkin(2) + t102 * pkin(9) + t152;
t143 = pkin(4) * t171 + t91 * t122 - t90 * t142 + t144;
t139 = cos(qJ(6));
t135 = sin(qJ(6));
t87 = t101 * t124 + t110 * t123;
t83 = t103 * t123 + t93 * t124;
t81 = t102 * t123 + t91 * t124;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t132 - g(2) * t129, t151, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t114 - g(2) * t112 - g(3) * t164, -g(1) * t113 - g(2) * t111 - g(3) * t163, -g(3) * t134 - t151 * t131, -g(1) * t159 - g(2) * t152 - g(3) * t158, 0, 0, 0, 0, 0, 0, -g(1) * t93 - g(2) * t91 - g(3) * t101, t149, -g(1) * t103 - g(2) * t102 - g(3) * t110, -g(1) * t148 - g(2) * t144 - g(3) * t147, 0, 0, 0, 0, 0, 0, -g(1) * (t93 * t140 + t170) - g(2) * (t91 * t140 + t171) - g(3) * (t101 * t140 + t169) -g(1) * (t103 * t140 - t93 * t136) - g(2) * (t102 * t140 - t91 * t136) - g(3) * (-t101 * t136 + t110 * t140) -t149, -g(1) * (t93 * pkin(3) + t92 * pkin(10) + t148) - g(2) * (t91 * pkin(3) + t90 * pkin(10) + t144) - g(3) * (t101 * pkin(3) + t100 * pkin(10) + t147) 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t87, t150, -t149, -g(1) * t146 - g(2) * t143 - g(3) * t145, 0, 0, 0, 0, 0, 0, -g(1) * (t92 * t135 + t83 * t139) - g(2) * (t90 * t135 + t81 * t139) - g(3) * (t100 * t135 + t87 * t139) -g(1) * (-t83 * t135 + t92 * t139) - g(2) * (-t81 * t135 + t90 * t139) - g(3) * (t100 * t139 - t87 * t135) -t150, -g(1) * (t83 * pkin(5) + t82 * pkin(12) + t146) - g(2) * (t81 * pkin(5) + t80 * pkin(12) + t143) - g(3) * (t87 * pkin(5) + t86 * pkin(12) + t145);];
U_reg  = t1;
