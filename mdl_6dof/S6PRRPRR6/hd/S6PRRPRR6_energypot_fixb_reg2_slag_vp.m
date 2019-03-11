% Calculate inertial parameters regressor of potential energy for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:27:37
% EndTime: 2019-03-08 22:27:37
% DurationCPUTime: 0.37s
% Computational Cost: add. (502->110), mult. (1179->171), div. (0->0), fcn. (1486->16), ass. (0->64)
t130 = sin(pkin(7));
t134 = cos(pkin(7));
t135 = cos(pkin(6));
t131 = sin(pkin(6));
t141 = cos(qJ(2));
t162 = t131 * t141;
t109 = -t130 * t162 + t135 * t134;
t129 = sin(pkin(12));
t133 = cos(pkin(12));
t139 = sin(qJ(2));
t159 = t135 * t141;
t112 = -t129 * t159 - t133 * t139;
t164 = t131 * t134;
t102 = -t112 * t130 + t129 * t164;
t171 = cos(qJ(3));
t110 = -t129 * t139 + t133 * t159;
t101 = -t110 * t130 - t133 * t164;
t128 = sin(pkin(13));
t170 = t101 * t128;
t169 = t102 * t128;
t168 = t109 * t128;
t166 = t129 * t131;
t165 = t131 * t133;
t163 = t131 * t139;
t160 = t135 * t139;
t158 = t133 * pkin(1) + pkin(8) * t166;
t157 = t135 * pkin(8) + qJ(1);
t154 = t130 * t171;
t153 = t134 * t171;
t152 = t131 * t154;
t151 = t129 * pkin(1) - pkin(8) * t165;
t150 = g(1) * t129 - g(2) * t133;
t127 = pkin(13) + qJ(5);
t122 = sin(t127);
t123 = cos(t127);
t111 = t129 * t141 + t133 * t160;
t138 = sin(qJ(3));
t90 = t111 * t171 + (t110 * t134 - t130 * t165) * t138;
t79 = -t101 * t123 + t90 * t122;
t113 = -t129 * t160 + t133 * t141;
t92 = t113 * t171 + (t112 * t134 + t130 * t166) * t138;
t81 = -t102 * t123 + t92 * t122;
t100 = t135 * t130 * t138 + (t134 * t138 * t141 + t171 * t139) * t131;
t85 = t100 * t122 - t109 * t123;
t149 = g(1) * t81 + g(2) * t79 + g(3) * t85;
t89 = -t110 * t153 + t111 * t138 + t133 * t152;
t91 = -t112 * t153 + t113 * t138 - t129 * t152;
t99 = -t135 * t154 + t138 * t163 - t153 * t162;
t148 = g(1) * t91 + g(2) * t89 + g(3) * t99;
t147 = t113 * pkin(2) + t102 * pkin(9) + t158;
t146 = pkin(2) * t163 + t109 * pkin(9) + t157;
t132 = cos(pkin(13));
t121 = t132 * pkin(4) + pkin(3);
t136 = -pkin(10) - qJ(4);
t145 = pkin(4) * t169 + t92 * t121 - t91 * t136 + t147;
t144 = pkin(4) * t168 + t100 * t121 - t99 * t136 + t146;
t143 = t111 * pkin(2) + t101 * pkin(9) + t151;
t142 = pkin(4) * t170 + t90 * t121 - t89 * t136 + t143;
t140 = cos(qJ(6));
t137 = sin(qJ(6));
t86 = t100 * t123 + t109 * t122;
t82 = t102 * t122 + t92 * t123;
t80 = t101 * t122 + t90 * t123;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t133 - g(2) * t129, t150, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t113 - g(2) * t111 - g(3) * t163, -g(1) * t112 - g(2) * t110 - g(3) * t162, -g(3) * t135 - t150 * t131, -g(1) * t158 - g(2) * t151 - g(3) * t157, 0, 0, 0, 0, 0, 0, -g(1) * t92 - g(2) * t90 - g(3) * t100, t148, -g(1) * t102 - g(2) * t101 - g(3) * t109, -g(1) * t147 - g(2) * t143 - g(3) * t146, 0, 0, 0, 0, 0, 0, -g(1) * (t92 * t132 + t169) - g(2) * (t90 * t132 + t170) - g(3) * (t100 * t132 + t168) -g(1) * (t102 * t132 - t92 * t128) - g(2) * (t101 * t132 - t90 * t128) - g(3) * (-t100 * t128 + t109 * t132) -t148, -g(1) * (t92 * pkin(3) + t91 * qJ(4) + t147) - g(2) * (t90 * pkin(3) + t89 * qJ(4) + t143) - g(3) * (t100 * pkin(3) + t99 * qJ(4) + t146) 0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t80 - g(3) * t86, t149, -t148, -g(1) * t145 - g(2) * t142 - g(3) * t144, 0, 0, 0, 0, 0, 0, -g(1) * (t91 * t137 + t82 * t140) - g(2) * (t89 * t137 + t80 * t140) - g(3) * (t99 * t137 + t86 * t140) -g(1) * (-t82 * t137 + t91 * t140) - g(2) * (-t80 * t137 + t89 * t140) - g(3) * (-t86 * t137 + t99 * t140) -t149, -g(1) * (t82 * pkin(5) + t81 * pkin(11) + t145) - g(2) * (t80 * pkin(5) + t79 * pkin(11) + t142) - g(3) * (t86 * pkin(5) + t85 * pkin(11) + t144);];
U_reg  = t1;
