% Calculate inertial parameters regressor of potential energy for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:34:24
% EndTime: 2019-03-09 00:34:24
% DurationCPUTime: 0.33s
% Computational Cost: add. (576->95), mult. (1520->149), div. (0->0), fcn. (1949->14), ass. (0->67)
t126 = sin(pkin(7));
t129 = cos(pkin(7));
t130 = cos(pkin(6));
t127 = sin(pkin(6));
t136 = cos(qJ(2));
t164 = t127 * t136;
t149 = t126 * t164 - t130 * t129;
t125 = sin(pkin(12));
t128 = cos(pkin(12));
t134 = sin(qJ(2));
t161 = t130 * t136;
t113 = -t125 * t161 - t128 * t134;
t166 = t127 * t129;
t150 = -t113 * t126 + t125 * t166;
t171 = cos(qJ(3));
t170 = cos(qJ(4));
t168 = t125 * t127;
t167 = t127 * t128;
t165 = t127 * t134;
t162 = t130 * t134;
t160 = t128 * pkin(1) + pkin(8) * t168;
t159 = t130 * pkin(8) + qJ(1);
t156 = t126 * t171;
t155 = t129 * t171;
t154 = t127 * t156;
t153 = t125 * pkin(1) - pkin(8) * t167;
t152 = g(1) * t125 - g(2) * t128;
t111 = -t125 * t134 + t128 * t161;
t151 = t111 * t126 + t128 * t166;
t131 = sin(qJ(5));
t135 = cos(qJ(5));
t132 = sin(qJ(4));
t112 = t125 * t136 + t128 * t162;
t133 = sin(qJ(3));
t95 = t112 * t171 + (t111 * t129 - t126 * t167) * t133;
t86 = -t132 * t151 + t170 * t95;
t94 = -t111 * t155 + t112 * t133 + t128 * t154;
t77 = t86 * t131 - t94 * t135;
t114 = -t125 * t162 + t128 * t136;
t97 = t114 * t171 + (t113 * t129 + t126 * t168) * t133;
t88 = t132 * t150 + t170 * t97;
t96 = -t113 * t155 + t114 * t133 - t125 * t154;
t79 = t88 * t131 - t96 * t135;
t104 = -t130 * t156 + t133 * t165 - t155 * t164;
t105 = t130 * t126 * t133 + (t129 * t133 * t136 + t134 * t171) * t127;
t99 = t105 * t170 - t132 * t149;
t83 = -t104 * t135 + t99 * t131;
t148 = g(1) * t79 + g(2) * t77 + g(3) * t83;
t85 = t95 * t132 + t151 * t170;
t87 = t97 * t132 - t150 * t170;
t98 = t105 * t132 + t149 * t170;
t147 = g(1) * t87 + g(2) * t85 + g(3) * t98;
t146 = g(1) * t96 + g(2) * t94 + g(3) * t104;
t145 = t114 * pkin(2) + t150 * pkin(9) + t160;
t144 = pkin(2) * t165 - t149 * pkin(9) + t159;
t143 = t97 * pkin(3) + t96 * pkin(10) + t145;
t142 = t105 * pkin(3) + t104 * pkin(10) + t144;
t141 = t112 * pkin(2) - pkin(9) * t151 + t153;
t140 = t88 * pkin(4) + t87 * pkin(11) + t143;
t139 = t99 * pkin(4) + t98 * pkin(11) + t142;
t138 = t95 * pkin(3) + t94 * pkin(10) + t141;
t137 = t86 * pkin(4) + t85 * pkin(11) + t138;
t84 = t104 * t131 + t99 * t135;
t80 = t96 * t131 + t88 * t135;
t78 = t94 * t131 + t86 * t135;
t75 = -g(1) * t80 - g(2) * t78 - g(3) * t84;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t128 - g(2) * t125, t152, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t114 - g(2) * t112 - g(3) * t165, -g(1) * t113 - g(2) * t111 - g(3) * t164, -g(3) * t130 - t127 * t152, -g(1) * t160 - g(2) * t153 - g(3) * t159, 0, 0, 0, 0, 0, 0, -g(1) * t97 - g(2) * t95 - g(3) * t105, t146, -g(1) * t150 + g(2) * t151 + g(3) * t149, -g(1) * t145 - g(2) * t141 - g(3) * t144, 0, 0, 0, 0, 0, 0, -g(1) * t88 - g(2) * t86 - g(3) * t99, t147, -t146, -g(1) * t143 - g(2) * t138 - g(3) * t142, 0, 0, 0, 0, 0, 0, t75, t148, -t147, -g(1) * t140 - g(2) * t137 - g(3) * t139, 0, 0, 0, 0, 0, 0, t75, -t147, -t148, -g(1) * (t80 * pkin(5) + t79 * qJ(6) + t140) - g(2) * (t78 * pkin(5) + t77 * qJ(6) + t137) - g(3) * (t84 * pkin(5) + t83 * qJ(6) + t139);];
U_reg  = t1;
