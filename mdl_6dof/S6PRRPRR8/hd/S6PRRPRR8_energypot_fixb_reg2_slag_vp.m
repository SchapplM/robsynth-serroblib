% Calculate inertial parameters regressor of potential energy for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:36
% EndTime: 2019-03-08 22:40:37
% DurationCPUTime: 0.32s
% Computational Cost: add. (453->95), mult. (1177->148), div. (0->0), fcn. (1483->14), ass. (0->63)
t123 = sin(pkin(7));
t126 = cos(pkin(7));
t127 = cos(pkin(6));
t124 = sin(pkin(6));
t134 = cos(qJ(2));
t160 = t124 * t134;
t106 = -t123 * t160 + t127 * t126;
t122 = sin(pkin(12));
t125 = cos(pkin(12));
t131 = sin(qJ(2));
t156 = t127 * t134;
t109 = -t122 * t156 - t125 * t131;
t162 = t124 * t126;
t100 = -t109 * t123 + t122 * t162;
t167 = cos(qJ(3));
t165 = t122 * t124;
t130 = sin(qJ(3));
t164 = t123 * t130;
t163 = t124 * t125;
t161 = t124 * t131;
t159 = t126 * t130;
t157 = t127 * t131;
t155 = t125 * pkin(1) + pkin(8) * t165;
t154 = t127 * pkin(8) + qJ(1);
t151 = t123 * t167;
t150 = t126 * t167;
t149 = t124 * t151;
t148 = t122 * pkin(1) - pkin(8) * t163;
t147 = g(1) * t122 - g(2) * t125;
t107 = -t122 * t131 + t125 * t156;
t99 = -t107 * t123 - t125 * t162;
t129 = sin(qJ(5));
t133 = cos(qJ(5));
t108 = t122 * t134 + t125 * t157;
t87 = -t107 * t150 + t108 * t130 + t125 * t149;
t78 = t99 * t129 - t87 * t133;
t110 = -t122 * t157 + t125 * t134;
t89 = -t109 * t150 + t110 * t130 - t122 * t149;
t80 = t100 * t129 - t89 * t133;
t97 = -t127 * t151 + t130 * t161 - t150 * t160;
t91 = t106 * t129 - t97 * t133;
t146 = g(1) * t80 + g(2) * t78 + g(3) * t91;
t145 = g(1) * t89 + g(2) * t87 + g(3) * t97;
t88 = t107 * t159 + t108 * t167 - t163 * t164;
t90 = t110 * t167 + (t109 * t126 + t123 * t165) * t130;
t98 = t127 * t164 + (t167 * t131 + t134 * t159) * t124;
t144 = g(1) * t90 + g(2) * t88 + g(3) * t98;
t143 = t110 * pkin(2) + t100 * pkin(9) + t155;
t142 = pkin(2) * t161 + t106 * pkin(9) + t154;
t141 = t90 * pkin(3) + t89 * qJ(4) + t143;
t140 = t98 * pkin(3) + t97 * qJ(4) + t142;
t139 = t108 * pkin(2) + t99 * pkin(9) + t148;
t138 = t100 * pkin(4) + t90 * pkin(10) + t141;
t137 = t106 * pkin(4) + t98 * pkin(10) + t140;
t136 = t88 * pkin(3) + t87 * qJ(4) + t139;
t135 = t99 * pkin(4) + t88 * pkin(10) + t136;
t132 = cos(qJ(6));
t128 = sin(qJ(6));
t92 = t106 * t133 + t97 * t129;
t82 = -g(1) * t100 - g(2) * t99 - g(3) * t106;
t81 = t100 * t133 + t89 * t129;
t79 = t87 * t129 + t99 * t133;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t125 - g(2) * t122, t147, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t110 - g(2) * t108 - g(3) * t161, -g(1) * t109 - g(2) * t107 - g(3) * t160, -g(3) * t127 - t147 * t124, -g(1) * t155 - g(2) * t148 - g(3) * t154, 0, 0, 0, 0, 0, 0, -t144, t145, t82, -g(1) * t143 - g(2) * t139 - g(3) * t142, 0, 0, 0, 0, 0, 0, t82, t144, -t145, -g(1) * t141 - g(2) * t136 - g(3) * t140, 0, 0, 0, 0, 0, 0, -g(1) * t81 - g(2) * t79 - g(3) * t92, t146, -t144, -g(1) * t138 - g(2) * t135 - g(3) * t137, 0, 0, 0, 0, 0, 0, -g(1) * (t90 * t128 + t81 * t132) - g(2) * (t88 * t128 + t79 * t132) - g(3) * (t98 * t128 + t92 * t132) -g(1) * (-t81 * t128 + t90 * t132) - g(2) * (-t79 * t128 + t88 * t132) - g(3) * (-t92 * t128 + t98 * t132) -t146, -g(1) * (t81 * pkin(5) + t80 * pkin(11) + t138) - g(2) * (t79 * pkin(5) + t78 * pkin(11) + t135) - g(3) * (t92 * pkin(5) + t91 * pkin(11) + t137);];
U_reg  = t1;
