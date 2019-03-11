% Calculate inertial parameters regressor of potential energy for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:13:57
% EndTime: 2019-03-09 14:13:58
% DurationCPUTime: 0.23s
% Computational Cost: add. (333->107), mult. (459->154), div. (0->0), fcn. (534->14), ass. (0->54)
t123 = sin(qJ(1));
t146 = g(1) * t123;
t126 = cos(qJ(1));
t145 = g(2) * t126;
t116 = sin(pkin(12));
t144 = t116 * pkin(3);
t120 = -pkin(9) - qJ(3);
t118 = cos(pkin(12));
t106 = t118 * pkin(3) + pkin(2);
t119 = cos(pkin(6));
t143 = t119 * pkin(8) + pkin(7);
t117 = sin(pkin(6));
t122 = sin(qJ(2));
t142 = t117 * t122;
t141 = t117 * t123;
t125 = cos(qJ(2));
t140 = t117 * t125;
t139 = t117 * t126;
t138 = t119 * t116;
t137 = t123 * t122;
t136 = t123 * t125;
t135 = t126 * t122;
t134 = t126 * t125;
t133 = t126 * pkin(1) + pkin(8) * t141;
t115 = pkin(12) + qJ(4);
t112 = t123 * pkin(1);
t132 = -pkin(8) * t139 + t112;
t114 = -pkin(10) + t120;
t108 = cos(t115);
t97 = pkin(4) * t108 + t106;
t107 = sin(t115);
t99 = pkin(4) * t107 + t144;
t131 = t114 * t140 + t119 * t99 + t97 * t142 + t143;
t94 = t119 * t136 + t135;
t95 = -t119 * t137 + t134;
t130 = -t94 * t114 + t99 * t141 + t95 * t97 + t133;
t129 = -t145 + t146;
t109 = qJ(5) + t115;
t104 = sin(t109);
t105 = cos(t109);
t93 = t119 * t135 + t136;
t80 = t93 * t104 + t105 * t139;
t82 = t95 * t104 - t105 * t141;
t86 = t104 * t142 - t119 * t105;
t128 = g(1) * t82 + g(2) * t80 + g(3) * t86;
t92 = -t119 * t134 + t137;
t127 = t112 + t93 * t97 - t92 * t114 + (-pkin(8) - t99) * t139;
t79 = -g(1) * t94 - g(2) * t92 + g(3) * t140;
t124 = cos(qJ(6));
t121 = sin(qJ(6));
t87 = t119 * t104 + t105 * t142;
t83 = t104 * t141 + t95 * t105;
t81 = -t104 * t139 + t93 * t105;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t126 - g(2) * t123, t129, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t142, -t79, -g(3) * t119 - t117 * t129, -g(1) * t133 - g(2) * t132 - g(3) * t143, 0, 0, 0, 0, 0, 0, -g(1) * (t116 * t141 + t95 * t118) - g(2) * (-t116 * t139 + t93 * t118) - g(3) * (t118 * t142 + t138) -g(1) * (-t95 * t116 + t118 * t141) - g(2) * (-t93 * t116 - t118 * t139) - g(3) * (-t116 * t142 + t119 * t118) t79, -g(1) * (t95 * pkin(2) + t94 * qJ(3) + t133) - g(2) * (t93 * pkin(2) + t92 * qJ(3) + t132) - g(3) * ((pkin(2) * t122 - qJ(3) * t125) * t117 + t143) 0, 0, 0, 0, 0, 0, -g(1) * (t107 * t141 + t95 * t108) - g(2) * (-t107 * t139 + t93 * t108) - g(3) * (t119 * t107 + t108 * t142) -g(1) * (-t95 * t107 + t108 * t141) - g(2) * (-t93 * t107 - t108 * t139) - g(3) * (-t107 * t142 + t119 * t108) t79, -g(1) * (t95 * t106 - t94 * t120 + t133) - g(2) * (t93 * t106 - t92 * t120 + t112) - g(3) * (pkin(3) * t138 + t143) + (-t144 * t146 - g(3) * (t106 * t122 + t120 * t125) - (-pkin(8) - t144) * t145) * t117, 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t87, t128, t79, -g(1) * t130 - g(2) * t127 - g(3) * t131, 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t121 + t83 * t124) - g(2) * (t92 * t121 + t81 * t124) - g(3) * (-t121 * t140 + t87 * t124) -g(1) * (-t83 * t121 + t94 * t124) - g(2) * (-t81 * t121 + t92 * t124) - g(3) * (-t87 * t121 - t124 * t140) -t128, -g(1) * (t83 * pkin(5) + t82 * pkin(11) + t130) - g(2) * (t81 * pkin(5) + t80 * pkin(11) + t127) - g(3) * (t87 * pkin(5) + t86 * pkin(11) + t131);];
U_reg  = t1;
