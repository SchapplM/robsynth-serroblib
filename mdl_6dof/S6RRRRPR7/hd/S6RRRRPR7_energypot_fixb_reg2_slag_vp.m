% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:34:07
% EndTime: 2019-03-09 22:34:08
% DurationCPUTime: 0.24s
% Computational Cost: add. (333->107), mult. (459->157), div. (0->0), fcn. (534->14), ass. (0->55)
t126 = -pkin(10) - pkin(9);
t125 = cos(qJ(1));
t147 = g(2) * t125;
t119 = sin(qJ(3));
t146 = t119 * pkin(3);
t117 = cos(pkin(6));
t145 = t117 * pkin(8) + pkin(7);
t123 = cos(qJ(3));
t106 = t123 * pkin(3) + pkin(2);
t116 = sin(pkin(6));
t120 = sin(qJ(2));
t144 = t116 * t120;
t121 = sin(qJ(1));
t143 = t116 * t121;
t142 = t116 * t123;
t124 = cos(qJ(2));
t141 = t116 * t124;
t140 = t116 * t125;
t139 = t117 * t119;
t138 = t119 * t121;
t137 = t121 * t120;
t136 = t121 * t124;
t135 = t125 * t120;
t134 = t125 * t124;
t115 = qJ(3) + qJ(4);
t133 = t125 * pkin(1) + pkin(8) * t143;
t111 = t121 * pkin(1);
t132 = -pkin(8) * t140 + t111;
t114 = -qJ(5) + t126;
t109 = cos(t115);
t98 = pkin(4) * t109 + t106;
t108 = sin(t115);
t99 = pkin(4) * t108 + t146;
t131 = t114 * t141 + t117 * t99 + t98 * t144 + t145;
t94 = t117 * t136 + t135;
t95 = -t117 * t137 + t134;
t130 = -t94 * t114 + t99 * t143 + t95 * t98 + t133;
t129 = g(1) * t121 - t147;
t107 = pkin(12) + t115;
t104 = sin(t107);
t105 = cos(t107);
t93 = t117 * t135 + t136;
t80 = t93 * t104 + t105 * t140;
t82 = t95 * t104 - t105 * t143;
t86 = t104 * t144 - t117 * t105;
t128 = g(1) * t82 + g(2) * t80 + g(3) * t86;
t92 = -t117 * t134 + t137;
t127 = t111 + t93 * t98 - t92 * t114 + (-pkin(8) - t99) * t140;
t79 = -g(1) * t94 - g(2) * t92 + g(3) * t141;
t122 = cos(qJ(6));
t118 = sin(qJ(6));
t87 = t117 * t104 + t105 * t144;
t83 = t104 * t143 + t95 * t105;
t81 = -t104 * t140 + t93 * t105;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t125 - g(2) * t121, t129, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t144, -t79, -g(3) * t117 - t129 * t116, -g(1) * t133 - g(2) * t132 - g(3) * t145, 0, 0, 0, 0, 0, 0, -g(1) * (t116 * t138 + t95 * t123) - g(2) * (-t119 * t140 + t93 * t123) - g(3) * (t120 * t142 + t139) -g(1) * (-t95 * t119 + t121 * t142) - g(2) * (-t93 * t119 - t123 * t140) - g(3) * (t117 * t123 - t119 * t144) t79, -g(1) * (t95 * pkin(2) + t94 * pkin(9) + t133) - g(2) * (t93 * pkin(2) + t92 * pkin(9) + t132) - g(3) * ((pkin(2) * t120 - pkin(9) * t124) * t116 + t145) 0, 0, 0, 0, 0, 0, -g(1) * (t108 * t143 + t95 * t109) - g(2) * (-t108 * t140 + t93 * t109) - g(3) * (t117 * t108 + t109 * t144) -g(1) * (-t95 * t108 + t109 * t143) - g(2) * (-t93 * t108 - t109 * t140) - g(3) * (-t108 * t144 + t117 * t109) t79, -g(1) * (t95 * t106 - t94 * t126 + t133) - g(2) * (t93 * t106 - t92 * t126 + t111) - g(3) * (pkin(3) * t139 + t145) + (-g(1) * pkin(3) * t138 - g(3) * (t106 * t120 + t124 * t126) - (-pkin(8) - t146) * t147) * t116, 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t87, t128, t79, -g(1) * t130 - g(2) * t127 - g(3) * t131, 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t118 + t83 * t122) - g(2) * (t92 * t118 + t81 * t122) - g(3) * (-t118 * t141 + t87 * t122) -g(1) * (-t83 * t118 + t94 * t122) - g(2) * (-t81 * t118 + t92 * t122) - g(3) * (-t87 * t118 - t122 * t141) -t128, -g(1) * (t83 * pkin(5) + t82 * pkin(11) + t130) - g(2) * (t81 * pkin(5) + t80 * pkin(11) + t127) - g(3) * (t87 * pkin(5) + t86 * pkin(11) + t131);];
U_reg  = t1;
