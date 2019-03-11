% Calculate inertial parameters regressor of potential energy for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:39:25
% EndTime: 2019-03-09 18:39:25
% DurationCPUTime: 0.24s
% Computational Cost: add. (333->107), mult. (459->157), div. (0->0), fcn. (534->14), ass. (0->55)
t126 = cos(qJ(1));
t147 = g(2) * t126;
t120 = sin(qJ(3));
t146 = t120 * pkin(3);
t118 = -qJ(4) - pkin(9);
t117 = cos(pkin(6));
t145 = t117 * pkin(8) + pkin(7);
t124 = cos(qJ(3));
t106 = t124 * pkin(3) + pkin(2);
t116 = sin(pkin(6));
t121 = sin(qJ(2));
t144 = t116 * t121;
t122 = sin(qJ(1));
t143 = t116 * t122;
t142 = t116 * t124;
t125 = cos(qJ(2));
t141 = t116 * t125;
t140 = t116 * t126;
t139 = t117 * t120;
t138 = t120 * t122;
t137 = t122 * t121;
t136 = t122 * t125;
t135 = t126 * t121;
t134 = t126 * t125;
t133 = t126 * pkin(1) + pkin(8) * t143;
t115 = qJ(3) + pkin(12);
t111 = t122 * pkin(1);
t132 = -pkin(8) * t140 + t111;
t114 = -pkin(10) + t118;
t108 = cos(t115);
t98 = pkin(4) * t108 + t106;
t107 = sin(t115);
t99 = pkin(4) * t107 + t146;
t131 = t114 * t141 + t117 * t99 + t98 * t144 + t145;
t94 = t117 * t136 + t135;
t95 = -t117 * t137 + t134;
t130 = -t94 * t114 + t99 * t143 + t95 * t98 + t133;
t129 = g(1) * t122 - t147;
t109 = qJ(5) + t115;
t104 = sin(t109);
t105 = cos(t109);
t93 = t117 * t135 + t136;
t80 = t93 * t104 + t105 * t140;
t82 = t95 * t104 - t105 * t143;
t86 = t104 * t144 - t117 * t105;
t128 = g(1) * t82 + g(2) * t80 + g(3) * t86;
t92 = -t117 * t134 + t137;
t127 = t111 + t93 * t98 - t92 * t114 + (-pkin(8) - t99) * t140;
t79 = -g(1) * t94 - g(2) * t92 + g(3) * t141;
t123 = cos(qJ(6));
t119 = sin(qJ(6));
t87 = t117 * t104 + t105 * t144;
t83 = t104 * t143 + t95 * t105;
t81 = -t104 * t140 + t93 * t105;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t126 - g(2) * t122, t129, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t144, -t79, -g(3) * t117 - t129 * t116, -g(1) * t133 - g(2) * t132 - g(3) * t145, 0, 0, 0, 0, 0, 0, -g(1) * (t116 * t138 + t95 * t124) - g(2) * (-t120 * t140 + t93 * t124) - g(3) * (t121 * t142 + t139) -g(1) * (-t95 * t120 + t122 * t142) - g(2) * (-t93 * t120 - t124 * t140) - g(3) * (t117 * t124 - t120 * t144) t79, -g(1) * (t95 * pkin(2) + t94 * pkin(9) + t133) - g(2) * (t93 * pkin(2) + t92 * pkin(9) + t132) - g(3) * ((pkin(2) * t121 - pkin(9) * t125) * t116 + t145) 0, 0, 0, 0, 0, 0, -g(1) * (t107 * t143 + t95 * t108) - g(2) * (-t107 * t140 + t93 * t108) - g(3) * (t117 * t107 + t108 * t144) -g(1) * (-t95 * t107 + t108 * t143) - g(2) * (-t93 * t107 - t108 * t140) - g(3) * (-t107 * t144 + t117 * t108) t79, -g(1) * (t95 * t106 - t94 * t118 + t133) - g(2) * (t93 * t106 - t92 * t118 + t111) - g(3) * (pkin(3) * t139 + t145) + (-g(1) * pkin(3) * t138 - g(3) * (t106 * t121 + t118 * t125) - (-pkin(8) - t146) * t147) * t116, 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t87, t128, t79, -g(1) * t130 - g(2) * t127 - g(3) * t131, 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t119 + t83 * t123) - g(2) * (t92 * t119 + t81 * t123) - g(3) * (-t119 * t141 + t87 * t123) -g(1) * (-t83 * t119 + t94 * t123) - g(2) * (-t81 * t119 + t92 * t123) - g(3) * (-t87 * t119 - t123 * t141) -t128, -g(1) * (t83 * pkin(5) + t82 * pkin(11) + t130) - g(2) * (t81 * pkin(5) + t80 * pkin(11) + t127) - g(3) * (t87 * pkin(5) + t86 * pkin(11) + t131);];
U_reg  = t1;
