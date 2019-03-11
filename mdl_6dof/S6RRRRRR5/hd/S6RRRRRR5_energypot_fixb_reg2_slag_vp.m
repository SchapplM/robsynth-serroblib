% Calculate inertial parameters regressor of potential energy for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:04:46
% EndTime: 2019-03-10 04:04:46
% DurationCPUTime: 0.24s
% Computational Cost: add. (333->107), mult. (459->157), div. (0->0), fcn. (534->14), ass. (0->55)
t127 = -pkin(10) - pkin(9);
t126 = cos(qJ(1));
t148 = g(2) * t126;
t120 = sin(qJ(3));
t147 = t120 * pkin(3);
t118 = cos(pkin(6));
t146 = t118 * pkin(8) + pkin(7);
t124 = cos(qJ(3));
t107 = t124 * pkin(3) + pkin(2);
t117 = sin(pkin(6));
t121 = sin(qJ(2));
t145 = t117 * t121;
t122 = sin(qJ(1));
t144 = t117 * t122;
t143 = t117 * t124;
t125 = cos(qJ(2));
t142 = t117 * t125;
t141 = t117 * t126;
t140 = t118 * t120;
t139 = t120 * t122;
t138 = t122 * t121;
t137 = t122 * t125;
t136 = t126 * t121;
t135 = t126 * t125;
t116 = qJ(3) + qJ(4);
t134 = t126 * pkin(1) + pkin(8) * t144;
t112 = t122 * pkin(1);
t133 = -pkin(8) * t141 + t112;
t108 = sin(t116);
t100 = pkin(4) * t108 + t147;
t115 = -pkin(11) + t127;
t109 = cos(t116);
t99 = pkin(4) * t109 + t107;
t132 = t118 * t100 + t115 * t142 + t99 * t145 + t146;
t95 = t118 * t137 + t136;
t96 = -t118 * t138 + t135;
t131 = t100 * t144 - t95 * t115 + t96 * t99 + t134;
t130 = g(1) * t122 - t148;
t111 = qJ(5) + t116;
t105 = sin(t111);
t106 = cos(t111);
t94 = t118 * t136 + t137;
t81 = t94 * t105 + t106 * t141;
t83 = t96 * t105 - t106 * t144;
t87 = t105 * t145 - t118 * t106;
t129 = g(1) * t83 + g(2) * t81 + g(3) * t87;
t93 = -t118 * t135 + t138;
t128 = t112 + t94 * t99 - t93 * t115 + (-pkin(8) - t100) * t141;
t80 = -g(1) * t95 - g(2) * t93 + g(3) * t142;
t123 = cos(qJ(6));
t119 = sin(qJ(6));
t88 = t118 * t105 + t106 * t145;
t84 = t105 * t144 + t96 * t106;
t82 = -t105 * t141 + t94 * t106;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t126 - g(2) * t122, t130, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t145, -t80, -g(3) * t118 - t130 * t117, -g(1) * t134 - g(2) * t133 - g(3) * t146, 0, 0, 0, 0, 0, 0, -g(1) * (t117 * t139 + t96 * t124) - g(2) * (-t120 * t141 + t94 * t124) - g(3) * (t121 * t143 + t140) -g(1) * (-t96 * t120 + t122 * t143) - g(2) * (-t94 * t120 - t124 * t141) - g(3) * (t118 * t124 - t120 * t145) t80, -g(1) * (t96 * pkin(2) + t95 * pkin(9) + t134) - g(2) * (t94 * pkin(2) + t93 * pkin(9) + t133) - g(3) * ((pkin(2) * t121 - pkin(9) * t125) * t117 + t146) 0, 0, 0, 0, 0, 0, -g(1) * (t108 * t144 + t96 * t109) - g(2) * (-t108 * t141 + t94 * t109) - g(3) * (t118 * t108 + t109 * t145) -g(1) * (-t96 * t108 + t109 * t144) - g(2) * (-t94 * t108 - t109 * t141) - g(3) * (-t108 * t145 + t118 * t109) t80, -g(1) * (t96 * t107 - t95 * t127 + t134) - g(2) * (t94 * t107 - t93 * t127 + t112) - g(3) * (pkin(3) * t140 + t146) + (-g(1) * pkin(3) * t139 - g(3) * (t107 * t121 + t125 * t127) - (-pkin(8) - t147) * t148) * t117, 0, 0, 0, 0, 0, 0, -g(1) * t84 - g(2) * t82 - g(3) * t88, t129, t80, -g(1) * t131 - g(2) * t128 - g(3) * t132, 0, 0, 0, 0, 0, 0, -g(1) * (t95 * t119 + t84 * t123) - g(2) * (t93 * t119 + t82 * t123) - g(3) * (-t119 * t142 + t88 * t123) -g(1) * (-t84 * t119 + t95 * t123) - g(2) * (-t82 * t119 + t93 * t123) - g(3) * (-t88 * t119 - t123 * t142) -t129, -g(1) * (t84 * pkin(5) + t83 * pkin(12) + t131) - g(2) * (t82 * pkin(5) + t81 * pkin(12) + t128) - g(3) * (t88 * pkin(5) + t87 * pkin(12) + t132);];
U_reg  = t1;
