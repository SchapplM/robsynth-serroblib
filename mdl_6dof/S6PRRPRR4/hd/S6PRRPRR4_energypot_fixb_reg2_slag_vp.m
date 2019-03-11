% Calculate inertial parameters regressor of potential energy for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:27
% EndTime: 2019-03-08 22:14:27
% DurationCPUTime: 0.23s
% Computational Cost: add. (310->87), mult. (725->126), div. (0->0), fcn. (900->12), ass. (0->56)
t142 = pkin(8) - pkin(9);
t106 = sin(pkin(11));
t108 = cos(pkin(11));
t112 = sin(qJ(2));
t115 = cos(qJ(2));
t136 = cos(pkin(6));
t128 = t115 * t136;
t91 = t106 * t112 - t108 * t128;
t141 = t91 * pkin(8);
t93 = t106 * t128 + t108 * t112;
t140 = t93 * pkin(8);
t139 = cos(qJ(3));
t107 = sin(pkin(6));
t138 = pkin(7) * t107;
t137 = t108 * pkin(1) + t106 * t138;
t111 = sin(qJ(3));
t135 = t107 * t111;
t134 = t107 * t112;
t133 = t107 * t115;
t132 = t136 * pkin(7) + qJ(1);
t129 = t112 * t136;
t94 = -t106 * t129 + t108 * t115;
t131 = t94 * pkin(2) + t137;
t130 = t107 * t139;
t127 = t106 * pkin(1) - t108 * t138;
t126 = g(1) * t106 - g(2) * t108;
t92 = t106 * t115 + t108 * t129;
t125 = t92 * pkin(2) + t127;
t84 = -t106 * t130 + t94 * t111;
t85 = t106 * t135 + t94 * t139;
t124 = t85 * pkin(3) + t84 * qJ(4) + t131;
t123 = pkin(2) * t134 - pkin(8) * t133 + t132;
t110 = sin(qJ(5));
t114 = cos(qJ(5));
t82 = t108 * t130 + t92 * t111;
t83 = -t108 * t135 + t92 * t139;
t68 = t83 * t110 - t82 * t114;
t70 = t85 * t110 - t84 * t114;
t95 = t111 * t134 - t136 * t139;
t96 = t136 * t111 + t112 * t130;
t74 = t96 * t110 - t95 * t114;
t122 = g(1) * t70 + g(2) * t68 + g(3) * t74;
t121 = g(1) * t84 + g(2) * t82 + g(3) * t95;
t72 = -g(1) * t93 - g(2) * t91 + g(3) * t133;
t120 = t83 * pkin(3) + t82 * qJ(4) + t125;
t119 = t96 * pkin(3) + t95 * qJ(4) + t123;
t118 = t85 * pkin(4) + t142 * t93 + t124;
t117 = t96 * pkin(4) + pkin(9) * t133 + t119;
t116 = t83 * pkin(4) + t142 * t91 + t120;
t113 = cos(qJ(6));
t109 = sin(qJ(6));
t75 = t95 * t110 + t96 * t114;
t71 = t84 * t110 + t85 * t114;
t69 = t82 * t110 + t83 * t114;
t67 = -g(1) * t85 - g(2) * t83 - g(3) * t96;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t108 - g(2) * t106, t126, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t94 - g(2) * t92 - g(3) * t134, -t72, -g(3) * t136 - t126 * t107, -g(1) * t137 - g(2) * t127 - g(3) * t132, 0, 0, 0, 0, 0, 0, t67, t121, t72, -g(1) * (t131 + t140) - g(2) * (t125 + t141) - g(3) * t123, 0, 0, 0, 0, 0, 0, t67, t72, -t121, -g(1) * (t124 + t140) - g(2) * (t120 + t141) - g(3) * t119, 0, 0, 0, 0, 0, 0, -g(1) * t71 - g(2) * t69 - g(3) * t75, t122, -t72, -g(1) * t118 - g(2) * t116 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (-t93 * t109 + t71 * t113) - g(2) * (-t91 * t109 + t69 * t113) - g(3) * (t109 * t133 + t75 * t113) -g(1) * (-t71 * t109 - t93 * t113) - g(2) * (-t69 * t109 - t91 * t113) - g(3) * (-t75 * t109 + t113 * t133) -t122, -g(1) * (t71 * pkin(5) + t70 * pkin(10) + t118) - g(2) * (t69 * pkin(5) + t68 * pkin(10) + t116) - g(3) * (t75 * pkin(5) + t74 * pkin(10) + t117);];
U_reg  = t1;
