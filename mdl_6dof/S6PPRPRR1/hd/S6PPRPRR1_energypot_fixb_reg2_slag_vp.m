% Calculate inertial parameters regressor of potential energy for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:37
% EndTime: 2019-03-08 18:43:37
% DurationCPUTime: 0.34s
% Computational Cost: add. (574->109), mult. (1452->172), div. (0->0), fcn. (1852->16), ass. (0->71)
t105 = sin(pkin(7));
t109 = cos(pkin(11));
t106 = sin(pkin(6));
t110 = cos(pkin(7));
t134 = t106 * t110;
t103 = sin(pkin(12));
t104 = sin(pkin(11));
t108 = cos(pkin(12));
t111 = cos(pkin(6));
t133 = t109 * t111;
t90 = -t104 * t103 + t108 * t133;
t75 = -t90 * t105 - t109 * t134;
t137 = t104 * t111;
t92 = -t109 * t103 - t108 * t137;
t76 = t104 * t134 - t92 * t105;
t136 = t106 * t108;
t89 = -t105 * t136 + t111 * t110;
t67 = -g(1) * t76 - g(2) * t75 - g(3) * t89;
t114 = sin(qJ(3));
t143 = pkin(3) * t114;
t142 = pkin(8) + qJ(4);
t138 = t104 * t106;
t141 = t109 * pkin(1) + qJ(2) * t138;
t140 = t111 * qJ(2) + qJ(1);
t139 = t103 * t106;
t135 = t106 * t109;
t87 = t105 * t143 + t142 * t110;
t88 = -t142 * t105 + t110 * t143;
t93 = -t103 * t137 + t109 * t108;
t117 = cos(qJ(3));
t98 = t117 * pkin(3) + pkin(2);
t132 = t87 * t138 + t92 * t88 + t93 * t98 + t141;
t131 = t111 * t87 + t88 * t136 + t98 * t139 + t140;
t130 = g(1) * t104 - g(2) * t109;
t100 = t104 * pkin(1);
t129 = -qJ(2) * t135 + t100;
t102 = sin(pkin(13));
t107 = cos(pkin(13));
t128 = t117 * t102 + t114 * t107;
t95 = -t114 * t102 + t117 * t107;
t113 = sin(qJ(5));
t116 = cos(qJ(5));
t84 = t128 * t105;
t86 = t128 * t110;
t91 = t103 * t133 + t104 * t108;
t62 = -t84 * t135 + t90 * t86 + t91 * t95;
t55 = t62 * t113 - t75 * t116;
t64 = t84 * t138 + t92 * t86 + t93 * t95;
t57 = t64 * t113 - t76 * t116;
t70 = t111 * t84 + (t103 * t95 + t108 * t86) * t106;
t65 = t70 * t113 - t89 * t116;
t127 = g(1) * t57 + g(2) * t55 + g(3) * t65;
t125 = t95 * t105;
t121 = t106 * t125;
t85 = t95 * t110;
t61 = -t109 * t121 - t128 * t91 + t85 * t90;
t63 = t104 * t121 - t128 * t93 + t92 * t85;
t69 = (-t103 * t128 + t108 * t85) * t106 + t111 * t125;
t126 = g(1) * t63 + g(2) * t61 + g(3) * t69;
t124 = t70 * pkin(4) - pkin(9) * t69 + t131;
t123 = t64 * pkin(4) - t63 * pkin(9) + t132;
t122 = g(1) * t93 + g(2) * t91 + g(3) * t139;
t120 = t100 + t90 * t88 + t91 * t98 + (-qJ(2) - t87) * t135;
t119 = t62 * pkin(4) - t61 * pkin(9) + t120;
t118 = -g(1) * (t105 * t138 + t110 * t92) - g(2) * (-t105 * t135 + t110 * t90) - g(3) * (t105 * t111 + t108 * t134);
t115 = cos(qJ(6));
t112 = sin(qJ(6));
t66 = t89 * t113 + t70 * t116;
t58 = t76 * t113 + t64 * t116;
t56 = t75 * t113 + t62 * t116;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t109 - g(2) * t104, t130, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t122, -g(1) * t92 - g(2) * t90 - g(3) * t136, -g(3) * t111 - t130 * t106, -g(1) * t141 - g(2) * t129 - g(3) * t140, 0, 0, 0, 0, 0, 0, t118 * t114 - t122 * t117, t122 * t114 + t118 * t117, t67, -g(1) * (t93 * pkin(2) + t141) - g(2) * (t91 * pkin(2) + t129) - g(3) * (pkin(2) * t139 + t140) + t67 * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t64 - g(2) * t62 - g(3) * t70, -t126, t67, -g(1) * t132 - g(2) * t120 - g(3) * t131, 0, 0, 0, 0, 0, 0, -g(1) * t58 - g(2) * t56 - g(3) * t66, t127, t126, -g(1) * t123 - g(2) * t119 - g(3) * t124, 0, 0, 0, 0, 0, 0, -g(1) * (-t63 * t112 + t58 * t115) - g(2) * (-t61 * t112 + t56 * t115) - g(3) * (-t69 * t112 + t66 * t115) -g(1) * (-t58 * t112 - t63 * t115) - g(2) * (-t56 * t112 - t61 * t115) - g(3) * (-t66 * t112 - t69 * t115) -t127, -g(1) * (t58 * pkin(5) + t57 * pkin(10) + t123) - g(2) * (t56 * pkin(5) + t55 * pkin(10) + t119) - g(3) * (pkin(5) * t66 + pkin(10) * t65 + t124);];
U_reg  = t1;
