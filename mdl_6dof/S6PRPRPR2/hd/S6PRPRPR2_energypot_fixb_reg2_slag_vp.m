% Calculate inertial parameters regressor of potential energy for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:32:52
% EndTime: 2019-03-08 19:32:52
% DurationCPUTime: 0.29s
% Computational Cost: add. (376->102), mult. (841->157), div. (0->0), fcn. (1059->14), ass. (0->58)
t100 = sin(pkin(6));
t131 = pkin(7) * t100;
t102 = cos(pkin(10));
t103 = cos(pkin(6));
t106 = sin(qJ(2));
t124 = t103 * t106;
t81 = pkin(2) * t124 + (-pkin(7) - qJ(3)) * t100;
t108 = cos(qJ(2));
t91 = pkin(2) * t108 + pkin(1);
t99 = sin(pkin(10));
t130 = t102 * t81 + t99 * t91;
t129 = t103 * pkin(7) + qJ(1);
t128 = cos(pkin(11));
t105 = sin(qJ(4));
t127 = t100 * t105;
t126 = t100 * t106;
t107 = cos(qJ(4));
t125 = t100 * t107;
t123 = t103 * t108;
t98 = sin(pkin(11));
t82 = -t106 * t128 - t108 * t98;
t80 = t82 * t103;
t119 = t108 * t128;
t83 = -t106 * t98 + t119;
t69 = -t102 * t80 + t83 * t99;
t122 = t69 * pkin(3) + t130;
t97 = sin(pkin(12));
t121 = pkin(5) * t97 + pkin(8);
t120 = t102 * t91 - t81 * t99;
t118 = pkin(2) * t126 + t103 * qJ(3) + t129;
t71 = t102 * t83 + t80 * t99;
t117 = t71 * pkin(3) + t120;
t79 = t82 * t100;
t116 = -t79 * pkin(3) + t118;
t115 = g(1) * t99 - g(2) * t102;
t109 = t103 * t83;
t68 = t102 * t109 + t99 * t82;
t114 = -pkin(8) * t68 + t122;
t70 = t102 * t82 - t99 * t109;
t113 = -pkin(8) * t70 + t117;
t62 = t102 * t125 + t105 * t69;
t64 = t105 * t71 - t99 * t125;
t72 = -t103 * t107 - t105 * t79;
t112 = g(1) * t64 + g(2) * t62 + g(3) * t72;
t78 = -t100 * t119 + t98 * t126;
t111 = g(1) * t70 + g(2) * t68 - g(3) * t78;
t110 = pkin(8) * t78 + t116;
t104 = -pkin(9) - qJ(5);
t101 = cos(pkin(12));
t96 = pkin(12) + qJ(6);
t93 = cos(t96);
t92 = sin(t96);
t90 = pkin(5) * t101 + pkin(4);
t77 = -g(3) * t103 - t115 * t100;
t73 = t103 * t105 - t107 * t79;
t65 = t107 * t71 + t99 * t127;
t63 = -t102 * t127 + t107 * t69;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t102 - g(2) * t99, t115, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * (t102 * t108 - t99 * t124) - g(2) * (t102 * t124 + t108 * t99) - g(3) * t126, -g(1) * (-t102 * t106 - t99 * t123) - g(2) * (t102 * t123 - t99 * t106) - g(3) * t100 * t108, t77, -g(1) * (pkin(1) * t102 + t99 * t131) - g(2) * (pkin(1) * t99 - t102 * t131) - g(3) * t129, 0, 0, 0, 0, 0, 0, -g(1) * t71 - g(2) * t69 + g(3) * t79, -t111, t77, -g(1) * t120 - g(2) * t130 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * t65 - g(2) * t63 - g(3) * t73, t112, t111, -g(1) * t113 - g(2) * t114 - g(3) * t110, 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t65 - t70 * t97) - g(2) * (t101 * t63 - t68 * t97) - g(3) * (t101 * t73 + t78 * t97) -g(1) * (-t101 * t70 - t65 * t97) - g(2) * (-t101 * t68 - t63 * t97) - g(3) * (t101 * t78 - t73 * t97) -t112, -g(1) * (pkin(4) * t65 + qJ(5) * t64 + t113) - g(2) * (pkin(4) * t63 + qJ(5) * t62 + t114) - g(3) * (pkin(4) * t73 + qJ(5) * t72 + t110) 0, 0, 0, 0, 0, 0, -g(1) * (t65 * t93 - t70 * t92) - g(2) * (t63 * t93 - t68 * t92) - g(3) * (t73 * t93 + t78 * t92) -g(1) * (-t65 * t92 - t70 * t93) - g(2) * (-t63 * t92 - t68 * t93) - g(3) * (-t73 * t92 + t78 * t93) -t112, -g(1) * (-t104 * t64 - t121 * t70 + t65 * t90 + t117) - g(2) * (-t104 * t62 - t121 * t68 + t63 * t90 + t122) - g(3) * (-t72 * t104 + t121 * t78 + t73 * t90 + t116);];
U_reg  = t1;
