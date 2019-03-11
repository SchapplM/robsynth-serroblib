% Calculate inertial parameters regressor of potential energy for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:28
% EndTime: 2019-03-08 19:58:28
% DurationCPUTime: 0.26s
% Computational Cost: add. (364->91), mult. (841->138), div. (0->0), fcn. (1059->12), ass. (0->59)
t95 = sin(pkin(6));
t129 = t95 * pkin(7);
t101 = sin(qJ(2));
t97 = cos(pkin(6));
t80 = t97 * t101 * pkin(2) + (-pkin(7) - qJ(3)) * t95;
t104 = cos(qJ(2));
t90 = t104 * pkin(2) + pkin(1);
t94 = sin(pkin(10));
t96 = cos(pkin(10));
t128 = t96 * t80 + t94 * t90;
t100 = sin(qJ(4));
t127 = t100 * t95;
t126 = t101 * t95;
t103 = cos(qJ(4));
t125 = t103 * t95;
t124 = t94 * t101;
t123 = t94 * t104;
t122 = t96 * t101;
t121 = t96 * t104;
t120 = t97 * pkin(7) + qJ(1);
t119 = cos(pkin(11));
t93 = sin(pkin(11));
t81 = -t101 * t119 - t104 * t93;
t79 = t81 * t97;
t115 = t104 * t119;
t82 = -t101 * t93 + t115;
t68 = -t96 * t79 + t94 * t82;
t118 = t68 * pkin(3) + t128;
t99 = sin(qJ(5));
t117 = pkin(5) * t99 + pkin(8);
t116 = -t94 * t80 + t96 * t90;
t114 = pkin(2) * t126 + t97 * qJ(3) + t120;
t70 = t94 * t79 + t96 * t82;
t113 = t70 * pkin(3) + t116;
t112 = g(1) * t94 - g(2) * t96;
t78 = t81 * t95;
t111 = -t78 * pkin(3) + t114;
t105 = t82 * t97;
t67 = t96 * t105 + t94 * t81;
t110 = -t67 * pkin(8) + t118;
t69 = -t94 * t105 + t96 * t81;
t109 = -t69 * pkin(8) + t113;
t61 = t68 * t100 + t96 * t125;
t63 = t70 * t100 - t94 * t125;
t71 = -t78 * t100 - t97 * t103;
t108 = g(1) * t63 + g(2) * t61 + g(3) * t71;
t77 = -t95 * t115 + t93 * t126;
t107 = g(1) * t69 + g(2) * t67 - g(3) * t77;
t106 = t77 * pkin(8) + t111;
t102 = cos(qJ(5));
t98 = -qJ(6) - pkin(9);
t89 = t102 * pkin(5) + pkin(4);
t76 = -g(3) * t97 - t112 * t95;
t72 = t97 * t100 - t78 * t103;
t64 = t70 * t103 + t94 * t127;
t62 = t68 * t103 - t96 * t127;
t59 = -g(1) * (t64 * t102 - t69 * t99) - g(2) * (t62 * t102 - t67 * t99) - g(3) * (t72 * t102 + t77 * t99);
t58 = -g(1) * (-t69 * t102 - t64 * t99) - g(2) * (-t67 * t102 - t62 * t99) - g(3) * (t77 * t102 - t72 * t99);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94, t112, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * (-t97 * t124 + t121) - g(2) * (t97 * t122 + t123) - g(3) * t126, -g(1) * (-t97 * t123 - t122) - g(2) * (t97 * t121 - t124) - g(3) * t95 * t104, t76, -g(1) * (t96 * pkin(1) + t94 * t129) - g(2) * (t94 * pkin(1) - t96 * t129) - g(3) * t120, 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68 + g(3) * t78, -t107, t76, -g(1) * t116 - g(2) * t128 - g(3) * t114, 0, 0, 0, 0, 0, 0, -g(1) * t64 - g(2) * t62 - g(3) * t72, t108, t107, -g(1) * t109 - g(2) * t110 - g(3) * t106, 0, 0, 0, 0, 0, 0, t59, t58, -t108, -g(1) * (t64 * pkin(4) + t63 * pkin(9) + t109) - g(2) * (t62 * pkin(4) + t61 * pkin(9) + t110) - g(3) * (t72 * pkin(4) + t71 * pkin(9) + t106) 0, 0, 0, 0, 0, 0, t59, t58, -t108, -g(1) * (-t117 * t69 - t63 * t98 + t64 * t89 + t113) - g(2) * (-t117 * t67 - t61 * t98 + t62 * t89 + t118) - g(3) * (t117 * t77 - t71 * t98 + t72 * t89 + t111);];
U_reg  = t1;
