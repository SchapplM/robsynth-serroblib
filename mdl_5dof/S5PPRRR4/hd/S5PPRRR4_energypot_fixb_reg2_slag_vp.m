% Calculate inertial parameters regressor of potential energy for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:42
% EndTime: 2019-12-05 15:19:42
% DurationCPUTime: 0.32s
% Computational Cost: add. (336->85), mult. (878->138), div. (0->0), fcn. (1108->14), ass. (0->57)
t94 = sin(pkin(5));
t95 = cos(pkin(11));
t125 = t94 * t95;
t93 = sin(pkin(6));
t97 = cos(pkin(6));
t98 = cos(pkin(5));
t76 = -t93 * t125 + t98 * t97;
t124 = t94 * t97;
t92 = sin(pkin(10));
t128 = t92 * t98;
t91 = sin(pkin(11));
t96 = cos(pkin(10));
t79 = -t95 * t128 - t96 * t91;
t70 = t92 * t124 - t79 * t93;
t131 = cos(qJ(3));
t129 = t91 * t94;
t127 = t93 * t94;
t126 = t93 * t98;
t123 = t96 * t98;
t120 = qJ(2) * t94;
t121 = t96 * pkin(1) + t92 * t120;
t119 = t98 * qJ(2) + qJ(1);
t116 = t94 * t131;
t115 = t97 * t131;
t114 = t93 * t116;
t113 = t92 * pkin(1) - t96 * t120;
t112 = g(1) * t92 - g(2) * t96;
t77 = t95 * t123 - t92 * t91;
t69 = -t96 * t124 - t77 * t93;
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t101 = sin(qJ(3));
t78 = t91 * t123 + t92 * t95;
t59 = t78 * t131 + (-t96 * t127 + t77 * t97) * t101;
t52 = t59 * t100 - t69 * t103;
t80 = -t91 * t128 + t96 * t95;
t61 = t80 * t131 + (t92 * t127 + t79 * t97) * t101;
t54 = t61 * t100 - t70 * t103;
t68 = t91 * t116 + (t95 * t124 + t126) * t101;
t62 = t68 * t100 - t76 * t103;
t111 = g(1) * t54 + g(2) * t52 + g(3) * t62;
t58 = t78 * t101 + t96 * t114 - t77 * t115;
t60 = t80 * t101 - t92 * t114 - t79 * t115;
t67 = t101 * t129 - t115 * t125 - t131 * t126;
t110 = g(1) * t60 + g(2) * t58 + g(3) * t67;
t109 = t80 * pkin(2) + t70 * pkin(7) + t121;
t108 = pkin(2) * t129 + t76 * pkin(7) + t119;
t107 = t61 * pkin(3) + t60 * pkin(8) + t109;
t106 = t68 * pkin(3) + t67 * pkin(8) + t108;
t105 = t78 * pkin(2) + t69 * pkin(7) + t113;
t104 = t59 * pkin(3) + t58 * pkin(8) + t105;
t102 = cos(qJ(5));
t99 = sin(qJ(5));
t63 = t76 * t100 + t68 * t103;
t55 = t70 * t100 + t61 * t103;
t53 = t69 * t100 + t59 * t103;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t92, t112, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t80 - g(2) * t78 - g(3) * t129, -g(1) * t79 - g(2) * t77 - g(3) * t125, -g(3) * t98 - t112 * t94, -g(1) * t121 - g(2) * t113 - g(3) * t119, 0, 0, 0, 0, 0, 0, -g(1) * t61 - g(2) * t59 - g(3) * t68, t110, -g(1) * t70 - g(2) * t69 - g(3) * t76, -g(1) * t109 - g(2) * t105 - g(3) * t108, 0, 0, 0, 0, 0, 0, -g(1) * t55 - g(2) * t53 - g(3) * t63, t111, -t110, -g(1) * t107 - g(2) * t104 - g(3) * t106, 0, 0, 0, 0, 0, 0, -g(1) * (t55 * t102 + t60 * t99) - g(2) * (t53 * t102 + t58 * t99) - g(3) * (t63 * t102 + t67 * t99), -g(1) * (t60 * t102 - t55 * t99) - g(2) * (t58 * t102 - t53 * t99) - g(3) * (t67 * t102 - t63 * t99), -t111, -g(1) * (t55 * pkin(4) + t54 * pkin(9) + t107) - g(2) * (t53 * pkin(4) + t52 * pkin(9) + t104) - g(3) * (t63 * pkin(4) + t62 * pkin(9) + t106);];
U_reg = t1;
