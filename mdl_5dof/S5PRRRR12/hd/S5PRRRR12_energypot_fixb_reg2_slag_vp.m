% Calculate inertial parameters regressor of potential energy for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:11
% EndTime: 2024-09-28 18:08:11
% DurationCPUTime: 0.08s
% Computational Cost: add. (296->135), mult. (387->217), div. (0->0), fcn. (399->26), ass. (0->82)
t103 = pkin(7) + pkin(8);
t90 = sin(pkin(6));
t130 = pkin(10) * t90;
t91 = sin(pkin(5));
t129 = g(3) * t91;
t102 = cos(qJ(2));
t78 = t102 * pkin(2) + pkin(1);
t88 = qJ(2) + qJ(3);
t85 = qJ(4) + t88;
t75 = sin(t85);
t89 = sin(pkin(11));
t128 = t89 * t75;
t94 = cos(pkin(5));
t127 = t89 * t94;
t126 = t90 * t91;
t125 = t91 * t89;
t92 = cos(pkin(11));
t124 = t91 * t92;
t98 = sin(qJ(2));
t123 = t91 * t98;
t122 = t92 * t75;
t121 = t92 * t94;
t93 = cos(pkin(6));
t95 = sin(qJ(5));
t120 = t93 * t95;
t99 = cos(qJ(5));
t119 = t93 * t99;
t118 = t94 * t95;
t117 = t94 * t98;
t116 = t94 * t99;
t115 = t89 * t102;
t114 = t92 * t102;
t87 = pkin(9) + t103;
t113 = t89 * t130;
t112 = t92 * t130;
t111 = t95 * t126;
t110 = t99 * t126;
t109 = t93 * t118;
t108 = t93 * t116;
t80 = pkin(5) - t88;
t79 = pkin(5) + t88;
t107 = g(1) * t89 - g(2) * t92;
t101 = cos(qJ(3));
t97 = sin(qJ(3));
t106 = pkin(3) * t102 * t97 + (t101 * pkin(3) + pkin(2)) * t98;
t100 = cos(qJ(4));
t53 = t92 * pkin(4) + t94 * t113;
t55 = pkin(4) * t127 - t112;
t96 = sin(qJ(4));
t105 = t92 * pkin(3) + t53 * t100 - t55 * t96;
t54 = -t89 * pkin(4) + t94 * t112;
t56 = pkin(4) * t121 + t113;
t104 = pkin(3) * t89 - t54 * t100 + t56 * t96;
t50 = -g(3) * t94 - t107 * t91;
t84 = t92 * pkin(1);
t83 = t89 * pkin(1);
t82 = cos(t88);
t81 = sin(t88);
t76 = cos(t85);
t74 = -qJ(4) + t80;
t73 = qJ(4) + t79;
t72 = cos(t79);
t71 = sin(t80);
t70 = cos(t73);
t69 = sin(t74);
t67 = cos(t80) / 0.2e1;
t66 = sin(t79) / 0.2e1;
t65 = t93 * pkin(10) + t87;
t64 = cos(t74) / 0.2e1;
t63 = sin(t73) / 0.2e1;
t62 = pkin(3) * t82 + t78;
t61 = -t96 * pkin(4) + t100 * t130;
t60 = pkin(4) * t100 + t96 * t130 + pkin(3);
t59 = pkin(2) * t117 - t91 * t103;
t58 = t67 + t72 / 0.2e1;
t57 = t66 - t71 / 0.2e1;
t52 = t64 + t70 / 0.2e1;
t51 = t63 - t69 / 0.2e1;
t49 = t106 * t94 - t91 * t87;
t48 = pkin(3) * t121 + t56 * t100 + t54 * t96;
t47 = -pkin(3) * t127 - t55 * t100 - t53 * t96;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t92 - g(2) * t89, t107, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * (-t89 * t117 + t114) - g(2) * (t92 * t117 + t115) - g(3) * t123, -g(1) * (-t94 * t115 - t92 * t98) - g(2) * (t94 * t114 - t89 * t98) - t102 * t129, t50, -g(1) * (pkin(7) * t125 + t84) - g(2) * (-pkin(7) * t124 + t83) - g(3) * (t94 * pkin(7) + qJ(1)), 0, 0, 0, 0, 0, 0, -g(1) * (-t89 * t57 + t92 * t82) - g(2) * (t92 * t57 + t89 * t82) - g(3) * (t67 - t72 / 0.2e1), -g(1) * (-t89 * t58 - t92 * t81) - g(2) * (t92 * t58 - t89 * t81) - g(3) * (t66 + t71 / 0.2e1), t50, -g(1) * (-t89 * t59 + t92 * t78) - g(2) * (t92 * t59 + t89 * t78) - g(3) * (pkin(2) * t123 + t94 * t103 + qJ(1)), 0, 0, 0, 0, 0, 0, -g(1) * (-t89 * t51 + t92 * t76) - g(2) * (t92 * t51 + t89 * t76) - g(3) * (t64 - t70 / 0.2e1), -g(1) * (-t89 * t52 - t122) - g(2) * (t92 * t52 - t128) - g(3) * (t63 + t69 / 0.2e1), t50, -g(1) * (-t89 * t49 + t92 * t62) - g(2) * (t92 * t49 + t89 * t62) - g(3) * (t106 * t91 + t94 * t87 + qJ(1)), 0, 0, 0, 0, 0, 0, -g(1) * ((-t89 * t109 + t92 * t99) * t76 + (-t89 * t116 - t92 * t120) * t75 + t89 * t111) - g(2) * ((t92 * t109 + t89 * t99) * t76 + (t92 * t116 - t89 * t120) * t75 - t92 * t111) - g(3) * (t90 * t118 + (t76 * t120 + t75 * t99) * t91), -g(1) * ((-t89 * t108 - t92 * t95) * t76 + (t89 * t118 - t92 * t119) * t75 + t89 * t110) - g(2) * ((t92 * t108 - t89 * t95) * t76 + (-t92 * t118 - t89 * t119) * t75 - t92 * t110) - g(3) * (t90 * t116 + (t76 * t119 - t75 * t95) * t91), t50 * t93 + (-g(1) * (t76 * t127 + t122) - g(2) * (-t76 * t121 + t128) + t76 * t129) * t90, -g(1) * ((t92 * pkin(2) + t105 * t101 + t47 * t97) * t102 + (-pkin(2) * t127 + t47 * t101 - t105 * t97) * t98 + t65 * t125 + t84) - g(2) * ((t89 * pkin(2) + t104 * t101 + t48 * t97) * t102 + (pkin(2) * t121 + t48 * t101 - t104 * t97) * t98 - t65 * t124 + t83) - g(3) * (t65 * t94 + qJ(1) + ((t60 * t101 + t61 * t97 + pkin(2)) * t98 - (t61 * t101 - t97 * t60) * t102) * t91);];
U_reg = t1;
