% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:11:16
% EndTime: 2019-05-07 20:11:21
% DurationCPUTime: 1.04s
% Computational Cost: add. (1659->136), mult. (3222->239), div. (0->0), fcn. (3893->10), ass. (0->113)
t101 = sin(qJ(4));
t105 = cos(qJ(4));
t98 = sin(pkin(11));
t99 = cos(pkin(11));
t73 = -t98 * t101 + t99 * t105;
t131 = t73 * pkin(5);
t125 = t105 * pkin(4);
t127 = cos(qJ(3));
t109 = t127 * pkin(2);
t92 = -t109 - pkin(3);
t83 = t92 - t125;
t57 = t83 - t131;
t138 = 0.2e1 * t57;
t93 = -pkin(3) - t125;
t58 = t93 - t131;
t137 = 0.2e1 * t58;
t102 = sin(qJ(3));
t103 = sin(qJ(2));
t128 = cos(qJ(2));
t80 = t102 * t103 - t127 * t128;
t136 = -0.2e1 * t80;
t135 = 0.2e1 * t80;
t94 = -t128 * pkin(2) - pkin(1);
t134 = 0.2e1 * t94;
t81 = t102 * t128 + t127 * t103;
t48 = t80 * pkin(3) - t81 * pkin(9) + t94;
t86 = (-pkin(8) - pkin(7)) * t103;
t110 = t128 * pkin(7);
t87 = t128 * pkin(8) + t110;
t56 = t102 * t86 + t127 * t87;
t27 = -t101 * t56 + t105 * t48;
t95 = t105 * qJ(5);
t18 = t80 * pkin(4) - t81 * t95 + t27;
t119 = t105 * t56;
t24 = t119 + (-qJ(5) * t81 + t48) * t101;
t10 = t98 * t18 + t99 * t24;
t74 = t99 * t101 + t98 * t105;
t9 = t99 * t18 - t98 * t24;
t133 = t10 * t73 - t9 * t74;
t132 = pkin(4) * t98;
t130 = t74 * pkin(10);
t129 = pkin(3) - t92;
t126 = t102 * pkin(2);
t124 = t105 * pkin(9);
t91 = pkin(9) + t126;
t71 = (-qJ(5) - t91) * t101;
t117 = t105 * t91;
t72 = t95 + t117;
t44 = t99 * t71 - t98 * t72;
t45 = t98 * t71 + t99 * t72;
t123 = -t44 * t74 + t45 * t73;
t84 = (-qJ(5) - pkin(9)) * t101;
t85 = t95 + t124;
t53 = t99 * t84 - t98 * t85;
t54 = t98 * t84 + t99 * t85;
t122 = -t53 * t74 + t54 * t73;
t121 = t57 + t58;
t120 = t101 * t81;
t118 = t105 * t81;
t55 = t102 * t87 - t127 * t86;
t52 = t55 * t101;
t116 = t55 * t105;
t113 = t101 * t105;
t112 = 0.2e1 * t128;
t111 = t81 * t136;
t100 = sin(qJ(6));
t104 = cos(qJ(6));
t42 = t73 * t81;
t4 = t80 * pkin(5) - t42 * pkin(10) + t9;
t41 = t74 * t81;
t6 = -t41 * pkin(10) + t10;
t1 = -t100 * t6 + t104 * t4;
t35 = pkin(4) * t120 + t55;
t108 = -pkin(3) * t81 - pkin(9) * t80;
t2 = t100 * t4 + t104 * t6;
t107 = -t80 * t91 + t81 * t92;
t97 = t105 ^ 2;
t96 = t101 ^ 2;
t89 = t99 * pkin(4) + pkin(5);
t88 = 0.2e1 * t113;
t76 = t81 ^ 2;
t75 = t80 ^ 2;
t70 = t73 * pkin(10);
t69 = t105 * t80;
t68 = t101 * t80;
t65 = t100 * t89 + t104 * t132;
t64 = -t100 * t132 + t104 * t89;
t60 = t81 * t113;
t51 = (-t96 + t97) * t81;
t50 = t100 * t73 + t104 * t74;
t49 = t100 * t74 - t104 * t73;
t47 = t50 ^ 2;
t46 = (t73 * t98 - t74 * t99) * pkin(4);
t37 = t70 + t54;
t36 = t53 - t130;
t34 = t50 * t80;
t33 = t49 * t80;
t30 = t70 + t45;
t29 = t44 - t130;
t28 = t101 * t48 + t119;
t26 = -0.2e1 * t50 * t49;
t25 = t41 * pkin(5) + t35;
t23 = -t100 * t41 + t104 * t42;
t22 = t100 * t42 + t104 * t41;
t20 = t100 * t36 + t104 * t37;
t19 = -t100 * t37 + t104 * t36;
t15 = t100 * t29 + t104 * t30;
t14 = -t100 * t30 + t104 * t29;
t13 = t25 * t50;
t12 = t25 * t49;
t11 = t23 * t50;
t5 = -t50 * t22 - t23 * t49;
t3 = [1, 0, 0, t103 ^ 2, t103 * t112, 0, 0, 0, pkin(1) * t112, -0.2e1 * pkin(1) * t103, t76, t111, 0, 0, 0, t80 * t134, t81 * t134, t97 * t76, -0.2e1 * t76 * t113, t118 * t135, t101 * t111, t75, 0.2e1 * t27 * t80 + 0.2e1 * t81 * t52, 0.2e1 * t81 * t116 - 0.2e1 * t28 * t80, -0.2e1 * t10 * t41 - 0.2e1 * t9 * t42, t10 ^ 2 + t35 ^ 2 + t9 ^ 2, t23 ^ 2, -0.2e1 * t23 * t22, t23 * t135, t22 * t136, t75, 0.2e1 * t1 * t80 + 0.2e1 * t25 * t22, -0.2e1 * t2 * t80 + 0.2e1 * t25 * t23; 0, 0, 0, 0, 0, t103, t128, 0, -t103 * pkin(7), -t110, 0, 0, t81, -t80, 0, -t55, -t56, t60, t51, t68, t69, 0, t107 * t101 - t116, t107 * t105 + t52, -t45 * t41 - t44 * t42 + t133, t10 * t45 + t35 * t83 + t9 * t44, t11, t5, t34, -t33, 0, t14 * t80 + t57 * t22 + t12, -t15 * t80 + t57 * t23 + t13; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t109, -0.2e1 * t126, t96, t88, 0, 0, 0, -0.2e1 * t92 * t105, 0.2e1 * t92 * t101, 0.2e1 * t123, t44 ^ 2 + t45 ^ 2 + t83 ^ 2, t47, t26, 0, 0, 0, t49 * t138, t50 * t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t80, 0, -t55, -t56, t60, t51, t68, t69, 0, t108 * t101 - t116, t108 * t105 + t52, -t54 * t41 - t53 * t42 + t133, t10 * t54 + t35 * t93 + t9 * t53, t11, t5, t34, -t33, 0, t19 * t80 + t58 * t22 + t12, -t20 * t80 + t58 * t23 + t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t109, -t126, t96, t88, 0, 0, 0, t129 * t105, -t129 * t101, t122 + t123, t44 * t53 + t45 * t54 + t83 * t93, t47, t26, 0, 0, 0, t121 * t49, t121 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t96, t88, 0, 0, 0, 0.2e1 * pkin(3) * t105, -0.2e1 * pkin(3) * t101, 0.2e1 * t122, t53 ^ 2 + t54 ^ 2 + t93 ^ 2, t47, t26, 0, 0, 0, t49 * t137, t50 * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, -t120, t80, t27, -t28 (-t41 * t98 - t42 * t99) * pkin(4) (t10 * t98 + t9 * t99) * pkin(4), 0, 0, t23, -t22, t80, t64 * t80 + t1, -t65 * t80 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t105, 0, -t101 * t91, -t117, t46 (t44 * t99 + t45 * t98) * pkin(4), 0, 0, t50, -t49, 0, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t105, 0, -t101 * pkin(9), -t124, t46 (t53 * t99 + t54 * t98) * pkin(4), 0, 0, t50, -t49, 0, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t98 ^ 2 + t99 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t64, -0.2e1 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, t22, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, 0, 0, 0, 0, t49, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, 0, 0, 0, 0, t49, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, t80, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t49, 0, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t49, 0, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t64, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
