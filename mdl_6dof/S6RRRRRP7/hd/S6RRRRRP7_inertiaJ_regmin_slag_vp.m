% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:31:30
% EndTime: 2019-05-08 05:31:35
% DurationCPUTime: 1.40s
% Computational Cost: add. (2052->179), mult. (4628->332), div. (0->0), fcn. (5378->10), ass. (0->115)
t115 = cos(qJ(4));
t79 = sin(pkin(6));
t84 = sin(qJ(2));
t112 = t79 * t84;
t80 = cos(pkin(6));
t83 = sin(qJ(3));
t86 = cos(qJ(3));
t48 = t83 * t112 - t80 * t86;
t49 = t86 * t112 + t80 * t83;
t82 = sin(qJ(4));
t32 = t115 * t48 + t82 * t49;
t127 = -0.2e1 * t32;
t74 = -t86 * pkin(3) - pkin(2);
t126 = 0.2e1 * t74;
t125 = 0.2e1 * t86;
t124 = -pkin(10) - pkin(9);
t123 = pkin(1) * t84;
t87 = cos(qJ(2));
t122 = pkin(1) * t87;
t111 = t79 * t87;
t97 = pkin(8) * t111;
t43 = t97 + (pkin(9) + t123) * t80;
t44 = (-pkin(2) * t87 - pkin(9) * t84 - pkin(1)) * t79;
t27 = -t83 * t43 + t86 * t44;
t98 = pkin(3) * t111;
t19 = -t49 * pkin(10) + t27 - t98;
t28 = t86 * t43 + t83 * t44;
t23 = -t48 * pkin(10) + t28;
t102 = -t115 * t19 + t82 * t23;
t8 = pkin(4) * t111 + t102;
t85 = cos(qJ(5));
t121 = t8 * t85;
t81 = sin(qJ(5));
t120 = t81 * pkin(5);
t119 = t82 * pkin(3);
t118 = t85 * pkin(5);
t117 = t85 * pkin(11);
t94 = t115 * pkin(3);
t72 = -t94 - pkin(4);
t116 = pkin(4) - t72;
t33 = t115 * t49 - t82 * t48;
t26 = -t81 * t111 + t85 * t33;
t24 = t26 * t81;
t64 = t124 * t86;
t91 = t115 * t83;
t39 = -t124 * t91 - t82 * t64;
t114 = t39 * t85;
t76 = t79 ^ 2;
t113 = t76 * t87;
t110 = t80 * t84;
t29 = t81 * t32;
t59 = t82 * t86 + t91;
t109 = t81 * t59;
t71 = pkin(11) + t119;
t108 = t81 * t71;
t107 = t81 * t85;
t106 = t82 * t83;
t30 = t85 * t32;
t40 = t124 * t106 - t115 * t64;
t105 = t85 * t40;
t104 = t85 * t59;
t103 = t85 * t71;
t75 = t85 * qJ(6);
t58 = -t115 * t86 + t106;
t101 = -0.2e1 * t59 * t58;
t100 = -0.2e1 * t111;
t99 = 0.2e1 * t111;
t96 = t83 * t111;
t95 = t86 * t111;
t66 = pkin(8) * t112;
t42 = t66 + (-pkin(2) - t122) * t80;
t35 = t48 * pkin(3) + t42;
t13 = t32 * pkin(4) - t33 * pkin(11) + t35;
t92 = t115 * t23;
t11 = t82 * t19 + t92;
t9 = -pkin(11) * t111 + t11;
t4 = t85 * t13 - t81 * t9;
t1 = t32 * pkin(5) - t26 * qJ(6) + t4;
t25 = t85 * t111 + t81 * t33;
t5 = t81 * t13 + t85 * t9;
t3 = -t25 * qJ(6) + t5;
t93 = -t1 * t81 + t3 * t85;
t36 = t58 * pkin(4) - t59 * pkin(11) + t74;
t21 = t85 * t36 - t81 * t40;
t15 = t58 * pkin(5) - t59 * t75 + t21;
t17 = t105 + (-qJ(6) * t59 + t36) * t81;
t90 = -t15 * t81 + t17 * t85;
t89 = -pkin(4) * t59 - pkin(11) * t58;
t88 = -t58 * t71 + t59 * t72;
t78 = t85 ^ 2;
t77 = t81 ^ 2;
t73 = -pkin(4) - t118;
t69 = t76 * t87 ^ 2;
t68 = 0.2e1 * t107;
t63 = t75 + t117;
t62 = (-qJ(6) - pkin(11)) * t81;
t61 = t72 - t118;
t57 = t63 * t85;
t56 = t59 ^ 2;
t55 = t75 + t103;
t54 = (-qJ(6) - t71) * t81;
t53 = pkin(1) * t110 + t97;
t52 = t80 * t122 - t66;
t51 = t85 * t58;
t50 = t81 * t58;
t47 = t55 * t85;
t46 = t81 * t104;
t38 = t39 * t81;
t37 = (-t77 + t78) * t59;
t31 = pkin(5) * t109 + t39;
t22 = t81 * t36 + t105;
t14 = -t81 * t25 + t26 * t85;
t7 = t8 * t81;
t6 = t25 * pkin(5) + t8;
t2 = [1, 0, 0, t76 * t84 ^ 2, 0.2e1 * t84 * t113, 0.2e1 * t79 * t110, t80 * t99, t80 ^ 2, 0.2e1 * pkin(1) * t113 + 0.2e1 * t52 * t80, -0.2e1 * t76 * t123 - 0.2e1 * t53 * t80, t49 ^ 2, -0.2e1 * t49 * t48, t49 * t100, t48 * t99, t69, -0.2e1 * t27 * t111 + 0.2e1 * t42 * t48, 0.2e1 * t28 * t111 + 0.2e1 * t42 * t49, t33 ^ 2, t33 * t127, t33 * t100, t32 * t99, t69, 0.2e1 * t102 * t111 + 0.2e1 * t35 * t32, 0.2e1 * t11 * t111 + 0.2e1 * t35 * t33, t26 ^ 2, -0.2e1 * t26 * t25, 0.2e1 * t26 * t32, t25 * t127, t32 ^ 2, 0.2e1 * t8 * t25 + 0.2e1 * t4 * t32, 0.2e1 * t8 * t26 - 0.2e1 * t5 * t32, -0.2e1 * t1 * t26 - 0.2e1 * t3 * t25, t1 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t112, t111, t80, t52, -t53, t49 * t83, -t83 * t48 + t49 * t86, -t96, -t95, 0, -pkin(2) * t48 + pkin(9) * t96 - t42 * t86, -pkin(2) * t49 + pkin(9) * t95 + t42 * t83, t33 * t59, -t59 * t32 - t33 * t58, -t59 * t111, t58 * t111, 0, t39 * t111 + t74 * t32 + t35 * t58, t40 * t111 + t74 * t33 + t35 * t59, t26 * t104 (-t25 * t85 - t24) * t59, t32 * t104 + t26 * t58, -t32 * t109 - t25 * t58, t32 * t58, t8 * t109 + t21 * t32 + t39 * t25 + t4 * t58, t8 * t104 - t22 * t32 + t39 * t26 - t5 * t58, -t15 * t26 - t17 * t25 + (-t1 * t85 - t3 * t81) * t59, t1 * t15 + t3 * t17 + t6 * t31; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t83 ^ 2, t83 * t125, 0, 0, 0, pkin(2) * t125, -0.2e1 * pkin(2) * t83, t56, t101, 0, 0, 0, t58 * t126, t59 * t126, t78 * t56, -0.2e1 * t56 * t107, 0.2e1 * t58 * t104, t81 * t101, t58 ^ 2, 0.2e1 * t39 * t109 + 0.2e1 * t21 * t58, 0.2e1 * t39 * t104 - 0.2e1 * t22 * t58, 0.2e1 * (-t15 * t85 - t17 * t81) * t59, t15 ^ 2 + t17 ^ 2 + t31 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t48, -t111, t27, -t28, 0, 0, t33, -t32, -t111, -t94 * t111 - t102, -t92 + (-t19 + t98) * t82, t24, t14, t29, t30, 0, -t32 * t108 + t72 * t25 - t121, -t32 * t103 + t72 * t26 + t7, -t55 * t25 - t54 * t26 + t93, t1 * t54 + t3 * t55 + t6 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t86, 0, -t83 * pkin(9), -t86 * pkin(9), 0, 0, t59, -t58, 0, -t39, -t40, t46, t37, t50, t51, 0, t88 * t81 - t114, t88 * t85 + t38 (-t54 * t85 - t55 * t81) * t59 + t90, t15 * t54 + t17 * t55 + t31 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t94, -0.2e1 * t119, t77, t68, 0, 0, 0, -0.2e1 * t72 * t85, 0.2e1 * t72 * t81, -0.2e1 * t54 * t81 + 0.2e1 * t47, t54 ^ 2 + t55 ^ 2 + t61 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, -t111, -t102, -t11, t24, t14, t29, t30, 0, -pkin(4) * t25 - pkin(11) * t29 - t121, -pkin(4) * t26 - pkin(11) * t30 + t7, -t63 * t25 - t62 * t26 + t93, t1 * t62 + t3 * t63 + t6 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t58, 0, -t39, -t40, t46, t37, t50, t51, 0, t89 * t81 - t114, t89 * t85 + t38 (-t62 * t85 - t63 * t81) * t59 + t90, t15 * t62 + t17 * t63 + t31 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t94, -t119, t77, t68, 0, 0, 0, t116 * t85, -t116 * t81, t47 + t57 + (-t54 - t62) * t81, t54 * t62 + t55 * t63 + t61 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t77, t68, 0, 0, 0, 0.2e1 * pkin(4) * t85, -0.2e1 * pkin(4) * t81, -0.2e1 * t62 * t81 + 0.2e1 * t57, t62 ^ 2 + t63 ^ 2 + t73 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t32, t4, -t5, -pkin(5) * t26, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, -t109, t58, t21, -t22, -pkin(5) * t104, t15 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t85, 0, -t108, -t103, -t120, t54 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t85, 0, -t81 * pkin(11), -t117, -t120, t62 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t2;
