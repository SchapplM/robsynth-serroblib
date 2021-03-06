% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:45:23
% EndTime: 2019-05-06 02:45:29
% DurationCPUTime: 2.12s
% Computational Cost: add. (1824->170), mult. (3379->289), div. (0->0), fcn. (3891->10), ass. (0->116)
t93 = sin(qJ(5));
t86 = t93 ^ 2;
t97 = cos(qJ(5));
t88 = t97 ^ 2;
t110 = t86 + t88;
t94 = sin(qJ(4));
t132 = t94 * pkin(3);
t80 = pkin(9) + t132;
t147 = t110 * t80;
t90 = sin(pkin(11));
t135 = t90 * pkin(1);
t76 = pkin(7) + t135;
t127 = pkin(8) + t76;
t95 = sin(qJ(3));
t107 = t127 * t95;
t99 = cos(qJ(3));
t53 = t127 * t99;
t98 = cos(qJ(4));
t31 = t98 * t107 + t53 * t94;
t146 = t31 ^ 2;
t63 = t94 * t95 - t98 * t99;
t145 = t63 ^ 2;
t144 = 0.2e1 * t63;
t91 = cos(pkin(11));
t134 = t91 * pkin(1);
t77 = -pkin(2) - t134;
t67 = -pkin(3) * t99 + t77;
t143 = 0.2e1 * t67;
t129 = t98 * pkin(3);
t82 = -pkin(5) * t97 - pkin(4);
t68 = t82 - t129;
t142 = 0.2e1 * t68;
t141 = 0.2e1 * t82;
t140 = 0.2e1 * t95;
t33 = -t94 * t107 + t98 * t53;
t118 = t97 * t33;
t136 = t63 * pkin(4);
t66 = t94 * t99 + t95 * t98;
t137 = pkin(9) * t66;
t25 = t136 + t67 - t137;
t10 = t118 + (-pkin(10) * t66 + t25) * t93;
t11 = t97 * t25 - t33 * t93;
t117 = t97 * t66;
t138 = pkin(5) * t63;
t7 = -pkin(10) * t117 + t11 + t138;
t92 = sin(qJ(6));
t96 = cos(qJ(6));
t3 = -t10 * t92 + t96 * t7;
t119 = t96 * t10;
t4 = t7 * t92 + t119;
t61 = t92 * t93 - t96 * t97;
t65 = t92 * t97 + t93 * t96;
t139 = -t3 * t65 - t4 * t61;
t133 = t92 * pkin(5);
t131 = t96 * pkin(5);
t130 = t97 * pkin(9);
t81 = -pkin(4) - t129;
t128 = pkin(4) - t81;
t28 = t65 * t66;
t126 = t28 * t65;
t125 = t31 * t63;
t124 = t31 * t97;
t123 = t63 * t61;
t122 = t66 * t86;
t121 = t93 * t66;
t120 = t93 * t97;
t116 = t97 * t80;
t54 = (-pkin(10) - t80) * t93;
t85 = t97 * pkin(10);
t55 = t85 + t116;
t34 = t54 * t96 - t55 * t92;
t35 = t54 * t92 + t55 * t96;
t115 = -t34 * t65 - t35 * t61;
t69 = (-pkin(10) - pkin(9)) * t93;
t70 = t85 + t130;
t43 = t69 * t96 - t70 * t92;
t44 = t69 * t92 + t70 * t96;
t114 = -t43 * t65 - t44 * t61;
t113 = t68 + t82;
t111 = pkin(9) * t110;
t87 = t95 ^ 2;
t89 = t99 ^ 2;
t109 = t87 + t89;
t108 = -0.2e1 * t66 * t63;
t105 = -pkin(4) * t66 - pkin(9) * t63;
t12 = t25 * t93 + t118;
t5 = -t11 * t93 + t12 * t97;
t104 = -t63 * t80 + t66 * t81;
t74 = 0.2e1 * t120;
t60 = t66 ^ 2;
t59 = t65 ^ 2;
t56 = t61 ^ 2;
t52 = t97 * t63;
t51 = t88 * t66;
t50 = t88 * t60;
t49 = t93 * t63;
t48 = t86 * t60;
t46 = t93 * t117;
t42 = t65 * t63;
t39 = -0.2e1 * t65 * t61;
t38 = t51 + t122;
t37 = t51 - t122;
t36 = (-t61 * t92 - t65 * t96) * pkin(5);
t30 = t96 * t117 - t92 * t121;
t27 = t30 ^ 2;
t26 = t28 ^ 2;
t22 = t31 * t93;
t18 = t30 * t65;
t17 = t30 * t61;
t16 = t28 * t61;
t15 = pkin(5) * t121 + t31;
t14 = t15 * t65;
t13 = t15 * t61;
t9 = -t17 - t126;
t8 = -t17 + t126;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t134, -0.2e1 * t135, 0 (t90 ^ 2 + t91 ^ 2) * pkin(1) ^ 2, t87, t99 * t140, 0, t89, 0, 0, -0.2e1 * t77 * t99, t77 * t140, 0.2e1 * t109 * t76, t109 * t76 ^ 2 + t77 ^ 2, t60, t108, 0, t145, 0, 0, t63 * t143, t66 * t143, 0.2e1 * t31 * t66 - 0.2e1 * t33 * t63, t33 ^ 2 + t67 ^ 2 + t146, t50, -0.2e1 * t60 * t120, t117 * t144, t48, t93 * t108, t145, 0.2e1 * t11 * t63 + 0.2e1 * t31 * t121, 0.2e1 * t31 * t117 - 0.2e1 * t12 * t63, 0.2e1 * (-t11 * t97 - t12 * t93) * t66, t11 ^ 2 + t12 ^ 2 + t146, t27, -0.2e1 * t30 * t28, t30 * t144, t26, -t28 * t144, t145, 0.2e1 * t15 * t28 + 0.2e1 * t3 * t63, 0.2e1 * t15 * t30 - 0.2e1 * t4 * t63, -0.2e1 * t28 * t4 - 0.2e1 * t3 * t30, t15 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 * t66 + t125, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t66 + t125, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t63 - t28 * t3 + t30 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 + t145, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 + t48 + t145, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 + t26 + t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, t99, 0, -t95 * t76, -t99 * t76, 0, 0, 0, 0, t66, 0, -t63, 0, -t31, -t33 (-t63 * t94 - t66 * t98) * pkin(3) (-t31 * t98 + t33 * t94) * pkin(3), t46, t37, t49, -t46, t52, 0, t104 * t93 - t124, t104 * t97 + t22, t5, t31 * t81 + t5 * t80, t18, t9, t42, t16, -t123, 0, t28 * t68 + t34 * t63 + t13, t30 * t68 - t35 * t63 + t14, -t28 * t35 - t30 * t34 + t139, t15 * t68 + t3 * t34 + t35 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, -t95, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t66, 0 (-t63 * t98 + t66 * t94) * pkin(3), 0, 0, 0, 0, 0, 0, -t52, t49, t38, t147 * t66 + t63 * t81, 0, 0, 0, 0, 0, 0, t123, t42, t8, -t28 * t34 + t30 * t35 + t63 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t129, -0.2e1 * t132, 0 (t94 ^ 2 + t98 ^ 2) * pkin(3) ^ 2, t86, t74, 0, t88, 0, 0, -0.2e1 * t81 * t97, 0.2e1 * t81 * t93, 0.2e1 * t147, t110 * t80 ^ 2 + t81 ^ 2, t59, t39, 0, t56, 0, 0, t61 * t142, t65 * t142, 0.2e1 * t115, t34 ^ 2 + t35 ^ 2 + t68 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, -t63, 0, -t31, -t33, 0, 0, t46, t37, t49, -t46, t52, 0, t105 * t93 - t124, t105 * t97 + t22, t5, -t31 * pkin(4) + t5 * pkin(9), t18, t9, t42, t16, -t123, 0, t28 * t82 + t43 * t63 + t13, t30 * t82 - t44 * t63 + t14, -t28 * t44 - t30 * t43 + t139, t15 * t82 + t3 * t43 + t4 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t66, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t49, t38, t110 * t137 - t136, 0, 0, 0, 0, 0, 0, t123, t42, t8, -t28 * t43 + t30 * t44 + t63 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t129, -t132, 0, 0, t86, t74, 0, t88, 0, 0, t128 * t97, -t128 * t93, t111 + t147, -t81 * pkin(4) + pkin(9) * t147, t59, t39, 0, t56, 0, 0, t113 * t61, t113 * t65, t114 + t115, t34 * t43 + t35 * t44 + t68 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t86, t74, 0, t88, 0, 0, 0.2e1 * pkin(4) * t97, -0.2e1 * pkin(4) * t93, 0.2e1 * t111, t110 * pkin(9) ^ 2 + pkin(4) ^ 2, t59, t39, 0, t56, 0, 0, t61 * t141, t65 * t141, 0.2e1 * t114, t43 ^ 2 + t44 ^ 2 + t82 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, -t121, t63, t11, -t12, 0, 0, 0, 0, t30, 0, -t28, t63, t63 * t131 + t3, -t119 + (-t7 - t138) * t92 (-t28 * t92 - t30 * t96) * pkin(5) (t3 * t96 + t4 * t92) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t117, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t30, 0 (-t28 * t96 + t30 * t92) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, t97, 0, -t93 * t80, -t116, 0, 0, 0, 0, t65, 0, -t61, 0, t34, -t35, t36 (t34 * t96 + t35 * t92) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, t97, 0, -t93 * pkin(9), -t130, 0, 0, 0, 0, t65, 0, -t61, 0, t43, -t44, t36 (t43 * t96 + t44 * t92) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t131, -0.2e1 * t133, 0 (t92 ^ 2 + t96 ^ 2) * pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, -t28, t63, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, -t61, 0, t34, -t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, -t61, 0, t43, -t44, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t131, -t133, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t1;
