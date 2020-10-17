% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:22
% EndTime: 2019-12-31 22:12:27
% DurationCPUTime: 1.37s
% Computational Cost: add. (1276->180), mult. (3103->338), div. (0->0), fcn. (3368->8), ass. (0->113)
t69 = sin(pkin(5));
t73 = sin(qJ(2));
t101 = t69 * t73;
t70 = cos(pkin(5));
t72 = sin(qJ(3));
t75 = cos(qJ(3));
t38 = t72 * t101 - t70 * t75;
t37 = t38 ^ 2;
t76 = cos(qJ(2));
t100 = t69 * t76;
t40 = t75 * t101 + t70 * t72;
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t25 = t74 * t100 + t40 * t71;
t120 = -0.2e1 * t25;
t119 = -0.2e1 * t40;
t118 = 0.2e1 * t69;
t117 = -0.2e1 * t72;
t116 = 0.2e1 * t72;
t115 = 0.2e1 * t75;
t114 = pkin(1) * t73;
t113 = pkin(1) * t76;
t112 = pkin(3) * t74;
t111 = pkin(8) * t71;
t110 = t38 * pkin(4);
t66 = t72 ^ 2;
t109 = t66 * pkin(8);
t108 = t71 * pkin(4);
t107 = t72 * pkin(8);
t87 = pkin(7) * t100;
t33 = t87 + (pkin(8) + t114) * t70;
t34 = (-pkin(2) * t76 - pkin(8) * t73 - pkin(1)) * t69;
t15 = -t72 * t33 + t75 * t34;
t13 = pkin(3) * t100 - t15;
t106 = t13 * t71;
t105 = t13 * t74;
t21 = t25 * t74;
t27 = -t71 * t100 + t40 * t74;
t22 = t27 * t71;
t104 = t38 * t75;
t103 = t40 * t72;
t64 = t69 ^ 2;
t102 = t64 * t76;
t99 = t70 * t73;
t35 = t71 * t38;
t98 = t71 * t72;
t97 = t71 * t74;
t96 = t71 * t75;
t95 = t72 * t38;
t36 = t74 * t38;
t61 = t74 * t72;
t94 = t74 * t75;
t93 = -qJ(5) - pkin(9);
t65 = t71 ^ 2;
t67 = t74 ^ 2;
t92 = t65 + t67;
t91 = qJ(5) * t72;
t90 = 0.2e1 * t100;
t89 = t72 * t115;
t88 = pkin(8) * t94;
t86 = t25 * t98;
t85 = t72 * t100;
t84 = t75 * t100;
t54 = pkin(7) * t101;
t32 = t54 + (-pkin(2) - t113) * t70;
t12 = t38 * pkin(3) - t40 * pkin(9) + t32;
t16 = t75 * t33 + t72 * t34;
t14 = -pkin(9) * t100 + t16;
t3 = t74 * t12 - t71 * t14;
t46 = -t75 * pkin(3) - t72 * pkin(9) - pkin(2);
t43 = t74 * t46;
t83 = -t74 * t91 + t43;
t4 = t71 * t12 + t74 * t14;
t82 = -t3 * t71 + t4 * t74;
t81 = -t15 * t72 + t16 * t75;
t30 = -pkin(8) * t96 + t43;
t31 = t71 * t46 + t88;
t80 = -t30 * t71 + t31 * t74;
t79 = -t27 * qJ(5) + t3;
t78 = pkin(8) ^ 2;
t68 = t75 ^ 2;
t63 = t66 * t78;
t62 = -t74 * pkin(4) - pkin(3);
t60 = t67 * t66;
t59 = t65 * t66;
t57 = t64 * t76 ^ 2;
t56 = 0.2e1 * t97;
t53 = t71 * t61;
t51 = t94 * t117;
t50 = -0.2e1 * t66 * t97;
t49 = t71 * t89;
t48 = t93 * t74;
t47 = t93 * t71;
t45 = (pkin(8) + t108) * t72;
t44 = (-t65 + t67) * t72;
t42 = pkin(1) * t99 + t87;
t41 = t70 * t113 - t54;
t28 = t88 + (t46 - t91) * t71;
t24 = t27 ^ 2;
t23 = t25 ^ 2;
t20 = (-pkin(4) - t111) * t75 + t83;
t19 = t27 * t61;
t18 = 0.2e1 * t27 * t38;
t17 = t38 * t120;
t11 = -t27 * t75 + t38 * t61;
t10 = t25 * t75 - t71 * t95;
t8 = t27 * t120;
t7 = -t71 * t25 + t27 * t74;
t6 = (-t22 - t21) * t72;
t5 = t25 * pkin(4) + t13;
t2 = -t25 * qJ(5) + t4;
t1 = t79 + t110;
t9 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t64 * t73 ^ 2, 0.2e1 * t73 * t102, t99 * t118, t57, t70 * t90, t70 ^ 2, 0.2e1 * pkin(1) * t102 + 0.2e1 * t41 * t70, -0.2e1 * t64 * t114 - 0.2e1 * t42 * t70, (-t41 * t73 + t42 * t76) * t118, t64 * pkin(1) ^ 2 + t41 ^ 2 + t42 ^ 2, t40 ^ 2, t38 * t119, t100 * t119, t37, t38 * t90, t57, -0.2e1 * t15 * t100 + 0.2e1 * t32 * t38, 0.2e1 * t16 * t100 + 0.2e1 * t32 * t40, -0.2e1 * t15 * t40 - 0.2e1 * t16 * t38, t15 ^ 2 + t16 ^ 2 + t32 ^ 2, t24, t8, t18, t23, t17, t37, 0.2e1 * t13 * t25 + 0.2e1 * t3 * t38, 0.2e1 * t13 * t27 - 0.2e1 * t4 * t38, -0.2e1 * t4 * t25 - 0.2e1 * t3 * t27, t13 ^ 2 + t3 ^ 2 + t4 ^ 2, t24, t8, t18, t23, t17, t37, 0.2e1 * t1 * t38 + 0.2e1 * t5 * t25, -0.2e1 * t2 * t38 + 0.2e1 * t5 * t27, -0.2e1 * t1 * t27 - 0.2e1 * t2 * t25, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, t100, t70, t41, -t42, 0, 0, t103, t40 * t75 - t95, -t85, -t104, -t84, 0, -pkin(2) * t38 + pkin(8) * t85 - t32 * t75, -pkin(2) * t40 + pkin(8) * t84 + t32 * t72, (t103 - t104) * pkin(8) + t81, -t32 * pkin(2) + t81 * pkin(8), t19, t6, t11, t86, t10, -t104, -t3 * t75 + t30 * t38 + (pkin(8) * t25 + t106) * t72, -t31 * t38 + t4 * t75 + (pkin(8) * t27 + t105) * t72, -t31 * t25 - t30 * t27 + (-t3 * t74 - t4 * t71) * t72, t107 * t13 + t3 * t30 + t4 * t31, t19, t6, t11, t86, t10, -t104, -t1 * t75 + t20 * t38 + t45 * t25 + t5 * t98, t2 * t75 + t45 * t27 - t28 * t38 + t5 * t61, -t20 * t27 - t28 * t25 + (-t1 * t74 - t2 * t71) * t72, t1 * t20 + t2 * t28 + t5 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t66, t89, 0, t68, 0, 0, pkin(2) * t115, pkin(2) * t117, 0.2e1 * (t66 + t68) * pkin(8), pkin(2) ^ 2 + t68 * t78 + t63, t60, t50, t51, t59, t49, t68, 0.2e1 * t71 * t109 - 0.2e1 * t30 * t75, 0.2e1 * t74 * t109 + 0.2e1 * t31 * t75, (-t30 * t74 - t31 * t71) * t116, t30 ^ 2 + t31 ^ 2 + t63, t60, t50, t51, t59, t49, t68, -0.2e1 * t20 * t75 + 0.2e1 * t45 * t98, 0.2e1 * t28 * t75 + 0.2e1 * t45 * t61, (-t20 * t74 - t28 * t71) * t116, t20 ^ 2 + t28 ^ 2 + t45 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t38, -t100, t15, -t16, 0, 0, t22, t7, t35, -t21, t36, 0, -pkin(3) * t25 - pkin(9) * t35 - t105, -pkin(3) * t27 - pkin(9) * t36 + t106, (t22 - t21) * pkin(9) + t82, -t13 * pkin(3) + pkin(9) * t82, t22, t7, t35, -t21, t36, 0, t62 * t25 + t47 * t38 - t5 * t74, t62 * t27 + t48 * t38 + t5 * t71, -t1 * t71 + t2 * t74 + t48 * t25 - t47 * t27, t1 * t47 - t2 * t48 + t5 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, t75, 0, -t107, -t75 * pkin(8), 0, 0, t53, t44, -t96, -t53, -t94, 0, -pkin(8) * t61 + (-pkin(3) * t72 + pkin(9) * t75) * t71, pkin(9) * t94 + (t111 - t112) * t72, t80, -pkin(3) * t107 + pkin(9) * t80, t53, t44, -t96, -t53, -t94, 0, -t45 * t74 - t47 * t75 + t62 * t98, t45 * t71 - t48 * t75 + t61 * t62, (-t47 * t72 + t28) * t74 + (t48 * t72 - t20) * t71, t20 * t47 - t28 * t48 + t45 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t65, t56, 0, t67, 0, 0, 0.2e1 * t112, -0.2e1 * pkin(3) * t71, 0.2e1 * t92 * pkin(9), pkin(9) ^ 2 * t92 + pkin(3) ^ 2, t65, t56, 0, t67, 0, 0, -0.2e1 * t62 * t74, 0.2e1 * t62 * t71, -0.2e1 * t47 * t71 - 0.2e1 * t48 * t74, t47 ^ 2 + t48 ^ 2 + t62 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t25, t38, t3, -t4, 0, 0, 0, 0, t27, 0, -t25, t38, t79 + 0.2e1 * t110, -t2, -t27 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, -t98, -t75, t30, -t31, 0, 0, 0, 0, t61, 0, -t98, -t75, (-0.2e1 * pkin(4) - t111) * t75 + t83, -t28, -pkin(4) * t61, t20 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, t74, 0, -t71 * pkin(9), -t74 * pkin(9), 0, 0, 0, 0, t71, 0, t74, 0, t47, t48, -t108, t47 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t27, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t61, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t71, 0, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t9;
