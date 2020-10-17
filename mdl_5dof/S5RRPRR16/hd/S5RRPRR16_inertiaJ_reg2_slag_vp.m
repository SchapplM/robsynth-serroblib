% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR16_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:26
% EndTime: 2019-12-31 20:47:30
% DurationCPUTime: 1.10s
% Computational Cost: add. (919->136), mult. (2137->283), div. (0->0), fcn. (2236->8), ass. (0->97)
t104 = -2 * pkin(2);
t49 = cos(pkin(5));
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t48 = sin(pkin(5));
t55 = cos(qJ(2));
t88 = t48 * t55;
t19 = t49 * t51 + t54 * t88;
t103 = t19 ^ 2;
t102 = -0.2e1 * t19;
t101 = 0.2e1 * t48;
t53 = cos(qJ(5));
t100 = 0.2e1 * t53;
t99 = 2 * qJ(3);
t56 = -pkin(2) - pkin(8);
t98 = pkin(1) * t55;
t52 = sin(qJ(2));
t37 = t48 * t52;
t64 = -qJ(3) * t52 - pkin(1);
t13 = (t55 * t56 + t64) * t48;
t31 = pkin(7) * t37;
t66 = -pkin(2) - t98;
t8 = pkin(3) * t37 + t31 + (-pkin(8) + t66) * t49;
t5 = -t13 * t51 + t54 * t8;
t3 = -pkin(4) * t37 - t5;
t50 = sin(qJ(5));
t97 = t3 * t50;
t96 = t3 * t53;
t95 = t54 * pkin(4);
t21 = t49 * t54 - t51 * t88;
t9 = t21 * t50 - t37 * t53;
t94 = t9 * t53;
t11 = t21 * t53 + t37 * t50;
t93 = t11 * t50;
t92 = t19 * t51;
t91 = t21 * t54;
t41 = t48 ^ 2;
t90 = t41 * t55;
t46 = t54 ^ 2;
t89 = t46 * t56;
t87 = t49 * t52;
t86 = t49 * t55;
t85 = t50 * t19;
t84 = t50 * t51;
t83 = t50 * t53;
t82 = t50 * t54;
t81 = t51 * t56;
t80 = t53 * t19;
t79 = t53 * t51;
t38 = t53 * t54;
t78 = t53 * t56;
t77 = t54 * t11;
t76 = t54 * t19;
t75 = t54 * t51;
t74 = t54 * t56;
t23 = pkin(1) * t87 + pkin(7) * t88;
t43 = t50 ^ 2;
t45 = t53 ^ 2;
t73 = t43 + t45;
t44 = t51 ^ 2;
t30 = t44 + t46;
t72 = -0.2e1 * t75;
t71 = t48 * t87;
t70 = t51 * t37;
t69 = t56 * t37;
t68 = t48 * t86;
t67 = t50 * t38;
t40 = t49 * qJ(3);
t16 = -t40 - t23;
t65 = t73 * t51;
t12 = pkin(3) * t88 - t16;
t63 = -pkin(9) * t51 - t95;
t6 = t13 * t54 + t51 * t8;
t4 = pkin(9) * t37 + t6;
t7 = pkin(4) * t19 - pkin(9) * t21 + t12;
t1 = -t4 * t50 + t53 * t7;
t2 = t4 * t53 + t50 * t7;
t62 = -t1 * t50 + t2 * t53;
t61 = t5 * t54 + t6 * t51;
t60 = t93 - t94;
t25 = pkin(4) * t51 - pkin(9) * t54 + qJ(3);
t14 = t25 * t53 - t50 * t81;
t15 = t25 * t50 + t51 * t78;
t59 = -t14 * t50 + t15 * t53;
t57 = qJ(3) ^ 2;
t47 = t56 ^ 2;
t42 = t49 ^ 2;
t39 = t46 * t47;
t36 = t41 * t55 ^ 2;
t35 = t41 * t52 ^ 2;
t28 = t54 * t37;
t26 = 0.2e1 * t52 * t90;
t24 = t30 * t56;
t22 = pkin(1) * t86 - t31;
t18 = t49 * t66 + t31;
t17 = (-pkin(2) * t55 + t64) * t48;
t10 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t35, t26, 0.2e1 * t71, t36, 0.2e1 * t68, t42, 0.2e1 * pkin(1) * t90 + 0.2e1 * t22 * t49, -0.2e1 * pkin(1) * t41 * t52 - 0.2e1 * t23 * t49, (-t22 * t52 + t23 * t55) * t101, pkin(1) ^ 2 * t41 + t22 ^ 2 + t23 ^ 2, t42, -0.2e1 * t71, -0.2e1 * t68, t35, t26, t36, (-t16 * t55 + t18 * t52) * t101, 0.2e1 * t17 * t88 + 0.2e1 * t18 * t49, -0.2e1 * t16 * t49 - 0.2e1 * t17 * t37, t16 ^ 2 + t17 ^ 2 + t18 ^ 2, t21 ^ 2, t21 * t102, 0.2e1 * t21 * t37, t103, t37 * t102, t35, 0.2e1 * t12 * t19 + 0.2e1 * t37 * t5, 0.2e1 * t12 * t21 - 0.2e1 * t37 * t6, -0.2e1 * t19 * t6 - 0.2e1 * t21 * t5, t12 ^ 2 + t5 ^ 2 + t6 ^ 2, t11 ^ 2, -0.2e1 * t11 * t9, 0.2e1 * t11 * t19, t9 ^ 2, t9 * t102, t103, 0.2e1 * t1 * t19 + 0.2e1 * t3 * t9, 0.2e1 * t11 * t3 - 0.2e1 * t19 * t2, -0.2e1 * t1 * t11 - 0.2e1 * t2 * t9, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t88, t49, t22, -t23, 0, 0, t49, -t37, -t88, 0, 0, 0, (-pkin(2) * t52 + qJ(3) * t55) * t48, t31 + (t104 - t98) * t49, 0.2e1 * t40 + t23, -pkin(2) * t18 - qJ(3) * t16, t91, -t21 * t51 - t76, t28, t92, -t70, 0, qJ(3) * t19 + t12 * t51 + t54 * t69, qJ(3) * t21 + t12 * t54 - t51 * t69, (-t21 * t56 - t5) * t54 + (-t19 * t56 - t6) * t51, t12 * qJ(3) + t56 * t61, t53 * t77, (-t93 - t94) * t54, t11 * t51 + t53 * t76, t9 * t82, -t50 * t76 - t51 * t9, t92, t1 * t51 + t14 * t19 + (-t56 * t9 + t97) * t54, -t15 * t19 - t2 * t51 + (-t11 * t56 + t96) * t54, -t14 * t11 - t15 * t9 + (-t1 * t53 - t2 * t50) * t54, t1 * t14 + t15 * t2 - t3 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, t104, t99, pkin(2) ^ 2 + t57, t46, t72, 0, t44, 0, 0, t51 * t99, t54 * t99, -0.2e1 * t24, t44 * t47 + t39 + t57, t45 * t46, -0.2e1 * t46 * t83, t75 * t100, t43 * t46, t50 * t72, t44, 0.2e1 * t14 * t51 - 0.2e1 * t50 * t89, -0.2e1 * t15 * t51 - 0.2e1 * t46 * t78, 0.2e1 * (-t14 * t53 - t15 * t50) * t54, t14 ^ 2 + t15 ^ 2 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t49, 0, t18, 0, 0, 0, 0, 0, 0, t28, -t70, -t91 - t92, t61, 0, 0, 0, 0, 0, 0, -t19 * t84 - t54 * t9, -t19 * t79 - t77, t60 * t51, -t3 * t54 + t51 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t30, t24, 0, 0, 0, 0, 0, 0, -t30 * t50, -t30 * t53, 0, t51 * t59 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t73 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, t37, t5, -t6, 0, 0, t93, t11 * t53 - t50 * t9, t85, -t94, t80, 0, -pkin(4) * t9 - pkin(9) * t85 - t96, -pkin(4) * t11 - pkin(9) * t80 + t97, pkin(9) * t60 + t62, -t3 * pkin(4) + pkin(9) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, -t51, 0, t74, -t81, 0, 0, t67, (-t43 + t45) * t54, t84, -t67, t79, 0, t50 * t63 + t53 * t74, -t50 * t74 + t53 * t63, t59, pkin(4) * t74 + pkin(9) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t51, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t82, t65, pkin(9) * t65 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t43, 0.2e1 * t83, 0, t45, 0, 0, pkin(4) * t100, -0.2e1 * pkin(4) * t50, 0.2e1 * t73 * pkin(9), pkin(9) ^ 2 * t73 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, t19, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, -t82, t51, t14, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t79, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t53, 0, -t50 * pkin(9), -t53 * pkin(9), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t10;
