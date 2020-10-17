% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:08
% EndTime: 2019-12-31 22:02:12
% DurationCPUTime: 0.96s
% Computational Cost: add. (795->132), mult. (1686->233), div. (0->0), fcn. (1778->6), ass. (0->83)
t56 = sin(qJ(4));
t57 = sin(qJ(3));
t59 = cos(qJ(4));
t60 = cos(qJ(3));
t33 = t56 * t57 - t59 * t60;
t46 = -t60 * pkin(3) - pkin(2);
t23 = t33 * pkin(4) + t46;
t90 = 0.2e1 * t23;
t89 = 0.2e1 * t46;
t58 = sin(qJ(2));
t88 = -0.2e1 * t58;
t61 = cos(qJ(2));
t87 = -0.2e1 * t61;
t86 = 0.2e1 * t61;
t85 = -pkin(8) - pkin(7);
t84 = pkin(2) * t60;
t83 = pkin(6) * t57;
t53 = t58 ^ 2;
t82 = t53 * pkin(6);
t81 = t56 * pkin(3);
t50 = t58 * pkin(6);
t51 = t59 * pkin(3);
t80 = t61 * pkin(3);
t79 = t61 * pkin(4);
t78 = t33 * t61;
t35 = t56 * t60 + t59 * t57;
t77 = t35 * t61;
t76 = t57 * t58;
t75 = t57 * t60;
t74 = t57 * t61;
t73 = t60 * t58;
t72 = t60 * t61;
t37 = pkin(3) * t76 + t50;
t52 = t57 ^ 2;
t54 = t60 ^ 2;
t71 = t52 + t54;
t70 = t58 * t86;
t69 = pkin(6) * t72;
t68 = t57 * t73;
t38 = -t61 * pkin(2) - t58 * pkin(7) - pkin(1);
t30 = t60 * t38;
t10 = -pkin(8) * t73 + t30 + (-pkin(3) - t83) * t61;
t14 = t69 + (-pkin(8) * t58 + t38) * t57;
t3 = t59 * t10 - t56 * t14;
t39 = t85 * t57;
t40 = t85 * t60;
t16 = t59 * t39 + t56 * t40;
t4 = t56 * t10 + t59 * t14;
t21 = -pkin(6) * t74 + t30;
t22 = t57 * t38 + t69;
t67 = -t21 * t57 + t22 * t60;
t17 = t56 * t39 - t59 * t40;
t28 = -t56 * t76 + t59 * t73;
t66 = -t28 * qJ(5) + t3;
t26 = t35 * t58;
t2 = -t26 * qJ(5) + t4;
t65 = pkin(3) ^ 2;
t64 = pkin(6) ^ 2;
t62 = 0.2e1 * pkin(4);
t55 = t61 ^ 2;
t49 = t53 * t64;
t48 = t56 ^ 2 * t65;
t47 = -0.2e1 * t81;
t45 = t51 + pkin(4);
t42 = t56 * t80;
t32 = t35 ^ 2;
t31 = t33 ^ 2;
t29 = t33 * t81;
t25 = t28 ^ 2;
t24 = t26 ^ 2;
t20 = t28 * t87;
t19 = t26 * t87;
t18 = t26 * t81;
t15 = -0.2e1 * t35 * t33;
t13 = t28 * t35;
t12 = t26 * t33;
t11 = t26 * pkin(4) + t37;
t8 = -0.2e1 * t28 * t26;
t7 = -t33 * qJ(5) + t17;
t6 = -t35 * qJ(5) + t16;
t5 = -t35 * t26 - t28 * t33;
t1 = t66 - t79;
t9 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t53, t70, 0, t55, 0, 0, pkin(1) * t86, pkin(1) * t88, 0.2e1 * (t53 + t55) * pkin(6), pkin(1) ^ 2 + t55 * t64 + t49, t54 * t53, -0.2e1 * t53 * t75, t72 * t88, t52 * t53, t57 * t70, t55, -0.2e1 * t21 * t61 + 0.2e1 * t57 * t82, 0.2e1 * t22 * t61 + 0.2e1 * t60 * t82, 0.2e1 * (-t21 * t60 - t22 * t57) * t58, t21 ^ 2 + t22 ^ 2 + t49, t25, t8, t20, t24, -t19, t55, 0.2e1 * t37 * t26 - 0.2e1 * t3 * t61, 0.2e1 * t37 * t28 + 0.2e1 * t4 * t61, -0.2e1 * t4 * t26 - 0.2e1 * t3 * t28, t3 ^ 2 + t37 ^ 2 + t4 ^ 2, t25, t8, t20, t24, -t19, t55, -0.2e1 * t1 * t61 + 0.2e1 * t11 * t26, 0.2e1 * t11 * t28 + 0.2e1 * t2 * t61, -0.2e1 * t1 * t28 - 0.2e1 * t2 * t26, t1 ^ 2 + t11 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t61, 0, -t50, -t61 * pkin(6), 0, 0, t68, (-t52 + t54) * t58, -t74, -t68, -t72, 0, -pkin(6) * t73 + (-pkin(2) * t58 + pkin(7) * t61) * t57, pkin(7) * t72 + (t83 - t84) * t58, t67, -pkin(2) * t50 + pkin(7) * t67, t13, t5, -t77, t12, t78, 0, -t16 * t61 + t46 * t26 + t37 * t33, t17 * t61 + t46 * t28 + t37 * t35, -t16 * t28 - t17 * t26 - t3 * t35 - t4 * t33, t3 * t16 + t4 * t17 + t37 * t46, t13, t5, -t77, t12, t78, 0, t11 * t33 + t23 * t26 - t6 * t61, t11 * t35 + t23 * t28 + t7 * t61, -t1 * t35 - t2 * t33 - t7 * t26 - t6 * t28, t1 * t6 + t11 * t23 + t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t52, 0.2e1 * t75, 0, t54, 0, 0, 0.2e1 * t84, -0.2e1 * pkin(2) * t57, 0.2e1 * t71 * pkin(7), pkin(7) ^ 2 * t71 + pkin(2) ^ 2, t32, t15, 0, t31, 0, 0, t33 * t89, t35 * t89, -0.2e1 * t16 * t35 - 0.2e1 * t17 * t33, t16 ^ 2 + t17 ^ 2 + t46 ^ 2, t32, t15, 0, t31, 0, 0, t33 * t90, t35 * t90, -0.2e1 * t7 * t33 - 0.2e1 * t6 * t35, t23 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, -t76, -t61, t21, -t22, 0, 0, 0, 0, t28, 0, -t26, -t61, -t59 * t80 + t3, -t4 + t42, -t28 * t51 - t18, (t3 * t59 + t4 * t56) * pkin(3), 0, 0, t28, 0, -t26, -t61, (-pkin(4) - t45) * t61 + t66, -t2 + t42, -t45 * t28 - t18, t1 * t45 + t2 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, t60, 0, -t57 * pkin(7), -t60 * pkin(7), 0, 0, 0, 0, t35, 0, -t33, 0, t16, -t17, -t35 * t51 - t29, (t16 * t59 + t17 * t56) * pkin(3), 0, 0, t35, 0, -t33, 0, t6, -t7, -t45 * t35 - t29, t6 * t45 + t7 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t51, t47, 0, t59 ^ 2 * t65 + t48, 0, 0, 0, 0, 0, 1, 0.2e1 * t45, t47, 0, t45 ^ 2 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, -t61, t3, -t4, 0, 0, 0, 0, t28, 0, -t26, -t61, t66 - 0.2e1 * t79, -t2, -t28 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, -t33, 0, t16, -t17, 0, 0, 0, 0, t35, 0, -t33, 0, t6, -t7, -t35 * pkin(4), t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t51, -t81, 0, 0, 0, 0, 0, 0, 0, 1, t62 + t51, -t81, 0, t45 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t62, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t28, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t35, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t9;
