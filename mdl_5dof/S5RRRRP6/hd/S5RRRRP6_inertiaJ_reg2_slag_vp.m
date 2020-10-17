% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRP6
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:34
% EndTime: 2019-12-31 21:54:38
% DurationCPUTime: 1.04s
% Computational Cost: add. (757->116), mult. (1503->209), div. (0->0), fcn. (1626->6), ass. (0->80)
t55 = sin(qJ(3));
t84 = t55 * pkin(2);
t44 = pkin(8) + t84;
t54 = sin(qJ(4));
t50 = t54 ^ 2;
t57 = cos(qJ(4));
t52 = t57 ^ 2;
t72 = t50 + t52;
t74 = t72 * t44;
t59 = cos(qJ(2));
t87 = -pkin(7) - pkin(6);
t38 = t87 * t59;
t58 = cos(qJ(3));
t56 = sin(qJ(2));
t67 = t87 * t56;
t15 = -t38 * t55 - t58 * t67;
t93 = t15 ^ 2;
t31 = t55 * t56 - t58 * t59;
t28 = t31 ^ 2;
t33 = t55 * t59 + t56 * t58;
t92 = 0.2e1 * t33;
t47 = -pkin(2) * t59 - pkin(1);
t91 = 0.2e1 * t47;
t90 = 0.2e1 * t54;
t89 = -0.2e1 * t57;
t88 = 0.2e1 * t59;
t86 = t31 * pkin(4);
t85 = t54 * pkin(4);
t83 = t58 * pkin(2);
t79 = t54 * t33;
t8 = pkin(4) * t79 + t15;
t82 = t8 * t57;
t45 = -pkin(3) - t83;
t81 = pkin(3) - t45;
t80 = t15 * t57;
t78 = t54 * t57;
t17 = -t58 * t38 + t55 * t67;
t77 = t57 * t17;
t25 = t57 * t33;
t76 = -qJ(5) - pkin(8);
t46 = -pkin(4) * t57 - pkin(3);
t35 = t46 - t83;
t75 = t35 + t46;
t73 = t72 * pkin(8);
t51 = t56 ^ 2;
t53 = t59 ^ 2;
t71 = t51 + t53;
t70 = qJ(5) * t33;
t69 = qJ(5) + t44;
t68 = -0.2e1 * t33 * t31;
t10 = pkin(3) * t31 - pkin(8) * t33 + t47;
t5 = t57 * t10 - t17 * t54;
t63 = -t57 * t70 + t5;
t2 = t63 + t86;
t4 = t77 + (t10 - t70) * t54;
t66 = -t2 * t54 + t4 * t57;
t65 = -pkin(3) * t33 - pkin(8) * t31;
t6 = t10 * t54 + t77;
t1 = -t5 * t54 + t57 * t6;
t64 = -t31 * t44 + t33 * t45;
t41 = 0.2e1 * t78;
t37 = t76 * t57;
t36 = t76 * t54;
t30 = t37 * t57;
t29 = t33 ^ 2;
t27 = t69 * t57;
t26 = t69 * t54;
t24 = t57 * t31;
t23 = t52 * t29;
t22 = t54 * t31;
t21 = t50 * t29;
t20 = t27 * t57;
t19 = t54 * t25;
t18 = -0.2e1 * t29 * t78;
t14 = t15 * t54;
t13 = 0.2e1 * t31 * t25;
t12 = t54 * t68;
t11 = (-t50 + t52) * t33;
t7 = t8 * t54;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t51, t56 * t88, 0, t53, 0, 0, pkin(1) * t88, -0.2e1 * pkin(1) * t56, 0.2e1 * t71 * pkin(6), pkin(6) ^ 2 * t71 + pkin(1) ^ 2, t29, t68, 0, t28, 0, 0, t31 * t91, t33 * t91, 0.2e1 * t15 * t33 - 0.2e1 * t17 * t31, t17 ^ 2 + t47 ^ 2 + t93, t23, t18, t13, t21, t12, t28, 0.2e1 * t15 * t79 + 0.2e1 * t31 * t5, 0.2e1 * t15 * t25 - 0.2e1 * t31 * t6, (-t5 * t57 - t54 * t6) * t92, t5 ^ 2 + t6 ^ 2 + t93, t23, t18, t13, t21, t12, t28, 0.2e1 * t2 * t31 + 0.2e1 * t79 * t8, 0.2e1 * t25 * t8 - 0.2e1 * t31 * t4, (-t2 * t57 - t4 * t54) * t92, t2 ^ 2 + t4 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, t59, 0, -t56 * pkin(6), -t59 * pkin(6), 0, 0, 0, 0, t33, 0, -t31, 0, -t15, -t17, (-t31 * t55 - t33 * t58) * pkin(2), (-t15 * t58 + t17 * t55) * pkin(2), t19, t11, t22, -t19, t24, 0, t54 * t64 - t80, t57 * t64 + t14, t1, t1 * t44 + t15 * t45, t19, t11, t22, -t19, t24, 0, -t26 * t31 + t35 * t79 - t82, t25 * t35 - t27 * t31 + t7, (t26 * t57 - t27 * t54) * t33 + t66, -t2 * t26 + t27 * t4 + t35 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t83, -0.2e1 * t84, 0, (t55 ^ 2 + t58 ^ 2) * pkin(2) ^ 2, t50, t41, 0, t52, 0, 0, t45 * t89, t45 * t90, 0.2e1 * t74, t44 ^ 2 * t72 + t45 ^ 2, t50, t41, 0, t52, 0, 0, t35 * t89, t35 * t90, 0.2e1 * t26 * t54 + 0.2e1 * t20, t26 ^ 2 + t27 ^ 2 + t35 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t31, 0, -t15, -t17, 0, 0, t19, t11, t22, -t19, t24, 0, t54 * t65 - t80, t57 * t65 + t14, t1, -t15 * pkin(3) + pkin(8) * t1, t19, t11, t22, -t19, t24, 0, t31 * t36 + t46 * t79 - t82, t25 * t46 + t31 * t37 + t7, (-t36 * t57 + t37 * t54) * t33 + t66, t2 * t36 - t37 * t4 + t46 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t83, -t84, 0, 0, t50, t41, 0, t52, 0, 0, t81 * t57, -t81 * t54, t73 + t74, -t45 * pkin(3) + pkin(8) * t74, t50, t41, 0, t52, 0, 0, -t75 * t57, t75 * t54, t20 - t30 + (t26 - t36) * t54, -t26 * t36 - t27 * t37 + t35 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t50, t41, 0, t52, 0, 0, 0.2e1 * pkin(3) * t57, -0.2e1 * pkin(3) * t54, 0.2e1 * t73, pkin(8) ^ 2 * t72 + pkin(3) ^ 2, t50, t41, 0, t52, 0, 0, t46 * t89, t46 * t90, -0.2e1 * t36 * t54 - 0.2e1 * t30, t36 ^ 2 + t37 ^ 2 + t46 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t79, t31, t5, -t6, 0, 0, 0, 0, t25, 0, -t79, t31, t63 + 0.2e1 * t86, -t4, -pkin(4) * t25, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t57, 0, -t54 * t44, -t57 * t44, 0, 0, 0, 0, t54, 0, t57, 0, -t26, -t27, -t85, -t26 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t57, 0, -t54 * pkin(8), -t57 * pkin(8), 0, 0, 0, 0, t54, 0, t57, 0, t36, t37, -t85, t36 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t25, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t54, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t54, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
