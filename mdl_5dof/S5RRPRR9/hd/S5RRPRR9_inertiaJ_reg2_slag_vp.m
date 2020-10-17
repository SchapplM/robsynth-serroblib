% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:58
% EndTime: 2019-12-31 20:22:01
% DurationCPUTime: 0.86s
% Computational Cost: add. (1019->107), mult. (1980->217), div. (0->0), fcn. (2286->8), ass. (0->74)
t48 = sin(pkin(9));
t49 = cos(pkin(9));
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t29 = t48 * t55 + t49 * t52;
t87 = -0.2e1 * t29;
t67 = -qJ(3) - pkin(6);
t37 = t67 * t55;
t62 = t67 * t52;
t17 = -t48 * t37 - t49 * t62;
t86 = t17 ^ 2;
t27 = t48 * t52 - t49 * t55;
t25 = t27 ^ 2;
t85 = 0.2e1 * t27;
t79 = t49 * pkin(2);
t41 = -pkin(3) - t79;
t54 = cos(qJ(4));
t36 = -t54 * pkin(4) + t41;
t84 = 0.2e1 * t36;
t43 = -t55 * pkin(2) - pkin(1);
t83 = 0.2e1 * t43;
t82 = 0.2e1 * t55;
t81 = t27 * pkin(4);
t80 = t48 * pkin(2);
t50 = sin(qJ(5));
t78 = t50 * pkin(4);
t53 = cos(qJ(5));
t77 = t53 * pkin(4);
t16 = t27 * pkin(3) - t29 * pkin(7) + t43;
t51 = sin(qJ(4));
t19 = -t49 * t37 + t48 * t62;
t69 = t54 * t19;
t5 = t69 + (-pkin(8) * t29 + t16) * t51;
t76 = t53 * t5;
t40 = pkin(7) + t80;
t75 = pkin(8) + t40;
t68 = t54 * t29;
t71 = t51 * t29;
t12 = -t50 * t71 + t53 * t68;
t33 = t50 * t51 - t53 * t54;
t74 = t12 * t33;
t35 = t50 * t54 + t53 * t51;
t73 = t35 * t27;
t72 = t51 * t27;
t70 = t51 * t54;
t44 = t51 ^ 2;
t46 = t54 ^ 2;
t66 = t44 + t46;
t45 = t52 ^ 2;
t47 = t55 ^ 2;
t65 = t45 + t47;
t64 = t27 * t87;
t63 = t51 * t68;
t6 = t54 * t16 - t51 * t19;
t4 = -pkin(8) * t68 + t6 + t81;
t1 = t53 * t4 - t50 * t5;
t7 = t51 * t16 + t69;
t61 = t7 * t51 + t6 * t54;
t60 = -t6 * t51 + t7 * t54;
t59 = -t27 * t40 + t29 * t41;
t31 = t35 ^ 2;
t30 = t33 ^ 2;
t26 = t29 ^ 2;
t24 = t75 * t54;
t23 = t75 * t51;
t22 = t54 * t27;
t20 = t33 * t27;
t15 = -t50 * t23 + t53 * t24;
t14 = -t53 * t23 - t50 * t24;
t10 = t35 * t29;
t9 = pkin(4) * t71 + t17;
t8 = t35 * t10;
t2 = t50 * t4 + t76;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t45, t52 * t82, 0, t47, 0, 0, pkin(1) * t82, -0.2e1 * pkin(1) * t52, 0.2e1 * t65 * pkin(6), t65 * pkin(6) ^ 2 + pkin(1) ^ 2, t26, t64, 0, t25, 0, 0, t27 * t83, t29 * t83, 0.2e1 * t17 * t29 - 0.2e1 * t19 * t27, t19 ^ 2 + t43 ^ 2 + t86, t46 * t26, -0.2e1 * t26 * t70, t68 * t85, t44 * t26, t51 * t64, t25, 0.2e1 * t17 * t71 + 0.2e1 * t6 * t27, 0.2e1 * t17 * t68 - 0.2e1 * t7 * t27, t61 * t87, t6 ^ 2 + t7 ^ 2 + t86, t12 ^ 2, -0.2e1 * t12 * t10, t12 * t85, t10 ^ 2, -t10 * t85, t25, 0.2e1 * t1 * t27 + 0.2e1 * t9 * t10, 0.2e1 * t9 * t12 - 0.2e1 * t2 * t27, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t10, t1 ^ 2 + t2 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t55, 0, -t52 * pkin(6), -t55 * pkin(6), 0, 0, 0, 0, t29, 0, -t27, 0, -t17, -t19, (-t27 * t48 - t29 * t49) * pkin(2), (-t17 * t49 + t19 * t48) * pkin(2), t63, (-t44 + t46) * t29, t72, -t63, t22, 0, -t17 * t54 + t59 * t51, t17 * t51 + t59 * t54, t60, t17 * t41 + t60 * t40, t12 * t35, -t8 - t74, t73, t10 * t33, -t20, 0, t36 * t10 + t14 * t27 + t9 * t33, t36 * t12 - t15 * t27 + t9 * t35, -t1 * t35 - t15 * t10 - t14 * t12 - t2 * t33, t1 * t14 + t2 * t15 + t9 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t79, -0.2e1 * t80, 0, (t48 ^ 2 + t49 ^ 2) * pkin(2) ^ 2, t44, 0.2e1 * t70, 0, t46, 0, 0, -0.2e1 * t41 * t54, 0.2e1 * t41 * t51, 0.2e1 * t66 * t40, t66 * t40 ^ 2 + t41 ^ 2, t31, -0.2e1 * t35 * t33, 0, t30, 0, 0, t33 * t84, t35 * t84, -0.2e1 * t14 * t35 - 0.2e1 * t15 * t33, t14 ^ 2 + t15 ^ 2 + t36 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t29, 0, t43, 0, 0, 0, 0, 0, 0, t22, -t72, -t66 * t29, t61, 0, 0, 0, 0, 0, 0, -t20, -t73, -t8 + t74, -t1 * t33 + t2 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t33 + t15 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, -t71, t27, t6, -t7, 0, 0, 0, 0, t12, 0, -t10, t27, t27 * t77 + t1, -t76 + (-t4 - t81) * t50, (-t10 * t50 - t12 * t53) * pkin(4), (t1 * t53 + t2 * t50) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t54, 0, -t51 * t40, -t54 * t40, 0, 0, 0, 0, t35, 0, -t33, 0, t14, -t15, (-t33 * t50 - t35 * t53) * pkin(4), (t14 * t53 + t15 * t50) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t51, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t35, 0, (-t33 * t53 + t35 * t50) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t77, -0.2e1 * t78, 0, (t50 ^ 2 + t53 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, t27, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, -t33, 0, t14, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t77, -t78, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
