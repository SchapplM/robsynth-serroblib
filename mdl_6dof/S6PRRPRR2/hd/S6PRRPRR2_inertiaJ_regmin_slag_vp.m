% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:46:37
% EndTime: 2021-01-16 03:46:41
% DurationCPUTime: 0.66s
% Computational Cost: add. (606->107), mult. (1336->193), div. (0->0), fcn. (1686->12), ass. (0->74)
t47 = sin(pkin(12));
t49 = cos(pkin(12));
t52 = sin(qJ(3));
t71 = cos(qJ(3));
t33 = t47 * t52 - t49 * t71;
t82 = -0.2e1 * t33;
t81 = 0.2e1 * t33;
t76 = t49 * pkin(3);
t42 = -pkin(4) - t76;
t54 = cos(qJ(5));
t38 = -t54 * pkin(5) + t42;
t80 = 0.2e1 * t38;
t44 = -t71 * pkin(3) - pkin(2);
t79 = 0.2e1 * t44;
t78 = t33 * pkin(5);
t77 = t47 * pkin(3);
t50 = sin(qJ(6));
t75 = t50 * pkin(5);
t53 = cos(qJ(6));
t74 = t53 * pkin(5);
t34 = t47 * t71 + t49 * t52;
t21 = t33 * pkin(4) - t34 * pkin(9) + t44;
t51 = sin(qJ(5));
t61 = t71 * pkin(8);
t39 = t71 * qJ(4) + t61;
t59 = (-pkin(8) - qJ(4)) * t52;
t24 = t49 * t39 + t47 * t59;
t65 = t54 * t24;
t7 = t65 + (-pkin(10) * t34 + t21) * t51;
t73 = t53 * t7;
t41 = pkin(9) + t77;
t72 = pkin(10) + t41;
t37 = t50 * t54 + t53 * t51;
t70 = t37 * t33;
t48 = sin(pkin(6));
t55 = cos(qJ(2));
t69 = t48 * t55;
t68 = t51 * t33;
t67 = t51 * t34;
t66 = t51 * t54;
t64 = t54 * t34;
t63 = cos(pkin(6));
t62 = 0.2e1 * t71;
t8 = t54 * t21 - t51 * t24;
t6 = -pkin(10) * t64 + t78 + t8;
t1 = -t50 * t7 + t53 * t6;
t60 = t48 * sin(qJ(2));
t22 = t47 * t39 - t49 * t59;
t58 = -t33 * t41 + t34 * t42;
t36 = t50 * t51 - t53 * t54;
t57 = -t52 * t60 + t63 * t71;
t46 = t54 ^ 2;
t45 = t51 ^ 2;
t32 = t34 ^ 2;
t31 = t33 ^ 2;
t30 = t72 * t54;
t29 = t72 * t51;
t28 = t63 * t52 + t71 * t60;
t27 = t54 * t33;
t25 = t36 * t33;
t20 = -t50 * t29 + t53 * t30;
t19 = -t53 * t29 - t50 * t30;
t17 = t36 * t34;
t16 = t37 * t34;
t15 = t49 * t28 + t47 * t57;
t13 = t47 * t28 - t49 * t57;
t12 = pkin(5) * t67 + t22;
t11 = t54 * t15 - t51 * t69;
t10 = -t51 * t15 - t54 * t69;
t9 = t51 * t21 + t65;
t4 = t50 * t10 + t53 * t11;
t3 = t53 * t10 - t50 * t11;
t2 = t50 * t6 + t73;
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 ^ 2 * t55 ^ 2 + t13 ^ 2 + t15 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t69, -t60, 0, 0, 0, 0, 0, t71 * t69, -t52 * t69, -t33 * t69, -t34 * t69, t13 * t34 - t15 * t33, t13 * t22 + t15 * t24 - t44 * t69, 0, 0, 0, 0, 0, t10 * t33 + t13 * t67, -t11 * t33 + t13 * t64, 0, 0, 0, 0, 0, t13 * t16 + t3 * t33, -t13 * t17 - t4 * t33; 0, 1, 0, 0, t52 ^ 2, t52 * t62, 0, 0, 0, pkin(2) * t62, -0.2e1 * pkin(2) * t52, t33 * t79, t34 * t79, 0.2e1 * t22 * t34 - 0.2e1 * t24 * t33, t22 ^ 2 + t24 ^ 2 + t44 ^ 2, t46 * t32, -0.2e1 * t32 * t66, t64 * t81, t67 * t82, t31, 0.2e1 * t22 * t67 + 0.2e1 * t8 * t33, 0.2e1 * t22 * t64 - 0.2e1 * t9 * t33, t17 ^ 2, 0.2e1 * t17 * t16, -t17 * t81, t16 * t82, t31, 0.2e1 * t1 * t33 + 0.2e1 * t12 * t16, -0.2e1 * t12 * t17 - 0.2e1 * t2 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t28, -t13, -t15, 0, (-t13 * t49 + t15 * t47) * pkin(3), 0, 0, 0, 0, 0, -t13 * t54, t13 * t51, 0, 0, 0, 0, 0, t13 * t36, t13 * t37; 0, 0, 0, 0, 0, 0, t52, t71, 0, -t52 * pkin(8), -t61, -t22, -t24, (-t33 * t47 - t34 * t49) * pkin(3), (-t22 * t49 + t24 * t47) * pkin(3), t51 * t64, (-t45 + t46) * t34, t68, t27, 0, -t22 * t54 + t58 * t51, t22 * t51 + t58 * t54, -t17 * t37, -t37 * t16 + t17 * t36, t70, -t25, 0, t12 * t36 + t38 * t16 + t19 * t33, t12 * t37 - t38 * t17 - t20 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t76, -0.2e1 * t77, 0, (t47 ^ 2 + t49 ^ 2) * pkin(3) ^ 2, t45, 0.2e1 * t66, 0, 0, 0, -0.2e1 * t42 * t54, 0.2e1 * t42 * t51, t37 ^ 2, -0.2e1 * t37 * t36, 0, 0, 0, t36 * t80, t37 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t34, 0, t44, 0, 0, 0, 0, 0, t27, -t68, 0, 0, 0, 0, 0, -t25, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t67, t33, t8, -t9, 0, 0, -t17, -t16, t33, t33 * t74 + t1, -t73 + (-t6 - t78) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t54, 0, -t51 * t41, -t54 * t41, 0, 0, t37, -t36, 0, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t51, 0, 0, 0, 0, 0, -t36, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t74, -0.2e1 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t16, t33, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, 0, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t74, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
