% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR14_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:38
% EndTime: 2019-12-31 20:38:40
% DurationCPUTime: 0.59s
% Computational Cost: add. (822->103), mult. (2003->232), div. (0->0), fcn. (2340->10), ass. (0->74)
t49 = sin(pkin(10));
t51 = cos(pkin(10));
t52 = cos(pkin(5));
t50 = sin(pkin(5));
t55 = sin(qJ(2));
t75 = t50 * t55;
t30 = t49 * t75 - t52 * t51;
t31 = t52 * t49 + t51 * t75;
t54 = sin(qJ(4));
t78 = cos(qJ(4));
t18 = -t54 * t30 + t78 * t31;
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t57 = cos(qJ(2));
t74 = t50 * t57;
t12 = t53 * t18 + t56 * t74;
t83 = -0.2e1 * t12;
t82 = -0.2e1 * t18;
t43 = -t51 * pkin(3) - pkin(2);
t81 = 0.2e1 * t43;
t80 = pkin(1) * t55;
t79 = pkin(1) * t57;
t13 = t56 * t18 - t53 * t74;
t77 = t13 * t53;
t45 = t50 ^ 2;
t76 = t45 * t57;
t73 = t52 * t55;
t17 = t78 * t30 + t54 * t31;
t72 = t53 * t17;
t37 = t54 * t49 - t78 * t51;
t71 = t53 * t37;
t38 = t78 * t49 + t54 * t51;
t70 = t53 * t38;
t69 = t53 * t56;
t16 = t56 * t17;
t68 = t56 * t38;
t67 = pkin(8) + qJ(3);
t63 = pkin(7) * t74;
t26 = t63 + (qJ(3) + t80) * t52;
t27 = (-pkin(2) * t57 - qJ(3) * t55 - pkin(1)) * t50;
t15 = t51 * t26 + t49 * t27;
t66 = t49 ^ 2 + t51 ^ 2;
t65 = -0.2e1 * t38 * t37;
t64 = 0.2e1 * t74;
t62 = qJ(3) * t74;
t61 = t67 * t49;
t14 = -t49 * t26 + t51 * t27;
t60 = -pkin(4) * t38 - pkin(9) * t37;
t59 = -t14 * t49 + t15 * t51;
t41 = pkin(7) * t75;
t29 = t41 + (-pkin(2) - t79) * t52;
t11 = -t30 * pkin(8) + t15;
t8 = -pkin(3) * t74 - t31 * pkin(8) + t14;
t5 = -t54 * t11 + t78 * t8;
t6 = t78 * t11 + t54 * t8;
t19 = t30 * pkin(3) + t29;
t48 = t56 ^ 2;
t47 = t53 ^ 2;
t39 = t67 * t51;
t35 = t38 ^ 2;
t34 = pkin(1) * t73 + t63;
t33 = t52 * t79 - t41;
t32 = t56 * t37;
t22 = t78 * t39 - t54 * t61;
t21 = t54 * t39 + t78 * t61;
t20 = t37 * pkin(4) - t38 * pkin(9) + t43;
t10 = t53 * t20 + t56 * t22;
t9 = t56 * t20 - t53 * t22;
t7 = t17 * pkin(4) - t18 * pkin(9) + t19;
t4 = -pkin(9) * t74 + t6;
t3 = pkin(4) * t74 - t5;
t2 = t56 * t4 + t53 * t7;
t1 = -t53 * t4 + t56 * t7;
t23 = [1, 0, 0, t45 * t55 ^ 2, 0.2e1 * t55 * t76, 0.2e1 * t50 * t73, t52 * t64, t52 ^ 2, 0.2e1 * pkin(1) * t76 + 0.2e1 * t33 * t52, -0.2e1 * t34 * t52 - 0.2e1 * t45 * t80, -0.2e1 * t14 * t74 + 0.2e1 * t29 * t30, 0.2e1 * t15 * t74 + 0.2e1 * t29 * t31, -0.2e1 * t14 * t31 - 0.2e1 * t15 * t30, t14 ^ 2 + t15 ^ 2 + t29 ^ 2, t18 ^ 2, t17 * t82, t74 * t82, t17 * t64, t45 * t57 ^ 2, 0.2e1 * t19 * t17 - 0.2e1 * t5 * t74, 0.2e1 * t19 * t18 + 0.2e1 * t6 * t74, t13 ^ 2, t13 * t83, 0.2e1 * t13 * t17, t17 * t83, t17 ^ 2, 0.2e1 * t1 * t17 + 0.2e1 * t3 * t12, 0.2e1 * t3 * t13 - 0.2e1 * t2 * t17; 0, 0, 0, 0, 0, t75, t74, t52, t33, -t34, -pkin(2) * t30 - t29 * t51 + t49 * t62, -pkin(2) * t31 + t29 * t49 + t51 * t62, (-t30 * t51 + t31 * t49) * qJ(3) + t59, -t29 * pkin(2) + t59 * qJ(3), t18 * t38, -t38 * t17 - t18 * t37, -t38 * t74, t37 * t74, 0, t43 * t17 + t19 * t37 + t21 * t74, t43 * t18 + t19 * t38 + t22 * t74, t13 * t68, (-t12 * t56 - t77) * t38, t13 * t37 + t17 * t68, -t12 * t37 - t17 * t70, t17 * t37, t1 * t37 + t21 * t12 + t9 * t17 + t3 * t70, -t10 * t17 + t21 * t13 - t2 * t37 + t3 * t68; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t51, -0.2e1 * pkin(2) * t49, 0.2e1 * t66 * qJ(3), t66 * qJ(3) ^ 2 + pkin(2) ^ 2, t35, t65, 0, 0, 0, t37 * t81, t38 * t81, t48 * t35, -0.2e1 * t35 * t69, 0.2e1 * t37 * t68, t53 * t65, t37 ^ 2, 0.2e1 * t21 * t70 + 0.2e1 * t9 * t37, -0.2e1 * t10 * t37 + 0.2e1 * t21 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t31, 0, t29, 0, 0, 0, 0, 0, t17, t18, 0, 0, 0, 0, 0, t16, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t49, 0, -pkin(2), 0, 0, 0, 0, 0, t37, t38, 0, 0, 0, 0, 0, t32, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, -t74, t5, -t6, t77, -t53 * t12 + t13 * t56, t72, t16, 0, -pkin(4) * t12 - pkin(9) * t72 - t3 * t56, -pkin(4) * t13 - pkin(9) * t16 + t3 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, 0, -t21, -t22, t53 * t68, (-t47 + t48) * t38, t71, t32, 0, -t21 * t56 + t60 * t53, t21 * t53 + t60 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t47, 0.2e1 * t69, 0, 0, 0, 0.2e1 * pkin(4) * t56, -0.2e1 * pkin(4) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, t17, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t70, t37, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t56, 0, -t53 * pkin(9), -t56 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t23;
