% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:31:25
% EndTime: 2019-05-05 09:31:27
% DurationCPUTime: 0.60s
% Computational Cost: add. (591->105), mult. (1234->179), div. (0->0), fcn. (1483->10), ass. (0->71)
t55 = cos(qJ(3));
t43 = -t55 * pkin(3) - pkin(2);
t77 = 0.2e1 * t43;
t76 = 0.2e1 * t55;
t75 = -pkin(9) - pkin(8);
t49 = sin(qJ(5));
t74 = t49 * pkin(5);
t50 = sin(qJ(4));
t73 = t50 * pkin(3);
t53 = cos(qJ(5));
t72 = t53 * pkin(10);
t54 = cos(qJ(4));
t71 = t54 * pkin(3);
t41 = -pkin(4) - t71;
t70 = pkin(4) - t41;
t48 = cos(pkin(6));
t51 = sin(qJ(3));
t47 = sin(pkin(6));
t67 = t47 * sin(qJ(2));
t22 = t48 * t55 - t51 * t67;
t23 = t48 * t51 + t55 * t67;
t11 = -t54 * t22 + t50 * t23;
t69 = t11 * t53;
t36 = t75 * t51;
t37 = t75 * t55;
t17 = -t54 * t36 - t50 * t37;
t68 = t17 * t53;
t66 = t47 * cos(qJ(2));
t31 = t50 * t55 + t54 * t51;
t65 = t49 * t31;
t64 = t49 * t53;
t18 = t50 * t36 - t54 * t37;
t63 = t53 * t18;
t62 = t53 * t31;
t40 = pkin(10) + t73;
t61 = t53 * t40;
t44 = t53 * qJ(6);
t30 = t50 * t51 - t54 * t55;
t60 = -0.2e1 * t31 * t30;
t42 = -t53 * pkin(5) - pkin(4);
t14 = t30 * pkin(4) - t31 * pkin(10) + t43;
t5 = t53 * t14 - t49 * t18;
t2 = t30 * pkin(5) - t31 * t44 + t5;
t4 = t63 + (-qJ(6) * t31 + t14) * t49;
t59 = -t2 * t49 + t4 * t53;
t58 = -pkin(4) * t31 - pkin(10) * t30;
t57 = -t30 * t40 + t31 * t41;
t46 = t53 ^ 2;
t45 = t49 ^ 2;
t38 = 0.2e1 * t64;
t35 = t44 + t72;
t34 = (-qJ(6) - pkin(10)) * t49;
t33 = t42 - t71;
t29 = t35 * t53;
t28 = t31 ^ 2;
t27 = t44 + t61;
t26 = (-qJ(6) - t40) * t49;
t25 = t53 * t30;
t24 = t49 * t30;
t21 = t27 * t53;
t20 = t49 * t62;
t16 = t17 * t49;
t15 = (-t45 + t46) * t31;
t12 = t50 * t22 + t54 * t23;
t10 = pkin(5) * t65 + t17;
t9 = t11 * t49;
t8 = t53 * t12 - t49 * t66;
t7 = -t49 * t12 - t53 * t66;
t6 = t49 * t14 + t63;
t1 = -t7 * t49 + t8 * t53;
t3 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 ^ 2 + t7 ^ 2 + t8 ^ 2; 0, 0, t66, -t67, 0, 0, 0, 0, 0, t55 * t66, -t51 * t66, 0, 0, 0, 0, 0, -t30 * t66, -t31 * t66, 0, 0, 0, 0, 0, t11 * t65 + t7 * t30, t11 * t62 - t8 * t30 (-t49 * t8 - t53 * t7) * t31, t11 * t10 + t7 * t2 + t8 * t4; 0, 1, 0, 0, t51 ^ 2, t51 * t76, 0, 0, 0, pkin(2) * t76, -0.2e1 * pkin(2) * t51, t28, t60, 0, 0, 0, t30 * t77, t31 * t77, t46 * t28, -0.2e1 * t28 * t64, 0.2e1 * t30 * t62, t49 * t60, t30 ^ 2, 0.2e1 * t17 * t65 + 0.2e1 * t5 * t30, 0.2e1 * t17 * t62 - 0.2e1 * t6 * t30, 0.2e1 * (-t2 * t53 - t4 * t49) * t31, t10 ^ 2 + t2 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t23, 0, 0, 0, 0, 0, -t11, -t12, 0, 0, 0, 0, 0, -t69, t9, t1, t11 * t33 + t7 * t26 + t8 * t27; 0, 0, 0, 0, 0, 0, t51, t55, 0, -t51 * pkin(8), -t55 * pkin(8), 0, 0, t31, -t30, 0, -t17, -t18, t20, t15, t24, t25, 0, t57 * t49 - t68, t57 * t53 + t16 (-t26 * t53 - t27 * t49) * t31 + t59, t10 * t33 + t2 * t26 + t4 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t71, -0.2e1 * t73, t45, t38, 0, 0, 0, -0.2e1 * t41 * t53, 0.2e1 * t41 * t49, -0.2e1 * t26 * t49 + 0.2e1 * t21, t26 ^ 2 + t27 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, 0, 0, 0, 0, -t69, t9, t1, t11 * t42 + t7 * t34 + t8 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, 0, -t17, -t18, t20, t15, t24, t25, 0, t58 * t49 - t68, t58 * t53 + t16 (-t34 * t53 - t35 * t49) * t31 + t59, t10 * t42 + t2 * t34 + t4 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t71, -t73, t45, t38, 0, 0, 0, t70 * t53, -t70 * t49, t21 + t29 + (-t26 - t34) * t49, t26 * t34 + t27 * t35 + t33 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t45, t38, 0, 0, 0, 0.2e1 * pkin(4) * t53, -0.2e1 * pkin(4) * t49, -0.2e1 * t34 * t49 + 0.2e1 * t29, t34 ^ 2 + t35 ^ 2 + t42 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t65, t30, t5, -t6, -pkin(5) * t62, t2 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t53, 0, -t49 * t40, -t61, -t74, t26 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t53, 0, -t49 * pkin(10), -t72, -t74, t34 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
