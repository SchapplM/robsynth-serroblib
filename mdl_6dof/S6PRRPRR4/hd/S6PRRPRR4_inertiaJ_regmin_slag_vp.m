% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:10:59
% EndTime: 2019-05-05 05:11:00
% DurationCPUTime: 0.55s
% Computational Cost: add. (291->86), mult. (625->144), div. (0->0), fcn. (770->10), ass. (0->70)
t38 = sin(qJ(5));
t39 = sin(qJ(3));
t42 = cos(qJ(5));
t43 = cos(qJ(3));
t15 = t38 * t39 + t42 * t43;
t76 = 0.2e1 * t15;
t37 = sin(qJ(6));
t75 = -0.2e1 * t37;
t74 = -0.2e1 * t39;
t41 = cos(qJ(6));
t73 = 0.2e1 * t41;
t72 = 0.2e1 * t43;
t56 = cos(pkin(6));
t36 = sin(pkin(6));
t66 = t36 * sin(qJ(2));
t11 = t39 * t66 - t43 * t56;
t12 = t39 * t56 + t43 * t66;
t6 = -t11 * t42 + t12 * t38;
t71 = t6 * t37;
t70 = t6 * t41;
t28 = t39 * pkin(8);
t22 = -pkin(9) * t39 + t28;
t29 = t43 * pkin(8);
t23 = -pkin(9) * t43 + t29;
t9 = -t22 * t42 + t23 * t38;
t69 = t9 * t37;
t68 = t9 * t41;
t45 = -pkin(3) - pkin(4);
t19 = qJ(4) * t38 - t42 * t45;
t17 = pkin(5) + t19;
t67 = pkin(5) + t17;
t44 = cos(qJ(2));
t65 = t36 * t44;
t64 = t37 * t15;
t16 = -t38 * t43 + t39 * t42;
t63 = t37 * t16;
t62 = t37 * t41;
t61 = t41 * t15;
t60 = t41 * t16;
t59 = t42 * t37;
t58 = t42 * t41;
t33 = t39 ^ 2;
t57 = t43 ^ 2 + t33;
t55 = -0.2e1 * t16 * t15;
t54 = -0.2e1 * t62;
t21 = -pkin(3) * t43 - qJ(4) * t39 - pkin(2);
t53 = t39 * t65;
t52 = t37 * t60;
t13 = pkin(4) * t43 - t21;
t51 = -pkin(5) * t16 - pkin(10) * t15;
t50 = -pkin(3) * t39 + qJ(4) * t43;
t49 = t11 * t39 + t12 * t43;
t20 = qJ(4) * t42 + t38 * t45;
t18 = -pkin(10) + t20;
t48 = -t15 * t18 + t16 * t17;
t47 = -t15 * t38 - t16 * t42;
t34 = t41 ^ 2;
t32 = t37 ^ 2;
t25 = 0.2e1 * t62;
t24 = t43 * t65;
t14 = t16 ^ 2;
t10 = t22 * t38 + t23 * t42;
t8 = (t32 - t34) * t16;
t7 = t11 * t38 + t12 * t42;
t5 = pkin(5) * t15 - pkin(10) * t16 + t13;
t4 = t37 * t65 + t41 * t7;
t3 = -t37 * t7 + t41 * t65;
t2 = t10 * t41 + t37 * t5;
t1 = -t10 * t37 + t41 * t5;
t26 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 ^ 2 * t44 ^ 2 + t11 ^ 2 + t12 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t65, -t66, 0, 0, 0, 0, 0, t24, -t53, t24, t49, t53, pkin(8) * t49 - t21 * t65, 0, 0, 0, 0, 0, t15 * t65, t16 * t65, 0, 0, 0, 0, 0, t15 * t3 + t6 * t63, -t15 * t4 + t6 * t60; 0, 1, 0, 0, t33, t39 * t72, 0, 0, 0, pkin(2) * t72, pkin(2) * t74, -0.2e1 * t21 * t43, 0.2e1 * t57 * pkin(8), t21 * t74, pkin(8) ^ 2 * t57 + t21 ^ 2, t14, t55, 0, 0, 0, t13 * t76, 0.2e1 * t13 * t16, t34 * t14, t14 * t54, t60 * t76, t37 * t55, t15 ^ 2, 0.2e1 * t1 * t15 + 0.2e1 * t63 * t9, -0.2e1 * t15 * t2 + 0.2e1 * t60 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -t11, 0, t12, -pkin(3) * t11 + qJ(4) * t12, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, t70, -t71; 0, 0, 0, 0, 0, 0, t39, t43, 0, -t28, -t29, -t28, t50, t29, t50 * pkin(8), 0, 0, -t16, t15, 0, t9, t10, -t52, t8, -t64, -t61, 0, t37 * t48 + t68, t41 * t48 - t69; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t19, 0.2e1 * t20, t32, t25, 0, 0, 0, t17 * t73, t17 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t37, t47 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, -t42, t38, 0, 0, 0, 0, 0, -t58, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t70, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, -t9, -t10, t52, -t8, t64, t61, 0, t37 * t51 - t68, t41 * t51 + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t19, -t20, -t32, t54, 0, 0, 0, -t67 * t41, t67 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t38, 0, 0, 0, 0, 0, t58, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t32, t25, 0, 0, 0, pkin(5) * t73, pkin(5) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t63, t15, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t41, 0, -t37 * t18, -t41 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t38, -t41 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t41, 0, -t37 * pkin(10), -t41 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t26;
