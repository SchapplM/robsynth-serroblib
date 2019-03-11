% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t37 = sin(qJ(3));
t40 = cos(qJ(3));
t49 = cos(pkin(6));
t34 = sin(pkin(6));
t60 = t34 * sin(qJ(2));
t18 = t49 * t37 + t40 * t60;
t33 = sin(pkin(11));
t35 = cos(pkin(11));
t43 = -t37 * t60 + t49 * t40;
t8 = t33 * t18 - t35 * t43;
t63 = t8 ^ 2;
t62 = 0.2e1 * t40;
t39 = cos(qJ(5));
t61 = t39 * pkin(5);
t41 = cos(qJ(2));
t59 = t34 * t41;
t22 = t33 * t37 - t35 * t40;
t36 = sin(qJ(5));
t58 = t36 * t22;
t23 = t33 * t40 + t35 * t37;
t57 = t36 * t23;
t56 = t36 * t39;
t53 = -qJ(4) - pkin(8);
t26 = t53 * t40;
t48 = t53 * t37;
t15 = -t35 * t26 + t33 * t48;
t55 = t39 * t15;
t54 = t39 * t23;
t31 = t36 ^ 2;
t32 = t39 ^ 2;
t52 = t31 + t32;
t51 = qJ(6) * t23;
t28 = t33 * pkin(3) + pkin(9);
t50 = qJ(6) + t28;
t30 = -t40 * pkin(3) - pkin(2);
t29 = -t35 * pkin(3) - pkin(4);
t12 = t22 * pkin(4) - t23 * pkin(9) + t30;
t3 = t39 * t12 - t36 * t15;
t13 = -t33 * t26 - t35 * t48;
t1 = t22 * pkin(5) - t39 * t51 + t3;
t2 = t55 + (t12 - t51) * t36;
t47 = t1 * t39 + t2 * t36;
t10 = t35 * t18 + t33 * t43;
t5 = -t36 * t10 - t39 * t59;
t6 = t39 * t10 - t36 * t59;
t46 = t6 * t36 + t5 * t39;
t19 = t50 * t36;
t20 = t50 * t39;
t45 = -t19 * t39 + t20 * t36;
t44 = -t22 * t28 + t23 * t29;
t25 = t29 - t61;
t21 = t23 ^ 2;
t17 = t39 * t22;
t7 = pkin(5) * t57 + t13;
t4 = t36 * t12 + t55;
t9 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 ^ 2 * t41 ^ 2 + t10 ^ 2 + t63, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t63; 0, 0, t59, -t60, 0, 0, 0, 0, 0, t40 * t59, -t37 * t59, -t10 * t22 + t8 * t23, t10 * t15 + t8 * t13 - t30 * t59, 0, 0, 0, 0, 0, t5 * t22 + t8 * t57, -t6 * t22 + t8 * t54, -t46 * t23, t5 * t1 + t6 * t2 + t8 * t7; 0, 1, 0, 0, t37 ^ 2, t37 * t62, 0, 0, 0, pkin(2) * t62, -0.2e1 * pkin(2) * t37, 0.2e1 * t13 * t23 - 0.2e1 * t15 * t22, t13 ^ 2 + t15 ^ 2 + t30 ^ 2, t32 * t21, -0.2e1 * t21 * t56, 0.2e1 * t22 * t54, -0.2e1 * t22 * t57, t22 ^ 2, 0.2e1 * t13 * t57 + 0.2e1 * t3 * t22, 0.2e1 * t13 * t54 - 0.2e1 * t4 * t22, -0.2e1 * t47 * t23, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t18, 0 (t10 * t33 - t35 * t8) * pkin(3), 0, 0, 0, 0, 0, -t8 * t39, t8 * t36, -t5 * t36 + t6 * t39, -t5 * t19 + t6 * t20 + t8 * t25; 0, 0, 0, 0, 0, 0, t37, t40, 0, -t37 * pkin(8), -t40 * pkin(8) (-t22 * t33 - t23 * t35) * pkin(3) (-t13 * t35 + t15 * t33) * pkin(3), t36 * t54 (-t31 + t32) * t23, t58, t17, 0, -t13 * t39 + t44 * t36, t13 * t36 + t44 * t39, -t1 * t36 + t2 * t39 - t45 * t23, -t1 * t19 + t2 * t20 + t7 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t33 ^ 2 + t35 ^ 2) * pkin(3) ^ 2, t31, 0.2e1 * t56, 0, 0, 0, -0.2e1 * t29 * t39, 0.2e1 * t29 * t36, 0.2e1 * t19 * t36 + 0.2e1 * t20 * t39, t19 ^ 2 + t20 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, t17, -t58, -t52 * t23, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t5 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t57, t22, t3, -t4, -pkin(5) * t54, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t39, 0, -t36 * t28, -t39 * t28, -t36 * pkin(5), -t19 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t36, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t9;
