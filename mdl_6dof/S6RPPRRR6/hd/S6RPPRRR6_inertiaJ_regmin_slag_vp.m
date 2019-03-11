% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t33 = sin(qJ(6));
t34 = sin(qJ(5));
t36 = cos(qJ(6));
t37 = cos(qJ(5));
t17 = t33 * t37 + t36 * t34;
t38 = cos(qJ(4));
t11 = t38 * t17;
t63 = -0.2e1 * t11;
t25 = -t37 * pkin(5) - pkin(4);
t62 = 0.2e1 * t25;
t31 = (pkin(1) + qJ(3));
t61 = 2 * t31;
t35 = sin(qJ(4));
t60 = 0.2e1 * t35;
t59 = 0.2e1 * t37;
t58 = pkin(8) + pkin(9);
t57 = t33 * pkin(5);
t56 = t36 * pkin(5);
t18 = t35 * pkin(4) - t38 * pkin(8) + t31;
t30 = -pkin(7) + qJ(2);
t48 = t37 * t35;
t42 = t30 * t48;
t5 = t42 + (-pkin(9) * t38 + t18) * t34;
t55 = t36 * t5;
t16 = t33 * t34 - t36 * t37;
t54 = t16 * t35;
t53 = t17 * t35;
t52 = t30 * t34;
t22 = t34 * t35;
t51 = t34 * t37;
t50 = t34 * t38;
t49 = t35 * t30;
t24 = t37 * t38;
t47 = t38 * t16;
t46 = t38 * t30;
t45 = t38 * t35;
t27 = t35 ^ 2;
t29 = t38 ^ 2;
t44 = -t27 - t29;
t43 = -0.2e1 * t45;
t14 = t37 * t18;
t4 = -pkin(9) * t24 + t14 + (pkin(5) - t52) * t35;
t1 = -t33 * t5 + t36 * t4;
t41 = -pkin(4) * t38 - pkin(8) * t35;
t40 = (qJ(2) ^ 2);
t39 = 2 * qJ(2);
t28 = t37 ^ 2;
t26 = t34 ^ 2;
t20 = t58 * t37;
t19 = t58 * t34;
t15 = (pkin(5) * t34 - t30) * t38;
t12 = -t33 * t22 + t36 * t48;
t9 = t34 * t18 + t42;
t8 = -t34 * t49 + t14;
t7 = -t33 * t19 + t36 * t20;
t6 = -t36 * t19 - t33 * t20;
t2 = t33 * t4 + t55;
t3 = [1, 0, 0, -2 * pkin(1), t39, pkin(1) ^ 2 + t40, t39, t61, t31 ^ 2 + t40, t29, t43, 0, 0, 0, t31 * t60, t38 * t61, t28 * t29, -0.2e1 * t29 * t51, t45 * t59, t34 * t43, t27, -0.2e1 * t29 * t52 + 0.2e1 * t8 * t35, -0.2e1 * t29 * t30 * t37 - 0.2e1 * t9 * t35, t47 ^ 2, -t47 * t63, -t47 * t60, t35 * t63, t27, 0.2e1 * t1 * t35 + 0.2e1 * t15 * t11, -0.2e1 * t15 * t47 - 0.2e1 * t2 * t35; 0, 0, 0, 1, 0, -pkin(1), 0, -1, -t31, 0, 0, 0, 0, 0, -t35, -t38, 0, 0, 0, 0, 0, -t48, t22, 0, 0, 0, 0, 0, t54, t53; 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, qJ(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t34, t44 * t37, 0, 0, 0, 0, 0, -t38 * t11 - t35 * t53, -t12 * t35 + t38 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, t46, -t49, t34 * t24 (-t26 + t28) * t38, t22, t48, 0, t41 * t34 + t37 * t46, -t34 * t46 + t41 * t37, -t47 * t17, -t17 * t11 + t16 * t47, t53, -t54, 0, t25 * t11 + t15 * t16 + t6 * t35, t15 * t17 - t25 * t47 - t7 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, 0, 0, 0, 0, t24, -t50, 0, 0, 0, 0, 0, -t47, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t26, 0.2e1 * t51, 0, 0, 0, pkin(4) * t59, -0.2e1 * pkin(4) * t34, t17 ^ 2, -0.2e1 * t17 * t16, 0, 0, 0, t16 * t62, t17 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t50, t35, t8, -t9, 0, 0, -t47, -t11, t35, t35 * t56 + t1, -t55 + (-t35 * pkin(5) - t4) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t34, 0, 0, 0, 0, 0, t16, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t48, 0, 0, 0, 0, 0, -t53, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t37, 0, -t34 * pkin(8), -t37 * pkin(8), 0, 0, t17, -t16, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t11, t35, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t56, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
