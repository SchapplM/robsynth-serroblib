% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t37 = cos(pkin(10));
t27 = -t37 * pkin(1) - pkin(2);
t36 = cos(pkin(11));
t22 = -t36 * pkin(3) + t27;
t34 = sin(pkin(11));
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t48 = -t40 * t34 + t42 * t36;
t16 = -t48 * pkin(4) + t22;
t62 = 0.2e1 * t16;
t21 = t42 * t34 + t40 * t36;
t61 = 0.2e1 * t21;
t39 = sin(qJ(5));
t60 = t39 * pkin(4);
t35 = sin(pkin(10));
t25 = t35 * pkin(1) + qJ(3);
t57 = pkin(7) + t25;
t19 = t57 * t34;
t20 = t57 * t36;
t49 = -t42 * t19 - t40 * t20;
t44 = -t21 * pkin(8) + t49;
t56 = cos(qJ(5));
t45 = t40 * t19 - t42 * t20;
t9 = t48 * pkin(8) - t45;
t4 = t39 * t9 - t56 * t44;
t41 = cos(qJ(6));
t59 = t4 * t41;
t50 = t56 * pkin(4);
t29 = -t50 - pkin(5);
t58 = pkin(5) - t29;
t14 = t39 * t21 - t56 * t48;
t38 = sin(qJ(6));
t11 = t38 * t14;
t15 = t56 * t21 + t39 * t48;
t55 = t38 * t15;
t54 = t38 * t41;
t53 = t41 * t15;
t52 = t34 ^ 2 + t36 ^ 2;
t51 = -0.2e1 * t15 * t14;
t47 = -pkin(5) * t15 - pkin(9) * t14;
t28 = pkin(9) + t60;
t46 = -t14 * t28 + t15 * t29;
t33 = t41 ^ 2;
t32 = t38 ^ 2;
t24 = 0.2e1 * t54;
t13 = t15 ^ 2;
t12 = t41 * t14;
t10 = t38 * t53;
t7 = (-t32 + t33) * t15;
t6 = t14 * pkin(5) - t15 * pkin(9) + t16;
t5 = t39 * t44 + t56 * t9;
t3 = t4 * t38;
t2 = t38 * t6 + t41 * t5;
t1 = -t38 * t5 + t41 * t6;
t8 = [1, 0, 0 (t35 ^ 2 + t37 ^ 2) * pkin(1) ^ 2, -0.2e1 * t27 * t36, 0.2e1 * t27 * t34, 0.2e1 * t52 * t25, t52 * t25 ^ 2 + t27 ^ 2, t21 ^ 2, t48 * t61, 0, 0, 0, -0.2e1 * t22 * t48, t22 * t61, t13, t51, 0, 0, 0, t14 * t62, t15 * t62, t33 * t13, -0.2e1 * t13 * t54, 0.2e1 * t14 * t53, t38 * t51, t14 ^ 2, 0.2e1 * t1 * t14 + 0.2e1 * t4 * t55, -0.2e1 * t2 * t14 + 0.2e1 * t4 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t36, t34, 0, t27, 0, 0, 0, 0, 0, -t48, t21, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t48, 0, t49, t45, 0, 0, t15, -t14, 0, -t4, -t5, t10, t7, t11, t12, 0, t46 * t38 - t59, t46 * t41 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t21, 0, 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, -t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t50, -0.2e1 * t60, t32, t24, 0, 0, 0, -0.2e1 * t29 * t41, 0.2e1 * t29 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, -t4, -t5, t10, t7, t11, t12, 0, t47 * t38 - t59, t47 * t41 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, -t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t50, -t60, t32, t24, 0, 0, 0, t58 * t41, -t58 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t32, t24, 0, 0, 0, 0.2e1 * pkin(5) * t41, -0.2e1 * pkin(5) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t55, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t41, 0, -t38 * t28, -t41 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t41, 0, -t38 * pkin(9), -t41 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
