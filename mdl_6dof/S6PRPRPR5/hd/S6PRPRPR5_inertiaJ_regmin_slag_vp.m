% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t31 = sin(pkin(11));
t33 = cos(pkin(11));
t36 = sin(qJ(4));
t39 = cos(qJ(4));
t20 = t36 * t31 - t39 * t33;
t21 = t39 * t31 + t36 * t33;
t26 = -t33 * pkin(3) - pkin(2);
t43 = -t21 * qJ(5) + t26;
t11 = t20 * pkin(4) + t43;
t62 = -0.2e1 * t11;
t61 = 0.2e1 * t26;
t60 = 0.2e1 * qJ(5);
t59 = pkin(4) + pkin(9);
t58 = t21 * t20;
t32 = sin(pkin(6));
t57 = t32 * sin(qJ(2));
t40 = cos(qJ(2));
t56 = t32 * t40;
t35 = sin(qJ(6));
t55 = t35 * t20;
t54 = t35 * t21;
t38 = cos(qJ(6));
t53 = t38 * t20;
t17 = t38 * t21;
t52 = t38 * t35;
t51 = pkin(8) + qJ(3);
t50 = t31 ^ 2 + t33 ^ 2;
t49 = qJ(5) * t20;
t48 = 0.2e1 * t58;
t47 = t20 * t56;
t46 = t21 * t56;
t34 = cos(pkin(6));
t15 = -t31 * t57 + t34 * t33;
t16 = t34 * t31 + t33 * t57;
t45 = -t15 * t31 + t16 * t33;
t22 = t51 * t31;
t23 = t51 * t33;
t12 = t39 * t22 + t36 * t23;
t13 = -t36 * t22 + t39 * t23;
t44 = t21 * t59 + t49;
t30 = t38 ^ 2;
t29 = t35 ^ 2;
t24 = t32 ^ 2 * t40 ^ 2;
t19 = t21 ^ 2;
t18 = t20 ^ 2;
t9 = t36 * t15 + t39 * t16;
t8 = -t39 * t15 + t36 * t16;
t7 = -t20 * pkin(5) + t13;
t6 = t21 * pkin(5) + t12;
t5 = t59 * t20 + t43;
t4 = -t35 * t8 + t38 * t56;
t3 = t35 * t56 + t38 * t8;
t2 = t35 * t6 + t38 * t5;
t1 = -t35 * t5 + t38 * t6;
t10 = [1, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t16 ^ 2 + t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 ^ 2 + t9 ^ 2 + t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, t56, -t57, t33 * t56, -t31 * t56, t45, pkin(2) * t56 + t45 * qJ(3), 0, 0, 0, 0, 0, -t47, -t46, -t9 * t20 + t8 * t21, t47, t46, -t11 * t56 + t8 * t12 + t9 * t13, 0, 0, 0, 0, 0, t3 * t21 - t9 * t53, t4 * t21 + t9 * t55; 0, 1, 0, 0, 0.2e1 * pkin(2) * t33, -0.2e1 * pkin(2) * t31, 0.2e1 * t50 * qJ(3), t50 * qJ(3) ^ 2 + pkin(2) ^ 2, t19, -0.2e1 * t58, 0, 0, 0, t20 * t61, t21 * t61, 0.2e1 * t12 * t21 - 0.2e1 * t13 * t20, t20 * t62, t21 * t62, t11 ^ 2 + t12 ^ 2 + t13 ^ 2, t29 * t18, 0.2e1 * t18 * t52, t35 * t48, t38 * t48, t19, 0.2e1 * t1 * t21 - 0.2e1 * t7 * t53, -0.2e1 * t2 * t21 + 0.2e1 * t7 * t55; 0, 0, 0, 0, 0, 0, 0, -t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t33, t31, 0, -pkin(2), 0, 0, 0, 0, 0, t20, t21, 0, -t20, -t21, t11, 0, 0, 0, 0, 0, -t54, -t17; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, 0, t8, t9, -t8 * pkin(4) + t9 * qJ(5), 0, 0, 0, 0, 0, t9 * t35, t9 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t12, -t13, -pkin(4) * t21 - t49, t12, t13, -t12 * pkin(4) + t13 * qJ(5), t20 * t52 (-t29 + t30) * t20, t17, -t54, 0, t7 * t35 - t44 * t38, t44 * t35 + t7 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t60, pkin(4) ^ 2 + qJ(5) ^ 2, t30, -0.2e1 * t52, 0, 0, 0, t35 * t60, t38 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, t12, 0, 0, 0, 0, 0, t17, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t53, t21, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, -t38 * t59, t35 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
