% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR13_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t32 = sin(pkin(8));
t33 = cos(pkin(8));
t35 = sin(qJ(3));
t58 = cos(qJ(3));
t17 = t35 * t32 - t58 * t33;
t19 = t58 * t32 + t35 * t33;
t26 = -t33 * pkin(2) - pkin(1);
t41 = -t19 * qJ(4) + t26;
t7 = t17 * pkin(3) + t41;
t63 = -0.2e1 * t7;
t14 = t17 ^ 2;
t15 = t19 ^ 2;
t62 = 0.2e1 * t26;
t61 = 0.2e1 * t33;
t60 = 0.2e1 * qJ(4);
t59 = pkin(3) + pkin(7);
t57 = t19 * t17;
t34 = sin(qJ(5));
t29 = t34 ^ 2;
t56 = t29 * t17;
t55 = t34 * t17;
t54 = t34 * t19;
t36 = cos(qJ(5));
t53 = t36 * t17;
t13 = t36 * t19;
t52 = t36 * t34;
t51 = pkin(6) + qJ(2);
t27 = t32 ^ 2;
t28 = t33 ^ 2;
t50 = t27 + t28;
t30 = t36 ^ 2;
t49 = t29 + t30;
t48 = qJ(4) * t17;
t47 = -0.2e1 * t57;
t46 = 0.2e1 * t57;
t45 = t17 * t52;
t22 = t51 * t32;
t23 = t51 * t33;
t10 = -t35 * t22 + t58 * t23;
t8 = t58 * t22 + t35 * t23;
t44 = t10 ^ 2 + t8 ^ 2;
t4 = t59 * t17 + t41;
t5 = t19 * pkin(4) + t8;
t2 = -t34 * t4 + t36 * t5;
t3 = t34 * t5 + t36 * t4;
t1 = t2 * t36 + t3 * t34;
t43 = -t2 * t34 + t3 * t36;
t42 = t19 * t59 + t48;
t40 = -0.2e1 * t10 * t17 + 0.2e1 * t8 * t19;
t38 = qJ(4) ^ 2;
t21 = t49 * t59;
t12 = t30 * t17;
t6 = -t17 * pkin(4) + t10;
t9 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t27, t32 * t61, 0, t28, 0, 0, pkin(1) * t61, -0.2e1 * pkin(1) * t32, 0.2e1 * t50 * qJ(2), t50 * qJ(2) ^ 2 + pkin(1) ^ 2, t15, t47, 0, t14, 0, 0, t17 * t62, t19 * t62, t40, t26 ^ 2 + t44, 0, 0, 0, t15, t47, t14, t40, t17 * t63, t19 * t63, t7 ^ 2 + t44, t29 * t14, 0.2e1 * t14 * t52, t34 * t46, t30 * t14, t36 * t46, t15, 0.2e1 * t2 * t19 - 0.2e1 * t6 * t53, -0.2e1 * t3 * t19 + 0.2e1 * t6 * t55, 0.2e1 * t43 * t17, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t17, t19, 0, t26, 0, 0, 0, 0, 0, 0, 0, -t17, -t19, t7, 0, 0, 0, 0, 0, 0, -t54, -t13, t12 + t56, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, 0, -t8, -t10, 0, 0, 0, -t19, t17, 0, 0, 0, -pkin(3) * t19 - t48, t8, t10, -t8 * pkin(3) + t10 * qJ(4), t45, t12 - t56, t13, -t45, -t54, 0, t6 * t34 - t36 * t42, t34 * t42 + t6 * t36, -t1, t6 * qJ(4) - t1 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(3), t60, pkin(3) ^ 2 + t38, t30, -0.2e1 * t52, 0, t29, 0, 0, t34 * t60, t36 * t60, 0.2e1 * t21, t49 * t59 ^ 2 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, t8, 0, 0, 0, 0, 0, 0, t13, -t54, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t53, t19, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, -t34, 0, -t36 * t59, t34 * t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t9;
