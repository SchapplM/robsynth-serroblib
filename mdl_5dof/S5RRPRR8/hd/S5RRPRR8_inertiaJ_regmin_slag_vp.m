% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR8
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
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t41 = cos(qJ(2));
t32 = -t41 * pkin(2) - pkin(1);
t35 = sin(pkin(9));
t36 = cos(pkin(9));
t39 = sin(qJ(2));
t44 = t35 * t39 - t36 * t41;
t18 = t44 * pkin(3) + t32;
t57 = 0.2e1 * t18;
t56 = 0.2e1 * t41;
t55 = t35 * pkin(2);
t40 = cos(qJ(5));
t38 = sin(qJ(4));
t48 = -qJ(3) - pkin(6);
t28 = t48 * t39;
t29 = t48 * t41;
t16 = t36 * t28 + t35 * t29;
t24 = t35 * t41 + t36 * t39;
t43 = -t24 * pkin(7) + t16;
t52 = cos(qJ(4));
t17 = t35 * t28 - t36 * t29;
t9 = -t44 * pkin(7) + t17;
t5 = t38 * t9 - t52 * t43;
t54 = t5 * t40;
t31 = t36 * pkin(2) + pkin(3);
t21 = t52 * t31 - t38 * t55;
t19 = -pkin(4) - t21;
t53 = pkin(4) - t19;
t14 = t38 * t24 + t52 * t44;
t37 = sin(qJ(5));
t11 = t37 * t14;
t15 = t52 * t24 - t38 * t44;
t51 = t37 * t15;
t50 = t37 * t40;
t49 = t40 * t15;
t47 = -0.2e1 * t15 * t14;
t46 = -pkin(4) * t15 - pkin(8) * t14;
t22 = -t38 * t31 - t52 * t55;
t20 = pkin(8) - t22;
t45 = -t14 * t20 + t15 * t19;
t34 = t40 ^ 2;
t33 = t37 ^ 2;
t30 = 0.2e1 * t50;
t13 = t15 ^ 2;
t12 = t40 * t14;
t10 = t37 * t49;
t7 = (-t33 + t34) * t15;
t6 = t38 * t43 + t52 * t9;
t4 = t14 * pkin(4) - t15 * pkin(8) + t18;
t3 = t5 * t37;
t2 = t37 * t4 + t40 * t6;
t1 = -t37 * t6 + t40 * t4;
t8 = [1, 0, 0, t39 ^ 2, t39 * t56, 0, 0, 0, pkin(1) * t56, -0.2e1 * pkin(1) * t39, -0.2e1 * t16 * t24 - 0.2e1 * t17 * t44, t16 ^ 2 + t17 ^ 2 + t32 ^ 2, t13, t47, 0, 0, 0, t14 * t57, t15 * t57, t34 * t13, -0.2e1 * t13 * t50, 0.2e1 * t14 * t49, t37 * t47, t14 ^ 2, 0.2e1 * t1 * t14 + 0.2e1 * t5 * t51, -0.2e1 * t2 * t14 + 0.2e1 * t5 * t49; 0, 0, 0, 0, 0, t39, t41, 0, -t39 * pkin(6), -t41 * pkin(6), (-t36 * t24 - t35 * t44) * pkin(2), (t16 * t36 + t17 * t35) * pkin(2), 0, 0, t15, -t14, 0, -t5, -t6, t10, t7, t11, t12, 0, t45 * t37 - t54, t45 * t40 + t3; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t35 ^ 2 + t36 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t21, 0.2e1 * t22, t33, t30, 0, 0, 0, -0.2e1 * t19 * t40, 0.2e1 * t19 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, -t5, -t6, t10, t7, t11, t12, 0, t46 * t37 - t54, t46 * t40 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t21, t22, t33, t30, 0, 0, 0, t53 * t40, -t53 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t33, t30, 0, 0, 0, 0.2e1 * pkin(4) * t40, -0.2e1 * pkin(4) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t51, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t40, 0, -t37 * t20, -t40 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t40, 0, -t37 * pkin(8), -t40 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
