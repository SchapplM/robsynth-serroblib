% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t40 = sin(qJ(2));
t41 = cos(qJ(2));
t37 = sin(pkin(5));
t47 = qJ(3) * t37;
t25 = -t41 * pkin(2) - t40 * t47 - pkin(1);
t39 = cos(pkin(5));
t44 = qJ(3) * t39 + pkin(7);
t26 = t44 * t40;
t38 = cos(pkin(8));
t58 = (t25 * t37 - t26 * t39) * t38;
t57 = 0.2e1 * t37;
t56 = 0.2e1 * t41;
t55 = pkin(2) * t38;
t54 = t37 * pkin(2);
t36 = sin(pkin(8));
t53 = t36 * t39;
t33 = t37 * t36;
t34 = t37 * t38;
t52 = t37 * t41;
t51 = t38 * t39;
t50 = t39 * t41;
t49 = pkin(3) + qJ(5);
t27 = t44 * t41;
t23 = t36 * t27;
t48 = pkin(3) * t52 + t23;
t22 = pkin(2) * t53 + t38 * t47;
t8 = t25 * t33 - t26 * t53 + t38 * t27;
t46 = -pkin(3) - t55;
t45 = -qJ(4) * t36 - pkin(2);
t9 = t39 * t25 + t37 * t26;
t13 = -t39 * qJ(4) - t22;
t20 = t36 * t50 + t40 * t38;
t42 = -t20 * qJ(4) + t9;
t5 = qJ(4) * t52 - t8;
t35 = t37 ^ 2;
t29 = t36 * t47;
t21 = pkin(2) * t51 - t29;
t19 = t40 * t36 - t38 * t50;
t15 = (-pkin(3) * t38 + t45) * t37;
t14 = t46 * t39 + t29;
t12 = (-t49 * t38 + t45) * t37;
t11 = pkin(4) * t34 - t13;
t10 = pkin(4) * t33 + t29 + (-qJ(5) + t46) * t39;
t7 = -t23 + t58;
t6 = t48 - t58;
t4 = t19 * pkin(3) + t42;
t3 = -t19 * pkin(4) - t5;
t2 = t49 * t19 + t42;
t1 = t26 * t51 + t20 * pkin(4) + (qJ(5) * t41 - t25 * t38) * t37 + t48;
t16 = [1, 0, 0, t40 ^ 2, t40 * t56, 0, 0, 0, pkin(1) * t56, -0.2e1 * pkin(1) * t40, 0.2e1 * t9 * t19 - 0.2e1 * t7 * t52, 0.2e1 * t9 * t20 + 0.2e1 * t8 * t52, -0.2e1 * t8 * t19 - 0.2e1 * t7 * t20, t7 ^ 2 + t8 ^ 2 + t9 ^ 2, 0.2e1 * t5 * t19 + 0.2e1 * t6 * t20, -0.2e1 * t4 * t19 - 0.2e1 * t6 * t52, -0.2e1 * t4 * t20 + 0.2e1 * t5 * t52, t4 ^ 2 + t5 ^ 2 + t6 ^ 2, 0.2e1 * t1 * t20 - 0.2e1 * t3 * t19, -0.2e1 * t2 * t20 - 0.2e1 * t3 * t52, 0.2e1 * t1 * t52 + 0.2e1 * t2 * t19, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t40, t41, 0, -t40 * pkin(7), -t41 * pkin(7), t7 * t39 + (-pkin(2) * t19 - t21 * t41 - t38 * t9) * t37, -t8 * t39 + (-pkin(2) * t20 + t22 * t41 + t36 * t9) * t37, -t22 * t19 - t21 * t20 + (-t36 * t7 + t38 * t8) * t37, t7 * t21 + t8 * t22 - t9 * t54, t13 * t19 + t14 * t20 + (t36 * t6 - t38 * t5) * t37, -t15 * t19 + t6 * t39 + (-t14 * t41 + t38 * t4) * t37, -t15 * t20 - t5 * t39 + (t13 * t41 - t36 * t4) * t37, t5 * t13 + t6 * t14 + t4 * t15, t10 * t20 - t11 * t19 + (t1 * t36 + t3 * t38) * t37, -t12 * t20 + t3 * t39 + (-t11 * t41 - t2 * t36) * t37, -t1 * t39 + t12 * t19 + (t10 * t41 - t2 * t38) * t37, t1 * t10 + t3 * t11 + t2 * t12; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t21 * t39 + 0.2e1 * t35 * t55, -0.2e1 * t35 * pkin(2) * t36 - 0.2e1 * t22 * t39, (-t21 * t36 + t22 * t38) * t57, t35 * pkin(2) ^ 2 + t21 ^ 2 + t22 ^ 2, (-t13 * t38 + t14 * t36) * t57, 0.2e1 * t14 * t39 + 0.2e1 * t15 * t34, -0.2e1 * t13 * t39 - 0.2e1 * t15 * t33, t13 ^ 2 + t14 ^ 2 + t15 ^ 2, (t10 * t36 + t11 * t38) * t57, 0.2e1 * t11 * t39 - 0.2e1 * t12 * t33, -0.2e1 * t10 * t39 - 0.2e1 * t12 * t34, t10 ^ 2 + t11 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, t9, 0, -t19, -t20, t4, 0, -t20, t19, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33, 0, -t54, 0, t34, -t33, t15, 0, -t33, -t34, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t52, 0, t6, t20, 0, t52, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t39, 0, t14, t33, 0, -t39, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t52, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t39, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t16;
