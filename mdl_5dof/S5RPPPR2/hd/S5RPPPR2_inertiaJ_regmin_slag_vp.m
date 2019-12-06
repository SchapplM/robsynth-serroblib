% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t38 = sin(pkin(7));
t59 = -0.2e1 * t38;
t36 = sin(pkin(9));
t39 = cos(pkin(9));
t41 = cos(pkin(7));
t40 = cos(pkin(8));
t52 = t40 * t38;
t18 = -t41 * t36 + t39 * t52;
t42 = sin(qJ(5));
t37 = sin(pkin(8));
t43 = cos(qJ(5));
t54 = t37 * t43;
t9 = t18 * t42 - t38 * t54;
t58 = -0.2e1 * t9;
t22 = -t41 * pkin(2) - t38 * qJ(3) - pkin(1);
t47 = qJ(2) * t41;
t15 = t37 * t22 + t40 * t47;
t11 = -t41 * qJ(4) + t15;
t28 = t38 * qJ(2);
t16 = t28 + (pkin(3) * t37 - qJ(4) * t40) * t38;
t7 = t39 * t11 + t36 * t16;
t57 = t36 * t37;
t56 = t37 * t38;
t55 = t37 * t42;
t53 = t38 * t39;
t51 = t42 * t36;
t50 = t43 * t36;
t49 = t36 ^ 2 + t39 ^ 2;
t31 = t37 ^ 2;
t34 = t40 ^ 2;
t48 = t31 + t34;
t32 = t38 ^ 2;
t46 = t32 * qJ(2);
t14 = t40 * t22 - t37 * t47;
t12 = t41 * pkin(3) - t14;
t6 = -t36 * t11 + t39 * t16;
t45 = t14 * t40 + t15 * t37;
t44 = qJ(2) ^ 2;
t35 = t41 ^ 2;
t27 = t32 * t44;
t20 = t39 * t54 - t42 * t40;
t19 = -t39 * t55 - t43 * t40;
t17 = t36 * t52 + t41 * t39;
t10 = t18 * t43 + t38 * t55;
t5 = pkin(6) * t56 + t7;
t4 = -pkin(4) * t56 - t6;
t3 = t17 * pkin(4) - t18 * pkin(6) + t12;
t2 = t42 * t3 + t43 * t5;
t1 = t43 * t3 - t42 * t5;
t8 = [1, 0, 0, 0.2e1 * pkin(1) * t41, pkin(1) * t59, 0.2e1 * (t32 + t35) * qJ(2), pkin(1) ^ 2 + t35 * t44 + t27, -0.2e1 * t14 * t41 + 0.2e1 * t37 * t46, 0.2e1 * t15 * t41 + 0.2e1 * t40 * t46, t45 * t59, t14 ^ 2 + t15 ^ 2 + t27, 0.2e1 * t12 * t17 + 0.2e1 * t6 * t56, 0.2e1 * t12 * t18 - 0.2e1 * t7 * t56, -0.2e1 * t7 * t17 - 0.2e1 * t6 * t18, t12 ^ 2 + t6 ^ 2 + t7 ^ 2, t10 ^ 2, t10 * t58, 0.2e1 * t10 * t17, t17 * t58, t17 ^ 2, 0.2e1 * t1 * t17 + 0.2e1 * t4 * t9, 0.2e1 * t4 * t10 - 0.2e1 * t2 * t17; 0, 0, 0, -t41, t38, 0, -pkin(1), -t40 * t41, t37 * t41, -t48 * t38, t45, -t36 * t31 * t38 - t40 * t17, -t40 * t18 - t31 * t53, (-t17 * t39 + t18 * t36) * t37, -t12 * t40 + (-t36 * t6 + t39 * t7) * t37, 0, 0, 0, 0, 0, t19 * t17 + t9 * t57, t10 * t57 - t20 * t17; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t48, 0, 0, 0, t49 * t31 + t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t56, t52, 0, t28, t37 * t53, -t36 * t56, -t36 * t17 - t39 * t18, t7 * t36 + t6 * t39, 0, 0, 0, 0, 0, -t17 * t51 - t39 * t9, -t39 * t10 - t17 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, 0, t12, 0, 0, 0, 0, 0, t43 * t17, -t42 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, t17, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
