% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t39 = sin(pkin(8));
t29 = t39 * pkin(2) + qJ(4);
t38 = sin(pkin(9));
t36 = t38 ^ 2;
t40 = cos(pkin(9));
t37 = t40 ^ 2;
t48 = t36 + t37;
t49 = t48 * t29;
t61 = 0.2e1 * t38;
t60 = -0.2e1 * t40;
t35 = cos(qJ(2)) * pkin(1);
t34 = t35 + pkin(2);
t41 = cos(pkin(8));
t57 = sin(qJ(2)) * pkin(1);
t16 = t39 * t34 + t41 * t57;
t13 = qJ(4) + t16;
t42 = sin(qJ(5));
t44 = cos(qJ(5));
t52 = t44 * t40;
t15 = t41 * t34 - t39 * t57;
t47 = -t40 * pkin(4) - t38 * pkin(7) - pkin(3);
t7 = -t15 + t47;
t3 = t13 * t52 + t42 * t7;
t54 = t36 * t44;
t59 = t13 * t54 + t3 * t40;
t58 = t41 * pkin(2);
t17 = t47 - t58;
t6 = t42 * t17 + t29 * t52;
t56 = t29 * t54 + t6 * t40;
t55 = t36 * t42;
t53 = t42 * t38;
t30 = t42 * t40;
t51 = t48 * t13;
t14 = -pkin(3) - t15;
t33 = -pkin(3) - t58;
t50 = t14 + t33;
t32 = t44 * t38;
t31 = t44 ^ 2 * t36;
t24 = -0.2e1 * t42 * t54;
t23 = -0.2e1 * t38 * t52;
t22 = t30 * t61;
t18 = t29 * t55;
t8 = t13 * t55;
t5 = t44 * t17 - t29 * t30;
t2 = -t13 * t30 + t44 * t7;
t1 = [1, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t57, t15 ^ 2 + t16 ^ 2, t14 * t60, t14 * t61, 0.2e1 * t51, t48 * t13 ^ 2 + t14 ^ 2, t31, t24, t23, t22, t37, -0.2e1 * t2 * t40 + 0.2e1 * t8, 0.2e1 * t59; 0, 0, 0, 1, t35, -t57, (t15 * t41 + t16 * t39) * pkin(2), -t50 * t40, t50 * t38, t49 + t51, t13 * t49 + t14 * t33, t31, t24, t23, t22, t37, t18 + t8 + (-t2 - t5) * t40, t56 + t59; 0, 0, 0, 1, 0, 0, (t39 ^ 2 + t41 ^ 2) * pkin(2) ^ 2, t33 * t60, t33 * t61, 0.2e1 * t49, t48 * t29 ^ 2 + t33 ^ 2, t31, t24, t23, t22, t37, -0.2e1 * t5 * t40 + 0.2e1 * t18, 0.2e1 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t48, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t40, t38, 0, t14, 0, 0, 0, 0, 0, -t52, t30; 0, 0, 0, 0, 0, 0, 0, -t40, t38, 0, t33, 0, 0, 0, 0, 0, -t52, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t53, -t40, t2, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t53, -t40, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t1;
