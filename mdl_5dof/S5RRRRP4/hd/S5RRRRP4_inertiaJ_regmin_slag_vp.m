% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t36 = cos(qJ(3));
t33 = sin(qJ(4));
t34 = sin(qJ(3));
t44 = t33 * t34;
t45 = cos(qJ(4));
t17 = -t45 * t36 + t44;
t59 = 0.2e1 * t17;
t40 = t45 * t34;
t18 = t33 * t36 + t40;
t58 = -0.2e1 * t18;
t29 = -t36 * pkin(3) - pkin(2);
t48 = cos(qJ(2)) * pkin(1);
t20 = t29 - t48;
t57 = 0.2e1 * t20;
t56 = 0.2e1 * t29;
t55 = 0.2e1 * t36;
t54 = -pkin(8) - pkin(7);
t31 = t36 * pkin(8);
t50 = sin(qJ(2)) * pkin(1);
t25 = pkin(7) + t50;
t43 = t36 * t25;
t15 = t31 + t43;
t46 = -pkin(8) - t25;
t7 = t33 * t15 - t46 * t40;
t8 = t45 * t15 + t46 * t44;
t53 = -t8 * t17 + t7 * t18;
t49 = t36 * pkin(7);
t21 = t31 + t49;
t12 = t33 * t21 - t54 * t40;
t13 = t45 * t21 + t54 * t44;
t52 = t12 * t18 - t13 * t17;
t9 = t17 * pkin(4) - t18 * qJ(5) + t29;
t6 = t9 - t48;
t51 = t6 + t9;
t28 = -pkin(2) - t48;
t47 = pkin(2) - t28;
t42 = t20 + t29;
t41 = t45 * pkin(3);
t39 = 0.2e1 * pkin(4);
t38 = 0.2e1 * qJ(5);
t32 = t34 ^ 2;
t30 = t33 * pkin(3);
t26 = t41 + pkin(4);
t23 = t30 + qJ(5);
t22 = t34 * t55;
t16 = t18 ^ 2;
t11 = t17 * t58;
t10 = -pkin(4) * t18 - t17 * qJ(5);
t3 = -t23 * t17 - t26 * t18;
t1 = [1, 0, 0, 1, 0.2e1 * t48, -0.2e1 * t50, t32, t22, 0, 0, 0, -0.2e1 * t28 * t36, 0.2e1 * t28 * t34, t16, t11, 0, 0, 0, t17 * t57, t18 * t57, t6 * t59, 0.2e1 * t53, t6 * t58, t6 ^ 2 + t7 ^ 2 + t8 ^ 2; 0, 0, 0, 1, t48, -t50, t32, t22, 0, 0, 0, t47 * t36, -t47 * t34, t16, t11, 0, 0, 0, t42 * t17, t42 * t18, t51 * t17, t52 + t53, -t51 * t18, t7 * t12 + t8 * t13 + t6 * t9; 0, 0, 0, 1, 0, 0, t32, t22, 0, 0, 0, pkin(2) * t55, -0.2e1 * pkin(2) * t34, t16, t11, 0, 0, 0, t17 * t56, t18 * t56, t9 * t59, 0.2e1 * t52, t9 * t58, t12 ^ 2 + t13 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, t34, t36, 0, -t34 * t25, -t43, 0, 0, t18, -t17, 0, -t7, -t8, -t7, t3, t8, t8 * t23 - t7 * t26; 0, 0, 0, 0, 0, 0, 0, 0, t34, t36, 0, -t34 * pkin(7), -t49, 0, 0, t18, -t17, 0, -t12, -t13, -t12, t3, t13, -t12 * t26 + t13 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t30, 0.2e1 * t26, 0, 0.2e1 * t23, t23 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, -t7, -t8, -t7, t10, t8, -t7 * pkin(4) + t8 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, -t12, -t13, -t12, t10, t13, -t12 * pkin(4) + t13 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t41, -t30, t39 + t41, 0, t38 + t30, t26 * pkin(4) + t23 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t39, 0, t38, pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
