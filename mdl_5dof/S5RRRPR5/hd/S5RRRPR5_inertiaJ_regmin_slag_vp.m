% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t46 = cos(qJ(2));
t35 = -t46 * pkin(2) - pkin(1);
t61 = 0.2e1 * t35;
t41 = sin(qJ(5));
t60 = 0.2e1 * t41;
t44 = cos(qJ(5));
t59 = -0.2e1 * t44;
t58 = 0.2e1 * t46;
t57 = pkin(6) + pkin(7);
t42 = sin(qJ(3));
t56 = t42 * pkin(2);
t43 = sin(qJ(2));
t28 = t57 * t43;
t29 = t57 * t46;
t45 = cos(qJ(3));
t18 = t42 * t28 - t45 * t29;
t25 = t42 * t43 - t45 * t46;
t10 = -t25 * qJ(4) - t18;
t39 = sin(pkin(9));
t40 = cos(pkin(9));
t17 = -t45 * t28 - t42 * t29;
t26 = t42 * t46 + t45 * t43;
t48 = -t26 * qJ(4) + t17;
t5 = t39 * t10 - t40 * t48;
t55 = t5 * t44;
t15 = t40 * t25 + t39 * t26;
t12 = t41 * t15;
t16 = -t39 * t25 + t40 * t26;
t54 = t41 * t16;
t53 = t41 * t44;
t52 = t44 * t16;
t36 = t45 * pkin(2);
t34 = t36 + pkin(3);
t22 = t40 * t34 - t39 * t56;
t20 = -pkin(4) - t22;
t33 = -t40 * pkin(3) - pkin(4);
t51 = t20 + t33;
t23 = t39 * t34 + t40 * t56;
t21 = pkin(8) + t23;
t50 = -t15 * t21 + t16 * t20;
t32 = t39 * pkin(3) + pkin(8);
t49 = -t15 * t32 + t16 * t33;
t19 = t25 * pkin(3) + t35;
t38 = t44 ^ 2;
t37 = t41 ^ 2;
t31 = 0.2e1 * t53;
t14 = t16 ^ 2;
t13 = t44 * t15;
t11 = t41 * t52;
t8 = (-t37 + t38) * t16;
t7 = t40 * t10 + t39 * t48;
t4 = t15 * pkin(4) - t16 * pkin(8) + t19;
t3 = t5 * t41;
t2 = t41 * t4 + t44 * t7;
t1 = t44 * t4 - t41 * t7;
t6 = [1, 0, 0, t43 ^ 2, t43 * t58, 0, 0, 0, pkin(1) * t58, -0.2e1 * pkin(1) * t43, t26 ^ 2, -0.2e1 * t26 * t25, 0, 0, 0, t25 * t61, t26 * t61, -0.2e1 * t7 * t15 + 0.2e1 * t5 * t16, t19 ^ 2 + t5 ^ 2 + t7 ^ 2, t38 * t14, -0.2e1 * t14 * t53, 0.2e1 * t15 * t52, -0.2e1 * t15 * t54, t15 ^ 2, 0.2e1 * t1 * t15 + 0.2e1 * t5 * t54, -0.2e1 * t2 * t15 + 0.2e1 * t5 * t52; 0, 0, 0, 0, 0, t43, t46, 0, -t43 * pkin(6), -t46 * pkin(6), 0, 0, t26, -t25, 0, t17, t18, -t23 * t15 - t22 * t16, -t5 * t22 + t7 * t23, t11, t8, t12, t13, 0, t50 * t41 - t55, t50 * t44 + t3; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t36, -0.2e1 * t56, 0, t22 ^ 2 + t23 ^ 2, t37, t31, 0, 0, 0, t20 * t59, t20 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t17, t18, (-t15 * t39 - t16 * t40) * pkin(3), (t39 * t7 - t40 * t5) * pkin(3), t11, t8, t12, t13, 0, t49 * t41 - t55, t49 * t44 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t36, -t56, 0, (t22 * t40 + t23 * t39) * pkin(3), t37, t31, 0, 0, 0, -t51 * t44, t51 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t39 ^ 2 + t40 ^ 2) * pkin(3) ^ 2, t37, t31, 0, 0, 0, t33 * t59, t33 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, t13, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t54, t15, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t44, 0, -t41 * t21, -t44 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t44, 0, -t41 * t32, -t44 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
