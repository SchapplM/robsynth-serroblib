% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t34 = cos(qJ(4));
t23 = -pkin(4) * t34 - pkin(3);
t55 = 0.2e1 * t23;
t35 = cos(qJ(3));
t54 = -0.2e1 * t35;
t53 = 0.2e1 * t35;
t52 = pkin(7) + pkin(8);
t51 = pkin(3) * t34;
t30 = sin(qJ(5));
t50 = t30 * pkin(4);
t33 = cos(qJ(5));
t49 = t33 * pkin(4);
t29 = cos(pkin(9));
t21 = -pkin(1) * t29 - pkin(2);
t32 = sin(qJ(3));
t14 = -pkin(3) * t35 - pkin(7) * t32 + t21;
t31 = sin(qJ(4));
t28 = sin(pkin(9));
t20 = pkin(1) * t28 + pkin(6);
t39 = t35 * t20;
t37 = t34 * t39;
t5 = t37 + (-pkin(8) * t32 + t14) * t31;
t48 = t33 * t5;
t16 = t30 * t34 + t31 * t33;
t47 = t16 * t35;
t46 = t20 * t31;
t45 = t31 * t32;
t44 = t31 * t34;
t43 = t31 * t35;
t42 = t34 * t32;
t41 = t34 * t35;
t15 = t30 * t31 - t33 * t34;
t40 = t35 * t15;
t38 = t32 * t53;
t12 = t34 * t14;
t4 = -pkin(8) * t42 + t12 + (-pkin(4) - t46) * t35;
t1 = -t30 * t5 + t33 * t4;
t27 = t35 ^ 2;
t26 = t34 ^ 2;
t25 = t32 ^ 2;
t24 = t31 ^ 2;
t18 = t52 * t34;
t17 = t52 * t31;
t13 = (pkin(4) * t31 + t20) * t32;
t11 = -t30 * t45 + t33 * t42;
t10 = t16 * t32;
t9 = -t17 * t30 + t18 * t33;
t8 = -t17 * t33 - t18 * t30;
t7 = t14 * t31 + t37;
t6 = -t31 * t39 + t12;
t2 = t30 * t4 + t48;
t3 = [1, 0, 0, (t28 ^ 2 + t29 ^ 2) * pkin(1) ^ 2, t25, t38, 0, 0, 0, t21 * t54, 0.2e1 * t21 * t32, t26 * t25, -0.2e1 * t25 * t44, -0.2e1 * t32 * t41, t31 * t38, t27, 0.2e1 * t25 * t46 - 0.2e1 * t35 * t6, 0.2e1 * t20 * t25 * t34 + 0.2e1 * t35 * t7, t11 ^ 2, -0.2e1 * t11 * t10, t11 * t54, t10 * t53, t27, -0.2e1 * t1 * t35 + 0.2e1 * t10 * t13, 0.2e1 * t11 * t13 + 0.2e1 * t2 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t32, t35, 0, -t32 * t20, -t39, t31 * t42, (-t24 + t26) * t32, -t43, -t41, 0, -t20 * t42 + (-pkin(3) * t32 + pkin(7) * t35) * t31, pkin(7) * t41 + (t46 - t51) * t32, t11 * t16, -t10 * t16 - t11 * t15, -t47, t40, 0, t10 * t23 + t13 * t15 - t35 * t8, t11 * t23 + t13 * t16 + t35 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t32, 0, 0, 0, 0, 0, t41, -t43, 0, 0, 0, 0, 0, -t40, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t24, 0.2e1 * t44, 0, 0, 0, 0.2e1 * t51, -0.2e1 * pkin(3) * t31, t16 ^ 2, -0.2e1 * t16 * t15, 0, 0, 0, t15 * t55, t16 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t45, -t35, t6, -t7, 0, 0, t11, -t10, -t35, -t35 * t49 + t1, -t48 + (pkin(4) * t35 - t4) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t42, 0, 0, 0, 0, 0, -t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t34, 0, -t31 * pkin(7), -t34 * pkin(7), 0, 0, t16, -t15, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t49, -0.2e1 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, -t35, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t49, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
