% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t36 = sin(pkin(8));
t56 = -0.2e1 * t36;
t38 = cos(pkin(8));
t55 = -0.2e1 * t38;
t54 = 0.2e1 * t38;
t53 = t38 * pkin(4);
t39 = sin(qJ(5));
t52 = t39 * pkin(4);
t41 = cos(qJ(5));
t51 = t41 * pkin(4);
t35 = sin(pkin(9));
t37 = cos(pkin(9));
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t24 = t42 * t35 + t40 * t37;
t19 = t24 * t36;
t26 = -t38 * pkin(2) - t36 * qJ(3) - pkin(1);
t22 = t37 * t26;
t48 = t37 * t36;
t11 = -pkin(6) * t48 + t22 + (-qJ(2) * t35 - pkin(3)) * t38;
t46 = qJ(2) * t38;
t17 = t35 * t26 + t37 * t46;
t49 = t35 * t36;
t15 = -pkin(6) * t49 + t17;
t7 = t40 * t11 + t42 * t15;
t5 = -t19 * pkin(7) + t7;
t50 = t41 * t5;
t30 = t36 * qJ(2);
t25 = pkin(3) * t49 + t30;
t47 = t35 ^ 2 + t37 ^ 2;
t32 = t36 ^ 2;
t45 = t32 * qJ(2);
t23 = -t40 * t35 + t42 * t37;
t20 = t23 * t36;
t6 = t42 * t11 - t40 * t15;
t4 = -t20 * pkin(7) - t53 + t6;
t1 = -t39 * t5 + t41 * t4;
t16 = -t35 * t46 + t22;
t44 = t16 * t37 + t17 * t35;
t43 = qJ(2) ^ 2;
t34 = t38 ^ 2;
t29 = t32 * t43;
t14 = t39 * t23 + t41 * t24;
t13 = t41 * t23 - t39 * t24;
t12 = t19 * pkin(4) + t25;
t9 = -t39 * t19 + t41 * t20;
t8 = t41 * t19 + t39 * t20;
t2 = t39 * t4 + t50;
t3 = [1, 0, 0, pkin(1) * t54, pkin(1) * t56, 0.2e1 * (t32 + t34) * qJ(2), pkin(1) ^ 2 + t34 * t43 + t29, -0.2e1 * t16 * t38 + 0.2e1 * t35 * t45, 0.2e1 * t17 * t38 + 0.2e1 * t37 * t45, t44 * t56, t16 ^ 2 + t17 ^ 2 + t29, t20 ^ 2, -0.2e1 * t20 * t19, t20 * t55, t19 * t54, t34, 0.2e1 * t25 * t19 - 0.2e1 * t6 * t38, 0.2e1 * t25 * t20 + 0.2e1 * t7 * t38, t9 ^ 2, -0.2e1 * t9 * t8, t9 * t55, t8 * t54, t34, -0.2e1 * t1 * t38 + 0.2e1 * t12 * t8, 0.2e1 * t12 * t9 + 0.2e1 * t2 * t38; 0, 0, 0, -t38, t36, 0, -pkin(1), -t37 * t38, t35 * t38, -t47 * t36, t44, 0, 0, 0, 0, 0, -t23 * t38, t24 * t38, 0, 0, 0, 0, 0, -t13 * t38, t14 * t38; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t49, t48, 0, t30, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, t8, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, -t38, t6, -t7, 0, 0, t9, -t8, -t38, -t38 * t51 + t1, -t50 + (-t4 + t53) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t24, 0, 0, 0, 0, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t51, -0.2e1 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, -t38, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t51, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
