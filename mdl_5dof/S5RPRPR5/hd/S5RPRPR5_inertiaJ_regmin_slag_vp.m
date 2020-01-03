% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t36 = sin(pkin(9));
t38 = cos(pkin(9));
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t27 = t36 * t43 + t38 * t41;
t37 = sin(pkin(8));
t21 = t27 * t37;
t52 = t43 * t37;
t54 = t41 * t37;
t22 = -t36 * t54 + t38 * t52;
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t9 = -t40 * t21 + t42 * t22;
t58 = -0.2e1 * t9;
t57 = -0.2e1 * t37;
t39 = cos(pkin(8));
t56 = 0.2e1 * t39;
t55 = pkin(3) * t36;
t53 = t41 * t39;
t51 = t43 * t39;
t29 = -t39 * pkin(2) - t37 * pkin(6) - pkin(1);
t25 = t43 * t29;
t48 = qJ(4) * t37;
t49 = qJ(2) * t41;
t12 = -t43 * t48 + t25 + (-pkin(3) - t49) * t39;
t46 = qJ(2) * t51;
t17 = t46 + (t29 - t48) * t41;
t7 = t36 * t12 + t38 * t17;
t28 = pkin(3) * t54 + t37 * qJ(2);
t34 = t37 ^ 2;
t35 = t39 ^ 2;
t50 = t34 + t35;
t47 = t34 * qJ(2);
t6 = t38 * t12 - t36 * t17;
t4 = -t39 * pkin(4) - t22 * pkin(7) + t6;
t5 = -t21 * pkin(7) + t7;
t1 = t42 * t4 - t40 * t5;
t2 = t40 * t4 + t42 * t5;
t32 = t38 * pkin(3) + pkin(4);
t26 = -t36 * t41 + t38 * t43;
t24 = t40 * t32 + t42 * t55;
t23 = t42 * t32 - t40 * t55;
t20 = t41 * t29 + t46;
t19 = -t39 * t49 + t25;
t15 = t21 * pkin(4) + t28;
t14 = t40 * t26 + t42 * t27;
t13 = t42 * t26 - t40 * t27;
t8 = t42 * t21 + t40 * t22;
t3 = [1, 0, 0, pkin(1) * t56, pkin(1) * t57, 0.2e1 * t50 * qJ(2), t50 * qJ(2) ^ 2 + pkin(1) ^ 2, t43 ^ 2 * t34, -0.2e1 * t43 * t34 * t41, t51 * t57, 0.2e1 * t37 * t53, t35, -0.2e1 * t19 * t39 + 0.2e1 * t41 * t47, 0.2e1 * t20 * t39 + 0.2e1 * t43 * t47, -0.2e1 * t7 * t21 - 0.2e1 * t6 * t22, t28 ^ 2 + t6 ^ 2 + t7 ^ 2, t9 ^ 2, t8 * t58, t39 * t58, t8 * t56, t35, -0.2e1 * t1 * t39 + 0.2e1 * t15 * t8, 0.2e1 * t15 * t9 + 0.2e1 * t2 * t39; 0, 0, 0, -t39, t37, 0, -pkin(1), 0, 0, 0, 0, 0, -t51, t53, -t27 * t21 - t26 * t22, t6 * t26 + t7 * t27, 0, 0, 0, 0, 0, -t13 * t39, t14 * t39; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t26 ^ 2 + t27 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t54, -t39, t19, -t20, (-t21 * t36 - t22 * t38) * pkin(3), (t36 * t7 + t38 * t6) * pkin(3), 0, 0, t9, -t8, -t39, -t23 * t39 + t1, t24 * t39 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t41, 0, (t26 * t38 + t27 * t36) * pkin(3), 0, 0, 0, 0, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t36 ^ 2 + t38 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t23, -0.2e1 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, t8, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, -t39, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t23, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
