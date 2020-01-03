% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t27 = sin(pkin(8));
t28 = cos(pkin(8));
t43 = sin(qJ(3));
t44 = cos(qJ(3));
t17 = t27 * t43 - t28 * t44;
t19 = t27 * t44 + t28 * t43;
t29 = sin(qJ(4));
t30 = cos(qJ(4));
t10 = t30 * t17 + t29 * t19;
t49 = t10 ^ 2;
t21 = -t28 * pkin(2) - pkin(1);
t15 = t17 * pkin(3) + t21;
t48 = 0.2e1 * t15;
t47 = 0.2e1 * t21;
t46 = 0.2e1 * t28;
t45 = t30 * pkin(3);
t12 = -t29 * t17 + t30 * t19;
t42 = t12 * t10;
t41 = pkin(6) + qJ(2);
t25 = t27 ^ 2;
t26 = t28 ^ 2;
t40 = t25 + t26;
t37 = t41 * t28;
t38 = t41 * t27;
t13 = -t37 * t43 - t38 * t44;
t35 = -t19 * pkin(7) + t13;
t14 = t37 * t44 - t38 * t43;
t8 = -t17 * pkin(7) + t14;
t4 = t29 * t8 - t30 * t35;
t6 = t29 * t35 + t30 * t8;
t39 = t4 ^ 2 + t6 ^ 2;
t36 = -0.2e1 * t6 * t10 + 0.2e1 * t4 * t12;
t33 = 2 * pkin(4);
t31 = 2 * qJ(5);
t24 = t29 * pkin(3);
t22 = pkin(4) + t45;
t20 = t24 + qJ(5);
t9 = t12 ^ 2;
t3 = t10 * pkin(4) - t12 * qJ(5) + t15;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t25, t27 * t46, 0, t26, 0, 0, pkin(1) * t46, -0.2e1 * pkin(1) * t27, 0.2e1 * t40 * qJ(2), qJ(2) ^ 2 * t40 + pkin(1) ^ 2, t19 ^ 2, -0.2e1 * t19 * t17, 0, t17 ^ 2, 0, 0, t17 * t47, t19 * t47, -0.2e1 * t13 * t19 - 0.2e1 * t14 * t17, t13 ^ 2 + t14 ^ 2 + t21 ^ 2, t9, -0.2e1 * t42, 0, t49, 0, 0, t10 * t48, t12 * t48, t36, t15 ^ 2 + t39, t9, 0, 0.2e1 * t42, 0, 0, t49, 0.2e1 * t3 * t10, t36, -0.2e1 * t3 * t12, t3 ^ 2 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t27, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t17, t19, 0, t21, 0, 0, 0, 0, 0, 0, t10, t12, 0, t15, 0, 0, 0, 0, 0, 0, t10, 0, -t12, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, 0, t13, -t14, 0, 0, 0, 0, t12, 0, -t10, 0, -t4, -t6, (-t10 * t29 - t12 * t30) * pkin(3), (t29 * t6 - t30 * t4) * pkin(3), 0, t12, 0, 0, t10, 0, -t4, -t20 * t10 - t22 * t12, t6, t6 * t20 - t4 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t45, -0.2e1 * t24, 0, (t29 ^ 2 + t30 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t22, 0, 0.2e1 * t20, t20 ^ 2 + t22 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, 0, -t4, -t6, 0, 0, 0, t12, 0, 0, t10, 0, -t4, -pkin(4) * t12 - t10 * qJ(5), t6, -t4 * pkin(4) + t6 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t45, -t24, 0, 0, 0, 0, 0, 1, 0, 0, t33 + t45, 0, t31 + t24, t22 * pkin(4) + t20 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t33, 0, t31, pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
