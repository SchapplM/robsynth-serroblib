% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t32 = sin(qJ(4));
t33 = sin(qJ(3));
t35 = cos(qJ(4));
t36 = cos(qJ(3));
t15 = t32 * t33 - t35 * t36;
t24 = -t36 * pkin(3) - pkin(2);
t9 = t15 * pkin(4) + t24;
t51 = 0.2e1 * t9;
t50 = 0.2e1 * t24;
t49 = 0.2e1 * t36;
t48 = -pkin(7) - pkin(6);
t47 = t32 * pkin(3);
t27 = t35 * pkin(3);
t34 = sin(qJ(2));
t46 = t33 * t34;
t45 = t36 * t34;
t37 = cos(qJ(2));
t44 = t37 * t15;
t17 = t32 * t36 + t35 * t33;
t43 = t37 * t17;
t28 = t33 ^ 2;
t30 = t36 ^ 2;
t42 = t28 + t30;
t41 = t42 * t34;
t19 = t48 * t33;
t20 = t48 * t36;
t6 = t35 * t19 + t32 * t20;
t7 = t32 * t19 - t35 * t20;
t40 = pkin(3) ^ 2;
t38 = 0.2e1 * pkin(4);
t31 = t37 ^ 2;
t29 = t34 ^ 2;
t26 = t32 ^ 2 * t40;
t25 = -0.2e1 * t47;
t23 = t27 + pkin(4);
t14 = t17 ^ 2;
t13 = t15 ^ 2;
t12 = t15 * t47;
t11 = -t32 * t46 + t35 * t45;
t10 = t17 * t34;
t8 = t11 * t47;
t5 = -0.2e1 * t17 * t15;
t4 = -t15 * qJ(5) + t7;
t3 = -t17 * qJ(5) + t6;
t2 = t10 ^ 2 + t11 ^ 2 + t31;
t1 = t10 * t17 - t11 * t15;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 + t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t29 + t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t34, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t36, -t37 * t33, t41, t37 * pkin(2) + pkin(6) * t41, 0, 0, 0, 0, 0, 0, -t44, -t43, t1, -t10 * t6 + t11 * t7 - t37 * t24, 0, 0, 0, 0, 0, 0, -t44, -t43, t1, -t10 * t3 + t11 * t4 - t37 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t28, t33 * t49, 0, t30, 0, 0, pkin(2) * t49, -0.2e1 * pkin(2) * t33, 0.2e1 * t42 * pkin(6), t42 * pkin(6) ^ 2 + pkin(2) ^ 2, t14, t5, 0, t13, 0, 0, t15 * t50, t17 * t50, -0.2e1 * t7 * t15 - 0.2e1 * t6 * t17, t24 ^ 2 + t6 ^ 2 + t7 ^ 2, t14, t5, 0, t13, 0, 0, t15 * t51, t17 * t51, -0.2e1 * t4 * t15 - 0.2e1 * t3 * t17, t3 ^ 2 + t4 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t45, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, -t10 * t27 + t8, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, -t10 * t23 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, t36, 0, -t33 * pkin(6), -t36 * pkin(6), 0, 0, 0, 0, t17, 0, -t15, 0, t6, -t7, -t17 * t27 - t12, (t32 * t7 + t35 * t6) * pkin(3), 0, 0, t17, 0, -t15, 0, t3, -t4, -t23 * t17 - t12, t3 * t23 + t4 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t27, t25, 0, t35 ^ 2 * t40 + t26, 0, 0, 0, 0, 0, 1, 0.2e1 * t23, t25, 0, t23 ^ 2 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, -t10 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t15, 0, t6, -t7, 0, 0, 0, 0, t17, 0, -t15, 0, t3, -t4, -t17 * pkin(4), t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, -t47, 0, 0, 0, 0, 0, 0, 0, 1, t38 + t27, -t47, 0, t23 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t17, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t16;
