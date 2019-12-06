% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t31 = sin(pkin(9));
t32 = cos(pkin(9));
t34 = sin(qJ(3));
t37 = cos(qJ(3));
t18 = t31 * t34 - t32 * t37;
t26 = -t37 * pkin(3) - pkin(2);
t12 = t18 * pkin(4) + t26;
t48 = 0.2e1 * t12;
t47 = 0.2e1 * t26;
t46 = 0.2e1 * t37;
t45 = t31 * pkin(3);
t44 = t32 * pkin(3);
t43 = -qJ(4) - pkin(6);
t27 = t34 ^ 2;
t29 = t37 ^ 2;
t42 = t27 + t29;
t35 = sin(qJ(2));
t41 = t42 * t35;
t22 = t43 * t34;
t23 = t43 * t37;
t10 = t32 * t22 + t31 * t23;
t11 = t31 * t22 - t32 * t23;
t20 = t31 * t37 + t32 * t34;
t38 = cos(qJ(2));
t36 = cos(qJ(5));
t33 = sin(qJ(5));
t30 = t38 ^ 2;
t28 = t35 ^ 2;
t25 = pkin(4) + t44;
t16 = t33 * t25 + t36 * t45;
t15 = t36 * t25 - t33 * t45;
t14 = t18 * t35;
t13 = t20 * t35;
t9 = -t33 * t18 + t36 * t20;
t7 = t36 * t18 + t33 * t20;
t6 = -t18 * pkin(7) + t11;
t5 = -t20 * pkin(7) + t10;
t4 = -t33 * t13 - t36 * t14;
t3 = -t36 * t13 + t33 * t14;
t2 = t33 * t5 + t36 * t6;
t1 = -t33 * t6 + t36 * t5;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 + t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t28 + t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 ^ 2 + t14 ^ 2 + t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t37, -t38 * t34, t41, t38 * pkin(2) + pkin(6) * t41, 0, 0, 0, 0, 0, 0, -t38 * t18, -t38 * t20, t13 * t20 + t14 * t18, -t13 * t10 - t14 * t11 - t38 * t26, 0, 0, 0, 0, 0, 0, -t38 * t7, -t38 * t9, -t3 * t9 - t4 * t7, t3 * t1 - t38 * t12 + t4 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t27, t34 * t46, 0, t29, 0, 0, pkin(2) * t46, -0.2e1 * pkin(2) * t34, 0.2e1 * t42 * pkin(6), t42 * pkin(6) ^ 2 + pkin(2) ^ 2, t20 ^ 2, -0.2e1 * t20 * t18, 0, t18 ^ 2, 0, 0, t18 * t47, t20 * t47, -0.2e1 * t10 * t20 - 0.2e1 * t11 * t18, t10 ^ 2 + t11 ^ 2 + t26 ^ 2, t9 ^ 2, -0.2e1 * t9 * t7, 0, t7 ^ 2, 0, 0, t7 * t48, t9 * t48, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t7, t1 ^ 2 + t12 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34 * t35, -t37 * t35, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t14, 0, (-t13 * t32 - t14 * t31) * pkin(3), 0, 0, 0, 0, 0, 0, t3, -t4, 0, t3 * t15 + t4 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t37, 0, -t34 * pkin(6), -t37 * pkin(6), 0, 0, 0, 0, t20, 0, -t18, 0, t10, -t11, (-t18 * t31 - t20 * t32) * pkin(3), (t10 * t32 + t11 * t31) * pkin(3), 0, 0, t9, 0, -t7, 0, t1, -t2, -t15 * t9 - t16 * t7, t1 * t15 + t2 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t44, -0.2e1 * t45, 0, (t31 ^ 2 + t32 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t15, -0.2e1 * t16, 0, t15 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t20, 0, t26, 0, 0, 0, 0, 0, 0, t7, t9, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t15, -t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
