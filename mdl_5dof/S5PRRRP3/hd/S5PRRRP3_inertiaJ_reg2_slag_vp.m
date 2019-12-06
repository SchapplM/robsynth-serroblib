% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRP3
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t30 = cos(qJ(3));
t21 = -t30 * pkin(3) - pkin(2);
t27 = sin(qJ(4));
t28 = sin(qJ(3));
t29 = cos(qJ(4));
t13 = t27 * t28 - t29 * t30;
t36 = t13 * pkin(4);
t8 = t21 + t36;
t40 = 0.2e1 * t8;
t39 = 0.2e1 * t21;
t38 = 0.2e1 * t30;
t37 = -pkin(7) - pkin(6);
t35 = t27 * pkin(3);
t24 = t29 * pkin(3);
t25 = t28 ^ 2;
t26 = t30 ^ 2;
t34 = t25 + t26;
t17 = t37 * t28;
t18 = t37 * t30;
t6 = t29 * t17 + t27 * t18;
t7 = t27 * t17 - t29 * t18;
t33 = pkin(3) ^ 2;
t31 = 0.2e1 * pkin(4);
t23 = t27 ^ 2 * t33;
t22 = -0.2e1 * t35;
t20 = t24 + pkin(4);
t15 = t27 * t30 + t29 * t28;
t12 = t15 ^ 2;
t11 = t13 ^ 2;
t10 = t15 * t35;
t9 = t13 * t35;
t5 = -0.2e1 * t15 * t13;
t4 = t12 + t11;
t3 = -t13 * qJ(5) + t7;
t2 = -t15 * qJ(5) + t6;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t6 + t15 * t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t2 + t15 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t25, t28 * t38, 0, t26, 0, 0, pkin(2) * t38, -0.2e1 * pkin(2) * t28, 0.2e1 * t34 * pkin(6), t34 * pkin(6) ^ 2 + pkin(2) ^ 2, t12, t5, 0, t11, 0, 0, t13 * t39, t15 * t39, -0.2e1 * t7 * t13 - 0.2e1 * t6 * t15, t21 ^ 2 + t6 ^ 2 + t7 ^ 2, t12, t5, 0, t11, 0, 0, t13 * t40, t15 * t40, -0.2e1 * t3 * t13 - 0.2e1 * t2 * t15, t2 ^ 2 + t3 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t28, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, -t13 * t24 + t10, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, -t13 * t20 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t30, 0, -t28 * pkin(6), -t30 * pkin(6), 0, 0, 0, 0, t15, 0, -t13, 0, t6, -t7, -t15 * t24 - t9, (t27 * t7 + t29 * t6) * pkin(3), 0, 0, t15, 0, -t13, 0, t2, -t3, -t20 * t15 - t9, t2 * t20 + t3 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t24, t22, 0, t29 ^ 2 * t33 + t23, 0, 0, 0, 0, 0, 1, 0.2e1 * t20, t22, 0, t20 ^ 2 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t13, 0, t6, -t7, 0, 0, 0, 0, t15, 0, -t13, 0, t2, -t3, -t15 * pkin(4), t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t24, -t35, 0, 0, 0, 0, 0, 0, 0, 1, t31 + t24, -t35, 0, t20 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t31, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
