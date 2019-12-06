% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t24 = cos(pkin(9));
t38 = t24 * pkin(2);
t18 = pkin(3) + t38;
t26 = sin(qJ(4));
t29 = cos(qJ(4));
t23 = sin(pkin(9));
t39 = t23 * pkin(2);
t13 = t26 * t18 + t29 * t39;
t11 = pkin(7) + t13;
t25 = sin(qJ(5));
t21 = t25 ^ 2;
t28 = cos(qJ(5));
t22 = t28 ^ 2;
t34 = t21 + t22;
t43 = t34 * t11;
t27 = sin(qJ(2));
t30 = cos(qJ(2));
t15 = -t23 * t27 + t24 * t30;
t16 = t23 * t30 + t24 * t27;
t4 = -t29 * t15 + t26 * t16;
t42 = t4 ^ 2;
t41 = 0.2e1 * t28;
t37 = t4 * t28;
t12 = t29 * t18 - t26 * t39;
t10 = -pkin(4) - t12;
t36 = pkin(4) - t10;
t35 = pkin(7) * t34;
t6 = t26 * t15 + t29 * t16;
t1 = t34 * t6;
t17 = t25 * t41;
t3 = t6 ^ 2;
t2 = t4 * t25;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 ^ 2 + t30 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t16 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 + t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t3 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t27, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, (t15 * t24 + t16 * t23) * pkin(2), 0, 0, 0, 0, 0, 0, -t4, -t6, 0, -t4 * t12 + t6 * t13, 0, 0, 0, 0, 0, 0, -t37, t2, t1, t4 * t10 + t43 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t39, 0, (t23 ^ 2 + t24 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t12, -0.2e1 * t13, 0, t12 ^ 2 + t13 ^ 2, t21, t17, 0, t22, 0, 0, -0.2e1 * t10 * t28, 0.2e1 * t10 * t25, 0.2e1 * t43, t34 * t11 ^ 2 + t10 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t6, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t2, t1, -t4 * pkin(4) + pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t12, -t13, 0, 0, t21, t17, 0, t22, 0, 0, t36 * t28, -t36 * t25, t35 + t43, -t10 * pkin(4) + pkin(7) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t21, t17, 0, t22, 0, 0, pkin(4) * t41, -0.2e1 * pkin(4) * t25, 0.2e1 * t35, t34 * pkin(7) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25 * t6, -t28 * t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t28, 0, -t25 * t11, -t28 * t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t28, 0, -t25 * pkin(7), -t28 * pkin(7), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
