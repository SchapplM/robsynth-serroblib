% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP2
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t26 = sin(qJ(4));
t22 = t26 ^ 2;
t28 = cos(qJ(4));
t23 = t28 ^ 2;
t15 = t22 + t23;
t25 = cos(pkin(8));
t38 = t25 * pkin(1);
t18 = pkin(2) + t38;
t27 = sin(qJ(3));
t29 = cos(qJ(3));
t24 = sin(pkin(8));
t39 = t24 * pkin(1);
t11 = t27 * t18 + t29 * t39;
t9 = pkin(7) + t11;
t40 = t15 * t9;
t44 = 0.2e1 * t26;
t43 = -0.2e1 * t28;
t42 = 0.2e1 * t28;
t10 = t29 * t18 - t27 * t39;
t8 = -pkin(3) - t10;
t41 = pkin(3) - t8;
t37 = t26 * pkin(4);
t36 = t28 * pkin(4);
t19 = -pkin(3) - t36;
t4 = t8 - t36;
t35 = t19 + t4;
t34 = -qJ(5) - pkin(7);
t33 = qJ(5) + t9;
t32 = t15 * pkin(7);
t17 = t26 * t42;
t14 = t34 * t28;
t13 = t34 * t26;
t12 = t14 * t28;
t3 = t33 * t28;
t2 = t33 * t26;
t1 = t3 * t28;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t39, 0, (t24 ^ 2 + t25 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t10, -0.2e1 * t11, 0, t10 ^ 2 + t11 ^ 2, t22, t17, 0, t23, 0, 0, t8 * t43, t8 * t44, 0.2e1 * t40, t15 * t9 ^ 2 + t8 ^ 2, t22, t17, 0, t23, 0, 0, t4 * t43, t4 * t44, 0.2e1 * t2 * t26 + 0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t28 + t3 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t10, -t11, 0, 0, t22, t17, 0, t23, 0, 0, t41 * t28, -t41 * t26, t32 + t40, -t8 * pkin(3) + pkin(7) * t40, t22, t17, 0, t23, 0, 0, -t35 * t28, t35 * t26, t1 - t12 + (-t13 + t2) * t26, -t2 * t13 - t3 * t14 + t4 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t13 - t26 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t22, t17, 0, t23, 0, 0, pkin(3) * t42, -0.2e1 * pkin(3) * t26, 0.2e1 * t32, pkin(7) ^ 2 * t15 + pkin(3) ^ 2, t22, t17, 0, t23, 0, 0, t19 * t43, t19 * t44, -0.2e1 * t13 * t26 - 0.2e1 * t12, t13 ^ 2 + t14 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t28, 0, -t26 * t9, -t28 * t9, 0, 0, 0, 0, t26, 0, t28, 0, -t2, -t3, -t37, -t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t28, 0, -t26 * pkin(7), -t28 * pkin(7), 0, 0, 0, 0, t26, 0, t28, 0, t13, t14, -t37, t13 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t26, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t26, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t5;
