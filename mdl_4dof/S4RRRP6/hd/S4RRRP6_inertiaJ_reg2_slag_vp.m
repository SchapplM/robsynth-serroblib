% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t26 = sin(qJ(2));
t49 = -0.2e1 * t26;
t48 = 0.2e1 * t26;
t28 = cos(qJ(2));
t47 = 0.2e1 * t28;
t27 = cos(qJ(3));
t46 = pkin(2) * t27;
t25 = sin(qJ(3));
t45 = pkin(5) * t25;
t22 = t26 ^ 2;
t44 = t22 * pkin(5);
t43 = t25 * pkin(3);
t42 = t26 * pkin(5);
t41 = t25 * t26;
t40 = t25 * t27;
t39 = t25 * t28;
t18 = t27 * t26;
t38 = t27 * t28;
t37 = -qJ(4) - pkin(6);
t21 = t25 ^ 2;
t23 = t27 ^ 2;
t36 = t21 + t23;
t35 = qJ(4) * t26;
t34 = t26 * t47;
t33 = pkin(5) * t38;
t8 = -t28 * pkin(2) - t26 * pkin(6) - pkin(1);
t5 = t27 * t8;
t32 = -t27 * t35 + t5;
t3 = -pkin(5) * t39 + t5;
t4 = t25 * t8 + t33;
t31 = -t3 * t25 + t4 * t27;
t30 = pkin(5) ^ 2;
t24 = t28 ^ 2;
t20 = t22 * t30;
t19 = -t27 * pkin(3) - pkin(2);
t17 = t23 * t22;
t16 = t21 * t22;
t15 = 0.2e1 * t40;
t14 = t25 * t18;
t13 = t38 * t49;
t12 = -0.2e1 * t22 * t40;
t11 = t25 * t34;
t10 = t37 * t27;
t9 = t37 * t25;
t7 = (pkin(5) + t43) * t26;
t6 = (-t21 + t23) * t26;
t2 = t33 + (t8 - t35) * t25;
t1 = (-pkin(3) - t45) * t28 + t32;
t29 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t22, t34, 0, t24, 0, 0, pkin(1) * t47, pkin(1) * t49, 0.2e1 * (t22 + t24) * pkin(5), pkin(1) ^ 2 + t24 * t30 + t20, t17, t12, t13, t16, t11, t24, 0.2e1 * t25 * t44 - 0.2e1 * t3 * t28, 0.2e1 * t27 * t44 + 0.2e1 * t4 * t28, (-t25 * t4 - t27 * t3) * t48, t3 ^ 2 + t4 ^ 2 + t20, t17, t12, t13, t16, t11, t24, -0.2e1 * t1 * t28 + 0.2e1 * t7 * t41, 0.2e1 * t7 * t18 + 0.2e1 * t2 * t28, (-t1 * t27 - t2 * t25) * t48, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t28, 0, -t42, -t28 * pkin(5), 0, 0, t14, t6, -t39, -t14, -t38, 0, -pkin(5) * t18 + (-pkin(2) * t26 + pkin(6) * t28) * t25, pkin(6) * t38 + (t45 - t46) * t26, t31, -pkin(2) * t42 + t31 * pkin(6), t14, t6, -t39, -t14, -t38, 0, t19 * t41 - t7 * t27 - t9 * t28, -t10 * t28 + t19 * t18 + t7 * t25, (-t26 * t9 + t2) * t27 + (t10 * t26 - t1) * t25, t1 * t9 - t2 * t10 + t7 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t21, t15, 0, t23, 0, 0, 0.2e1 * t46, -0.2e1 * pkin(2) * t25, 0.2e1 * t36 * pkin(6), t36 * pkin(6) ^ 2 + pkin(2) ^ 2, t21, t15, 0, t23, 0, 0, -0.2e1 * t19 * t27, 0.2e1 * t19 * t25, -0.2e1 * t10 * t27 - 0.2e1 * t9 * t25, t10 ^ 2 + t19 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t41, -t28, t3, -t4, 0, 0, 0, 0, t18, 0, -t41, -t28, (-0.2e1 * pkin(3) - t45) * t28 + t32, -t2, -pkin(3) * t18, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t27, 0, -t25 * pkin(6), -t27 * pkin(6), 0, 0, 0, 0, t25, 0, t27, 0, t9, t10, -t43, t9 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(3), 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t18, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t25, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t29;
