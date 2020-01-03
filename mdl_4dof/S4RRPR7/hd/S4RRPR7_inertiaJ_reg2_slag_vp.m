% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t23 = sin(pkin(7));
t24 = cos(pkin(7));
t26 = sin(qJ(2));
t28 = cos(qJ(2));
t11 = t23 * t28 + t24 * t26;
t50 = -0.2e1 * t11;
t39 = -qJ(3) - pkin(5);
t13 = t39 * t28;
t34 = t39 * t26;
t4 = -t23 * t13 - t24 * t34;
t49 = t4 ^ 2;
t9 = t23 * t26 - t24 * t28;
t48 = t9 ^ 2;
t18 = -t28 * pkin(2) - pkin(1);
t47 = 0.2e1 * t18;
t46 = 0.2e1 * t28;
t45 = t23 * pkin(2);
t44 = t24 * pkin(2);
t25 = sin(qJ(4));
t43 = t25 * t9;
t42 = t25 * t11;
t27 = cos(qJ(4));
t41 = t25 * t27;
t40 = t27 * t11;
t19 = t25 ^ 2;
t21 = t27 ^ 2;
t38 = t19 + t21;
t20 = t26 ^ 2;
t22 = t28 ^ 2;
t37 = t20 + t22;
t36 = t9 * t50;
t35 = t25 * t40;
t3 = t9 * pkin(3) - t11 * pkin(6) + t18;
t6 = -t24 * t13 + t23 * t34;
t1 = -t25 * t6 + t27 * t3;
t2 = t25 * t3 + t27 * t6;
t33 = t1 * t27 + t2 * t25;
t32 = -t1 * t25 + t2 * t27;
t16 = pkin(6) + t45;
t17 = -pkin(3) - t44;
t31 = t11 * t17 - t16 * t9;
t8 = t11 ^ 2;
t7 = t27 * t9;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t20, t26 * t46, 0, t22, 0, 0, pkin(1) * t46, -0.2e1 * pkin(1) * t26, 0.2e1 * t37 * pkin(5), t37 * pkin(5) ^ 2 + pkin(1) ^ 2, t8, t36, 0, t48, 0, 0, t9 * t47, t11 * t47, 0.2e1 * t4 * t11 - 0.2e1 * t6 * t9, t18 ^ 2 + t6 ^ 2 + t49, t21 * t8, -0.2e1 * t8 * t41, 0.2e1 * t9 * t40, t19 * t8, t25 * t36, t48, 0.2e1 * t1 * t9 + 0.2e1 * t4 * t42, -0.2e1 * t2 * t9 + 0.2e1 * t4 * t40, t33 * t50, t1 ^ 2 + t2 ^ 2 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t28, 0, -t26 * pkin(5), -t28 * pkin(5), 0, 0, 0, 0, t11, 0, -t9, 0, -t4, -t6, (-t11 * t24 - t23 * t9) * pkin(2), (t23 * t6 - t24 * t4) * pkin(2), t35, (-t19 + t21) * t11, t43, -t35, t7, 0, t31 * t25 - t4 * t27, t4 * t25 + t31 * t27, t32, t32 * t16 + t4 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t44, -0.2e1 * t45, 0, (t23 ^ 2 + t24 ^ 2) * pkin(2) ^ 2, t19, 0.2e1 * t41, 0, t21, 0, 0, -0.2e1 * t17 * t27, 0.2e1 * t17 * t25, 0.2e1 * t38 * t16, t38 * t16 ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t11, 0, t18, 0, 0, 0, 0, 0, 0, t7, -t43, -t38 * t11, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t42, t9, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t27, 0, -t25 * t16, -t27 * t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
