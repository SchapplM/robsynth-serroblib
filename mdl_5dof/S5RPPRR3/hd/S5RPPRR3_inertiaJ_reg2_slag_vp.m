% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRR3
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
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t28 = sin(pkin(9));
t30 = cos(pkin(9));
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t19 = t33 * t28 - t35 * t30;
t31 = cos(pkin(8));
t42 = t31 * pkin(1);
t25 = -pkin(2) - t42;
t22 = -t30 * pkin(3) + t25;
t12 = t19 * pkin(4) + t22;
t46 = 0.2e1 * t12;
t45 = 0.2e1 * t22;
t44 = 0.2e1 * t28;
t29 = sin(pkin(8));
t43 = t29 * pkin(1);
t32 = sin(qJ(5));
t41 = t32 * pkin(4);
t34 = cos(qJ(5));
t40 = t34 * pkin(4);
t24 = qJ(3) + t43;
t39 = pkin(6) + t24;
t26 = t28 ^ 2;
t27 = t30 ^ 2;
t38 = t26 + t27;
t15 = t39 * t28;
t16 = t39 * t30;
t5 = -t35 * t15 - t33 * t16;
t6 = -t33 * t15 + t35 * t16;
t21 = t35 * t28 + t33 * t30;
t18 = t21 ^ 2;
t17 = t19 ^ 2;
t11 = -t32 * t19 + t34 * t21;
t9 = t34 * t19 + t32 * t21;
t8 = t11 ^ 2;
t7 = t9 ^ 2;
t4 = -t19 * pkin(7) + t6;
t3 = -t21 * pkin(7) + t5;
t2 = t32 * t3 + t34 * t4;
t1 = t34 * t3 - t32 * t4;
t10 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t42, -0.2e1 * t43, 0, (t29 ^ 2 + t31 ^ 2) * pkin(1) ^ 2, t26, t30 * t44, 0, t27, 0, 0, -0.2e1 * t25 * t30, t25 * t44, 0.2e1 * t38 * t24, t38 * t24 ^ 2 + t25 ^ 2, t18, -0.2e1 * t21 * t19, 0, t17, 0, 0, t19 * t45, t21 * t45, -0.2e1 * t6 * t19 - 0.2e1 * t5 * t21, t22 ^ 2 + t5 ^ 2 + t6 ^ 2, t8, -0.2e1 * t11 * t9, 0, t7, 0, 0, t9 * t46, t11 * t46, -0.2e1 * t1 * t11 - 0.2e1 * t2 * t9, t1 ^ 2 + t12 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t19 + t6 * t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t9 + t2 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 + t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t28, 0, t25, 0, 0, 0, 0, 0, 0, t19, t21, 0, t22, 0, 0, 0, 0, 0, 0, t9, t11, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, 0, t5, -t6, 0, 0, 0, 0, t11, 0, -t9, 0, t1, -t2, (-t11 * t34 - t32 * t9) * pkin(4), (t1 * t34 + t2 * t32) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t21, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t11, 0, (t11 * t32 - t34 * t9) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t40, -0.2e1 * t41, 0, (t32 ^ 2 + t34 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t40, -t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t10;
