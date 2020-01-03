% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t16 = (pkin(1) + qJ(3));
t44 = (t16 ^ 2);
t43 = 2 * t16;
t20 = cos(qJ(5));
t42 = 0.2e1 * t20;
t21 = cos(qJ(4));
t41 = 0.2e1 * t21;
t40 = t21 * pkin(4);
t18 = sin(qJ(5));
t11 = t18 ^ 2;
t39 = t11 * t21;
t14 = t21 ^ 2;
t15 = -pkin(6) + qJ(2);
t38 = t14 * t15;
t37 = t18 * t20;
t36 = t18 * t21;
t19 = sin(qJ(4));
t35 = t19 * t15;
t34 = t20 * t19;
t9 = t20 * t21;
t33 = t21 * t15;
t32 = t21 * t19;
t13 = t20 ^ 2;
t31 = t11 + t13;
t12 = t19 ^ 2;
t5 = t12 + t14;
t30 = -0.2e1 * t32;
t29 = t18 * t9;
t28 = t31 * t19;
t27 = -pkin(7) * t19 - t40;
t4 = t19 * pkin(4) - t21 * pkin(7) + t16;
t1 = -t18 * t35 + t20 * t4;
t2 = t15 * t34 + t18 * t4;
t26 = -t1 * t20 - t2 * t18;
t25 = -t1 * t18 + t2 * t20;
t23 = (qJ(2) ^ 2);
t22 = 2 * qJ(2);
t10 = t15 ^ 2;
t8 = t13 * t21;
t7 = t18 * t19;
t6 = t14 * t10;
t3 = t5 * t15;
t17 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t22, pkin(1) ^ 2 + t23, 1, 0, 0, 0, 0, 0, 0, t22, t43, t23 + t44, t14, t30, 0, t12, 0, 0, t19 * t43, t16 * t41, -0.2e1 * t3, t12 * t10 + t44 + t6, t13 * t14, -0.2e1 * t14 * t37, t32 * t42, t11 * t14, t18 * t30, t12, 0.2e1 * t1 * t19 - 0.2e1 * t18 * t38, -0.2e1 * t2 * t19 - 0.2e1 * t20 * t38, t26 * t41, t1 ^ 2 + t2 ^ 2 + t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -1, -t16, 0, 0, 0, 0, 0, 0, -t19, -t21, 0, -t16, 0, 0, 0, 0, 0, 0, -t34, t7, t8 + t39, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(2), 0, 0, 0, 0, 0, 0, 0, 0, -t5, t3, 0, 0, 0, 0, 0, 0, -t5 * t18, -t5 * t20, 0, t25 * t19 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * t12 + t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, 0, t33, -t35, 0, 0, t29, t8 - t39, t7, -t29, t34, 0, t27 * t18 + t20 * t33, -t18 * t33 + t27 * t20, t25, pkin(4) * t33 + t25 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t36, t28, pkin(7) * t28 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t11, 0.2e1 * t37, 0, t13, 0, 0, pkin(4) * t42, -0.2e1 * pkin(4) * t18, 0.2e1 * t31 * pkin(7), t31 * pkin(7) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t36, t19, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t20, 0, -t18 * pkin(7), -t20 * pkin(7), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t17;
