% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t21 = sin(qJ(3));
t41 = -0.2e1 * t21;
t24 = cos(qJ(3));
t40 = 0.2e1 * t24;
t23 = cos(qJ(4));
t39 = pkin(3) * t23;
t20 = sin(qJ(4));
t38 = pkin(7) * t20;
t37 = t20 * pkin(4);
t18 = sin(pkin(5));
t36 = t18 * sin(qJ(2));
t35 = t18 * cos(qJ(2));
t34 = t20 * t21;
t33 = t20 * t23;
t32 = t20 * t24;
t31 = t23 * t21;
t30 = t23 * t24;
t29 = -qJ(5) - pkin(8);
t28 = qJ(5) * t21;
t27 = t21 * t40;
t26 = pkin(7) * t30;
t19 = cos(pkin(5));
t17 = t23 ^ 2;
t16 = t21 ^ 2;
t15 = t20 ^ 2;
t14 = -t23 * pkin(4) - pkin(3);
t13 = t29 * t23;
t12 = t29 * t20;
t11 = -t24 * pkin(3) - t21 * pkin(8) - pkin(2);
t10 = (pkin(7) + t37) * t21;
t9 = t23 * t11;
t8 = t19 * t21 + t24 * t36;
t7 = -t19 * t24 + t21 * t36;
t6 = t20 * t11 + t26;
t5 = -pkin(7) * t32 + t9;
t4 = t26 + (t11 - t28) * t20;
t3 = -t20 * t35 + t8 * t23;
t2 = -t8 * t20 - t23 * t35;
t1 = -t23 * t28 + t9 + (-pkin(4) - t38) * t24;
t22 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 ^ 2 + t3 ^ 2 + t7 ^ 2; 0, 0, t35, -t36, 0, 0, 0, 0, 0, t24 * t35, -t21 * t35, 0, 0, 0, 0, 0, -t2 * t24 + t7 * t34, t3 * t24 + t7 * t31, (-t2 * t23 - t20 * t3) * t21, t2 * t1 + t7 * t10 + t3 * t4; 0, 1, 0, 0, t16, t27, 0, 0, 0, pkin(2) * t40, pkin(2) * t41, t17 * t16, -0.2e1 * t16 * t33, t30 * t41, t20 * t27, t24 ^ 2, 0.2e1 * t16 * t38 - 0.2e1 * t5 * t24, 0.2e1 * t16 * pkin(7) * t23 + 0.2e1 * t6 * t24, 0.2e1 * (-t1 * t23 - t20 * t4) * t21, t1 ^ 2 + t10 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, 0, 0, 0, 0, -t7 * t23, t7 * t20, -t2 * t20 + t3 * t23, t2 * t12 - t3 * t13 + t7 * t14; 0, 0, 0, 0, 0, 0, t21, t24, 0, -t21 * pkin(7), -t24 * pkin(7), t20 * t31, (-t15 + t17) * t21, -t32, -t30, 0, -pkin(7) * t31 + (-pkin(3) * t21 + pkin(8) * t24) * t20, pkin(8) * t30 + (t38 - t39) * t21, (-t12 * t21 + t4) * t23 + (t13 * t21 - t1) * t20, t1 * t12 + t10 * t14 - t4 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t15, 0.2e1 * t33, 0, 0, 0, 0.2e1 * t39, -0.2e1 * pkin(3) * t20, -0.2e1 * t12 * t20 - 0.2e1 * t13 * t23, t12 ^ 2 + t13 ^ 2 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t3, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t34, -t24, t5, -t6, -pkin(4) * t31, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t23, 0, -t20 * pkin(8), -t23 * pkin(8), -t37, t12 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t22;
