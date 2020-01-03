% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR6
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
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t25 = cos(pkin(8));
t17 = -t25 * pkin(1) - pkin(2);
t24 = cos(pkin(9));
t13 = -t24 * pkin(3) + t17;
t38 = 0.2e1 * t13;
t23 = sin(pkin(8));
t15 = t23 * pkin(1) + qJ(3);
t37 = pkin(6) + t15;
t36 = cos(qJ(4));
t22 = sin(pkin(9));
t27 = sin(qJ(4));
t11 = t27 * t22 - t36 * t24;
t26 = sin(qJ(5));
t6 = t26 * t11;
t12 = t36 * t22 + t27 * t24;
t35 = t26 * t12;
t28 = cos(qJ(5));
t34 = t26 * t28;
t33 = t28 * t12;
t32 = t22 ^ 2 + t24 ^ 2;
t31 = -0.2e1 * t12 * t11;
t30 = -pkin(4) * t12 - pkin(7) * t11;
t21 = t28 ^ 2;
t20 = t26 ^ 2;
t10 = t12 ^ 2;
t9 = t37 * t24;
t8 = t37 * t22;
t7 = t28 * t11;
t5 = -t27 * t8 + t36 * t9;
t4 = t27 * t9 + t36 * t8;
t3 = t11 * pkin(4) - t12 * pkin(7) + t13;
t2 = t26 * t3 + t28 * t5;
t1 = -t26 * t5 + t28 * t3;
t14 = [1, 0, 0, (t23 ^ 2 + t25 ^ 2) * pkin(1) ^ 2, -0.2e1 * t17 * t24, 0.2e1 * t17 * t22, 0.2e1 * t32 * t15, t32 * t15 ^ 2 + t17 ^ 2, t10, t31, 0, 0, 0, t11 * t38, t12 * t38, t21 * t10, -0.2e1 * t10 * t34, 0.2e1 * t11 * t33, t26 * t31, t11 ^ 2, 0.2e1 * t1 * t11 + 0.2e1 * t4 * t35, -0.2e1 * t2 * t11 + 0.2e1 * t4 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t24, t22, 0, t17, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, t7, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, -t4, -t5, t26 * t33, (-t20 + t21) * t12, t6, t7, 0, t30 * t26 - t4 * t28, t4 * t26 + t30 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, 0, 0, 0, 0, -t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t20, 0.2e1 * t34, 0, 0, 0, 0.2e1 * pkin(4) * t28, -0.2e1 * pkin(4) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t35, t11, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t28, 0, -t26 * pkin(7), -t28 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t14;
