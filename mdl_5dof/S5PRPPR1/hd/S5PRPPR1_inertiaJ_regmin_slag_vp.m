% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t23 = sin(pkin(8));
t38 = -0.2e1 * t23;
t22 = sin(pkin(9));
t24 = cos(pkin(9));
t26 = sin(qJ(5));
t27 = cos(qJ(5));
t11 = -t26 * t22 + t27 * t24;
t8 = t11 * t23;
t37 = -0.2e1 * t8;
t25 = cos(pkin(8));
t36 = 0.2e1 * t25;
t14 = -t25 * pkin(3) - t23 * qJ(4) - pkin(2);
t31 = qJ(3) * t25;
t6 = t22 * t14 + t24 * t31;
t35 = t22 * t23;
t34 = t24 * t23;
t33 = t22 ^ 2 + t24 ^ 2;
t19 = t23 ^ 2;
t21 = t25 ^ 2;
t32 = t19 + t21;
t30 = t19 * qJ(3);
t10 = t24 * t14;
t5 = -t22 * t31 + t10;
t29 = t6 * t22 + t5 * t24;
t12 = t27 * t22 + t26 * t24;
t28 = qJ(3) ^ 2;
t17 = t23 * qJ(3);
t16 = t19 * t28;
t13 = pkin(4) * t35 + t17;
t7 = t12 * t23;
t4 = -pkin(6) * t35 + t6;
t3 = -pkin(6) * t34 + t10 + (-qJ(3) * t22 - pkin(4)) * t25;
t2 = t26 * t3 + t27 * t4;
t1 = -t26 * t4 + t27 * t3;
t9 = [1, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, t33 * t19 + t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t22 * t5 + t24 * t6 - t31) * t23, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, pkin(2) * t36, pkin(2) * t38, 0.2e1 * t32 * qJ(3), pkin(2) ^ 2 + t21 * t28 + t16, 0.2e1 * t22 * t30 - 0.2e1 * t5 * t25, 0.2e1 * t24 * t30 + 0.2e1 * t6 * t25, t29 * t38, t5 ^ 2 + t6 ^ 2 + t16, t8 ^ 2, t7 * t37, t25 * t37, t7 * t36, t21, -0.2e1 * t1 * t25 + 0.2e1 * t13 * t7, 0.2e1 * t13 * t8 + 0.2e1 * t2 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t25, t23, 0, -pkin(2), -t24 * t25, t22 * t25, -t33 * t23, t29, 0, 0, 0, 0, 0, -t11 * t25, t12 * t25; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t35, t34, 0, t17, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t25, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
