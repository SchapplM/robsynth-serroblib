% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRR1
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
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t18 = sin(pkin(9));
t19 = cos(pkin(9));
t21 = sin(qJ(4));
t23 = cos(qJ(4));
t10 = t21 * t18 - t23 * t19;
t15 = -t19 * pkin(3) - pkin(2);
t32 = 0.2e1 * t10 * pkin(4) + 0.2e1 * t15;
t31 = 0.2e1 * t15;
t20 = sin(qJ(5));
t30 = t20 * pkin(4);
t22 = cos(qJ(5));
t29 = t22 * pkin(4);
t28 = pkin(6) + qJ(3);
t27 = t18 ^ 2 + t19 ^ 2;
t12 = t28 * t18;
t13 = t28 * t19;
t26 = -t23 * t12 - t21 * t13;
t25 = t21 * t12 - t23 * t13;
t11 = t23 * t18 + t21 * t19;
t6 = -t20 * t10 + t22 * t11;
t5 = t22 * t10 + t20 * t11;
t4 = -t10 * pkin(7) - t25;
t3 = -t11 * pkin(7) + t26;
t2 = -t20 * t3 - t22 * t4;
t1 = -t20 * t4 + t22 * t3;
t7 = [1, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0.2e1 * pkin(2) * t19, -0.2e1 * pkin(2) * t18, 0.2e1 * t27 * qJ(3), t27 * qJ(3) ^ 2 + pkin(2) ^ 2, t11 ^ 2, -0.2e1 * t11 * t10, 0, 0, 0, t10 * t31, t11 * t31, t6 ^ 2, -0.2e1 * t6 * t5, 0, 0, 0, t5 * t32, t6 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t19, t18, 0, -pkin(2), 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, 0, t26, t25, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t29, -0.2e1 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t7;
