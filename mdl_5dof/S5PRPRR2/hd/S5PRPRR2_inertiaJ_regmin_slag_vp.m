% Calculate minimal parameter regressor of joint inertia matrix for
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
% MM_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t19 = cos(qJ(5));
t26 = 0.2e1 * t19;
t15 = cos(pkin(9));
t12 = t15 * pkin(2) + pkin(3);
t17 = sin(qJ(4));
t20 = cos(qJ(4));
t14 = sin(pkin(9));
t24 = pkin(2) * t14;
t6 = t20 * t12 - t17 * t24;
t4 = -pkin(4) - t6;
t25 = pkin(4) - t4;
t18 = sin(qJ(2));
t21 = cos(qJ(2));
t10 = t14 * t21 + t15 * t18;
t9 = -t14 * t18 + t15 * t21;
t2 = t17 * t10 - t20 * t9;
t23 = t2 * t19;
t7 = -t17 * t12 - t20 * t24;
t16 = sin(qJ(5));
t13 = t16 ^ 2;
t11 = t16 * t26;
t5 = pkin(7) - t7;
t3 = t20 * t10 + t17 * t9;
t1 = t2 * t16;
t8 = [1, 0, 0, 0, t10 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t21, -t18, (t10 * t14 + t15 * t9) * pkin(2), 0, -t2, -t3, 0, 0, 0, 0, 0, -t23, t1; 0, 1, 0, 0, (t14 ^ 2 + t15 ^ 2) * pkin(2) ^ 2, 1, 0.2e1 * t6, 0.2e1 * t7, t13, t11, 0, 0, 0, -0.2e1 * t4 * t19, 0.2e1 * t4 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t2, -t3, 0, 0, 0, 0, 0, -t23, t1; 0, 0, 0, 0, 0, 1, t6, t7, t13, t11, 0, 0, 0, t25 * t19, -t25 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0, t13, t11, 0, 0, 0, pkin(4) * t26, -0.2e1 * pkin(4) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t3, -t19 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t19, 0, -t16 * t5, -t19 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t19, 0, -t16 * pkin(7), -t19 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
