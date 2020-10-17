% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:26
% EndTime: 2019-12-05 16:04:27
% DurationCPUTime: 0.23s
% Computational Cost: add. (65->40), mult. (172->83), div. (0->0), fcn. (201->8), ass. (0->37)
t19 = cos(qJ(5));
t37 = 0.2e1 * t19;
t36 = 2 * qJ(3);
t14 = sin(pkin(5));
t18 = sin(qJ(2));
t35 = t14 * t18;
t21 = cos(qJ(2));
t34 = t14 * t21;
t16 = sin(qJ(5));
t17 = sin(qJ(4));
t33 = t16 * t17;
t32 = t16 * t19;
t20 = cos(qJ(4));
t31 = t16 * t20;
t22 = -pkin(2) - pkin(7);
t30 = t17 * t22;
t29 = t19 * t17;
t8 = t19 * t20;
t28 = t19 * t22;
t27 = t20 * t17;
t26 = t20 * t22;
t11 = t17 ^ 2;
t13 = t20 ^ 2;
t25 = -t11 - t13;
t24 = -0.2e1 * t27;
t23 = -pkin(4) * t20 - pkin(8) * t17;
t15 = cos(pkin(5));
t12 = t19 ^ 2;
t10 = t16 ^ 2;
t7 = t17 * pkin(4) - t20 * pkin(8) + qJ(3);
t6 = t15 * t20 - t17 * t34;
t5 = t15 * t17 + t20 * t34;
t4 = t16 * t7 + t17 * t28;
t3 = -t16 * t30 + t19 * t7;
t2 = t16 * t35 + t6 * t19;
t1 = -t6 * t16 + t19 * t35;
t9 = [1, 0, 0, 0, 0, 0, t15 ^ 2 + (t18 ^ 2 + t21 ^ 2) * t14 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t34, -t35, -t34, t35, (pkin(2) * t21 + qJ(3) * t18) * t14, 0, 0, 0, 0, 0, t17 * t35, t20 * t35, 0, 0, 0, 0, 0, t1 * t17 + t5 * t31, -t2 * t17 + t5 * t8; 0, 1, 0, 0, -0.2e1 * pkin(2), t36, pkin(2) ^ 2 + (qJ(3) ^ 2), t13, t24, 0, 0, 0, t17 * t36, t20 * t36, t12 * t13, -0.2e1 * t13 * t32, t27 * t37, t16 * t24, t11, -0.2e1 * t13 * t22 * t16 + 0.2e1 * t3 * t17, -0.2e1 * t13 * t28 - 0.2e1 * t4 * t17; 0, 0, 0, 0, 0, 0, -t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 * t16, t25 * t19; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t5 * t19, t5 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t17, 0, t26, -t30, t16 * t8, (-t10 + t12) * t20, t33, t29, 0, t23 * t16 + t19 * t26, -t16 * t26 + t23 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t17, 0, 0, 0, 0, 0, t8, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t10, 0.2e1 * t32, 0, 0, 0, pkin(4) * t37, -0.2e1 * pkin(4) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t31, t17, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t19, 0, -t16 * pkin(8), -t19 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
