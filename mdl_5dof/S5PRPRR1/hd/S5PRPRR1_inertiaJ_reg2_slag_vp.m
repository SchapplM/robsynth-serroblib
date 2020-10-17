% Calculate inertial parameters regressor of joint inertia matrix for
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
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:04
% EndTime: 2019-12-05 15:43:06
% DurationCPUTime: 0.38s
% Computational Cost: add. (307->44), mult. (621->83), div. (0->0), fcn. (750->6), ass. (0->34)
t25 = sin(pkin(9));
t26 = cos(pkin(9));
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t17 = t28 * t25 - t30 * t26;
t22 = -t26 * pkin(3) - pkin(2);
t12 = t17 * pkin(4) + t22;
t39 = 0.2e1 * t12;
t38 = 0.2e1 * t22;
t37 = 0.2e1 * t26;
t27 = sin(qJ(5));
t36 = t27 * pkin(4);
t29 = cos(qJ(5));
t35 = t29 * pkin(4);
t34 = pkin(6) + qJ(3);
t23 = t25 ^ 2;
t24 = t26 ^ 2;
t33 = t23 + t24;
t20 = t34 * t25;
t21 = t34 * t26;
t10 = -t30 * t20 - t28 * t21;
t11 = -t28 * t20 + t30 * t21;
t19 = t30 * t25 + t28 * t26;
t15 = t19 ^ 2;
t14 = t17 ^ 2;
t9 = -t27 * t17 + t29 * t19;
t7 = t29 * t17 + t27 * t19;
t6 = t9 ^ 2;
t5 = t7 ^ 2;
t4 = -t17 * pkin(7) + t11;
t3 = -t19 * pkin(7) + t10;
t2 = t27 * t3 + t29 * t4;
t1 = -t27 * t4 + t29 * t3;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 + t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * t10 + t19 * t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t1 + t9 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t23, t25 * t37, 0, t24, 0, 0, pkin(2) * t37, -0.2e1 * pkin(2) * t25, 0.2e1 * t33 * qJ(3), t33 * qJ(3) ^ 2 + pkin(2) ^ 2, t15, -0.2e1 * t19 * t17, 0, t14, 0, 0, t17 * t38, t19 * t38, -0.2e1 * t10 * t19 - 0.2e1 * t11 * t17, t10 ^ 2 + t11 ^ 2 + t22 ^ 2, t6, -0.2e1 * t9 * t7, 0, t5, 0, 0, t7 * t39, t9 * t39, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t7, t1 ^ 2 + t12 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t25, 0, -pkin(2), 0, 0, 0, 0, 0, 0, t17, t19, 0, t22, 0, 0, 0, 0, 0, 0, t7, t9, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t19, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t9, 0, (t27 * t9 - t29 * t7) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, 0, t10, -t11, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, (-t27 * t7 - t29 * t9) * pkin(4), (t1 * t29 + t2 * t27) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t36, 0, (t27 ^ 2 + t29 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
