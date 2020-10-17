% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:04
% EndTime: 2019-12-05 15:31:07
% DurationCPUTime: 0.41s
% Computational Cost: add. (145->52), mult. (336->82), div. (0->0), fcn. (313->4), ass. (0->42)
t23 = sin(pkin(8));
t41 = -0.2e1 * t23;
t24 = cos(pkin(8));
t40 = 0.2e1 * t24;
t25 = sin(qJ(4));
t39 = t25 * t23;
t26 = cos(qJ(4));
t17 = t26 * t23;
t38 = t26 * t24;
t19 = t23 ^ 2;
t20 = t24 ^ 2;
t37 = t19 + t20;
t21 = t25 ^ 2;
t22 = t26 ^ 2;
t13 = t21 + t22;
t36 = qJ(3) * t24;
t35 = qJ(3) * t25;
t34 = qJ(5) * t23;
t33 = t19 * qJ(3);
t32 = t23 * t40;
t31 = t26 * t36;
t9 = -t24 * pkin(3) - t23 * pkin(6) - pkin(2);
t6 = t26 * t9;
t30 = -t26 * t34 + t6;
t1 = (-pkin(4) - t35) * t24 + t30;
t2 = t31 + (t9 - t34) * t25;
t29 = t1 * t26 + t2 * t25;
t3 = -t24 * t35 + t6;
t4 = t25 * t9 + t31;
t28 = t4 * t25 + t3 * t26;
t27 = qJ(3) ^ 2;
t18 = t19 * t27;
t16 = t22 * t19;
t15 = t25 * t24;
t14 = t21 * t19;
t12 = -0.2e1 * t26 * t19 * t25;
t11 = t38 * t41;
t10 = t25 * t32;
t8 = (pkin(4) * t25 + qJ(3)) * t23;
t7 = t13 * t23;
t5 = t16 + t14 + t20;
t42 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t25 * t3 + t26 * t4 - t36) * t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t8 + (-t1 * t25 + t2 * t26) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t19, t32, 0, t20, 0, 0, pkin(2) * t40, pkin(2) * t41, 0.2e1 * t37 * qJ(3), pkin(2) ^ 2 + t20 * t27 + t18, t16, t12, t11, t14, t10, t20, -0.2e1 * t3 * t24 + 0.2e1 * t25 * t33, 0.2e1 * t4 * t24 + 0.2e1 * t26 * t33, t28 * t41, t3 ^ 2 + t4 ^ 2 + t18, t16, t12, t11, t14, t10, t20, -0.2e1 * t1 * t24 + 0.2e1 * t8 * t39, 0.2e1 * t8 * t17 + 0.2e1 * t2 * t24, t29 * t41, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t23, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t38, t15, -t7, t28, 0, 0, 0, 0, 0, 0, -t38, t15, -t7, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t17, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t17, 0, -pkin(4) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t39, -t24, t3, -t4, 0, 0, 0, 0, t17, 0, -t39, -t24, (-0.2e1 * pkin(4) - t35) * t24 + t30, -t2, -pkin(4) * t17, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t26 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t17, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t42;
