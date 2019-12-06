% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPPRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:07
% EndTime: 2019-12-05 14:58:09
% DurationCPUTime: 0.34s
% Computational Cost: add. (512->68), mult. (1333->127), div. (0->0), fcn. (1150->8), ass. (0->59)
t29 = sin(pkin(9));
t31 = cos(pkin(9));
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t62 = -t34 * t29 + t36 * t31;
t30 = sin(pkin(8));
t51 = qJD(1) * t30;
t21 = t31 * qJD(2) - t29 * t51;
t22 = t29 * qJD(2) + t31 * t51;
t11 = t36 * t21 - t34 * t22;
t9 = t11 * qJD(4);
t12 = t34 * t21 + t36 * t22;
t10 = t12 * qJD(4);
t24 = t36 * t29 + t34 * t31;
t17 = t24 * t30;
t61 = t10 * t17;
t60 = t10 * t62;
t33 = sin(qJ(5));
t35 = cos(qJ(5));
t59 = t33 * t35;
t37 = qJD(5) ^ 2;
t56 = t37 * t33;
t55 = t37 * t35;
t27 = t33 ^ 2;
t28 = t35 ^ 2;
t54 = t27 - t28;
t53 = t27 + t28;
t52 = qJD(4) * pkin(4);
t50 = qJD(4) * t30;
t49 = qJD(5) * t33;
t48 = qJD(5) * t35;
t19 = t62 * qJD(4);
t47 = t19 * qJD(4);
t46 = qJD(4) * qJD(5);
t38 = qJD(4) ^ 2;
t45 = t38 * t59;
t7 = -t11 - t52;
t44 = -qJD(4) * t7 - t9;
t43 = t46 * t59;
t32 = cos(pkin(8));
t26 = -t32 * qJD(1) + qJD(3);
t8 = qJD(4) * pkin(6) + t12;
t5 = t35 * t26 - t33 * t8;
t6 = t33 * t26 + t35 * t8;
t42 = t33 * t5 - t35 * t6;
t18 = t62 * t30;
t14 = t35 * t18 - t32 * t33;
t13 = -t33 * t18 - t32 * t35;
t41 = pkin(6) * t37;
t40 = qJD(5) * (t11 + t7 - t52);
t1 = t5 * qJD(5) + t35 * t9;
t2 = -t6 * qJD(5) - t33 * t9;
t39 = t1 * t35 - t2 * t33 + (-t33 * t6 - t35 * t5) * qJD(5);
t20 = t24 * qJD(4);
t16 = t62 * t50;
t15 = t24 * t50;
t4 = -qJD(5) * t14 + t33 * t15;
t3 = qJD(5) * t13 - t35 * t15;
t23 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * qJD(4), t15 * qJD(4), 0, -t11 * t16 - t12 * t15 + t9 * t18 + t61, 0, 0, 0, 0, 0, 0, t4 * qJD(5) + (-t16 * t35 + t17 * t49) * qJD(4), -t3 * qJD(5) + (t16 * t33 + t17 * t48) * qJD(4), (t3 * t35 - t33 * t4 + (-t13 * t35 - t14 * t33) * qJD(5)) * qJD(4), t1 * t14 + t2 * t13 + t7 * t16 + t6 * t3 + t5 * t4 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20 * qJD(4), -t47, 0, -t11 * t20 + t12 * t19 + t9 * t24 - t60, 0, 0, 0, 0, 0, 0, -t19 * t49 - t24 * t55 + (-t20 * t35 - t49 * t62) * qJD(4), -t19 * t48 + t24 * t56 + (t20 * t33 - t48 * t62) * qJD(4), t53 * t47, -t42 * t19 + t7 * t20 + t39 * t24 - t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t55, 0, -t42 * qJD(5) + t1 * t33 + t2 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t43, -0.2e1 * t54 * t46, t55, -0.2e1 * t43, -t56, 0, t33 * t40 - t41 * t35, t41 * t33 + t35 * t40, -t53 * t9 + t39, -t10 * pkin(4) + t39 * pkin(6) + t42 * t11 - t7 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t54 * t38, 0, t45, 0, 0, t44 * t33, t44 * t35, 0, 0;];
tauc_reg = t23;
