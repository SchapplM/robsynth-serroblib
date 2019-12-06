% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPPRR2
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPPRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:38
% EndTime: 2019-12-05 14:59:40
% DurationCPUTime: 0.40s
% Computational Cost: add. (533->82), mult. (1298->153), div. (0->0), fcn. (1088->8), ass. (0->69)
t29 = sin(pkin(9));
t31 = cos(pkin(9));
t30 = sin(pkin(8));
t56 = qJD(1) * t30;
t24 = t29 * qJD(2) + t31 * t56;
t32 = cos(pkin(8));
t26 = -t32 * qJD(1) + qJD(3);
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t42 = t34 * t24 - t36 * t26;
t10 = t42 * qJD(4);
t53 = qJD(4) * qJD(5);
t71 = -0.2e1 * t53;
t37 = qJD(5) ^ 2;
t38 = qJD(4) ^ 2;
t70 = (t37 + t38) * t34;
t16 = t36 * t24 + t34 * t26;
t9 = t16 * qJD(4);
t69 = t34 * t9;
t65 = t30 * t31;
t19 = t32 * t36 + t34 * t65;
t68 = t9 * t19;
t67 = t29 * t30;
t66 = t29 * t36;
t33 = sin(qJ(5));
t35 = cos(qJ(5));
t64 = t33 * t35;
t63 = t38 * t34;
t62 = t38 * t36;
t27 = t33 ^ 2;
t28 = t35 ^ 2;
t61 = t27 - t28;
t60 = t27 + t28;
t58 = qJD(4) * pkin(4);
t7 = t42 - t58;
t57 = qJD(4) * t7;
t55 = qJD(4) * t34;
t54 = qJD(5) * t19;
t52 = t36 * t65;
t51 = t38 * t64;
t50 = t29 * t62;
t49 = t29 * t55;
t48 = t10 - t57;
t47 = t33 * t49;
t46 = t35 * t49;
t45 = t36 * t71;
t44 = t53 * t64;
t23 = -t31 * qJD(2) + t29 * t56;
t8 = qJD(4) * pkin(6) + t16;
t5 = t35 * t23 - t33 * t8;
t6 = t33 * t23 + t35 * t8;
t43 = t33 * t5 - t35 * t6;
t20 = -t32 * t34 + t52;
t11 = -t20 * t33 + t35 * t67;
t12 = t20 * t35 + t33 * t67;
t22 = -t33 * t31 + t35 * t66;
t21 = -t35 * t31 - t33 * t66;
t41 = pkin(6) * t37;
t40 = qJD(5) * (-t42 + t7 - t58);
t1 = t5 * qJD(5) - t35 * t10;
t2 = -t6 * qJD(5) + t33 * t10;
t39 = t1 * t35 - t2 * t33 + (-t33 * t6 - t35 * t5) * qJD(5);
t18 = t19 * qJD(4);
t17 = qJD(4) * t52 - t32 * t55;
t14 = -qJD(5) * t22 + t47;
t13 = qJD(5) * t21 - t46;
t4 = qJD(5) * t11 - t18 * t35;
t3 = -qJD(5) * t12 + t18 * t33;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * qJD(4), t18 * qJD(4), 0, -t10 * t20 - t16 * t18 + t17 * t42 + t68, 0, 0, 0, 0, 0, 0, t3 * qJD(5) + (-t17 * t35 + t33 * t54) * qJD(4), -t4 * qJD(5) + (t17 * t33 + t35 * t54) * qJD(4), (-t3 * t33 + t35 * t4 + (-t11 * t35 - t12 * t33) * qJD(5)) * qJD(4), t1 * t12 + t2 * t11 + t7 * t17 + t5 * t3 + t6 * t4 + t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t29 * t63, 0, (-t10 * t36 + t69 + (-t16 * t34 + t36 * t42) * qJD(4)) * t29, 0, 0, 0, 0, 0, 0, -t35 * t50 + (t14 + t47) * qJD(5), t33 * t50 + (-t13 + t46) * qJD(5), (t13 * t35 - t14 * t33 + (-t21 * t35 - t22 * t33) * qJD(5)) * qJD(4), t1 * t22 + t6 * t13 + t5 * t14 + t2 * t21 + (t36 * t57 + t69) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t62, 0, -t10 * t34 - t9 * t36 + (t16 * t36 + t34 * t42) * qJD(4), 0, 0, 0, 0, 0, 0, t33 * t45 - t35 * t70, t33 * t70 + t35 * t45, t60 * t62, (-t43 * qJD(4) - t9) * t36 + (t39 + t57) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t44, t61 * t71, t37 * t35, -0.2e1 * t44, -t37 * t33, 0, t33 * t40 - t41 * t35, t41 * t33 + t35 * t40, t60 * t10 + t39, -t9 * pkin(4) + t39 * pkin(6) - t7 * t16 - t42 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t61 * t38, 0, t51, 0, 0, t48 * t33, t48 * t35, 0, 0;];
tauc_reg = t15;
