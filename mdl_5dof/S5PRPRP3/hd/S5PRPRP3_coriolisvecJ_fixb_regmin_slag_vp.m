% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [5x14]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:41
% EndTime: 2019-12-05 15:33:43
% DurationCPUTime: 0.32s
% Computational Cost: add. (394->82), mult. (953->130), div. (0->0), fcn. (635->6), ass. (0->61)
t36 = sin(qJ(4));
t38 = cos(qJ(4));
t39 = cos(qJ(2));
t54 = t39 * qJD(1);
t25 = qJD(2) * pkin(2) + t54;
t35 = cos(pkin(8));
t37 = sin(qJ(2));
t58 = qJD(1) * t37;
t27 = t35 * t58;
t34 = sin(pkin(8));
t13 = t34 * t25 + t27;
t49 = t13 + (pkin(6) + qJ(5)) * qJD(2);
t4 = t38 * qJD(3) - t49 * t36;
t5 = t36 * qJD(3) + t49 * t38;
t60 = qJD(4) * pkin(4);
t3 = t4 + t60;
t66 = t3 - t4;
t65 = t35 * pkin(2);
t40 = qJD(4) ^ 2;
t64 = t40 * t36;
t63 = t40 * t38;
t32 = t36 ^ 2;
t33 = t38 ^ 2;
t62 = t32 - t33;
t61 = t32 + t33;
t29 = t34 * pkin(2) + pkin(6);
t59 = qJ(5) + t29;
t56 = t36 * qJD(4);
t55 = t38 * qJD(4);
t53 = qJD(2) * qJD(4);
t52 = -t38 * pkin(4) - pkin(3);
t51 = t36 * t53;
t26 = t34 * t58;
t12 = t35 * t25 - t26;
t22 = t34 * t37 - t35 * t39;
t18 = t22 * qJD(2);
t15 = qJD(1) * t18;
t8 = -qJD(2) * pkin(3) - t12;
t50 = -t8 * qJD(2) + t15;
t48 = qJD(4) * t59;
t47 = qJD(2) * qJD(5) - t15;
t46 = t3 * t36 - t5 * t38;
t23 = t34 * t39 + t35 * t37;
t16 = t23 * qJD(2);
t14 = qJD(1) * t16;
t17 = t34 * t54 + t27;
t44 = qJD(2) * t17 - t29 * t40 - t14;
t19 = t35 * t54 - t26;
t43 = qJD(4) * (qJD(2) * (-pkin(3) - t65) + t19 + t8);
t1 = t4 * qJD(4) + t47 * t38;
t2 = -qJD(4) * t5 - t47 * t36;
t42 = t1 * t38 - t2 * t36 + (-t3 * t38 - t36 * t5) * qJD(4);
t41 = qJD(2) ^ 2;
t28 = pkin(4) * t51;
t21 = t59 * t38;
t20 = t59 * t36;
t11 = -t36 * qJD(5) - t38 * t48;
t10 = t38 * qJD(5) - t36 * t48;
t7 = t28 + t14;
t6 = t52 * qJD(2) + qJD(5) - t12;
t9 = [0, 0, -t41 * t37, -t41 * t39, -t12 * t16 - t13 * t18 + t14 * t22 - t15 * t23, 0, 0, 0, 0, 0, t18 * t56 - t23 * t63 + (-t16 * t38 + t22 * t56) * qJD(2), t18 * t55 + t23 * t64 + (t16 * t36 + t22 * t55) * qJD(2), -t61 * t18 * qJD(2), t6 * t16 + t46 * t18 + t7 * t22 + t42 * t23; 0, 0, 0, 0, t12 * t17 - t13 * t19 + (-t14 * t35 - t15 * t34) * pkin(2), 0.2e1 * t38 * t51, -0.2e1 * t62 * t53, t63, -t64, 0, t36 * t43 + t44 * t38, -t44 * t36 + t38 * t43, (t10 * t38 - t11 * t36 - t61 * t19 + (t20 * t38 - t21 * t36) * qJD(4)) * qJD(2) + t42, t1 * t21 + t5 * t10 - t2 * t20 + t3 * t11 + t7 * (t52 - t65) + (pkin(4) * t56 - t17) * t6 + t46 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63, 0, -t46 * qJD(4) + t1 * t36 + t2 * t38; 0, 0, 0, 0, 0, -t36 * t41 * t38, t62 * t41, 0, 0, 0, t50 * t36, t50 * t38, (-t60 + t66) * t38 * qJD(2), t66 * t5 + (-qJD(2) * t36 * t6 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 * t41, t28 + (t23 * qJD(1) + t46) * qJD(2);];
tauc_reg = t9;
