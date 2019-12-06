% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:41
% EndTime: 2019-12-05 15:24:43
% DurationCPUTime: 0.31s
% Computational Cost: add. (321->78), mult. (876->122), div. (0->0), fcn. (671->8), ass. (0->60)
t45 = sin(pkin(9));
t47 = cos(pkin(9));
t67 = t45 ^ 2 + t47 ^ 2;
t72 = qJD(2) * t67;
t49 = sin(qJ(5));
t51 = cos(qJ(5));
t29 = t49 * t45 - t51 * t47;
t24 = t29 * qJD(5);
t31 = t51 * t45 + t49 * t47;
t23 = t31 * qJD(2);
t46 = sin(pkin(8));
t50 = sin(qJ(2));
t66 = qJD(1) * t50;
t37 = t46 * t66;
t48 = cos(pkin(8));
t52 = cos(qJ(2));
t62 = t52 * qJD(1);
t61 = t48 * t62;
t22 = -t37 + t61;
t71 = t22 - qJD(4);
t70 = t48 * pkin(2);
t40 = t46 * pkin(2) + qJ(4);
t69 = pkin(6) + t40;
t30 = t46 * t52 + t48 * t50;
t17 = t30 * qJD(2);
t12 = qJD(1) * t17;
t28 = t46 * t50 - t48 * t52;
t68 = t12 * t28;
t36 = qJD(2) * pkin(2) + t62;
t38 = t48 * t66;
t11 = t46 * t36 + t38;
t65 = qJD(2) * t45;
t64 = qJD(2) * t47;
t63 = t24 * qJD(5);
t60 = t49 * t65;
t59 = t51 * t64;
t58 = -t47 * pkin(4) - pkin(3);
t34 = qJD(2) * t61;
t9 = t34 + (qJD(4) - t37) * qJD(2);
t57 = t67 * t9;
t10 = t48 * t36 - t37;
t18 = t46 * t62 + t38;
t56 = qJD(2) * t18 - t12;
t55 = t67 * (qJD(2) * qJ(4) + t11);
t54 = qJD(4) - t10;
t25 = t31 * qJD(5);
t53 = qJD(2) ^ 2;
t35 = qJD(5) * t59;
t33 = t58 - t70;
t27 = t69 * t47;
t26 = t69 * t45;
t21 = t28 * qJD(2);
t19 = -t59 + t60;
t16 = t25 * qJD(5);
t15 = qJD(2) * t25;
t14 = -qJD(5) * t60 + t35;
t13 = -qJD(2) * t37 + t34;
t7 = -qJD(2) * pkin(3) + t54;
t5 = t58 * qJD(2) + t54;
t1 = [0, 0, -t53 * t50, -t53 * t52, -t10 * t17 - t11 * t21 + t13 * t30 + t68, -t17 * t64, t17 * t65, -t21 * t72, t7 * t17 - t55 * t21 + t30 * t57 + t68, 0, 0, 0, 0, 0, t28 * t15 + t17 * t19 + (t31 * t21 + t30 * t24) * qJD(5), t28 * t14 + t17 * t23 + (-t29 * t21 + t30 * t25) * qJD(5); 0, 0, 0, 0, t10 * t18 - t11 * t22 + (-t12 * t48 + t13 * t46) * pkin(2), t56 * t47, -t56 * t45, -t71 * t72 + t57, t12 * (-pkin(3) - t70) - t7 * t18 + t40 * t57 - t71 * t55, t14 * t31 - t23 * t24, -t14 * t29 - t31 * t15 + t24 * t19 - t23 * t25, -t63, -t16, 0, t12 * t29 + t33 * t15 - t18 * t19 + t5 * t25 + ((t26 * t49 - t27 * t51) * qJD(5) + t71 * t31) * qJD(5), t12 * t31 + t33 * t14 - t18 * t23 - t5 * t24 + ((t26 * t51 + t27 * t49) * qJD(5) - t71 * t29) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t63; 0, 0, 0, 0, 0, 0, 0, -t67 * t53, (t30 * qJD(1) - t55) * qJD(2), 0, 0, 0, 0, 0, 0.2e1 * t23 * qJD(5), t35 + (-t19 - t60) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t19, -t19 ^ 2 + t23 ^ 2, t35 + (t19 - t60) * qJD(5), 0, 0, -t5 * t23 - t31 * t9, t5 * t19 + t29 * t9;];
tauc_reg = t1;
