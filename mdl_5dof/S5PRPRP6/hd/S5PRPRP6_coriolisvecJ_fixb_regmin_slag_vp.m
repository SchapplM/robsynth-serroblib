% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:19
% EndTime: 2019-12-05 15:41:21
% DurationCPUTime: 0.32s
% Computational Cost: add. (318->84), mult. (662->113), div. (0->0), fcn. (318->4), ass. (0->64)
t24 = sin(qJ(4));
t25 = sin(qJ(2));
t53 = t25 * qJD(1);
t45 = qJD(2) * t53;
t19 = t24 * t45;
t28 = -pkin(2) - pkin(6);
t27 = cos(qJ(2));
t52 = t27 * qJD(1);
t40 = qJD(3) - t52;
t11 = t28 * qJD(2) + t40;
t26 = cos(qJ(4));
t65 = t26 * t11;
t2 = t19 + (qJD(5) + t65) * qJD(4);
t20 = t26 * t45;
t55 = qJD(4) * t24;
t3 = t11 * t55 - t20;
t42 = qJD(4) * pkin(4) - qJD(5);
t6 = -t42 - t65;
t50 = qJD(4) * qJ(5);
t8 = t24 * t11 + t50;
t37 = t24 * t6 + t26 * t8;
t31 = qJD(4) * t37 + t2 * t24 - t3 * t26;
t13 = t24 * pkin(4) - t26 * qJ(5) + qJ(3);
t56 = qJD(2) * t13;
t7 = t53 + t56;
t57 = t7 * qJD(2);
t68 = t31 - t57;
t36 = pkin(4) * t26 + qJ(5) * t24;
t9 = qJD(4) * t36 - t26 * qJD(5) + qJD(3);
t1 = (t9 + t52) * qJD(2);
t29 = qJD(4) ^ 2;
t64 = t28 * t29;
t67 = (-t9 + t52) * qJD(2) - t1 + t64;
t51 = qJD(2) * qJ(3);
t18 = t51 + t53;
t66 = t18 * t27;
t30 = qJD(2) ^ 2;
t63 = t30 * t25;
t62 = t30 * t27;
t22 = t24 ^ 2;
t23 = t26 ^ 2;
t61 = t22 - t23;
t60 = t22 + t23;
t59 = t29 + t30;
t58 = qJD(2) * pkin(2);
t54 = t18 * qJD(2);
t49 = qJD(2) * qJD(4);
t48 = 0.2e1 * t49;
t46 = t59 * t27;
t43 = -t18 + t53;
t41 = -0.2e1 * t24 * t49;
t38 = t24 * t8 - t26 * t6;
t35 = -t43 + t51;
t34 = qJD(4) * (-t53 + t7 + t56);
t14 = (qJD(3) + t52) * qJD(2);
t33 = t40 * qJD(2) + t14 - t64;
t21 = t26 * t30 * t24;
t17 = t40 - t58;
t16 = t59 * t26;
t15 = t59 * t24;
t12 = t36 * qJD(2);
t5 = t25 * t41 + t26 * t46;
t4 = t26 * t25 * t48 + t24 * t46;
t10 = [0, 0, -t63, -t62, t63, t62, t14 * t25 + (t66 + (t17 - t52) * t25) * qJD(2), 0, 0, 0, 0, 0, t4, t5, t4, -t60 * t63, -t5, (t38 * qJD(2) + t1) * t25 - t68 * t27; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t14 * qJ(3) + t18 * qJD(3) + (-t66 + (-t17 - t58) * t25) * qJD(1), t26 * t41, t61 * t48, -t29 * t24, -t29 * t26, 0, qJD(4) * t26 * t35 + t24 * t33, t26 * t33 - t35 * t55, -t67 * t24 + t26 * t34, t60 * t45 - t31, t24 * t34 + t67 * t26, t1 * t13 + t7 * t9 + (-t38 * t25 - t27 * t7) * qJD(1) + t31 * t28; 0, 0, 0, 0, 0, -t30, t43 * qJD(2), 0, 0, 0, 0, 0, -t15, -t16, -t15, 0, t16, t68; 0, 0, 0, 0, 0, 0, 0, t21, -t61 * t30, 0, 0, 0, -t26 * t54 + t20, t24 * t54 - t19, t20 + (-t12 * t24 - t26 * t7) * qJD(2), ((t8 - t50) * t26 + (t42 + t6) * t24) * qJD(2), 0.2e1 * qJD(4) * qJD(5) + t19 + (t12 * t26 - t24 * t7) * qJD(2), -t3 * pkin(4) + t2 * qJ(5) + t8 * qJD(5) - t11 * t37 - t7 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t23 * t30 - t29, -t8 * qJD(4) + t26 * t57 + t3;];
tauc_reg = t10;
