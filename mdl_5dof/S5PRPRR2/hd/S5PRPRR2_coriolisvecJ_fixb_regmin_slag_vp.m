% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:13
% EndTime: 2019-12-05 15:45:14
% DurationCPUTime: 0.30s
% Computational Cost: add. (426->72), mult. (983->109), div. (0->0), fcn. (746->8), ass. (0->61)
t51 = qJD(5) ^ 2;
t44 = cos(pkin(9));
t37 = t44 * pkin(2) + pkin(3);
t46 = sin(qJ(4));
t49 = cos(qJ(4));
t43 = sin(pkin(9));
t81 = pkin(2) * t43;
t57 = t46 * t37 + t49 * t81;
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t33 = t43 * t50 + t44 * t47;
t29 = t33 * qJD(1);
t32 = -t43 * t47 + t44 * t50;
t31 = t32 * qJD(1);
t40 = qJD(2) + qJD(4);
t66 = (-t57 * qJD(4) + t49 * t29 + t46 * t31) * t40;
t83 = (pkin(7) + t57) * t51 - t66;
t36 = qJD(2) * pkin(2) + t50 * qJD(1);
t69 = qJD(1) * t47;
t19 = t44 * t36 - t43 * t69;
t17 = qJD(2) * pkin(3) + t19;
t28 = t33 * qJD(2);
t24 = qJD(1) * t28;
t30 = t32 * qJD(2);
t25 = qJD(1) * t30;
t20 = t43 * t36 + t44 * t69;
t75 = t46 * t20;
t54 = -(qJD(4) * t17 + t25) * t49 + qJD(4) * t75 + t46 * t24;
t11 = t46 * t17 + t49 * t20;
t3 = t11 * qJD(4) + t49 * t24 + t46 * t25;
t45 = sin(qJ(5));
t48 = cos(qJ(5));
t10 = t49 * t17 - t75;
t80 = t40 * pkin(4);
t8 = -t10 - t80;
t70 = qJD(5) * t8;
t82 = t3 * t45 + t48 * t70;
t15 = t46 * t32 + t49 * t33;
t79 = (t15 * qJD(4) + t49 * t28 + t46 * t30) * t40;
t78 = t11 * t40;
t76 = t45 * t48;
t74 = t51 * t45;
t56 = t49 * t37 - t46 * t81;
t72 = -t56 * qJD(4) - t46 * t29 + t49 * t31;
t71 = t45 ^ 2 - t48 ^ 2;
t68 = 0.2e1 * qJD(5) * t40;
t67 = -t8 * t40 + t54;
t62 = pkin(7) * t51 - t78;
t61 = t15 * t51 + t79;
t60 = t49 * t32 - t46 * t33;
t59 = qJD(5) * (t10 - t80);
t4 = t60 * qJD(4) - t46 * t28 + t49 * t30;
t58 = qJD(5) * (-t40 * t60 - t4);
t55 = qJD(5) * ((-pkin(4) - t56) * t40 + t72);
t52 = qJD(2) ^ 2;
t39 = t40 ^ 2;
t38 = t51 * t48;
t35 = t68 * t76;
t23 = t71 * t68;
t6 = t45 * t70;
t1 = [0, 0, -t52 * t47, -t52 * t50, -t19 * t28 + t20 * t30 - t24 * t32 + t25 * t33, 0, -t79, -t4 * t40, 0, 0, 0, 0, 0, t45 * t58 - t61 * t48, t61 * t45 + t48 * t58; 0, 0, 0, 0, t19 * t29 - t20 * t31 + (-t24 * t44 + t25 * t43) * pkin(2), 0, -t3 + t66, t72 * t40 + t54, t35, -t23, t38, -t74, 0, t6 + t45 * t55 + (-t3 - t83) * t48, t83 * t45 + t48 * t55 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t38; 0, 0, 0, 0, 0, 0, -t3 + t78, t10 * t40 + t54, t35, -t23, t38, -t74, 0, t6 + t45 * t59 + (-t3 - t62) * t48, t62 * t45 + t48 * t59 + t82; 0, 0, 0, 0, 0, 0, 0, 0, -t39 * t76, t71 * t39, 0, 0, 0, t67 * t45, t67 * t48;];
tauc_reg = t1;
