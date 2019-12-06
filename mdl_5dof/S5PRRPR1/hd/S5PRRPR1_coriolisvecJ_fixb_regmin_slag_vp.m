% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:50
% EndTime: 2019-12-05 16:15:52
% DurationCPUTime: 0.35s
% Computational Cost: add. (392->75), mult. (757->120), div. (0->0), fcn. (480->6), ass. (0->57)
t47 = sin(pkin(9));
t48 = cos(pkin(9));
t69 = t47 ^ 2 + t48 ^ 2;
t46 = qJD(2) + qJD(3);
t52 = cos(qJ(3));
t67 = pkin(2) * qJD(3);
t62 = qJD(2) * t67;
t26 = qJD(4) * t46 + t52 * t62;
t78 = t69 * t26;
t49 = sin(qJ(5));
t51 = cos(qJ(5));
t28 = t47 * t51 + t48 * t49;
t14 = t28 * t46;
t77 = t46 * t69;
t22 = t28 * qJD(5);
t68 = pkin(2) * qJD(2);
t56 = -t52 * t68 + qJD(4);
t50 = sin(qJ(3));
t76 = t50 * pkin(2);
t75 = t52 * pkin(2);
t39 = -pkin(4) * t48 - pkin(3);
t11 = t39 * t46 + t56;
t71 = t51 * t48;
t72 = t49 * t47;
t27 = -t71 + t72;
t59 = t50 * t62;
t74 = t11 * t22 + t27 * t59;
t21 = t27 * qJD(5);
t73 = -t11 * t21 + t28 * t59;
t66 = t46 * t72;
t65 = t46 * t71;
t64 = t50 * t67;
t60 = t46 * t47 * t76;
t58 = (-qJD(3) + t46) * t68;
t57 = (-qJD(2) - t46) * t67;
t55 = t69 * (qJ(4) * t46 + t50 * t68);
t54 = t50 * t57;
t53 = t50 * t58;
t42 = t48 * pkin(7);
t38 = qJ(4) + t76;
t37 = t52 * t67 + qJD(4);
t35 = t47 * t59;
t34 = qJ(4) * t48 + t42;
t33 = (-pkin(7) - qJ(4)) * t47;
t32 = t39 - t75;
t31 = qJD(5) * t65;
t29 = -pkin(3) * t46 + t56;
t25 = t38 * t48 + t42;
t24 = (-pkin(7) - t38) * t47;
t18 = t22 * qJD(5);
t17 = t21 * qJD(5);
t12 = -t65 + t66;
t8 = t46 * t22;
t7 = -qJD(5) * t66 + t31;
t2 = -t14 * t21 + t28 * t7;
t1 = t12 * t21 - t14 * t22 - t27 * t7 - t28 * t8;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17; 0, 0, 0, 0, 0, t54, t52 * t57, t48 * t54, qJD(3) * t60 + t35, t37 * t77 + t78, t55 * t37 + t38 * t78 + (t29 + (-pkin(3) - t75) * qJD(2)) * t64, t2, t1, -t17, -t18, 0, t12 * t64 + t32 * t8 + ((-t24 * t49 - t25 * t51) * qJD(5) - t28 * t37) * qJD(5) + t74, t14 * t64 + t32 * t7 + ((-t24 * t51 + t25 * t49) * qJD(5) + t27 * t37) * qJD(5) + t73; 0, 0, 0, 0, 0, t53, t52 * t58, t48 * t53, -qJD(2) * t60 + t35, t56 * t77 + t78, t55 * qJD(4) + qJ(4) * t78 + (-t55 * t52 + (-pkin(3) * qJD(3) - t29) * t50) * t68, t2, t1, -t17, -t18, 0, t39 * t8 + ((-t33 * t49 - t34 * t51) * qJD(5) - t28 * qJD(4)) * qJD(5) + (-t50 * t12 + t22 * t52) * t68 + t74, t39 * t7 + ((-t33 * t51 + t34 * t49) * qJD(5) + t27 * qJD(4)) * qJD(5) + (-t50 * t14 - t21 * t52) * t68 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 * t46 ^ 2, -t46 * t55 + t59, 0, 0, 0, 0, 0, 0.2e1 * t14 * qJD(5), t31 + (-t12 - t66) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t12, -t12 ^ 2 + t14 ^ 2, t31 + (t12 - t66) * qJD(5), 0, 0, -t11 * t14 - t26 * t28, t11 * t12 + t26 * t27;];
tauc_reg = t3;
