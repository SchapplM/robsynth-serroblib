% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:15
% EndTime: 2019-12-05 17:38:17
% DurationCPUTime: 0.36s
% Computational Cost: add. (346->95), mult. (763->141), div. (0->0), fcn. (450->4), ass. (0->74)
t47 = sin(qJ(5));
t48 = sin(qJ(4));
t49 = cos(qJ(5));
t50 = cos(qJ(4));
t22 = t47 * t50 + t49 * t48;
t19 = t22 * qJD(1);
t41 = qJD(4) + qJD(5);
t81 = t19 * t41;
t46 = pkin(1) + qJ(3);
t87 = qJD(1) * t46;
t86 = qJD(5) - t41;
t79 = t49 * t50;
t59 = t41 * t79;
t71 = qJD(5) * t47;
t73 = qJD(4) * t48;
t9 = -t47 * t73 - t48 * t71 + t59;
t85 = t9 * t41;
t45 = -pkin(6) + qJ(2);
t84 = pkin(7) - t45;
t74 = qJD(1) * t50;
t64 = t49 * t74;
t75 = qJD(1) * t48;
t65 = t47 * t75;
t18 = -t64 + t65;
t83 = t18 * t19;
t82 = t18 * t41;
t39 = qJD(1) * qJ(2) + qJD(3);
t30 = -pkin(6) * qJD(1) + t39;
t16 = -pkin(7) * t75 + t48 * t30;
t80 = t49 * t16;
t67 = qJD(1) * qJD(4);
t61 = t48 * t67;
t78 = -qJD(5) * t65 - t47 * t61;
t77 = t48 ^ 2 - t50 ^ 2;
t51 = qJD(4) ^ 2;
t52 = qJD(1) ^ 2;
t76 = -t51 - t52;
t72 = qJD(4) * t50;
t70 = t48 * qJD(2);
t31 = -qJD(2) + t87;
t69 = qJD(2) - t31;
t68 = qJD(1) * qJD(2);
t66 = pkin(4) * t74;
t40 = 0.2e1 * t68;
t63 = 0.2e1 * qJD(3) * qJD(1);
t27 = t84 * t50;
t17 = -pkin(7) * t74 + t50 * t30;
t13 = qJD(4) * pkin(4) + t17;
t62 = -pkin(4) * t41 - t13;
t60 = -0.2e1 * t50 * t67;
t36 = t48 * pkin(4) + t46;
t32 = pkin(4) * t72 + qJD(3);
t57 = qJD(2) + t31 + t87;
t56 = -t45 * t51 + t63;
t37 = t50 * t68;
t10 = t37 + (pkin(7) * qJD(1) - t30) * t73;
t11 = t30 * t72 + (-pkin(7) * t72 + t70) * qJD(1);
t21 = t36 * qJD(1) - qJD(2);
t55 = t49 * t10 - t47 * t11 + t21 * t18;
t54 = -t41 * t64 - t78;
t53 = t16 * t71 + (-t16 * t41 - t10) * t47 + t21 * t19;
t8 = t41 * t22;
t4 = t8 * qJD(1);
t26 = t84 * t48;
t25 = t32 * qJD(1);
t23 = -t47 * t48 + t79;
t15 = -qJD(4) * t27 + t70;
t14 = t50 * qJD(2) + t84 * t73;
t6 = t8 * t41;
t5 = qJD(1) * t59 + t78;
t3 = t18 ^ 2 - t19 ^ 2;
t2 = t54 - t82;
t1 = -t4 + t81;
t7 = [0, 0, 0, 0, t40, qJ(2) * t40, t40, t63, t39 * qJD(2) + t31 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t46) * qJD(1), t48 * t60, 0.2e1 * t77 * t67, -t51 * t48, -t51 * t50, 0, t56 * t48 + t57 * t72, t56 * t50 - t57 * t73, t18 * t8 - t4 * t23, t18 * t9 + t8 * t19 + t4 * t22 - t23 * t5, -t6, -t85, 0, t32 * t19 + t36 * t5 + t25 * t22 + t21 * t9 + (t49 * t14 - t47 * t15 + (t26 * t49 + t27 * t47) * qJD(5)) * t41, -t32 * t18 - t36 * t4 + t25 * t23 - t21 * t8 - (t47 * t14 + t49 * t15 + (t26 * t47 - t27 * t49) * qJD(5)) * t41; 0, 0, 0, 0, -t52, -t52 * qJ(2), -t52, 0, (-qJD(3) - t39) * qJD(1), 0, 0, 0, 0, 0, t60, 0.2e1 * t61, 0, 0, 0, 0, 0, t54 + t82, 0.2e1 * t81; 0, 0, 0, 0, 0, 0, 0, -t52, t69 * qJD(1), 0, 0, 0, 0, 0, t76 * t48, t76 * t50, 0, 0, 0, 0, 0, -qJD(1) * t19 - t6, qJD(1) * t18 - t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t52 * t48, -t77 * t52, 0, 0, 0, -t31 * t74 + t37, -t69 * t75, -t83, t3, t1, t2, 0, -t19 * t66 - (-t47 * t17 - t80) * t41 + (t62 * t47 - t80) * qJD(5) + t55, t18 * t66 + (t62 * qJD(5) + t17 * t41 - t11) * t49 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t3, t1, t2, 0, t55 + t86 * (-t47 * t13 - t80), (-t86 * t13 - t11) * t49 + t53;];
tauc_reg = t7;
