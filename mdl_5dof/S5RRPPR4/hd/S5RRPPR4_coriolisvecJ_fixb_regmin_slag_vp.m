% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:48
% EndTime: 2019-12-31 19:27:49
% DurationCPUTime: 0.31s
% Computational Cost: add. (439->86), mult. (738->122), div. (0->0), fcn. (357->6), ass. (0->65)
t46 = sin(pkin(8));
t47 = cos(pkin(8));
t51 = cos(qJ(2));
t75 = pkin(1) * qJD(1);
t69 = t51 * t75;
t49 = sin(qJ(2));
t70 = t49 * t75;
t15 = t46 * t69 - t47 * t70;
t42 = qJD(1) + qJD(2);
t57 = (t46 * qJD(3) - t15) * t42;
t74 = pkin(1) * qJD(2);
t66 = qJD(1) * t74;
t34 = t51 * t66;
t38 = t42 * qJD(3);
t21 = t34 + t38;
t63 = t49 * t66;
t6 = t46 * t21 - t47 * t63;
t84 = t57 + t6;
t41 = t42 ^ 2;
t83 = 2 * qJD(5);
t52 = -pkin(2) - pkin(3);
t71 = t51 * t74;
t33 = qJD(3) + t71;
t72 = t49 * t74;
t11 = t46 * t33 - t47 * t72;
t82 = t11 * t42;
t53 = qJD(5) ^ 2;
t77 = t47 * qJ(3) + t46 * t52;
t81 = (-pkin(7) + t77) * t53;
t48 = sin(qJ(5));
t80 = t42 * t48;
t40 = t53 * t48;
t50 = cos(qJ(5));
t79 = t53 * t50;
t7 = t47 * t21 + t46 * t63;
t68 = -t51 * pkin(1) - pkin(2);
t35 = -pkin(3) + t68;
t36 = t49 * pkin(1) + qJ(3);
t78 = t46 * t35 + t47 * t36;
t76 = t48 ^ 2 - t50 ^ 2;
t73 = t50 * t83;
t61 = qJD(3) - t69;
t14 = t52 * t42 + t61;
t25 = t42 * qJ(3) + t70;
t3 = t47 * t14 - t46 * t25;
t1 = t42 * pkin(4) - t3;
t67 = t1 * t42 - t7;
t16 = (t46 * t49 + t47 * t51) * t75;
t64 = t47 * qJD(3) - t16;
t4 = t46 * t14 + t47 * t25;
t62 = t3 * t46 - t4 * t47;
t60 = -t53 * (-pkin(7) + t78) + t82;
t59 = t47 * t35 - t46 * t36;
t58 = -t46 * qJ(3) + t47 * t52;
t56 = t42 * t69 - t34;
t12 = t47 * t33 + t46 * t72;
t55 = qJD(5) * (-t42 * (pkin(4) - t59) - t1 - t12);
t54 = qJD(5) * (-(pkin(4) - t58) * t42 - t1 - t64);
t26 = t73 * t80;
t22 = -t42 * pkin(2) + t61;
t19 = (-qJD(1) - t42) * t72;
t18 = (-qJD(2) + t42) * t70;
t13 = -0.2e1 * t76 * t42 * qJD(5);
t5 = t6 * t50;
t2 = [0, 0, 0, 0, t19, -t42 * t71 - t34, t19, t33 * t42 + t21, t21 * t36 + t25 * t33 + (t68 * qJD(1) + t22) * t72, t6 + t82, t12 * t42 + t7, -t3 * t11 + t4 * t12 - t6 * t59 + t7 * t78, t26, t13, -t79, t40, 0, t48 * t55 + t60 * t50 + t5, (-t6 - t60) * t48 + t50 * t55; 0, 0, 0, 0, t18, t56, t18, 0.2e1 * t38 - t56, t21 * qJ(3) + t25 * qJD(3) + (-t25 * t51 + (-pkin(2) * qJD(2) - t22) * t49) * t75, t84, t64 * t42 + t7, -t62 * qJD(3) + t3 * t15 - t4 * t16 - t6 * t58 + t7 * t77, t26, t13, -t79, t40, 0, t5 + (t57 - t81) * t50 + t48 * t54, (t81 - t84) * t48 + t50 * t54; 0, 0, 0, 0, 0, 0, 0, -t41, -t25 * t42 + t63, -t46 * t41, -t47 * t41, t62 * t42 + t7 * t46 - t6 * t47, 0, 0, 0, 0, 0, -t46 * t79 + (-t46 * t50 * t42 + t47 * t48 * t83) * t42, t46 * t40 + (t46 * t80 + t47 * t73) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t41 * t50, t76 * t41, 0, 0, 0, t67 * t48, t67 * t50;];
tauc_reg = t2;
