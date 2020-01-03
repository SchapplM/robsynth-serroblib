% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tau_reg [4x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:45
% EndTime: 2019-12-31 16:40:46
% DurationCPUTime: 0.43s
% Computational Cost: add. (294->114), mult. (695->153), div. (0->0), fcn. (487->6), ass. (0->74)
t78 = (qJ(2) * qJDD(1));
t77 = (qJD(1) * qJD(2));
t73 = 2 * t77;
t94 = (t73 + t78) * qJ(2);
t93 = t77 + t78;
t54 = cos(qJ(1));
t46 = g(2) * t54;
t52 = sin(qJ(1));
t92 = g(1) * t52;
t71 = -t46 + t92;
t49 = sin(pkin(6));
t53 = cos(qJ(4));
t50 = cos(pkin(6));
t51 = sin(qJ(4));
t89 = t50 * t51;
t22 = t49 * t53 - t89;
t21 = t49 * t51 + t50 * t53;
t15 = t21 * qJD(1);
t75 = qJD(1) * t89;
t59 = qJD(4) * t75 - t21 * qJDD(1);
t91 = t50 * pkin(2);
t88 = t50 * t52;
t87 = -pkin(5) + qJ(2);
t48 = t50 ^ 2;
t86 = t94 * t48;
t85 = t54 * pkin(1) + t52 * qJ(2);
t84 = qJD(1) * t49;
t83 = qJD(3) * t49;
t82 = qJDD(1) * pkin(1);
t69 = t49 * qJ(3) + pkin(1);
t24 = -t69 - t91;
t81 = qJDD(1) * t24;
t42 = t49 * qJDD(1);
t79 = t50 * qJDD(1);
t76 = qJD(1) * qJD(3);
t74 = t53 * t84;
t19 = t93 * t49 + qJDD(3);
t72 = -t52 * pkin(1) + t54 * qJ(2);
t47 = t49 ^ 2;
t56 = qJD(1) ^ 2;
t25 = (-t47 - t48) * t56;
t70 = t49 * t76;
t41 = qJDD(2) - t82;
t67 = g(1) * t54 + g(2) * t52;
t66 = t53 * t42 - t51 * t79;
t26 = t87 * t49;
t27 = t87 * t50;
t65 = t53 * t26 - t51 * t27;
t64 = t51 * t26 + t53 * t27;
t63 = -t41 + t82 - t46;
t58 = qJDD(2) + t81;
t4 = t58 - t70;
t62 = -t4 - t81 - t46;
t61 = -t67 + (t73 + 2 * t78) * t48;
t13 = t21 * qJD(4);
t18 = (pkin(2) + pkin(3)) * t50 + t69;
t60 = t93 * t47;
t55 = qJD(4) ^ 2;
t36 = g(1) * t88;
t29 = qJ(2) * t84 + qJD(3);
t17 = t74 - t75;
t14 = t22 * qJD(4);
t12 = t24 * qJD(1) + qJD(2);
t11 = t21 * t54;
t10 = t22 * t54;
t9 = t21 * t52;
t8 = t22 * t52;
t7 = (t87 * qJDD(1) + t77) * t50;
t6 = -pkin(5) * t42 + t19;
t5 = t18 * qJD(1) - qJD(2);
t3 = t18 * qJDD(1) - qJDD(2) + t70;
t2 = qJD(4) * t74 - t59;
t1 = -qJD(1) * t13 + t66;
t16 = [qJDD(1), t71, t67, t63 * t50 + t36, (-t63 - t92) * t49, 0.2e1 * t60 + t61, -t41 * pkin(1) - g(1) * t72 - g(2) * t85 + t94 * t47 + t86, t36 + (t62 + t70) * t50, t19 * t49 + t60 + t61, t47 * t76 + (t62 + t92) * t49, t4 * t24 - g(1) * (-pkin(2) * t88 + t72) - g(2) * (t54 * t91 + t85) + (t19 * qJ(2) + t71 * qJ(3) + t29 * qJD(2) - t12 * qJD(3)) * t49 + t86, t1 * t22 - t17 * t13, -t1 * t21 + t13 * t15 - t17 * t14 - t22 * t2, -t13 * qJD(4) + t22 * qJDD(4), -t14 * qJD(4) - t21 * qJDD(4), 0, t15 * t83 + t18 * t2 + t3 * t21 + t5 * t14 + t65 * qJDD(4) + g(1) * t9 - g(2) * t11 + (t22 * qJD(2) - t64 * qJD(4)) * qJD(4), t17 * t83 + t18 * t1 + t3 * t22 - t5 * t13 - t64 * qJDD(4) + g(1) * t8 - g(2) * t10 + (-t21 * qJD(2) - t65 * qJD(4)) * qJD(4); 0, 0, 0, -t79, t42, t25, qJ(2) * t25 + t41 - t71, -t79, t25, -t42, -t48 * t56 * qJ(2) + (-qJD(3) - t29) * t84 + t58 - t71, 0, 0, 0, 0, 0, (-t17 - t74) * qJD(4) + t59, 0.2e1 * t15 * qJD(4) - t66; 0, 0, 0, 0, 0, 0, 0, -t49 * t56 * t50, t42, -t47 * t56, g(3) * t50 + (qJD(1) * t12 - t67) * t49 + t19, 0, 0, 0, 0, 0, qJDD(4) * t53 - t15 * t84 - t55 * t51, -qJDD(4) * t51 - t17 * t84 - t55 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t15, -t15 ^ 2 + t17 ^ 2, t66, (t17 - t74) * qJD(4) + t59, qJDD(4), -g(1) * t10 - g(2) * t8 + g(3) * t21 - t5 * t17 - t51 * t7 + t53 * t6, g(1) * t11 + g(2) * t9 + g(3) * t22 + t5 * t15 - t51 * t6 - t53 * t7;];
tau_reg = t16;
