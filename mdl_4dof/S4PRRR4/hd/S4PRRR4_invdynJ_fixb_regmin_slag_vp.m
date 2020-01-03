% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:41
% EndTime: 2019-12-31 16:32:42
% DurationCPUTime: 0.38s
% Computational Cost: add. (398->104), mult. (848->157), div. (0->0), fcn. (585->8), ass. (0->75)
t51 = pkin(7) + qJ(2);
t44 = sin(t51);
t45 = cos(t51);
t73 = g(1) * t45 + g(2) * t44;
t57 = sin(qJ(3));
t59 = cos(qJ(3));
t96 = pkin(5) + pkin(6);
t77 = qJD(2) * t96;
t14 = t59 * qJD(1) - t57 * t77;
t52 = qJD(3) + qJD(4);
t58 = cos(qJ(4));
t88 = qJD(2) * t59;
t78 = t58 * t88;
t56 = sin(qJ(4));
t89 = qJD(2) * t57;
t79 = t56 * t89;
t17 = -t78 + t79;
t19 = -t56 * t88 - t58 * t89;
t93 = t19 * t17;
t32 = t96 * t59;
t85 = t57 * qJD(1);
t15 = qJD(2) * t32 + t85;
t92 = t58 * t15;
t53 = t57 ^ 2;
t91 = -t59 ^ 2 + t53;
t90 = qJD(3) * pkin(3);
t87 = qJD(4) * t56;
t86 = (qJDD(2) * pkin(2));
t84 = qJDD(1) - g(3);
t83 = t57 * qJDD(2);
t82 = t59 * qJDD(2);
t81 = qJD(2) * qJD(3);
t80 = t57 * t90;
t31 = t96 * t57;
t43 = -t59 * pkin(3) - pkin(2);
t76 = qJD(3) * t96;
t75 = t59 * t81;
t72 = g(1) * t44 - g(2) * t45;
t71 = t56 * t83 - t58 * t82;
t20 = t56 * t57 - t58 * t59;
t10 = t52 * t20;
t21 = t56 * t59 + t58 * t57;
t50 = qJDD(3) + qJDD(4);
t70 = -t10 * t52 + t21 * t50;
t13 = t14 + t90;
t69 = -t56 * t13 - t92;
t68 = -t58 * t31 - t56 * t32;
t67 = -t56 * t31 + t58 * t32;
t66 = -0.2e1 * pkin(2) * t81 - pkin(5) * qJDD(3);
t60 = qJD(3) ^ 2;
t65 = -pkin(5) * t60 + t72 + (2 * t86);
t61 = qJD(2) ^ 2;
t64 = t61 * pkin(2) - qJDD(2) * pkin(5) + t73;
t4 = qJD(4) * t78 - t52 * t79 + t56 * t82 + (t75 + t83) * t58;
t30 = t43 * qJD(2);
t55 = qJ(3) + qJ(4);
t48 = sin(t55);
t49 = cos(t55);
t46 = t59 * qJDD(1);
t8 = qJDD(3) * pkin(3) + t46 - qJDD(2) * t31 + (-t59 * t77 - t85) * qJD(3);
t63 = t30 * t17 + t15 * t87 + g(3) * t48 + (-t15 * t52 - t8) * t56 + t73 * t49;
t11 = t52 * t21;
t9 = t14 * qJD(3) + t57 * qJDD(1) + qJDD(2) * t32;
t62 = -g(3) * t49 + t69 * qJD(4) + t30 * t19 + t73 * t48 - t56 * t9 + t58 * t8;
t29 = qJDD(3) * t59 - t60 * t57;
t28 = qJDD(3) * t57 + t60 * t59;
t23 = t59 * t76;
t22 = t57 * t76;
t16 = -t86 + (t57 * t81 - t82) * pkin(3);
t6 = -t17 ^ 2 + t19 ^ 2;
t5 = t11 * qJD(2) + t71;
t3 = -t11 * t52 - t20 * t50;
t2 = -t71 + (-qJD(2) * t21 - t19) * t52;
t1 = t17 * t52 + t4;
t7 = [t84, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, 0, 0, 0, 0, t3, -t70; 0, qJDD(2), t72, t73, t53 * qJDD(2) + 0.2e1 * t57 * t75, 0.2e1 * t57 * t82 - 0.2e1 * t91 * t81, t28, t29, 0, t66 * t57 + t65 * t59, -t65 * t57 + t66 * t59, t19 * t10 + t4 * t21, t10 * t17 + t19 * t11 - t4 * t20 - t21 * t5, t70, t3, 0, t17 * t80 + t43 * t5 + t16 * t20 + t30 * t11 + (-t67 * qJD(4) + t56 * t22 - t58 * t23) * t52 + t68 * t50 + t72 * t49, -t19 * t80 + t43 * t4 + t16 * t21 - t30 * t10 - (t68 * qJD(4) - t58 * t22 - t56 * t23) * t52 - t67 * t50 - t72 * t48; 0, 0, 0, 0, -t57 * t61 * t59, t91 * t61, t83, t82, qJDD(3), -g(3) * t59 + t64 * t57 + t46, -t84 * t57 + t64 * t59, -t93, t6, t1, t2, t50, -(-t56 * t14 - t92) * t52 + (-t17 * t89 + t58 * t50 - t52 * t87) * pkin(3) + t62, (-qJD(4) * t13 + t14 * t52 - t9) * t58 + (-qJD(4) * t58 * t52 + t19 * t89 - t56 * t50) * pkin(3) + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t6, t1, t2, t50, -t69 * t52 + t62, (-t9 + (-qJD(4) + t52) * t13) * t58 + t63;];
tau_reg = t7;
