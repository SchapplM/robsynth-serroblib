% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:38
% EndTime: 2019-12-05 18:14:41
% DurationCPUTime: 0.55s
% Computational Cost: add. (1051->125), mult. (1860->158), div. (0->0), fcn. (1095->14), ass. (0->94)
t58 = qJ(1) + pkin(9) + qJ(3);
t52 = qJ(4) + t58;
t47 = sin(t52);
t48 = cos(t52);
t126 = g(2) * t48 + g(3) * t47;
t59 = qJDD(1) + qJDD(3);
t63 = sin(pkin(9));
t119 = pkin(1) * t63;
t67 = sin(qJ(3));
t64 = cos(pkin(9));
t51 = t64 * pkin(1) + pkin(2);
t39 = t51 * qJD(1);
t106 = qJD(3) * t39;
t37 = t51 * qJDD(1);
t71 = cos(qJ(3));
t89 = -t67 * t106 + t71 * t37;
t98 = qJD(1) * qJD(3) * t71;
t76 = (-qJDD(1) * t67 - t98) * t119 + t89;
t13 = t59 * pkin(3) + t76;
t100 = qJD(1) * t119;
t94 = t67 * t100;
t16 = (qJDD(1) * t119 + t106) * t71 - qJD(3) * t94 + t67 * t37;
t22 = t71 * t39 - t94;
t60 = qJD(1) + qJD(3);
t20 = t60 * pkin(3) + t22;
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t105 = qJD(4) * t66;
t23 = t71 * t100 + t67 * t39;
t93 = -g(2) * t47 + g(3) * t48 + t23 * t105;
t78 = -(qJD(4) * t20 + t16) * t70 - t66 * t13 + t93;
t40 = t71 * t51;
t124 = -t67 * t119 + t40;
t29 = pkin(3) + t124;
t30 = t71 * t119 + t67 * t51;
t108 = t66 * t29 + t70 * t30;
t49 = sin(t58);
t50 = cos(t58);
t125 = g(2) * t50 + g(3) * t49;
t109 = t70 * t23;
t11 = t66 * t20 + t109;
t123 = t11 * qJD(4);
t14 = t66 * t22 + t109;
t53 = t66 * pkin(3) + pkin(8);
t114 = t70 * pkin(3);
t54 = -pkin(4) - t114;
t56 = qJDD(4) + t59;
t57 = qJD(4) + t60;
t73 = qJD(5) ^ 2;
t122 = -t53 * t73 - t54 * t56 - (pkin(3) * t105 - t14) * t57;
t69 = cos(qJ(5));
t104 = t69 * qJD(5);
t116 = t56 * pkin(4);
t97 = -t70 * t13 + t66 * t16;
t3 = -t116 + t97 + t123;
t65 = sin(qJ(5));
t111 = t66 * t23;
t10 = t70 * t20 - t111;
t115 = t57 * pkin(4);
t8 = -t10 - t115;
t120 = t8 * t104 + t3 * t65;
t27 = t124 * qJD(3);
t28 = t30 * qJD(3);
t117 = (t108 * qJD(4) + t66 * t27 + t70 * t28) * t57;
t113 = t11 * t57;
t110 = t69 * t56;
t61 = t65 ^ 2;
t107 = -t69 ^ 2 + t61;
t103 = qJDD(2) - g(1);
t101 = t8 * qJD(5) * t65 + t126 * t69;
t99 = -pkin(3) * t57 - t20;
t68 = sin(qJ(1));
t72 = cos(qJ(1));
t90 = g(2) * t72 + g(3) * t68;
t88 = t70 * t29 - t66 * t30;
t86 = -t97 + t126;
t84 = -pkin(8) * t73 + t113 + t116;
t17 = -pkin(4) - t88;
t18 = pkin(8) + t108;
t83 = -t17 * t56 - t18 * t73 - t117;
t82 = -t56 * pkin(8) - t8 * t57 + t78;
t81 = -pkin(8) * qJDD(5) + (t10 - t115) * qJD(5);
t4 = t88 * qJD(4) + t70 * t27 - t66 * t28;
t80 = -qJDD(5) * t18 + (t17 * t57 - t4) * qJD(5);
t15 = t70 * t22 - t111;
t79 = -qJDD(5) * t53 + (-qJD(4) * t114 + t54 * t57 + t15) * qJD(5);
t77 = t86 - t123;
t75 = -g(2) * t49 + g(3) * t50 - t16;
t55 = t57 ^ 2;
t36 = qJDD(5) * t69 - t73 * t65;
t35 = qJDD(5) * t65 + t73 * t69;
t24 = 0.2e1 * t65 * t57 * t104 + t61 * t56;
t19 = -0.2e1 * t107 * t57 * qJD(5) + 0.2e1 * t65 * t110;
t1 = [qJDD(1), t90, -g(2) * t68 + g(3) * t72, (t90 + (t63 ^ 2 + t64 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t59, -t28 * t60 + t40 * t59 + (-t98 + (-qJDD(1) - t59) * t67) * t119 + t89 + t125, -t27 * t60 - t30 * t59 + t75, t56, t88 * t56 - t117 + t77, -t108 * t56 - t4 * t57 + t78, t24, t19, t35, t36, 0, t80 * t65 + (-t3 + t83) * t69 + t101, t80 * t69 + (-t83 - t126) * t65 + t120; 0, 0, 0, t103, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35; 0, 0, 0, 0, t59, t23 * t60 + t125 + t76, t22 * t60 + t75, t56, t56 * t114 + t14 * t57 + (t99 * t66 - t109) * qJD(4) + t86, t15 * t57 + (-pkin(3) * t56 - t13) * t66 + (t99 * qJD(4) - t16) * t70 + t93, t24, t19, t35, t36, 0, t79 * t65 + (t122 - t3) * t69 + t101, t79 * t69 + (-t122 - t126) * t65 + t120; 0, 0, 0, 0, 0, 0, 0, t56, t77 + t113, t10 * t57 + t78, t24, t19, t35, t36, 0, t81 * t65 + (-t3 + t84) * t69 + t101, t81 * t69 + (-t84 - t126) * t65 + t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 * t55 * t69, t107 * t55, t65 * t56, t110, qJDD(5), t103 * t69 + t82 * t65, -t103 * t65 + t82 * t69;];
tau_reg = t1;
