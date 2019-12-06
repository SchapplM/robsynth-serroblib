% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:36:27
% EndTime: 2019-12-05 17:36:31
% DurationCPUTime: 0.90s
% Computational Cost: add. (857->191), mult. (1699->264), div. (0->0), fcn. (1071->10), ass. (0->109)
t58 = cos(pkin(7));
t44 = -pkin(1) * t58 - pkin(2);
t55 = sin(pkin(8));
t57 = cos(pkin(8));
t30 = -pkin(3) * t57 - pkin(6) * t55 + t44;
t20 = t30 * qJD(1) + qJD(3);
t56 = sin(pkin(7));
t43 = pkin(1) * t56 + qJ(3);
t35 = t43 * qJD(1);
t27 = qJD(2) * t55 + t35 * t57;
t60 = sin(qJ(4));
t62 = cos(qJ(4));
t105 = qJD(1) * t55;
t84 = qJ(5) * t105;
t7 = (t20 - t84) * t60 + t27 * t62;
t128 = qJD(4) * t7;
t95 = qJDD(1) * t44;
t34 = qJDD(3) + t95;
t52 = qJ(1) + pkin(7);
t48 = sin(t52);
t49 = cos(t52);
t75 = g(2) * t49 + g(3) * t48;
t126 = -t34 + t75;
t90 = qJD(1) * qJD(3);
t31 = qJDD(1) * t43 + t90;
t113 = t57 * t60;
t21 = t48 * t113 + t49 * t62;
t23 = t49 * t113 - t48 * t62;
t125 = -g(2) * t21 + g(3) * t23;
t104 = qJD(1) * t57;
t41 = -qJD(4) + t104;
t81 = t62 * t20 - t27 * t60;
t6 = -t62 * t84 + t81;
t3 = -pkin(4) * t41 + t6;
t124 = -t6 + t3;
t123 = pkin(4) * t60;
t118 = t31 * t55;
t46 = t57 * qJDD(2);
t18 = -t46 + t118;
t120 = t18 * t55;
t91 = t57 * qJDD(1);
t40 = -qJDD(4) + t91;
t117 = t40 * t57;
t116 = (pkin(4) * t62 + pkin(3)) * t57;
t115 = t49 * t60;
t50 = t55 ^ 2;
t64 = qJD(1) ^ 2;
t114 = t50 * t64;
t112 = t57 * t62;
t101 = qJD(3) * t62;
t99 = qJD(4) * t62;
t111 = t57 * t101 + t30 * t99;
t32 = t43 * t112;
t110 = t60 * t30 + t32;
t109 = t57 ^ 2 + t50;
t53 = t60 ^ 2;
t54 = t62 ^ 2;
t108 = -t53 - t54;
t107 = t53 - t54;
t106 = qJ(5) * t55;
t103 = qJD(1) * t60;
t102 = qJD(3) * t57;
t100 = qJD(4) * t60;
t98 = qJD(5) * t55;
t97 = -qJD(4) - t41;
t94 = qJDD(1) * t60;
t93 = qJDD(1) * t62;
t92 = t55 * qJDD(1);
t89 = qJD(1) * qJD(4);
t88 = qJD(1) * qJD(5);
t87 = t62 * t106;
t86 = t41 * t100;
t85 = t57 * t100;
t61 = sin(qJ(1));
t83 = -pkin(1) * t61 + t49 * qJ(3);
t82 = t62 * t89;
t80 = t40 - t91;
t79 = t40 + t91;
t78 = qJD(1) * t97;
t77 = t55 * pkin(4) * t82 + t92 * t123 + qJDD(5) - t46;
t63 = cos(qJ(1));
t76 = -pkin(1) * t63 - pkin(2) * t49;
t74 = g(2) * t48 - g(3) * t49;
t73 = g(2) * t63 + g(3) * t61;
t72 = t3 * t62 + t60 * t7;
t71 = t3 * t60 - t62 * t7;
t19 = qJDD(2) * t55 + t31 * t57;
t70 = t19 * t57 + t120;
t47 = t57 * qJD(2);
t26 = t35 * t55 - t47;
t69 = t26 * t55 + t27 * t57;
t17 = t30 * qJDD(1) + qJDD(3);
t68 = t27 * t100 - t60 * t17 - t62 * t19 - t20 * t99;
t67 = -t95 + t126;
t66 = -t41 ^ 2 - t114;
t33 = t85 * t105;
t29 = t62 * t30;
t24 = -t49 * t112 - t48 * t60;
t22 = t48 * t112 - t115;
t14 = t62 * t17;
t11 = qJD(5) - t47 + (pkin(4) * t103 + t35) * t55;
t10 = -t60 * t106 + t110;
t9 = t77 + t118;
t8 = -t87 + t29 + (-t43 * t60 - pkin(4)) * t57;
t5 = -t60 * t102 - t62 * t98 + (-t32 + (-t30 + t106) * t60) * qJD(4);
t4 = -t60 * t98 + (-t43 * t113 - t87) * qJD(4) + t111;
t2 = (-t60 * t88 + (-t82 - t94) * qJ(5)) * t55 - t68;
t1 = -pkin(4) * t40 - t60 * t19 + t14 + (-qJ(5) * qJDD(1) - t88) * t62 * t55 - t128;
t12 = [qJDD(1), t73, -g(2) * t61 + g(3) * t63, (t73 + (t56 ^ 2 + t58 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t67 * t57, -t67 * t55, t31 * t109 + t70 + t74, t34 * t44 - g(2) * (-qJ(3) * t48 + t76) - g(3) * (-pkin(2) * t48 + t83) + t70 * t43 + t69 * qJD(3), (qJDD(1) * t54 - 0.2e1 * t60 * t82) * t50, 0.2e1 * (t107 * t89 - t60 * t93) * t50, t33 + (-t62 * t79 + t86) * t55, (t79 * t60 + (t41 + t104) * t99) * t55, t117, -g(2) * t24 + g(3) * t22 - t14 * t57 - t29 * t40 + ((qJD(1) * t50 + t41 * t57) * t43 + t69) * t99 + (-(-qJD(4) * t30 - t102) * t41 - (-qJD(4) * t20 - t19) * t57 + t50 * t90 + t120 + (t50 * qJDD(1) + t117) * t43) * t60, (-t43 * t85 + t111) * t41 + t110 * t40 - t68 * t57 - g(2) * t23 - g(3) * t21 + (-t26 * t100 + t18 * t62) * t55 + (t43 * t93 + (-t43 * t100 + t101) * qJD(1)) * t50, ((-t128 - qJDD(1) * t8 - t1 + (-qJD(4) * t10 - t5) * qJD(1)) * t62 + (qJD(4) * t3 - qJDD(1) * t10 - t2 + (qJD(4) * t8 - t4) * qJD(1)) * t60 + t75) * t55, t2 * t10 + t7 * t4 + t1 * t8 + t3 * t5 - g(2) * (-t49 * t116 + t76) - g(3) * (pkin(4) * t115 + t83) + (-g(2) * (-qJ(3) - t123) - g(3) * (-pkin(2) - t116)) * t48 + (t9 * (t43 + t123) + t11 * (pkin(4) * t99 + qJD(3)) - t75 * (-qJ(5) - pkin(6))) * t55; 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, -t18 * t57 + t19 * t55 - g(1), 0, 0, 0, 0, 0, (t80 * t60 + (t41 - t104) * t99) * t55, t33 + (t62 * t80 - t86) * t55, 0, -t57 * t9 - g(1) + (-qJD(4) * t72 - t1 * t60 + t2 * t62) * t55; 0, 0, 0, 0, -t91, t92, -t109 * t64, -qJD(1) * t69 - t126, 0, 0, 0, 0, 0, -t40 * t62 + t60 * t66, t40 * t60 + t62 * t66, t108 * t92, t1 * t62 + t2 * t60 - t71 * qJD(4) + (-t11 * t55 + t57 * t71) * qJD(1) - t75; 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60 * t114, -t107 * t114, (t60 * t78 + t93) * t55, (t62 * t78 - t94) * t55, -t40, t14 + (-t26 * t105 + t97 * t27) * t62 + (g(1) * t55 + t97 * t20 - t19) * t60 + t125, -t81 * t41 - g(2) * t22 - g(3) * t24 + (g(1) * t62 + t26 * t103) * t55 + t68, (-pkin(4) * t93 + (pkin(4) * qJD(4) - t124) * t103) * t55, t124 * t7 + (t1 + (-t11 * t62 * qJD(1) + g(1) * t60) * t55 + t125) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108 * t114, g(1) * t57 + (qJD(1) * t72 + t31 + t74) * t55 + t77;];
tau_reg = t12;
