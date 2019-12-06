% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:12
% EndTime: 2019-12-05 15:09:15
% DurationCPUTime: 0.67s
% Computational Cost: add. (683->142), mult. (1453->178), div. (0->0), fcn. (1089->10), ass. (0->100)
t55 = sin(pkin(8));
t57 = cos(pkin(8));
t60 = sin(qJ(3));
t62 = cos(qJ(3));
t132 = -t55 * t60 + t57 * t62;
t25 = t132 * qJD(1);
t101 = qJD(3) * t132;
t30 = t55 * t62 + t57 * t60;
t128 = t30 * qJDD(1);
t135 = qJDD(3) * pkin(6) + qJD(1) * t101 + qJD(2) * qJD(4) + t128;
t28 = t30 * qJD(3);
t59 = sin(qJ(4));
t61 = cos(qJ(4));
t75 = pkin(4) * t61 + qJ(5) * t59 + pkin(3);
t134 = t75 * qJD(3);
t26 = t30 * qJD(1);
t19 = qJD(3) * pkin(6) + t26;
t120 = t19 * t59;
t13 = qJD(2) * t61 - t120;
t133 = qJD(5) - t13;
t52 = pkin(8) + qJ(3);
t47 = cos(t52);
t121 = g(3) * t47;
t46 = sin(t52);
t56 = sin(pkin(7));
t58 = cos(pkin(7));
t79 = g(1) * t58 + g(2) * t56;
t129 = t79 * t46;
t69 = -t121 + t129;
t131 = t75 * qJDD(3);
t78 = pkin(4) * t59 - qJ(5) * t61;
t24 = t78 * qJD(4) - t59 * qJD(5);
t104 = qJD(3) * t24;
t74 = -qJD(1) * t28 + t132 * qJDD(1);
t5 = t104 - t74 - t131;
t130 = -t5 + t131;
t10 = -qJD(4) * pkin(4) + t133;
t119 = t19 * t61;
t14 = qJD(2) * t59 + t119;
t11 = qJD(4) * qJ(5) + t14;
t98 = qJDD(4) * pkin(4);
t127 = qJDD(5) - t98;
t122 = g(3) * t46;
t126 = t79 * t47 + t122;
t63 = qJD(4) ^ 2;
t125 = pkin(6) * t63;
t116 = t56 * t59;
t115 = t56 * t61;
t112 = t58 * t59;
t111 = t58 * t61;
t110 = t59 * t61;
t109 = t24 - t26;
t53 = t59 ^ 2;
t54 = t61 ^ 2;
t108 = t53 - t54;
t107 = t53 + t54;
t106 = qJD(3) * pkin(3);
t105 = pkin(6) * qJDD(4);
t103 = qJD(3) * t26;
t102 = qJD(3) * t101;
t100 = qJD(3) * t59;
t95 = t61 * qJDD(3);
t92 = qJD(3) * qJD(4);
t91 = qJDD(4) * qJ(5);
t90 = t59 * qJDD(2) + t135 * t61;
t64 = qJD(3) ^ 2;
t89 = t64 * t110;
t88 = -g(1) * t56 + g(2) * t58;
t87 = t13 + t120;
t86 = qJD(4) * t119 - t61 * qJDD(2) + t135 * t59;
t18 = -t25 - t106;
t85 = t18 - t106;
t84 = t59 * t25 * qJD(4) + t61 * t103 + (g(1) * t111 + g(2) * t115) * t46;
t12 = -t25 - t134;
t83 = t12 - t134;
t82 = t107 * qJDD(3);
t80 = t121 + t125;
t77 = t10 * t59 + t11 * t61;
t76 = -qJD(3) * t28 + qJDD(3) * t132;
t21 = t47 * t115 - t112;
t23 = t47 * t111 + t116;
t73 = g(1) * t23 + g(2) * t21 - t90;
t72 = t30 * t63 - t76;
t71 = 0.2e1 * qJDD(3) * pkin(3) + t74 - t80;
t70 = -0.2e1 * t101 * qJD(4) - qJDD(4) * t30;
t20 = t47 * t116 + t111;
t22 = t47 * t112 - t115;
t68 = g(1) * t22 + g(2) * t20 + t59 * t122 - t86;
t67 = qJD(4) * t14 + t68;
t3 = t91 + (qJD(5) - t120) * qJD(4) + t90;
t4 = t86 + t127;
t66 = t3 * t61 + t4 * t59 + (t10 * t61 - t11 * t59) * qJD(4);
t65 = -t126 + t66;
t48 = t59 * qJDD(3);
t36 = qJDD(4) * t61 - t59 * t63;
t35 = qJDD(4) * t59 + t61 * t63;
t31 = t78 * qJD(3);
t2 = t70 * t59 - t72 * t61;
t1 = t72 * t59 + t70 * t61;
t6 = [qJDD(1) - g(3), -g(3) + (t55 ^ 2 + t57 ^ 2) * qJDD(1), 0, t76, -qJDD(3) * t30 - t102, 0, 0, 0, 0, 0, t2, t1, t2, t107 * t102 + t30 * t82, -t1, t101 * t77 + t12 * t28 - t132 * t5 + t66 * t30 - g(3); 0, qJDD(2) + t88, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, t36, 0, t35, t77 * qJD(4) + t3 * t59 - t4 * t61 + t88; 0, 0, qJDD(3), t69 + t74 + t103, -t128 + t126, qJDD(3) * t53 + 0.2e1 * t92 * t110, -0.2e1 * t108 * t92 + 0.2e1 * t59 * t95, t35, t36, 0, (t85 * qJD(4) - t105) * t59 + t71 * t61 + t84, (-t105 + (t25 + t85) * qJD(4)) * t61 + (-t103 - t71 - t129) * t59, (t83 * qJD(4) - t105) * t59 + (-t80 - t104 + t130) * t61 + t84, -t107 * t25 * qJD(3) + pkin(6) * t82 + t65, (t105 + (-t25 - t83) * qJD(4)) * t61 + (-t109 * qJD(3) - t125 + t130 + t69) * t59, t65 * pkin(6) + t109 * t12 - t77 * t25 + (-t5 + t69) * t75; 0, 0, 0, 0, 0, -t89, t108 * t64, t48, t95, qJDD(4), -t18 * t100 + t67, (-qJD(3) * t18 + t122) * t61 + t87 * qJD(4) + t73, 0.2e1 * t98 - qJDD(5) + (-t12 * t59 + t31 * t61) * qJD(3) + t67, -t78 * qJDD(3), -t61 * t122 + 0.2e1 * t91 + (t12 * t61 + t31 * t59) * qJD(3) + (0.2e1 * qJD(5) - t87) * qJD(4) - t73, t3 * qJ(5) - t4 * pkin(4) - t12 * t31 - t10 * t14 - g(1) * (-pkin(4) * t22 + qJ(5) * t23) - g(2) * (-pkin(4) * t20 + qJ(5) * t21) + t78 * t122 + t133 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t89, t48, -t53 * t64 - t63, -qJD(4) * t11 + t12 * t100 + t127 - t68;];
tau_reg = t6;
