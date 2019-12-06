% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPR1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:53
% EndTime: 2019-12-05 16:15:56
% DurationCPUTime: 0.58s
% Computational Cost: add. (787->146), mult. (1092->181), div. (0->0), fcn. (739->12), ass. (0->95)
t80 = pkin(8) + qJ(2);
t73 = qJ(3) + t80;
t60 = cos(t73);
t132 = g(2) * t60;
t59 = sin(t73);
t55 = g(1) * t59;
t136 = t132 - t55;
t118 = pkin(2) * qJD(2);
t85 = sin(qJ(3));
t109 = t85 * t118;
t87 = cos(qJ(3));
t130 = t87 * pkin(2);
t120 = -qJD(3) * t109 + qJDD(2) * t130;
t106 = qJDD(4) - t120;
t76 = qJDD(2) + qJDD(3);
t131 = t76 * pkin(3);
t25 = t106 - t131;
t104 = -t25 - t132;
t83 = cos(pkin(9));
t86 = cos(qJ(5));
t124 = t86 * t83;
t82 = sin(pkin(9));
t84 = sin(qJ(5));
t127 = t84 * t82;
t34 = -t124 + t127;
t35 = t86 * t82 + t84 * t83;
t119 = t82 ^ 2 + t83 ^ 2;
t81 = qJD(2) + qJD(3);
t28 = t35 * t81;
t135 = g(1) * t60 + g(2) * t59;
t114 = qJDD(2) * t85;
t115 = qJD(3) * t87;
t20 = t76 * qJ(4) + t81 * qJD(4) + (qJD(2) * t115 + t114) * pkin(2);
t66 = t83 * qJDD(1);
t12 = -t82 * t20 + t66;
t13 = t82 * qJDD(1) + t83 * t20;
t98 = -t12 * t82 + t13 * t83;
t134 = t119 * t81;
t133 = t120 + t55;
t30 = t35 * qJD(5);
t117 = qJD(2) * t87;
t99 = -pkin(2) * t117 + qJD(4);
t128 = t83 * t76;
t123 = t104 * t82;
t122 = t60 * pkin(3) + t59 * qJ(4);
t116 = qJD(3) * t85;
t113 = t81 * t127;
t112 = t81 * t124;
t111 = qJD(5) * t112 + t35 * t76;
t110 = pkin(2) * t116;
t107 = t81 * t116;
t62 = -t83 * pkin(4) - pkin(3);
t105 = -t59 * pkin(3) + t60 * qJ(4);
t103 = t119 * t76;
t102 = t81 * t109;
t101 = t34 * t76;
t97 = t119 * (t81 * qJ(4) + t109);
t61 = t85 * pkin(2) + qJ(4);
t32 = (-pkin(7) - t61) * t82;
t74 = t83 * pkin(7);
t33 = t83 * t61 + t74;
t96 = t86 * t32 - t84 * t33;
t95 = t84 * t32 + t86 * t33;
t44 = (-pkin(7) - qJ(4)) * t82;
t45 = t83 * qJ(4) + t74;
t94 = t86 * t44 - t84 * t45;
t93 = t84 * t44 + t86 * t45;
t92 = -t135 + t98;
t18 = t62 * t76 + t106;
t24 = t62 * t81 + t99;
t29 = t34 * qJD(5);
t79 = pkin(9) + qJ(5);
t69 = sin(t79);
t91 = t136 * t69 + t18 * t35 - t24 * t29;
t71 = cos(t79);
t90 = -t136 * t71 + t18 * t34 + t24 * t30;
t89 = t102 - t132;
t63 = -pkin(3) - t130;
t88 = pkin(2) * t107 + t63 * t76;
t72 = cos(t80);
t70 = sin(t80);
t52 = pkin(2) * t115 + qJD(4);
t43 = t62 - t130;
t42 = t83 * t55;
t36 = -t81 * pkin(3) + t99;
t26 = -t112 + t113;
t15 = -t30 * qJD(5) - t34 * qJDD(5);
t14 = -t29 * qJD(5) + t35 * qJDD(5);
t6 = pkin(7) * t128 + t13;
t5 = t66 + (-pkin(7) * t76 - t20) * t82;
t4 = t81 * t30 + t101;
t3 = -qJD(5) * t113 + t111;
t2 = -t28 * t29 + t3 * t35;
t1 = t29 * t26 - t28 * t30 - t3 * t34 - t35 * t4;
t7 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t83 + t13 * t82 - g(3), 0, 0, 0, 0, 0, t15, -t14; 0, qJDD(2), g(1) * t70 - g(2) * t72, g(1) * t72 + g(2) * t70, t76, -t132 + (t76 * t87 - t107) * pkin(2) + t133, ((-qJDD(2) - t76) * t85 + (-qJD(2) - t81) * t115) * pkin(2) + t135, t42 + (t104 - t88) * t83, (t88 - t55) * t82 - t123, t61 * t103 + t52 * t134 + t92, t25 * t63 + t36 * t110 - g(1) * (-pkin(2) * t70 + t105) - g(2) * (pkin(2) * t72 + t122) + t98 * t61 + t97 * t52, t2, t1, t14, t15, 0, t26 * t110 + t43 * t4 + t96 * qJDD(5) + (-t95 * qJD(5) - t35 * t52) * qJD(5) + t90, t28 * t110 + t43 * t3 - t95 * qJDD(5) + (-t96 * qJD(5) + t34 * t52) * qJD(5) + t91; 0, 0, 0, 0, t76, t89 + t133, (-t114 + (-qJD(3) + t81) * t117) * pkin(2) + t135, t42 + (-t25 + t89 + t131) * t83, (-t102 - t131 - t55) * t82 - t123, qJ(4) * t103 + t99 * t134 + t92, -t25 * pkin(3) - g(1) * t105 - g(2) * t122 + t97 * qJD(4) + t98 * qJ(4) + (-t36 * t85 - t97 * t87) * t118, t2, t1, t14, t15, 0, t62 * t4 + t94 * qJDD(5) + (-t35 * qJD(4) - t93 * qJD(5)) * qJD(5) + (-t85 * t26 + t87 * t30) * t118 + t90, t62 * t3 - t93 * qJDD(5) + (t34 * qJD(4) - t94 * qJD(5)) * qJD(5) + (-t85 * t28 - t87 * t29) * t118 + t91; 0, 0, 0, 0, 0, 0, 0, -t128, t82 * t76, -t119 * t81 ^ 2, -t97 * t81 - t104 - t55, 0, 0, 0, 0, 0, 0.2e1 * t28 * qJD(5) + t101, (-t26 - t113) * qJD(5) + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t26, -t26 ^ 2 + t28 ^ 2, (t26 - t113) * qJD(5) + t111, -t101, qJDD(5), -g(3) * t71 + t135 * t69 - t24 * t28 + t86 * t5 - t84 * t6, g(3) * t69 + t135 * t71 + t24 * t26 - t84 * t5 - t86 * t6;];
tau_reg = t7;
