% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRR1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:18
% EndTime: 2019-12-05 17:38:21
% DurationCPUTime: 0.60s
% Computational Cost: add. (565->158), mult. (1010->195), div. (0->0), fcn. (625->8), ass. (0->102)
t74 = cos(qJ(4));
t106 = t74 * qJDD(1);
t104 = qJD(1) * qJD(4);
t71 = sin(qJ(4));
t99 = t71 * t104;
t131 = t99 - t106;
t72 = sin(qJ(1));
t75 = cos(qJ(1));
t94 = g(1) * t75 + g(2) * t72;
t70 = sin(qJ(5));
t73 = cos(qJ(5));
t25 = t70 * t74 + t73 * t71;
t22 = t25 * qJD(1);
t59 = qJD(4) + qJD(5);
t121 = t22 * t59;
t130 = g(1) * t72 - g(2) * t75;
t69 = pkin(1) + qJ(3);
t129 = qJD(1) * t69;
t37 = -qJD(2) + t129;
t68 = -pkin(6) + qJ(2);
t128 = (qJD(2) + t37 + t129) * qJD(4) + qJDD(4) * t68;
t62 = qJD(1) * qJD(2);
t54 = 0.2e1 * t62;
t125 = pkin(7) - t68;
t10 = t59 * t25;
t26 = -t70 * t71 + t73 * t74;
t58 = qJDD(4) + qJDD(5);
t124 = -t10 * t59 + t26 * t58;
t115 = qJD(1) * t71;
t102 = t70 * t115;
t114 = qJD(1) * t74;
t21 = -t73 * t114 + t102;
t123 = t21 * t22;
t122 = t21 * t59;
t47 = qJ(2) * qJD(1) + qJD(3);
t36 = -pkin(6) * qJD(1) + t47;
t17 = -pkin(7) * t115 + t71 * t36;
t120 = t73 * t17;
t119 = t75 * pkin(1) + t72 * qJ(2);
t65 = t74 ^ 2;
t117 = t71 ^ 2 - t65;
t76 = qJD(4) ^ 2;
t77 = qJD(1) ^ 2;
t116 = -t76 - t77;
t113 = qJD(4) * t71;
t112 = qJD(4) * t74;
t111 = qJD(5) * t70;
t109 = qJDD(4) * t71;
t108 = t69 * qJDD(1);
t107 = t71 * qJDD(1);
t67 = qJDD(1) * pkin(1);
t105 = t67 - qJDD(2);
t103 = qJD(3) * qJD(1);
t61 = qJDD(1) * qJ(2);
t31 = t125 * t74;
t101 = qJDD(2) - t130;
t100 = qJDD(3) + t61 + t62;
t98 = t74 * t104;
t97 = t59 * t74;
t96 = -0.2e1 * t98;
t95 = -t67 + t101;
t43 = t71 * pkin(4) + t69;
t18 = -pkin(7) * t114 + t74 * t36;
t92 = t73 * t106 - t70 * t107;
t11 = -t71 * t111 - t70 * t113 + t73 * t97;
t91 = -t11 * t59 - t25 * t58;
t14 = qJD(4) * pkin(4) + t18;
t90 = -t70 * t14 - t120;
t30 = t125 * t71;
t89 = -t73 * t30 - t70 * t31;
t88 = t70 * t30 - t73 * t31;
t60 = qJDD(1) * qJ(3);
t86 = -t60 + t95;
t85 = -qJD(5) * t102 - t131 * t70;
t28 = t60 + t103 + t105;
t84 = 0.2e1 * t61 + t54 - t94;
t83 = t98 + t107;
t82 = t37 * qJD(1) + t94;
t24 = t43 * qJD(1) - qJD(2);
t66 = qJ(4) + qJ(5);
t49 = sin(t66);
t50 = cos(t66);
t27 = -pkin(6) * qJDD(1) + t100;
t20 = t74 * t27;
t7 = qJDD(4) * pkin(4) + t131 * pkin(7) - t36 * t113 + t20;
t81 = t17 * t111 + g(3) * t50 + (-t17 * t59 - t7) * t70 + t24 * t22 + t94 * t49;
t80 = -t68 * t76 + t103 + t108 + t130 + t28;
t79 = (-t59 * t114 - t107) * t73 - t85;
t8 = -t83 * pkin(7) + t36 * t112 + t71 * t27;
t78 = g(3) * t49 + t90 * qJD(5) + t24 * t21 - t94 * t50 + t73 * t7 - t70 * t8;
t3 = -t10 * qJD(1) + t92;
t52 = t75 * qJ(2);
t48 = qJDD(4) * t74;
t38 = pkin(4) * t112 + qJD(3);
t16 = t71 * qJD(2) - qJD(4) * t31;
t15 = t74 * qJD(2) + t125 * t113;
t12 = t83 * pkin(4) + t28;
t6 = t21 ^ 2 - t22 ^ 2;
t4 = (qJD(1) * t97 + t107) * t73 + t85;
t2 = t79 - t122;
t1 = t3 + t121;
t5 = [qJDD(1), t130, t94, -0.2e1 * t67 + t101, t84, t105 * pkin(1) - g(1) * (-t72 * pkin(1) + t52) - g(2) * t119 + (t61 + t54) * qJ(2), qJDD(3) + t84, 0.2e1 * t103 - t86 + t108, t28 * t69 + t37 * qJD(3) + t100 * qJ(2) + t47 * qJD(2) - g(1) * (-t69 * t72 + t52) - g(2) * (t75 * qJ(3) + t119), t65 * qJDD(1) + t71 * t96, 0.2e1 * t117 * t104 - 0.2e1 * t71 * t106, -t76 * t71 + t48, -t76 * t74 - t109, 0, t128 * t74 + t80 * t71, -t128 * t71 + t80 * t74, t21 * t10 + t3 * t26, t10 * t22 + t21 * t11 - t3 * t25 - t26 * t4, t124, t91, 0, t38 * t22 + t43 * t4 + t12 * t25 + t24 * t11 + (-t89 * qJD(5) + t73 * t15 - t70 * t16) * t59 + t88 * t58 + t130 * t49, -t38 * t21 + t43 * t3 + t12 * t26 - t24 * t10 - (t88 * qJD(5) + t70 * t15 + t73 * t16) * t59 - t89 * t58 + t130 * t50; 0, 0, 0, qJDD(1), -t77, -t77 * qJ(2) + t95, -t77, -qJDD(1), (-qJD(3) - t47) * qJD(1) + t86, 0, 0, 0, 0, 0, t96 - t107, 0.2e1 * t99 - t106, 0, 0, 0, 0, 0, t79 + t122, 0.2e1 * t121 - t92; 0, 0, 0, 0, 0, 0, qJDD(1), -t77, -t82 + t100, 0, 0, 0, 0, 0, t116 * t71 + t48, t116 * t74 - t109, 0, 0, 0, 0, 0, -qJD(1) * t22 + t124, qJD(1) * t21 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t77 * t71, -t117 * t77, t106, -t107, qJDD(4), g(3) * t71 - t82 * t74 + t20, g(3) * t74 + (-t27 + t82) * t71, -t123, t6, t1, t2, t58, -(-t70 * t18 - t120) * t59 + (-t59 * t111 - t22 * t114 + t73 * t58) * pkin(4) + t78, (-qJD(5) * t14 + t18 * t59 - t8) * t73 + (-qJD(5) * t73 * t59 + t21 * t114 - t70 * t58) * pkin(4) + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, t6, t1, t2, t58, -t90 * t59 + t78, (-t8 + (-qJD(5) + t59) * t14) * t73 + t81;];
tau_reg = t5;
