% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:13:15
% EndTime: 2021-01-15 17:13:19
% DurationCPUTime: 0.74s
% Computational Cost: add. (815->192), mult. (1386->219), div. (0->0), fcn. (777->6), ass. (0->115)
t70 = sin(pkin(7));
t71 = cos(pkin(7));
t75 = sin(qJ(1));
t77 = cos(qJ(1));
t31 = -t75 * t70 - t77 * t71;
t32 = t77 * t70 - t75 * t71;
t146 = -g(1) * t31 - g(2) * t32;
t109 = qJD(1) * qJD(4);
t102 = qJ(5) * t109;
t108 = qJD(1) * qJD(5);
t111 = qJ(5) * qJDD(1);
t74 = sin(qJ(4));
t76 = cos(qJ(4));
t145 = t76 * t102 + (t108 + t111) * t74;
t140 = t76 * pkin(4);
t105 = pkin(3) + t140;
t73 = qJ(5) + pkin(6);
t144 = t105 * t70 - t73 * t71 + qJ(2);
t78 = pkin(1) + pkin(2);
t128 = qJD(4) * pkin(4);
t118 = qJ(5) * qJD(1);
t117 = qJD(1) * qJ(2);
t41 = -t78 * qJD(1) + qJD(2);
t23 = t71 * t117 + t70 * t41;
t18 = -qJD(1) * pkin(6) + t23;
t61 = t76 * qJD(3);
t9 = t61 + (t118 - t18) * t74;
t6 = t9 + t128;
t143 = -t9 + t6;
t64 = g(3) * t76;
t136 = t76 * t18;
t89 = -t74 * qJD(3) - t136;
t10 = -t76 * t118 - t89;
t139 = t10 * t76;
t138 = t31 * t74;
t137 = t32 * t74;
t80 = qJD(1) ^ 2;
t135 = t76 * t80;
t134 = g(1) * t137 - g(2) * t138;
t133 = t71 * qJ(2) - t70 * t78;
t132 = g(1) * t75 - g(2) * t77;
t67 = t74 ^ 2;
t68 = t76 ^ 2;
t131 = t67 - t68;
t130 = t67 + t68;
t79 = qJD(4) ^ 2;
t129 = t79 + t80;
t34 = -pkin(6) + t133;
t127 = qJ(5) - t34;
t126 = qJD(1) * t76;
t125 = qJD(2) * t71;
t124 = qJDD(1) * pkin(1);
t123 = qJDD(4) * pkin(4);
t22 = -t70 * t117 + t71 * t41;
t17 = qJD(1) * pkin(3) - t22;
t13 = pkin(4) * t126 + qJD(5) + t17;
t122 = t13 * qJD(1);
t121 = t17 * qJD(1);
t120 = t74 * qJD(4);
t119 = qJDD(3) + g(3);
t116 = qJDD(4) * t74;
t115 = qJDD(4) * t76;
t114 = t74 * qJDD(1);
t113 = t76 * qJDD(1);
t112 = qJ(2) * qJDD(1);
t110 = qJD(1) * qJD(2);
t14 = t18 * t120;
t107 = t146 * t76 + t14;
t40 = -t78 * qJDD(1) + qJDD(2);
t16 = t70 * t40 + (t110 + t112) * t71;
t106 = 0.2e1 * t110;
t44 = t70 * t110;
t104 = t74 * t109;
t103 = -t70 * qJ(2) - t71 * t78;
t60 = t76 * qJDD(3);
t100 = -g(1) * t138 - g(2) * t137 + t60 + t64;
t33 = pkin(3) - t103;
t30 = t33 + t140;
t99 = -qJD(1) * t30 - t13;
t98 = qJD(4) * t127;
t97 = 0.2e1 * t76 * t109;
t96 = -qJD(5) + t125;
t12 = -qJDD(1) * pkin(6) + t16;
t95 = qJD(4) * qJD(3) + t12;
t94 = qJDD(2) - t124;
t15 = -t70 * t112 + t71 * t40 - t44;
t93 = -g(1) * t32 + g(2) * t31;
t92 = g(1) * t77 + g(2) * t75;
t91 = -t6 * t74 + t139;
t90 = t22 * t70 - t23 * t71;
t11 = qJDD(1) * pkin(3) - t15;
t88 = qJDD(3) + t102;
t83 = pkin(4) * t113 + qJDD(5) + t11;
t3 = -pkin(4) * t104 + t83;
t37 = -pkin(4) * t120 + t70 * qJD(2);
t87 = -qJD(1) * t37 - qJDD(1) * t30 - t3;
t86 = t97 + t114;
t85 = 0.2e1 * t104 - t113;
t84 = -t95 + t111;
t82 = -qJDD(1) * t33 + t34 * t79 - t11 - t44;
t81 = -qJDD(4) * t34 + (-qJD(1) * t33 - t125 - t17) * qJD(4);
t63 = t77 * qJ(2);
t62 = t75 * qJ(2);
t39 = -t79 * t74 + t115;
t38 = -t79 * t76 - t116;
t21 = t105 * t71 + t73 * t70 + t78;
t20 = t127 * t76;
t19 = t127 * t74;
t8 = -t96 * t74 + t76 * t98;
t7 = t74 * t98 + t96 * t76;
t5 = t86 * t71 + (t129 * t74 - t115) * t70;
t4 = t85 * t71 + (-t129 * t76 - t116) * t70;
t2 = -t14 + t88 * t74 + (-t84 - t108) * t76;
t1 = qJD(4) * t89 - t74 * t12 + t123 + t145 + t60;
t24 = [qJDD(1), t132, t92, -qJDD(2) + 0.2e1 * t124 + t132, t106 - t92 + 0.2e1 * t112, -t94 * pkin(1) - g(1) * (-t75 * pkin(1) + t63) - g(2) * (t77 * pkin(1) + t62) + (t106 + t112) * qJ(2), t16 * t133 + t15 * t103 - g(1) * (-t78 * t75 + t63) - g(2) * (t78 * t77 + t62) - t90 * qJD(2), t67 * qJDD(1) + t74 * t97, -0.2e1 * t131 * t109 + 0.2e1 * t74 * t113, t38, -t39, 0, t81 * t74 + (-t82 + t93) * t76, t74 * t82 + t76 * t81 + t134, t19 * qJDD(4) + (t74 * t99 + t8) * qJD(4) + (-t87 + t93) * t76, t20 * qJDD(4) + t87 * t74 + (t76 * t99 - t7) * qJD(4) + t134, (qJD(4) * t6 + qJDD(1) * t20 - t2 + (qJD(4) * t19 - t7) * qJD(1)) * t76 + (t10 * qJD(4) + qJDD(1) * t19 + t1 + (-qJD(4) * t20 + t8) * qJD(1)) * t74 + t146, -t2 * t20 + t10 * t7 + t1 * t19 + t6 * t8 + t3 * t30 + t13 * t37 - g(1) * (t144 * t77 - t21 * t75) - g(2) * (t144 * t75 + t21 * t77); 0, 0, 0, -qJDD(1), -t80, -t80 * qJ(2) - t132 + t94, qJD(1) * t90 + t15 * t71 + t16 * t70 - t132, 0, 0, 0, 0, 0, t4, t5, t4, t5, (-qJDD(1) * t70 + t71 * t80) * t130, (-qJD(1) * t91 - t3) * t71 + (-t122 - t1 * t74 + t2 * t76 + (-t10 * t74 - t6 * t76) * qJD(4)) * t70 - t132; 0, 0, 0, 0, 0, 0, t119, 0, 0, 0, 0, 0, t39, t38, t39, t38, 0, qJD(4) * t91 + t1 * t76 + t2 * t74 + g(3); 0, 0, 0, 0, 0, 0, 0, -t74 * t135, t131 * t80, -t114, -t113, qJDD(4), (-t12 + t121) * t74 + t100, t61 * qJD(4) + (-t18 * qJD(4) - t119) * t74 + (-t95 + t121) * t76 + t107, 0.2e1 * t123 + (t10 - t136) * qJD(4) + (pkin(4) * t135 + t122 - t95) * t74 + t100 + t145, -t67 * t80 * pkin(4) + t9 * qJD(4) + (-g(3) - t88) * t74 + ((qJD(5) + t13) * qJD(1) + t84) * t76 + t107, pkin(4) * t114 + (t128 - t143) * t126, t143 * t10 + (t64 + t1 + (t122 + t146) * t74) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t86, -t130 * t80, (t139 + (-t6 - t128) * t74) * qJD(1) + t83 + t93;];
tau_reg = t24;
