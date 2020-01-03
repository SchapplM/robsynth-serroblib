% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:50
% EndTime: 2019-12-31 19:27:52
% DurationCPUTime: 0.61s
% Computational Cost: add. (901->160), mult. (1117->190), div. (0->0), fcn. (618->10), ass. (0->106)
t77 = qJ(1) + qJ(2);
t68 = sin(t77);
t69 = cos(t77);
t78 = sin(pkin(8));
t79 = cos(pkin(8));
t30 = -t68 * t79 + t69 * t78;
t31 = t68 * t78 + t69 * t79;
t140 = g(1) * t30 + g(2) * t31;
t124 = pkin(1) * qJD(1);
t109 = qJD(2) * t124;
t123 = pkin(1) * qJDD(1);
t81 = sin(qJ(2));
t84 = cos(qJ(2));
t126 = -t81 * t109 + t84 * t123;
t108 = -qJDD(3) + t126;
t72 = qJDD(1) + qJDD(2);
t86 = -pkin(2) - pkin(3);
t14 = t86 * t72 - t108;
t129 = t84 * t109 + t81 * t123;
t64 = t72 * qJ(3);
t73 = qJD(1) + qJD(2);
t65 = t73 * qJD(3);
t15 = t64 + t65 + t129;
t5 = t78 * t14 + t79 * t15;
t94 = -g(1) * t31 + g(2) * t30 + t5;
t139 = -g(1) * t69 - g(2) * t68;
t111 = t84 * t124;
t112 = t81 * t124;
t28 = t78 * t111 - t79 * t112;
t138 = (t78 * qJD(3) - t28) * t73;
t71 = t73 ^ 2;
t135 = t72 * pkin(2);
t82 = sin(qJ(1));
t134 = t82 * pkin(1);
t122 = qJD(2) * t81;
t113 = pkin(1) * t122;
t121 = qJD(2) * t84;
t47 = pkin(1) * t121 + qJD(3);
t23 = -t79 * t113 + t78 * t47;
t133 = t23 * t73;
t80 = sin(qJ(5));
t132 = t80 * t72;
t83 = cos(qJ(5));
t131 = t83 * t72;
t130 = t140 * t80;
t61 = -t84 * pkin(1) - pkin(2);
t50 = -pkin(3) + t61;
t53 = t81 * pkin(1) + qJ(3);
t22 = t78 * t50 + t79 * t53;
t128 = t69 * pkin(2) + t68 * qJ(3);
t127 = g(1) * t68 - g(2) * t69;
t38 = t79 * qJ(3) + t78 * t86;
t75 = t80 ^ 2;
t125 = -t83 ^ 2 + t75;
t120 = qJD(5) * t73;
t119 = qJDD(4) + g(3);
t118 = qJDD(5) * t80;
t117 = qJDD(5) * t83;
t87 = qJD(5) ^ 2;
t116 = t87 + t71;
t115 = 0.2e1 * t120;
t85 = cos(qJ(1));
t114 = t85 * pkin(1) + t128;
t110 = t73 * t122;
t52 = t69 * qJ(3);
t107 = -t68 * pkin(2) + t52;
t4 = t79 * t14 - t78 * t15;
t106 = t83 * t115;
t105 = t126 + t127;
t104 = t129 + t139;
t29 = (t78 * t81 + t79 * t84) * t124;
t102 = t79 * qJD(3) - t29;
t101 = t86 * t68 + t52;
t97 = qJD(3) - t111;
t27 = t86 * t73 + t97;
t36 = t73 * qJ(3) + t112;
t8 = t79 * t27 - t78 * t36;
t9 = t78 * t27 + t79 * t36;
t98 = t8 * t78 - t9 * t79;
t37 = -t78 * qJ(3) + t79 * t86;
t34 = pkin(4) - t37;
t35 = -pkin(7) + t38;
t96 = -t34 * t72 + t35 * t87;
t21 = t79 * t50 - t78 * t53;
t95 = -qJDD(3) + t105;
t25 = -t108 - t135;
t16 = pkin(4) - t21;
t17 = -pkin(7) + t22;
t93 = -t16 * t72 + t17 * t87 - t133;
t6 = t73 * pkin(4) - t8;
t92 = t72 * pkin(7) + t6 * t73 - t94;
t91 = t73 * t111 - t104;
t24 = t78 * t113 + t79 * t47;
t90 = -qJDD(5) * t17 + (-t16 * t73 - t24 - t6) * qJD(5);
t89 = -t140 + t138;
t88 = -qJDD(5) * t35 + (-t34 * t73 - t102 - t6) * qJD(5);
t55 = t69 * pkin(3);
t43 = t73 * t112;
t40 = -t87 * t80 + t117;
t39 = -t87 * t83 - t118;
t33 = -t73 * pkin(2) + t97;
t26 = t80 * t106 + t75 * t72;
t11 = -0.2e1 * t125 * t120 + 0.2e1 * t80 * t131;
t2 = t72 * pkin(4) - t4;
t1 = t2 * t83;
t3 = [qJDD(1), g(1) * t82 - g(2) * t85, g(1) * t85 + g(2) * t82, t72, (t72 * t84 - t110) * pkin(1) + t105, (-t73 * t121 - t72 * t81) * pkin(1) - t104, -pkin(1) * t110 + (pkin(2) - t61) * t72 + t95, t47 * t73 + t53 * t72 + t139 + t15, t15 * t53 + t36 * t47 + t25 * t61 + t33 * t113 - g(1) * (t107 - t134) - g(2) * t114, -t21 * t72 + t133 - t140 - t4, t22 * t72 + t24 * t73 + t94, t5 * t22 + t9 * t24 + t4 * t21 - t8 * t23 - g(1) * (t101 - t134) - g(2) * (t55 + t114), t26, t11, t39, -t40, 0, t1 + t90 * t80 + (-t93 - t140) * t83, t90 * t83 + (-t2 + t93) * t80 + t130; 0, 0, 0, t72, t43 + t105, t91, t43 + t95 + 0.2e1 * t135, 0.2e1 * t64 + 0.2e1 * t65 - t91, t15 * qJ(3) + t36 * qJD(3) - t25 * pkin(2) - g(1) * t107 - g(2) * t128 + (-t33 * t81 - t36 * t84) * t124, -t37 * t72 - t4 + t89, t102 * t73 + t38 * t72 + t94, t5 * t38 + t4 * t37 - t9 * t29 + t8 * t28 - g(1) * t101 - g(2) * (t55 + t128) - t98 * qJD(3), t26, t11, t39, -t40, 0, t1 + t88 * t80 + (t89 - t96) * t83, t88 * t83 + (-t2 + t96 - t138) * t80 + t130; 0, 0, 0, 0, 0, 0, -t72, -t71, -t36 * t73 - t127 + t25, -t78 * t71 - t79 * t72, -t79 * t71 + t78 * t72, t4 * t79 + t5 * t78 + t98 * t73 - t127, 0, 0, 0, 0, 0, (t80 * t115 - t131) * t79 + (-t116 * t83 - t118) * t78, (t106 + t132) * t79 + (t116 * t80 - t117) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, 0, 0, 0, 0, 0, t40, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 * t71 * t83, t125 * t71, -t132, -t131, qJDD(5), t119 * t83 + t92 * t80, -t119 * t80 + t92 * t83;];
tau_reg = t3;
