% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPPP1
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
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% tau_reg [4x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:31
% EndTime: 2019-03-08 18:26:31
% DurationCPUTime: 0.44s
% Computational Cost: add. (440->169), mult. (1252->220), div. (0->0), fcn. (939->6), ass. (0->127)
t72 = cos(pkin(4));
t146 = pkin(1) * t72;
t73 = sin(qJ(1));
t145 = g(1) * t73;
t69 = sin(pkin(6));
t65 = t69 ^ 2;
t70 = sin(pkin(4));
t66 = t70 ^ 2;
t144 = t65 * t66;
t75 = qJD(1) ^ 2;
t143 = t66 * t75;
t142 = t69 * t70;
t71 = cos(pkin(6));
t141 = t70 * t71;
t74 = cos(qJ(1));
t140 = t70 * t74;
t139 = t72 * t75;
t138 = t73 * t69;
t137 = t73 * t71;
t136 = t74 * t69;
t135 = t74 * t71;
t122 = qJDD(1) * t70;
t107 = qJ(2) * t122;
t120 = qJD(1) * qJD(2);
t105 = t70 * t120;
t43 = t69 * t105;
t134 = t69 * t107 + t43;
t127 = qJD(1) * t70;
t111 = qJ(2) * t127;
t126 = qJD(1) * t72;
t113 = t69 * t126;
t25 = pkin(1) * t113 + t71 * t111;
t132 = qJ(2) * t70;
t36 = t71 * t132 + t69 * t146;
t133 = t74 * pkin(1) + t73 * t132;
t131 = t72 * qJ(3);
t88 = pkin(3) * t141 + t131;
t8 = qJD(1) * t88 + qJD(4) + t25;
t130 = -qJD(4) - t8;
t129 = pkin(1) * qJDD(1);
t125 = qJD(2) * t70;
t39 = t72 * qJD(3) + t71 * t125;
t128 = qJD(1) * t39;
t46 = t69 * t111;
t124 = qJD(3) + t46;
t123 = qJDD(1) * t69;
t121 = t72 * qJDD(1);
t119 = qJD(1) * qJD(3);
t118 = pkin(3) * t142;
t117 = t71 * t146;
t116 = t71 * t139;
t115 = pkin(1) * t121;
t45 = t71 * t105;
t12 = t71 * t107 + t69 * t115 + t45;
t114 = t69 * t127;
t112 = -pkin(1) * t71 - pkin(2);
t110 = t71 * t122;
t109 = qJDD(3) + t134;
t108 = -t73 * pkin(1) + t74 * t132;
t106 = t66 * t120;
t104 = t69 * t119;
t103 = -qJ(3) * t69 - pkin(1);
t102 = t69 * t71 * t143;
t55 = t69 * t132;
t98 = t112 * t72;
t19 = t55 + t98;
t92 = t112 * qJDD(1);
t6 = t72 * t92 + t109;
t101 = qJDD(1) * t19 + t6;
t89 = -pkin(2) * t71 + t103;
t20 = t89 * t70;
t80 = t89 * qJDD(1);
t7 = qJDD(2) + (t80 - t104) * t70;
t100 = qJDD(1) * t20 + t7;
t99 = -qJ(4) + t112;
t4 = -qJ(3) * t121 - t72 * t119 - t12;
t31 = -t72 * t135 + t138;
t33 = t72 * t137 + t136;
t97 = g(1) * t31 - g(2) * t33;
t32 = t72 * t136 + t137;
t34 = -t72 * t138 + t135;
t96 = g(1) * t32 - g(2) * t34;
t95 = g(1) * t74 + g(2) * t73;
t94 = g(2) * t140 - g(3) * t72 + qJDD(2);
t93 = (qJD(1) * t117 - t46) * t69 - t25 * t71;
t91 = t34 * pkin(2) + t33 * qJ(3) + t133;
t90 = -qJD(3) * t69 - qJD(4) * t71;
t87 = t99 * qJDD(1);
t54 = -pkin(1) * t122 + qJDD(2);
t86 = t66 * t129 - t54 * t70;
t50 = qJDD(1) * t118;
t1 = t50 + (-qJD(1) * qJD(4) + t87) * t72 + t109;
t76 = t72 * t99 + t118;
t10 = t55 + t76;
t38 = -t72 * qJD(4) + t69 * t125;
t85 = qJD(1) * t38 + qJDD(1) * t10 + t1;
t13 = t88 + t36;
t3 = pkin(3) * t110 + qJDD(4) - t4;
t84 = qJDD(1) * t13 + t128 + t3;
t18 = -t131 - t36;
t83 = -qJDD(1) * t18 + t128 - t4;
t82 = -t32 * pkin(2) - t31 * qJ(3) + t108;
t81 = (-pkin(2) - qJ(4)) * t71 + t103;
t14 = t81 * t70;
t77 = t81 * qJDD(1);
t2 = qJDD(2) + (qJD(1) * t90 + t77) * t70;
t30 = t90 * t70;
t79 = (-qJD(1) * t30 - qJDD(1) * t14 - t2) * t70;
t78 = -g(1) * t33 - g(2) * t31 + g(3) * t141 + t109;
t68 = t72 ^ 2;
t67 = t71 ^ 2;
t41 = t65 * t106;
t40 = t139 * t142;
t37 = (-t68 - t144) * t75;
t35 = -t55 + t117;
t29 = (-t65 - t67) * t143;
t28 = t102 + t121;
t23 = -t40 + t110;
t22 = (-t116 + t123) * t70;
t21 = (t116 + t123) * t70;
t17 = -qJ(3) * t126 - t25;
t16 = qJD(1) * t20 + qJD(2);
t15 = qJD(1) * t98 + t124;
t11 = t71 * t115 - t134;
t9 = qJD(1) * t14 + qJD(2);
t5 = qJD(1) * t76 + t124;
t24 = [qJDD(1), -g(2) * t74 + t145, t95, t86 * t71 + (qJDD(1) * t35 + t11 - t43) * t72 + t96, -t86 * t69 + (-qJDD(1) * t36 - t12 - t45) * t72 - t97, t67 * t106 + t41 + (-t11 * t69 + t12 * t71 + (-t35 * t69 + t36 * t71) * qJDD(1) - t95) * t70, t12 * t36 + t11 * t35 - g(1) * t108 - g(2) * t133 + (-t54 * pkin(1) - qJD(2) * t93) * t70, t41 + (t101 * t69 + t83 * t71 - t95) * t70, -t66 * t71 * t104 + t101 * t72 + (qJD(2) * t113 + t100 * t71) * t70 - t96, -t100 * t142 + t119 * t144 + t83 * t72 + t97, t7 * t20 + t4 * t18 - t17 * t39 + t6 * t19 - g(1) * t82 - g(2) * t91 + (qJD(2) * t15 - qJD(3) * t16) * t142 (t69 * t85 + t71 * t84 - t95) * t70, t69 * t79 + t72 * t84 + t97, t71 * t79 - t72 * t85 + t96, t2 * t14 + t9 * t30 + t1 * t10 + t5 * t38 + t3 * t13 + t8 * t39 - g(1) * (pkin(3) * t140 - t32 * qJ(4) + t82) - g(2) * (t73 * t70 * pkin(3) + t34 * qJ(4) + t91); 0, 0, 0, -t23, t21, t29 (qJD(1) * t93 - t129 - t145) * t70 + t94, t29, t23, -t21 (-t145 + t80 + (t17 * t71 + (-qJD(3) - t15) * t69) * qJD(1)) * t70 + t94, t29, -t21, -t23 (-t145 + t77 + (t130 * t71 + (-qJD(3) - t5) * t69) * qJD(1)) * t70 + t94; 0, 0, 0, 0, 0, 0, 0, t22, t28, t37, t16 * t114 + (qJD(1) * t17 + t92) * t72 + t78, t22, t37, -t28, t9 * t114 + t50 + (t130 * qJD(1) + t87) * t72 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 + t110, -t102 + t121 (-t66 * t67 - t68) * t75, -g(3) * t142 - g(1) * t34 - g(2) * t32 + (t9 * t141 + t5 * t72) * qJD(1) + t3;];
tau_reg  = t24;
