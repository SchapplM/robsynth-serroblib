% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:30
% EndTime: 2019-12-05 16:04:36
% DurationCPUTime: 1.42s
% Computational Cost: add. (846->230), mult. (1863->339), div. (0->0), fcn. (1452->10), ass. (0->135)
t63 = sin(qJ(4));
t50 = t63 * qJD(2) + qJD(5);
t161 = t50 - qJD(5);
t111 = qJD(2) * qJD(4);
t66 = cos(qJ(4));
t112 = t66 * qJDD(2);
t160 = -t63 * t111 + t112;
t61 = cos(pkin(5));
t116 = qJDD(1) * t61;
t102 = t66 * t116;
t68 = -pkin(2) - pkin(7);
t59 = sin(pkin(5));
t117 = qJDD(1) * t59;
t67 = cos(qJ(2));
t104 = t67 * t117;
t132 = qJD(1) * t59;
t64 = sin(qJ(2));
t107 = t64 * t132;
t46 = qJD(2) * t107;
t82 = qJDD(3) + t46 - t104;
t21 = t68 * qJDD(2) + t82;
t131 = qJD(1) * t61;
t106 = t67 * t132;
t92 = qJD(3) - t106;
t34 = t68 * qJD(2) + t92;
t85 = t63 * t131 - t66 * t34;
t1 = qJDD(4) * pkin(8) - qJD(4) * t85 + t63 * t21 + t102;
t126 = qJD(4) * t66;
t103 = t63 * t116;
t20 = t66 * t131 + t63 * t34;
t2 = -qJDD(4) * pkin(4) + t20 * qJD(4) - t66 * t21 + t103;
t43 = t63 * pkin(4) - t66 * pkin(8) + qJ(3);
t23 = qJD(2) * t43 + t107;
t100 = t66 * t111;
t113 = t63 * qJDD(2);
t35 = qJDD(5) + t100 + t113;
t9 = -qJD(4) * pkin(4) + t85;
t144 = t61 * t67;
t58 = sin(pkin(9));
t60 = cos(pkin(9));
t28 = -t60 * t144 + t58 * t64;
t30 = t58 * t144 + t60 * t64;
t94 = g(1) * t30 + g(2) * t28;
t159 = -t50 * (qJD(5) * t43 + t68 * t126) + t2 * t66 + (-t9 * qJD(4) - qJD(5) * t23 - t68 * t35 - t1) * t63 + t94;
t149 = t59 * t63;
t146 = t59 * t67;
t32 = t66 * t146 + t61 * t63;
t80 = g(1) * (-t58 * t149 + t30 * t66) + g(2) * (t60 * t149 + t28 * t66) - g(3) * t32;
t95 = pkin(4) * t66 + pkin(8) * t63;
t158 = t50 * (pkin(8) * qJD(5) + t95 * qJD(2)) + t2 + t80;
t62 = sin(qJ(5));
t121 = t62 * qJD(4);
t129 = qJD(2) * t66;
t65 = cos(qJ(5));
t40 = t65 * t129 + t121;
t8 = qJD(5) * t40 - t65 * qJDD(4) + t160 * t62;
t10 = qJD(4) * pkin(8) + t20;
t145 = t61 * t64;
t29 = t60 * t145 + t58 * t67;
t31 = -t58 * t145 + t60 * t67;
t93 = g(1) * t31 + g(2) * t29;
t157 = -(t50 * t68 + t10) * qJD(5) - t93;
t119 = qJDD(1) - g(3);
t148 = t59 * t64;
t156 = -t119 * t148 + t93;
t118 = qJD(2) * qJ(3);
t42 = t107 + t118;
t155 = qJD(4) * (t107 - t42 - t118) - qJDD(4) * t68;
t120 = t65 * qJD(4);
t124 = qJD(5) * t66;
t79 = -t63 * t120 - t62 * t124;
t7 = qJD(2) * t79 + qJD(5) * t120 + t62 * qJDD(4) + t65 * t112;
t154 = t66 * t7;
t153 = t7 * t62;
t38 = t62 * t129 - t120;
t152 = t38 * t50;
t151 = t40 * t50;
t150 = t50 * t65;
t147 = t59 * t66;
t143 = t62 * t35;
t142 = t62 * t64;
t141 = t63 * t68;
t140 = t64 * t65;
t139 = t65 * t35;
t138 = t66 * t38;
t137 = t66 * t40;
t136 = t66 * t68;
t70 = qJD(2) ^ 2;
t135 = t67 * t70;
t57 = t66 ^ 2;
t134 = t63 ^ 2 - t57;
t69 = qJD(4) ^ 2;
t133 = -t69 - t70;
t130 = qJD(2) * t59;
t128 = qJD(4) * t38;
t127 = qJD(4) * t40;
t125 = qJD(5) * t50;
t123 = qJDD(2) * pkin(2);
t122 = t42 * qJD(2);
t115 = qJDD(4) * t63;
t110 = qJDD(2) * qJ(3);
t109 = t64 * t130;
t108 = t67 * t130;
t105 = t64 * t117;
t97 = -t21 + t122;
t96 = qJD(5) * t63 + qJD(2);
t4 = t65 * t10 + t62 * t23;
t91 = t62 * t10 - t65 * t23;
t90 = (-qJD(2) * pkin(2) + t92) * t64 + t42 * t67;
t89 = qJDD(2) * t64 + t135;
t33 = -t63 * t146 + t61 * t66;
t17 = t59 * t140 - t33 * t62;
t18 = t59 * t142 + t33 * t65;
t87 = t63 * t140 + t62 * t67;
t86 = t63 * t142 - t65 * t67;
t84 = t65 * t125 + t143;
t83 = -t62 * t125 + t139;
t36 = qJD(4) * t95 + qJD(3);
t76 = -t9 * t124 - t43 * t35 - t36 * t50;
t75 = -g(3) * t146 + t104 + t94;
t74 = -pkin(8) * t35 + (-t85 + t9) * t50;
t73 = qJDD(3) - t75;
t22 = t105 + t110 + (qJD(3) + t106) * qJD(2);
t71 = -g(3) * t148 + t92 * qJD(2) - t68 * t69 + t110 + t22 - t93;
t54 = qJDD(4) * t66;
t27 = t89 * t59;
t26 = (-qJDD(2) * t67 + t64 * t70) * t59;
t24 = t82 - t123;
t16 = qJD(4) * t33 - t109 * t66;
t15 = -qJD(4) * t32 + t109 * t63;
t14 = t60 * t147 - t28 * t63;
t12 = t58 * t147 + t30 * t63;
t6 = t105 + t43 * qJDD(2) + (t36 + t106) * qJD(2);
t5 = t65 * t6;
t3 = [t119, 0, -t26, -t27, t26, t27, t61 ^ 2 * qJDD(1) - g(3) + (qJD(2) * t90 + t22 * t64 - t24 * t67) * t59, 0, 0, 0, 0, 0, -t16 * qJD(4) - t32 * qJDD(4) + (t100 * t64 + t63 * t89) * t59, -t15 * qJD(4) - t33 * qJDD(4) + (t66 * t135 + t160 * t64) * t59, 0, 0, 0, 0, 0, (-qJD(5) * t18 + t108 * t65 - t15 * t62) * t50 + t17 * t35 + t16 * t38 + t32 * t8, -(qJD(5) * t17 + t108 * t62 + t15 * t65) * t50 - t18 * t35 + t16 * t40 + t32 * t7; 0, qJDD(2), t75, t156, t73 - 0.2e1 * t123, 0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t110 - t156, t22 * qJ(3) + t42 * qJD(3) - t24 * pkin(2) - g(1) * (-t30 * pkin(2) + t31 * qJ(3)) - g(2) * (-t28 * pkin(2) + t29 * qJ(3)) + (-g(3) * (pkin(2) * t67 + qJ(3) * t64) - t90 * qJD(1)) * t59, t57 * qJDD(2) - 0.2e1 * t100 * t63, 0.2e1 * t134 * t111 - 0.2e1 * t63 * t112, -t69 * t63 + t54, -t69 * t66 - t115, 0, -t155 * t66 + t71 * t63, t155 * t63 + t71 * t66, t65 * t154 + t40 * t79, (t38 * t65 + t40 * t62) * t63 * qJD(4) + (-t153 - t65 * t8 + (t38 * t62 - t40 * t65) * qJD(5)) * t66, (-t50 * t120 + t7) * t63 + (t83 + t127) * t66, (t50 * t121 - t8) * t63 + (-t84 - t128) * t66, t50 * t126 + t35 * t63, -t8 * t136 + t5 * t63 + (t38 * t141 - t66 * t91) * qJD(4) + (t157 * t63 - t76) * t65 + t159 * t62 + (-g(3) * t87 + (t64 * t138 + t50 * t86) * qJD(1)) * t59, -t7 * t136 + (t40 * t141 - t4 * t66) * qJD(4) + ((-t157 - t6) * t63 + t76) * t62 + t159 * t65 + (g(3) * t86 + (t64 * t137 + t50 * t87) * qJD(1)) * t59; 0, 0, 0, 0, qJDD(2), -t70, -t122 + t46 + t73 - t123, 0, 0, 0, 0, 0, t133 * t63 + t54, t133 * t66 - t115, 0, 0, 0, 0, 0, -t66 * t8 + (t128 - t143) * t63 + (-t66 * t121 - t65 * t96) * t50, -t154 + (t127 - t139) * t63 + (-t66 * t120 + t62 * t96) * t50; 0, 0, 0, 0, 0, 0, 0, t66 * t70 * t63, -t134 * t70, t112, -t113, qJDD(4), -t66 * t97 - t103 - t80, g(1) * t12 - g(2) * t14 + g(3) * t33 + t97 * t63 - t102, t40 * t150 + t153, (t7 - t152) * t65 + (-t8 - t151) * t62, (t63 * t150 - t137) * qJD(2) + t84, (-t62 * t63 * t50 + t138) * qJD(2) + t83, -t50 * t129, -pkin(4) * t8 + t129 * t91 - t158 * t65 - t20 * t38 + t74 * t62, -pkin(4) * t7 + t4 * t129 + t158 * t62 - t20 * t40 + t74 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t38, -t38 ^ 2 + t40 ^ 2, t7 + t152, t151 - t8, t35, -t62 * t1 + t5 - t9 * t40 - g(1) * (-t12 * t62 + t31 * t65) - g(2) * (t14 * t62 + t29 * t65) - g(3) * t17 + t161 * t4, -t65 * t1 - t62 * t6 + t9 * t38 - g(1) * (-t12 * t65 - t31 * t62) - g(2) * (t14 * t65 - t29 * t62) + g(3) * t18 - t161 * t91;];
tau_reg = t3;
