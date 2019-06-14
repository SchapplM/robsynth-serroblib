% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPPRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:52:47
% EndTime: 2019-05-04 21:52:50
% DurationCPUTime: 1.39s
% Computational Cost: add. (4205->220), mult. (7646->318), div. (0->0), fcn. (5452->12), ass. (0->149)
t118 = sin(qJ(6));
t121 = cos(qJ(6));
t122 = cos(qJ(5));
t146 = qJD(2) * t122;
t78 = -t121 * qJD(5) + t118 * t146;
t80 = t118 * qJD(5) + t121 * t146;
t66 = t80 * t78;
t145 = qJD(2) * qJD(5);
t140 = t122 * t145;
t119 = sin(qJ(5));
t144 = t119 * qJDD(2);
t83 = -t140 - t144;
t77 = qJDD(6) - t83;
t165 = -t66 + t77;
t167 = t118 * t165;
t166 = t121 * t165;
t113 = sin(pkin(6));
t120 = sin(qJ(2));
t123 = cos(qJ(2));
t111 = sin(pkin(11));
t114 = cos(pkin(11));
t125 = qJD(2) ^ 2;
t86 = t114 * qJDD(2) - t111 * t125;
t87 = t111 * qJDD(2) + t114 * t125;
t164 = (t120 * t86 + t123 * t87) * t113;
t163 = (t120 * t87 - t123 * t86) * t113;
t101 = t122 * qJDD(2);
t141 = t119 * t145;
t84 = t101 - t141;
t139 = -t121 * qJDD(5) + t118 * t84;
t97 = t119 * qJD(2) + qJD(6);
t41 = (qJD(6) - t97) * t80 + t139;
t75 = t78 ^ 2;
t76 = t80 ^ 2;
t96 = t97 ^ 2;
t162 = pkin(3) + pkin(8);
t108 = -g(3) + qJDD(1);
t116 = cos(pkin(6));
t112 = sin(pkin(10));
t115 = cos(pkin(10));
t90 = t112 * g(1) - t115 * g(2);
t161 = t116 * t90;
t131 = t108 * t113 + t161;
t91 = -t115 * g(1) - t112 * g(2);
t57 = -t120 * t91 + t131 * t123;
t126 = qJDD(2) * pkin(2) + t57;
t58 = t131 * t120 + t123 * t91;
t52 = -t125 * pkin(2) + t58;
t36 = t111 * t126 + t114 * t52;
t124 = qJD(5) ^ 2;
t110 = qJDD(2) * pkin(3);
t35 = -t111 * t52 + t114 * t126;
t34 = -t125 * qJ(4) + qJDD(4) - t110 - t35;
t127 = -qJDD(2) * pkin(8) + t34;
t98 = t116 * t108;
t71 = -t113 * t90 + qJDD(3) + t98;
t27 = t119 * t71 - t122 * t127;
t136 = t119 * pkin(5) - t122 * pkin(9);
t81 = t136 * qJD(2);
t19 = -qJDD(5) * pkin(5) - t124 * pkin(9) + t81 * t146 + t27;
t160 = t118 * t19;
t55 = t66 + t77;
t159 = t118 * t55;
t158 = t118 * t97;
t142 = t122 * t125 * t119;
t92 = qJDD(5) + t142;
t157 = t119 * t92;
t156 = t121 * t19;
t155 = t121 * t55;
t154 = t121 * t97;
t153 = t122 * t71;
t93 = qJDD(5) - t142;
t152 = t122 * t93;
t106 = t119 ^ 2;
t151 = t106 * t125;
t107 = t122 ^ 2;
t150 = t107 * t125;
t148 = qJD(6) + t97;
t147 = t106 + t107;
t143 = t119 * t66;
t20 = -t124 * pkin(5) + qJDD(5) * pkin(9) + t153 + (-qJD(2) * t81 + t127) * t119;
t134 = -t83 + t140;
t135 = -t84 + t141;
t103 = 0.2e1 * qJD(4) * qJD(2);
t104 = qJDD(2) * qJ(4);
t138 = t103 + t104 + t36;
t32 = -t162 * t125 + t138;
t22 = t134 * pkin(5) + t135 * pkin(9) + t32;
t8 = t118 * t20 - t121 * t22;
t9 = t118 * t22 + t121 * t20;
t5 = t118 * t8 + t121 * t9;
t137 = -pkin(2) * t87 - t36;
t4 = t118 * t9 - t121 * t8;
t28 = t119 * t127 + t153;
t11 = t119 * t28 - t122 * t27;
t130 = -t118 * qJDD(5) - t121 * t84;
t129 = pkin(2) * t86 + t35;
t128 = qJ(4) + t136;
t60 = -t78 * qJD(6) - t130;
t95 = -t124 - t150;
t94 = -t124 - t151;
t89 = t147 * t125;
t88 = t147 * qJDD(2);
t85 = t101 - 0.2e1 * t141;
t82 = 0.2e1 * t140 + t144;
t74 = t97 * t78;
t73 = -t76 + t96;
t72 = t75 - t96;
t70 = -t119 * t95 - t122 * t92;
t69 = -t119 * t93 + t122 * t94;
t68 = t122 * t95 - t157;
t67 = t119 * t94 + t152;
t65 = t116 * t71;
t64 = t76 - t75;
t63 = -t76 - t96;
t62 = -t111 * t89 + t114 * t88;
t61 = -t96 - t75;
t59 = -t80 * qJD(6) - t139;
t53 = t75 + t76;
t51 = t111 * t85 - t114 * t68;
t50 = t111 * t82 - t114 * t67;
t46 = t148 * t78 + t130;
t45 = t60 + t74;
t44 = t60 - t74;
t42 = -t148 * t80 - t139;
t40 = -t118 * t63 - t155;
t39 = t121 * t63 - t159;
t38 = t121 * t61 - t167;
t37 = t118 * t61 + t166;
t33 = -t125 * pkin(3) + t138;
t30 = t118 * t45 - t121 * t41;
t29 = -t118 * t41 - t121 * t45;
t26 = -t119 * t46 + t122 * t40;
t25 = t119 * t40 + t122 * t46;
t24 = -t119 * t42 + t122 * t38;
t23 = t119 * t38 + t122 * t42;
t18 = -t119 * t53 + t122 * t30;
t17 = t119 * t30 + t122 * t53;
t16 = t111 * t36 + t114 * t35;
t15 = t111 * t33 - t114 * t34;
t14 = t111 * t39 - t114 * t25;
t13 = t111 * t37 - t114 * t23;
t12 = t119 * t27 + t122 * t28;
t10 = t111 * t29 - t114 * t17;
t6 = -t114 * t11 + t111 * t32;
t3 = t119 * t19 + t122 * t5;
t2 = t119 * t5 - t122 * t19;
t1 = t111 * t4 - t114 * t2;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t108, 0, 0, 0, 0, 0, 0, (qJDD(2) * t123 - t120 * t125) * t113, (-qJDD(2) * t120 - t123 * t125) * t113, 0, t116 * t98 + (t120 * t58 + t123 * t57 - t161) * t113, 0, 0, 0, 0, 0, 0, -t163, -t164, 0, t65 + (t120 * (-t111 * t35 + t114 * t36) + t123 * t16) * t113, 0, 0, 0, 0, 0, 0, 0, t163, t164, t65 + (t120 * (t111 * t34 + t114 * t33) + t123 * t15) * t113, 0, 0, 0, 0, 0, 0, t116 * t69 + (t120 * (t111 * t67 + t114 * t82) + t123 * t50) * t113, t116 * t70 + (t120 * (t111 * t68 + t114 * t85) + t123 * t51) * t113, (t120 * (-t111 * t88 - t114 * t89) + t123 * t62) * t113, t116 * t12 + (t120 * (t111 * t11 + t114 * t32) + t123 * t6) * t113, 0, 0, 0, 0, 0, 0, t116 * t24 + (t120 * (t111 * t23 + t114 * t37) + t123 * t13) * t113, t116 * t26 + (t120 * (t111 * t25 + t114 * t39) + t123 * t14) * t113, t116 * t18 + (t120 * (t111 * t17 + t114 * t29) + t123 * t10) * t113, t116 * t3 + (t120 * (t111 * t2 + t114 * t4) + t123 * t1) * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t57, -t58, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t129, t137, 0, pkin(2) * t16, qJDD(2), 0, 0, 0, 0, 0, 0, qJDD(4) - 0.2e1 * t110 - t129, t103 + 0.2e1 * t104 - t137, pkin(2) * t15 - pkin(3) * t34 + qJ(4) * t33, -t135 * t122, -t119 * t85 - t122 * t82, t152 - t119 * (t124 - t150), t134 * t119, t122 * (-t124 + t151) - t157, 0, pkin(2) * t50 + qJ(4) * t82 + t119 * t32 - t162 * t67, pkin(2) * t51 + qJ(4) * t85 + t122 * t32 - t162 * t68, pkin(2) * t62 - qJ(4) * t89 + t162 * t88 - t11, pkin(2) * t6 + qJ(4) * t32 - t162 * t11, t122 * (t121 * t60 - t80 * t158) + t143, t122 * (-t118 * t44 + t121 * t42) + t119 * t64, t122 * (-t118 * t73 + t166) + t119 * t45, t122 * (-t118 * t59 + t78 * t154) - t143, t122 * (t121 * t72 - t159) - t119 * t41, t119 * t77 + t122 * (t118 * t80 - t121 * t78) * t97, pkin(2) * t13 + qJ(4) * t37 + t122 * (-pkin(9) * t37 + t160) - t119 * (-pkin(5) * t37 + t8) - t162 * t23, pkin(2) * t14 + qJ(4) * t39 + t122 * (-pkin(9) * t39 + t156) - t119 * (-pkin(5) * t39 + t9) - t162 * t25, pkin(2) * t10 - t122 * t4 + t128 * t29 - t162 * t17, pkin(2) * t1 + t128 * t4 - t162 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, 0, 0, 0, 0, t69, t70, 0, t12, 0, 0, 0, 0, 0, 0, t24, t26, t18, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t125, t34, 0, 0, 0, 0, 0, 0, t67, t68, -t88, t11, 0, 0, 0, 0, 0, 0, t23, t25, t17, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, (-t106 + t107) * t125, t101, -t142, -t144, qJDD(5), -t27, -t28, 0, 0, t118 * t60 + t80 * t154, t118 * t42 + t121 * t44, t121 * t73 + t167, t121 * t59 + t78 * t158, t118 * t72 + t155, (-t118 * t78 - t121 * t80) * t97, pkin(5) * t42 + pkin(9) * t38 - t156, pkin(5) * t46 + pkin(9) * t40 + t160, pkin(5) * t53 + pkin(9) * t30 + t5, -pkin(5) * t19 + pkin(9) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t64, t45, -t66, -t41, t77, -t8, -t9, 0, 0;];
tauJ_reg  = t7;
