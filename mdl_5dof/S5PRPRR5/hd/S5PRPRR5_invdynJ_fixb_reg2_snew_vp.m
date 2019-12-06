% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRR5
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:43
% EndTime: 2019-12-05 15:54:50
% DurationCPUTime: 2.16s
% Computational Cost: add. (5968->254), mult. (13922->381), div. (0->0), fcn. (10402->10), ass. (0->165)
t151 = sin(qJ(4));
t146 = sin(pkin(9));
t148 = cos(pkin(9));
t154 = cos(qJ(4));
t121 = (-t146 * t151 + t148 * t154) * qJD(2);
t164 = t146 * t154 + t148 * t151;
t123 = t164 * qJD(2);
t180 = t123 * t121;
t200 = qJDD(4) + t180;
t205 = t151 * t200;
t204 = t154 * t200;
t147 = sin(pkin(8));
t149 = cos(pkin(8));
t131 = -t149 * g(1) - t147 * g(2);
t152 = sin(qJ(2));
t155 = cos(qJ(2));
t177 = -g(3) + qJDD(1);
t110 = t155 * t131 + t152 * t177;
t193 = qJD(2) ^ 2;
t201 = -t193 * pkin(2) + qJDD(2) * qJ(3) + 0.2e1 * qJD(2) * qJD(3) + t110;
t150 = sin(qJ(5));
t142 = qJDD(4) + qJDD(5);
t153 = cos(qJ(5));
t91 = -t153 * t121 + t150 * t123;
t93 = t150 * t121 + t153 * t123;
t70 = t93 * t91;
t198 = -t70 + t142;
t203 = t150 * t198;
t202 = t153 * t198;
t174 = t123 * qJD(4);
t170 = t148 * qJDD(2);
t171 = t146 * qJDD(2);
t80 = -t151 * t171 + t154 * t170;
t101 = t80 - t174;
t120 = t164 * qJDD(2);
t175 = t121 * qJD(4);
t103 = t120 + t175;
t61 = -t91 * qJD(5) + t150 * t101 + t153 * t103;
t143 = qJD(4) + qJD(5);
t88 = t143 * t91;
t199 = t61 - t88;
t197 = t148 * t193;
t157 = t146 ^ 2;
t159 = t148 ^ 2;
t196 = t157 + t159;
t130 = -t147 * g(1) + t149 * g(2);
t127 = t148 * t130;
t195 = t127 + (pkin(3) * t197 - pkin(6) * qJDD(2) - t201) * t146;
t129 = t196 * t193;
t167 = t146 * t130 + t201 * t148;
t173 = t159 * t193;
t76 = -pkin(3) * t173 + pkin(6) * t170 + t167;
t52 = t151 * t76 - t154 * t195;
t194 = -t52 + (-t103 + t175) * pkin(7);
t89 = t91 ^ 2;
t90 = t93 ^ 2;
t118 = t121 ^ 2;
t119 = t123 ^ 2;
t141 = t143 ^ 2;
t161 = pkin(4) * t200 + t194;
t108 = qJD(4) * pkin(4) - t123 * pkin(7);
t53 = t195 * t151 + t154 * t76;
t33 = -t118 * pkin(4) + t101 * pkin(7) - qJD(4) * t108 + t53;
t17 = t150 * t33 - t153 * t161;
t185 = t153 * t33;
t18 = t150 * t161 + t185;
t8 = t150 * t18 - t153 * t17;
t192 = t151 * t8;
t191 = t154 * t8;
t30 = t151 * t53 - t154 * t52;
t190 = t146 * t30;
t109 = -t152 * t131 + t155 * t177;
t145 = qJDD(2) * pkin(2);
t106 = -t193 * qJ(3) + qJDD(3) - t109 - t145;
t84 = -pkin(3) * t170 + t106 + (-t157 * t193 - t173) * pkin(6);
t51 = -t101 * pkin(4) - t118 * pkin(7) + t123 * t108 + t84;
t189 = t150 * t51;
t67 = t70 + t142;
t188 = t150 * t67;
t187 = t151 * t84;
t98 = qJDD(4) - t180;
t186 = t151 * t98;
t184 = t153 * t51;
t183 = t153 * t67;
t182 = t154 * t84;
t181 = t154 * t98;
t179 = t143 * t150;
t178 = t143 * t153;
t169 = t155 * qJDD(2);
t62 = t146 * (t201 * t146 - t127) + t148 * t167;
t9 = t150 * t17 + t153 * t18;
t31 = t151 * t52 + t154 * t53;
t166 = -t153 * t101 + t150 * t103;
t165 = -t106 + t145;
t163 = (-qJD(5) + t143) * t93 - t166;
t156 = qJD(4) ^ 2;
t140 = t159 * qJDD(2);
t139 = t157 * qJDD(2);
t128 = t140 + t139;
t125 = t196 * t197;
t124 = t146 * t129;
t113 = -t119 - t156;
t112 = -t119 + t156;
t111 = t118 - t156;
t102 = t120 + 0.2e1 * t175;
t100 = -t80 + 0.2e1 * t174;
t96 = -t156 - t118;
t86 = -t90 + t141;
t85 = t89 - t141;
t83 = -t90 - t141;
t82 = -t118 - t119;
t78 = -t151 * t113 - t181;
t77 = t154 * t113 - t186;
t74 = t151 * t120 + t154 * t80;
t73 = -t154 * t120 + t151 * t80;
t72 = t154 * t96 - t205;
t71 = t151 * t96 + t204;
t69 = t90 - t89;
t65 = -t141 - t89;
t64 = (t150 * t93 - t153 * t91) * t143;
t63 = (-t150 * t91 - t153 * t93) * t143;
t60 = -t93 * qJD(5) - t166;
t59 = -t89 - t90;
t58 = -t146 * t77 + t148 * t78;
t57 = t153 * t85 - t188;
t56 = -t150 * t86 + t202;
t55 = t150 * t85 + t183;
t54 = t153 * t86 + t203;
t50 = -t150 * t83 - t183;
t49 = t153 * t83 - t188;
t47 = -t146 * t73 + t148 * t74;
t46 = t61 + t88;
t41 = (qJD(5) + t143) * t93 + t166;
t40 = -t146 * t71 + t148 * t72;
t39 = t153 * t61 - t93 * t179;
t38 = t150 * t61 + t93 * t178;
t37 = -t150 * t60 + t91 * t178;
t36 = t153 * t60 + t91 * t179;
t35 = t153 * t65 - t203;
t34 = t150 * t65 + t202;
t29 = -pkin(7) * t49 + t184;
t28 = -t151 * t49 + t154 * t50;
t27 = t151 * t50 + t154 * t49;
t26 = -pkin(7) * t34 + t189;
t25 = t150 * t46 + t153 * t163;
t24 = -t150 * t199 - t153 * t41;
t23 = t150 * t163 - t153 * t46;
t22 = -t150 * t41 + t153 * t199;
t21 = -t151 * t34 + t154 * t35;
t20 = t151 * t35 + t154 * t34;
t19 = -pkin(4) * t199 + pkin(7) * t50 + t189;
t15 = -pkin(4) * t41 + pkin(7) * t35 - t184;
t14 = t148 * t31 - t190;
t13 = -t146 * t27 + t148 * t28;
t12 = -t151 * t23 + t154 * t25;
t11 = t151 * t25 + t154 * t23;
t10 = -t146 * t20 + t148 * t21;
t7 = -pkin(4) * t51 + pkin(7) * t9;
t6 = -pkin(7) * t23 - t8;
t5 = -t146 * t11 + t148 * t12;
t4 = -pkin(4) * t59 + pkin(7) * t25 + t9;
t3 = t154 * t9 - t192;
t2 = t151 * t9 + t191;
t1 = -t146 * t2 + t148 * t3;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t177, 0, 0, 0, 0, 0, 0, -t152 * t193 + t169, -t152 * qJDD(2) - t155 * t193, 0, t155 * t109 + t152 * t110, 0, 0, 0, 0, 0, 0, -t152 * t125 + t148 * t169, t152 * t124 - t146 * t169, t152 * t128 + t155 * t129, -t155 * t106 + t152 * t62, 0, 0, 0, 0, 0, 0, -t155 * t100 + t152 * t40, -t155 * t102 + t152 * t58, t152 * t47 - t155 * t82, t152 * t14 - t155 * t84, 0, 0, 0, 0, 0, 0, t152 * t10 - t155 * t41, t152 * t13 - t155 * t199, t152 * t5 - t155 * t59, t152 * t1 - t155 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t109, -t110, 0, 0, t139, 0.2e1 * t146 * t170, 0, t140, 0, 0, -qJ(3) * t125 + t165 * t148, qJ(3) * t124 - t146 * t165, pkin(2) * t129 + qJ(3) * t128 + t62, -pkin(2) * t106 + qJ(3) * t62, t146 * (t154 * t103 - t151 * t174) + t148 * (t151 * t103 + t154 * t174), t146 * (-t154 * t100 - t151 * t102) + t148 * (-t151 * t100 + t154 * t102), t146 * (-t151 * t112 + t204) + t148 * (t154 * t112 + t205), t146 * (-t151 * t101 - t154 * t175) + t148 * (t154 * t101 - t151 * t175), t146 * (t154 * t111 - t186) + t148 * (t151 * t111 + t181), (t146 * (t121 * t154 + t123 * t151) + t148 * (t121 * t151 - t123 * t154)) * qJD(4), t146 * (-pkin(6) * t71 + t187) + t148 * (-pkin(3) * t100 + pkin(6) * t72 - t182) - pkin(2) * t100 + qJ(3) * t40, t146 * (-pkin(6) * t77 + t182) + t148 * (-pkin(3) * t102 + pkin(6) * t78 + t187) - pkin(2) * t102 + qJ(3) * t58, t146 * (-pkin(6) * t73 - t30) + t148 * (-pkin(3) * t82 + pkin(6) * t74 + t31) - pkin(2) * t82 + qJ(3) * t47, -pkin(6) * t190 + t148 * (-pkin(3) * t84 + pkin(6) * t31) - pkin(2) * t84 + qJ(3) * t14, t146 * (-t151 * t38 + t154 * t39) + t148 * (t151 * t39 + t154 * t38), t146 * (-t151 * t22 + t154 * t24) + t148 * (t151 * t24 + t154 * t22), t146 * (-t151 * t54 + t154 * t56) + t148 * (t151 * t56 + t154 * t54), t146 * (-t151 * t36 + t154 * t37) + t148 * (t151 * t37 + t154 * t36), t146 * (-t151 * t55 + t154 * t57) + t148 * (t151 * t57 + t154 * t55), t146 * (-t151 * t63 + t154 * t64) + t148 * (t151 * t64 + t154 * t63), t146 * (-pkin(6) * t20 - t151 * t15 + t154 * t26) + t148 * (-pkin(3) * t41 + pkin(6) * t21 + t154 * t15 + t151 * t26) - pkin(2) * t41 + qJ(3) * t10, t146 * (-pkin(6) * t27 - t151 * t19 + t154 * t29) + t148 * (-pkin(3) * t199 + pkin(6) * t28 + t151 * t29 + t154 * t19) - pkin(2) * t199 + qJ(3) * t13, t146 * (-pkin(6) * t11 - t151 * t4 + t154 * t6) + t148 * (-pkin(3) * t59 + pkin(6) * t12 + t151 * t6 + t154 * t4) - pkin(2) * t59 + qJ(3) * t5, t146 * (-pkin(6) * t2 - pkin(7) * t191 - t151 * t7) + t148 * (-pkin(3) * t51 + pkin(6) * t3 - pkin(7) * t192 + t154 * t7) - pkin(2) * t51 + qJ(3) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, t171, -t129, t106, 0, 0, 0, 0, 0, 0, t100, t102, t82, t84, 0, 0, 0, 0, 0, 0, t41, t199, t59, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t119 - t118, t120, t180, t80, qJDD(4), -t52, -t53, 0, 0, t70, t69, t46, -t70, t163, t142, pkin(4) * t34 - t17, -t185 - t150 * t194 + (-t150 * t200 + t49) * pkin(4), pkin(4) * t23, pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t69, t46, -t70, t163, t142, -t17, -t18, 0, 0;];
tauJ_reg = t16;
