% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 23:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRPR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:21:41
% EndTime: 2019-05-04 23:21:47
% DurationCPUTime: 1.75s
% Computational Cost: add. (3631->244), mult. (7096->300), div. (0->0), fcn. (4358->10), ass. (0->159)
t120 = sin(pkin(6));
t122 = cos(pkin(6));
t127 = sin(qJ(4));
t128 = sin(qJ(2));
t130 = cos(qJ(4));
t131 = cos(qJ(2));
t133 = qJD(2) ^ 2;
t167 = t130 * t133;
t99 = t127 * t167;
t87 = qJDD(4) + t99;
t179 = t127 * t87;
t189 = t130 ^ 2;
t108 = t189 * t133;
t132 = qJD(4) ^ 2;
t95 = -t108 - t132;
t54 = -t130 * t95 + t179;
t163 = qJD(2) * qJD(4);
t105 = t127 * t163;
t159 = t130 * qJDD(2);
t81 = -0.2e1 * t105 + t159;
t206 = t120 * (t128 * t81 + t131 * t54) - t122 * (t127 * t95 + t130 * t87);
t88 = qJDD(4) - t99;
t172 = t130 * t88;
t115 = t127 ^ 2;
t107 = t115 * t133;
t93 = -t107 - t132;
t52 = t127 * t93 + t172;
t156 = t130 * t163;
t160 = t127 * qJDD(2);
t78 = 0.2e1 * t156 + t160;
t205 = t120 * (-t128 * t78 + t131 * t52) + t122 * (t127 * t88 - t130 * t93);
t204 = pkin(8) + pkin(2);
t203 = t204 * t52;
t202 = t204 * t54;
t80 = -t105 + t159;
t201 = t80 - t105;
t126 = sin(qJ(6));
t129 = cos(qJ(6));
t166 = qJD(2) * t127;
t72 = qJD(4) * t126 - t129 * t166;
t74 = qJD(4) * t129 + t126 * t166;
t49 = t74 * t72;
t69 = qJDD(6) + t80;
t195 = -t49 + t69;
t200 = t126 * t195;
t199 = t129 * t195;
t116 = -g(3) + qJDD(1);
t119 = sin(pkin(10));
t121 = cos(pkin(10));
t85 = g(1) * t119 - g(2) * t121;
t198 = t116 * t120 + t122 * t85;
t168 = t130 * qJ(5);
t152 = t127 * pkin(4) - t168;
t75 = t152 * qJD(2);
t194 = -t132 * pkin(4) - t75 * t166;
t162 = qJD(3) * qJD(2);
t110 = 0.2e1 * t162;
t165 = qJD(2) * t130;
t157 = qJD(5) * t165;
t191 = -0.2e1 * t157 + t110;
t118 = qJDD(2) * pkin(2);
t86 = -g(1) * t121 - g(2) * t119;
t40 = -t128 * t86 + t131 * t198;
t146 = qJDD(3) - t40;
t35 = -t133 * qJ(3) - t118 + t146;
t33 = -qJDD(2) * pkin(8) + t35;
t31 = t130 * t33;
t141 = -qJDD(4) * pkin(4) - t132 * qJ(5) + t75 * t165 + qJDD(5) - t31;
t136 = -qJDD(4) * pkin(9) + t141;
t59 = t116 * t122 - t120 * t85;
t190 = pkin(5) * t80 + t136 + (pkin(5) * t163 + pkin(9) * t167 + t59) * t127;
t67 = t72 ^ 2;
t68 = t74 ^ 2;
t101 = qJD(6) + t165;
t96 = t101 ^ 2;
t188 = 0.2e1 * qJD(5);
t187 = pkin(4) + pkin(9);
t161 = qJDD(4) * qJ(5);
t173 = t130 * t59;
t181 = t127 * t33;
t22 = t173 + t181;
t79 = t156 + t160;
t91 = pkin(5) * t165 - qJD(4) * pkin(9);
t12 = t161 - pkin(9) * t107 - t79 * pkin(5) + (t188 + t91) * qJD(4) + t22 + t194;
t185 = t12 * t129;
t183 = t126 * t12;
t38 = t49 + t69;
t182 = t126 * t38;
t180 = t127 * t59;
t112 = qJDD(2) * qJ(3);
t41 = t128 * t198 + t131 * t86;
t153 = -t133 * pkin(2) + t112 + t41;
t145 = -t133 * pkin(8) + t153;
t138 = pkin(4) * t156 - qJ(5) * t201 + t145;
t13 = -pkin(5) * t107 - t91 * t165 + t187 * t79 + t138 + t191;
t176 = t129 * t13;
t175 = t129 * t38;
t155 = qJDD(4) * t126 - t129 * t79;
t140 = (-qJD(6) + t101) * t74 - t155;
t147 = qJDD(4) * t129 + t126 * t79;
t43 = -qJD(6) * t72 + t147;
t60 = t101 * t72;
t30 = t43 + t60;
t14 = t126 * t140 - t129 * t30;
t174 = t130 * t14;
t171 = t101 * t126;
t170 = t101 * t129;
t164 = qJD(6) + t101;
t158 = t130 * t49;
t4 = t126 * t13 - t129 * t190;
t82 = (t115 + t189) * qJDD(2);
t84 = -t108 - t107;
t154 = qJ(3) * t84 + t204 * t82;
t5 = t190 * t126 + t176;
t2 = t126 * t5 - t129 * t4;
t3 = t126 * t4 + t129 * t5;
t1 = t127 * t12 - t130 * t2;
t21 = -t31 + t180;
t8 = t127 * t22 - t130 * t21;
t151 = -t130 * (t107 - t132) + t179;
t148 = t127 * (-t108 + t132) - t172;
t144 = qJ(3) + t152;
t143 = t127 * t187 + qJ(3) - t168;
t142 = qJD(4) * t188 + t181 + t194;
t137 = -t142 - t161;
t17 = t141 + t180;
t135 = pkin(4) * t79 + t138;
t18 = t135 + t191;
t83 = t108 - t107;
t63 = (-qJDD(2) * t131 + t128 * t133) * t120;
t62 = (qJDD(2) * t128 + t131 * t133) * t120;
t58 = -t68 + t96;
t57 = t67 - t96;
t56 = t201 * t130;
t51 = (t79 + t156) * t127;
t50 = t122 * t59;
t48 = t68 - t67;
t47 = -t68 - t96;
t46 = -t127 * t81 - t130 * t78;
t45 = (t128 * t84 + t131 * t82) * t120;
t44 = -t96 - t67;
t42 = -qJD(6) * t74 - t155;
t36 = -t67 - t68;
t34 = t110 + t153;
t32 = t110 + t145;
t29 = t43 - t60;
t28 = -t164 * t72 + t147;
t25 = t164 * t74 + t155;
t24 = -t126 * t47 - t175;
t23 = t129 * t47 - t182;
t20 = t129 * t44 - t200;
t19 = t126 * t44 + t199;
t16 = -t137 + t173;
t15 = t126 * t30 + t129 * t140;
t11 = t127 * t28 - t130 * t23;
t9 = t127 * t25 - t130 * t19;
t7 = t127 * t36 - t174;
t6 = t127 * t16 - t130 * t17;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t116, 0, 0, 0, 0, 0, 0, -t63, -t62, 0, t50 + (t128 * t41 + t131 * t40) * t120, 0, 0, 0, 0, 0, 0, 0, t63, t62, t50 + (t128 * t34 - t131 * t35) * t120, 0, 0, 0, 0, 0, 0, -t205, t206, t45, t122 * (t127 * t21 + t130 * t22) + (t128 * t32 - t131 * t8) * t120, 0, 0, 0, 0, 0, 0, t45, t205, -t206, t122 * (t127 * t17 + t130 * t16) + (t128 * t18 - t131 * t6) * t120, 0, 0, 0, 0, 0, 0, t122 * (t127 * t19 + t130 * t25) + (t128 * t20 - t131 * t9) * t120, t122 * (t127 * t23 + t130 * t28) + (-t11 * t131 + t128 * t24) * t120, t122 * (t127 * t14 + t130 * t36) + (t128 * t15 - t131 * t7) * t120, t122 * (t12 * t130 + t127 * t2) + (-t1 * t131 + t128 * t3) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t40, -t41, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, -0.2e1 * t118 + t146, t110 + 0.2e1 * t112 + t41, -pkin(2) * t35 + qJ(3) * t34, t56, t46, -t148, t51, -t151, 0, qJ(3) * t78 + t127 * t32 - t203, qJ(3) * t81 + t130 * t32 + t202, t154 - t8, qJ(3) * t32 - t204 * t8, 0, t148, t151, t56, t46, t51, t130 * (-qJ(5) * t84 + t141) + (pkin(4) * t84 + t137) * t127 + t154, -t127 * t18 - t144 * t78 + t203, t130 * (-t135 + 0.2e1 * t157 - 0.2e1 * t162) - t144 * t81 - t202, t144 * t18 - t204 * t6, t158 - t127 * (-t126 * t43 - t74 * t170), t130 * t48 - t127 * (t126 * t25 - t129 * t29), t130 * t30 - t127 * (-t129 * t58 - t200), -t158 - t127 * (-t129 * t42 - t72 * t171), t130 * t140 - t127 * (-t126 * t57 - t175), t130 * t69 - t127 * (t126 * t72 + t129 * t74) * t101, t130 * (pkin(5) * t19 - t4) - t127 * (pkin(5) * t25 + t185) - t204 * t9 + t143 * t20, t130 * (-t176 - t126 * (pkin(9) * t99 + t136 + t180) - qJ(5) * t24 + (-t126 * (t80 + t105) + t23) * pkin(5)) - t127 * (pkin(5) * t28 - t187 * t24 - t183) + qJ(3) * t24 - t204 * t11, pkin(5) * t174 - t127 * (pkin(5) * t36 - t3) - t204 * t7 + t143 * t15, t143 * t3 + (-pkin(5) - t204) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t133, t35, 0, 0, 0, 0, 0, 0, t52, -t54, -t82, t8, 0, 0, 0, 0, 0, 0, -t82, -t52, t54, t6, 0, 0, 0, 0, 0, 0, t9, t11, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t83, t159, -t99, -t160, qJDD(4), -t21, -t22, 0, 0, qJDD(4), -t159, t160, t99, t83, -t99, (-pkin(4) * t130 - qJ(5) * t127) * qJDD(2), -pkin(4) * t88 - qJ(5) * t93 + t17, -pkin(4) * t95 + t173 + (qJDD(4) + t87) * qJ(5) + t142, -pkin(4) * t17 + qJ(5) * t16, t129 * t43 - t74 * t171, -t126 * t29 - t129 * t25, -t126 * t58 + t199, -t126 * t42 + t72 * t170, t129 * t57 - t182, (t126 * t74 - t129 * t72) * t101, qJ(5) * t25 - t187 * t19 + t183, qJ(5) * t28 - t187 * t23 + t185, qJ(5) * t36 - t187 * t14 - t2, qJ(5) * t12 - t187 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t88, t95, t17, 0, 0, 0, 0, 0, 0, t19, t23, t14, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t48, t30, -t49, t140, t69, -t4, -t5, 0, 0;];
tauJ_reg  = t10;
