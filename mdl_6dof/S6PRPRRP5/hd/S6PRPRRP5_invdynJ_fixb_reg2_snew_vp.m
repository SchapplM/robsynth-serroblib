% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 23:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:57:03
% EndTime: 2019-05-04 23:57:09
% DurationCPUTime: 1.83s
% Computational Cost: add. (5053->241), mult. (9600->320), div. (0->0), fcn. (6482->10), ass. (0->170)
t148 = sin(qJ(5));
t151 = cos(qJ(5));
t152 = cos(qJ(4));
t181 = qJD(2) * t152;
t113 = -t151 * qJD(4) + t148 * t181;
t133 = t152 * qJDD(2);
t149 = sin(qJ(4));
t177 = qJD(2) * qJD(4);
t173 = t149 * t177;
t119 = t133 - t173;
t84 = -t113 * qJD(5) + t148 * qJDD(4) + t151 * t119;
t131 = t149 * qJD(2) + qJD(5);
t99 = t131 * t113;
t67 = t84 + t99;
t219 = qJ(6) * t67;
t209 = -pkin(8) - pkin(2);
t172 = t152 * t177;
t176 = t149 * qJDD(2);
t118 = -t172 - t176;
t111 = qJDD(5) - t118;
t115 = t148 * qJD(4) + t151 * t181;
t89 = t115 * t113;
t211 = t111 - t89;
t218 = pkin(5) * t211;
t217 = t148 * t211;
t216 = t151 * t211;
t140 = -g(3) + qJDD(1);
t144 = sin(pkin(6));
t146 = cos(pkin(6));
t143 = sin(pkin(10));
t145 = cos(pkin(10));
t165 = t143 * g(1) - t145 * g(2);
t215 = t140 * t144 + t146 * t165;
t109 = t113 ^ 2;
t201 = t149 * pkin(4);
t166 = -t152 * pkin(9) + t201;
t116 = t166 * qJD(2);
t142 = qJDD(2) * pkin(2);
t154 = qJD(2) ^ 2;
t123 = -t145 * g(1) - t143 * g(2);
t150 = sin(qJ(2));
t153 = cos(qJ(2));
t78 = -t150 * t123 + t215 * t153;
t160 = qJDD(3) - t78;
t71 = -t154 * qJ(3) - t142 + t160;
t155 = -qJDD(2) * pkin(8) + t71;
t97 = t146 * t140 - t144 * t165;
t193 = t152 * t97;
t210 = qJD(4) ^ 2;
t38 = -t210 * pkin(4) + qJDD(4) * pkin(9) + t193 + (-qJD(2) * t116 + t155) * t149;
t163 = -t119 + t173;
t164 = -t118 + t172;
t136 = qJDD(2) * qJ(3);
t79 = t153 * t123 + t215 * t150;
t167 = 0.2e1 * qJD(3) * qJD(2) + t79;
t162 = t136 + t167;
t69 = t209 * t154 + t162;
t44 = t164 * pkin(4) + t163 * pkin(9) + t69;
t19 = t148 * t44 + t151 * t38;
t171 = -t151 * qJDD(4) + t148 * t119;
t83 = -t115 * qJD(5) - t171;
t93 = t131 * pkin(5) - t115 * qJ(6);
t159 = t83 * qJ(6) - 0.2e1 * qJD(6) * t113 - t131 * t93 + t19;
t110 = t115 ^ 2;
t130 = t131 ^ 2;
t86 = -t110 - t130;
t214 = -t159 + (t109 + t86) * pkin(5);
t212 = t84 - t99;
t63 = (qJD(5) - t131) * t115 + t171;
t85 = -t130 - t109;
t47 = t148 * t85 + t216;
t208 = pkin(4) * t47;
t76 = t111 + t89;
t196 = t148 * t76;
t53 = t151 * t86 - t196;
t207 = pkin(4) * t53;
t179 = qJD(6) * t115;
t105 = -0.2e1 * t179;
t18 = t148 * t38 - t151 * t44;
t157 = -t18 + t218 - t219;
t13 = t105 + t157;
t206 = pkin(5) * t13;
t205 = pkin(5) * t67;
t34 = -t148 * t63 - t151 * t67;
t204 = pkin(9) * t34;
t203 = pkin(9) * t47;
t202 = pkin(9) * t53;
t35 = t148 * t67 - t151 * t63;
t74 = -t109 - t110;
t200 = -pkin(4) * t74 + pkin(9) * t35;
t48 = t151 * t85 - t217;
t62 = (qJD(5) + t131) * t115 + t171;
t199 = -pkin(4) * t62 + pkin(9) * t48;
t194 = t151 * t76;
t54 = -t148 * t86 - t194;
t198 = -pkin(4) * t212 + pkin(9) * t54;
t49 = t149 * t97 - t152 * t155;
t37 = -qJDD(4) * pkin(4) - t210 * pkin(9) + t116 * t181 + t49;
t197 = t148 * t37;
t195 = t151 * t37;
t191 = qJ(6) * t148;
t190 = qJ(6) * t151;
t189 = t131 * t148;
t188 = t131 * t151;
t138 = t149 ^ 2;
t187 = t138 * t154;
t139 = t152 ^ 2;
t186 = t139 * t154;
t174 = t152 * t154 * t149;
t124 = qJDD(4) + t174;
t184 = t149 * t124;
t125 = qJDD(4) - t174;
t183 = t152 * t125;
t182 = t138 + t139;
t175 = t149 * t89;
t8 = t148 * t18 + t151 * t19;
t22 = t149 * t35 - t152 * t74;
t170 = qJ(3) * t34 + t209 * t22;
t26 = t149 * t48 - t152 * t62;
t169 = qJ(3) * t47 + t209 * t26;
t28 = t149 * t54 - t152 * t212;
t168 = qJ(3) * t53 + t209 * t28;
t7 = t148 * t19 - t151 * t18;
t50 = t149 * t155 + t193;
t23 = t149 * t50 - t152 * t49;
t156 = t157 + t218;
t20 = -t83 * pkin(5) - t109 * qJ(6) + t115 * t93 + qJDD(6) + t37;
t129 = -t186 - t210;
t128 = -t187 - t210;
t122 = t182 * t154;
t121 = t182 * qJDD(2);
t120 = t133 - 0.2e1 * t173;
t117 = 0.2e1 * t172 + t176;
t106 = 0.2e1 * t179;
t101 = (-qJDD(2) * t153 + t150 * t154) * t144;
t100 = (qJDD(2) * t150 + t153 * t154) * t144;
t95 = -t110 + t130;
t94 = t109 - t130;
t92 = t152 * t129 - t184;
t91 = t149 * t128 + t183;
t90 = t146 * t97;
t87 = t110 - t109;
t72 = (-t113 * t148 - t115 * t151) * t131;
t70 = -t154 * pkin(2) + t162;
t59 = t115 * t188 + t148 * t84;
t58 = t113 * t189 + t151 * t83;
t57 = t149 * t111 + t152 * (-t113 * t151 + t115 * t148) * t131;
t56 = t148 * t94 + t194;
t55 = t151 * t95 + t217;
t43 = t152 * (-t115 * t189 + t151 * t84) + t175;
t42 = t152 * (t113 * t188 - t148 * t83) - t175;
t39 = -pkin(5) * t212 - qJ(6) * t76;
t33 = -t148 * t62 + t151 * t212;
t30 = t152 * (t151 * t94 - t196) - t149 * t63;
t29 = t152 * (-t148 * t95 + t216) + t149 * t67;
t24 = t152 * (-t148 * t212 - t151 * t62) + t149 * t87;
t16 = -qJ(6) * t86 + t20;
t15 = -pkin(5) * t62 + qJ(6) * t85 - t20;
t14 = -t109 * pkin(5) + t159;
t12 = t146 * (t149 * t212 + t152 * t54) + (t150 * t53 - t153 * t28) * t144;
t11 = t106 - t157 + t219;
t10 = t146 * (t149 * t62 + t152 * t48) + (t150 * t47 - t153 * t26) * t144;
t9 = -qJ(6) * t63 + (-t109 - t74) * pkin(5) + t159;
t6 = t146 * (t149 * t74 + t152 * t35) + (t150 * t34 - t153 * t22) * t144;
t5 = -pkin(5) * t20 + qJ(6) * t14;
t4 = t149 * t8 - t152 * t37;
t3 = -t148 * t13 + t151 * t14;
t2 = t151 * t13 + t148 * t14;
t1 = t149 * t3 - t152 * t20;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t140, 0, 0, 0, 0, 0, 0, -t101, -t100, 0, t90 + (t150 * t79 + t153 * t78) * t144, 0, 0, 0, 0, 0, 0, 0, t101, t100, t90 + (t150 * t70 - t153 * t71) * t144, 0, 0, 0, 0, 0, 0, t146 * (-t149 * t125 + t152 * t128) + (t150 * t117 - t153 * t91) * t144, t146 * (-t152 * t124 - t149 * t129) + (t150 * t120 - t153 * t92) * t144, (t121 * t153 - t122 * t150) * t144, t146 * (t149 * t49 + t152 * t50) + (t150 * t69 - t153 * t23) * t144, 0, 0, 0, 0, 0, 0, t10, t12, t6, t146 * (t149 * t37 + t152 * t8) + (t150 * t7 - t153 * t4) * t144, 0, 0, 0, 0, 0, 0, t10, t12, t6, t146 * (t149 * t20 + t152 * t3) + (-t153 * t1 + t150 * t2) * t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t78, -t79, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, -0.2e1 * t142 + t160, 0.2e1 * t136 + t167, -pkin(2) * t71 + qJ(3) * t70, -t163 * t152, -t152 * t117 - t149 * t120, t183 - t149 * (-t186 + t210), t164 * t149, t152 * (t187 - t210) - t184, 0, qJ(3) * t117 + t149 * t69 + t209 * t91, qJ(3) * t120 + t152 * t69 + t209 * t92, -qJ(3) * t122 - t209 * t121 - t23, qJ(3) * t69 + t209 * t23, t43, t24, t29, t42, t30, t57, t152 * (t197 - t203) - t149 * (t18 - t208) + t169, t152 * (t195 - t202) - t149 * (t19 - t207) + t168, t152 * (-t7 - t204) + t34 * t201 + t170, t209 * t4 + (qJ(3) + t166) * t7, t43, t24, t29, t42, t30, t57, t152 * (-t148 * t15 - t190 * t211 - t203) - t149 * (t106 - t156 - t208) + t169, t152 * (-t148 * t39 + t151 * t16 - t202) - t149 * (-t207 - t214) + t168, t152 * (t151 * t11 - t148 * t9 - t204) - t149 * (-pkin(4) * t34 + t205) + t170, t152 * (-pkin(9) * t2 - t13 * t190 - t148 * t5) - t149 * (-pkin(4) * t2 - t206) + qJ(3) * t2 + t209 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t154, t71, 0, 0, 0, 0, 0, 0, t91, t92, -t121, t23, 0, 0, 0, 0, 0, 0, t26, t28, t22, t4, 0, 0, 0, 0, 0, 0, t26, t28, t22, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, (-t138 + t139) * t154, t133, -t174, -t176, qJDD(4), -t49, -t50, 0, 0, t59, t33, t55, t58, t56, t72, -t195 + t199, t197 + t198, t8 + t200, -pkin(4) * t37 + pkin(9) * t8, t59, t33, t55, t58, t56, t72, t151 * t15 - t191 * t211 + t199, t148 * t16 + t151 * t39 + t198, t148 * t11 + t151 * t9 + t200, -pkin(4) * t20 + pkin(9) * t3 - t13 * t191 + t151 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t87, t67, -t89, -t63, t111, -t18, -t19, 0, 0, t89, t87, t67, -t89, -t63, t111, t105 + t156, t214, -t205, t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t212, t74, t20;];
tauJ_reg  = t17;
