% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 01:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRRR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:34:29
% EndTime: 2019-05-05 01:34:37
% DurationCPUTime: 2.73s
% Computational Cost: add. (12030->325), mult. (22941->466), div. (0->0), fcn. (16056->12), ass. (0->207)
t172 = sin(qJ(6));
t178 = cos(qJ(4));
t203 = qJD(2) * qJD(4);
t157 = t178 * t203;
t174 = sin(qJ(4));
t158 = t174 * qJDD(2);
t142 = -t158 - t157;
t135 = qJDD(5) - t142;
t132 = qJDD(6) + t135;
t173 = sin(qJ(5));
t177 = cos(qJ(5));
t207 = qJD(2) * t178;
t137 = -t177 * qJD(4) + t173 * t207;
t139 = qJD(4) * t173 + t177 * t207;
t176 = cos(qJ(6));
t113 = t176 * t137 + t139 * t172;
t115 = -t137 * t172 + t139 * t176;
t85 = t115 * t113;
t237 = -t85 + t132;
t243 = t172 * t237;
t118 = t139 * t137;
t236 = -t118 + t135;
t242 = t173 * t236;
t241 = t176 * t237;
t240 = t177 * t236;
t169 = sin(pkin(6));
t170 = cos(pkin(6));
t221 = sin(pkin(11));
t222 = cos(pkin(11));
t184 = t221 * g(1) - t222 * g(2);
t209 = -g(3) + qJDD(1);
t239 = t169 * t209 + t170 * t184;
t160 = t178 * qJDD(2);
t199 = t174 * t203;
t143 = t160 - t199;
t190 = -t173 * qJDD(4) - t177 * t143;
t109 = -qJD(5) * t137 - t190;
t196 = -t177 * qJDD(4) + t173 * t143;
t187 = qJD(5) * t139 + t196;
t64 = -t113 * qJD(6) + t176 * t109 - t172 * t187;
t156 = qJD(2) * t174 + qJD(5);
t152 = qJD(6) + t156;
t99 = t152 * t113;
t238 = -t99 + t64;
t126 = t156 * t137;
t90 = t109 + t126;
t197 = t172 * t109 + t176 * t187;
t48 = (qJD(6) - t152) * t115 + t197;
t86 = (qJD(5) - t156) * t139 + t196;
t111 = t113 ^ 2;
t112 = t115 ^ 2;
t235 = t137 ^ 2;
t134 = t139 ^ 2;
t151 = t152 ^ 2;
t155 = t156 ^ 2;
t234 = qJD(4) ^ 2;
t194 = pkin(4) * t174 - pkin(9) * t178;
t140 = t194 * qJD(2);
t168 = qJDD(2) * pkin(2);
t180 = qJD(2) ^ 2;
t175 = sin(qJ(2));
t179 = cos(qJ(2));
t185 = -t222 * g(1) - t221 * g(2);
t105 = -t175 * t185 + t239 * t179;
t182 = qJDD(3) - t105;
t96 = -t180 * qJ(3) - t168 + t182;
t181 = -qJDD(2) * pkin(8) + t96;
t125 = -t169 * t184 + t170 * t209;
t211 = t178 * t125;
t67 = -t234 * pkin(4) + qJDD(4) * pkin(9) + t211 + (-qJD(2) * t140 + t181) * t174;
t192 = -t143 + t199;
t193 = -t142 + t157;
t163 = qJDD(2) * qJ(3);
t106 = t239 * t175 + t179 * t185;
t195 = 0.2e1 * qJD(3) * qJD(2) + t106;
t191 = t163 + t195;
t232 = -pkin(8) - pkin(2);
t93 = t232 * t180 + t191;
t71 = t193 * pkin(4) + t192 * pkin(9) + t93;
t36 = t173 * t67 - t177 * t71;
t31 = t236 * pkin(5) - t90 * pkin(10) - t36;
t122 = pkin(5) * t156 - pkin(10) * t139;
t37 = t173 * t71 + t177 * t67;
t34 = -t235 * pkin(5) - t187 * pkin(10) - t156 * t122 + t37;
t13 = t172 * t34 - t176 * t31;
t14 = t172 * t31 + t176 * t34;
t7 = -t13 * t176 + t14 * t172;
t233 = pkin(5) * t7;
t51 = t99 + t64;
t26 = -t172 * t48 - t176 * t51;
t231 = pkin(5) * t26;
t230 = t173 * t7;
t229 = t177 * t7;
t79 = t125 * t174 - t178 * t181;
t66 = -qJDD(4) * pkin(4) - t234 * pkin(9) + t140 * t207 + t79;
t38 = t187 * pkin(5) - t235 * pkin(10) + t122 * t139 + t66;
t228 = t172 * t38;
t75 = t85 + t132;
t227 = t172 * t75;
t226 = t173 * t66;
t225 = t176 * t38;
t224 = t176 * t75;
t223 = t177 * t66;
t103 = t118 + t135;
t220 = t103 * t173;
t219 = t103 * t177;
t218 = t152 * t172;
t217 = t152 * t176;
t216 = t156 * t173;
t215 = t156 * t177;
t165 = t174 ^ 2;
t214 = t165 * t180;
t166 = t178 ^ 2;
t213 = t166 * t180;
t200 = t178 * t180 * t174;
t147 = qJDD(4) + t200;
t212 = t174 * t147;
t148 = qJDD(4) - t200;
t210 = t178 * t148;
t208 = t165 + t166;
t205 = qJD(5) + t156;
t202 = t174 * t85;
t201 = t174 * t118;
t8 = t13 * t172 + t176 * t14;
t21 = t173 * t36 + t177 * t37;
t20 = t173 * t37 - t177 * t36;
t80 = t174 * t181 + t211;
t40 = t174 * t80 - t178 * t79;
t81 = -t151 - t111;
t41 = t172 * t81 + t241;
t189 = pkin(5) * t41 - t13;
t188 = qJ(3) + t194;
t94 = -t112 - t151;
t54 = t176 * t94 - t227;
t186 = pkin(5) * t54 - t14;
t154 = -t213 - t234;
t153 = -t214 - t234;
t146 = t208 * t180;
t145 = t208 * qJDD(2);
t144 = t160 - 0.2e1 * t199;
t141 = t158 + 0.2e1 * t157;
t128 = (-qJDD(2) * t179 + t175 * t180) * t169;
t127 = (qJDD(2) * t175 + t179 * t180) * t169;
t124 = -t134 + t155;
t123 = -t155 + t235;
t121 = t154 * t178 - t212;
t120 = t153 * t174 + t210;
t119 = t170 * t125;
t117 = t134 - t235;
t116 = -t134 - t155;
t110 = -t155 - t235;
t101 = t134 + t235;
t98 = -t112 + t151;
t97 = t111 - t151;
t95 = -t180 * pkin(2) + t191;
t91 = t205 * t137 + t190;
t89 = t109 - t126;
t87 = -t205 * t139 - t196;
t84 = t112 - t111;
t83 = -t116 * t173 - t219;
t82 = t116 * t177 - t220;
t78 = t110 * t177 - t242;
t77 = t110 * t173 + t240;
t73 = (-t113 * t176 + t115 * t172) * t152;
t72 = (-t113 * t172 - t115 * t176) * t152;
t68 = -t111 - t112;
t63 = -qJD(6) * t115 - t197;
t62 = t173 * t90 - t177 * t86;
t61 = -t173 * t86 - t177 * t90;
t60 = t176 * t97 - t227;
t59 = -t172 * t98 + t241;
t58 = t172 * t97 + t224;
t57 = t176 * t98 + t243;
t56 = t174 * t83 + t178 * t91;
t55 = -t172 * t94 - t224;
t53 = t174 * t78 + t178 * t87;
t47 = (qJD(6) + t152) * t115 + t197;
t46 = -t115 * t218 + t176 * t64;
t45 = t115 * t217 + t172 * t64;
t44 = t113 * t217 - t172 * t63;
t43 = t113 * t218 + t176 * t63;
t42 = t176 * t81 - t243;
t39 = t101 * t178 + t174 * t62;
t33 = -t173 * t54 + t177 * t55;
t32 = t173 * t55 + t177 * t54;
t29 = -t172 * t238 - t176 * t47;
t28 = t172 * t51 - t176 * t48;
t27 = -t172 * t47 + t176 * t238;
t25 = -pkin(10) * t54 + t225;
t24 = -t173 * t41 + t177 * t42;
t23 = t173 * t42 + t177 * t41;
t22 = -pkin(10) * t41 + t228;
t19 = t174 * t33 - t178 * t238;
t18 = -pkin(5) * t238 + pkin(10) * t55 + t228;
t17 = t174 * t24 - t178 * t47;
t16 = -pkin(5) * t47 + pkin(10) * t42 - t225;
t15 = t174 * t21 - t178 * t66;
t11 = -t173 * t26 + t177 * t28;
t10 = t173 * t28 + t177 * t26;
t9 = t11 * t174 - t178 * t68;
t6 = -pkin(5) * t38 + pkin(10) * t8;
t5 = -pkin(10) * t26 - t7;
t4 = -pkin(5) * t68 + pkin(10) * t28 + t8;
t3 = t177 * t8 - t230;
t2 = t173 * t8 + t229;
t1 = t174 * t3 - t178 * t38;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t209, 0, 0, 0, 0, 0, 0, -t128, -t127, 0, t119 + (t105 * t179 + t106 * t175) * t169, 0, 0, 0, 0, 0, 0, 0, t128, t127, t119 + (t175 * t95 - t179 * t96) * t169, 0, 0, 0, 0, 0, 0, t170 * (-t148 * t174 + t153 * t178) + (-t120 * t179 + t141 * t175) * t169, t170 * (-t147 * t178 - t154 * t174) + (-t121 * t179 + t144 * t175) * t169, (t145 * t179 - t146 * t175) * t169, t170 * (t174 * t79 + t178 * t80) + (t175 * t93 - t179 * t40) * t169, 0, 0, 0, 0, 0, 0, t170 * (-t174 * t87 + t178 * t78) + (t175 * t77 - t179 * t53) * t169, t170 * (-t174 * t91 + t178 * t83) + (t175 * t82 - t179 * t56) * t169, t170 * (-t101 * t174 + t178 * t62) + (t175 * t61 - t179 * t39) * t169, t170 * (t174 * t66 + t178 * t21) + (-t15 * t179 + t175 * t20) * t169, 0, 0, 0, 0, 0, 0, t170 * (t174 * t47 + t178 * t24) + (-t17 * t179 + t175 * t23) * t169, t170 * (t174 * t238 + t178 * t33) + (t175 * t32 - t179 * t19) * t169, t170 * (t11 * t178 + t174 * t68) + (t10 * t175 - t179 * t9) * t169, t170 * (t174 * t38 + t178 * t3) + (-t1 * t179 + t175 * t2) * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t105, -t106, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, -0.2e1 * t168 + t182, 0.2e1 * t163 + t195, -pkin(2) * t96 + qJ(3) * t95, -t192 * t178, -t141 * t178 - t144 * t174, t210 - t174 * (-t213 + t234), t193 * t174, t178 * (t214 - t234) - t212, 0, qJ(3) * t141 + t232 * t120 + t174 * t93, qJ(3) * t144 + t232 * t121 + t178 * t93, -qJ(3) * t146 - t232 * t145 - t40, qJ(3) * t93 + t232 * t40, t178 * (t109 * t177 - t139 * t216) + t201, t178 * (-t173 * t89 + t177 * t87) + t174 * t117, t178 * (-t124 * t173 + t240) + t174 * t90, t178 * (t137 * t215 + t173 * t187) - t201, t178 * (t123 * t177 - t220) - t174 * t86, t174 * t135 + t178 * (-t137 * t177 + t139 * t173) * t156, t178 * (-pkin(9) * t77 + t226) - t174 * (-pkin(4) * t77 + t36) + qJ(3) * t77 + t232 * t53, t178 * (-pkin(9) * t82 + t223) - t174 * (-pkin(4) * t82 + t37) + qJ(3) * t82 + t232 * t56, -t178 * t20 + t188 * t61 + t232 * t39, t232 * t15 + t188 * t20, t178 * (-t173 * t45 + t177 * t46) + t202, t178 * (-t173 * t27 + t177 * t29) + t174 * t84, t178 * (-t173 * t57 + t177 * t59) + t174 * t51, t178 * (-t173 * t43 + t177 * t44) - t202, t178 * (-t173 * t58 + t177 * t60) - t174 * t48, t178 * (-t173 * t72 + t177 * t73) + t174 * t132, t178 * (-pkin(9) * t23 - t16 * t173 + t177 * t22) - t174 * (-pkin(4) * t23 - t189) + qJ(3) * t23 + t232 * t17, t178 * (-pkin(9) * t32 - t173 * t18 + t177 * t25) - t174 * (-pkin(4) * t32 - t186) + qJ(3) * t32 + t232 * t19, t178 * (-pkin(9) * t10 - t173 * t4 + t177 * t5) - t174 * (-pkin(4) * t10 - t231) + qJ(3) * t10 + t232 * t9, t178 * (-pkin(9) * t2 - pkin(10) * t229 - t173 * t6) - t174 * (-pkin(4) * t2 - t233) + qJ(3) * t2 + t232 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t180, t96, 0, 0, 0, 0, 0, 0, t120, t121, -t145, t40, 0, 0, 0, 0, 0, 0, t53, t56, t39, t15, 0, 0, 0, 0, 0, 0, t17, t19, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, (-t165 + t166) * t180, t160, -t200, -t158, qJDD(4), -t79, -t80, 0, 0, t109 * t173 + t139 * t215, t173 * t87 + t177 * t89, t124 * t177 + t242, t137 * t216 - t177 * t187, t123 * t173 + t219, (-t137 * t173 - t139 * t177) * t156, pkin(4) * t87 + pkin(9) * t78 - t223, pkin(4) * t91 + pkin(9) * t83 + t226, pkin(4) * t101 + pkin(9) * t62 + t21, -pkin(4) * t66 + pkin(9) * t21, t173 * t46 + t177 * t45, t173 * t29 + t177 * t27, t173 * t59 + t177 * t57, t173 * t44 + t177 * t43, t173 * t60 + t177 * t58, t173 * t73 + t177 * t72, -pkin(4) * t47 + pkin(9) * t24 + t16 * t177 + t173 * t22, -pkin(4) * t238 + pkin(9) * t33 + t173 * t25 + t177 * t18, -pkin(4) * t68 + pkin(9) * t11 + t173 * t5 + t177 * t4, -pkin(4) * t38 + pkin(9) * t3 - pkin(10) * t230 + t177 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t117, t90, -t118, -t86, t135, -t36, -t37, 0, 0, t85, t84, t51, -t85, -t48, t132, t189, t186, t231, t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t84, t51, -t85, -t48, t132, -t13, -t14, 0, 0;];
tauJ_reg  = t12;
