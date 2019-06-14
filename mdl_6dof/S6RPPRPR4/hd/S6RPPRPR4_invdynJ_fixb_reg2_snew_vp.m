% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:15:53
% EndTime: 2019-05-05 14:16:02
% DurationCPUTime: 3.57s
% Computational Cost: add. (12206->350), mult. (24702->487), div. (0->0), fcn. (14477->10), ass. (0->201)
t165 = sin(pkin(10));
t167 = cos(pkin(10));
t175 = cos(qJ(4));
t172 = sin(qJ(4));
t214 = qJD(1) * t172;
t133 = -t167 * t175 * qJD(1) + t165 * t214;
t134 = (-t175 * t165 - t172 * t167) * qJD(1);
t218 = t133 * t134;
t235 = qJDD(4) + t218;
t240 = t165 * t235;
t239 = t167 * t235;
t238 = 2 * qJD(5);
t171 = sin(qJ(6));
t206 = qJD(1) * qJD(4);
t199 = t175 * t206;
t205 = t172 * qJDD(1);
t140 = -t199 - t205;
t200 = t172 * t206;
t204 = t175 * qJDD(1);
t184 = t200 - t204;
t192 = t140 * t165 - t167 * t184;
t106 = qJDD(6) + t192;
t174 = cos(qJ(6));
t114 = -t174 * qJD(4) + t134 * t171;
t116 = qJD(4) * t171 + t134 * t174;
t91 = t116 * t114;
t234 = t106 - t91;
t237 = t171 * t234;
t236 = t174 * t234;
t212 = qJD(4) * t134;
t92 = t192 + t212;
t162 = t175 ^ 2;
t178 = qJD(1) ^ 2;
t157 = t162 * t178;
t177 = qJD(4) ^ 2;
t150 = -t157 - t177;
t163 = g(3) + qJDD(3);
t160 = qJDD(1) * qJ(2);
t173 = sin(qJ(1));
t176 = cos(qJ(1));
t190 = g(1) * t176 + g(2) * t173;
t186 = 0.2e1 * qJD(2) * qJD(1) - t190;
t185 = t160 + t186;
t233 = pkin(1) + pkin(2);
t121 = -t233 * t178 + t185;
t166 = sin(pkin(9));
t168 = cos(pkin(9));
t196 = t173 * g(1) - g(2) * t176;
t188 = qJDD(2) - t196;
t182 = -qJ(2) * t178 + t188;
t180 = -t233 * qJDD(1) + t182;
t89 = t168 * t121 + t166 * t180;
t183 = -pkin(3) * t178 - qJDD(1) * pkin(7) + t89;
t181 = t172 * t183;
t216 = t172 * t178;
t179 = -t181 - t140 * qJ(5) + qJDD(4) * pkin(4) + (pkin(4) * t216 - qJ(5) * t206 + t163) * t175;
t82 = t172 * t163 + t175 * t183;
t69 = pkin(4) * t150 - qJ(5) * t204 + t82;
t39 = t133 * t238 + t165 * t179 + t167 * t69;
t128 = -qJD(6) + t133;
t109 = t167 * t140 + t165 * t184;
t194 = -t174 * qJDD(4) + t109 * t171;
t57 = (qJD(6) + t128) * t116 + t194;
t112 = t114 ^ 2;
t113 = t116 ^ 2;
t127 = t128 ^ 2;
t131 = t133 ^ 2;
t132 = t134 ^ 2;
t232 = pkin(5) * t165;
t193 = t121 * t166 - t168 * t180;
t86 = qJDD(1) * pkin(3) - t178 * pkin(7) + t193;
t71 = qJDD(5) - (qJD(4) * pkin(4) + qJ(5) * t214) * t214 - t184 * pkin(4) - qJ(5) * t157 + t86;
t230 = t165 * t71;
t229 = t167 * t71;
t101 = -pkin(5) * t133 - pkin(8) * t134;
t195 = t165 * t69 - t167 * t179;
t27 = -qJDD(4) * pkin(5) - t177 * pkin(8) - (-(2 * qJD(5)) - t101) * t134 + t195;
t228 = t171 * t27;
t73 = t106 + t91;
t227 = t171 * t73;
t38 = t134 * t238 + t195;
t21 = t165 * t39 - t167 * t38;
t226 = t172 * t21;
t225 = t174 * t27;
t224 = t174 * t73;
t223 = qJDD(1) * pkin(1);
t104 = qJDD(4) - t218;
t222 = t104 * t165;
t221 = t104 * t167;
t220 = t128 * t171;
t219 = t128 * t174;
t151 = t175 * t216;
t146 = qJDD(4) + t151;
t217 = t172 * t146;
t147 = qJDD(4) - t151;
t215 = t175 * t147;
t213 = qJD(4) * t133;
t211 = qJD(4) * t165;
t210 = qJD(4) * t167;
t208 = qJD(6) - t128;
t203 = t165 * t91;
t202 = t167 * t91;
t201 = -pkin(5) * t167 - pkin(4);
t198 = qJ(2) * t166 + pkin(3);
t197 = qJ(2) * t168 - pkin(7);
t22 = t165 * t38 + t167 * t39;
t28 = -pkin(5) * t177 + qJDD(4) * pkin(8) + t101 * t133 + t39;
t191 = -t109 - t213;
t44 = pkin(5) * t92 + t191 * pkin(8) + t71;
t19 = t171 * t28 - t174 * t44;
t20 = t171 * t44 + t174 * t28;
t8 = t171 * t19 + t174 * t20;
t3 = t165 * t8 - t167 * t27;
t4 = t165 * t27 + t167 * t8;
t189 = t172 * t3 - t175 * t4;
t7 = t171 * t20 - t174 * t19;
t81 = -t175 * t163 + t181;
t51 = t172 * t81 + t175 * t82;
t187 = -qJDD(4) * t171 - t109 * t174;
t93 = -t192 + t212;
t141 = 0.2e1 * t200 - t204;
t80 = -qJD(6) * t114 - t187;
t161 = t172 ^ 2;
t155 = t161 * t178;
t149 = -t155 - t177;
t145 = t155 + t157;
t144 = (-t161 - t162) * qJDD(1);
t143 = qJDD(1) * t168 + t166 * t178;
t142 = -t166 * qJDD(1) + t168 * t178;
t139 = 0.2e1 * t199 + t205;
t130 = -t182 + t223;
t124 = -t132 - t177;
t123 = -t132 + t177;
t122 = t131 - t177;
t118 = -t149 * t172 - t215;
t117 = t150 * t175 - t217;
t110 = t144 * t166 + t145 * t168;
t102 = -t177 - t131;
t100 = t114 * t128;
t99 = t118 * t166 + t139 * t168;
t98 = t117 * t166 + t141 * t168;
t97 = -t113 + t127;
t96 = t112 - t127;
t95 = t109 - t213;
t90 = -t131 - t132;
t87 = t113 - t112;
t85 = -t113 - t127;
t84 = -t124 * t165 - t221;
t83 = t124 * t167 - t222;
t79 = -qJD(6) * t116 - t194;
t78 = -t127 - t112;
t77 = t112 + t113;
t76 = t102 * t167 - t240;
t75 = t102 * t165 + t239;
t70 = (t114 * t174 - t116 * t171) * t128;
t65 = t165 * t95 + t167 * t93;
t64 = t165 * t93 - t167 * t95;
t63 = t166 * t89 - t168 * t193;
t62 = t208 * t114 + t187;
t61 = -t100 + t80;
t60 = t100 + t80;
t58 = -t208 * t116 - t194;
t56 = t116 * t220 + t174 * t80;
t55 = -t114 * t219 - t171 * t79;
t54 = -t172 * t83 + t175 * t84;
t53 = t174 * t96 - t227;
t52 = -t171 * t97 + t236;
t50 = -t171 * t85 - t224;
t49 = t174 * t85 - t227;
t48 = t174 * t78 - t237;
t47 = t171 * t78 + t236;
t46 = -t172 * t75 + t175 * t76;
t45 = t166 * t54 + t168 * t191;
t42 = t166 * t51 - t168 * t86;
t41 = -t172 * t64 + t175 * t65;
t40 = t166 * t46 - t168 * t92;
t36 = t171 * t61 - t174 * t57;
t35 = -t171 * t60 + t174 * t58;
t34 = -t171 * t57 - t174 * t61;
t33 = -t165 * t62 + t167 * t50;
t32 = t165 * t50 + t167 * t62;
t31 = -t165 * t58 + t167 * t48;
t30 = t165 * t48 + t167 * t58;
t29 = t166 * t41 - t168 * t90;
t26 = -t165 * t77 + t167 * t36;
t25 = t165 * t36 + t167 * t77;
t24 = -pkin(8) * t49 + t225;
t23 = -pkin(8) * t47 + t228;
t17 = -t172 * t32 + t175 * t33;
t16 = -t172 * t30 + t175 * t31;
t15 = -t172 * t25 + t175 * t26;
t14 = -pkin(5) * t49 + t20;
t13 = -pkin(5) * t47 + t19;
t12 = t166 * t17 - t168 * t49;
t11 = t16 * t166 - t168 * t47;
t10 = t15 * t166 - t168 * t34;
t9 = t175 * t22 - t226;
t6 = t166 * t9 - t168 * t71;
t5 = -pkin(8) * t34 - t7;
t1 = -t166 * t189 - t168 * t7;
t2 = [0, 0, 0, 0, 0, qJDD(1), t196, t190, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -t188 + 0.2e1 * t223, 0, 0.2e1 * t160 + t186, pkin(1) * t130 + qJ(2) * (-pkin(1) * t178 + t185), 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t142 + t233 * t143 + t193, qJ(2) * t143 + t233 * t142 + t89, 0, qJ(2) * (t166 * t193 + t168 * t89) - t233 * t63, (-t140 + t199) * t172, t139 * t175 - t141 * t172, -t217 - t175 * (-t155 + t177), -t141 * t175, -t172 * (t157 - t177) - t215, 0, qJ(2) * (t117 * t168 - t141 * t166) + t175 * t86 - pkin(3) * t141 - pkin(7) * t117 - t233 * t98, qJ(2) * (t118 * t168 - t139 * t166) - t172 * t86 - pkin(3) * t139 - pkin(7) * t118 - t233 * t99, qJ(2) * (t144 * t168 - t145 * t166) - pkin(3) * t145 - pkin(7) * t144 - t233 * t110 - t51, qJ(2) * (t166 * t86 + t168 * t51) + pkin(3) * t86 - pkin(7) * t51 - t233 * t42, -t172 * (t109 * t167 - t134 * t211) - t175 * (t109 * t165 + t134 * t210), -t172 * (t165 * t191 - t167 * t92) - t175 * (-t165 * t92 - t167 * t191), -t172 * (-t123 * t165 + t239) - t175 * (t123 * t167 + t240), -t172 * (-t133 * t210 + t165 * t192) - t175 * (-t133 * t211 - t167 * t192), -t172 * (t122 * t167 - t222) - t175 * (t122 * t165 + t221), (-t172 * (t133 * t167 + t134 * t165) - t175 * (t133 * t165 - t134 * t167)) * qJD(4), qJ(2) * (t166 * t92 + t168 * t46) - t172 * (-qJ(5) * t75 + t230) - t175 * (-pkin(4) * t92 + qJ(5) * t76 - t229) + pkin(3) * t92 - pkin(7) * t46 - t233 * t40, qJ(2) * (-t166 * t191 + t168 * t54) - t172 * (-qJ(5) * t83 + t229) - t175 * (pkin(4) * t191 + qJ(5) * t84 + t230) - pkin(3) * t191 - pkin(7) * t54 - t233 * t45, qJ(2) * (t166 * t90 + t168 * t41) - t172 * (-qJ(5) * t64 - t21) - t175 * (-pkin(4) * t90 + qJ(5) * t65 + t22) + pkin(3) * t90 - pkin(7) * t41 - t233 * t29, qJ(2) * (t166 * t71 + t168 * t9) + qJ(5) * t226 - t175 * (-pkin(4) * t71 + qJ(5) * t22) + pkin(3) * t71 - pkin(7) * t9 - t233 * t6, -t172 * (t167 * t56 + t203) - t175 * (t165 * t56 - t202), -t172 * (t165 * t87 + t167 * t35) - t175 * (t165 * t35 - t167 * t87), -t172 * (t165 * t61 + t167 * t52) - t175 * (t165 * t52 - t167 * t61), -t172 * (t167 * t55 - t203) - t175 * (t165 * t55 + t202), -t172 * (-t165 * t57 + t167 * t53) - t175 * (t165 * t53 + t167 * t57), -t172 * (t106 * t165 + t167 * t70) - t175 * (-t106 * t167 + t165 * t70), qJ(2) * (t16 * t168 + t166 * t47) - t172 * (-qJ(5) * t30 - t13 * t165 + t167 * t23) - t175 * (-pkin(4) * t47 + qJ(5) * t31 + t13 * t167 + t165 * t23) + pkin(3) * t47 - pkin(7) * t16 - t233 * t11, qJ(2) * (t166 * t49 + t168 * t17) - t172 * (-qJ(5) * t32 - t14 * t165 + t167 * t24) - t175 * (-pkin(4) * t49 + qJ(5) * t33 + t14 * t167 + t165 * t24) + pkin(3) * t49 - pkin(7) * t17 - t233 * t12, -t172 * (-qJ(5) * t25 + t167 * t5) - t175 * (qJ(5) * t26 + t165 * t5) + t197 * t15 + (-t172 * t232 - t175 * t201 + t198) * t34 - t233 * t10, -t233 * t1 + (-t172 * (-pkin(8) * t167 + t232) - t175 * (-pkin(8) * t165 + t201) + t198) * t7 + (-t197 + qJ(5)) * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t178, -t130, 0, 0, 0, 0, 0, 0, -t143, -t142, 0, t63, 0, 0, 0, 0, 0, 0, t98, t99, t110, t42, 0, 0, 0, 0, 0, 0, t40, t45, t29, t6, 0, 0, 0, 0, 0, 0, t11, t12, t10, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, 0, 0, 0, 0, 0, 0, t146 * t175 + t150 * t172, -t147 * t172 + t149 * t175, 0, t172 * t82 - t175 * t81, 0, 0, 0, 0, 0, 0, t172 * t76 + t175 * t75, t172 * t84 + t175 * t83, t172 * t65 + t175 * t64, t172 * t22 + t175 * t21, 0, 0, 0, 0, 0, 0, t172 * t31 + t175 * t30, t172 * t33 + t175 * t32, t172 * t26 + t175 * t25, t172 * t4 + t175 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, t155 - t157, -t205, t151, -t204, qJDD(4), -t81, -t82, 0, 0, -t218, t132 - t131, t95, t218, t93, qJDD(4), pkin(4) * t75 - t38, pkin(4) * t83 - t39, pkin(4) * t64, pkin(4) * t21, -t116 * t219 + t171 * t80, t171 * t58 + t174 * t60, t174 * t97 + t237, -t114 * t220 + t174 * t79, t171 * t96 + t224, (t114 * t171 + t116 * t174) * t128, pkin(4) * t30 + pkin(5) * t58 + pkin(8) * t48 - t225, pkin(4) * t32 + pkin(5) * t62 + pkin(8) * t50 + t228, pkin(4) * t25 + pkin(5) * t77 + pkin(8) * t36 + t8, pkin(4) * t3 - pkin(5) * t27 + pkin(8) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, -t191, t90, t71, 0, 0, 0, 0, 0, 0, t47, t49, t34, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t87, t61, -t91, -t57, t106, -t19, -t20, 0, 0;];
tauJ_reg  = t2;
