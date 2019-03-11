% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tau_reg [6x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:27:59
% EndTime: 2019-03-08 19:28:06
% DurationCPUTime: 2.59s
% Computational Cost: add. (2809->335), mult. (6665->478), div. (0->0), fcn. (5621->16), ass. (0->191)
t145 = sin(pkin(12));
t155 = sin(qJ(4));
t212 = qJD(2) * t155;
t149 = cos(pkin(12));
t158 = cos(qJ(4));
t222 = t149 * t158;
t96 = qJD(2) * t222 - t145 * t212;
t95 = qJD(6) - t96;
t256 = t95 - qJD(6);
t255 = qJ(5) + pkin(8);
t156 = sin(qJ(2));
t159 = cos(qJ(2));
t148 = sin(pkin(6));
t213 = qJD(2) * t148;
t198 = qJD(1) * t213;
t207 = qJDD(1) * t148;
t254 = t156 * t207 + t159 * t198;
t146 = sin(pkin(11));
t150 = cos(pkin(11));
t218 = t159 * t150;
t104 = t146 * t156 - t218;
t147 = sin(pkin(10));
t151 = cos(pkin(10));
t152 = cos(pkin(6));
t181 = t146 * t159 + t150 * t156;
t93 = t181 * t152;
t52 = t104 * t147 - t151 * t93;
t55 = t104 * t151 + t147 * t93;
t253 = g(1) * t55 + g(2) * t52;
t157 = cos(qJ(6));
t154 = sin(qJ(6));
t211 = qJD(6) * t154;
t205 = t158 * qJDD(2);
t206 = t155 * qJDD(2);
t184 = -t145 * t206 + t149 * t205;
t105 = t145 * t158 + t149 * t155;
t97 = t105 * qJD(4);
t65 = qJD(2) * t97 + qJDD(6) - t184;
t177 = -t157 * t65 + t95 * t211;
t92 = t181 * t148;
t71 = t152 * t158 - t155 * t92;
t72 = t152 * t155 + t158 * t92;
t26 = t145 * t71 + t149 * t72;
t225 = t148 * t156;
t91 = t146 * t225 - t148 * t218;
t84 = t91 * t157;
t252 = -t154 * t26 + t84;
t210 = t155 * qJD(4);
t214 = qJD(1) * t148;
t200 = t156 * t214;
t114 = t150 * t200;
t199 = t159 * t214;
t86 = t146 * t199 + t114;
t251 = pkin(4) * t210 - t86;
t124 = t152 * qJDD(1) + qJDD(3);
t112 = t158 * t124;
t126 = qJD(1) * t152 + qJD(3);
t121 = t159 * t207;
t90 = qJDD(2) * pkin(2) - t156 * t198 + t121;
t48 = t146 * t90 + t254 * t150;
t44 = qJDD(2) * pkin(8) + t48;
t167 = qJ(5) * qJDD(2) + qJD(2) * qJD(5) + qJD(4) * t126 + t44;
t109 = qJD(2) * pkin(2) + t199;
t76 = t146 * t109 + t114;
t190 = t255 * qJD(2) + t76;
t180 = t190 * qJD(4);
t11 = qJDD(4) * pkin(4) - t167 * t155 - t158 * t180 + t112;
t12 = (t124 - t180) * t155 + t167 * t158;
t3 = t11 * t149 - t12 * t145;
t1 = -qJDD(4) * pkin(5) - t3;
t130 = pkin(4) * t145 + pkin(9);
t142 = qJ(4) + pkin(12);
t137 = sin(t142);
t138 = cos(t142);
t226 = t148 * t151;
t228 = t147 * t148;
t98 = t105 * qJD(2);
t250 = (pkin(4) * t212 + pkin(5) * t98 - pkin(9) * t96 + qJD(6) * t130) * t95 + g(1) * (t137 * t55 + t138 * t228) + g(2) * (t137 * t52 - t138 * t226) + g(3) * (-t137 * t92 + t138 * t152) + t1;
t208 = qJD(2) * qJD(4);
t196 = t158 * t208;
t197 = t155 * t208;
t67 = t105 * qJDD(2) - t145 * t197 + t149 * t196;
t80 = qJD(4) * t154 + t157 * t98;
t24 = t80 * qJD(6) - t157 * qJDD(4) + t154 * t67;
t103 = t145 * t155 - t222;
t100 = t103 * qJD(4);
t46 = t126 * t155 + t190 * t158;
t235 = t145 * t46;
t45 = t158 * t126 - t190 * t155;
t40 = qJD(4) * pkin(4) + t45;
t17 = t149 * t40 - t235;
t15 = -qJD(4) * pkin(5) - t17;
t175 = -g(3) * t92 + t253;
t4 = t145 * t11 + t149 * t12;
t2 = qJDD(4) * pkin(9) + t4;
t135 = pkin(4) * t158 + pkin(3);
t113 = t146 * t200;
t75 = t109 * t150 - t113;
t64 = -t135 * qJD(2) + qJD(5) - t75;
t22 = -pkin(5) * t96 - pkin(9) * t98 + t64;
t131 = pkin(2) * t146 + pkin(8);
t217 = qJ(5) + t131;
t188 = qJD(4) * t217;
t169 = -qJD(5) * t155 - t158 * t188;
t85 = qJD(5) * t158 - t155 * t188;
t89 = t150 * t199 - t113;
t237 = t103 * t89 + t145 * t169 + t149 * t85;
t245 = pkin(2) * t150;
t179 = -t135 - t245;
t58 = pkin(5) * t103 - pkin(9) * t105 + t179;
t102 = t217 * t158;
t191 = t217 * t155;
t62 = t149 * t102 - t145 * t191;
t249 = -(qJD(6) * t22 + t2) * t103 + t1 * t105 - t15 * t100 + (-qJD(6) * t58 - t237) * t95 - t62 * t65 + t175;
t219 = t157 * t100;
t248 = -t177 * t105 - t95 * t219;
t172 = t104 * t152;
t53 = -t147 * t181 - t151 * t172;
t56 = t147 * t172 - t151 * t181;
t174 = g(1) * t56 + g(2) * t53 - g(3) * t91;
t247 = -t174 * t138 + t58 * t65;
t246 = g(3) * t71;
t209 = t157 * qJD(4);
t78 = t154 * t98 - t209;
t243 = t78 * t95;
t242 = t78 * t98;
t241 = t80 * t95;
t240 = t80 * t98;
t23 = qJD(6) * t209 + t154 * qJDD(4) + t157 * t67 - t98 * t211;
t239 = t23 * t103 + t80 * t97;
t38 = t149 * t46;
t18 = t145 * t40 + t38;
t238 = -t105 * t89 + t145 * t85 - t149 * t169;
t236 = pkin(5) * t97 + pkin(9) * t100 + t251;
t234 = t15 * t105;
t233 = t154 * t23;
t231 = t154 * t65;
t230 = t154 * t91;
t229 = t154 * t95;
t192 = t157 * t95;
t227 = t147 * t156;
t224 = t148 * t158;
t223 = t148 * t159;
t221 = t152 * t156;
t220 = t152 * t159;
t216 = qJDD(1) - g(3);
t143 = t155 ^ 2;
t215 = -t158 ^ 2 + t143;
t201 = t151 * t220;
t47 = -t254 * t146 + t150 * t90;
t186 = g(1) * t147 - g(2) * t151;
t185 = -t103 * t24 - t78 * t97;
t16 = qJD(4) * pkin(9) + t18;
t6 = t154 * t22 + t157 * t16;
t183 = t154 * t16 - t157 * t22;
t182 = t157 * t26 + t230;
t178 = t96 * t229 - t177;
t176 = -t147 * t220 - t151 * t156;
t173 = -g(1) * t228 + g(2) * t226 - g(3) * t152;
t73 = -qJD(2) * pkin(3) - t75;
t170 = -t73 * qJD(2) - t253 - t44;
t20 = t149 * t45 - t235;
t168 = -t130 * t65 + (t15 + t20) * t95;
t133 = -pkin(3) - t245;
t166 = -qJDD(4) * t131 + (qJD(2) * t133 + t73 + t89) * qJD(4);
t165 = -g(1) * t176 - g(3) * t223;
t164 = (-qJD(6) * t192 - t231) * t105 + t100 * t229;
t27 = pkin(4) * t197 - t135 * qJDD(2) + qJDD(5) - t47;
t160 = qJD(4) ^ 2;
t163 = -qJD(2) * t86 + t131 * t160 + t174 - t47 + (-pkin(3) + t133) * qJDD(2);
t161 = qJD(2) ^ 2;
t132 = -pkin(4) * t149 - pkin(5);
t118 = pkin(2) * t201;
t117 = qJDD(4) * t158 - t155 * t160;
t116 = qJDD(4) * t155 + t158 * t160;
t88 = t104 * t213;
t87 = qJD(2) * t92;
t69 = t137 * t152 + t138 * t92;
t66 = -qJD(4) * t98 + t184;
t61 = t102 * t145 + t149 * t191;
t34 = t71 * qJD(4) - t158 * t88;
t33 = -t72 * qJD(4) + t155 * t88;
t32 = t137 * t228 - t138 * t55;
t30 = -t137 * t226 - t138 * t52;
t25 = t145 * t72 - t149 * t71;
t19 = t145 * t45 + t38;
t14 = t145 * t33 + t149 * t34;
t13 = t145 * t34 - t149 * t33;
t9 = -pkin(5) * t66 - pkin(9) * t67 + t27;
t7 = t157 * t9;
t5 = [t216, 0 (qJDD(2) * t159 - t156 * t161) * t148 (-qJDD(2) * t156 - t159 * t161) * t148, t124 * t152 - t47 * t91 + t48 * t92 - t75 * t87 - t76 * t88 - g(3), 0, 0, 0, 0, 0, -t91 * t205 + qJD(4) * t33 + qJDD(4) * t71 + (-t158 * t87 + t91 * t210) * qJD(2), t91 * t206 - qJD(4) * t34 - qJDD(4) * t72 + (qJD(4) * t158 * t91 + t155 * t87) * qJD(2), t13 * t98 + t14 * t96 + t25 * t67 + t26 * t66, -t13 * t17 + t14 * t18 - t25 * t3 + t26 * t4 + t27 * t91 + t64 * t87 - g(3), 0, 0, 0, 0, 0 (-t182 * qJD(6) - t14 * t154 + t157 * t87) * t95 + t252 * t65 + t13 * t78 + t25 * t24 -(t252 * qJD(6) + t14 * t157 + t154 * t87) * t95 - t182 * t65 + t13 * t80 + t25 * t23; 0, qJDD(2), t121 - g(2) * (t201 - t227) + t165, -g(1) * (t147 * t221 - t151 * t159) - g(2) * (-t147 * t159 - t151 * t221) - t216 * t225, -g(2) * t118 + t75 * t86 - t76 * t89 + (g(2) * t227 + t48 * t146 + t47 * t150 + t165) * pkin(2), qJDD(2) * t143 + 0.2e1 * t155 * t196, 0.2e1 * t155 * t205 - 0.2e1 * t215 * t208, t116, t117, 0, t155 * t166 - t158 * t163, t155 * t163 + t158 * t166, t100 * t17 - t103 * t4 - t105 * t3 - t18 * t97 + t237 * t96 + t238 * t98 + t61 * t67 + t62 * t66 + t175, t4 * t62 - t3 * t61 + t27 * t179 - g(1) * (pkin(2) * t176 + t135 * t56 - t255 * t55) - g(2) * (-pkin(2) * t227 + t135 * t53 - t255 * t52 + t118) - g(3) * (pkin(2) * t223 - t135 * t91 + t255 * t92) + t251 * t64 + t237 * t18 - t238 * t17, -t80 * t219 + (t157 * t23 - t80 * t211) * t105 -(-t154 * t80 - t157 * t78) * t100 + (-t233 - t157 * t24 + (t154 * t78 - t157 * t80) * qJD(6)) * t105, t239 + t248, t164 + t185, t103 * t65 + t95 * t97, t7 * t103 + t61 * t24 - t183 * t97 + t238 * t78 + (t236 * t95 + (-t16 * t103 - t62 * t95 + t234) * qJD(6) + t247) * t157 + t249 * t154, t61 * t23 - t6 * t97 + t238 * t80 + (-(-qJD(6) * t16 + t9) * t103 - qJD(6) * t234 + (qJD(6) * t62 - t236) * t95 - t247) * t154 + t249 * t157; 0, 0, 0, 0, t173 + t124, 0, 0, 0, 0, 0, t117, -t116, -t100 * t96 + t103 * t67 + t105 * t66 + t97 * t98, -t100 * t18 - t103 * t3 + t105 * t4 - t17 * t97 + t173, 0, 0, 0, 0, 0, t164 - t185, t239 - t248; 0, 0, 0, 0, 0, -t155 * t161 * t158, t215 * t161, t206, t205, qJDD(4), t170 * t155 - t186 * t224 + t112 - t246, g(3) * t72 + (t186 * t148 - t124) * t155 + t170 * t158 (t18 - t19) * t98 + (t17 - t20) * t96 + (t145 * t66 - t149 * t67) * pkin(4), t17 * t19 - t18 * t20 + (t4 * t145 + t3 * t149 - t64 * t212 - g(1) * (t147 * t224 + t155 * t55) - g(2) * (-t151 * t224 + t155 * t52) - t246) * pkin(4), t80 * t192 + t233 (t23 - t243) * t157 + (-t24 - t241) * t154, t95 * t192 + t231 - t240, t178 + t242, -t95 * t98, t132 * t24 + t168 * t154 - t250 * t157 + t183 * t98 - t19 * t78, t132 * t23 + t250 * t154 + t168 * t157 - t19 * t80 + t6 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96 ^ 2 - t98 ^ 2, t17 * t98 - t18 * t96 + t174 + t27, 0, 0, 0, 0, 0, t178 - t242, -t95 ^ 2 * t157 - t231 - t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t78, -t78 ^ 2 + t80 ^ 2, t23 + t243, -t24 + t241, t65, -t154 * t2 + t7 - t15 * t80 - g(1) * (-t154 * t32 - t157 * t56) - g(2) * (-t154 * t30 - t157 * t53) - g(3) * (-t154 * t69 + t84) + t256 * t6, -t157 * t2 - t154 * t9 + t15 * t78 - g(1) * (t154 * t56 - t157 * t32) - g(2) * (t154 * t53 - t157 * t30) - g(3) * (-t157 * t69 - t230) - t256 * t183;];
tau_reg  = t5;
