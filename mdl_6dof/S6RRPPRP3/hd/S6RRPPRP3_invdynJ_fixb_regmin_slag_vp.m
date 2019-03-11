% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:37
% EndTime: 2019-03-09 08:35:45
% DurationCPUTime: 3.36s
% Computational Cost: add. (2559->439), mult. (5025->514), div. (0->0), fcn. (2755->6), ass. (0->233)
t277 = pkin(2) + pkin(3);
t206 = t277 * qJD(2);
t143 = sin(qJ(2));
t231 = qJD(1) * t143;
t108 = pkin(7) * t231;
t68 = -qJ(4) * t231 + t108;
t210 = qJD(3) + t68;
t46 = -t206 + t210;
t142 = sin(qJ(5));
t145 = cos(qJ(5));
t223 = t145 * qJD(2);
t146 = cos(qJ(2));
t230 = qJD(1) * t146;
t65 = t142 * t230 - t223;
t88 = qJD(5) + t231;
t269 = t65 * t88;
t226 = qJD(5) * t142;
t205 = t146 * t226;
t219 = t146 * qJDD(1);
t27 = -t142 * qJDD(2) - t145 * t219 - qJD(5) * t223 + (t143 * t223 + t205) * qJD(1);
t286 = t27 - t269;
t229 = qJD(2) * t142;
t66 = t145 * t230 + t229;
t268 = t66 * t88;
t222 = qJD(1) * qJD(2);
t202 = t143 * t222;
t28 = qJD(5) * t66 - t145 * qJDD(2) + (-t202 + t219) * t142;
t165 = t28 - t268;
t290 = t143 * t277;
t201 = t146 * t222;
t220 = t143 * qJDD(1);
t289 = t201 + t220;
t203 = t277 * qJDD(2);
t182 = pkin(4) * t143 + pkin(8) * t146;
t57 = -qJD(1) * pkin(1) - pkin(2) * t230 - qJ(3) * t231;
t43 = pkin(3) * t230 + qJD(4) - t57;
t33 = qJD(1) * t182 + t43;
t132 = -pkin(8) - t277;
t41 = qJD(2) * t132 + t210;
t13 = t142 * t33 + t145 * t41;
t104 = pkin(7) * t220;
t87 = pkin(7) * t201;
t204 = qJDD(3) + t104 + t87;
t221 = qJD(1) * qJD(4);
t151 = -qJ(4) * t289 - t143 * t221 + t204;
t20 = qJDD(2) * t132 + t151;
t64 = qJDD(5) + t289;
t166 = pkin(4) * t146 + t132 * t143;
t158 = t166 * qJD(2);
t114 = t143 * qJD(3);
t137 = qJDD(1) * pkin(1);
t184 = pkin(2) * t219 + qJ(3) * t289 + qJD(1) * t114 + t137;
t162 = pkin(3) * t219 + qJDD(4) + t184;
t11 = qJD(1) * t158 + qJDD(1) * t182 + t162;
t8 = t145 * t11;
t1 = pkin(5) * t64 - qJ(6) * t27 - qJD(5) * t13 + qJD(6) * t66 - t142 * t20 + t8;
t10 = qJ(6) * t65 + t13;
t225 = qJD(5) * t145;
t169 = -t11 * t142 - t145 * t20 - t225 * t33 + t226 * t41;
t2 = qJ(6) * t28 + qJD(6) * t65 - t169;
t126 = g(3) * t146;
t147 = cos(qJ(1));
t246 = t143 * t147;
t144 = sin(qJ(1));
t248 = t143 * t144;
t216 = g(1) * t246 + g(2) * t248 - t126;
t12 = -t142 * t41 + t145 * t33;
t9 = qJ(6) * t66 + t12;
t5 = pkin(5) * t88 + t9;
t288 = -(t10 * t88 + t1) * t142 - (t5 * t88 - t2) * t145 - t216;
t287 = t88 ^ 2;
t133 = qJD(2) * qJ(3);
t109 = pkin(7) * t230;
t70 = -qJ(4) * t230 + t109;
t55 = -t133 - t70;
t127 = g(2) * t144;
t285 = g(1) * t147 + t127;
t241 = t145 * t147;
t58 = t142 * t248 - t241;
t244 = t144 * t145;
t60 = -t142 * t246 - t244;
t284 = -g(1) * t60 + g(2) * t58;
t283 = -pkin(5) * t28 + qJDD(6);
t136 = qJDD(2) * pkin(4);
t228 = qJD(2) * t143;
t215 = pkin(7) * t228;
t167 = qJD(4) * t146 + t215;
t105 = pkin(7) * t219;
t130 = qJDD(2) * qJ(3);
t131 = qJD(2) * qJD(3);
t208 = t105 + t130 + t131;
t82 = qJ(4) * t202;
t187 = -t82 - t208;
t23 = qJ(4) * t219 + qJD(1) * t167 + t187;
t21 = t136 - t23;
t270 = g(3) * t143;
t280 = qJD(5) * t132 * t88 + t146 * t285 - t21 + t270;
t256 = pkin(7) * qJDD(2);
t116 = t143 * qJ(3);
t123 = t146 * pkin(2);
t235 = t123 + t116;
t74 = -pkin(1) - t235;
t279 = (qJD(1) * t74 + t57) * qJD(2) - t256;
t278 = t66 ^ 2;
t276 = -t9 + t5;
t272 = pkin(5) * t142;
t128 = g(1) * t144;
t271 = g(2) * t147;
t122 = t146 * pkin(3);
t267 = pkin(7) - qJ(4);
t103 = pkin(5) * t145 + pkin(4);
t239 = qJ(6) - t132;
t191 = qJD(5) * t239;
t100 = qJ(3) * t230;
t36 = qJD(1) * t166 + t100;
t264 = t142 * t36 + t145 * t70;
t266 = -qJD(6) * t145 - t264 + (qJ(6) * t231 + t191) * t142;
t247 = t143 * t145;
t172 = pkin(5) * t146 - qJ(6) * t247;
t35 = t145 * t36;
t265 = -qJD(1) * t172 + t145 * t191 - t35 + (qJD(6) + t70) * t142;
t207 = t122 + t235;
t62 = pkin(1) + t207;
t40 = t182 + t62;
t77 = t267 * t143;
t67 = t145 * t77;
t263 = t142 * t40 + t67;
t262 = t142 * t27;
t261 = t142 * t64;
t260 = t145 * t64;
t259 = t145 * t66;
t258 = t146 * t65;
t257 = t146 * t66;
t255 = qJ(3) * t146;
t254 = qJD(2) * t55;
t253 = qJD(2) * t65;
t252 = qJD(2) * t66;
t251 = qJDD(2) * pkin(2);
t250 = t142 * t146;
t249 = t142 * t147;
t150 = qJD(1) ^ 2;
t245 = t143 * t150;
t243 = t144 * t146;
t242 = t145 * t146;
t240 = t146 * t147;
t237 = -qJD(4) - t43;
t227 = qJD(2) * t146;
t236 = qJ(3) * t227 + t114;
t134 = t143 ^ 2;
t135 = t146 ^ 2;
t233 = t134 - t135;
t232 = t134 + t135;
t224 = qJD(6) * t146;
t140 = -qJ(6) - pkin(8);
t218 = -t140 + t277;
t32 = t158 + t236;
t53 = -qJD(4) * t143 + t227 * t267;
t217 = t142 * t32 + t145 * t53 + t225 * t40;
t214 = t142 * t143 * t88;
t213 = t88 * t247;
t212 = qJ(6) * t242;
t211 = t146 * t245;
t209 = t105 + 0.2e1 * t130 + 0.2e1 * t131;
t200 = -pkin(1) - t116;
t199 = t128 - t271;
t198 = qJD(1) * t62 + t43;
t196 = -qJD(5) * t33 - t20;
t195 = t237 * t143;
t194 = -qJD(2) * pkin(2) + qJD(3);
t124 = t147 * pkin(7);
t193 = -qJ(4) * t147 + t124;
t192 = g(1) * t218;
t190 = pkin(1) * t147 + pkin(2) * t240 + pkin(7) * t144 + qJ(3) * t246;
t189 = -t104 + t216;
t188 = 0.2e1 * t201;
t186 = t143 * t206;
t185 = -t225 * t41 + t8;
t183 = pkin(3) * t240 + t190;
t149 = qJD(2) ^ 2;
t181 = pkin(7) * t149 + t271;
t179 = t10 * t145 - t142 * t5;
t73 = t108 + t194;
t75 = t109 + t133;
t178 = t143 * t75 - t146 * t73;
t177 = -qJDD(3) + t189;
t174 = t200 - t123;
t173 = t87 - t177;
t171 = -t225 * t88 - t261;
t170 = t226 * t88 - t260;
t168 = -0.2e1 * pkin(1) * t222 - t256;
t49 = qJD(2) * pkin(4) - t55;
t161 = -t181 + 0.2e1 * t137;
t160 = -qJ(4) * qJDD(1) - t285;
t157 = -t132 * t64 - t49 * t88;
t17 = -qJD(1) * t186 + t162;
t42 = -t186 + t236;
t156 = -qJD(1) * t42 - qJDD(1) * t62 - t17 + t271;
t155 = -qJ(4) * t220 + t173;
t30 = pkin(2) * t202 - t184;
t54 = pkin(2) * t228 - t236;
t154 = -qJD(1) * t54 - qJDD(1) * t74 - t181 - t30;
t44 = -pkin(7) * t202 + t208;
t50 = t204 - t251;
t153 = -qJD(2) * t178 + t50 * t143 + t44 * t146;
t152 = qJD(5) * t179 + t1 * t145 + t142 * t2 - t271;
t141 = qJ(3) + pkin(4);
t117 = t146 * qJ(4);
t101 = qJ(4) * t228;
t93 = g(1) * t243;
t92 = g(1) * t248;
t86 = qJ(3) * t240;
t84 = qJ(3) * t243;
t80 = -t134 * t150 - t149;
t78 = pkin(7) * t146 - t117;
t76 = qJDD(2) + t211;
t72 = t239 * t145;
t71 = t239 * t142;
t69 = pkin(2) * t231 - t100;
t63 = t65 ^ 2;
t61 = -t142 * t144 + t143 * t241;
t59 = -t143 * t244 - t249;
t52 = -t101 + t167;
t51 = -t231 * t277 + t100;
t39 = t145 * t40;
t31 = -pkin(5) * t65 + qJD(6) + t49;
t26 = t145 * t32;
t22 = -t203 + t151;
t19 = qJ(6) * t250 + t263;
t16 = pkin(5) * t143 - t142 * t77 + t212 + t39;
t6 = t21 + t283;
t4 = qJD(5) * t212 + (-qJ(6) * t228 - qJD(5) * t77 + t224) * t142 + t217;
t3 = t145 * t224 - t142 * t53 + t26 + t172 * qJD(2) + (-t67 + (-qJ(6) * t146 - t40) * t142) * qJD(5);
t7 = [qJDD(1), t199, t285, qJDD(1) * t134 + t143 * t188, 0.2e1 * t143 * t219 - 0.2e1 * t222 * t233, qJDD(2) * t143 + t146 * t149, qJDD(2) * t146 - t143 * t149, 0, t143 * t168 + t146 * t161 + t93, -t143 * t161 + t146 * t168 - t92, t143 * t279 + t146 * t154 + t93, pkin(7) * qJDD(1) * t232 + t153 - t285, t143 * t154 - t146 * t279 + t92, pkin(7) * t153 - g(1) * t124 - g(2) * t190 - t128 * t174 + t30 * t74 + t57 * t54, qJDD(2) * t78 + t92 + (t146 * t198 - t52) * qJD(2) - t156 * t143, qJDD(2) * t77 - t93 + (t143 * t198 + t53) * qJD(2) + t156 * t146 (-qJD(2) * t46 - qJDD(1) * t78 + t23 + (-qJD(2) * t77 + t52) * qJD(1)) * t146 + (-t254 - qJDD(1) * t77 - t22 + (qJD(2) * t78 - t53) * qJD(1)) * t143 + t285, t22 * t77 + t46 * t53 - t23 * t78 + t55 * t52 + t17 * t62 + t43 * t42 - g(1) * t193 - g(2) * t183 + (-g(1) * (t174 - t122) + g(2) * qJ(4)) * t144, -t66 * t205 + (-t146 * t27 - t228 * t66) * t145 (t142 * t66 + t145 * t65) * t228 + (t262 - t145 * t28 + (t142 * t65 - t259) * qJD(5)) * t146 (t223 * t88 + t27) * t143 + (t170 - t252) * t146 (-t229 * t88 + t28) * t143 + (-t171 + t253) * t146, t143 * t64 + t227 * t88 (-t225 * t77 + t26) * t88 + t39 * t64 + t185 * t143 + t52 * t65 - t78 * t28 - g(1) * t59 - g(2) * t61 + (qJD(2) * t12 - t225 * t49) * t146 + ((-qJD(5) * t40 - t53) * t88 - t77 * t64 - t21 * t146 + (qJD(2) * t49 + t196) * t143) * t142 -(-t226 * t77 + t217) * t88 - t263 * t64 + t52 * t66 + t78 * t27 - g(1) * t58 - g(2) * t60 + (t223 * t49 + t169) * t143 + (-qJD(2) * t13 - t21 * t145 + t226 * t49) * t146, -t16 * t27 + t19 * t28 + t3 * t66 + t4 * t65 + t93 + (-t10 * t142 - t145 * t5) * t228 + t152 * t146, t2 * t19 + t10 * t4 + t1 * t16 + t5 * t3 - t6 * t117 + t31 * (t228 * t272 + t101 - t215) - g(1) * (-pkin(5) * t249 + t193) - g(2) * (t103 * t246 + t183) + (t6 * (pkin(7) - t272) + t31 * (-pkin(5) * t225 - qJD(4)) + t140 * t271) * t146 + (-g(1) * (-t103 * t143 + t200) - g(2) * (-qJ(4) - t272) + t146 * t192) * t144; 0, 0, 0, -t211, t233 * t150, t220, t219, qJDD(2), pkin(1) * t245 + t189, t270 - t105 + (pkin(1) * t150 + t285) * t146, 0.2e1 * t251 + (-t143 * t57 + t146 * t69) * qJD(1) + t177 (-pkin(2) * t143 + t255) * qJDD(1) + ((t75 - t133) * t143 + (t194 - t73) * t146) * qJD(1) (qJD(1) * t69 - g(3)) * t143 + (qJD(1) * t57 - t285) * t146 + t209, t44 * qJ(3) + t75 * qJD(3) - t50 * pkin(2) - t57 * t69 - g(1) * (-pkin(2) * t246 + t86) - g(2) * (-pkin(2) * t248 + t84) - g(3) * t235 + t178 * qJD(1) * pkin(7), qJD(2) * t68 + t82 + (-g(3) + (-pkin(7) * qJD(2) - t51) * qJD(1)) * t143 + (qJD(1) * t237 + t160) * t146 + t209, -qJD(2) * t70 - 0.2e1 * t203 + ((-qJ(4) * qJD(2) + t51) * t146 + t195) * qJD(1) + t155 (-t255 + t290) * qJDD(1), -g(1) * t86 - g(2) * t84 - g(3) * t207 - t23 * qJ(3) - t210 * t55 - t22 * t277 + t285 * t290 - t43 * t51 - t46 * t70, t259 * t88 - t262 (-t27 - t269) * t145 + (-t28 - t268) * t142 (-t213 + t257) * qJD(1) + t171 (t214 - t258) * qJD(1) + t170, -t88 * t230, -t12 * t230 - t141 * t28 - t35 * t88 - t210 * t65 + (t70 * t88 + t157) * t142 - t280 * t145, t13 * t230 + t141 * t27 + t142 * t280 + t145 * t157 - t210 * t66 + t264 * t88, t265 * t66 + t266 * t65 - t27 * t71 - t28 * t72 - t288, -t2 * t72 + t1 * t71 + t6 * (qJ(3) + t103) - g(1) * (t103 * t240 + t86) - g(2) * (t103 * t243 + t84) - g(3) * (-t140 * t146 + t207) + t265 * t5 + (-pkin(5) * t226 + t210) * t31 + t266 * t10 + (-qJD(1) * t272 * t31 - g(3) * t103 + t127 * t218 + t147 * t192) * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t220, t80, -qJD(2) * t75 + t231 * t57 + t173 - t251, t80, t76, -t220, t254 - t203 + (-qJ(4) * t227 + t195) * qJD(1) + t155, 0, 0, 0, 0, 0, -t145 * t287 + t253 - t261, t142 * t287 + t252 - t260, t142 * t286 + t145 * t165, -qJD(2) * t31 + t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188 + t220, 0.2e1 * t202 - t219, -t232 * t150 (-t146 * t55 + (t46 - t206) * t143) * qJD(1) + t162 + t199, 0, 0, 0, 0, 0 (-t214 - t258) * qJD(1) - t170 (-t213 - t257) * qJD(1) + t171, t142 * t165 - t145 * t286, t128 + (t143 * t179 + t146 * t31) * qJD(1) + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t65, -t63 + t278, t286, t165, t64, t13 * t88 + t49 * t66 + (t196 - t126) * t142 + t185 + t284, g(1) * t61 - g(2) * t59 - g(3) * t242 + t12 * t88 - t49 * t65 + t169, -pkin(5) * t27 + t276 * t65, t276 * t10 + (-g(3) * t250 + t31 * t66 + t1 + t284) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63 - t278, -t10 * t65 - t5 * t66 + t136 + (-pkin(7) * t222 - g(3)) * t143 + (t160 - t221) * t146 - t187 + t283;];
tau_reg  = t7;
