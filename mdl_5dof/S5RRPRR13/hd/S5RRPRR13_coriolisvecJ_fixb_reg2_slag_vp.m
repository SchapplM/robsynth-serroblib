% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR13_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:33:52
% EndTime: 2019-12-31 20:34:05
% DurationCPUTime: 4.75s
% Computational Cost: add. (7713->409), mult. (19306->595), div. (0->0), fcn. (13975->8), ass. (0->203)
t212 = sin(pkin(9));
t216 = sin(qJ(2));
t264 = qJD(1) * t216;
t247 = t212 * t264;
t213 = cos(pkin(9));
t257 = t213 * qJD(2);
t163 = t247 - t257;
t246 = t213 * t264;
t258 = t212 * qJD(2);
t165 = t246 + t258;
t215 = sin(qJ(4));
t217 = cos(qJ(4));
t103 = t215 * t163 - t217 * t165;
t104 = t217 * t163 + t215 * t165;
t214 = sin(qJ(5));
t291 = cos(qJ(5));
t44 = t291 * t103 + t214 * t104;
t311 = t44 ^ 2;
t48 = t214 * t103 - t291 * t104;
t310 = t48 ^ 2;
t218 = cos(qJ(2));
t256 = t218 * qJD(1);
t199 = -qJD(4) + t256;
t192 = -qJD(5) + t199;
t309 = t44 * t192;
t308 = t48 * t192;
t171 = t217 * t212 + t215 * t213;
t153 = t171 * qJD(4);
t226 = t171 * t218;
t269 = qJD(1) * t226 - t153;
t273 = t217 * t213;
t170 = t215 * t212 - t273;
t260 = qJD(4) * t217;
t261 = qJD(4) * t215;
t295 = -t212 * t261 + t213 * t260;
t301 = t170 * t256 + t295;
t233 = pkin(2) * t216 - qJ(3) * t218;
t173 = t233 * qJD(1);
t154 = t212 * t173;
t275 = t213 * t216;
t276 = t212 * t218;
t227 = -pkin(6) * t275 - pkin(7) * t276;
t110 = qJD(1) * t227 + t154;
t289 = pkin(7) + qJ(3);
t183 = t289 * t212;
t184 = t289 * t213;
t120 = -t215 * t183 + t217 * t184;
t122 = pkin(6) * t247 + t213 * t173;
t274 = t213 * t218;
t230 = pkin(3) * t216 - pkin(7) * t274;
t93 = qJD(1) * t230 + t122;
t286 = t171 * qJD(3) + t120 * qJD(4) - t215 * t110 + t217 * t93;
t285 = -qJD(3) * t273 + t217 * t110 + t183 * t260 + (qJD(3) * t212 + qJD(4) * t184 + t93) * t215;
t307 = t103 ^ 2;
t306 = t104 ^ 2;
t305 = -pkin(4) * t264 - t301 * pkin(8) - t286;
t304 = -t269 * pkin(8) + t285;
t290 = t48 * t44;
t303 = t103 * t199;
t302 = t104 * t199;
t300 = -t310 + t311;
t244 = qJD(5) * t291;
t259 = qJD(5) * t214;
t255 = qJD(1) * qJD(2);
t243 = t218 * t255;
t236 = t212 * t243;
t237 = t213 * t243;
t62 = t215 * (qJD(4) * t165 + t236) + t163 * t260 - t217 * t237;
t63 = -t163 * t261 + t165 * t260 + t215 * t237 + t217 * t236;
t19 = -t103 * t259 + t104 * t244 + t214 * t63 + t291 * t62;
t299 = -t19 + t308;
t180 = -t218 * pkin(2) - t216 * qJ(3) - pkin(1);
t156 = t180 * qJD(1);
t204 = pkin(6) * t256;
t186 = qJD(2) * qJ(3) + t204;
t112 = t213 * t156 - t212 * t186;
t69 = -pkin(3) * t256 - t165 * pkin(7) + t112;
t113 = t212 * t156 + t213 * t186;
t73 = -t163 * pkin(7) + t113;
t34 = t215 * t69 + t217 * t73;
t225 = t230 * qJD(2);
t149 = qJD(2) * t233 - t216 * qJD(3);
t133 = t149 * qJD(1);
t203 = pkin(6) * t264;
t178 = (qJD(3) - t203) * qJD(2);
t91 = t213 * t133 - t212 * t178;
t66 = qJD(1) * t225 + t91;
t92 = t212 * t133 + t213 * t178;
t70 = -pkin(7) * t236 + t92;
t15 = -qJD(4) * t34 - t215 * t70 + t217 * t66;
t202 = t216 * t255;
t10 = pkin(4) * t202 + t62 * pkin(8) + t15;
t14 = t215 * t66 + t217 * t70 + t69 * t260 - t261 * t73;
t11 = -t63 * pkin(8) + t14;
t33 = -t215 * t73 + t217 * t69;
t28 = pkin(8) * t103 + t33;
t26 = -t199 * pkin(4) + t28;
t29 = -t104 * pkin(8) + t34;
t223 = -t214 * t10 - t291 * t11 - t26 * t244 + t29 * t259;
t282 = qJD(2) * pkin(2);
t240 = qJD(3) - t282;
t179 = t203 + t240;
t121 = t163 * pkin(3) + t179;
t59 = t104 * pkin(4) + t121;
t298 = -t48 * t59 + t223;
t296 = -0.2e1 * t255;
t162 = t213 * t180;
t111 = -pkin(7) * t275 + t162 + (-pkin(6) * t212 - pkin(3)) * t218;
t197 = pkin(6) * t274;
t129 = t212 * t180 + t197;
t277 = t212 * t216;
t118 = -pkin(7) * t277 + t129;
t53 = t215 * t111 + t217 * t118;
t210 = t216 ^ 2;
t211 = t218 ^ 2;
t294 = qJD(1) * (t210 - 0.2e1 * t211);
t252 = t291 * t29;
t6 = t214 * t26 + t252;
t2 = -qJD(5) * t6 + t291 * t10 - t214 * t11;
t293 = t44 * t59 + t2;
t20 = -qJD(5) * t44 - t214 * t62 + t291 * t63;
t292 = -t20 + t309;
t119 = -t217 * t183 - t215 * t184;
t86 = -t171 * pkin(8) + t119;
t87 = -t170 * pkin(8) + t120;
t37 = -t214 * t87 + t291 * t86;
t288 = t37 * qJD(5) + t305 * t214 - t304 * t291;
t38 = t214 * t86 + t291 * t87;
t287 = -t38 * qJD(5) + t304 * t214 + t305 * t291;
t284 = t170 * t244 + t171 * t259 - t269 * t214 - t301 * t291;
t109 = -t214 * t170 + t291 * t171;
t283 = t109 * qJD(5) + t301 * t214 - t269 * t291;
t281 = t214 * t29;
t280 = t103 * t104;
t279 = t163 * t213;
t220 = qJD(1) ^ 2;
t278 = t211 * t220;
t272 = t218 * t220;
t219 = qJD(2) ^ 2;
t271 = t219 * t216;
t270 = t219 * t218;
t263 = qJD(2) * t216;
t253 = pkin(6) * t263;
t116 = t213 * t149 + t212 * t253;
t198 = pkin(6) * t243;
t148 = pkin(3) * t236 + t198;
t251 = t212 * t256;
t157 = pkin(3) * t251 + t204;
t262 = qJD(2) * t218;
t205 = pkin(6) * t262;
t250 = t218 * t258;
t158 = pkin(3) * t250 + t205;
t174 = pkin(3) * t277 + t216 * pkin(6);
t265 = t210 - t211;
t254 = pkin(6) * t276;
t201 = -t213 * pkin(3) - pkin(2);
t245 = -t269 * pkin(4) - t157;
t241 = pkin(1) * t296;
t52 = t217 * t111 - t215 * t118;
t239 = t163 + t257;
t238 = -t165 + t258;
t42 = t63 * pkin(4) + t148;
t235 = t218 * t202;
t234 = -t179 + t240;
t232 = qJD(1) * t239;
t231 = qJD(1) * t238;
t143 = t170 * t216;
t36 = -t218 * pkin(4) + t143 * pkin(8) + t52;
t142 = t171 * t216;
t39 = -t142 * pkin(8) + t53;
t22 = -t214 * t39 + t291 * t36;
t23 = t214 * t36 + t291 * t39;
t82 = -t214 * t142 - t291 * t143;
t83 = t225 + t116;
t135 = t212 * t149;
t95 = qJD(2) * t227 + t135;
t24 = t111 * t260 - t118 * t261 + t215 * t83 + t217 * t95;
t228 = t218 * t231;
t25 = -t53 * qJD(4) - t215 * t95 + t217 * t83;
t209 = t213 ^ 2;
t208 = t212 ^ 2;
t195 = t216 * t272;
t187 = -0.2e1 * t235;
t134 = t170 * pkin(4) + t201;
t128 = t162 - t254;
t123 = -pkin(6) * t246 + t154;
t117 = -t213 * t253 + t135;
t115 = t142 * pkin(4) + t174;
t108 = t291 * t170 + t214 * t171;
t89 = qJD(2) * t226 + t295 * t216;
t88 = -t217 * t218 * t257 + t153 * t216 + t215 * t250;
t81 = t291 * t142 - t214 * t143;
t64 = t89 * pkin(4) + t158;
t31 = t82 * qJD(5) - t214 * t88 + t291 * t89;
t30 = t142 * t244 - t143 * t259 + t214 * t89 + t291 * t88;
t21 = -t89 * pkin(8) + t24;
t18 = pkin(4) * t263 + t88 * pkin(8) + t25;
t9 = t291 * t28 - t281;
t8 = -t214 * t28 - t252;
t5 = t291 * t26 - t281;
t4 = -t23 * qJD(5) + t291 * t18 - t214 * t21;
t3 = t22 * qJD(5) + t214 * t18 + t291 * t21;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t235, t265 * t296, t270, t187, -t271, 0, -pkin(6) * t270 + t216 * t241, pkin(6) * t271 + t218 * t241, 0, 0, (t165 * t213 + t209 * t264) * t262, (-t279 + (-t165 - 0.2e1 * t246) * t212) * t262, (t165 * t216 + t213 * t294) * qJD(2), (t163 * t212 + t208 * t264) * t262, (-t163 * t216 - t212 * t294) * qJD(2), t187, (-qJD(1) * t116 - t91) * t218 + ((pkin(6) * t163 + t179 * t212) * t218 + (t112 + (t128 + 0.2e1 * t254) * qJD(1)) * t216) * qJD(2), (qJD(1) * t117 + t92) * t218 + ((pkin(6) * t165 + t179 * t213) * t218 + (-t113 + (-t129 + 0.2e1 * t197) * qJD(1)) * t216) * qJD(2), -t116 * t165 - t117 * t163 + (-t212 * t92 - t213 * t91) * t216 + (-t112 * t213 - t113 * t212 + (-t128 * t213 - t129 * t212) * qJD(1)) * t262, t112 * t116 + t113 * t117 + t91 * t128 + t92 * t129 + (t179 + t203) * t205, t103 * t88 + t62 * t143, t103 * t89 + t88 * t104 + t62 * t142 + t143 * t63, t88 * t199 + t62 * t218 + (-qJD(1) * t143 - t103) * t263, t104 * t89 + t63 * t142, t89 * t199 + t63 * t218 + (-qJD(1) * t142 - t104) * t263, (-t199 - t256) * t263, t158 * t104 + t121 * t89 + t148 * t142 - t15 * t218 + t174 * t63 - t25 * t199 + (qJD(1) * t52 + t33) * t263, -t158 * t103 - t121 * t88 + t14 * t218 - t148 * t143 - t174 * t62 + t24 * t199 + (-qJD(1) * t53 - t34) * t263, t103 * t25 - t24 * t104 - t14 * t142 + t15 * t143 + t33 * t88 - t34 * t89 + t52 * t62 - t53 * t63, t121 * t158 + t14 * t53 + t148 * t174 + t15 * t52 + t34 * t24 + t33 * t25, -t19 * t82 + t30 * t44, t19 * t81 - t82 * t20 - t30 * t48 + t31 * t44, t19 * t218 + t30 * t192 + (qJD(1) * t82 - t44) * t263, t20 * t81 - t31 * t48, t31 * t192 + t20 * t218 + (-qJD(1) * t81 + t48) * t263, (-t192 - t256) * t263, t115 * t20 - t4 * t192 - t2 * t218 + t59 * t31 + t42 * t81 - t64 * t48 + (qJD(1) * t22 + t5) * t263, -t223 * t218 - t115 * t19 + t3 * t192 - t59 * t30 + t42 * t82 - t64 * t44 + (-qJD(1) * t23 - t6) * t263, t22 * t19 - t2 * t82 - t23 * t20 + t223 * t81 + t3 * t48 + t5 * t30 - t6 * t31 + t4 * t44, t42 * t115 + t2 * t22 - t223 * t23 + t6 * t3 + t5 * t4 + t59 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t195, t265 * t220, 0, t195, 0, 0, t220 * pkin(1) * t216, pkin(1) * t272, 0, 0, t213 * t228, (t279 + t165 * t212 + (-t208 + t209) * qJD(2)) * t256, t213 * t278 + t216 * t231, -t239 * t251, -t212 * t278 + t216 * t232, t195, ((-qJ(3) * t258 - t112) * t216 + (-pkin(6) * t239 + t212 * t234 + t122) * t218) * qJD(1), ((-qJ(3) * t257 + t113) * t216 + (pkin(6) * t238 + t213 * t234 - t123) * t218) * qJD(1), t122 * t165 + t123 * t163 + (-qJD(3) * t163 + t112 * t256 + t92) * t213 + (qJD(3) * t165 + t113 * t256 - t91) * t212, -t112 * t122 - t113 * t123 + (-t112 * t212 + t113 * t213) * qJD(3) + (-t91 * t212 + t92 * t213) * qJ(3) + (-t179 - t282) * t204, -t103 * t301 - t62 * t171, -t103 * t269 - t104 * t301 + t62 * t170 - t171 * t63, -t301 * t199 + (qJD(2) * t171 + t103) * t264, -t269 * t104 + t63 * t170, -t269 * t199 + (-qJD(2) * t170 + t104) * t264, t199 * t264, -t157 * t104 + t148 * t170 + t201 * t63 + t286 * t199 - t269 * t121 + (qJD(2) * t119 - t33) * t264, t157 * t103 + t148 * t171 - t201 * t62 - t285 * t199 + t301 * t121 + (-qJD(2) * t120 + t34) * t264, -t103 * t286 + t285 * t104 + t119 * t62 - t120 * t63 - t14 * t170 - t15 * t171 + t269 * t34 - t301 * t33, t15 * t119 + t14 * t120 - t121 * t157 + t148 * t201 - t285 * t34 - t286 * t33, -t19 * t109 + t284 * t44, t19 * t108 - t109 * t20 + t283 * t44 - t284 * t48, t284 * t192 + (qJD(2) * t109 + t44) * t264, t20 * t108 - t283 * t48, t283 * t192 + (-qJD(2) * t108 - t48) * t264, t192 * t264, t42 * t108 + t134 * t20 + t283 * t59 - t245 * t48 - t287 * t192 + (qJD(2) * t37 - t5) * t264, t42 * t109 - t134 * t19 - t284 * t59 - t245 * t44 + t288 * t192 + (-qJD(2) * t38 + t6) * t264, t108 * t223 - t2 * t109 + t37 * t19 - t38 * t20 - t283 * t6 + t284 * t5 + t287 * t44 + t288 * t48, t42 * t134 + t2 * t37 - t223 * t38 + t245 * t59 + t287 * t5 + t288 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t218 * t232, -t163 ^ 2 - t165 ^ 2, t112 * t165 + t113 * t163 + t198, 0, 0, 0, 0, 0, 0, t63 + t303, -t62 + t302, -t306 - t307, -t33 * t103 + t104 * t34 + t148, 0, 0, 0, 0, 0, 0, t20 + t309, -t19 - t308, -t310 - t311, -t5 * t44 - t6 * t48 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t280, -t306 + t307, -t62 - t302, t280, -t63 + t303, t202, t103 * t121 - t34 * t199 + t15, t121 * t104 - t33 * t199 - t14, 0, 0, t290, t300, t299, -t290, t292, t202, t8 * t192 + (-t103 * t48 + t192 * t259 + t291 * t202) * pkin(4) + t293, -t9 * t192 + (-t103 * t44 + t192 * t244 - t202 * t214) * pkin(4) + t298, -t6 * t44 - t9 * t48 + t5 * t48 - t8 * t44 + (t291 * t19 - t20 * t214 + (-t214 * t44 + t291 * t48) * qJD(5)) * pkin(4), -t5 * t8 - t6 * t9 + (t291 * t2 - t223 * t214 + t103 * t59 + (-t214 * t5 + t291 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, t300, t299, -t290, t292, t202, -t6 * t192 + t293, -t5 * t192 + t298, 0, 0;];
tauc_reg = t1;
