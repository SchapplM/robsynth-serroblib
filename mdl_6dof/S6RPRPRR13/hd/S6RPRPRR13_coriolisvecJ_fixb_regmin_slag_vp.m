% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR13_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:25:21
% EndTime: 2019-03-09 04:25:35
% DurationCPUTime: 5.26s
% Computational Cost: add. (8658->450), mult. (28915->641), div. (0->0), fcn. (24963->12), ass. (0->224)
t186 = sin(pkin(12));
t188 = sin(pkin(6));
t189 = cos(pkin(12));
t192 = sin(qJ(3));
t187 = sin(pkin(7));
t301 = cos(pkin(6));
t262 = t301 * t187;
t300 = cos(pkin(7));
t263 = t192 * t300;
t311 = cos(qJ(3));
t136 = (t311 * t186 + t189 * t263) * t188 + t192 * t262;
t319 = qJD(1) * t136;
t327 = qJD(5) + t319;
t246 = t300 * t311;
t293 = t188 * t189;
t214 = t246 * t293;
t210 = qJD(1) * t214;
t230 = t311 * t262;
t287 = qJD(1) * t188;
t270 = t192 * t287;
t123 = -qJD(1) * t230 + t186 * t270 - t210;
t241 = t301 * t300;
t271 = t189 * t287;
t149 = -qJD(1) * t241 + t187 * t271 - qJD(3);
t191 = sin(qJ(5));
t194 = cos(qJ(5));
t84 = -t194 * t123 - t149 * t191;
t330 = t327 * t84;
t295 = t186 * t188;
t329 = t192 * t295 - t230;
t202 = t188 * (t186 * t246 + t189 * t192);
t144 = qJD(1) * t202;
t285 = qJD(3) * t192;
t269 = t187 * t285;
t328 = t144 - t269;
t83 = qJD(6) + t84;
t154 = t187 * t293 - t241;
t264 = t188 * t300;
t201 = (t189 * t264 + t262) * pkin(9);
t274 = pkin(1) * t301;
t289 = qJ(2) * t293 + t186 * t274;
t133 = t201 + t289;
t181 = t189 * t274;
t199 = t301 * pkin(2) + (-t300 * pkin(9) - qJ(2)) * t295;
t137 = t181 + t199;
t296 = t186 * t187;
t146 = (-pkin(2) * t189 - pkin(9) * t296 - pkin(1)) * t188;
t254 = t311 * t293;
t171 = qJD(2) * t254;
t228 = qJD(3) * t246;
t249 = t186 * t264;
t229 = qJD(2) * t249;
t268 = qJD(3) * t311;
t250 = t187 * t268;
t197 = -t192 * (qJD(3) * t133 + t229) + t137 * t228 + t146 * t250 + t171;
t53 = t154 * qJD(4) - t197;
t326 = t329 * qJD(3);
t190 = sin(qJ(6));
t193 = cos(qJ(6));
t86 = t123 * t191 - t149 * t194;
t71 = t190 * t86 - t193 * t327;
t325 = t327 * t71;
t324 = t188 ^ 2 * (t186 ^ 2 + t189 ^ 2);
t323 = t194 * t327;
t260 = qJD(1) * t301;
t251 = pkin(1) * t260;
t152 = qJ(2) * t271 + t186 * t251;
t110 = qJD(1) * t201 + t152;
t177 = t189 * t251;
t120 = qJD(1) * t199 + t177;
t142 = qJD(1) * t146 + qJD(2);
t273 = t187 * t311;
t63 = t192 * t110 - t120 * t246 - t142 * t273;
t322 = qJD(4) + t63;
t207 = -t191 * t300 - t194 * t273;
t272 = t186 * t287;
t253 = t187 * t272;
t321 = -qJD(5) * t207 + t328 * t191 + t194 * t253;
t156 = -t191 * t273 + t194 * t300;
t320 = qJD(5) * t156 - t191 * t253 + t328 * t194;
t145 = (-t192 * t249 + t254) * qJD(1);
t227 = t250 - t145;
t318 = -qJD(5) + t327;
t108 = t326 * qJD(1) - qJD(3) * t210;
t290 = pkin(4) * t319 + t322;
t312 = pkin(3) + pkin(10);
t29 = t312 * t149 + t290;
t81 = -t120 * t187 + t300 * t142;
t221 = -qJ(4) * t319 + t81;
t34 = t312 * t123 + t221;
t14 = t191 * t29 + t194 * t34;
t212 = qJD(1) * t229;
t99 = t120 * t263;
t52 = t189 * qJD(2) * t270 + qJD(3) * t99 + t110 * t268 + t142 * t269 + t311 * t212;
t32 = -pkin(4) * t108 + t52;
t126 = t136 * qJD(3);
t109 = qJD(1) * t126;
t286 = qJD(2) * t188;
t169 = t286 * t296;
t164 = qJD(1) * t169;
t216 = qJ(4) * t108 - qJD(4) * t319 + t164;
t40 = t312 * t109 + t216;
t265 = t191 * t40 - t194 * t32;
t4 = pkin(5) * t108 + t14 * qJD(5) + t265;
t317 = (pkin(5) * t86 + t83 * pkin(11)) * t83 + t4;
t148 = t149 * qJ(4);
t294 = t187 * t192;
t64 = t311 * t110 + t142 * t294 + t99;
t49 = -pkin(4) * t123 + t64;
t35 = -t148 + t49;
t316 = -t108 * t312 - t327 * t35;
t281 = qJD(5) * t194;
t283 = qJD(5) * t191;
t69 = t191 * t109 + t123 * t281 + t149 * t283;
t73 = t190 * t327 + t193 * t86;
t23 = t73 * qJD(6) + t193 * t108 + t190 * t69;
t315 = t300 * t109 - t123 * t253 - t149 * t328;
t314 = -t300 * t108 + t149 * t227 - t253 * t319;
t244 = qJD(1) * t171 - t110 * t285 + t120 * t228 + t142 * t250 - t192 * t212;
t39 = qJD(4) * t149 - t244;
t27 = -pkin(4) * t109 - t39;
t98 = t194 * t109;
t70 = t86 * qJD(5) - t98;
t10 = pkin(5) * t70 - pkin(11) * t69 + t27;
t12 = pkin(11) * t327 + t14;
t20 = pkin(5) * t84 - pkin(11) * t86 + t35;
t240 = t12 * t190 - t193 * t20;
t219 = t191 * t32 + t194 * t40 + t29 * t281 - t34 * t283;
t3 = -pkin(11) * t108 + t219;
t1 = -t240 * qJD(6) + t10 * t190 + t193 * t3;
t313 = t319 ^ 2;
t310 = t71 * t83;
t309 = t73 * t83;
t198 = t192 * t133 - t137 * t246 - t146 * t273;
t41 = t136 * pkin(4) + t312 * t154 + t198;
t135 = -t214 + t329;
t87 = -t137 * t187 + t300 * t146;
t220 = -qJ(4) * t136 + t87;
t46 = t312 * t135 + t220;
t235 = t191 * t41 + t194 * t46;
t299 = qJ(4) * t123;
t68 = t312 * t319 + t299;
t308 = t191 * t49 + t194 * t68;
t54 = pkin(3) * t123 + t221;
t307 = t319 * t54;
t278 = qJD(6) * t193;
t279 = qJD(6) * t190;
t22 = -t190 * t108 + t193 * t69 + t278 * t327 - t86 * t279;
t306 = t190 * t22;
t305 = t190 * t70;
t304 = t193 * t70;
t259 = t193 * t83;
t303 = t312 * t83;
t245 = pkin(5) * t194 + pkin(11) * t191;
t302 = t245 * qJD(5) - (-pkin(4) - t245) * t319 + t322;
t298 = t319 * t123;
t297 = t319 * t191;
t284 = qJD(5) * t190;
t282 = qJD(5) * t193;
t280 = qJD(5) * t312;
t118 = t311 * t133;
t261 = t300 * t137;
t258 = t83 * t327;
t257 = t327 * t73;
t173 = pkin(5) * t191 - pkin(11) * t194 + qJ(4);
t256 = -pkin(11) * t123 - qJD(6) * t173 + t308;
t255 = t191 * t327;
t61 = t154 * qJ(4) - t146 * t294 - t192 * t261 - t118;
t196 = qJD(1) ^ 2;
t248 = t188 * t196 * t301;
t6 = t12 * t193 + t190 * t20;
t58 = qJD(2) * t202 + (t118 + (t146 * t187 + t261) * t192) * qJD(3);
t239 = t149 * t58 + t154 * t52;
t17 = pkin(11) * t136 + t235;
t232 = t135 * t194 + t154 * t191;
t51 = -pkin(4) * t135 - t61;
t89 = t135 * t191 - t154 * t194;
t21 = -pkin(5) * t232 - pkin(11) * t89 + t51;
t238 = t17 * t193 + t190 * t21;
t237 = -t17 * t190 + t193 * t21;
t13 = -t191 * t34 + t194 * t29;
t236 = -t191 * t46 + t194 * t41;
t125 = -qJD(3) * t214 + t326;
t43 = -t125 * pkin(4) + t58;
t215 = qJ(4) * t125 - qJD(4) * t136 + t169;
t50 = t312 * t126 + t215;
t234 = -t191 * t50 + t194 * t43;
t233 = t136 * t193 - t190 * t89;
t77 = t136 * t190 + t193 * t89;
t231 = (-qJ(2) * t272 + t177) * t186 - t152 * t189;
t226 = -t83 * t278 - t305;
t225 = t83 * t279 - t304;
t224 = -t194 * t108 - t255 * t327;
t223 = -t156 * t190 + t193 * t294;
t222 = t156 * t193 + t190 * t294;
t218 = t191 * t43 + t194 * t50 + t41 * t281 - t46 * t283;
t217 = -0.2e1 * t260 * t286;
t211 = -t149 * t64 - t52;
t209 = t191 * t108 - t323 * t327;
t11 = -pkin(5) * t327 - t13;
t208 = -pkin(11) * t70 + (t11 + t13) * t83;
t2 = -t6 * qJD(6) + t193 * t10 - t190 * t3;
t200 = -t123 * t149 - t108;
t33 = -t126 * pkin(4) - t53;
t82 = t108 * t136;
t80 = pkin(3) * t319 + t299;
t79 = -t123 * t190 + t193 * t297;
t78 = t193 * t123 + t190 * t297;
t75 = t89 * qJD(5) - t126 * t194;
t74 = t232 * qJD(5) + t126 * t191;
t65 = pkin(3) * t126 + t215;
t62 = t154 * pkin(3) + t198;
t60 = pkin(3) * t135 + t220;
t57 = pkin(3) * t109 + t216;
t56 = t148 - t64;
t55 = pkin(3) * t149 + t322;
t26 = t233 * qJD(6) - t125 * t190 + t193 * t74;
t25 = t77 * qJD(6) + t125 * t193 + t190 * t74;
t18 = pkin(5) * t123 + t191 * t68 - t194 * t49;
t16 = -pkin(5) * t136 - t236;
t15 = t75 * pkin(5) - t74 * pkin(11) + t33;
t8 = pkin(5) * t125 + t235 * qJD(5) - t234;
t7 = -pkin(11) * t125 + t218;
t5 = [0, 0, 0, t186 * t217, t189 * t217, 0.2e1 * qJD(2) * qJD(1) * t324 ((t189 * t289 + (qJ(2) * t295 - t181) * t186) * qJD(1) - t231) * t286, -t125 * t319 - t82, t108 * t135 - t109 * t136 + t123 * t125 - t126 * t319, t108 * t154 + t125 * t149, t109 * t154 + t126 * t149, 0, t109 * t87 + t126 * t81 + (qJD(1) * t135 + t123) * t169 + t239, -t87 * t108 - t81 * t125 + t136 * t164 + t149 * t197 + t154 * t244 + t169 * t319, -t108 * t62 + t109 * t61 + t123 * t53 - t125 * t55 + t126 * t56 + t135 * t39 + t136 * t52 + t319 * t58, -t109 * t60 - t123 * t65 - t126 * t54 - t135 * t57 - t239, t108 * t60 + t125 * t54 - t136 * t57 + t149 * t53 + t154 * t39 - t319 * t65, t39 * t61 + t52 * t62 + t53 * t56 + t54 * t65 + t55 * t58 + t57 * t60, t69 * t89 + t74 * t86, t232 * t69 - t70 * t89 - t74 * t84 - t75 * t86, -t108 * t89 - t125 * t86 + t136 * t69 + t327 * t74, -t108 * t232 + t125 * t84 - t136 * t70 - t327 * t75, -t125 * t327 - t82, t234 * t327 - t236 * t108 - t265 * t136 - t13 * t125 + t33 * t84 + t51 * t70 - t27 * t232 + t35 * t75 + (-t136 * t14 - t235 * t327) * qJD(5), t235 * t108 + t14 * t125 - t219 * t136 - t218 * t327 + t27 * t89 + t33 * t86 + t35 * t74 + t51 * t69, t22 * t77 + t26 * t73, t22 * t233 - t23 * t77 - t25 * t73 - t26 * t71, -t22 * t232 + t26 * t83 + t70 * t77 + t73 * t75, t23 * t232 + t233 * t70 - t25 * t83 - t71 * t75, -t232 * t70 + t75 * t83 (-qJD(6) * t238 + t15 * t193 - t190 * t7) * t83 + t237 * t70 - t2 * t232 - t240 * t75 + t8 * t71 + t16 * t23 - t4 * t233 + t11 * t25 -(qJD(6) * t237 + t15 * t190 + t193 * t7) * t83 - t238 * t70 + t1 * t232 - t6 * t75 + t8 * t73 + t16 * t22 + t4 * t77 + t11 * t26; 0, 0, 0, t186 * t248, t189 * t248, -t196 * t324, t231 * t287, 0, 0, 0, 0, 0, t315, t314, t145 * t123 - t144 * t319 + (t311 * t108 - t109 * t192 + (-t311 * t123 + t192 * t319) * qJD(3)) * t187, -t315, -t314, t57 * t300 - t55 * t144 + t56 * t145 + (-t54 * t272 - t311 * t52 - t192 * t39 + (t192 * t55 - t311 * t56) * qJD(3)) * t187, 0, 0, 0, 0, 0, -t207 * t108 - t145 * t84 + (t192 * t70 + t268 * t84) * t187 - t320 * t327, t156 * t108 - t145 * t86 + (t192 * t69 + t268 * t86) * t187 + t321 * t327, 0, 0, 0, 0, 0, -t207 * t23 + t223 * t70 + (-qJD(6) * t222 + t190 * t321 + t193 * t227) * t83 + t320 * t71, -t207 * t22 - t222 * t70 + (-qJD(6) * t223 - t190 * t227 + t193 * t321) * t83 + t320 * t73; 0, 0, 0, 0, 0, 0, 0, t298, -t123 ^ 2 + t313, t200 -(qJD(3) + t149) * t319, 0, -t319 * t81 + t211, t123 * t81 + t149 * t63 - t244, pkin(3) * t108 - qJ(4) * t109 + (-t56 - t64) * t319 + (t55 - t322) * t123, t123 * t80 - t211 + t307, -t123 * t54 - t149 * t322 + t319 * t80 - t39, -pkin(3) * t52 - qJ(4) * t39 - t322 * t56 - t54 * t80 - t55 * t64, t194 * t69 - t255 * t86 (-t327 * t86 - t70) * t194 + (-t69 + t330) * t191, t123 * t86 + t224, -t123 * t84 + t209, t327 * t123, qJ(4) * t70 + t13 * t123 + t290 * t84 + (t27 + (t68 + t280) * t327) * t191 + (-t327 * t49 - t316) * t194, qJ(4) * t69 - t14 * t123 + t27 * t194 + t290 * t86 + (t194 * t280 + t308) * t327 + t316 * t191, t193 * t194 * t22 + (-t191 * t282 - t194 * t279 - t79) * t73, t71 * t79 + t73 * t78 + (t190 * t73 + t193 * t71) * t283 + (-t306 - t193 * t23 + (t190 * t71 - t193 * t73) * qJD(6)) * t194, -t79 * t83 + (-t282 * t83 + t22) * t191 + (t257 - t225) * t194, t78 * t83 + (t284 * t83 - t23) * t191 + (t226 - t325) * t194, t191 * t70 + t83 * t323, t173 * t304 - t11 * t78 - t18 * t71 + (t190 * t256 + t193 * t302) * t83 + (-t11 * t284 + t2 - (qJD(5) * t71 + t226) * t312) * t191 + (t11 * t278 - t240 * t319 + t4 * t190 + t312 * t23 + (t190 * t303 - t240) * qJD(5)) * t194, -t173 * t305 - t11 * t79 - t18 * t73 + (-t190 * t302 + t193 * t256) * t83 + (-t11 * t282 - t1 - (qJD(5) * t73 + t225) * t312) * t191 + (-t11 * t279 - t6 * t319 + t4 * t193 + t312 * t22 + (t193 * t303 - t6) * qJD(5)) * t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, -t298, -t149 ^ 2 - t313, -t149 * t56 + t307 + t52, 0, 0, 0, 0, 0, t149 * t84 + t224, t149 * t86 + t209, 0, 0, 0, 0, 0, t149 * t259 + (-t190 * t258 - t23) * t194 + (t226 + t325) * t191, -t149 * t190 * t83 + (-t193 * t258 - t22) * t194 + (t257 + t225) * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86 * t84, -t84 ^ 2 + t86 ^ 2, t69 + t330, t318 * t86 + t98, -t108, t14 * t318 - t35 * t86 - t265, t13 * t327 + t35 * t84 - t219, t259 * t73 + t306 (t22 - t310) * t193 + (-t23 - t309) * t190, t259 * t83 - t73 * t86 + t305, -t190 * t83 ^ 2 + t71 * t86 + t304, -t83 * t86, -pkin(5) * t23 - t14 * t71 + t208 * t190 - t193 * t317 + t240 * t86, -pkin(5) * t22 - t14 * t73 + t190 * t317 + t208 * t193 + t6 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, t22 + t310, -t23 + t309, t70, -t11 * t73 + t6 * t83 + t2, t11 * t71 - t240 * t83 - t1;];
tauc_reg  = t5;
