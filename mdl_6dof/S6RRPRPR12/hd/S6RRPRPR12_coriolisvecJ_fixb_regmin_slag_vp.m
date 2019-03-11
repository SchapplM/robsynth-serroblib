% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:21:47
% EndTime: 2019-03-09 11:22:01
% DurationCPUTime: 4.74s
% Computational Cost: add. (7139->442), mult. (18472->626), div. (0->0), fcn. (14017->10), ass. (0->234)
t216 = sin(qJ(4));
t219 = cos(qJ(4));
t213 = sin(pkin(6));
t217 = sin(qJ(2));
t300 = qJD(1) * t217;
t275 = t213 * t300;
t257 = t219 * t275;
t221 = -pkin(2) - pkin(9);
t309 = qJ(5) - t221;
t267 = t309 * t219;
t190 = pkin(2) * t275;
t220 = cos(qJ(2));
t250 = pkin(9) * t217 - qJ(3) * t220;
t301 = qJD(1) * t213;
t130 = t250 * t301 + t190;
t282 = t220 * t301;
t214 = cos(pkin(6));
t292 = t214 * qJD(1);
t288 = pkin(1) * t292;
t151 = pkin(8) * t282 + t217 * t288;
t132 = pkin(3) * t282 + t151;
t306 = t219 * t130 + t216 * t132;
t346 = -qJ(5) * t257 - qJD(4) * t267 - t216 * qJD(5) - t306;
t264 = -t216 * t130 + t219 * t132;
t296 = qJD(4) * t216;
t311 = t216 * t217;
t345 = -t219 * qJD(5) + t309 * t296 - (pkin(4) * t220 - qJ(5) * t311) * t301 - t264;
t194 = qJD(2) + t292;
t140 = t216 * t194 + t219 * t282;
t142 = t219 * t194 - t216 * t282;
t212 = sin(pkin(11));
t320 = cos(pkin(11));
t263 = -t320 * t140 - t212 * t142;
t338 = qJD(6) - t263;
t178 = qJD(4) + t275;
t215 = sin(qJ(6));
t218 = cos(qJ(6));
t235 = -t212 * t140 + t320 * t142;
t68 = -t218 * t178 + t215 * t235;
t344 = t338 * t68;
t269 = qJD(4) * t320;
t314 = t212 * t216;
t305 = -t212 * t296 + t219 * t269 + t320 * t257 - t275 * t314;
t234 = t212 * t219 + t320 * t216;
t126 = t234 * t275;
t295 = qJD(4) * t219;
t154 = -t212 * t295 - t216 * t269;
t304 = -t154 + t126;
t268 = t218 * t338;
t289 = qJD(1) * qJD(2);
t272 = t213 * t289;
t256 = t217 * t272;
t96 = -qJD(4) * t140 + t216 * t256;
t299 = qJD(2) * t217;
t231 = t219 * t299 + t220 * t296;
t97 = t194 * t295 - t231 * t301;
t63 = t212 * t96 + t320 * t97;
t323 = t215 * t63;
t343 = -t268 * t338 - t323;
t150 = pkin(8) * t275 - t220 * t288;
t339 = qJD(3) + t150;
t342 = t339 + pkin(4) * t295 + (pkin(4) * t219 + pkin(3)) * t275;
t328 = -t346 * t212 + t345 * t320;
t326 = t345 * t212 + t346 * t320;
t313 = t213 * t217;
t196 = pkin(8) * t313;
t284 = -pkin(1) * t220 - pkin(2);
t110 = pkin(3) * t313 + t196 + (-pkin(9) + t284) * t214;
t271 = -qJ(3) * t217 - pkin(1);
t128 = (t221 * t220 + t271) * t213;
t307 = t216 * t110 + t219 * t128;
t332 = pkin(1) * t217;
t199 = t214 * t332;
t312 = t213 * t220;
t340 = pkin(8) * t312 + t199;
t291 = pkin(3) * t275 + t339;
t203 = t212 * pkin(4) + pkin(10);
t180 = t220 * t272;
t287 = pkin(1) * qJD(2) * t214;
t259 = qJD(1) * t287;
t143 = pkin(8) * t180 + t217 * t259;
t111 = pkin(3) * t180 + t143;
t106 = qJD(1) * t128;
t83 = t221 * t194 + t291;
t57 = t219 * t106 + t216 * t83;
t173 = pkin(2) * t256;
t297 = qJD(3) * t217;
t223 = (qJD(2) * t250 - t297) * t213;
t90 = qJD(1) * t223 + t173;
t226 = -qJD(4) * t57 + t219 * t111 - t216 * t90;
t18 = pkin(4) * t180 - t96 * qJ(5) - t142 * qJD(5) + t226;
t237 = -t106 * t296 + t216 * t111 + t219 * t90 + t83 * t295;
t22 = -t97 * qJ(5) - t140 * qJD(5) + t237;
t5 = t320 * t18 - t212 * t22;
t3 = -pkin(5) * t180 - t5;
t337 = (t142 * pkin(4) + pkin(5) * t235 - pkin(10) * t263 + qJD(6) * t203) * t338 + t3;
t167 = t309 * t216;
t115 = -t320 * t167 - t212 * t267;
t159 = t320 * t219 - t314;
t310 = t216 * pkin(4) + qJ(3);
t99 = pkin(5) * t234 - t159 * pkin(10) + t310;
t336 = (-t305 * pkin(5) - t304 * pkin(10) + qJD(6) * t115 - t342) * t338 - t99 * t63;
t64 = -t212 * t97 + t320 * t96;
t303 = pkin(8) * t256 - t220 * t259;
t117 = -t194 * qJD(3) + t303;
t91 = -pkin(3) * t256 - t117;
t65 = t97 * pkin(4) + t91;
t16 = t63 * pkin(5) - t64 * pkin(10) + t65;
t56 = -t216 * t106 + t219 * t83;
t49 = -t142 * qJ(5) + t56;
t46 = t178 * pkin(4) + t49;
t50 = -t140 * qJ(5) + t57;
t47 = t320 * t50;
t20 = t212 * t46 + t47;
t14 = t178 * pkin(10) + t20;
t181 = t194 * qJ(3);
t100 = t181 + t132;
t73 = t140 * pkin(4) + qJD(5) + t100;
t38 = -pkin(5) * t263 - pkin(10) * t235 + t73;
t249 = t215 * t14 - t218 * t38;
t6 = t212 * t18 + t320 * t22;
t4 = pkin(10) * t180 + t6;
t1 = -qJD(6) * t249 + t215 * t16 + t218 * t4;
t324 = t212 * t50;
t19 = t320 * t46 - t324;
t334 = t5 * t159 - t304 * t19 + t305 * t20 + t234 * t6;
t333 = pkin(3) + pkin(8);
t331 = t68 * t235;
t70 = t215 * t178 + t218 * t235;
t330 = t70 * t235;
t155 = t214 * t216 + t219 * t312;
t281 = t213 * t299;
t122 = -qJD(4) * t155 + t216 * t281;
t156 = t214 * t219 - t216 * t312;
t185 = pkin(2) * t281;
t103 = t185 + t223;
t133 = (t333 * t312 + t199) * qJD(2);
t224 = -t307 * qJD(4) - t216 * t103 + t219 * t133;
t298 = qJD(2) * t220;
t280 = t213 * t298;
t27 = pkin(4) * t280 - t122 * qJ(5) - t156 * qJD(5) + t224;
t123 = -t231 * t213 + t214 * t295;
t236 = t219 * t103 + t110 * t295 - t128 * t296 + t216 * t133;
t31 = -t123 * qJ(5) - t155 * qJD(5) + t236;
t12 = t212 * t27 + t320 * t31;
t265 = t219 * t110 - t216 * t128;
t54 = pkin(4) * t313 - t156 * qJ(5) + t265;
t59 = -t155 * qJ(5) + t307;
t33 = t212 * t54 + t320 * t59;
t327 = pkin(5) * t282 - t328;
t325 = t115 * t63;
t293 = qJD(6) * t218;
t294 = qJD(6) * t215;
t36 = t178 * t293 + t215 * t180 + t218 * t64 - t235 * t294;
t322 = t36 * t215;
t321 = t63 * t234;
t319 = t140 * t178;
t318 = t142 * t178;
t317 = t159 * t218;
t316 = t178 * t221;
t209 = t213 ^ 2;
t315 = t209 * qJD(1) ^ 2;
t210 = t217 ^ 2;
t302 = -t220 ^ 2 + t210;
t286 = t178 * t217 * t219;
t285 = t220 * t315;
t144 = -t214 * qJ(3) - t340;
t277 = t178 * t295;
t273 = t209 * t289;
t270 = -t218 * t180 + t215 * t64;
t262 = t194 + t292;
t261 = -qJD(6) * t234 - t194;
t260 = 0.2e1 * t273;
t258 = t217 * t285;
t127 = pkin(3) * t312 - t144;
t254 = -0.2e1 * pkin(1) * t273;
t252 = t151 * t194 - t143;
t8 = t218 * t14 + t215 * t38;
t29 = pkin(10) * t313 + t33;
t240 = t155 * pkin(4) + t127;
t92 = t320 * t155 + t212 * t156;
t93 = -t212 * t155 + t320 * t156;
t44 = t92 * pkin(5) - t93 * pkin(10) + t240;
t248 = t215 * t44 + t218 * t29;
t247 = -t215 * t29 + t218 * t44;
t152 = t340 * qJD(2);
t246 = t143 * t214 + t152 * t194;
t245 = t218 * t63 + (t215 * t263 - t294) * t338;
t189 = t220 * t287;
t244 = -pkin(8) * t281 + t189;
t243 = -t194 * t282 + t180;
t242 = -t215 * t93 + t218 * t313;
t75 = t215 * t313 + t218 * t93;
t11 = -t212 * t31 + t320 * t27;
t32 = -t212 * t59 + t320 * t54;
t241 = t100 * t217 + t221 * t298;
t239 = t178 * t216;
t145 = (-pkin(2) * t220 + t271) * t213;
t88 = t215 * t126 - t218 * t282;
t233 = t154 * t215 + t159 * t293 - t88;
t89 = t218 * t126 + t215 * t282;
t232 = t154 * t218 - t159 * t294 - t89;
t206 = t214 * qJD(3);
t109 = -t333 * t281 + t189 + t206;
t230 = (-qJ(3) * t298 - t297) * t213;
t13 = -t178 * pkin(5) - t19;
t24 = t320 * t49 - t324;
t228 = -t203 * t63 + (t13 + t24) * t338;
t227 = t123 * pkin(4) + t109;
t2 = -qJD(6) * t8 + t218 * t16 - t215 * t4;
t225 = -t325 + t3 * t159 + (pkin(10) * t282 - qJD(6) * t99 - t326) * t338;
t204 = -t320 * pkin(4) - pkin(5);
t166 = t219 * t180;
t149 = -qJ(3) * t282 + t190;
t146 = t214 * t284 + t196;
t138 = -t206 - t244;
t137 = qJD(1) * t145;
t134 = t185 + t230;
t129 = -t181 - t151;
t124 = -t194 * pkin(2) + t339;
t114 = -t212 * t167 + t320 * t267;
t113 = qJD(1) * t230 + t173;
t108 = t137 * t275;
t72 = t320 * t122 - t212 * t123;
t71 = t212 * t122 + t320 * t123;
t42 = qJD(6) * t75 + t215 * t72 - t218 * t280;
t41 = qJD(6) * t242 + t215 * t280 + t218 * t72;
t37 = qJD(6) * t70 + t270;
t28 = -pkin(5) * t313 - t32;
t25 = t71 * pkin(5) - t72 * pkin(10) + t227;
t23 = t212 * t49 + t47;
t10 = pkin(10) * t280 + t12;
t9 = -pkin(5) * t280 - t11;
t7 = [0, 0, 0, t217 * t220 * t260, -t302 * t260, t262 * t280, -t262 * t281, 0, t217 * t254 - t246, -t244 * t194 + t303 * t214 + t220 * t254 (-t117 * t220 + t143 * t217 + (t124 * t220 + t129 * t217) * qJD(2) + (-t138 * t220 + t152 * t217 + (t144 * t217 + t146 * t220) * qJD(2)) * qJD(1)) * t213 (-t137 * t299 + t113 * t220 + (t134 * t220 - t145 * t299) * qJD(1)) * t213 + t246, -t117 * t214 - t138 * t194 + (-t137 * t298 - t113 * t217 + (-t134 * t217 - t145 * t298) * qJD(1)) * t213, t113 * t145 + t117 * t144 + t124 * t152 + t129 * t138 + t137 * t134 + t143 * t146, t142 * t122 + t96 * t156, -t122 * t140 - t142 * t123 - t96 * t155 - t156 * t97, t122 * t178 + (t217 * t96 + (qJD(1) * t156 + t142) * t298) * t213, -t123 * t178 + (-t217 * t97 + (-qJD(1) * t155 - t140) * t298) * t213 (t178 * t213 + t209 * t300) * t298, t224 * t178 + t109 * t140 + t127 * t97 + t91 * t155 + t100 * t123 + (t226 * t217 + (qJD(1) * t265 + t56) * t298) * t213, -t236 * t178 + t109 * t142 + t127 * t96 + t91 * t156 + t100 * t122 + (-t237 * t217 + (-t307 * qJD(1) - t57) * t298) * t213, -t11 * t235 + t12 * t263 - t19 * t72 - t20 * t71 - t32 * t64 - t33 * t63 - t5 * t93 - t6 * t92, t19 * t11 + t20 * t12 + t227 * t73 + t240 * t65 + t5 * t32 + t6 * t33, t36 * t75 + t70 * t41, t242 * t36 - t75 * t37 - t41 * t68 - t70 * t42, t338 * t41 + t36 * t92 + t75 * t63 + t70 * t71, t242 * t63 - t338 * t42 - t37 * t92 - t68 * t71, t338 * t71 + t63 * t92 (-qJD(6) * t248 - t215 * t10 + t218 * t25) * t338 + t247 * t63 + t2 * t92 - t249 * t71 + t9 * t68 + t28 * t37 - t3 * t242 + t13 * t42 -(qJD(6) * t247 + t218 * t10 + t215 * t25) * t338 - t248 * t63 - t1 * t92 - t8 * t71 + t9 * t70 + t28 * t36 + t3 * t75 + t13 * t41; 0, 0, 0, -t258, t302 * t315, t243 (-qJD(2) + t194) * t275, 0, t315 * t332 + t252, pkin(1) * t285 - t150 * t194 + t303 ((-qJ(3) * qJD(2) - t129 - t151) * t217 + (-pkin(2) * qJD(2) - t124 + t339) * t220) * t301, -t149 * t282 + t108 - t252, t339 * t194 + (t137 * t220 + t149 * t217) * t301 - t117, -t143 * pkin(2) - t117 * qJ(3) - t124 * t151 - t129 * t339 - t137 * t149, -t142 * t239 + t96 * t219 (-t97 - t318) * t219 + (-t96 + t319) * t216, -t178 * t296 + t166 + (-t142 * t220 - t178 * t311) * t301, -t277 + (-t286 + (-qJD(2) * t216 + t140) * t220) * t301, -t178 * t282, qJ(3) * t97 + t91 * t216 - t264 * t178 + t291 * t140 + (t100 * t219 - t216 * t316) * qJD(4) + (t219 * t241 - t56 * t220) * t301, qJ(3) * t96 + t91 * t219 + t306 * t178 + t291 * t142 + (-t100 * t216 - t219 * t316) * qJD(4) + (-t216 * t241 + t57 * t220) * t301, t114 * t64 - t235 * t328 + t263 * t326 - t325 - t334, -t5 * t114 + t6 * t115 + t328 * t19 + t326 * t20 + t65 * t310 + t342 * t73, t232 * t70 + t36 * t317, t89 * t68 + t70 * t88 + (-t215 * t70 - t218 * t68) * t154 + (-t322 - t218 * t37 + (t215 * t68 - t218 * t70) * qJD(6)) * t159, t232 * t338 + t234 * t36 + t305 * t70 + t63 * t317, -t159 * t323 - t233 * t338 - t234 * t37 - t305 * t68, t305 * t338 + t321, t114 * t37 + t233 * t13 + t2 * t234 + t225 * t215 - t218 * t336 - t249 * t305 + t327 * t68, -t1 * t234 + t114 * t36 + t232 * t13 + t215 * t336 + t225 * t218 - t305 * t8 + t327 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, t258, -t194 ^ 2 - t210 * t315, t129 * t194 + t108 + t143, 0, 0, 0, 0, 0, -t194 * t140 - t178 * t239 + t166, -t277 - t194 * t142 + (-t216 * t298 - t286) * t301, -t159 * t64 + t235 * t304 + t263 * t305 - t321, -t73 * t194 + t334, 0, 0, 0, 0, 0, -t215 * t321 - t159 * t37 + t304 * t68 + (-t305 * t215 + t218 * t261) * t338, -t218 * t321 - t159 * t36 + t304 * t70 + (-t261 * t215 - t218 * t305) * t338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142 * t140, -t140 ^ 2 + t142 ^ 2, t96 + t319, t318 - t97, t180, -t100 * t142 + t57 * t178 + t226, t100 * t140 + t178 * t56 - t237 (-t212 * t63 - t320 * t64) * pkin(4) + (t19 - t24) * t263 + (t20 - t23) * t235, t19 * t23 - t20 * t24 + (-t142 * t73 + t212 * t6 + t320 * t5) * pkin(4), t268 * t70 + t322 (t36 - t344) * t218 + (-t338 * t70 - t37) * t215, -t330 - t343, t245 + t331, -t338 * t235, t204 * t37 + t228 * t215 - t218 * t337 - t23 * t68 + t235 * t249, t204 * t36 + t215 * t337 + t228 * t218 - t23 * t70 + t8 * t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235 ^ 2 - t263 ^ 2, t19 * t235 - t20 * t263 + t65, 0, 0, 0, 0, 0, t245 - t331, -t330 + t343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t68, -t68 ^ 2 + t70 ^ 2, t36 + t344, -t270 + (-qJD(6) + t338) * t70, t63, -t13 * t70 + t338 * t8 + t2, t13 * t68 - t249 * t338 - t1;];
tauc_reg  = t7;
