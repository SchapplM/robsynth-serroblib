% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:33:05
% EndTime: 2019-03-09 16:33:16
% DurationCPUTime: 4.54s
% Computational Cost: add. (9795->470), mult. (23386->588), div. (0->0), fcn. (17480->14), ass. (0->266)
t234 = cos(qJ(2));
t371 = pkin(7) + pkin(8);
t182 = t371 * t234;
t175 = qJD(1) * t182;
t233 = cos(qJ(3));
t157 = t233 * t175;
t230 = sin(qJ(2));
t181 = t371 * t230;
t173 = qJD(1) * t181;
t229 = sin(qJ(3));
t311 = qJD(1) * t234;
t293 = t233 * t311;
t312 = qJD(1) * t230;
t296 = t229 * t312;
t150 = -t293 + t296;
t337 = qJ(4) * t150;
t107 = t173 * t229 - t157 + t337;
t152 = -t229 * t311 - t233 * t312;
t146 = t152 * qJ(4);
t153 = t229 * t175;
t316 = -t233 * t173 - t153;
t108 = t146 + t316;
t225 = sin(pkin(10));
t226 = cos(pkin(10));
t326 = t225 * t229;
t349 = pkin(2) * qJD(3);
t338 = -t107 * t225 - t108 * t226 + (t226 * t233 - t326) * t349;
t224 = qJ(2) + qJ(3);
t216 = sin(t224);
t217 = cos(t224);
t231 = sin(qJ(1));
t235 = cos(qJ(1));
t273 = g(1) * t235 + g(2) * t231;
t385 = -g(3) * t217 + t216 * t273;
t305 = qJD(1) * qJD(2);
t290 = t234 * t305;
t304 = t230 * qJDD(1);
t383 = t290 + t304;
t228 = sin(qJ(5));
t308 = qJD(5) * t228;
t284 = -t226 * t150 + t152 * t225;
t331 = t284 * t228;
t382 = t308 - t331;
t232 = cos(qJ(5));
t374 = qJD(5) - t284;
t279 = t374 * t232;
t302 = qJD(2) + qJD(3);
t303 = t234 * qJDD(1);
t100 = qJD(3) * t293 + t229 * t303 + t233 * t383 - t302 * t296;
t167 = t229 * t234 + t230 * t233;
t125 = t302 * t167;
t269 = t229 * t304 - t233 * t303;
t101 = qJD(1) * t125 + t269;
t62 = -t100 * t225 - t226 * t101;
t59 = qJDD(5) - t62;
t381 = -t228 * t59 - t374 * t279;
t219 = t234 * pkin(2);
t357 = pkin(1) + t219;
t262 = -t150 * t225 - t226 * t152;
t105 = t228 * t302 + t232 * t262;
t333 = t105 * t228;
t379 = t374 * t333;
t348 = qJD(2) * pkin(2);
t160 = -t173 + t348;
t283 = t233 * t160 - t153;
t98 = t146 + t283;
t214 = pkin(10) + t224;
t202 = sin(t214);
t203 = cos(t214);
t367 = pkin(5) * t232;
t209 = pkin(4) + t367;
t227 = -qJ(6) - pkin(9);
t260 = -t202 * t227 + t203 * t209;
t315 = -t229 * t181 + t233 * t182;
t378 = qJ(6) * t331 + t232 * qJD(6);
t377 = g(1) * t231 - g(2) * t235;
t212 = pkin(2) * t312;
t370 = pkin(3) * t152;
t75 = pkin(4) * t262 - pkin(9) * t284 - t370;
t72 = t212 + t75;
t376 = t228 * t72 - t232 * t338;
t319 = t232 * t235;
t324 = t228 * t231;
t135 = t203 * t324 + t319;
t321 = t231 * t232;
t323 = t228 * t235;
t137 = -t203 * t323 + t321;
t375 = -g(1) * t137 + g(2) * t135;
t126 = qJDD(2) * pkin(2) - t371 * t383;
t291 = t230 * t305;
t130 = t371 * (-t291 + t303);
t278 = qJD(3) * t160 + t130;
t310 = qJD(3) * t229;
t373 = t229 * t126 - t175 * t310 + t233 * t278;
t372 = t105 ^ 2;
t369 = pkin(3) * t216;
t368 = pkin(3) * t226;
t361 = g(3) * t202;
t360 = g(3) * t203;
t358 = g(3) * t228;
t261 = -t229 * t160 - t157;
t99 = -t261 - t337;
t346 = t226 * t99;
t89 = t302 * pkin(3) + t98;
t53 = t225 * t89 + t346;
t51 = t302 * pkin(9) + t53;
t180 = t357 * qJD(1);
t127 = pkin(3) * t150 + qJD(4) - t180;
t64 = -pkin(4) * t284 - pkin(9) * t262 + t127;
t30 = -t228 * t51 + t232 * t64;
t25 = -qJ(6) * t105 + t30;
t18 = pkin(5) * t374 + t25;
t356 = -t25 + t18;
t277 = t232 * t302;
t103 = t228 * t262 - t277;
t307 = qJD(5) * t232;
t263 = t226 * t100 - t225 * t101;
t301 = qJDD(2) + qJDD(3);
t244 = t228 * t263 - t232 * t301;
t46 = qJD(5) * t105 + t244;
t355 = -t103 * t307 - t228 * t46;
t123 = t233 * t126;
t322 = t229 * t130;
t41 = t301 * pkin(3) - t100 * qJ(4) + t261 * qJD(3) + t152 * qJD(4) + t123 - t322;
t47 = -qJ(4) * t101 - qJD(4) * t150 + t373;
t14 = t225 * t41 + t226 * t47;
t92 = t225 * t99;
t61 = t226 * t98 - t92;
t354 = t228 * t75 + t232 * t61;
t282 = -t233 * t181 - t182 * t229;
t114 = -qJ(4) * t167 + t282;
t166 = t229 * t230 - t233 * t234;
t115 = -qJ(4) * t166 + t315;
t82 = t114 * t225 + t115 * t226;
t79 = t232 * t82;
t121 = t226 * t166 + t167 * t225;
t122 = -t166 * t225 + t167 * t226;
t276 = pkin(3) * t166 - t357;
t80 = pkin(4) * t121 - pkin(9) * t122 + t276;
t352 = t228 * t80 + t79;
t210 = t233 * pkin(2) + pkin(3);
t325 = t226 * t229;
t145 = pkin(2) * t325 + t225 * t210;
t140 = pkin(9) + t145;
t318 = -qJ(6) - t140;
t281 = qJD(5) * t318;
t351 = t228 * t281 - t376 + t378;
t218 = t232 * qJ(6);
t271 = pkin(5) * t262 - t218 * t284;
t67 = t232 * t72;
t350 = t232 * t281 - t271 - t67 + (-qJD(6) - t338) * t228;
t347 = t18 * t232;
t45 = -qJD(5) * t277 - t228 * t301 - t232 * t263 + t262 * t308;
t345 = t228 * t45;
t124 = t302 * t166;
t87 = -t124 * t226 - t125 * t225;
t343 = t232 * t87;
t52 = t226 * t89 - t92;
t50 = -t302 * pkin(4) - t52;
t342 = t50 * t284;
t204 = pkin(3) * t225 + pkin(9);
t317 = -qJ(6) - t204;
t280 = qJD(5) * t317;
t341 = t228 * t280 - t354 + t378;
t74 = t232 * t75;
t340 = t232 * t280 - t271 - t74 + (-qJD(6) + t61) * t228;
t339 = t226 * t107 - t108 * t225 + (t225 * t233 + t325) * t349;
t336 = t103 * t262;
t335 = t103 * t284;
t334 = t105 * t262;
t332 = t374 * t262;
t330 = t122 * t228;
t329 = t122 * t232;
t328 = t152 * t150;
t207 = pkin(3) * t217;
t314 = t207 + t219;
t222 = t230 ^ 2;
t313 = -t234 ^ 2 + t222;
t309 = qJD(3) * t233;
t297 = qJD(2) * t371;
t174 = t230 * t297;
t176 = t234 * t297;
t247 = -t233 * t174 - t229 * t176 - t181 * t309 - t182 * t310;
t68 = -qJ(4) * t125 - qJD(4) * t166 + t247;
t240 = -qJD(3) * t315 + t229 * t174 - t233 * t176;
t69 = qJ(4) * t124 - qJD(4) * t167 + t240;
t35 = t225 * t69 + t226 * t68;
t213 = t230 * t348;
t289 = pkin(3) * t125 + t213;
t86 = -t124 * t225 + t226 * t125;
t39 = pkin(4) * t86 - pkin(9) * t87 + t289;
t300 = t228 * t39 + t232 * t35 + t80 * t307;
t48 = t50 * t308;
t295 = t122 * t307;
t294 = qJD(5) * t204 * t374;
t13 = -t225 * t47 + t226 * t41;
t11 = -t301 * pkin(4) - t13;
t292 = -t11 - t360;
t221 = -qJ(4) - t371;
t288 = pkin(5) * t228 - t221;
t34 = t225 * t68 - t226 * t69;
t60 = t225 * t98 + t346;
t12 = t301 * pkin(9) + t14;
t286 = -qJD(5) * t64 - t12;
t81 = -t226 * t114 + t115 * t225;
t144 = -pkin(2) * t326 + t210 * t226;
t139 = -pkin(4) - t144;
t275 = t382 * pkin(5);
t147 = pkin(2) * t291 - qJDD(1) * t357;
t241 = t101 * pkin(3) + qJDD(4) + t147;
t24 = -t62 * pkin(4) - pkin(9) * t263 + t241;
t23 = t232 * t24;
t274 = -t51 * t307 + t23;
t270 = t207 + t260;
t268 = t262 * t53 + t284 * t52;
t267 = -t140 * t59 - t342;
t266 = -t204 * t59 - t342;
t31 = t228 * t64 + t232 * t51;
t26 = -qJ(6) * t103 + t31;
t265 = -t228 * t26 - t347;
t264 = -qJ(6) * t87 - qJD(6) * t122;
t259 = -t202 * t209 - t203 * t227;
t258 = t11 * t228 + t203 * t358 + t262 * t31 + t50 * t307;
t257 = -t262 * t30 + t48 + (g(1) * t319 + g(2) * t321) * t202;
t256 = t232 * t59 - t374 * t382;
t255 = t232 * t12 + t228 * t24 + t64 * t307 - t51 * t308;
t254 = t273 * t202;
t253 = -0.2e1 * pkin(1) * t305 - pkin(7) * qJDD(2);
t252 = t228 * t87 + t295;
t251 = -t122 * t308 + t343;
t245 = -t180 * t152 + t123 + t385;
t236 = qJD(2) ^ 2;
t243 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t236 + t377;
t237 = qJD(1) ^ 2;
t242 = pkin(1) * t237 - pkin(7) * qJDD(1) + t273;
t7 = t46 * pkin(5) + qJDD(6) + t11;
t239 = g(3) * t216 - t180 * t150 + t217 * t273 - t373;
t1 = pkin(5) * t59 + qJ(6) * t45 - t31 * qJD(5) - qJD(6) * t105 - t12 * t228 + t23;
t3 = -qJ(6) * t46 - qJD(6) * t103 + t255;
t238 = t265 * qJD(5) - t1 * t228 - t203 * t273 + t3 * t232 + t26 * t331 + t284 * t347 - t361;
t205 = -pkin(4) - t368;
t177 = -pkin(2) * t230 - t369;
t172 = pkin(1) + t314;
t162 = t204 * t232 + t218;
t161 = t317 * t228;
t159 = t235 * t172;
t138 = t203 * t319 + t324;
t136 = -t203 * t321 + t323;
t129 = t140 * t232 + t218;
t128 = t318 * t228;
t109 = -t150 ^ 2 + t152 ^ 2;
t102 = t103 ^ 2;
t85 = -t269 + (-qJD(1) * t167 - t152) * t302;
t84 = t150 * t302 + t100;
t78 = t232 * t80;
t42 = t103 * pkin(5) + qJD(6) + t50;
t37 = t232 * t39;
t32 = -qJ(6) * t330 + t352;
t29 = pkin(5) * t121 - t122 * t218 - t228 * t82 + t78;
t19 = t105 * t279 - t345;
t16 = -t334 - t381;
t15 = t256 + t336;
t6 = (-t45 + t335) * t232 - t379 + t355;
t5 = -qJ(6) * t295 + (-qJD(5) * t82 + t264) * t228 + t300;
t4 = pkin(5) * t86 - t228 * t35 + t37 + t264 * t232 + (-t79 + (qJ(6) * t122 - t80) * t228) * qJD(5);
t2 = [qJDD(1), t377, t273, qJDD(1) * t222 + 0.2e1 * t230 * t290, 0.2e1 * t230 * t303 - 0.2e1 * t313 * t305, qJDD(2) * t230 + t234 * t236, qJDD(2) * t234 - t230 * t236, 0, t230 * t253 + t234 * t243, -t230 * t243 + t234 * t253, t100 * t167 + t124 * t152, -t100 * t166 - t101 * t167 + t124 * t150 + t125 * t152, -t124 * t302 + t167 * t301, -t125 * t302 - t166 * t301, 0, -t101 * t357 - t180 * t125 + t147 * t166 + t150 * t213 + t217 * t377 + t240 * t302 + t282 * t301, -t100 * t357 + t180 * t124 + t147 * t167 - t152 * t213 - t216 * t377 - t247 * t302 - t315 * t301, -t14 * t121 - t13 * t122 + t262 * t34 + t263 * t81 + t284 * t35 - t52 * t87 - t53 * t86 + t82 * t62 - t273, t14 * t82 + t53 * t35 - t13 * t81 - t52 * t34 + t241 * t276 + t127 * t289 - g(1) * (-t172 * t231 - t221 * t235) - g(2) * (-t221 * t231 + t159) t105 * t251 - t329 * t45 (-t103 * t232 - t333) * t87 + (t345 - t232 * t46 + (t103 * t228 - t105 * t232) * qJD(5)) * t122, t105 * t86 - t121 * t45 + t251 * t374 + t329 * t59, -t103 * t86 - t121 * t46 - t252 * t374 - t330 * t59, t121 * t59 + t374 * t86 (-t307 * t82 + t37) * t374 + t78 * t59 + t274 * t121 + t30 * t86 + t34 * t103 + t81 * t46 + t50 * t295 - g(1) * t136 - g(2) * t138 + ((-qJD(5) * t80 - t35) * t374 - t82 * t59 + t286 * t121 + t11 * t122 + t50 * t87) * t228 -(-t308 * t82 + t300) * t374 - t352 * t59 - t255 * t121 - t31 * t86 + t34 * t105 - t81 * t45 + t50 * t343 - g(1) * t135 - g(2) * t137 + (t11 * t232 - t48) * t122, -t103 * t5 - t105 * t4 + t29 * t45 - t32 * t46 + t265 * t87 + t377 * t202 + (-t1 * t232 - t228 * t3 + (t18 * t228 - t232 * t26) * qJD(5)) * t122, t3 * t32 + t26 * t5 + t1 * t29 + t18 * t4 + t7 * (pkin(5) * t330 + t81) + t42 * (pkin(5) * t252 + t34) - g(2) * t159 + (-g(1) * t288 - g(2) * t260) * t235 + (-g(1) * (-t172 - t260) - g(2) * t288) * t231; 0, 0, 0, -t230 * t237 * t234, t313 * t237, t304, t303, qJDD(2), -g(3) * t234 + t230 * t242, g(3) * t230 + t234 * t242, -t328, t109, t84, t85, t301, -t175 * t309 + t157 * t302 + (-t173 * t302 - t278) * t229 + (-t150 * t312 + t233 * t301 - t302 * t310) * pkin(2) + t245, t316 * t302 + (t152 * t312 - t229 * t301 - t302 * t309) * pkin(2) + t239, -t144 * t263 + t145 * t62 + t262 * t339 + t284 * t338 + t268, t14 * t145 + t13 * t144 - t127 * (t212 - t370) - g(3) * t314 + t338 * t53 - t339 * t52 - t273 * t177, t19, t6, t16, t15, -t332, t139 * t46 + t292 * t232 + t267 * t228 + t339 * t103 + (-t140 * t307 - t228 * t338 - t67) * t374 + t257, -t139 * t45 + t267 * t232 - t228 * t254 + t339 * t105 + (t140 * t308 + t376) * t374 + t258, -t103 * t351 - t105 * t350 + t128 * t45 - t129 * t46 + t238, t3 * t129 + t1 * t128 + t7 * (t139 - t367) - g(3) * (t219 + t270) + (t275 + t339) * t42 + t351 * t26 + t350 * t18 + t273 * (-t177 - t259); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t328, t109, t84, t85, t301, -qJD(2) * t261 + t245 - t322, t283 * t302 + t239, -t61 * t284 - t60 * t262 + (t225 * t62 - t226 * t263) * pkin(3) + t268, t52 * t60 - t53 * t61 + (t127 * t152 + t13 * t226 + t14 * t225 + t385) * pkin(3), t19, t6, t16, t15, -t332, -t60 * t103 - t74 * t374 + t205 * t46 + (t374 * t61 + t266) * t228 + (t292 - t294) * t232 + t257, -t205 * t45 + t354 * t374 - t60 * t105 + t266 * t232 + (-t254 + t294) * t228 + t258, -t103 * t341 - t105 * t340 + t161 * t45 - t162 * t46 + t238, t3 * t162 + t1 * t161 + t7 * (-t209 - t368) - g(3) * t270 + (t275 - t60) * t42 + t341 * t26 + t340 * t18 + t273 * (-t259 + t369); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262 ^ 2 - t284 ^ 2, t262 * t52 - t284 * t53 + t241 - t377, 0, 0, 0, 0, 0, t256 - t336, -t334 + t381 (t45 + t335) * t232 + t379 + t355, -t262 * t42 + (t26 * t374 + t1) * t232 + (-t18 * t374 + t3) * t228 - t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t103, -t102 + t372, t103 * t374 - t45, -t244 + (-qJD(5) + t374) * t105, t59, -t105 * t50 + t374 * t31 + (t286 + t361) * t228 + t274 + t375, g(1) * t138 - g(2) * t136 + t103 * t50 + t232 * t361 + t30 * t374 - t255, pkin(5) * t45 - t103 * t356, t356 * t26 + (-t42 * t105 + t202 * t358 + t1 + t375) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102 - t372, t26 * t103 + t18 * t105 - t254 + t360 + t7;];
tau_reg  = t2;
