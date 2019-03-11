% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:21:03
% EndTime: 2019-03-09 09:21:19
% DurationCPUTime: 7.12s
% Computational Cost: add. (4648->567), mult. (11143->739), div. (0->0), fcn. (8359->10), ass. (0->285)
t193 = sin(pkin(6));
t202 = cos(qJ(2));
t316 = qJD(1) * qJD(2);
t290 = t202 * t316;
t198 = sin(qJ(2));
t313 = qJDD(1) * t198;
t228 = t290 + t313;
t391 = t193 * t228;
t328 = qJD(1) * t193;
t297 = t198 * t328;
t194 = cos(pkin(6));
t319 = t194 * qJD(1);
t309 = pkin(1) * t319;
t333 = pkin(8) * t297 - t202 * t309;
t389 = qJD(3) + t333;
t393 = -qJ(4) * t297 + t389;
t170 = qJD(2) + t319;
t201 = cos(qJ(5));
t197 = sin(qJ(5));
t327 = qJD(1) * t202;
t296 = t193 * t327;
t272 = t197 * t296;
t97 = -t170 * t201 + t272;
t335 = qJD(6) - t97;
t196 = sin(qJ(6));
t200 = cos(qJ(6));
t388 = -t297 - qJD(5);
t98 = t197 * t170 + t201 * t296;
t60 = -t196 * t98 + t200 * t388;
t392 = t335 * t60;
t69 = ((qJD(2) - t170) * t327 + t313) * t193;
t291 = t198 * t316;
t270 = t193 * t291;
t312 = qJDD(1) * t202;
t289 = t193 * t312;
t390 = t270 - t289;
t199 = sin(qJ(1));
t341 = t199 * t202;
t203 = cos(qJ(1));
t342 = t198 * t203;
t115 = t194 * t342 + t341;
t338 = t202 * t203;
t344 = t198 * t199;
t114 = -t194 * t338 + t344;
t346 = t193 * t203;
t78 = t114 * t201 + t197 * t346;
t387 = -t115 * t200 + t196 * t78;
t386 = t115 * t196 + t200 * t78;
t189 = t193 ^ 2;
t311 = 0.2e1 * t189;
t62 = -t196 * t388 - t200 * t98;
t385 = t388 * t62;
t314 = qJDD(1) * t194;
t168 = qJDD(2) + t314;
t204 = -pkin(2) - pkin(3);
t307 = t168 * t204;
t227 = t291 - t312;
t384 = t193 * t227;
t116 = t194 * t341 + t342;
t347 = t193 * t202;
t383 = g(1) * t116 + g(2) * t114 - g(3) * t347;
t382 = qJD(2) * t204;
t112 = -t194 * t201 + t197 * t347;
t339 = t201 * t203;
t356 = t114 * t197;
t348 = t193 * t199;
t80 = -t116 * t197 - t201 * t348;
t229 = g(1) * t80 + g(2) * (t193 * t339 - t356) + g(3) * t112;
t102 = qJDD(5) + t391;
t190 = -pkin(9) + t204;
t42 = t170 * t190 + t393;
t268 = pkin(4) * t198 + pkin(9) * t202;
t242 = t268 * t193;
t90 = -pkin(1) * t328 - pkin(2) * t296 - qJ(3) * t297;
t68 = pkin(3) * t296 + qJD(4) - t90;
t47 = qJD(1) * t242 + t68;
t18 = t197 * t47 + t201 * t42;
t230 = pkin(4) * t202 + t190 * t198;
t220 = qJD(2) * t230;
t349 = t193 * t198;
t166 = qJD(3) * t349;
t186 = t193 * pkin(1);
t252 = pkin(2) * t289 + qJ(3) * t391 + qJD(1) * t166 + qJDD(1) * t186;
t223 = pkin(3) * t289 + qJDD(4) + t252;
t19 = (qJD(1) * t220 + qJDD(1) * t268) * t193 + t223;
t308 = pkin(1) * qJD(2) * t194;
t276 = qJD(1) * t308;
t306 = pkin(1) * t314;
t273 = pkin(8) * t391 + t198 * t276 - t202 * t306;
t248 = -qJDD(3) - t273;
t315 = qJD(1) * qJD(4);
t206 = (-qJ(4) * t228 - t198 * t315) * t193 - t248;
t22 = t168 * t190 + t206;
t214 = -qJD(5) * t18 + t201 * t19 - t197 * t22;
t4 = -t102 * pkin(5) - t214;
t381 = (-pkin(5) * t98 + pkin(10) * t335) * t335 + t229 + t4;
t195 = qJ(3) + pkin(4);
t127 = pkin(5) * t201 + pkin(10) * t197 + t195;
t322 = qJD(5) * t201;
t323 = qJD(5) * t197;
t38 = -t201 * t168 + t170 * t323 - t197 * t270 + (qJD(1) * t322 + qJDD(1) * t197) * t347;
t36 = -qJDD(6) + t38;
t380 = (t393 + t388 * (pkin(5) * t197 - pkin(10) * t201)) * t335 - t127 * t36;
t165 = t202 * t308;
t183 = t194 * qJD(3);
t326 = qJD(2) * t198;
t295 = t193 * t326;
t58 = t193 * (pkin(8) * t326 + qJD(4) * t202) - qJ(4) * t295 - t165 - t183;
t130 = qJ(4) * t270;
t140 = t168 * qJ(3);
t141 = t170 * qJD(3);
t274 = pkin(8) * t390 - t198 * t306 - t202 * t276;
t39 = t140 + t141 - t274;
t30 = (qJ(4) * qJDD(1) + t315) * t347 - t130 - t39;
t332 = pkin(2) * t347 + qJ(3) * t349;
t100 = -t186 - t332;
t174 = pkin(3) * t347;
t84 = t174 - t100;
t57 = t242 + t84;
t172 = pkin(8) * t349;
t260 = -qJ(4) * t349 + t172;
t298 = -pkin(1) * t202 - pkin(2);
t279 = -pkin(3) + t298;
t64 = (-pkin(9) + t279) * t194 + t260;
t253 = t197 * t57 + t201 * t64;
t325 = qJD(2) * t202;
t294 = t193 * t325;
t334 = qJ(3) * t294 + t166;
t44 = t193 * t220 + t334;
t376 = pkin(1) * t198;
t178 = t194 * t376;
t67 = -qJD(4) * t349 + (t178 + (pkin(8) - qJ(4)) * t347) * qJD(2);
t378 = -qJD(5) * t253 - t197 * t67 + t201 * t44;
t14 = -pkin(10) * t388 + t18;
t144 = t170 * qJ(3);
t106 = pkin(8) * t296 + t198 * t309;
t88 = -qJ(4) * t296 + t106;
t65 = -t144 - t88;
t55 = pkin(4) * t170 - t65;
t26 = -pkin(5) * t97 + pkin(10) * t98 + t55;
t258 = t14 * t196 - t200 * t26;
t238 = t197 * t19 + t201 * t22 + t47 * t322 - t323 * t42;
t3 = pkin(10) * t102 + t238;
t23 = t168 * pkin(4) - t30;
t37 = qJD(5) * t272 - t197 * t168 - t170 * t322 + t201 * t390;
t8 = -pkin(5) * t38 - pkin(10) * t37 + t23;
t1 = -t258 * qJD(6) + t196 * t8 + t200 * t3;
t377 = 0.2e1 * t140;
t375 = g(3) * t198;
t374 = t168 * pkin(2);
t91 = t196 * t201 * t297 - t200 * t296;
t373 = t91 * t335;
t343 = t198 * t201;
t92 = (t196 * t202 + t200 * t343) * t328;
t372 = t92 * t335;
t147 = qJ(3) * t296;
t56 = t230 * t328 + t147;
t371 = t197 * t56 + t201 * t88;
t320 = qJD(6) * t200;
t321 = qJD(6) * t196;
t11 = t196 * t102 + t200 * t37 - t320 * t388 + t321 * t98;
t370 = t11 * t196;
t368 = t388 * t97;
t367 = t196 * t36;
t366 = t196 * t335;
t365 = t200 * t36;
t364 = t200 * t62;
t363 = t200 * t335;
t362 = t202 * t97;
t361 = t202 * t98;
t360 = t98 * t388;
t358 = qJ(3) * t202;
t357 = qJD(5) * t335;
t231 = t388 * t197;
t353 = t168 * t194;
t352 = t189 * qJD(1) ^ 2;
t351 = t190 * t196;
t350 = t190 * t200;
t345 = t197 * t102;
t340 = t201 * t102;
t336 = -qJD(4) - t68;
t331 = pkin(8) * t347 + t178;
t191 = t198 ^ 2;
t192 = t202 ^ 2;
t330 = t191 - t192;
t329 = qJ(3) * qJD(2);
t324 = qJD(5) * t190;
t310 = g(3) * t349;
t305 = t196 * t357;
t304 = t200 * t357;
t303 = t198 * t231;
t302 = t388 * t343;
t301 = t202 * t352;
t300 = t201 * t347;
t99 = t194 * qJ(3) + t331;
t292 = qJ(4) * t313;
t285 = -t200 * t102 + t196 * t37;
t284 = t336 * t198;
t283 = -t114 * pkin(2) + qJ(3) * t115;
t117 = -t194 * t344 + t338;
t282 = -t116 * pkin(2) + qJ(3) * t117;
t281 = qJD(1) * t100 + t90;
t280 = t170 + t319;
t278 = t168 + t314;
t277 = t204 * t349;
t275 = t198 * t301;
t265 = g(1) * t114 - g(2) * t116;
t264 = -g(1) * t117 - g(2) * t115;
t263 = -g(1) * t115 + g(2) * t117;
t262 = g(1) * t203 + g(2) * t199;
t261 = g(1) * t199 - g(2) * t203;
t6 = t14 * t200 + t196 * t26;
t25 = pkin(10) * t349 + t253;
t113 = t194 * t197 + t300;
t83 = qJ(4) * t347 - t99;
t70 = t194 * pkin(4) - t83;
t34 = -pkin(5) * t112 + pkin(10) * t113 + t70;
t257 = t196 * t34 + t200 * t25;
t256 = -t196 * t25 + t200 * t34;
t17 = -t197 * t42 + t201 * t47;
t254 = -t197 * t64 + t201 * t57;
t251 = qJD(2) * (qJD(1) * t84 + t68);
t250 = qJD(2) * t277;
t249 = -t170 ^ 2 - t191 * t352;
t247 = -pkin(8) * t295 + t165;
t246 = t203 * pkin(1) + t117 * pkin(2) + pkin(8) * t348 + qJ(3) * t116;
t244 = t320 * t335 - t367;
t243 = t321 * t335 + t365;
t240 = t113 * t196 + t200 * t349;
t76 = -t113 * t200 + t196 * t349;
t31 = qJD(1) * t250 + t223;
t66 = t250 + t334;
t239 = qJD(1) * t66 + qJDD(1) * t84 + t31;
t237 = t197 * t44 + t201 * t67 + t57 * t322 - t323 * t64;
t236 = t335 * t388;
t235 = t388 * t60;
t234 = t322 * t388 - t345;
t233 = t323 * t388 + t340;
t43 = pkin(2) * t270 - t252;
t89 = pkin(2) * t295 - t334;
t232 = -qJD(1) * t89 - qJDD(1) * t100 - t43;
t226 = -pkin(1) * t199 - t115 * pkin(2) + pkin(8) * t346 - qJ(3) * t114;
t107 = t331 * qJD(2);
t224 = -t107 * t170 - t263;
t221 = -t264 + t310;
t13 = pkin(5) * t388 - t17;
t219 = pkin(10) * t36 + (t13 + t17) * t335;
t218 = t264 - t274;
t217 = -t221 + t23;
t216 = -t190 * t102 + t388 * t55;
t215 = -t273 + t383;
t2 = -qJD(6) * t6 - t196 * t3 + t200 * t8;
t213 = -qJDD(3) + t215;
t212 = t170 * t333 + t218;
t211 = qJD(6) * t190 * t335 + t221;
t210 = t106 * t170 + t215;
t209 = t243 - t385;
t208 = -t235 - t244;
t207 = -(-pkin(10) * t296 + qJD(6) * t127 - t371) * t335 + t383;
t118 = t170 * t297;
t104 = t168 + t275;
t103 = pkin(2) * t297 - t147;
t101 = t194 * t298 + t172;
t94 = t183 + t247;
t87 = qJD(1) * t277 + t147;
t85 = t144 + t106;
t82 = -pkin(2) * t170 + t389;
t81 = t116 * t201 - t197 * t348;
t74 = qJD(5) * t112 + t201 * t295;
t73 = -qJD(5) * t300 - t194 * t323 + t197 * t295;
t71 = t194 * t279 + t260;
t53 = t170 * t204 + t393;
t49 = -t248 - t374;
t46 = t117 * t196 + t200 * t81;
t45 = t117 * t200 - t196 * t81;
t33 = qJD(6) * t240 + t196 * t294 + t74 * t200;
t32 = qJD(6) * t76 + t74 * t196 - t200 * t294;
t29 = t307 + t206;
t27 = -pkin(5) * t296 + t197 * t88 - t201 * t56;
t24 = -pkin(5) * t349 - t254;
t20 = t73 * pkin(5) - t74 * pkin(10) - t58;
t12 = qJD(6) * t62 + t285;
t10 = -pkin(5) * t294 - t378;
t9 = pkin(10) * t294 + t237;
t5 = [qJDD(1), t261, t262 (qJDD(1) * t191 + 0.2e1 * t198 * t290) * t189 (t198 * t312 - t316 * t330) * t311 (t198 * t278 + t280 * t325) * t193 (t202 * t278 - t280 * t326) * t193, t353, -t172 * t168 - t273 * t194 + (t202 * t353 - t227 * t311) * pkin(1) + t224, -pkin(1) * t228 * t311 - t168 * t331 - t170 * t247 + t194 * t274 - t265, -t101 * t168 - t49 * t194 + (t202 * t232 + t281 * t326) * t193 + t224 ((qJD(2) * t82 + qJDD(1) * t99 + t39 + (qJD(2) * t101 + t94) * qJD(1)) * t202 + (-qJD(2) * t85 + qJDD(1) * t101 + t49 + (-qJD(2) * t99 + t107) * qJD(1)) * t198 - t262) * t193, t99 * t168 + t94 * t170 + t39 * t194 + (t198 * t232 - t281 * t325) * t193 + t265, -g(1) * t226 - g(2) * t246 + t43 * t100 + t49 * t101 + t82 * t107 + t39 * t99 + t85 * t94 + t90 * t89, -t83 * t168 - t58 * t170 - t30 * t194 + (t198 * t239 + t202 * t251) * t193 + t265, t71 * t168 + t67 * t170 + t29 * t194 + (t198 * t251 - t202 * t239) * t193 + t263 ((-qJD(2) * t53 + qJDD(1) * t83 + t30 + (-qJD(2) * t71 + t58) * qJD(1)) * t202 + (-qJD(2) * t65 - qJDD(1) * t71 - t29 + (-qJD(2) * t83 - t67) * qJD(1)) * t198 + t262) * t193, t29 * t71 + t53 * t67 + t30 * t83 + t65 * t58 + t31 * t84 + t68 * t66 - g(1) * (-pkin(3) * t115 - qJ(4) * t346 + t226) - g(2) * (pkin(3) * t117 - qJ(4) * t348 + t246) -t113 * t37 - t74 * t98, t112 * t37 - t113 * t38 + t73 * t98 + t74 * t97, -t113 * t102 - t74 * t388 + (t198 * t37 - t325 * t98) * t193, t112 * t102 + t73 * t388 + (t198 * t38 + t325 * t97) * t193 (t102 * t198 - t325 * t388) * t193, -t378 * t388 + t254 * t102 + t58 * t97 - t70 * t38 - t23 * t112 + t55 * t73 + g(1) * t78 - g(2) * t81 + (t17 * t325 + t198 * t214) * t193, t237 * t388 - t253 * t102 + t58 * t98 + t70 * t37 - t23 * t113 + t55 * t74 - g(1) * t356 - g(2) * t80 + (g(1) * t339 - t18 * t325 - t198 * t238) * t193, t11 * t76 + t33 * t62, t11 * t240 - t12 * t76 - t32 * t62 - t33 * t60, -t11 * t112 + t33 * t335 - t36 * t76 + t62 * t73, t112 * t12 - t240 * t36 - t32 * t335 - t60 * t73, t112 * t36 + t335 * t73 (-qJD(6) * t257 - t196 * t9 + t200 * t20) * t335 - t256 * t36 - t2 * t112 - t258 * t73 + t10 * t60 + t24 * t12 - t4 * t240 + t13 * t32 + g(1) * t386 - g(2) * t46 -(qJD(6) * t256 + t196 * t20 + t200 * t9) * t335 + t257 * t36 + t1 * t112 - t6 * t73 + t10 * t62 + t24 * t11 + t4 * t76 + t13 * t33 - g(1) * t387 - g(2) * t45; 0, 0, 0, -t275, t330 * t352, t69, t118 - t384, t168, t352 * t376 + t210, pkin(1) * t301 - t212 + t310, 0.2e1 * t374 - qJDD(3) + (t103 * t202 - t198 * t90) * t328 + t210 ((-pkin(2) * t198 + t358) * qJDD(1) + ((-t106 + t85 - t329) * t198 + (-pkin(2) * qJD(2) + t389 - t82) * t202) * qJD(1)) * t193, t377 + 0.2e1 * t141 + (-t375 + (t103 * t198 + t202 * t90) * qJD(1)) * t193 + t212, -t49 * pkin(2) - g(1) * t282 - g(2) * t283 - g(3) * t332 + t39 * qJ(3) - t90 * t103 - t82 * t106 + t389 * t85, t377 + t130 + t141 + t393 * t170 + (-qJ(4) * t312 - t375 + (-t198 * t87 + t202 * t336) * qJD(1)) * t193 + t218, -t170 * t88 + 0.2e1 * t307 + (-t292 + ((-qJ(4) * qJD(2) + t87) * t202 + t284) * qJD(1)) * t193 - t213 ((-t198 * t204 - t358) * qJDD(1) + ((t65 + t88 + t329) * t198 + (-t393 + t53 - t382) * t202) * qJD(1)) * t193, t29 * t204 - t30 * qJ(3) - t53 * t88 - t68 * t87 - g(1) * (-pkin(3) * t116 + t282) - g(2) * (-pkin(3) * t114 + t283) - g(3) * (t174 + t332) - t393 * t65, -t37 * t197 - t201 * t360 (-t37 + t368) * t201 + (-t38 + t360) * t197 (t302 + t361) * t328 + t234 (-t303 - t362) * t328 - t233, t388 * t296, -t17 * t296 - t195 * t38 - t393 * t97 + (-t388 * t88 + t216) * t197 + (-(-t56 - t324) * t388 + t217) * t201, t195 * t37 - t371 * t388 + t18 * t296 - t393 * t98 + t216 * t201 + (-t324 * t388 - t217) * t197, -t11 * t200 * t197 + (t197 * t321 - t200 * t322 - t92) * t62, t92 * t60 + t62 * t91 + (t196 * t62 + t200 * t60) * t322 + (t370 + t12 * t200 + (-t196 * t60 + t364) * qJD(6)) * t197, -t372 + (t11 - t304) * t201 + (t243 + t385) * t197, t373 + (-t12 + t305) * t201 + (-t235 + t244) * t197, -t36 * t201 + t231 * t335, -t13 * t91 - t27 * t60 + t380 * t200 + t207 * t196 + (t36 * t351 + t2 + (-t13 * t196 + t190 * t60) * qJD(5) - t211 * t200) * t201 + (t258 * t297 - t13 * t320 + t190 * t12 - t4 * t196 + (t335 * t351 + t258) * qJD(5)) * t197, -t13 * t92 - t27 * t62 - t380 * t196 + t207 * t200 + (t36 * t350 - t1 + (-t13 * t200 + t190 * t62) * qJD(5) + t211 * t196) * t201 + (t6 * t297 + t13 * t321 + t190 * t11 - t4 * t200 + (t335 * t350 + t6) * qJD(5)) * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t69, t249, -t170 * t85 + t297 * t90 - t213 - t374, t249, t104, -t69, t170 * t65 + t307 + (-t292 + (-qJ(4) * t325 + t284) * qJD(1)) * t193 - t213, 0, 0, 0, 0, 0, -t201 * t388 ^ 2 + t170 * t97 - t345, t170 * t98 + t231 * t388 - t340, 0, 0, 0, 0, 0, -t170 * t363 + (-t196 * t236 + t12) * t197 + t208 * t201, t170 * t366 + (-t200 * t236 + t11) * t197 + t209 * t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t313 + (qJD(2) + t170) * t327) * t193, t118 + t384 (-t191 - t192) * t352, g(3) * t194 + ((-t202 * t65 + (t53 + t382) * t198) * qJD(1) + t261) * t193 + t223, 0, 0, 0, 0, 0 (t303 - t362) * t328 + t233 (t302 - t361) * t328 + t234, 0, 0, 0, 0, 0, -t373 + (-t12 - t305) * t201 + t208 * t197, -t372 + (-t11 - t304) * t201 + t209 * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98 * t97, -t97 ^ 2 + t98 ^ 2, t37 + t368, t38 + t360, t102, -t18 * t388 + t55 * t98 + t214 - t229, g(1) * t81 + g(2) * t78 - g(3) * t113 - t17 * t388 - t55 * t97 - t238, t335 * t364 + t370 (t11 - t392) * t200 + (-t335 * t62 - t12) * t196, t335 * t363 + t62 * t98 - t367, -t335 * t366 - t60 * t98 - t365, t335 * t98, -pkin(5) * t12 - t18 * t60 + t219 * t196 - t200 * t381 - t258 * t98, -pkin(5) * t11 - t18 * t62 + t196 * t381 + t219 * t200 - t6 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60, -t60 ^ 2 + t62 ^ 2, t11 + t392, -t285 + (-qJD(6) + t335) * t62, -t36, -g(1) * t45 + g(2) * t387 - g(3) * t240 - t13 * t62 + t6 * t335 + t2, g(1) * t46 + g(2) * t386 + g(3) * t76 + t13 * t60 - t258 * t335 - t1;];
tau_reg  = t5;
