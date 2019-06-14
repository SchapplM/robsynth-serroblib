% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 09:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:38:00
% EndTime: 2019-05-05 09:38:17
% DurationCPUTime: 7.63s
% Computational Cost: add. (20046->441), mult. (39688->593), div. (0->0), fcn. (29331->12), ass. (0->275)
t262 = sin(qJ(4));
t266 = cos(qJ(4));
t267 = cos(qJ(3));
t263 = sin(qJ(3));
t316 = qJD(2) * t263;
t225 = -t266 * t267 * qJD(2) + t262 * t316;
t221 = qJD(5) + t225;
t360 = t221 ^ 2;
t320 = t267 * t262;
t227 = (t263 * t266 + t320) * qJD(2);
t254 = qJD(3) + qJD(4);
t261 = sin(qJ(5));
t265 = cos(qJ(5));
t208 = t227 * t261 - t265 * t254;
t361 = t208 ^ 2;
t187 = t361 - t360;
t314 = qJD(2) * qJD(3);
t304 = t267 * t314;
t313 = t263 * qJDD(2);
t231 = t304 + t313;
t251 = t267 * qJDD(2);
t305 = t263 * t314;
t232 = t251 - t305;
t298 = t231 * t262 - t266 * t232;
t181 = -qJD(4) * t227 - t298;
t179 = qJDD(5) - t181;
t210 = t227 * t265 + t254 * t261;
t328 = t210 * t208;
t138 = -t328 - t179;
t339 = t138 * t261;
t105 = -t187 * t265 - t339;
t192 = t210 * t221;
t182 = -qJD(4) * t225 + t231 * t266 + t232 * t262;
t312 = qJDD(3) + qJDD(4);
t299 = -t182 * t261 + t265 * t312;
t288 = qJD(5) * t210 - t299;
t124 = -t192 + t288;
t427 = t263 * (t105 * t266 + t124 * t262) + t267 * (t105 * t262 - t124 * t266);
t207 = t210 ^ 2;
t371 = -t207 - t360;
t94 = t265 * t371 + t339;
t426 = pkin(2) * t94;
t425 = pkin(3) * t94;
t424 = pkin(4) * t94;
t423 = pkin(10) * t94;
t338 = t138 * t265;
t96 = -t261 * t371 + t338;
t422 = pkin(10) * t96;
t421 = t262 * t96;
t420 = t266 * t96;
t268 = cos(qJ(2));
t419 = t268 * t94;
t370 = t207 - t361;
t284 = -t265 * t182 - t261 * t312;
t279 = -t208 * qJD(5) - t284;
t329 = t208 * t221;
t377 = t329 - t279;
t341 = t377 * t261;
t373 = t192 + t288;
t76 = t373 * t265 - t341;
t416 = t263 * (-t262 * t370 + t266 * t76) + t267 * (t262 * t76 + t266 * t370);
t368 = -t328 + t179;
t336 = t368 * t265;
t365 = -t360 - t361;
t375 = t261 * t365 + t336;
t337 = t368 * t261;
t374 = t265 * t365 - t337;
t393 = t262 * t373 + t266 * t374;
t394 = t262 * t374 - t266 * t373;
t407 = -t263 * t394 + t267 * t393;
t415 = -pkin(2) * t375 + pkin(8) * t407;
t258 = sin(pkin(6));
t259 = cos(pkin(6));
t264 = sin(qJ(2));
t412 = t259 * (t263 * t393 + t267 * t394) + (t264 * t407 - t268 * t375) * t258;
t411 = pkin(3) * t394;
t410 = pkin(9) * t394;
t409 = t377 * qJ(6);
t101 = -t187 * t261 + t338;
t408 = -pkin(3) * t375 + pkin(9) * t393;
t366 = t329 + t279;
t188 = -t207 + t360;
t395 = -t188 * t261 + t336;
t405 = t263 * (t262 * t366 + t266 * t395) + t267 * (t262 * t395 - t266 * t366);
t402 = pkin(4) * t375;
t401 = pkin(10) * t374;
t400 = pkin(10) * t375;
t396 = t265 * t188 + t337;
t340 = t377 * t265;
t75 = -t373 * t261 - t340;
t369 = t207 + t361;
t392 = pkin(4) * t369;
t391 = -qJ(6) * t261 - pkin(4);
t388 = t262 * t369;
t203 = t227 * t225;
t372 = -t203 + t312;
t386 = t262 * t372;
t382 = t266 * t369;
t380 = t266 * t372;
t352 = sin(pkin(11));
t353 = cos(pkin(11));
t283 = t352 * g(1) - t353 * g(2);
t280 = t259 * t283;
t318 = -g(3) + qJDD(1);
t378 = t258 * t318 + t280;
t236 = -t353 * g(1) - t352 * g(2);
t195 = t268 * t236 + t378 * t264;
t376 = qJDD(2) * pkin(8) + t195;
t200 = pkin(4) * t225 - pkin(10) * t227;
t359 = t254 ^ 2;
t358 = qJD(2) ^ 2;
t272 = -t358 * pkin(2) + t376;
t271 = t263 * t272;
t275 = -t258 * t283 + t259 * t318;
t240 = qJD(3) * pkin(3) - pkin(9) * t316;
t256 = t267 ^ 2;
t273 = t263 * t275;
t133 = t267 * t376 + t273 + t232 * pkin(9) - qJD(3) * t240 + (-t267 * pkin(2) - t256 * pkin(3)) * t358;
t322 = t266 * t133;
t348 = qJDD(3) * pkin(3);
t356 = t231 * pkin(9);
t85 = t322 + t262 * (-t271 + t348 - t356) + ((pkin(3) * t316 + qJD(3) * pkin(9)) * qJD(2) + t275) * t320;
t69 = -t359 * pkin(4) + t312 * pkin(10) - t225 * t200 + t85;
t293 = t236 * t264 - t378 * t268;
t185 = -qJDD(2) * pkin(2) - t358 * pkin(8) + t293;
t253 = t256 * t358;
t154 = -t232 * pkin(3) - pkin(9) * t253 + t240 * t316 + t185;
t218 = t254 * t225;
t164 = t182 - t218;
t83 = -t164 * pkin(10) + (t227 * t254 - t181) * pkin(4) + t154;
t45 = t261 * t69 - t265 * t83;
t46 = t261 * t83 + t265 * t69;
t19 = t261 * t45 + t265 * t46;
t170 = pkin(5) * t208 - qJ(6) * t210;
t297 = t179 * qJ(6) - t208 * t170 + t46;
t364 = -pkin(5) * (t371 + t360) - qJ(6) * t138 + t297;
t326 = t221 * t265;
t307 = t208 * t326;
t289 = t261 * t288 + t307;
t308 = t266 * t328;
t309 = t262 * t328;
t363 = t263 * (t266 * t289 - t309) + t267 * (t262 * t289 + t308);
t327 = t221 * t261;
t184 = t210 * t327;
t294 = t184 - t307;
t362 = t263 * (t179 * t262 + t266 * t294) + t267 * (-t266 * t179 + t262 * t294);
t223 = t225 ^ 2;
t224 = t227 ^ 2;
t357 = pkin(4) * t262;
t243 = t263 * t358 * t267;
t156 = -t267 * t275 + t271;
t270 = pkin(9) * t304 - t156 - t356;
t84 = t133 * t262 - t266 * (pkin(3) * t243 + t270 + t348);
t68 = -t312 * pkin(4) - t359 * pkin(10) + t200 * t227 + t84;
t355 = -pkin(4) * t68 + pkin(10) * t19;
t64 = t261 * t68;
t49 = t262 * t85 - t266 * t84;
t354 = t263 * t49;
t65 = t265 * t68;
t350 = qJ(6) * t265;
t343 = t366 * t261;
t342 = t366 * t265;
t335 = t154 * t262;
t334 = t154 * t266;
t198 = t203 + t312;
t331 = t198 * t262;
t330 = t198 * t266;
t325 = t254 * t262;
t324 = t254 * t266;
t237 = qJDD(3) + t243;
t323 = t263 * t237;
t238 = qJDD(3) - t243;
t321 = t267 * t238;
t315 = qJD(6) * t221;
t131 = (qJD(5) + t221) * t208 + t284;
t311 = pkin(4) * t131 + t422 + t64;
t310 = -pkin(4) * t373 + t401 - t65;
t306 = -pkin(4) * t266 - pkin(3);
t50 = t262 * t84 + t266 * t85;
t216 = 0.2e1 * t315;
t292 = t216 + t297;
t23 = (t369 - t360) * pkin(5) + t292;
t35 = -t179 * pkin(5) - qJ(6) * t360 + t170 * t210 + qJDD(6) + t45;
t25 = qJ(6) * t369 + t35;
t78 = -t124 * t265 + t343;
t302 = pkin(10) * t78 + t265 * t23 + t261 * t25 + t392;
t126 = (-qJD(5) + t221) * t210 + t299;
t77 = t126 * t265 + t343;
t301 = pkin(10) * t77 + t19 + t392;
t157 = t267 * t272 + t273;
t107 = t156 * t263 + t267 * t157;
t34 = -pkin(5) * t360 + t292;
t296 = -pkin(5) * t35 + qJ(6) * t34;
t295 = t208 * t327 - t265 * t288;
t291 = -pkin(5) * t366 - qJ(6) * t124;
t11 = t19 * t262 - t266 * t68;
t12 = t19 * t266 + t262 * t68;
t3 = -t263 * t11 + t267 * t12;
t18 = t261 * t46 - t265 * t45;
t278 = t288 * pkin(5) + t409 + t68;
t277 = 0.2e1 * qJD(6) * t210 - t278;
t30 = -pkin(5) * t192 + t277 - t409;
t290 = -pkin(4) * t377 - pkin(5) * t340 + t261 * t30 - t422;
t31 = (-t373 - t192) * pkin(5) + t277;
t287 = t265 * t31 + t391 * t373 + t401;
t285 = (-t208 * t261 - t210 * t265) * t221;
t282 = (-qJD(4) + t254) * t227 - t298;
t10 = t261 * t35 + t265 * t34;
t41 = (pkin(5) * t221 - 0.2e1 * qJD(6)) * t210 + t278;
t281 = pkin(10) * t10 + (-pkin(5) * t265 + t391) * t41;
t118 = t265 * t279 - t184;
t276 = t263 * (t266 * t118 + t309) + t267 * (t262 * t118 - t308);
t274 = pkin(5) * t368 + qJ(6) * t365 - t35;
t269 = qJD(3) ^ 2;
t255 = t263 ^ 2;
t252 = t255 * t358;
t242 = -t253 - t269;
t241 = -t252 - t269;
t235 = t252 + t253;
t234 = (t255 + t256) * qJDD(2);
t233 = t251 - 0.2e1 * t305;
t230 = 0.2e1 * t304 + t313;
t215 = -t224 + t359;
t214 = t223 - t359;
t213 = -t224 - t359;
t212 = -t241 * t263 - t321;
t211 = t242 * t267 - t323;
t202 = t224 - t223;
t196 = -t359 - t223;
t183 = -t223 - t224;
t167 = -t213 * t262 - t330;
t166 = t213 * t266 - t331;
t165 = t182 + t218;
t160 = (qJD(4) + t254) * t227 + t298;
t153 = t196 * t266 - t386;
t152 = t196 * t262 + t380;
t117 = t210 * t326 + t261 * t279;
t110 = -t263 * t166 + t267 * t167;
t109 = t165 * t262 + t266 * t282;
t108 = -t165 * t266 + t262 * t282;
t98 = -t263 * t152 + t267 * t153;
t74 = -t124 * t261 - t342;
t73 = t126 * t261 - t342;
t63 = -t263 * t108 + t267 * t109;
t61 = -t131 * t262 + t420;
t59 = t131 * t266 + t421;
t57 = t262 * t377 - t420;
t55 = -t266 * t377 - t421;
t54 = t266 * t78 - t388;
t53 = t266 * t77 - t388;
t52 = t262 * t78 + t382;
t51 = t262 * t77 + t382;
t48 = t65 - t423;
t47 = t64 - t400;
t42 = -pkin(4) * t74 - t291;
t38 = -t263 * t59 + t267 * t61;
t36 = -t263 * t55 + t267 * t57;
t33 = t46 - t424;
t32 = t45 - t402;
t27 = -t263 * t52 + t267 * t54;
t26 = -t263 * t51 + t267 * t53;
t21 = t267 * t50 - t354;
t20 = -t274 - t402;
t16 = -t261 * t31 - t350 * t373 - t400;
t15 = pkin(5) * t341 + t265 * t30 + t423;
t14 = -0.2e1 * t315 - t364 + t424;
t13 = -pkin(10) * t73 - t18;
t9 = t261 * t34 - t265 * t35;
t7 = -pkin(10) * t74 - t23 * t261 + t25 * t265;
t6 = t10 * t266 + t262 * t41;
t5 = t10 * t262 - t266 * t41;
t4 = -pkin(10) * t9 + (pkin(5) * t261 - t350) * t41;
t2 = -pkin(4) * t9 - t296;
t1 = -t263 * t5 + t267 * t6;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t318, 0, 0, 0, 0, 0, 0, (qJDD(2) * t268 - t264 * t358) * t258, (-qJDD(2) * t264 - t268 * t358) * t258, 0, t259 ^ 2 * t318 + (t264 * t195 - t268 * t293 - t280) * t258, 0, 0, 0, 0, 0, 0, t259 * (t237 * t267 + t242 * t263) + (t264 * t211 + t268 * t233) * t258, t259 * (-t238 * t263 + t241 * t267) + (t264 * t212 - t268 * t230) * t258, (t234 * t264 + t235 * t268) * t258, t259 * (-t156 * t267 + t157 * t263) + (t264 * t107 - t268 * t185) * t258, 0, 0, 0, 0, 0, 0, t259 * (t152 * t267 + t153 * t263) + (-t268 * t160 + t264 * t98) * t258, t259 * (t166 * t267 + t167 * t263) + (t264 * t110 - t164 * t268) * t258, t259 * (t108 * t267 + t109 * t263) + (-t268 * t183 + t264 * t63) * t258, t259 * (t263 * t50 + t267 * t49) + (-t268 * t154 + t264 * t21) * t258, 0, 0, 0, 0, 0, 0, t412, t259 * (t263 * t61 + t267 * t59) + (t264 * t38 - t419) * t258, t259 * (t263 * t53 + t267 * t51) + (t264 * t26 - t268 * t73) * t258, t259 * (t11 * t267 + t12 * t263) + (-t268 * t18 + t264 * t3) * t258, 0, 0, 0, 0, 0, 0, t412, t259 * (t263 * t54 + t267 * t52) + (t264 * t27 - t268 * t74) * t258, t259 * (t263 * t57 + t267 * t55) + (t264 * t36 + t419) * t258, t259 * (t263 * t6 + t267 * t5) + (t264 * t1 - t268 * t9) * t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t293, -t195, 0, 0, (t231 + t304) * t263, t230 * t267 + t233 * t263, t323 + t267 * (-t252 + t269), (t232 - t305) * t267, t263 * (t253 - t269) + t321, 0, pkin(2) * t233 + pkin(8) * t211 - t185 * t267, -pkin(2) * t230 + pkin(8) * t212 + t185 * t263, pkin(2) * t235 + pkin(8) * t234 + t107, -pkin(2) * t185 + pkin(8) * t107, t263 * (t182 * t266 - t227 * t325) + t267 * (t182 * t262 + t227 * t324), t263 * (-t160 * t266 - t164 * t262) + t267 * (-t160 * t262 + t164 * t266), t263 * (-t215 * t262 + t380) + t267 * (t215 * t266 + t386), t263 * (-t181 * t262 + t225 * t324) + t267 * (t181 * t266 + t225 * t325), t263 * (t214 * t266 - t331) + t267 * (t214 * t262 + t330), (t263 * (-t225 * t266 + t227 * t262) + t267 * (-t225 * t262 - t227 * t266)) * t254, t263 * (-pkin(9) * t152 + t335) + t267 * (-pkin(3) * t160 + pkin(9) * t153 - t334) - pkin(2) * t160 + pkin(8) * t98, t263 * (-pkin(9) * t166 + t334) + t267 * (-pkin(3) * t164 + pkin(9) * t167 + t335) - pkin(2) * t164 + pkin(8) * t110, t263 * (-pkin(9) * t108 - t49) + t267 * (-pkin(3) * t183 + pkin(9) * t109 + t50) - pkin(2) * t183 + pkin(8) * t63, -pkin(9) * t354 + t267 * (-pkin(3) * t154 + pkin(9) * t50) - pkin(2) * t154 + pkin(8) * t21, t276, -t416, t405, t363, -t427, t362, t263 * (-t262 * t32 + t266 * t47 - t410) + t267 * (t262 * t47 + t266 * t32 + t408) + t415, t263 * (-pkin(9) * t59 - t262 * t33 + t266 * t48) + t267 * (pkin(9) * t61 + t262 * t48 + t266 * t33 - t425) - t426 + pkin(8) * t38, t263 * (-pkin(9) * t51 + t13 * t266) + t267 * (pkin(9) * t53 + t13 * t262) + pkin(8) * t26 + (t263 * t357 + t267 * t306 - pkin(2)) * t73, (t263 * (-pkin(10) * t266 + t357) + t267 * (-pkin(10) * t262 + t306) - pkin(2)) * t18 + (pkin(8) + pkin(9)) * t3, t276, t405, t416, t362, t427, t363, t263 * (t16 * t266 - t20 * t262 - t410) + t267 * (t16 * t262 + t20 * t266 + t408) + t415, t263 * (-pkin(9) * t52 - t262 * t42 + t266 * t7) + t267 * (-pkin(3) * t74 + pkin(9) * t54 + t262 * t7 + t266 * t42) - pkin(2) * t74 + pkin(8) * t27, t263 * (-pkin(9) * t55 - t14 * t262 + t15 * t266) + t267 * (pkin(9) * t57 + t14 * t266 + t15 * t262 + t425) + t426 + pkin(8) * t36, t263 * (-pkin(9) * t5 - t2 * t262 + t266 * t4) + t267 * (-pkin(3) * t9 + pkin(9) * t6 + t2 * t266 + t262 * t4) - pkin(2) * t9 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, t252 - t253, t313, t243, t251, qJDD(3), -t156, -t157, 0, 0, t203, t202, t165, -t203, t282, t312, pkin(3) * t152 - t84, -t322 - t262 * t270 + (-t237 * t262 + t166) * pkin(3), pkin(3) * t108, pkin(3) * t49, t117, t75, t396, t295, -t101, t285, t310 + t411, pkin(3) * t59 + t311, pkin(3) * t51 + t301, pkin(3) * t11 + t355, t117, t396, -t75, t285, t101, t295, t287 + t411, pkin(3) * t52 + t302, pkin(3) * t55 + t290, pkin(3) * t5 + t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, t202, t165, -t203, t282, t312, -t84, -t85, 0, 0, t117, t75, t396, t295, -t101, t285, t310, t311, t301, t355, t117, t396, -t75, t285, t101, t295, t287, t302, t290, t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t328, t370, t366, -t328, -t124, t179, -t45, -t46, 0, 0, t328, t366, -t370, t179, t124, -t328, t274, t291, t216 + t364, t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t368, t366, t371, t35;];
tauJ_reg  = t8;
