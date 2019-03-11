% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPR11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:07
% EndTime: 2019-03-09 11:15:23
% DurationCPUTime: 6.86s
% Computational Cost: add. (6291->525), mult. (13327->682), div. (0->0), fcn. (9011->12), ass. (0->293)
t252 = cos(qJ(4));
t248 = sin(qJ(4));
t346 = qJD(2) * t248;
t253 = cos(qJ(2));
t349 = qJD(1) * t253;
t174 = t252 * t349 + t346;
t312 = t248 * t349;
t344 = qJD(2) * t252;
t176 = -t312 + t344;
t244 = sin(pkin(10));
t245 = cos(pkin(10));
t109 = t174 * t244 - t176 * t245;
t247 = sin(qJ(6));
t251 = cos(qJ(6));
t288 = -t174 * t245 - t176 * t244;
t336 = qJD(6) * t251;
t337 = qJD(6) * t247;
t249 = sin(qJ(2));
t331 = qJD(1) * qJD(2);
t311 = t249 * t331;
t329 = t253 * qJDD(1);
t419 = t311 - t329;
t102 = -qJD(4) * t174 + t252 * qJDD(2) + t248 * t419;
t103 = -qJD(4) * t312 + t248 * qJDD(2) + (qJD(2) * qJD(4) - t419) * t252;
t47 = -t102 * t244 - t103 * t245;
t48 = t102 * t245 - t103 * t244;
t10 = t109 * t337 + t247 * t47 + t251 * t48 + t288 * t336;
t350 = qJD(1) * t249;
t207 = qJD(4) + t350;
t199 = qJD(6) + t207;
t54 = t109 * t247 + t251 * t288;
t387 = t199 * t54;
t427 = t10 - t387;
t415 = -t251 * t109 + t247 * t288;
t426 = t415 * t54;
t425 = t415 ^ 2 - t54 ^ 2;
t229 = qJ(4) + pkin(10) + qJ(6);
t215 = sin(t229);
t216 = cos(t229);
t250 = sin(qJ(1));
t254 = cos(qJ(1));
t367 = t249 * t254;
t125 = t215 * t367 + t216 * t250;
t369 = t249 * t250;
t127 = -t215 * t369 + t216 * t254;
t310 = t253 * t331;
t330 = t249 * qJDD(1);
t272 = t310 + t330;
t171 = qJDD(4) + t272;
t255 = -pkin(2) - pkin(8);
t205 = pkin(7) * t310;
t219 = pkin(7) * t330;
t309 = qJDD(3) + t205 + t219;
t106 = t272 * pkin(3) + qJDD(2) * t255 + t309;
t206 = pkin(2) * t311;
t381 = qJ(3) * t253;
t291 = pkin(8) * t249 - t381;
t335 = t249 * qJD(3);
t264 = qJD(2) * t291 - t335;
t230 = t249 * qJ(3);
t308 = -pkin(1) - t230;
t270 = t253 * t255 + t308;
t78 = qJD(1) * t264 + qJDD(1) * t270 + t206;
t305 = t252 * t106 - t248 * t78;
t136 = t270 * qJD(1);
t222 = pkin(7) * t350;
t416 = qJD(3) + t222;
t333 = pkin(3) * t350 + t416;
t141 = qJD(2) * t255 + t333;
t80 = t136 * t252 + t141 * t248;
t261 = -qJD(4) * t80 + t305;
t16 = t171 * pkin(4) - t102 * qJ(5) - t176 * qJD(5) + t261;
t340 = qJD(4) * t252;
t326 = -t248 * t106 - t141 * t340 - t252 * t78;
t341 = qJD(4) * t248;
t18 = -qJ(5) * t103 - qJD(5) * t174 - t136 * t341 - t326;
t4 = t245 * t16 - t18 * t244;
t2 = pkin(5) * t171 - pkin(9) * t48 + t4;
t69 = -qJ(5) * t174 + t80;
t386 = t245 * t69;
t79 = -t136 * t248 + t252 * t141;
t68 = -qJ(5) * t176 + t79;
t60 = pkin(4) * t207 + t68;
t26 = t244 * t60 + t386;
t414 = pkin(9) * t288;
t21 = t26 + t414;
t20 = t21 * t337;
t238 = g(3) * t253;
t223 = pkin(7) * t349;
t182 = pkin(3) * t349 + t223;
t241 = qJD(2) * qJ(3);
t155 = t241 + t182;
t114 = pkin(4) * t174 + qJD(5) + t155;
t61 = -pkin(5) * t288 + t114;
t424 = g(1) * t125 - g(2) * t127 - t247 * t2 - t215 * t238 - t61 * t54 + t20;
t11 = qJD(6) * t415 + t247 * t48 - t251 * t47;
t384 = t415 * t199;
t422 = -t11 + t384;
t227 = pkin(2) * t350;
t147 = qJD(1) * t291 + t227;
t163 = t252 * t182;
t334 = t252 * qJD(5);
t357 = qJ(5) - t255;
t371 = t248 * t249;
t421 = t341 * t357 - t334 + t248 * t147 - t163 - (pkin(4) * t253 - qJ(5) * t371) * qJD(1);
t185 = t357 * t252;
t318 = t252 * t350;
t354 = t252 * t147 + t248 * t182;
t420 = -qJ(5) * t318 - qJD(4) * t185 - t248 * qJD(5) - t354;
t286 = t244 * t252 + t245 * t248;
t355 = t207 * t286;
t124 = -t215 * t250 + t216 * t367;
t126 = t215 * t254 + t216 * t369;
t5 = t244 * t16 + t245 * t18;
t3 = pkin(9) * t47 + t5;
t320 = t251 * t2 - t247 * t3;
t418 = -g(1) * t124 - g(2) * t126 + t216 * t238 - t61 * t415 + t320;
t402 = pkin(3) + pkin(7);
t417 = pkin(9) * t109;
t313 = t245 * t340;
t373 = t244 * t248;
t356 = -t244 * t341 + t245 * t318 - t350 * t373 + t313;
t389 = -t244 * t420 + t245 * t421;
t388 = t244 * t421 + t245 * t420;
t411 = t355 * t251;
t235 = t253 * pkin(2);
t352 = t235 + t230;
t187 = -pkin(1) - t352;
t169 = -pkin(8) * t253 + t187;
t194 = t402 * t249;
t353 = t252 * t169 + t248 * t194;
t285 = -t245 * t252 + t373;
t108 = -t247 * t286 - t251 * t285;
t152 = t252 * t171;
t410 = -t207 * t341 + t152;
t409 = pkin(4) * t340 + t416;
t218 = pkin(4) * t252 + pkin(3);
t339 = qJD(4) * t253;
t314 = t248 * t339;
t345 = qJD(2) * t249;
t408 = -pkin(4) * t314 + (-pkin(7) - t218) * t345;
t407 = t249 * t344 + t314;
t294 = g(1) * t254 + g(2) * t250;
t406 = t103 * pkin(4) + qJDD(5);
t362 = t252 * t254;
t157 = -t248 * t250 + t249 * t362;
t365 = t250 * t252;
t159 = t248 * t254 + t249 * t365;
t363 = t252 * t253;
t405 = -g(1) * t157 - g(2) * t159 + g(3) * t363;
t404 = t155 * t207 + t255 * t171;
t62 = t244 * t69;
t25 = t245 * t60 - t62;
t321 = -g(1) * t367 - g(2) * t369 + t238;
t403 = -t25 * t355 + t26 * t356 - t285 * t4 + t286 * t5 + t321;
t401 = pkin(4) * t244;
t399 = g(1) * t250;
t395 = g(2) * t254;
t394 = g(3) * t249;
t231 = t248 * pkin(4);
t226 = pkin(2) * t345;
t122 = t226 + t264;
t343 = qJD(2) * t253;
t183 = t402 * t343;
t164 = t252 * t183;
t303 = qJ(5) * t253 - t169;
t37 = pkin(4) * t343 + t164 + t303 * t340 + (-qJ(5) * t345 - qJD(4) * t194 + qJD(5) * t253 - t122) * t248;
t271 = t252 * t122 - t169 * t341 + t248 * t183 + t194 * t340;
t41 = qJ(5) * t407 - t253 * t334 + t271;
t13 = t244 * t37 + t245 * t41;
t107 = -t247 * t285 + t251 * t286;
t391 = -qJD(6) * t107 - t247 * t356 - t411;
t390 = qJD(6) * t108 - t247 * t355 + t251 * t356;
t31 = t245 * t68 - t62;
t178 = t252 * t194;
t92 = t249 * pkin(4) + t248 * t303 + t178;
t97 = -qJ(5) * t363 + t353;
t43 = t244 * t92 + t245 * t97;
t19 = pkin(5) * t207 + t25 + t417;
t385 = t251 * t19;
t293 = qJD(1) * t218;
t383 = pkin(5) * t356 + t249 * t293 + t409;
t382 = pkin(7) * qJDD(2);
t380 = qJDD(2) * pkin(2);
t379 = t102 * t252;
t165 = qJDD(6) + t171;
t378 = t107 * t165;
t377 = t108 * t165;
t375 = t174 * t207;
t374 = t176 * t207;
t372 = t248 * t171;
t370 = t248 * t253;
t368 = t249 * t252;
t257 = qJD(1) ^ 2;
t366 = t249 * t257;
t364 = t250 * t253;
t246 = -qJ(5) - pkin(8);
t361 = t253 * t246;
t360 = t253 * t254;
t358 = qJ(3) + t231;
t184 = t357 * t248;
t116 = -t245 * t184 - t244 * t185;
t195 = t402 * t253;
t242 = t249 ^ 2;
t243 = t253 ^ 2;
t351 = t242 - t243;
t348 = qJD(2) * t174;
t347 = qJD(2) * t176;
t342 = qJD(4) * t136;
t338 = qJD(4) * t255;
t328 = pkin(4) * t370;
t325 = t248 * t367;
t324 = t253 * t366;
t220 = pkin(7) * t329;
t239 = qJDD(2) * qJ(3);
t240 = qJD(2) * qJD(3);
t323 = t220 + t239 + t240;
t322 = pkin(4) * t363 + t195;
t319 = t402 * qJD(2);
t317 = t248 * t345;
t307 = qJD(6) * t19 + t3;
t12 = -t244 * t41 + t245 * t37;
t30 = -t244 * t68 - t386;
t42 = -t244 * t97 + t245 * t92;
t304 = -qJD(2) * pkin(2) + qJD(3);
t115 = t184 * t244 - t245 * t185;
t301 = t254 * pkin(1) + pkin(2) * t360 + t250 * pkin(7) + qJ(3) * t367;
t300 = -t219 - t321;
t299 = pkin(3) * t329 + t323;
t298 = qJD(1) * t319;
t88 = -pkin(9) * t286 + t116;
t297 = pkin(5) * t349 - pkin(9) * t355 + qJD(6) * t88 - t389;
t87 = pkin(9) * t285 + t115;
t296 = pkin(9) * t356 - qJD(6) * t87 - t388;
t256 = qJD(2) ^ 2;
t295 = pkin(7) * t256 + t395;
t292 = qJD(6) * t285 - t356;
t7 = t247 * t19 + t251 * t21;
t144 = t285 * t253;
t145 = t286 * t253;
t289 = t251 * t144 + t145 * t247;
t90 = t144 * t247 - t145 * t251;
t186 = t222 + t304;
t193 = -t223 - t241;
t287 = t186 * t253 + t193 * t249;
t284 = t207 ^ 2;
t282 = t308 - t235;
t217 = pkin(4) * t245 + pkin(5);
t280 = t217 * t247 + t251 * t401;
t279 = t217 * t251 - t247 * t401;
t278 = t294 * t253;
t277 = -0.2e1 * pkin(1) * t331 - t382;
t156 = t282 * qJD(1);
t276 = t156 * t350 + qJDD(3) - t300;
t275 = -t207 * t340 - t372;
t273 = -qJ(3) * t343 - t335;
t269 = 0.2e1 * qJDD(1) * pkin(1) - t295;
t267 = t382 + (-qJD(1) * t187 - t156) * qJD(2);
t265 = -t278 - t394;
t113 = -t249 * t298 + t299;
t263 = t113 + t265;
t104 = qJD(1) * t273 + qJDD(1) * t282 + t206;
t149 = t226 + t273;
t259 = qJD(1) * t149 + qJDD(1) * t187 + t104 + t295;
t137 = pkin(7) * t311 - t323;
t146 = t309 - t380;
t258 = qJD(2) * t287 - t137 * t253 + t146 * t249;
t58 = t113 + t406;
t236 = t254 * pkin(7);
t211 = g(1) * t364;
t204 = qJ(3) * t360;
t202 = qJ(3) * t364;
t181 = t249 * t319;
t179 = -qJ(3) * t349 + t227;
t160 = -t248 * t369 + t362;
t158 = t325 + t365;
t135 = pkin(5) * t286 + t358;
t105 = -pkin(5) * t144 + t322;
t94 = -t244 * t407 - t245 * t317 + t253 * t313;
t93 = qJD(4) * t145 - t285 * t345;
t74 = pkin(4) * t176 - pkin(5) * t109;
t59 = -t93 * pkin(5) + t408;
t34 = pkin(9) * t144 + t43;
t33 = pkin(5) * t249 + pkin(9) * t145 + t42;
t28 = qJD(6) * t90 - t247 * t94 - t251 * t93;
t27 = qJD(6) * t289 + t247 * t93 - t251 * t94;
t24 = -t47 * pkin(5) + t58;
t23 = t31 + t417;
t22 = t30 - t414;
t9 = pkin(9) * t93 + t13;
t8 = pkin(5) * t343 + pkin(9) * t94 + t12;
t6 = -t21 * t247 + t385;
t1 = [qJDD(1), -t395 + t399, t294, qJDD(1) * t242 + 0.2e1 * t249 * t310, 0.2e1 * t249 * t329 - 0.2e1 * t331 * t351, qJDD(2) * t249 + t253 * t256, qJDD(2) * t253 - t249 * t256, 0, t249 * t277 + t253 * t269 + t211, t277 * t253 + (-t269 - t399) * t249 (t242 + t243) * qJDD(1) * pkin(7) + t258 - t294, t249 * t267 + t253 * t259 - t211, t267 * t253 + (-t259 + t399) * t249, pkin(7) * t258 - g(1) * t236 - g(2) * t301 + t104 * t187 + t156 * t149 - t282 * t399, -t102 * t370 + (-t252 * t339 + t317) * t176 (-t174 * t248 + t176 * t252) * t345 + (-t379 + t103 * t248 + (t174 * t252 + t176 * t248) * qJD(4)) * t253 (t207 * t346 + t102) * t249 + (t275 + t347) * t253 (t207 * t344 - t103) * t249 + (-t348 - t410) * t253, t171 * t249 + t207 * t343 (-t248 * t122 + t164) * t207 + (-t169 * t248 + t178) * t171 + t305 * t249 - t181 * t174 + t195 * t103 + t113 * t363 - g(1) * t160 - g(2) * t158 + (-t155 * t368 + t253 * t79) * qJD(2) + (-t155 * t370 - t207 * t353 - t80 * t249) * qJD(4), -t271 * t207 - t353 * t171 - t181 * t176 + t195 * t102 + g(1) * t159 - g(2) * t157 + ((qJD(2) * t155 + t342) * t248 + t326) * t249 + (-qJD(2) * t80 - t113 * t248 - t155 * t340) * t253, -g(2) * t360 + t109 * t12 + t13 * t288 + t144 * t5 + t145 * t4 + t25 * t94 + t26 * t93 - t42 * t48 + t43 * t47 + t211, t5 * t43 + t26 * t13 + t4 * t42 + t25 * t12 + t58 * t322 - g(1) * (t254 * t218 + t236) - g(2) * (pkin(4) * t325 - t246 * t360 + t301) + (-g(1) * (-pkin(4) * t371 + t282 + t361) - g(2) * t218) * t250 + t408 * t114, t10 * t90 + t27 * t415, t10 * t289 - t11 * t90 + t27 * t54 - t28 * t415, t10 * t249 + t165 * t90 + t199 * t27 + t343 * t415, -t11 * t249 + t165 * t289 - t199 * t28 + t343 * t54, t165 * t249 + t199 * t343 (-t247 * t9 + t251 * t8) * t199 + (-t247 * t34 + t251 * t33) * t165 + t320 * t249 + t6 * t343 - t59 * t54 + t105 * t11 - t24 * t289 + t61 * t28 - g(1) * t127 - g(2) * t125 + ((-t247 * t33 - t251 * t34) * t199 - t7 * t249) * qJD(6), -t7 * t343 + g(1) * t126 - g(2) * t124 + t105 * t10 + t20 * t249 + t24 * t90 + t61 * t27 + t59 * t415 + (-(-qJD(6) * t34 + t8) * t199 - t33 * t165 - t2 * t249) * t247 + (-(qJD(6) * t33 + t9) * t199 - t34 * t165 - t307 * t249) * t251; 0, 0, 0, -t324, t351 * t257, t330, t329, qJDD(2), pkin(1) * t366 + t300, t394 - t220 + (pkin(1) * t257 + t294) * t253 (-pkin(2) * t249 + t381) * qJDD(1) + ((-t193 - t241) * t249 + (-t186 + t304) * t253) * qJD(1), -t179 * t349 + t276 - 0.2e1 * t380, t220 + 0.2e1 * t239 + 0.2e1 * t240 + (qJD(1) * t179 - g(3)) * t249 + (qJD(1) * t156 - t294) * t253, -t137 * qJ(3) - t193 * qJD(3) - t146 * pkin(2) - t156 * t179 - g(1) * (-pkin(2) * t367 + t204) - g(2) * (-pkin(2) * t369 + t202) - g(3) * t352 - t287 * qJD(1) * pkin(7), -t248 * t374 + t379 (-t103 - t374) * t252 + (-t102 + t375) * t248 (-t176 * t253 - t207 * t371) * qJD(1) + t410 (t174 * t253 - t207 * t368) * qJD(1) + t275, -t207 * t349, -t79 * t349 + qJ(3) * t103 - t163 * t207 + t333 * t174 + t404 * t252 + ((t147 - t338) * t207 + t263) * t248, qJ(3) * t102 + t354 * t207 + t80 * t349 + t333 * t176 - t404 * t248 + (-t207 * t338 + t263) * t252, t109 * t389 - t115 * t48 + t116 * t47 + t288 * t388 - t403, t5 * t116 + t4 * t115 + t58 * t358 - g(1) * (t254 * t328 + t204) - g(2) * (t250 * t328 + t202) - g(3) * (t352 - t361) + t388 * t26 + t389 * t25 + t409 * t114 + (-g(3) * t231 + t114 * t293 + t294 * (pkin(2) - t246)) * t249, t10 * t108 + t391 * t415, -t10 * t107 - t108 * t11 - t390 * t415 + t391 * t54, t199 * t391 - t349 * t415 + t377, -t199 * t390 - t349 * t54 - t378, -t199 * t349 (-t247 * t88 + t251 * t87) * t165 + t135 * t11 + t24 * t107 - t6 * t349 + t390 * t61 - t383 * t54 + (t247 * t296 - t251 * t297) * t199 + t265 * t215 -(t247 * t87 + t251 * t88) * t165 + t135 * t10 + t24 * t108 + t7 * t349 + t391 * t61 + t383 * t415 + (t247 * t297 + t251 * t296) * t199 + t265 * t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t330, qJDD(2) + t324, -t242 * t257 - t256, qJD(2) * t193 + t205 + t276 - t380, 0, 0, 0, 0, 0, -t248 * t284 + t152 - t348, -t252 * t284 - t347 - t372, -t109 * t355 + t285 * t48 + t286 * t47 + t288 * t356, -t114 * qJD(2) + t403, 0, 0, 0, 0, 0, t377 + qJD(2) * t54 + (t247 * t292 - t286 * t336 - t411) * t199, -t378 - qJD(2) * t415 + (t292 * t251 + (qJD(6) * t286 + t355) * t247) * t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176 * t174, -t174 ^ 2 + t176 ^ 2, t102 + t375, -t103 + t374, t171, -t155 * t176 + t80 * t207 + t261 + t405, g(1) * t158 - g(2) * t160 + t155 * t174 + t79 * t207 + (t342 - t238) * t248 + t326 (t244 * t47 - t245 * t48) * pkin(4) + (-t31 + t25) * t288 + (-t26 - t30) * t109, -t25 * t30 - t26 * t31 + (-t114 * t176 + t5 * t244 + t4 * t245 + t405) * pkin(4), -t426, t425, t427, t422, t165, t279 * t165 - (t22 * t251 - t23 * t247) * t199 + t74 * t54 + (-t199 * t280 - t7) * qJD(6) + t418, -t280 * t165 - t251 * t3 + (t22 * t247 + t23 * t251) * t199 - t74 * t415 + (-t199 * t279 - t385) * qJD(6) + t424; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109 ^ 2 - t288 ^ 2, -t25 * t109 - t26 * t288 - t278 + (-g(3) - t298) * t249 + t299 + t406, 0, 0, 0, 0, 0, t11 + t384, t10 + t387; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, t425, t427, t422, t165 (-qJD(6) + t199) * t7 + t418, t6 * t199 - t251 * t307 + t424;];
tau_reg  = t1;
