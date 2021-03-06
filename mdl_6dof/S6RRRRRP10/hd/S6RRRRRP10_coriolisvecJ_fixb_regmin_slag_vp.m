% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRP10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:27:26
% EndTime: 2019-03-10 02:28:01
% DurationCPUTime: 13.88s
% Computational Cost: add. (18470->635), mult. (47008->854), div. (0->0), fcn. (37015->10), ass. (0->279)
t266 = sin(pkin(6));
t273 = cos(qJ(2));
t367 = qJD(1) * t273;
t351 = t266 * t367;
t240 = -qJD(3) + t351;
t268 = sin(qJ(4));
t271 = cos(qJ(4));
t272 = cos(qJ(3));
t324 = t272 * t351;
t270 = sin(qJ(2));
t368 = qJD(1) * t266;
t352 = t270 * t368;
t185 = t268 * t324 - t271 * t352;
t364 = qJD(3) * t272;
t449 = -t268 * t364 + t185;
t402 = cos(pkin(6));
t332 = t402 * qJD(1);
t323 = pkin(1) * t332;
t245 = t270 * t323;
t215 = pkin(8) * t351 + t245;
t269 = sin(qJ(3));
t453 = t215 + t240 * (pkin(3) * t269 - pkin(10) * t272);
t362 = qJD(4) * t271;
t452 = t269 * t362 - t449;
t417 = cos(qJ(5));
t339 = t417 * qJD(5);
t451 = t271 * (t417 * qJD(4) + t339);
t365 = qJD(3) * t269;
t416 = pkin(9) * t268;
t450 = -t271 * t453 + t365 * t416;
t212 = -pkin(8) * t352 + t273 * t323;
t302 = (pkin(2) * t270 - pkin(9) * t273) * t266;
t213 = qJD(1) * t302;
t372 = t272 * t212 + t269 * t213;
t135 = pkin(10) * t352 + t372;
t236 = -pkin(3) * t272 - pkin(10) * t269 - pkin(2);
t363 = qJD(4) * t268;
t448 = -t236 * t362 - (-t271 * t365 - t272 * t363) * pkin(9) + t271 * t135 + t453 * t268;
t383 = t272 * t273;
t186 = (t268 * t270 + t271 * t383) * t368;
t424 = -t269 * t363 + t271 * t364 - t186;
t384 = t271 * t272;
t255 = pkin(9) * t384;
t325 = t269 * t351;
t447 = -pkin(4) * t325 + pkin(11) * t186 + t135 * t268 + (pkin(4) * t269 - pkin(11) * t384) * qJD(3) + (-t255 + (pkin(11) * t269 - t236) * t268) * qJD(4) + t450;
t446 = pkin(11) * t452 + t448;
t418 = -pkin(11) - pkin(10);
t353 = qJD(4) * t418;
t311 = t332 + qJD(2);
t174 = pkin(9) * t311 + t215;
t211 = (-pkin(2) * t273 - pkin(9) * t270 - pkin(1)) * t266;
t191 = qJD(1) * t211;
t122 = -t269 * t174 + t191 * t272;
t197 = t269 * t352 - t272 * t311;
t199 = t269 * t311 + t272 * t352;
t139 = pkin(3) * t199 + pkin(10) * t197;
t381 = t271 * t122 + t268 * t139;
t396 = t197 * t268;
t445 = pkin(11) * t396 - t268 * t353 + t381;
t138 = t271 * t139;
t395 = t197 * t271;
t444 = -pkin(4) * t199 - pkin(11) * t395 + t122 * t268 + t271 * t353 - t138;
t392 = t266 * t270;
t221 = t269 * t402 + t272 * t392;
t391 = t266 * t273;
t300 = -t221 * t271 + t268 * t391;
t443 = pkin(11) * t300;
t147 = -t268 * t199 - t240 * t271;
t148 = t199 * t271 - t240 * t268;
t267 = sin(qJ(5));
t296 = t267 * t147 + t148 * t417;
t108 = pkin(3) * t240 - t122;
t74 = -pkin(4) * t147 + t108;
t86 = -t417 * t147 + t148 * t267;
t32 = pkin(5) * t86 - qJ(6) * t296 + t74;
t442 = t32 * t86;
t441 = t74 * t86;
t414 = t296 * t86;
t229 = t267 * t271 + t268 * t417;
t425 = qJD(4) + qJD(5);
t161 = t425 * t229;
t322 = t417 * t364;
t380 = t161 * t269 + t186 * t417 - t449 * t267 - t271 * t322;
t361 = qJD(5) * t267;
t387 = t268 * t269;
t379 = -t185 * t417 + t424 * t267 + t268 * t322 + t269 * t451 - t361 * t387;
t390 = t267 * t268;
t293 = t271 * t417 - t390;
t375 = t293 * t197 - t390 * t425 + t451;
t374 = t229 * t197 + t161;
t426 = t365 - t325;
t419 = t296 ^ 2;
t439 = -t86 ^ 2 + t419;
t192 = qJD(4) + t197;
t181 = qJD(5) + t192;
t358 = qJD(1) * qJD(2);
t337 = t266 * t358;
t320 = t273 * t337;
t278 = -qJD(3) * t197 + t272 * t320;
t321 = t270 * t337;
t277 = t199 * t363 + t240 * t362 - t268 * t321 - t271 * t278;
t283 = -t199 * t362 + t240 * t363 - t268 * t278 + t271 * t321;
t34 = -t147 * t339 + t148 * t361 - t267 * t283 + t417 * t277;
t21 = t181 * t86 - t34;
t54 = pkin(5) * t296 + qJ(6) * t86;
t438 = t192 ^ 2;
t123 = t272 * t174 + t269 * t191;
t319 = -t123 + (t363 + t396) * pkin(4);
t152 = t199 * qJD(3) + t269 * t320;
t388 = t268 * t152;
t436 = t270 * t273;
t241 = t418 * t268;
t242 = t418 * t271;
t294 = t241 * t417 + t267 * t242;
t435 = -qJD(5) * t294 - t444 * t267 + t445 * t417;
t179 = t267 * t241 - t242 * t417;
t434 = -qJD(5) * t179 + t445 * t267 + t444 * t417;
t354 = pkin(1) * t402;
t428 = -pkin(8) * t392 + t273 * t354;
t209 = -pkin(2) * t402 - t428;
t220 = t269 * t392 - t272 * t402;
t130 = t220 * pkin(3) - t221 * pkin(10) + t209;
t284 = pkin(8) * t391 + t270 * t354;
t210 = pkin(9) * t402 + t284;
t373 = t272 * t210 + t269 * t211;
t132 = -pkin(10) * t391 + t373;
t378 = t268 * t130 + t271 * t132;
t227 = t271 * t236;
t386 = t269 * t271;
t157 = -pkin(11) * t386 + t227 + (-pkin(4) - t416) * t272;
t370 = t268 * t236 + t255;
t170 = -pkin(11) * t387 + t370;
t433 = t267 * t157 + t417 * t170;
t329 = -t269 * t212 + t213 * t272;
t134 = -pkin(3) * t352 - t329;
t432 = pkin(4) * t452 + pkin(9) * t364 - t134;
t360 = t148 * qJD(3);
t342 = t272 * t360;
t431 = t269 * t277 - t342;
t430 = -t157 * t339 + t170 * t361 - t447 * t267 + t446 * t417;
t429 = t324 - t364;
t150 = t152 * pkin(5);
t173 = -pkin(2) * t311 - t212;
t103 = t197 * pkin(3) - t199 * pkin(10) + t173;
t109 = -pkin(10) * t240 + t123;
t214 = qJD(2) * t302;
t206 = qJD(1) * t214;
t216 = t428 * qJD(2);
t207 = qJD(1) * t216;
t289 = -t174 * t365 + t191 * t364 + t269 * t206 + t272 * t207;
t72 = pkin(10) * t321 + t289;
t314 = pkin(8) * t320;
t80 = t152 * pkin(3) - pkin(10) * t278 + qJD(2) * t245 + t314;
t335 = -t268 * t72 + t271 * t80;
t11 = t152 * pkin(4) + pkin(11) * t277 - t103 * t363 - t109 * t362 + t335;
t26 = t103 * t362 - t109 * t363 + t268 * t80 + t271 * t72;
t14 = pkin(11) * t283 + t26;
t62 = t271 * t103 - t109 * t268;
t51 = -pkin(11) * t148 + t62;
t44 = pkin(4) * t192 + t51;
t63 = t268 * t103 + t271 * t109;
t52 = pkin(11) * t147 + t63;
t333 = -t417 * t11 + t267 * t14 + t52 * t339 + t44 * t361;
t2 = -t150 + t333;
t286 = t296 * t32 + t2;
t35 = t147 * t361 + t148 * t339 - t267 * t277 - t417 * t283;
t423 = t181 * t296 - t35;
t422 = -t74 * t296 - t333;
t349 = qJD(2) * t391;
t165 = -qJD(3) * t220 + t272 * t349;
t166 = t221 * t268 + t271 * t391;
t366 = qJD(2) * t270;
t350 = t266 * t366;
t102 = -qJD(4) * t166 + t165 * t271 + t268 * t350;
t164 = qJD(3) * t221 + t269 * t349;
t288 = -t210 * t365 + t211 * t364 + t269 * t214 + t272 * t216;
t76 = pkin(10) * t350 + t288;
t217 = t284 * qJD(2);
t95 = t164 * pkin(3) - t165 * pkin(10) + t217;
t334 = -t268 * t76 + t271 * t95;
t18 = pkin(4) * t164 - pkin(11) * t102 - qJD(4) * t378 + t334;
t348 = t271 * t366;
t399 = t165 * t268;
t285 = t266 * t348 - t399;
t357 = t130 * t362 + t268 * t95 + t271 * t76;
t389 = t268 * t132;
t25 = t285 * pkin(11) + (-t389 + t443) * qJD(4) + t357;
t331 = t271 * t130 - t389;
t59 = pkin(4) * t220 + t331 + t443;
t64 = -pkin(11) * t166 + t378;
t303 = t267 * t59 + t417 * t64;
t421 = -qJD(5) * t303 + t18 * t417 - t267 * t25;
t420 = -qJD(5) * t433 + t446 * t267 + t447 * t417;
t274 = qJD(1) ^ 2;
t415 = pkin(9) * t272;
t413 = -qJ(6) * t426 + qJD(6) * t272 + t430;
t412 = pkin(5) * t426 + t420;
t411 = -pkin(5) * t374 + qJ(6) * t375 + qJD(6) * t229 - t319;
t219 = t293 * t269;
t410 = -pkin(5) * t379 - qJ(6) * t380 + qJD(6) * t219 - t432;
t355 = t417 * t52;
t20 = t267 * t44 + t355;
t407 = t181 * t20;
t406 = t267 * t52;
t405 = -qJ(6) * t199 - t435;
t404 = t199 * pkin(5) - t434;
t23 = t417 * t51 - t406;
t403 = -pkin(4) * t339 - qJD(6) + t23;
t401 = t148 * t192;
t400 = t152 * t272;
t398 = t294 * t152;
t397 = t179 * t152;
t394 = t199 * t240;
t290 = t240 * t269;
t263 = t266 ^ 2;
t393 = t263 * t274;
t385 = t271 * t152;
t19 = t417 * t44 - t406;
t382 = qJD(6) - t19;
t233 = pkin(4) * t387 + t269 * pkin(9);
t369 = t270 ^ 2 - t273 ^ 2;
t359 = t148 * qJD(4);
t260 = -pkin(4) * t271 - pkin(3);
t346 = t240 * t365;
t341 = t269 * t359;
t338 = t263 * t358;
t336 = -t267 * t11 - t417 * t14 - t44 * t339 + t52 * t361;
t330 = -t269 * t210 + t211 * t272;
t328 = 0.2e1 * t338;
t327 = t174 * t364 + t191 * t365 - t272 * t206 + t269 * t207;
t22 = t267 * t51 + t355;
t318 = pkin(4) * t361 - t22;
t317 = t266 * t274 * t402;
t316 = -0.2e1 * pkin(1) * t338;
t131 = pkin(3) * t391 - t330;
t310 = 0.2e1 * t332 + qJD(2);
t149 = t152 * qJ(6);
t175 = t181 * qJD(6);
t1 = t149 + t175 - t336;
t309 = -t210 * t364 - t211 * t365 + t272 * t214 - t269 * t216;
t308 = t181 * t19 + t336;
t307 = -t267 * t64 + t417 * t59;
t73 = -pkin(3) * t321 + t327;
t301 = t108 * t362 + t73 * t268;
t299 = t267 * t18 + t417 * t25 + t59 * t339 - t361 * t64;
t298 = t157 * t417 - t267 * t170;
t112 = -t267 * t166 - t300 * t417;
t292 = -t192 * t362 - t388;
t91 = pkin(4) * t166 + t131;
t282 = pkin(1) * (-qJD(2) * t332 + t393);
t27 = -qJD(4) * t63 + t335;
t280 = qJD(4) * t300 + t285;
t46 = -pkin(4) * t283 + t73;
t53 = -(-t221 * t362 - t399) * pkin(4) + (-pkin(3) * t366 - (t273 * t363 + t348) * pkin(4)) * t266 - t309;
t275 = t277 * t271;
t259 = -pkin(4) * t417 - pkin(5);
t257 = pkin(4) * t267 + qJ(6);
t218 = t229 * t269;
t208 = qJD(1) * t217;
t154 = -pkin(5) * t293 - qJ(6) * t229 + t260;
t136 = pkin(5) * t218 - qJ(6) * t219 + t233;
t133 = t152 * t220;
t111 = t166 * t417 - t267 * t300;
t106 = t272 * pkin(5) - t298;
t105 = -qJ(6) * t272 + t433;
t77 = -pkin(3) * t350 - t309;
t45 = pkin(4) * t148 + t54;
t42 = pkin(5) * t111 - qJ(6) * t112 + t91;
t38 = qJD(5) * t112 + t267 * t102 - t280 * t417;
t37 = -t102 * t417 + t166 * t339 - t267 * t280 - t300 * t361;
t29 = -t220 * pkin(5) - t307;
t28 = qJ(6) * t220 + t303;
t16 = t181 * qJ(6) + t20;
t15 = -t181 * pkin(5) + t382;
t8 = t38 * pkin(5) + t37 * qJ(6) - t112 * qJD(6) + t53;
t7 = t35 * pkin(5) + t34 * qJ(6) - qJD(6) * t296 + t46;
t6 = -t164 * pkin(5) - t421;
t5 = qJ(6) * t164 + qJD(6) * t220 + t299;
t3 = [0, 0, 0, t328 * t436, -t369 * t328, t310 * t349, -t310 * t350, 0, -t208 * t402 - t217 * t311 + t270 * t316, -t207 * t402 - t216 * t311 + t273 * t316, t199 * t165 + t221 * t278, -t221 * t152 - t199 * t164 - t165 * t197 - t220 * t278, -t165 * t240 + t199 * t350 + t221 * t321 - t278 * t391, t164 * t240 + (t152 * t273 + (-qJD(1) * t220 - t197) * t366) * t266 (-t240 * t266 - t263 * t367) * t366, -t309 * t240 + t217 * t197 + t209 * t152 + t208 * t220 + t173 * t164 + (t327 * t273 + (qJD(1) * t330 + t122) * t366) * t266, -t123 * t350 + t173 * t165 + t217 * t199 + t208 * t221 + t209 * t278 + t240 * t288 + t289 * t391 - t321 * t373, t148 * t102 + t277 * t300, t102 * t147 + t148 * t280 + t166 * t277 - t283 * t300, t102 * t192 + t148 * t164 - t152 * t300 - t220 * t277, t147 * t164 - t166 * t152 + t192 * t280 + t220 * t283, t164 * t192 + t133, t334 * t192 + t331 * t152 + t27 * t220 + t62 * t164 - t77 * t147 - t131 * t283 + t73 * t166 - t108 * t285 + (-t108 * t300 - t192 * t378) * qJD(4) -(-t132 * t363 + t357) * t192 - t378 * t152 - t26 * t220 - t63 * t164 + t77 * t148 - t131 * t277 - t73 * t300 + t108 * t102, -t112 * t34 - t296 * t37, t111 * t34 - t112 * t35 - t296 * t38 + t37 * t86, t112 * t152 + t164 * t296 - t181 * t37 - t220 * t34, -t111 * t152 - t164 * t86 - t181 * t38 - t220 * t35, t164 * t181 + t133, t46 * t111 + t307 * t152 + t19 * t164 + t181 * t421 - t220 * t333 + t91 * t35 + t74 * t38 + t53 * t86, t46 * t112 - t152 * t303 - t20 * t164 - t181 * t299 + t220 * t336 + t296 * t53 - t91 * t34 - t74 * t37, t111 * t7 - t15 * t164 - t152 * t29 - t181 * t6 - t2 * t220 + t32 * t38 + t35 * t42 + t8 * t86, -t1 * t111 + t112 * t2 - t15 * t37 - t16 * t38 - t28 * t35 - t29 * t34 + t296 * t6 - t5 * t86, t1 * t220 - t112 * t7 + t152 * t28 + t16 * t164 + t181 * t5 - t296 * t8 + t32 * t37 + t34 * t42, t1 * t28 + t15 * t6 + t16 * t5 + t2 * t29 + t32 * t8 + t42 * t7; 0, 0, 0, -t393 * t436, t369 * t393, -t273 * t317, t270 * t317, 0, t215 * t311 + t270 * t282 - t314, pkin(8) * t321 + t212 * t311 + t273 * t282, -qJD(3) * t269 ^ 2 * t352 + ((qJD(3) * t311 + t320) * t269 - t394) * t272, -t269 * t152 + t197 * t429 - t199 * t426 + t272 * t278, -t240 * t364 + (t240 * t383 + (qJD(2) * t269 - t199) * t270) * t368, t346 + (-t273 * t290 + (qJD(2) * t272 + t197) * t270) * t368, t240 * t352, -pkin(2) * t152 - t208 * t272 + t329 * t240 - t215 * t197 + (t173 * t269 + t240 * t415) * qJD(3) + (-t122 * t270 + (-pkin(9) * t366 - t173 * t273) * t269) * t368, -pkin(2) * t278 - pkin(9) * t346 + t123 * t352 - t173 * t429 - t215 * t199 + t208 * t269 - t240 * t372 - t321 * t415, -t148 * t186 - t268 * t341 - t269 * t275 + t271 * t342, t147 * t424 + t148 * t185 + t268 * t431 - t271 * t341 + t283 * t386, -t148 * t325 + t272 * t277 + (t360 + t385) * t269 + t424 * t192, -t283 * t272 + t449 * t192 + (-t147 * t240 + t292) * t269, -t192 * t290 - t400, -t108 * t185 + t134 * t147 + t227 * t152 + ((-qJD(4) * t236 + t135) * t268 + t450) * t192 + (t108 * t268 * qJD(3) - t27 + (-qJD(3) * t147 + t292) * pkin(9)) * t272 + (-pkin(9) * t283 - t240 * t62 + t301) * t269, -t431 * pkin(9) + t424 * t108 - t134 * t148 - t152 * t370 + t448 * t192 + t26 * t272 + t386 * t73 - t426 * t63, -t219 * t34 - t296 * t380, t218 * t34 - t219 * t35 - t296 * t379 + t380 * t86, t152 * t219 - t181 * t380 + t272 * t34 - t290 * t296, -t152 * t218 - t181 * t379 + t272 * t35 + t290 * t86, -t181 * t290 - t400, t298 * t152 + t181 * t420 + t19 * t426 + t46 * t218 + t233 * t35 + t272 * t333 + t379 * t74 + t432 * t86, -t152 * t433 + t181 * t430 + t20 * t290 + t46 * t219 - t233 * t34 - t272 * t336 + t296 * t432 - t380 * t74, -t106 * t152 + t136 * t35 + t15 * t290 + t181 * t412 + t2 * t272 + t218 * t7 + t32 * t379 - t410 * t86, -t1 * t218 - t105 * t35 - t106 * t34 - t15 * t380 - t16 * t379 + t2 * t219 - t296 * t412 + t413 * t86, -t1 * t272 + t105 * t152 + t136 * t34 - t16 * t290 - t181 * t413 - t219 * t7 + t296 * t410 + t32 * t380, t1 * t105 + t106 * t2 + t136 * t7 - t15 * t412 - t16 * t413 - t32 * t410; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199 * t197, -t197 ^ 2 + t199 ^ 2, -t197 * t240 + t278, -t152 - t394, t321, -t123 * t240 - t173 * t199 - t327, -t122 * t240 + t173 * t197 - t289 (-qJD(4) * t199 + t321) * t268 ^ 2 + ((-qJD(4) * t240 + t278) * t268 + t401) * t271, -t148 * t396 - t275 + (t283 - t359) * t268 + (t362 + t395) * t147, -t148 * t199 + t271 * t438 + t388, -t147 * t199 - t268 * t438 + t385, -t192 * t199, -pkin(10) * t388 + pkin(3) * t283 + t123 * t147 - t62 * t199 - t73 * t271 + (-pkin(10) * t362 - t138 + (t108 + t122) * t268) * t192, pkin(3) * t277 + t108 * t395 - t123 * t148 + t192 * t381 + t63 * t199 + t301 + (t192 * t363 - t385) * pkin(10), -t229 * t34 + t296 * t375, -t229 * t35 - t293 * t34 - t296 * t374 - t375 * t86, t152 * t229 + t181 * t375 - t199 * t296, t152 * t293 - t181 * t374 + t199 * t86, -t181 * t199, t181 * t434 - t19 * t199 + t260 * t35 - t293 * t46 + t319 * t86 + t374 * t74 + t398, t181 * t435 + t20 * t199 + t46 * t229 - t260 * t34 + t319 * t296 + t375 * t74 - t397, t15 * t199 + t154 * t35 - t181 * t404 - t293 * t7 + t32 * t374 - t411 * t86 + t398, t1 * t293 + t15 * t375 - t16 * t374 - t179 * t35 + t2 * t229 + t294 * t34 + t296 * t404 - t405 * t86, t154 * t34 - t16 * t199 + t181 * t405 - t229 * t7 + t296 * t411 - t32 * t375 + t397, t1 * t179 + t15 * t404 + t154 * t7 + t16 * t405 - t2 * t294 - t32 * t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148 * t147, -t147 ^ 2 + t148 ^ 2, -t147 * t192 - t277, t283 + t401, t152, -t108 * t148 + t192 * t63 + t27, -t108 * t147 + t192 * t62 - t26, t414, t439, t21, t423, t152, t22 * t181 + (-t148 * t86 + t152 * t417 - t181 * t361) * pkin(4) + t422, t23 * t181 + t441 + (-t148 * t296 - t152 * t267 - t181 * t339) * pkin(4) + t336, -t152 * t259 - t181 * t318 - t45 * t86 - t286, -t257 * t35 - t259 * t34 + (t16 + t318) * t296 + (t15 + t403) * t86, t152 * t257 - t181 * t403 + t296 * t45 + t1 - t442, t1 * t257 + t15 * t318 - t16 * t403 + t2 * t259 - t32 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t414, t439, t21, t423, t152, t407 + t422, t308 + t441, -t54 * t86 + t150 - t286 + t407, pkin(5) * t34 - qJ(6) * t35 + (t16 - t20) * t296 + (t15 - t382) * t86, t296 * t54 + 0.2e1 * t149 + 0.2e1 * t175 - t308 - t442, -pkin(5) * t2 + qJ(6) * t1 - t15 * t20 + t16 * t382 - t32 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152 + t414, t21, -t181 ^ 2 - t419, -t16 * t181 + t286;];
tauc_reg  = t3;
