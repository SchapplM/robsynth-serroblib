% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR11_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR11_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:49
% EndTime: 2019-12-31 22:44:22
% DurationCPUTime: 15.79s
% Computational Cost: add. (17740->766), mult. (43785->1063), div. (0->0), fcn. (34529->14), ass. (0->315)
t298 = sin(qJ(2));
t450 = cos(pkin(5));
t368 = t298 * t450;
t282 = pkin(1) * t368;
t297 = sin(qJ(3));
t301 = cos(qJ(3));
t347 = pkin(3) * t297 - pkin(9) * t301;
t294 = sin(pkin(5));
t302 = cos(qJ(2));
t427 = t294 * t302;
t502 = -(t282 + (pkin(7) + t347) * t427) * qJD(1) + t347 * qJD(3);
t300 = cos(qJ(4));
t409 = qJD(1) * t302;
t386 = t294 * t409;
t410 = qJD(1) * t294;
t387 = t298 * t410;
t296 = sin(qJ(4));
t422 = t296 * t301;
t188 = -t300 * t387 + t386 * t422;
t405 = qJD(3) * t301;
t494 = -t296 * t405 + t188;
t403 = qJD(4) * t300;
t495 = t297 * t403 - t494;
t366 = t450 * qJD(1);
t358 = pkin(1) * t366;
t212 = -pkin(7) * t387 + t302 * t358;
t349 = pkin(2) * t298 - pkin(8) * t302;
t213 = t349 * t410;
t146 = t301 * t212 + t297 * t213;
t125 = pkin(9) * t387 + t146;
t348 = pkin(3) * t301 + pkin(9) * t297;
t250 = -pkin(2) - t348;
t404 = qJD(4) * t296;
t406 = qJD(3) * t297;
t454 = -t300 * t125 + t250 * t403 + (-t300 * t406 - t301 * t404) * pkin(8) + t502 * t296;
t468 = pkin(8) * t296;
t501 = t125 * t296 + t502 * t300 + t406 * t468;
t419 = t301 * t302;
t423 = t296 * t298;
t189 = (t300 * t419 + t423) * t410;
t420 = t300 * t301;
t283 = pkin(8) * t420;
t385 = t297 * t409;
t500 = -pkin(4) * t294 * t385 + pkin(10) * t189 + (pkin(4) * t297 - pkin(10) * t420) * qJD(3) + (-t283 + (pkin(10) * t297 - t250) * t296) * qJD(4) + t501;
t499 = t495 * pkin(10) - t454;
t339 = t366 + qJD(2);
t199 = t297 * t339 + t301 * t387;
t257 = -qJD(3) + t386;
t150 = t296 * t199 + t257 * t300;
t152 = t199 * t300 - t257 * t296;
t295 = sin(qJ(5));
t470 = cos(qJ(5));
t329 = -t295 * t150 + t470 * t152;
t69 = t470 * t150 + t152 * t295;
t463 = t69 * t329;
t303 = -pkin(10) - pkin(9);
t390 = qJD(4) * t303;
t197 = t297 * t387 - t301 * t339;
t441 = t197 * t296;
t413 = pkin(7) * t427 + t282;
t209 = t450 * pkin(8) + t413;
t179 = qJD(2) * pkin(8) + t209 * qJD(1);
t337 = -pkin(2) * t302 - pkin(8) * t298 - pkin(1);
t190 = t337 * t410;
t107 = -t297 * t179 + t190 * t301;
t138 = pkin(3) * t199 + pkin(9) * t197;
t59 = t300 * t107 + t296 * t138;
t498 = -pkin(10) * t441 + t296 * t390 - t59;
t58 = -t107 * t296 + t300 * t138;
t497 = -pkin(4) * t199 - t58 + (-pkin(10) * t197 + t390) * t300;
t377 = t470 * qJD(5);
t496 = (t470 * qJD(4) + t377) * t300;
t191 = qJD(4) + t197;
t400 = qJDD(1) * t302;
t272 = t294 * t400;
t401 = qJD(1) * qJD(2);
t375 = t298 * t401;
t355 = t294 * t375;
t211 = qJDD(3) - t272 + t355;
t335 = qJD(2) * t358;
t363 = t450 * qJDD(1);
t352 = pkin(1) * t363;
t392 = pkin(7) * t272 + t298 * t352 + t302 * t335;
t153 = -pkin(7) * t355 + t392;
t333 = t363 + qJDD(2);
t134 = pkin(8) * t333 + t153;
t330 = qJD(2) * t349;
t141 = (qJD(1) * t330 + qJDD(1) * t337) * t294;
t44 = t301 * t134 + t297 * t141 - t179 * t406 + t190 * t405;
t36 = pkin(9) * t211 + t44;
t321 = qJD(3) * t339;
t407 = qJD(2) * t302;
t382 = t301 * t407;
t398 = t298 * qJDD(1);
t109 = -(qJD(1) * (-t298 * t406 + t382) + t301 * t398) * t294 - t297 * t333 - t301 * t321;
t383 = t297 * t407;
t110 = (qJD(1) * (t298 * t405 + t383) + t297 * t398) * t294 + t297 * t321 - t301 * t333;
t373 = t294 * t398;
t374 = t302 * t401;
t360 = t298 * t335 - t302 * t352 + (t294 * t374 + t373) * pkin(7);
t135 = -pkin(2) * t333 + t360;
t41 = t110 * pkin(3) + t109 * pkin(9) + t135;
t178 = -pkin(2) * t339 - t212;
t86 = t197 * pkin(3) - t199 * pkin(9) + t178;
t108 = t179 * t301 + t190 * t297;
t89 = -pkin(9) * t257 + t108;
t331 = -t296 * t41 - t300 * t36 - t86 * t403 + t404 * t89;
t46 = -t296 * t89 + t300 * t86;
t493 = -t46 * t191 - t331;
t343 = t300 * t405 - t189;
t492 = -t297 * t404 + t343;
t361 = t297 * t134 - t301 * t141 + t179 * t405 + t190 * t406;
t491 = t108 * t257 + t361;
t367 = t302 * t450;
t430 = t294 * t298;
t232 = pkin(1) * t367 - pkin(7) * t430;
t216 = qJD(2) * t232;
t490 = t329 ^ 2 - t69 ^ 2;
t397 = -t296 * t109 + t199 * t403 - t257 * t404;
t334 = t300 * t211 - t397;
t402 = qJD(5) * t295;
t54 = t300 * t109 + t199 * t404 - t296 * t211 + t257 * t403;
t16 = t150 * t377 + t152 * t402 - t295 * t334 + t470 * t54;
t184 = qJD(5) + t191;
t489 = t184 * t69 - t16;
t33 = -pkin(10) * t152 + t46;
t29 = pkin(4) * t191 + t33;
t47 = t296 * t86 + t300 * t89;
t34 = -pkin(10) * t150 + t47;
t102 = qJDD(4) + t110;
t9 = -t47 * qJD(4) - t296 * t36 + t300 * t41;
t6 = pkin(4) * t102 + pkin(10) * t54 + t9;
t7 = pkin(10) * t334 - t331;
t1 = t29 * t377 + t295 * t6 - t34 * t402 + t470 * t7;
t299 = sin(qJ(1));
t471 = cos(qJ(1));
t229 = -t299 * t368 + t471 * t302;
t429 = t294 * t299;
t173 = t229 * t301 + t297 * t429;
t228 = t471 * t298 + t299 * t367;
t293 = qJ(4) + qJ(5);
t287 = sin(t293);
t288 = cos(t293);
t119 = t173 * t288 + t228 * t287;
t428 = t294 * t301;
t225 = t297 * t450 + t298 * t428;
t351 = t450 * t471;
t227 = t298 * t351 + t299 * t302;
t389 = t294 * t471;
t169 = t227 * t301 - t297 * t389;
t226 = t298 * t299 - t302 * t351;
t486 = t169 * t288 + t226 * t287;
t88 = pkin(3) * t257 - t107;
t60 = pkin(4) * t150 + t88;
t488 = t60 * t69 + g(1) * t119 + g(2) * t486 - g(3) * (-t225 * t288 + t287 * t427) - t1;
t487 = t169 * t287 - t226 * t288;
t485 = t169 * t296 - t226 * t300;
t484 = t169 * t300 + t226 * t296;
t290 = t294 ^ 2;
t482 = 0.2e1 * t290;
t481 = -t47 * t191 - t9;
t480 = t107 * t257 + t44;
t208 = -pkin(2) * t450 - t232;
t224 = t297 * t430 - t301 * t450;
t121 = t224 * pkin(3) - t225 * pkin(9) + t208;
t414 = pkin(2) * t427 + pkin(8) * t430;
t210 = -pkin(1) * t294 - t414;
t140 = t301 * t209 + t297 * t210;
t123 = -pkin(9) * t427 + t140;
t57 = t296 * t121 + t300 * t123;
t207 = t296 * t250 + t283;
t479 = (qJDD(2) + 0.2e1 * t363) * t294;
t217 = t413 * qJD(2);
t477 = qJD(4) + qJD(5);
t126 = -t173 * t296 + t228 * t300;
t166 = t225 * t296 + t300 * t427;
t476 = -g(1) * t126 + g(2) * t485 + g(3) * t166;
t118 = -t173 * t287 + t228 * t288;
t396 = t470 * t34;
t12 = t295 * t29 + t396;
t2 = -qJD(5) * t12 - t295 * t7 + t470 * t6;
t475 = -t60 * t329 - g(1) * t118 + g(2) * t487 - g(3) * (-t225 * t287 - t288 * t427) + t2;
t17 = qJD(5) * t329 - t295 * t54 - t470 * t334;
t474 = t184 * t329 - t17;
t304 = qJD(1) ^ 2;
t464 = t211 * pkin(3);
t238 = t300 * t250;
t421 = t297 * t300;
t159 = -pkin(10) * t421 + t238 + (-pkin(4) - t468) * t301;
t424 = t296 * t297;
t176 = -pkin(10) * t424 + t207;
t90 = t470 * t159 - t295 * t176;
t462 = t90 * qJD(5) + t500 * t295 - t499 * t470;
t91 = t295 * t159 + t470 * t176;
t461 = -t91 * qJD(5) + t499 * t295 + t500 * t470;
t460 = t295 * t34;
t457 = t54 * t296;
t258 = t303 * t296;
t259 = t303 * t300;
t181 = t470 * t258 + t295 * t259;
t456 = t181 * qJD(5) + t497 * t295 + t470 * t498;
t182 = t295 * t258 - t470 * t259;
t455 = -t182 * qJD(5) - t295 * t498 + t497 * t470;
t453 = -t207 * qJD(4) + t501;
t145 = -t297 * t212 + t213 * t301;
t124 = -pkin(3) * t387 - t145;
t452 = pkin(4) * t495 + pkin(8) * t405 - t124;
t240 = t295 * t300 + t470 * t296;
t162 = t477 * t240;
t357 = t470 * t405;
t451 = -t162 * t297 - t470 * t189 + t295 * t494 + t300 * t357;
t449 = t102 * t300;
t446 = t150 * t191;
t445 = t150 * t296;
t444 = t152 * t150;
t443 = t152 * t191;
t442 = t197 * t257;
t440 = t199 * t197;
t439 = t199 * t257;
t434 = t257 * t297;
t433 = t287 * t301;
t432 = t288 * t301;
t431 = t290 * t304;
t426 = t295 * t296;
t425 = t296 * t102;
t418 = -t470 * t188 + t295 * t492 + t296 * t357 + t297 * t496 - t402 * t424;
t417 = -t240 * t197 - t162;
t328 = t470 * t300 - t426;
t416 = -t328 * t197 + t477 * t426 - t496;
t412 = t471 * pkin(1) + pkin(7) * t429;
t291 = t298 ^ 2;
t292 = t302 ^ 2;
t411 = t291 - t292;
t408 = qJD(2) * t298;
t395 = t302 * t431;
t394 = t296 * t427;
t393 = t229 * pkin(2) + t412;
t391 = pkin(4) * t296 + pkin(8);
t384 = t294 * t408;
t379 = g(3) * t414;
t376 = pkin(1) * t482;
t370 = -pkin(1) * t299 + pkin(7) * t389;
t56 = t300 * t121 - t123 * t296;
t139 = -t297 * t209 + t210 * t301;
t365 = -t227 * t297 - t301 * t389;
t364 = t191 * t300;
t362 = t298 * t395;
t356 = t298 * t374;
t354 = -t108 + (t404 + t441) * pkin(4);
t353 = -t227 * pkin(2) + t370;
t350 = t294 * t304 * t450;
t172 = t229 * t297 - t299 * t428;
t346 = g(1) * t365 + g(2) * t172;
t345 = -g(1) * t226 + g(2) * t228;
t344 = g(1) * t229 + g(2) * t227;
t122 = pkin(3) * t427 - t139;
t342 = -t47 * t296 - t46 * t300;
t286 = pkin(4) * t300 + pkin(3);
t341 = t286 * t301 - t297 * t303;
t338 = 0.2e1 * t366 + qJD(2);
t336 = pkin(8) * t228 + t393;
t214 = t294 * t330;
t67 = -t209 * t405 - t210 * t406 + t214 * t301 - t297 * t216;
t332 = g(1) * t471 + g(2) * t299;
t167 = t225 * t300 - t394;
t42 = pkin(4) * t224 - pkin(10) * t167 + t56;
t48 = -pkin(10) * t166 + t57;
t18 = -t295 * t48 + t470 * t42;
t19 = t295 * t42 + t470 * t48;
t93 = -t295 * t166 + t470 * t167;
t37 = t361 - t464;
t325 = -pkin(9) * t102 + t191 * t88;
t66 = -t209 * t406 + t210 * t405 + t297 * t214 + t301 * t216;
t62 = pkin(9) * t384 + t66;
t164 = qJD(3) * t225 + t294 * t383;
t165 = -qJD(3) * t224 + t294 * t382;
t77 = t164 * pkin(3) - t165 * pkin(9) + t217;
t23 = t121 * t403 - t123 * t404 + t296 * t77 + t300 * t62;
t324 = -pkin(8) * t226 + t353;
t323 = g(1) * t172 - g(2) * t365 + g(3) * t224;
t322 = -g(1) * t173 - g(2) * t169 - g(3) * t225;
t320 = t334 * t300;
t316 = t323 - t37;
t315 = -g(1) * t228 - g(2) * t226 + g(3) * t427;
t314 = -g(3) * t430 - t344;
t312 = -pkin(8) * t211 - t178 * t257;
t63 = -pkin(3) * t384 - t67;
t24 = -t57 * qJD(4) - t296 * t62 + t300 * t77;
t310 = qJD(4) * pkin(9) * t191 - t316;
t308 = pkin(8) * qJD(3) * t257 - t135 - t315;
t244 = t391 * t297;
t222 = t228 * pkin(2);
t220 = t226 * pkin(2);
t219 = t328 * t297;
t218 = t240 * t297;
t215 = t413 * qJD(1);
t206 = -pkin(8) * t422 + t238;
t127 = t173 * t300 + t228 * t296;
t98 = qJDD(5) + t102;
t92 = t470 * t166 + t167 * t295;
t85 = -qJD(4) * t166 + t165 * t300 + t296 * t384;
t84 = -qJD(4) * t394 + t165 * t296 + t225 * t403 - t300 * t384;
t73 = pkin(4) * t166 + t122;
t38 = pkin(4) * t84 + t63;
t28 = t93 * qJD(5) + t295 * t85 + t470 * t84;
t27 = t166 * t377 + t167 * t402 + t295 * t84 - t470 * t85;
t20 = -pkin(4) * t334 + t37;
t15 = -pkin(10) * t84 + t23;
t14 = t470 * t33 - t460;
t13 = -t295 * t33 - t396;
t11 = t470 * t29 - t460;
t10 = pkin(4) * t164 - pkin(10) * t85 + t24;
t4 = -t19 * qJD(5) + t470 * t10 - t295 * t15;
t3 = t18 * qJD(5) + t295 * t10 + t470 * t15;
t5 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t299 - g(2) * t471, t332, 0, 0, (qJDD(1) * t291 + 0.2e1 * t356) * t290, (t302 * t398 - t401 * t411) * t482, t294 * t338 * t407 + t479 * t298, (qJDD(1) * t292 - 0.2e1 * t356) * t290, t479 * t302 - t338 * t384, t333 * t450, -t217 * t339 + t232 * t333 - t360 * t450 + g(1) * t227 - g(2) * t229 + (-t375 + t400) * t376, -t216 * t339 - t413 * t333 - t153 * t450 + (-t374 - t398) * t376 + t345, ((-t212 * qJD(2) + qJDD(1) * t413 + t153) * t302 + (-qJD(2) * t215 - qJDD(1) * t232 + t360) * t298 - t332) * t294, t290 * qJDD(1) * pkin(1) ^ 2 - g(1) * t370 - g(2) * t412 + t153 * t413 - t212 * t217 + t215 * t216 - t232 * t360, -t109 * t225 + t165 * t199, t109 * t224 - t110 * t225 - t164 * t199 - t165 * t197, -t165 * t257 + t211 * t225 + (t109 * t302 + t199 * t408) * t294, t110 * t224 + t164 * t197, t164 * t257 - t211 * t224 + (t110 * t302 - t197 * t408) * t294, (-t211 * t302 - t257 * t408) * t294, g(1) * t169 - g(2) * t173 + t110 * t208 + t135 * t224 + t139 * t211 + t164 * t178 + t197 * t217 - t257 * t67 + (t107 * t408 + t302 * t361) * t294, -t109 * t208 + t135 * t225 - t140 * t211 + t165 * t178 + t199 * t217 + t257 * t66 + (-t108 * t408 + t302 * t44) * t294 + t346, -t107 * t165 - t108 * t164 + t109 * t139 - t110 * t140 - t197 * t66 - t199 * t67 - t224 * t44 + t225 * t361 - t345, -g(1) * t324 - g(2) * t336 + t107 * t67 + t108 * t66 + t135 * t208 - t139 * t361 + t44 * t140 + t178 * t217, t152 * t85 - t167 * t54, -t85 * t150 - t152 * t84 + t54 * t166 + t167 * t334, t102 * t167 + t152 * t164 + t191 * t85 - t224 * t54, t150 * t84 - t166 * t334, -t166 * t102 - t150 * t164 - t84 * t191 + t224 * t334, t102 * t224 + t164 * t191, g(1) * t484 - g(2) * t127 + t56 * t102 - t122 * t334 + t63 * t150 + t46 * t164 + t37 * t166 + t24 * t191 + t9 * t224 + t88 * t84, -g(1) * t485 - g(2) * t126 - t57 * t102 - t122 * t54 + t63 * t152 - t47 * t164 + t37 * t167 - t23 * t191 + t331 * t224 + t88 * t85, -t23 * t150 - t24 * t152 + t166 * t331 - t9 * t167 + t334 * t57 - t46 * t85 - t47 * t84 + t56 * t54 - t346, -t331 * t57 + t47 * t23 + t9 * t56 + t46 * t24 + t37 * t122 + t88 * t63 - g(1) * (-pkin(3) * t169 + pkin(9) * t365 + t324) - g(2) * (pkin(3) * t173 + pkin(9) * t172 + t336), -t16 * t93 - t27 * t329, t16 * t92 - t17 * t93 + t27 * t69 - t28 * t329, -t16 * t224 + t164 * t329 - t184 * t27 + t93 * t98, t17 * t92 + t28 * t69, -t164 * t69 - t17 * t224 - t184 * t28 - t92 * t98, t164 * t184 + t224 * t98, g(1) * t486 - g(2) * t119 + t11 * t164 + t73 * t17 + t18 * t98 + t4 * t184 + t2 * t224 + t20 * t92 + t60 * t28 + t38 * t69, -g(1) * t487 - g(2) * t118 - t1 * t224 - t12 * t164 - t73 * t16 - t3 * t184 - t19 * t98 + t20 * t93 - t60 * t27 + t38 * t329, -t1 * t92 + t11 * t27 - t12 * t28 + t16 * t18 - t17 * t19 - t2 * t93 - t3 * t69 - t329 * t4 - t346, t1 * t19 + t12 * t3 + t2 * t18 + t11 * t4 + t20 * t73 + t60 * t38 - g(1) * (-t169 * t286 - t226 * t391 - t303 * t365 + t353) - g(2) * (-t172 * t303 + t173 * t286 + t228 * t391 + t393); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t362, t411 * t431, -t302 * t350 + t373, t362, t298 * t350 + t272, t333, pkin(1) * t298 * t431 + t215 * t339 - t315 - t360, pkin(1) * t395 + t212 * t339 + (pkin(7) * t401 + g(3)) * t430 + t344 - t392, 0, 0, -t109 * t297 - t301 * t439, (-t109 + t442) * t301 + (-t110 + t439) * t297, -t257 * t405 + t211 * t297 + (-t199 * t298 + t257 * t419) * t410, -t110 * t301 - t197 * t434, t257 * t406 + t211 * t301 + (t197 * t298 - t302 * t434) * t410, t257 * t387, -pkin(2) * t110 - t107 * t387 + t145 * t257 - t197 * t215 + t297 * t312 + t301 * t308, pkin(2) * t109 + t108 * t387 - t146 * t257 - t199 * t215 - t297 * t308 + t301 * t312, t145 * t199 + t146 * t197 + ((qJD(3) * t199 - t110) * pkin(8) + t480) * t301 + ((qJD(3) * t197 - t109) * pkin(8) + t491) * t297 + t314, -t135 * pkin(2) - t108 * t146 - t107 * t145 - t178 * t215 + g(1) * t222 + g(2) * t220 - t379 + (t361 * t297 + t44 * t301 + (-t107 * t301 - t108 * t297) * qJD(3) - t344) * pkin(8), t152 * t492 - t54 * t421, t189 * t150 + t152 * t188 + (-t150 * t300 - t152 * t296) * t405 + (t320 + t457 + (-t152 * t300 + t445) * qJD(4)) * t297, t301 * t54 + t343 * t191 + (-t152 * t257 - t191 * t404 + t449) * t297, t150 * t495 - t334 * t424, -t334 * t301 + t494 * t191 + (t150 * t257 - t191 * t403 - t425) * t297, -t102 * t301 - t191 * t434, t206 * t102 - t124 * t150 - t88 * t188 + t453 * t191 + t314 * t296 + (-t9 + (pkin(8) * t150 + t88 * t296) * qJD(3) - t315 * t300) * t301 + (-pkin(8) * t334 - t257 * t46 + t37 * t296 + t403 * t88) * t297, -t207 * t102 - t124 * t152 - t88 * t189 - t454 * t191 + t314 * t300 + (-t331 + (pkin(8) * t152 + t300 * t88) * qJD(3) + t315 * t296) * t301 + (-pkin(8) * t54 + t257 * t47 + t37 * t300 - t404 * t88) * t297, t207 * t334 + t206 * t54 + t47 * t188 + t46 * t189 - t453 * t152 - t454 * t150 + t342 * t405 + (t331 * t296 - t9 * t300 + (t296 * t46 - t300 * t47) * qJD(4) - t315) * t297, -t331 * t207 + t9 * t206 - t88 * t124 - g(1) * (-t348 * t228 - t222) - g(2) * (-t348 * t226 - t220) - g(3) * (t348 * t427 + t414) + t454 * t47 + t453 * t46 + (t37 * t297 + t405 * t88 - t344) * pkin(8), -t16 * t219 + t329 * t451, t16 * t218 - t17 * t219 - t329 * t418 - t451 * t69, t16 * t301 + t184 * t451 + t219 * t98 - t329 * t434, t17 * t218 + t418 * t69, t17 * t301 - t184 * t418 - t218 * t98 + t434 * t69, -t184 * t434 - t301 * t98, t244 * t17 + t20 * t218 + t90 * t98 - t2 * t301 + t11 * t406 - g(1) * (-t228 * t432 + t229 * t287) - g(2) * (-t226 * t432 + t227 * t287) + t452 * t69 + t418 * t60 + t461 * t184 + (-t11 * t385 - g(3) * (t287 * t298 + t288 * t419)) * t294, -t244 * t16 + t20 * t219 - t91 * t98 + t1 * t301 - t12 * t406 - g(1) * (t228 * t433 + t229 * t288) - g(2) * (t226 * t433 + t227 * t288) + t452 * t329 + t451 * t60 - t462 * t184 + (t12 * t385 - g(3) * (-t287 * t419 + t288 * t298)) * t294, -t1 * t218 - t451 * t11 - t418 * t12 + t16 * t90 - t17 * t91 - t2 * t219 - t315 * t297 - t329 * t461 - t462 * t69, t1 * t91 + t2 * t90 + t20 * t244 - g(1) * (-t228 * t341 + t229 * t391 - t222) - g(2) * (-t226 * t341 + t227 * t391 - t220) - t379 + t452 * t60 - g(3) * (pkin(4) * t423 + t302 * t341) * t294 + t462 * t12 + t461 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t440, -t197 ^ 2 + t199 ^ 2, -t109 - t442, -t440, -t110 - t439, t211, -t178 * t199 + t323 - t491, t178 * t197 - t322 - t480, 0, 0, t152 * t364 - t457, (-t54 - t446) * t300 + (t334 - t443) * t296, -t152 * t199 + t191 * t364 + t425, t191 * t445 + t320, -t191 ^ 2 * t296 + t150 * t199 + t449, -t191 * t199, -pkin(3) * t397 - t58 * t191 - t46 * t199 - t108 * t150 + t325 * t296 + (-t310 + t464) * t300, pkin(3) * t54 - t108 * t152 + t191 * t59 + t199 * t47 + t296 * t310 + t300 * t325, t59 * t150 + t58 * t152 + ((qJD(4) * t152 + t334) * pkin(9) + t493) * t300 + ((qJD(4) * t150 - t54) * pkin(9) + t481) * t296 + t322, -t88 * t108 - t46 * t58 - t47 * t59 + t316 * pkin(3) + (qJD(4) * t342 - t9 * t296 - t300 * t331 + t322) * pkin(9), -t16 * t240 - t329 * t416, -t16 * t328 - t17 * t240 + t329 * t417 + t416 * t69, -t184 * t416 - t199 * t329 + t240 * t98, -t17 * t328 - t417 * t69, t184 * t417 + t199 * t69 + t328 * t98, -t184 * t199, -t11 * t199 - t17 * t286 + t181 * t98 + t455 * t184 - t20 * t328 + t323 * t288 + t354 * t69 - t417 * t60, t12 * t199 + t16 * t286 - t182 * t98 - t456 * t184 + t20 * t240 - t323 * t287 + t329 * t354 - t416 * t60, t1 * t328 + t416 * t11 + t417 * t12 + t16 * t181 - t17 * t182 - t2 * t240 - t329 * t455 - t456 * t69 + t322, t1 * t182 + t2 * t181 - t20 * t286 - g(1) * (-t172 * t286 - t173 * t303) - g(2) * (-t169 * t303 + t286 * t365) - g(3) * (-t224 * t286 - t225 * t303) + t354 * t60 + t456 * t12 + t455 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t444, -t150 ^ 2 + t152 ^ 2, -t54 + t446, -t444, t334 + t443, t102, -t88 * t152 + t476 - t481, g(1) * t127 + g(2) * t484 + g(3) * t167 + t88 * t150 - t493, 0, 0, t463, t490, t489, -t463, t474, t98, -t13 * t184 + (-t152 * t69 - t184 * t402 + t470 * t98) * pkin(4) + t475, t14 * t184 + (-t152 * t329 - t184 * t377 - t295 * t98) * pkin(4) + t488, -t11 * t69 + t12 * t329 + t13 * t329 + t14 * t69 + (t470 * t16 - t17 * t295 + (t295 * t329 - t470 * t69) * qJD(5)) * pkin(4), -t11 * t13 - t12 * t14 + (t1 * t295 + t2 * t470 - t60 * t152 + (-t11 * t295 + t12 * t470) * qJD(5) + t476) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t463, t490, t489, -t463, t474, t98, t12 * t184 + t475, t11 * t184 + t488, 0, 0;];
tau_reg = t5;
