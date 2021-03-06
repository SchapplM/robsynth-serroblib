% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:16:17
% EndTime: 2019-03-09 22:16:39
% DurationCPUTime: 9.67s
% Computational Cost: add. (15576->627), mult. (36991->771), div. (0->0), fcn. (26413->8), ass. (0->323)
t284 = sin(qJ(4));
t385 = qJD(4) * t284;
t285 = sin(qJ(3));
t286 = sin(qJ(2));
t387 = qJD(1) * t286;
t363 = t285 * t387;
t289 = cos(qJ(2));
t452 = cos(qJ(3));
t364 = t452 * t289;
t226 = -qJD(1) * t364 + t363;
t414 = t226 * t284;
t495 = t385 + t414;
t451 = pkin(2) * t285;
t272 = pkin(9) + t451;
t288 = cos(qJ(4));
t355 = t452 * qJD(3);
t343 = pkin(2) * t355;
t333 = t288 * t343;
t302 = -t272 * t385 + t333;
t402 = t285 * t289;
t244 = t286 * t452 + t402;
t228 = t244 * qJD(1);
t173 = pkin(3) * t228 + pkin(9) * t226;
t372 = pkin(2) * t387;
t154 = t173 + t372;
t454 = -pkin(8) - pkin(7);
t258 = t454 * t289;
t248 = qJD(1) * t258;
t229 = t285 * t248;
t256 = t454 * t286;
t246 = qJD(1) * t256;
t182 = t246 * t452 + t229;
t104 = t284 * t154 + t288 * t182;
t222 = t228 * qJ(5);
t86 = t222 + t104;
t494 = -t86 + t302;
t274 = -pkin(2) * t289 - pkin(1);
t254 = qJD(1) * t274;
t152 = t226 * pkin(3) - t228 * pkin(9) + t254;
t230 = t452 * t248;
t443 = qJD(2) * pkin(2);
t232 = t246 + t443;
t179 = t285 * t232 - t230;
t377 = qJD(2) + qJD(3);
t156 = pkin(9) * t377 + t179;
t95 = t288 * t152 - t284 * t156;
t396 = qJD(5) - t95;
t493 = t495 * pkin(10);
t283 = sin(qJ(6));
t287 = cos(qJ(6));
t224 = qJD(4) + t226;
t455 = pkin(4) + pkin(5);
t195 = t288 * t228 + t284 * t377;
t486 = -t195 * pkin(10) + t396;
t55 = -t224 * t455 + t486;
t215 = t224 * qJ(5);
t193 = t228 * t284 - t288 * t377;
t96 = t152 * t284 + t156 * t288;
t72 = pkin(10) * t193 + t96;
t66 = t215 + t72;
t24 = -t283 * t66 + t287 * t55;
t379 = qJD(6) - t224;
t492 = t24 * t379;
t132 = t193 * t283 + t195 * t287;
t178 = t452 * t232 + t229;
t155 = -pkin(3) * t377 - t178;
t307 = t195 * qJ(5) - t155;
t70 = -t193 * t455 + t307;
t491 = t70 * t132;
t490 = t493 + t494;
t175 = t284 * t182;
t447 = -pkin(10) + t272;
t240 = t447 * t288;
t334 = t284 * t343;
t369 = t455 * t228;
t449 = pkin(10) * t226;
t489 = -qJD(4) * t240 - t334 + t175 + (-t154 + t449) * t288 - t369;
t111 = t284 * t173 + t288 * t178;
t91 = t222 + t111;
t488 = -pkin(9) * t385 + t493 - t91;
t165 = t284 * t178;
t453 = pkin(9) - pkin(10);
t257 = t453 * t288;
t487 = qJD(4) * t257 - t165 - (-t173 + t449) * t288 + t369;
t329 = -t287 * t193 + t195 * t283;
t427 = t329 * t132;
t367 = qJD(2) * t454;
t342 = qJD(1) * t367;
t234 = t289 * t342;
t221 = t452 * t234;
t233 = t286 * t342;
t350 = t285 * t233 - t221;
t114 = qJD(3) * t179 + t350;
t344 = qJD(4) * t377;
t356 = t452 * qJD(2);
t304 = t289 * (t356 + t355);
t378 = qJD(1) * qJD(2);
t354 = t286 * t378;
t389 = qJD(3) * t363 + t285 * t354;
t479 = -qJD(1) * t304 + t389;
t122 = t228 * t385 + (t479 - t344) * t288;
t485 = t122 * qJ(5) - t195 * qJD(5) + t114;
t243 = -t283 * t288 + t284 * t287;
t381 = qJD(6) * t287;
t382 = qJD(6) * t283;
t384 = qJD(4) * t288;
t463 = -t283 * t384 - t284 * t381 + t287 * t385 + t288 * t382;
t392 = t243 * t226 + t463;
t241 = t283 * t284 + t287 * t288;
t185 = t241 * qJD(6) - t283 * t385 - t287 * t384;
t484 = t241 * t226 - t185;
t181 = t246 * t285 - t230;
t386 = qJD(3) * t285;
t339 = pkin(2) * t386 - t181;
t483 = t132 ^ 2 - t329 ^ 2;
t25 = t283 * t55 + t287 * t66;
t188 = t377 * t244;
t170 = t188 * qJD(1);
t370 = t286 * t443;
t100 = t389 * pkin(9) + t170 * pkin(3) + (-pkin(9) * t304 + t370) * qJD(1);
t113 = t232 * t355 + t233 * t452 + t285 * t234 + t248 * t386;
t345 = -t288 * t100 + t284 * t113 + t152 * t385 + t156 * t384;
t12 = t122 * pkin(10) - t170 * t455 + t345;
t296 = t284 * t479;
t123 = qJD(4) * t195 - t296;
t163 = t170 * qJ(5);
t213 = t224 * qJD(5);
t318 = -t284 * t100 - t288 * t113 - t152 * t384 + t156 * t385;
t21 = t163 + t213 - t318;
t14 = pkin(10) * t123 + t21;
t4 = -qJD(6) * t25 + t287 * t12 - t283 * t14;
t482 = t25 * t379 + t4;
t351 = -t283 * t122 - t287 * t123;
t49 = qJD(6) * t132 + t351;
t481 = t132 * t379 - t49;
t48 = t287 * t122 - t283 * t123 - t193 * t381 + t195 * t382;
t461 = t329 * t379 - t48;
t341 = t289 * t356;
t403 = t285 * t286;
t187 = -t289 * t355 + t377 * t403 - t341;
t417 = t195 * t288;
t420 = t193 * t284;
t428 = t123 * t288;
t429 = t122 * t284;
t480 = t244 * (qJD(4) * (-t417 + t420) - t428 + t429) + (t193 * t288 + t195 * t284) * t187;
t418 = t195 * t224;
t422 = t193 * t224;
t26 = t284 * (t123 + t418) + t288 * (t122 + t422);
t477 = -0.2e1 * t378;
t400 = t288 * t170;
t476 = pkin(9) * (t224 * t385 - t400);
t450 = pkin(4) * t170;
t29 = t345 - t450;
t475 = t21 * t288 + t29 * t284;
t413 = t226 * t288;
t80 = t215 + t96;
t88 = t193 * pkin(4) - t307;
t474 = -t80 * t228 - t88 * t413;
t41 = t123 * pkin(4) + t485;
t473 = -t41 * t288 + t88 * t385;
t390 = pkin(4) * t414 - qJ(5) * t413;
t470 = -t390 - t339;
t223 = pkin(4) * t385 - qJ(5) * t384 - t284 * qJD(5);
t467 = -pkin(5) * t495 - t223;
t466 = t452 * t256 + t285 * t258;
t278 = t284 * qJ(5);
t465 = t288 * pkin(4) + t278;
t3 = t283 * t12 + t287 * t14 + t55 * t381 - t382 * t66;
t462 = t329 * t70 - t3;
t242 = -t364 + t403;
t425 = t187 * t284;
t321 = t244 * t384 - t425;
t406 = t284 * t170;
t460 = t242 * t123 + t188 * t193 + t224 * t321 + t244 * t406;
t432 = qJ(5) * t288;
t459 = t284 * t455 - t432;
t200 = t285 * t256 - t258 * t452;
t247 = t286 * t367;
t137 = t200 * qJD(3) + t285 * t247 - t454 * t341;
t457 = t195 ^ 2;
t456 = t224 ^ 2;
t448 = t228 * pkin(4);
t239 = t447 * t284;
t172 = t239 * t283 + t240 * t287;
t446 = qJD(6) * t172 + t283 * t490 + t287 * t489;
t171 = t239 * t287 - t240 * t283;
t445 = -qJD(6) * t171 + t283 * t489 - t287 * t490;
t444 = pkin(2) * qJD(3);
t442 = t195 * t88;
t31 = t318 * t288;
t441 = t41 * t284;
t255 = t453 * t284;
t197 = t255 * t287 - t257 * t283;
t438 = qJD(6) * t197 + t283 * t487 + t287 * t488;
t199 = t255 * t283 + t257 * t287;
t437 = -qJD(6) * t199 - t283 * t488 + t287 * t487;
t436 = t467 + t470;
t112 = t179 - t390;
t435 = t112 + t467;
t251 = -qJ(5) * t283 - t287 * t455;
t434 = qJD(6) * t251 - t283 * t72 + t287 * t486;
t252 = qJ(5) * t287 - t283 * t455;
t433 = -qJD(6) * t252 - t283 * t486 - t287 * t72;
t431 = t114 * t466;
t141 = t170 * t242;
t426 = t170 * t272;
t424 = t187 * t288;
t423 = t193 * qJ(5);
t421 = t193 * t272;
t419 = t195 * t193;
t416 = t379 * t228;
t415 = t224 * t228;
t412 = t228 * t226;
t411 = t244 * t284;
t410 = t244 * t288;
t409 = t254 * t228;
t292 = qJD(1) ^ 2;
t399 = t289 * t292;
t291 = qJD(2) ^ 2;
t398 = t291 * t286;
t397 = t291 * t289;
t395 = -t193 * t333 - t272 * t428;
t394 = t112 - t223;
t393 = -t223 + t470;
t177 = pkin(3) * t242 - pkin(9) * t244 + t274;
t125 = t284 * t177 + t288 * t200;
t388 = t286 ^ 2 - t289 ^ 2;
t383 = qJD(5) * t288;
t376 = -t95 * t413 - t96 * t414 - t31;
t374 = pkin(3) + t465;
t373 = t452 * pkin(2);
t368 = t286 * t399;
t101 = t242 * qJ(5) + t125;
t366 = t284 * t452;
t365 = t288 * t452;
t361 = t272 * t384;
t357 = t193 ^ 2 - t457;
t352 = pkin(1) * t477;
t103 = t288 * t154 - t175;
t110 = t288 * t173 - t165;
t190 = t284 * t200;
t124 = t288 * t177 - t190;
t348 = t379 ^ 2;
t347 = t224 * t284;
t346 = t224 * t288;
t273 = -t373 - pkin(3);
t340 = t289 * t354;
t338 = t114 * t284 + t155 * t384 + t96 * t228;
t337 = pkin(4) * t284 - t432;
t77 = t190 + (-pkin(10) * t244 - t177) * t288 - t455 * t242;
t84 = pkin(10) * t411 + t101;
t46 = -t283 * t84 + t287 * t77;
t47 = t283 * t77 + t287 * t84;
t79 = -pkin(4) * t224 + t396;
t336 = t284 * t80 - t288 * t79;
t335 = t284 * t96 + t288 * t95;
t331 = t155 * t226 - t426;
t330 = -t228 * t193 - t400;
t326 = t417 + t420;
t121 = pkin(3) * t188 + pkin(9) * t187 + t370;
t136 = qJD(3) * t466 + t452 * t247 + t367 * t402;
t51 = t288 * t121 - t284 * t136 - t177 * t385 - t200 * t384;
t325 = t79 * t228 + t473;
t324 = t384 * t88 + t441;
t323 = -t114 * t288 + t155 * t385 - t95 * t228;
t236 = t273 - t465;
t320 = t244 * t385 + t424;
t319 = t224 * t96 - t345;
t22 = -t123 * t455 - t485;
t317 = t22 * t241 + t228 * t24 - t392 * t70;
t316 = t22 * t243 - t228 * t25 + t484 * t70;
t50 = t284 * t121 + t288 * t136 + t177 * t384 - t200 * t385;
t314 = -t80 * t414 + t475 + (t384 + t413) * t79;
t312 = (-t224 * t384 - t406) * pkin(9);
t310 = t193 * t347 - t428;
t309 = -t24 * t484 - t3 * t241 - t4 * t243 + t25 * t392;
t308 = -t122 * t272 + t195 * t343;
t36 = t188 * qJ(5) + t242 * qJD(5) + t50;
t305 = t224 * t95 + t318;
t303 = -t334 - t361;
t301 = t123 * t411 + t193 * t321;
t300 = -qJD(4) * t336 + t475;
t299 = -qJD(4) * t335 + t284 * t345 - t31;
t295 = t254 * t226 - t113;
t293 = t228 * t384 + t284 * t344 - t296 - t418;
t279 = t288 * pkin(5);
t237 = t279 + t374;
t214 = t279 - t236;
t162 = t241 * t244;
t161 = t243 * t244;
t142 = -t226 ^ 2 + t228 ^ 2;
t139 = t226 * t377 - t479;
t138 = pkin(4) * t195 + t423;
t133 = t244 * t337 - t466;
t116 = pkin(9) * t428;
t106 = -t244 * t459 + t466;
t102 = -t242 * pkin(4) - t124;
t99 = -t195 * t455 - t423;
t93 = -t110 - t448;
t87 = -t103 - t448;
t85 = t188 * t224 + t141;
t78 = -t122 + t422;
t65 = -t195 * t228 + t224 * t346 + t406;
t64 = -t284 * t456 - t330;
t63 = t224 * t347 + t330;
t60 = t195 * t346 - t429;
t59 = t185 * t244 + t187 * t243;
t58 = t187 * t241 + t244 * t463;
t53 = -t337 * t187 + (qJD(4) * t465 - t383) * t244 + t137;
t52 = -t122 * t410 - t195 * t320;
t45 = t459 * t187 + (t383 + (-t288 * t455 - t278) * qJD(4)) * t244 - t137;
t44 = -t188 * pkin(4) - t51;
t43 = t241 * t170 - t228 * t329 + t379 * t392;
t42 = t228 * t132 - t170 * t243 + t379 * t484;
t30 = -t122 * t242 + t195 * t188 - t224 * t320 + t244 * t400;
t23 = pkin(10) * t321 + t36;
t19 = pkin(10) * t320 - t188 * t455 - t51;
t9 = t49 * t241 - t329 * t392;
t8 = t132 * t484 - t48 * t243;
t7 = t132 * t392 + t241 * t48 - t49 * t243 - t329 * t484;
t6 = -qJD(6) * t47 + t287 * t19 - t283 * t23;
t5 = qJD(6) * t46 + t283 * t19 + t287 * t23;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t340, t388 * t477, t397, -0.2e1 * t340, -t398, 0, -pkin(7) * t397 + t286 * t352, pkin(7) * t398 + t289 * t352, 0, 0, -t228 * t187 - t244 * t479, -t244 * t170 + t187 * t226 - t228 * t188 + t242 * t479, -t187 * t377, t188 * t226 + t141, -t188 * t377, 0, t274 * t170 + t254 * t188 - t137 * t377 + (qJD(1) * t242 + t226) * t370, pkin(2) * t244 * t354 - t136 * t377 - t254 * t187 + t228 * t370 - t274 * t479, -t113 * t242 + t114 * t244 - t136 * t226 + t137 * t228 - t200 * t170 + t178 * t187 - t179 * t188 + t466 * t479, t113 * t200 + t179 * t136 - t178 * t137 + 0.2e1 * t254 * t370 - t431, t52, t480, t30, t301, -t460, t85, t114 * t411 - t123 * t466 + t124 * t170 + t137 * t193 + t155 * t321 + t188 * t95 + t224 * t51 - t242 * t345, t114 * t410 + t122 * t466 - t125 * t170 + t137 * t195 - t155 * t320 - t188 * t96 - t224 * t50 + t242 * t318, t122 * t124 - t123 * t125 - t193 * t50 - t195 * t51 + t335 * t187 + (t284 * t318 + t288 * t345 + (t284 * t95 - t288 * t96) * qJD(4)) * t244, -t124 * t345 - t125 * t318 + t137 * t155 + t50 * t96 + t51 * t95 - t431, t52, t30, -t480, t85, t460, t301, -t102 * t170 + t123 * t133 - t188 * t79 + t193 * t53 - t224 * t44 - t242 * t29 + t244 * t324 - t425 * t88, -t101 * t123 - t102 * t122 - t193 * t36 + t195 * t44 + t336 * t187 + (-t21 * t284 + t288 * t29 + (-t284 * t79 - t288 * t80) * qJD(4)) * t244, t101 * t170 + t122 * t133 + t188 * t80 - t195 * t53 + t21 * t242 + t224 * t36 + t244 * t473 + t88 * t424, t101 * t21 + t102 * t29 + t133 * t41 + t36 * t80 + t44 * t79 + t53 * t88, -t132 * t58 - t162 * t48, -t132 * t59 - t161 * t48 - t162 * t49 + t329 * t58, -t132 * t188 - t162 * t170 + t242 * t48 - t379 * t58, -t161 * t49 + t329 * t59, -t161 * t170 + t188 * t329 + t242 * t49 - t379 * t59, -t188 * t379 + t141, t106 * t49 - t161 * t22 - t170 * t46 - t188 * t24 - t242 * t4 + t329 * t45 + t379 * t6 + t59 * t70, -t106 * t48 + t132 * t45 + t162 * t22 + t170 * t47 + t188 * t25 + t242 * t3 - t379 * t5 - t58 * t70, -t132 * t6 + t161 * t3 - t162 * t4 + t24 * t58 - t25 * t59 - t329 * t5 + t46 * t48 - t47 * t49, t106 * t22 + t24 * t6 + t25 * t5 + t3 * t47 + t4 * t46 + t45 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t368, t388 * t292, 0, t368, 0, 0, t292 * pkin(1) * t286, pkin(1) * t399, 0, 0, t412, t142, t139, -t412, 0, 0, t248 * t355 + t221 - t226 * t372 - t409 + t181 * t377 + (-qJD(3) * t232 - t377 * t444 - t233) * t285, t182 * t377 + (-t228 * t387 - t355 * t377) * pkin(2) + t295, -t170 * t451 + t479 * t373 + (t179 + t339) * t228 + (-t178 + t182 - t343) * t226, t178 * t181 - t179 * t182 + (-t254 * t387 - t452 * t114 + t113 * t285 + (-t178 * t285 + t179 * t452) * qJD(3)) * pkin(2), t60, -t26, t65, t310, t64, -t415, t273 * t123 + t331 * t284 + t339 * t193 + (-t103 + t303) * t224 + t323, -t273 * t122 + t331 * t288 + t339 * t195 + (t104 - t302) * t224 + t338, t103 * t195 + t104 * t193 + (t195 * t272 - t95) * t384 + (t345 + (-t96 + t421) * qJD(4) + t308) * t284 + t376 + t395, -t95 * t103 - t96 * t104 + t114 * t273 - t155 * t181 + (t155 * t285 + t365 * t96 - t366 * t95) * t444 + t299 * t272, t60, t65, t26, -t415, t63, t310, t236 * t123 + (t226 * t88 - t426) * t284 - t393 * t193 + (t87 + t303) * t224 + t325, t86 * t193 + (-t87 + t361) * t195 + ((-t80 + t421) * qJD(4) + t308) * t284 + t314 + t395, t236 * t122 - t441 + (-qJD(4) * t88 + t426) * t288 + t393 * t195 + t494 * t224 + t474, t41 * t236 - t79 * t87 - t80 * t86 - t393 * t88 + (t365 * t80 + t366 * t79) * t444 + t300 * t272, t8, t7, t42, t9, t43, t416, -t170 * t171 + t214 * t49 + t329 * t436 - t379 * t446 + t317, t132 * t436 + t170 * t172 - t214 * t48 + t379 * t445 + t316, t132 * t446 + t171 * t48 - t172 * t49 + t329 * t445 + t309, t171 * t4 + t172 * t3 + t214 * t22 - t24 * t446 - t25 * t445 + t436 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t412, t142, t139, -t412, 0, 0, t179 * qJD(2) - t350 - t409, t178 * t377 + t295, 0, 0, t60, -t26, t65, t310, t64, -t415, -pkin(3) * t123 - t110 * t224 + t155 * t414 - t179 * t193 + t312 + t323, pkin(3) * t122 + t111 * t224 + t155 * t413 - t179 * t195 + t338 + t476, t110 * t195 + t111 * t193 - t116 + (-pkin(9) * t122 + t345) * t284 + (pkin(9) * t326 - t335) * qJD(4) + t376, -t114 * pkin(3) + pkin(9) * t299 - t95 * t110 - t96 * t111 - t155 * t179, t60, t65, t26, -t415, t63, t310, -t123 * t374 - t193 * t394 + t224 * t93 + t414 * t88 + t312 + t325, -t80 * t385 + t193 * t91 - t195 * t93 - t116 + (qJD(4) * t326 - t429) * pkin(9) + t314, -t122 * t374 + t195 * t394 - t224 * t91 - t324 + t474 - t476, pkin(9) * t300 - t374 * t41 - t394 * t88 - t79 * t93 - t80 * t91, t8, t7, t42, t9, t43, t416, -t170 * t197 + t237 * t49 + t329 * t435 + t379 * t437 + t317, t132 * t435 + t170 * t199 - t237 * t48 - t379 * t438 + t316, -t132 * t437 + t197 * t48 - t199 * t49 - t329 * t438 + t309, t197 * t4 + t199 * t3 + t22 * t237 + t24 * t437 + t25 * t438 + t435 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t419, -t357, t78, -t419, -t293, t170, -t155 * t195 + t319, t155 * t193 + t305, 0, 0, t419, t78, t357, t170, t293, -t419, -t138 * t193 + t319 - t442 + 0.2e1 * t450, pkin(4) * t122 - t123 * qJ(5) + (t80 - t96) * t195 + (t79 - t396) * t193, t138 * t195 - t193 * t88 + 0.2e1 * t163 + 0.2e1 * t213 - t305, -t29 * pkin(4) + t21 * qJ(5) - t88 * t138 + t396 * t80 - t79 * t96, -t427, -t483, -t461, t427, -t481, t170, -t170 * t251 - t329 * t99 + t379 * t433 - t4 + t491, -t132 * t99 + t170 * t252 - t379 * t434 - t462, t251 * t48 - t252 * t49 + (t24 - t434) * t329 + (-t25 - t433) * t132, t24 * t433 + t25 * t434 + t4 * t251 + t3 * t252 - t70 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t228 * t377 + t419, t78, -t456 - t457, -t224 * t80 + t29 + t442, 0, 0, 0, 0, 0, 0, -t287 * t170 - t195 * t329 - t283 * t348, -t195 * t132 + t283 * t170 - t287 * t348, t283 * t481 - t287 * t461, -t70 * t195 + t482 * t287 + (t3 - t492) * t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t427, t483, t461, -t427, -t351 + (-qJD(6) + t379) * t132, -t170, t482 - t491, t462 + t492, 0, 0;];
tauc_reg  = t1;
