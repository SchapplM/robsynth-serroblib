% Calculate time derivative of joint inertia matrix for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:20:19
% EndTime: 2019-03-09 02:20:39
% DurationCPUTime: 11.55s
% Computational Cost: add. (48372->789), mult. (38206->1105), div. (0->0), fcn. (36774->12), ass. (0->398)
t294 = pkin(11) + qJ(4);
t287 = sin(t294);
t503 = -qJD(1) * t287 / 0.2e1;
t289 = cos(t294);
t297 = qJ(5) + qJ(6);
t291 = sin(t297);
t292 = cos(t297);
t339 = Icges(7,5) * t292 - Icges(7,6) * t291;
t295 = qJD(5) + qJD(6);
t447 = t287 * t295;
t171 = (-Icges(7,5) * t291 - Icges(7,6) * t292) * t447 + (Icges(7,3) * t287 + t289 * t339) * qJD(4);
t466 = Icges(7,4) * t292;
t342 = -Icges(7,2) * t291 + t466;
t230 = -Icges(7,6) * t289 + t287 * t342;
t418 = qJD(4) * t291;
t502 = -t230 * t418 - t171;
t301 = sin(qJ(5));
t303 = cos(qJ(5));
t340 = Icges(6,5) * t303 - Icges(6,6) * t301;
t414 = qJD(5) * t287;
t202 = (-Icges(6,5) * t301 - Icges(6,6) * t303) * t414 + (Icges(6,3) * t287 + t289 * t340) * qJD(4);
t468 = Icges(6,4) * t303;
t343 = -Icges(6,2) * t301 + t468;
t238 = -Icges(6,6) * t289 + t287 * t343;
t416 = qJD(4) * t301;
t501 = -t238 * t416 - t202;
t296 = qJ(1) + pkin(10);
t290 = cos(t296);
t420 = qJD(4) * t289;
t382 = t420 / 0.2e1;
t288 = sin(t296);
t425 = qJD(1) * t288;
t393 = t287 * t425;
t500 = t290 * t382 - t393 / 0.2e1;
t423 = qJD(1) * t290;
t383 = t423 / 0.2e1;
t499 = t287 * t383 + t288 * t382;
t388 = t288 * t420;
t318 = t287 * t423 + t388;
t492 = 2 * m(5);
t476 = rSges(5,1) * t289;
t364 = -rSges(5,2) * t287 + t476;
t226 = -rSges(5,3) * t290 + t288 * t364;
t442 = t289 * t290;
t284 = t288 * rSges(5,3);
t448 = t287 * t290;
t495 = -rSges(5,2) * t448 + t284;
t227 = rSges(5,1) * t442 + t495;
t264 = rSges(5,1) * t287 + rSges(5,2) * t289;
t325 = qJD(4) * t264;
t309 = rSges(5,2) * t393 + rSges(5,3) * t423 - t290 * t325;
t93 = (qJD(1) * t226 + t309) * t290 + (-t288 * t325 + (-t227 + t495) * qJD(1)) * t288;
t498 = t492 * t93;
t305 = -pkin(9) - pkin(8);
t477 = pkin(8) + t305;
t497 = t289 * t477;
t471 = Icges(5,4) * t287;
t349 = Icges(5,1) * t289 - t471;
t225 = Icges(5,5) * t288 + t290 * t349;
t451 = t225 * t289;
t470 = Icges(5,4) * t289;
t345 = -Icges(5,2) * t287 + t470;
t223 = Icges(5,6) * t288 + t290 * t345;
t456 = t223 * t287;
t329 = -t451 + t456;
t496 = t329 * t290;
t300 = -pkin(7) - qJ(3);
t299 = cos(pkin(11));
t285 = pkin(3) * t299 + pkin(2);
t326 = -t285 - t364;
t480 = sin(qJ(1)) * pkin(1);
t200 = -t480 + (rSges(5,3) - t300) * t290 + t326 * t288;
t293 = cos(qJ(1)) * pkin(1);
t367 = t290 * t285 - t288 * t300 + t293;
t201 = t227 + t367;
t494 = t200 * t290 + t201 * t288;
t341 = Icges(5,5) * t289 - Icges(5,6) * t287;
t220 = -Icges(5,3) * t290 + t288 * t341;
t222 = -Icges(5,6) * t290 + t288 * t345;
t224 = -Icges(5,5) * t290 + t288 * t349;
t441 = t290 * t291;
t445 = t288 * t292;
t243 = -t289 * t441 + t445;
t440 = t290 * t292;
t446 = t288 * t291;
t244 = t289 * t440 + t446;
t168 = t244 * rSges(7,1) + t243 * rSges(7,2) + rSges(7,3) * t448;
t276 = pkin(4) * t442;
t248 = pkin(8) * t448 + t276;
t444 = t288 * t301;
t280 = pkin(5) * t444;
t286 = pkin(5) * t303 + pkin(4);
t327 = t286 * t442 - t305 * t448 + t280;
t180 = t327 - t248;
t434 = t168 + t180;
t241 = -t289 * t446 - t440;
t242 = t289 * t445 - t441;
t359 = -t242 * rSges(7,1) - t241 * rSges(7,2);
t449 = t287 * t288;
t167 = rSges(7,3) * t449 - t359;
t478 = pkin(4) - t286;
t315 = -t287 * t477 - t289 * t478;
t439 = t290 * t301;
t411 = pkin(5) * t439;
t179 = t288 * t315 - t411;
t435 = t167 + t179;
t493 = -t288 * t435 - t290 * t434;
t328 = rSges(4,1) * t299 - rSges(4,2) * sin(pkin(11)) + pkin(2);
t474 = rSges(4,3) + qJ(3);
t210 = t288 * t474 + t290 * t328 + t293;
t491 = 2 * m(6);
t490 = 2 * m(7);
t489 = t288 ^ 2;
t488 = t290 ^ 2;
t487 = t288 / 0.2e1;
t486 = -t289 / 0.2e1;
t485 = t290 / 0.2e1;
t484 = -rSges(6,3) - pkin(8);
t483 = m(5) * t264;
t482 = pkin(4) * t289;
t481 = t287 * pkin(4);
t479 = qJD(1) / 0.2e1;
t475 = rSges(7,3) * t287;
t473 = -rSges(7,3) + t305;
t379 = -t289 * t295 + qJD(1);
t314 = t287 * t418 + t292 * t379;
t424 = qJD(1) * t289;
t378 = -t295 + t424;
t155 = t288 * t314 - t378 * t441;
t417 = qJD(4) * t292;
t313 = -t287 * t417 + t291 * t379;
t156 = t288 * t313 + t378 * t440;
t360 = t156 * rSges(7,1) + t155 * rSges(7,2);
t106 = rSges(7,3) * t318 + t360;
t419 = qJD(4) * t290;
t387 = t289 * t419;
t472 = t106 * t448 + t167 * t387;
t469 = Icges(6,4) * t301;
t467 = Icges(7,4) * t291;
t438 = t290 * t303;
t249 = -t289 * t444 - t438;
t443 = t288 * t303;
t250 = t289 * t443 - t439;
t362 = -rSges(6,1) * t250 - rSges(6,2) * t249;
t194 = rSges(6,3) * t449 - t362;
t462 = t194 * t290;
t203 = (-Icges(6,2) * t303 - t469) * t414 + (Icges(6,6) * t287 + t289 * t343) * qJD(4);
t459 = t203 * t301;
t458 = t222 * t287;
t457 = t222 * t289;
t455 = t223 * t289;
t454 = t224 * t287;
t453 = t224 * t289;
t452 = t225 * t287;
t450 = t287 * t286;
t153 = t290 * t314 + t378 * t446;
t154 = t290 * t313 - t378 * t445;
t402 = t154 * rSges(7,1) + t153 * rSges(7,2) + rSges(7,3) * t387;
t105 = -rSges(7,3) * t393 + t402;
t270 = pkin(8) * t387;
t412 = qJD(5) * t303;
t408 = pkin(5) * t412;
t395 = qJD(1) * t411 + t288 * t408 + t305 * t393;
t413 = qJD(5) * t301;
t409 = pkin(5) * t413;
t437 = t105 - t270 + (pkin(8) * t425 + t419 * t478) * t287 + ((-qJD(4) * t305 - t409) * t290 + t478 * t425) * t289 + t395;
t358 = rSges(7,1) * t292 - rSges(7,2) * t291;
t232 = -rSges(7,3) * t289 + t287 * t358;
t422 = qJD(4) * t287;
t436 = t168 * t422 + t232 * t393;
t178 = (-rSges(7,1) * t291 - rSges(7,2) * t292) * t447 + (t289 * t358 + t475) * qJD(4);
t386 = t287 * t413;
t208 = -pkin(5) * t386 + qJD(4) * t315;
t433 = -t178 - t208;
t251 = -t289 * t439 + t443;
t252 = t289 * t438 + t444;
t195 = t252 * rSges(6,1) + t251 * rSges(6,2) + rSges(6,3) * t448;
t432 = -t195 - t248;
t361 = rSges(6,1) * t303 - rSges(6,2) * t301;
t205 = (-rSges(6,1) * t301 - rSges(6,2) * t303) * t414 + (rSges(6,3) * t287 + t289 * t361) * qJD(4);
t368 = pkin(8) * t287 + t482;
t258 = t368 * qJD(4);
t431 = -t205 - t258;
t137 = t289 * t167 + t232 * t449;
t219 = -t287 * t478 + t497;
t430 = t219 + t232;
t247 = t368 * t288;
t429 = t288 * t247 + t290 * t248;
t240 = -rSges(6,3) * t289 + t287 * t361;
t265 = -t289 * pkin(8) + t481;
t428 = -t240 - t265;
t283 = qJD(3) * t290;
t427 = t300 * t425 + t283;
t221 = Icges(5,3) * t288 + t290 * t341;
t426 = qJD(1) * t221;
t421 = qJD(4) * t288;
t415 = qJD(4) * t303;
t374 = -qJD(5) * t289 + qJD(1);
t312 = t287 * t416 + t303 * t374;
t373 = -qJD(5) + t424;
t174 = t290 * t312 + t373 * t444;
t311 = -t287 * t415 + t301 * t374;
t175 = t290 * t311 - t373 * t443;
t317 = t387 - t393;
t111 = Icges(6,5) * t175 + Icges(6,6) * t174 + Icges(6,3) * t317;
t113 = Icges(6,4) * t175 + Icges(6,2) * t174 + Icges(6,6) * t317;
t115 = Icges(6,1) * t175 + Icges(6,4) * t174 + Icges(6,5) * t317;
t189 = Icges(6,5) * t252 + Icges(6,6) * t251 + Icges(6,3) * t448;
t191 = Icges(6,4) * t252 + Icges(6,2) * t251 + Icges(6,6) * t448;
t193 = Icges(6,1) * t252 + Icges(6,4) * t251 + Icges(6,5) * t448;
t333 = -t191 * t301 + t193 * t303;
t32 = (qJD(4) * t333 - t111) * t289 + (qJD(4) * t189 - t113 * t301 + t115 * t303 + (-t191 * t303 - t193 * t301) * qJD(5)) * t287;
t347 = Icges(6,1) * t303 - t469;
t204 = (-Icges(6,1) * t301 - t468) * t414 + (Icges(6,5) * t287 + t289 * t347) * qJD(4);
t237 = -Icges(6,3) * t289 + t287 * t340;
t239 = -Icges(6,5) * t289 + t287 * t347;
t61 = t174 * t238 + t175 * t239 + t202 * t448 + t203 * t251 + t204 * t252 + t237 * t317;
t407 = t32 / 0.2e1 + t61 / 0.2e1;
t176 = t288 * t312 - t373 * t439;
t177 = t288 * t311 + t373 * t438;
t112 = Icges(6,5) * t177 + Icges(6,6) * t176 + Icges(6,3) * t318;
t114 = Icges(6,4) * t177 + Icges(6,2) * t176 + Icges(6,6) * t318;
t116 = Icges(6,1) * t177 + Icges(6,4) * t176 + Icges(6,5) * t318;
t188 = Icges(6,5) * t250 + Icges(6,6) * t249 + Icges(6,3) * t449;
t190 = Icges(6,4) * t250 + Icges(6,2) * t249 + Icges(6,6) * t449;
t192 = Icges(6,1) * t250 + Icges(6,4) * t249 + Icges(6,5) * t449;
t334 = -t190 * t301 + t192 * t303;
t31 = (qJD(4) * t334 - t112) * t289 + (qJD(4) * t188 - t114 * t301 + t116 * t303 + (-t190 * t303 - t192 * t301) * qJD(5)) * t287;
t62 = t176 * t238 + t177 * t239 + t202 * t449 + t203 * t249 + t204 * t250 + t237 * t318;
t406 = t62 / 0.2e1 + t31 / 0.2e1;
t405 = t230 * t292 * t295;
t125 = t237 * t449 + t238 * t249 + t239 * t250;
t94 = -t188 * t289 + t287 * t334;
t404 = t94 / 0.2e1 + t125 / 0.2e1;
t126 = t237 * t448 + t238 * t251 + t239 * t252;
t95 = -t189 * t289 + t287 * t333;
t403 = t95 / 0.2e1 + t126 / 0.2e1;
t346 = Icges(7,1) * t292 - t467;
t173 = (-Icges(7,1) * t291 - t466) * t447 + (Icges(7,5) * t287 + t289 * t346) * qJD(4);
t229 = -Icges(7,3) * t289 + t287 * t339;
t231 = -Icges(7,5) * t289 + t287 * t346;
t401 = t287 * t292 * t173 + t289 * t231 * t417 + t229 * t422;
t400 = t175 * rSges(6,1) + t174 * rSges(6,2) + rSges(6,3) * t387;
t399 = -t258 + t433;
t398 = t287 * t303 * t204 + t289 * t239 * t415 + t237 * t422;
t269 = t421 * t481;
t389 = t287 * t419;
t397 = t288 * (pkin(8) * t318 + qJD(1) * t276 - t269) + t290 * (-pkin(8) * t393 + t270 + (-t288 * t424 - t389) * pkin(4)) + t247 * t423;
t396 = -t265 - t430;
t394 = t240 * t425;
t385 = t449 / 0.2e1;
t384 = t448 / 0.2e1;
t381 = t430 * t290;
t207 = t428 * t290;
t380 = -t286 * t289 - t285;
t377 = t289 * t409;
t376 = t290 * t408;
t375 = t289 * t106 + t178 * t449 + t232 * t318;
t148 = t396 * t290;
t136 = -t229 * t289 + (-t230 * t291 + t231 * t292) * t287;
t161 = Icges(7,5) * t242 + Icges(7,6) * t241 + Icges(7,3) * t449;
t163 = Icges(7,4) * t242 + Icges(7,2) * t241 + Icges(7,6) * t449;
t165 = Icges(7,1) * t242 + Icges(7,4) * t241 + Icges(7,5) * t449;
t338 = -t163 * t291 + t165 * t292;
t89 = -t161 * t289 + t287 * t338;
t162 = Icges(7,5) * t244 + Icges(7,6) * t243 + Icges(7,3) * t448;
t164 = Icges(7,4) * t244 + Icges(7,2) * t243 + Icges(7,6) * t448;
t166 = Icges(7,1) * t244 + Icges(7,4) * t243 + Icges(7,5) * t448;
t337 = -t164 * t291 + t166 * t292;
t90 = -t162 * t289 + t287 * t337;
t351 = t288 * t89 + t290 * t90;
t120 = t229 * t449 + t230 * t241 + t231 * t242;
t73 = t161 * t449 + t163 * t241 + t165 * t242;
t74 = t162 * t449 + t164 * t241 + t166 * t242;
t356 = t288 * t73 + t290 * t74;
t40 = -t120 * t289 + t287 * t356;
t121 = t229 * t448 + t230 * t243 + t231 * t244;
t75 = t161 * t448 + t163 * t243 + t165 * t244;
t76 = t162 * t448 + t164 * t243 + t166 * t244;
t355 = t288 * t75 + t290 * t76;
t41 = -t121 * t289 + t287 * t355;
t100 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t318;
t102 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t318;
t104 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t318;
t19 = t100 * t448 + t102 * t243 + t104 * t244 + t153 * t163 + t154 * t165 + t161 * t317;
t101 = Icges(7,4) * t154 + Icges(7,2) * t153 + Icges(7,6) * t317;
t103 = Icges(7,1) * t154 + Icges(7,4) * t153 + Icges(7,5) * t317;
t99 = Icges(7,5) * t154 + Icges(7,6) * t153 + Icges(7,3) * t317;
t20 = t101 * t243 + t103 * t244 + t153 * t164 + t154 * t166 + t162 * t317 + t448 * t99;
t172 = (-Icges(7,2) * t292 - t467) * t447 + (Icges(7,6) * t287 + t289 * t342) * qJD(4);
t54 = t153 * t230 + t154 * t231 + t171 * t448 + t172 * t243 + t173 * t244 + t229 * t317;
t57 = t288 * t76 - t290 * t75;
t5 = (qJD(4) * t355 - t54) * t289 + (-qJD(1) * t57 + qJD(4) * t121 + t19 * t288 + t20 * t290) * t287;
t21 = t100 * t449 + t102 * t241 + t104 * t242 + t155 * t163 + t156 * t165 + t161 * t318;
t22 = t101 * t241 + t103 * t242 + t155 * t164 + t156 * t166 + t162 * t318 + t449 * t99;
t55 = t155 * t230 + t156 * t231 + t171 * t449 + t172 * t241 + t173 * t242 + t229 * t318;
t56 = t288 * t74 - t290 * t73;
t6 = (qJD(4) * t356 - t55) * t289 + (-qJD(1) * t56 + qJD(4) * t120 + t21 * t288 + t22 * t290) * t287;
t366 = t41 * t387 + t5 * t448 + t6 * t449 + (-t136 * t289 + t287 * t351) * t422 + t318 * t40;
t363 = rSges(6,1) * t177 + rSges(6,2) * t176;
t357 = -t290 * t300 - t480;
t85 = t188 * t449 + t190 * t249 + t192 * t250;
t86 = t189 * t449 + t191 * t249 + t193 * t250;
t58 = t288 * t86 - t290 * t85;
t354 = t288 * t85 + t290 * t86;
t87 = t188 * t448 + t190 * t251 + t192 * t252;
t88 = t189 * t448 + t191 * t251 + t193 * t252;
t59 = t288 * t88 - t290 * t87;
t353 = t288 * t87 + t290 * t88;
t352 = t288 * t90 - t290 * t89;
t350 = t288 * t94 + t290 * t95;
t344 = Icges(5,2) * t289 + t471;
t332 = -t195 * t288 + t462;
t331 = -t194 * t288 - t195 * t290;
t330 = -t453 + t458;
t324 = t179 * t290 - t288 * t434;
t322 = qJD(4) * t344;
t321 = qJD(4) * (-Icges(5,5) * t287 - Icges(5,6) * t289);
t320 = t287 * t484 - t285 - t482;
t316 = t287 * t473 + t380;
t12 = qJD(1) * t355 - t19 * t290 + t20 * t288;
t13 = qJD(1) * t356 - t21 * t290 + t22 * t288;
t25 = (qJD(4) * t338 - t100) * t289 + (qJD(4) * t161 + (-t163 * t295 + t104) * t292 + (-t165 * t295 - t102) * t291) * t287;
t26 = (qJD(4) * t337 - t99) * t289 + (qJD(4) * t162 + (-t164 * t295 + t103) * t292 + (-t166 * t295 - t101) * t291) * t287;
t310 = t5 * t487 + t13 * t385 + t12 * t384 + (qJD(1) * t351 - t25 * t290 + t26 * t288) * t486 + t40 * t425 / 0.2e1 + t41 * t383 - t290 * t6 / 0.2e1 + t352 * t422 / 0.2e1 + t500 * t57 + t499 * t56;
t131 = t136 * t422;
t63 = t502 * t289 + (-t405 + (-t231 * t295 - t172) * t291) * t287 + t401;
t7 = t131 + (qJD(4) * t351 - t63) * t289 + (-qJD(1) * t352 + t25 * t288 + t26 * t290) * t287;
t308 = -t289 * t7 - t393 * t41 + t366;
t307 = t131 + (t25 + t55) * t385 + (t26 + t54) * t384 + (t121 + t90) * t500 + (t120 + t89) * t499;
t306 = t288 * t320 + t357;
t209 = -t288 * t328 + t290 * t474 - t480;
t282 = qJD(3) * t288;
t257 = t364 * qJD(4);
t253 = t265 * t425;
t206 = t428 * t288;
t197 = -qJD(1) * t210 + t283;
t196 = qJD(1) * t209 + t282;
t182 = t288 * t321 + t426;
t181 = -qJD(1) * t220 + t290 * t321;
t152 = t167 * t448;
t147 = t396 * t288;
t146 = t264 * t421 + (t290 * t326 - t284 - t293) * qJD(1) + t427;
t145 = t282 + ((-t285 - t476) * t288 + t357) * qJD(1) + t309;
t144 = -t237 * t289 + (-t238 * t301 + t239 * t303) * t287;
t143 = -t195 * t289 - t240 * t448;
t142 = t194 * t289 + t240 * t449;
t141 = t367 - t432;
t140 = t306 + t362;
t139 = t144 * t422;
t138 = -t168 * t289 - t232 * t448;
t135 = t221 * t288 - t496;
t134 = t220 * t288 - t290 * t330;
t133 = -t221 * t290 - t329 * t288;
t132 = -t220 * t290 - t288 * t330;
t130 = t327 + t367 + t168;
t129 = -t480 + (pkin(5) * t301 - t300) * t290 + t316 * t288 + t359;
t128 = qJD(1) * t207 + t288 * t431;
t127 = t290 * t431 + t253 + t394;
t124 = t332 * t287;
t123 = -t376 + t269 + (-t377 + (-t450 - t497) * qJD(4)) * t288 + (t290 * t315 + t280) * qJD(1);
t119 = -t168 * t449 + t152;
t118 = rSges(6,3) * t318 + t363;
t117 = -rSges(6,3) * t393 + t400;
t98 = -t331 + t429;
t92 = -t287 * t381 - t289 * t434;
t91 = t179 * t289 + t219 * t449 + t137;
t84 = t269 + t484 * t388 + (t290 * t320 - t293) * qJD(1) - t363 + t427;
t83 = -pkin(4) * t389 + qJD(1) * t306 + t270 + t282 + t400;
t82 = qJD(1) * t148 + t288 * t399;
t81 = t290 * t399 + t425 * t430 + t253;
t72 = t287 * t324 + t152;
t71 = t376 + (t377 + (t289 * t473 + t450) * qJD(4)) * t288 + (t290 * t316 - t280 - t293) * qJD(1) - t360 + t427;
t70 = t282 + (-t377 + (-t289 * t305 - t450) * qJD(4)) * t290 + ((t380 - t475) * t288 + t357) * qJD(1) + t395 + t402;
t69 = t429 - t493;
t68 = t501 * t289 + (-t459 + (-t238 * t303 - t239 * t301) * qJD(5)) * t287 + t398;
t67 = (t240 * t421 + t118) * t289 + (-qJD(4) * t194 + t205 * t288 + t240 * t423) * t287;
t66 = (-t240 * t419 - t117) * t289 + (qJD(4) * t195 - t205 * t290 + t394) * t287;
t65 = -t167 * t422 + t375;
t64 = -t178 * t448 + (-t232 * t419 - t105) * t289 + t436;
t46 = t332 * t420 + (qJD(1) * t331 - t117 * t288 + t118 * t290) * t287;
t45 = t117 * t290 + t118 * t288 + (t288 * t432 + t462) * qJD(1) + t397;
t44 = -t126 * t289 + t287 * t353;
t43 = -t125 * t289 + t287 * t354;
t42 = -t168 * t388 + (-t105 * t288 + (-t167 * t288 - t168 * t290) * qJD(1)) * t287 + t472;
t34 = (t219 * t421 + t123) * t289 + (-qJD(4) * t435 + t208 * t288 + t219 * t423) * t287 + t375;
t33 = (-qJD(4) * t381 - t437) * t289 + (qJD(4) * t180 + t219 * t425 + t290 * t433) * t287 + t436;
t30 = t111 * t449 + t113 * t249 + t115 * t250 + t176 * t191 + t177 * t193 + t189 * t318;
t29 = t112 * t449 + t114 * t249 + t116 * t250 + t176 * t190 + t177 * t192 + t188 * t318;
t28 = t111 * t448 + t113 * t251 + t115 * t252 + t174 * t191 + t175 * t193 + t189 * t317;
t27 = t112 * t448 + t114 * t251 + t116 * t252 + t174 * t190 + t175 * t192 + t188 * t317;
t18 = t437 * t290 + (t106 + t123) * t288 + (t435 * t290 + (-t248 - t434) * t288) * qJD(1) + t397;
t17 = t324 * t420 + (qJD(1) * t493 + t123 * t290 - t437 * t288) * t287 + t472;
t16 = qJD(1) * t354 + t288 * t30 - t29 * t290;
t15 = qJD(1) * t353 - t27 * t290 + t28 * t288;
t9 = (qJD(4) * t354 - t62) * t289 + (-qJD(1) * t58 + qJD(4) * t125 + t288 * t29 + t290 * t30) * t287;
t8 = (qJD(4) * t353 - t61) * t289 + (-t59 * qJD(1) + qJD(4) * t126 + t27 * t288 + t28 * t290) * t287;
t1 = [t401 + (t145 * t201 + t146 * t200) * t492 + 0.2e1 * m(4) * (t196 * t210 + t197 * t209) + (t140 * t84 + t141 * t83) * t491 + (t129 * t71 + t130 * t70) * t490 - t291 * t231 * t447 + t398 - t239 * t386 + (-t344 + t349) * t422 + (Icges(5,1) * t287 + t345 + t470) * t420 + (t501 + t502) * t289 + (-t172 * t291 - t238 * t412 - t405 - t459) * t287; 0; 0; m(7) * (t288 * t71 - t290 * t70 + (t129 * t290 + t130 * t288) * qJD(1)) + m(6) * (t288 * t84 - t290 * t83 + (t140 * t290 + t141 * t288) * qJD(1)) + m(5) * (qJD(1) * t494 - t145 * t290 + t146 * t288) + m(4) * (-t196 * t290 + t197 * t288 + (t209 * t290 + t210 * t288) * qJD(1)); 0; 0; ((qJD(1) * t223 - t288 * t322) * t486 + t225 * t503 - t55 / 0.2e1 - t25 / 0.2e1 + (t458 / 0.2e1 - t453 / 0.2e1) * qJD(4) - t406) * t290 + ((-qJD(1) * t222 - t290 * t322) * t289 / 0.2e1 + t224 * t503 + t54 / 0.2e1 + t26 / 0.2e1 + (-t456 / 0.2e1 + t451 / 0.2e1) * qJD(4) + t407) * t288 + m(5) * ((-t145 * t288 - t146 * t290) * t264 - t494 * t257) + m(6) * (t127 * t140 + t128 * t141 + t206 * t83 + t207 * t84) + m(7) * (t129 * t81 + t130 * t82 + t147 * t70 + t148 * t71) + (t489 / 0.2e1 + t488 / 0.2e1) * t341 * qJD(4) + ((t455 / 0.2e1 + t452 / 0.2e1 - t201 * t483 + t121 / 0.2e1 + t90 / 0.2e1 + t403) * t290 + (t200 * t483 + t457 / 0.2e1 + t454 / 0.2e1 + t120 / 0.2e1 + t89 / 0.2e1 + t404) * t288) * qJD(1); m(5) * t93 + m(6) * t45 + m(7) * t18; m(6) * (t127 * t288 - t128 * t290 + (t206 * t288 + t207 * t290) * qJD(1)) + m(7) * (t288 * t81 - t290 * t82 + (t147 * t288 + t148 * t290) * qJD(1)); (t147 * t82 + t148 * t81 + t69 * t18) * t490 + (t127 * t207 + t128 * t206 + t45 * t98) * t491 + (t488 + t489) * t264 * t257 * t492 + (t227 * t498 - t13 - t16 + (-t133 * qJD(1) + (-qJD(1) * t330 - t182) * t290) * t290) * t290 + (t12 + t15 + t226 * t498 + (t134 * qJD(1) + (t329 * qJD(1) + t181) * t288) * t288 + ((-t182 + (-t452 - t455) * qJD(4) + t223 * t420 + t225 * t422 - t426) * t288 + (t222 * t420 + t224 * t422 + t181 - (t454 + t457) * qJD(4)) * t290 + (t135 - t132 + (t221 - t330) * t288 + t496) * qJD(1)) * t290) * t288 + (-t132 * t290 + t133 * t288 + t56 + t58) * t425 + (-t134 * t290 + t135 * t288 + t57 + t59) * t423; m(6) * (t140 * t67 + t141 * t66 + t142 * t84 + t143 * t83) + m(7) * (t129 * t34 + t130 * t33 + t70 * t92 + t71 * t91) + t139 + (t407 * t290 + t406 * t288 + (-t288 * t403 + t290 * t404) * qJD(1)) * t287 + t307 + (-t63 - t68 + (t288 * t404 + t290 * t403) * qJD(4)) * t289; m(6) * t46 + m(7) * t17; m(6) * (t288 * t67 - t290 * t66 + (t142 * t290 + t143 * t288) * qJD(1)) + m(7) * (t288 * t34 - t290 * t33 + (t288 * t92 + t290 * t91) * qJD(1)); (qJD(4) * (t288 * t95 - t290 * t94) / 0.2e1 + t15 * t485 + t16 * t487 + (t58 * t485 - t288 * t59 / 0.2e1) * qJD(1)) * t287 + (t43 * t479 + t58 * t382 + (qJD(1) * t94 + t32) * t486 + t8 / 0.2e1) * t288 + m(6) * (t124 * t45 + t127 * t142 + t128 * t143 + t206 * t66 + t207 * t67 + t46 * t98) + m(7) * (t147 * t33 + t148 * t34 + t17 * t69 + t18 * t72 + t81 * t91 + t82 * t92) + (t44 * t479 + t59 * t382 + (qJD(1) * t95 - t31) * t486 - t9 / 0.2e1) * t290 + t310; (t17 * t72 + t33 * t92 + t34 * t91) * t490 + (t124 * t46 + t142 * t67 + t143 * t66) * t491 + (t68 * t289 - t139 - t7 + (t288 * t43 - t289 * t350 + t290 * t44) * qJD(4)) * t289 + (t288 * t9 + t290 * t8 - t289 * (t31 * t288 + t32 * t290) + (-t144 * t289 + t287 * t350) * qJD(4) + ((-t289 * t94 + t43) * t290 + (t289 * t95 - t41 - t44) * t288) * qJD(1)) * t287 + t366; m(7) * (t129 * t65 + t130 * t64 + t137 * t71 + t138 * t70) + t307 - t63 * t289; m(7) * t42; m(7) * (t288 * t65 - t290 * t64 + (t137 * t290 + t138 * t288) * qJD(1)); m(7) * (t119 * t18 + t137 * t81 + t138 * t82 + t147 * t64 + t148 * t65 + t42 * t69) + t310; m(7) * (t119 * t17 + t137 * t34 + t138 * t33 + t42 * t72 + t64 * t92 + t65 * t91) + t308; (t119 * t42 + t137 * t65 + t138 * t64) * t490 + t308;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
