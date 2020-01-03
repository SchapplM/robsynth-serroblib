% Calculate vector of inverse dynamics joint torques for
% S5RRPRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:48
% EndTime: 2020-01-03 12:02:02
% DurationCPUTime: 10.09s
% Computational Cost: add. (19139->642), mult. (11926->816), div. (0->0), fcn. (9239->10), ass. (0->369)
t307 = qJ(1) + qJ(2);
t294 = sin(t307);
t283 = pkin(2) * t294;
t292 = pkin(9) + t307;
t281 = sin(t292);
t282 = cos(t292);
t367 = -rSges(4,1) * t281 - rSges(4,2) * t282;
t546 = t367 - t283;
t305 = qJD(1) + qJD(2);
t296 = cos(t307);
t525 = rSges(3,1) * t294 + rSges(3,2) * t296;
t198 = t525 * t305;
t309 = sin(qJ(1));
t486 = pkin(1) * qJD(1);
t404 = t309 * t486;
t177 = t404 + t198;
t306 = qJ(4) + qJ(5);
t293 = sin(t306);
t295 = cos(t306);
t477 = Icges(6,4) * t293;
t220 = Icges(6,1) * t295 - t477;
t148 = -Icges(6,5) * t282 + t220 * t281;
t449 = t282 * t293;
t230 = Icges(6,4) * t449;
t448 = t282 * t295;
t149 = Icges(6,1) * t448 + Icges(6,5) * t281 - t230;
t304 = qJD(4) + qJD(5);
t204 = t281 * t304;
t205 = t282 * t304;
t217 = Icges(6,2) * t295 + t477;
t276 = Icges(6,4) * t295;
t359 = -Icges(6,2) * t293 + t276;
t526 = Icges(6,1) * t293 + t276;
t541 = t526 + t359;
t319 = t204 * (-Icges(6,2) * t448 + t149 - t230) - t205 * (-t217 * t281 + t148) + t305 * t541;
t146 = -Icges(6,6) * t282 + t281 * t359;
t147 = Icges(6,4) * t448 - Icges(6,2) * t449 + Icges(6,6) * t281;
t506 = t204 * (t282 * t526 + t147) - t205 * (t281 * t526 + t146) + t305 * (t217 - t220);
t545 = t293 * t319 + t295 * t506;
t308 = sin(qJ(4));
t310 = cos(qJ(4));
t297 = Icges(5,4) * t310;
t360 = -Icges(5,2) * t308 + t297;
t523 = Icges(5,1) * t308 + t297;
t418 = t523 + t360;
t478 = Icges(5,4) * t308;
t260 = Icges(5,2) * t310 + t478;
t263 = Icges(5,1) * t310 - t478;
t419 = t260 - t263;
t544 = (t308 * t418 + t310 * t419) * t305;
t489 = pkin(3) * t282;
t209 = pkin(7) * t281 + t489;
t440 = t296 * t305;
t266 = pkin(2) * t440;
t543 = -t209 * t305 + t266;
t267 = rSges(5,1) * t308 + rSges(5,2) * t310;
t188 = t267 * t281;
t189 = t267 * t282;
t450 = t281 * t310;
t253 = rSges(5,1) * t450;
t451 = t281 * t308;
t406 = rSges(5,2) * t451;
t370 = t253 - t406;
t163 = -rSges(5,3) * t282 + t370;
t275 = t281 * pkin(3);
t208 = -pkin(7) * t282 + t275;
t387 = t208 + t283;
t372 = t163 + t387;
t414 = qJD(4) * t282;
t399 = t267 * t414;
t78 = t305 * t372 + t399 + t404;
t311 = cos(qJ(1));
t291 = t311 * t486;
t415 = qJD(4) * t281;
t369 = -t267 * t415 + t266;
t446 = t282 * t308;
t254 = rSges(5,2) * t446;
t445 = t282 * t310;
t409 = rSges(5,1) * t445;
t164 = rSges(5,3) * t281 - t254 + t409;
t428 = t164 + t209;
t79 = t305 * t428 + t291 + t369;
t542 = -t188 * t78 - t189 * t79;
t410 = -qJDD(4) - qJDD(5);
t133 = -t205 * t305 + t281 * t410;
t447 = t282 * t305;
t248 = pkin(7) * t447;
t452 = t281 * t305;
t179 = pkin(3) * t452 - t248;
t413 = qJD(4) * t305;
t190 = -qJDD(4) * t281 - t282 * t413;
t279 = t295 * rSges(6,1);
t223 = -rSges(6,2) * t293 + t279;
t197 = t223 * t304;
t221 = rSges(6,1) * t293 + rSges(6,2) * t295;
t303 = qJDD(1) + qJDD(2);
t299 = t309 * pkin(1);
t314 = qJD(1) ^ 2;
t468 = pkin(1) * qJDD(1);
t373 = -t299 * t314 + t311 * t468;
t490 = pkin(2) * t303;
t302 = t305 ^ 2;
t491 = pkin(2) * t302;
t324 = -t294 * t491 + t296 * t490 + t373;
t300 = t310 * pkin(4);
t285 = t300 + pkin(3);
t232 = t282 * t285;
t312 = -pkin(8) - pkin(7);
t140 = t489 - t232 + (pkin(7) + t312) * t281;
t408 = rSges(6,1) * t448;
t371 = -rSges(6,2) * t449 + t408;
t152 = rSges(6,3) * t281 + t371;
t437 = -t140 + t152;
t401 = t209 + t437;
t439 = t310 * qJD(4) ^ 2;
t271 = t282 * t312;
t412 = qJD(4) * t308;
t396 = t282 * t412;
t376 = pkin(4) * t396;
t112 = t376 + t248 + (t271 + (-pkin(3) + t285) * t281) * t305;
t236 = rSges(6,2) * t448;
t443 = t293 * t305;
t405 = rSges(6,2) * t443;
t427 = -rSges(6,3) * t447 - t281 * t405;
t444 = t293 * t304;
t94 = t304 * t236 + (t282 * t444 + t295 * t452) * rSges(6,1) + t427;
t481 = -t112 - t94;
t18 = t133 * t221 - t197 * t204 + (t190 * t308 - t281 * t439) * pkin(4) + (-t179 + t481) * t305 + t401 * t303 + t324;
t540 = -g(2) + t18;
t411 = qJD(4) * t310;
t425 = rSges(5,3) * t447 + t305 * t406;
t110 = rSges(5,2) * t282 * t411 + (t305 * t450 + t396) * rSges(5,1) - t425;
t488 = rSges(5,1) * t310;
t269 = -rSges(5,2) * t308 + t488;
t237 = t269 * qJD(4);
t46 = -t237 * t415 + t190 * t267 + (-t110 - t179) * t305 + t428 * t303 + t324;
t539 = -g(2) + t46;
t272 = t281 * rSges(4,2);
t207 = rSges(4,1) * t282 - t272;
t538 = t207 * t303 + t302 * t367 - g(2) + t324;
t210 = t285 * t447;
t403 = pkin(4) * t412;
t346 = -t305 * t312 - t403;
t420 = pkin(3) * t447 + pkin(7) * t452;
t113 = t281 * t346 + t210 - t420;
t241 = t281 * t413;
t134 = qJD(5) * t452 + t282 * t410 + t241;
t191 = -qJDD(4) * t282 + t241;
t301 = t311 * pkin(1);
t416 = t301 * t314 + t309 * t468;
t374 = t294 * t490 + t296 * t491 + t416;
t347 = t208 * t303 + t305 * t420 + t374;
t421 = t281 * t285 + t271;
t139 = -t208 + t421;
t453 = t281 * t295;
t454 = t281 * t293;
t151 = rSges(6,1) * t453 - rSges(6,2) * t454 - rSges(6,3) * t282;
t438 = t139 + t151;
t407 = rSges(6,1) * t444;
t426 = rSges(6,3) * t452 + t305 * t408;
t441 = t295 * t304;
t95 = -t281 * t407 + (-t281 * t441 - t282 * t443) * rSges(6,2) + t426;
t17 = -t134 * t221 + t197 * t205 + (t113 + t95) * t305 + t438 * t303 + (-t191 * t308 + t282 * t439) * pkin(4) + t347;
t537 = -g(3) + t17;
t397 = t281 * t411;
t398 = t281 * t412;
t424 = rSges(5,3) * t452 + t305 * t409;
t111 = -rSges(5,1) * t398 + (-t305 * t446 - t397) * rSges(5,2) + t424;
t45 = t111 * t305 + t163 * t303 - t191 * t267 + t237 * t414 + t347;
t536 = -g(3) + t45;
t246 = rSges(4,1) * t447;
t535 = -g(3) - t303 * t367 + t305 * (-rSges(4,2) * t452 + t246) + t374;
t224 = rSges(3,1) * t296 - rSges(3,2) * t294;
t534 = -t198 * t305 + t224 * t303 - g(2) + t373;
t442 = t294 * t305;
t199 = rSges(3,1) * t440 - rSges(3,2) * t442;
t533 = t199 * t305 + t303 * t525 - g(3) + t416;
t216 = Icges(6,5) * t295 - Icges(6,6) * t293;
t144 = -Icges(6,3) * t282 + t216 * t281;
t463 = t147 * t293;
t358 = -t149 * t295 + t463;
t345 = -t144 + t358;
t531 = t205 * t345;
t530 = t305 * t546;
t196 = t305 * t207;
t528 = t246 - t196;
t527 = -t275 - t283;
t141 = t305 * t152;
t333 = -pkin(4) * t398 - t204 * t221 + t266;
t522 = t140 * t305 - t141 + t210 - t333 + t426 + t543;
t150 = t305 * t164;
t521 = -t150 - t369 + t420 + t424 + t543;
t456 = t217 * t304;
t520 = -Icges(6,6) * t305 + t456;
t174 = t221 * t281;
t175 = rSges(6,1) * t449 + t236;
t53 = t151 * t204 + t152 * t205 + qJD(3) + (t139 * t281 - t140 * t282) * qJD(4);
t348 = -t205 * t221 - t376;
t349 = t387 + t438;
t59 = t305 * t349 - t348 + t404;
t60 = t305 * t401 + t291 + t333;
t519 = -t59 * (-t174 * t305 + t205 * t223) - t53 * (-t174 * t204 - t175 * t205) - t60 * (-t175 * t305 - t204 * t223);
t514 = -Icges(5,6) * t305 + qJD(4) * t260;
t107 = -t281 * t514 + t360 * t447;
t344 = t263 * t305;
t511 = -Icges(5,5) * t305 + qJD(4) * t523;
t109 = -t281 * t511 + t282 * t344;
t259 = Icges(5,5) * t310 - Icges(5,6) * t308;
t156 = -Icges(5,3) * t282 + t259 * t281;
t342 = t360 * t281;
t158 = -Icges(5,6) * t282 + t342;
t160 = -Icges(5,5) * t282 + t263 * t281;
t99 = t158 * t310 + t160 * t308;
t518 = qJD(4) * t99 + t107 * t308 - t109 * t310 - t156 * t305;
t215 = Icges(6,5) * t293 + Icges(6,6) * t295;
t517 = -Icges(6,3) * t305 + t215 * t304;
t516 = -Icges(6,5) * t305 + t304 * t526;
t258 = Icges(5,5) * t308 + Icges(5,6) * t310;
t515 = -Icges(5,3) * t305 + qJD(4) * t258;
t226 = t360 * qJD(4);
t227 = t263 * qJD(4);
t352 = t260 * t310 + t308 * t523;
t513 = qJD(4) * t352 + t226 * t308 - t227 * t310 - t258 * t305;
t106 = t282 * t514 + t305 * t342;
t108 = t281 * t344 + t282 * t511;
t157 = Icges(5,5) * t445 - Icges(5,6) * t446 + Icges(5,3) * t281;
t159 = Icges(5,4) * t445 - Icges(5,2) * t446 + Icges(5,6) * t281;
t252 = Icges(5,4) * t446;
t161 = Icges(5,1) * t445 + Icges(5,5) * t281 - t252;
t356 = t159 * t310 + t161 * t308;
t512 = qJD(4) * t356 - t106 * t308 + t108 * t310 - t157 * t305;
t377 = -t220 * t304 + t456;
t378 = t541 * t304;
t509 = -t215 * t305 + t293 * t378 + t295 * t377;
t145 = Icges(6,5) * t448 - Icges(6,6) * t449 + Icges(6,3) * t281;
t341 = t359 * t305;
t382 = t149 * t304 - t281 * t341 - t282 * t520;
t343 = t220 * t305;
t384 = t147 * t304 + t281 * t343 + t282 * t516;
t508 = -t145 * t305 + t293 * t382 + t295 * t384;
t383 = t148 * t304 - t281 * t520 + t282 * t341;
t385 = t146 * t304 + t281 * t516 - t282 * t343;
t507 = -t144 * t305 + t293 * t383 + t295 * t385;
t505 = t133 / 0.2e1;
t504 = t134 / 0.2e1;
t503 = t190 / 0.2e1;
t502 = t191 / 0.2e1;
t501 = -t204 / 0.2e1;
t500 = t204 / 0.2e1;
t499 = t205 / 0.2e1;
t498 = -t205 / 0.2e1;
t497 = -t281 / 0.2e1;
t496 = -t282 / 0.2e1;
t495 = t303 / 0.2e1;
t494 = -t305 / 0.2e1;
t493 = t305 / 0.2e1;
t492 = rSges(5,3) + pkin(7);
t484 = t305 * t60;
t483 = t305 * t79;
t458 = t215 * t281;
t98 = -t217 * t449 + t448 * t526 + t458;
t482 = t98 * t305;
t455 = t258 * t281;
t118 = -t260 * t446 + t445 * t523 + t455;
t467 = t118 * t305;
t464 = t146 * t293;
t462 = t148 * t295;
t461 = t158 * t308;
t460 = t159 * t308;
t459 = t160 * t310;
t169 = t215 * t282;
t339 = t216 * t305;
t183 = t258 * t282;
t340 = t259 * t305;
t432 = t281 * t523 + t158;
t431 = t282 * t523 + t159;
t430 = -t260 * t281 + t160;
t429 = -Icges(5,2) * t445 + t161 - t252;
t402 = t151 * t447 + (-t141 + t95) * t281;
t400 = t248 + t425;
t394 = t452 / 0.2e1;
t393 = -t447 / 0.2e1;
t392 = -t415 / 0.2e1;
t391 = t415 / 0.2e1;
t390 = -t414 / 0.2e1;
t389 = t414 / 0.2e1;
t386 = -pkin(4) * t308 - t221;
t381 = -t145 - t464;
t380 = -t145 + t462;
t178 = t224 * t305 + t291;
t284 = pkin(2) * t296;
t181 = t207 + t284;
t270 = rSges(2,1) * t311 - rSges(2,2) * t309;
t268 = rSges(2,1) * t309 + rSges(2,2) * t311;
t136 = t160 * t450;
t69 = -t156 * t282 - t158 * t451 + t136;
t137 = t161 * t450;
t70 = t157 * t282 + t159 * t451 - t137;
t366 = -t281 * t70 - t282 * t69;
t138 = t158 * t446;
t71 = -t156 * t281 - t160 * t445 + t138;
t355 = -t161 * t310 + t460;
t72 = t157 * t281 - t355 * t282;
t365 = -t281 * t72 - t282 * t71;
t364 = -t281 * t79 + t282 * t78;
t84 = -t147 * t295 - t149 * t293;
t357 = t459 - t461;
t354 = t163 * t281 + t164 * t282;
t353 = -t217 * t293 + t295 * t526;
t351 = -t260 * t308 + t310 * t523;
t350 = -pkin(2) * t442 - t427;
t335 = t169 * t204 - t205 * t458 - t339;
t332 = -t281 * t339 - t282 * t517 + t305 * t358;
t331 = -t282 * t339 + t517 * t281 + (t462 - t464) * t305;
t330 = -t281 * t340 - t282 * t515 + t305 * t355;
t329 = t281 * t515 - t282 * t340 + t305 * t357;
t328 = -t216 * t304 + t305 * t353;
t327 = -qJD(4) * t259 + t305 * t351;
t326 = t308 * t430 + t310 * t432;
t325 = t308 * t429 + t310 * t431;
t114 = t151 + t283 + t421;
t123 = -t282 * t492 + t370 - t527;
t115 = t232 + t284 + (rSges(6,3) - t312) * t281 + t371;
t13 = t281 * t331 + t282 * t507;
t14 = t281 * t332 - t282 * t508;
t15 = -t281 * t507 + t282 * t331;
t16 = t281 * t508 + t282 * t332;
t125 = t148 * t453;
t63 = -t144 * t282 - t146 * t454 + t125;
t126 = t149 * t453;
t64 = t145 * t282 + t147 * t454 - t126;
t97 = t281 * t353 - t169;
t96 = t97 * t305;
t30 = -t204 * t64 - t205 * t63 + t96;
t127 = t146 * t449;
t65 = -t144 * t281 - t148 * t448 + t127;
t66 = t145 * t281 - t282 * t358;
t31 = -t204 * t66 - t205 * t65 - t482;
t40 = -t293 * t385 + t295 * t383;
t41 = t293 * t384 - t295 * t382;
t47 = t281 * t328 + t282 * t509;
t48 = -t281 * t509 + t282 * t328;
t83 = t146 * t295 + t148 * t293;
t323 = (-t13 * t205 + t133 * t66 + t134 * t65 - t14 * t204 - t303 * t98 + t305 * t47) * t497 + (t335 * t281 + t282 * t545) * t500 + (-t281 * t545 + t335 * t282) * t499 + (t133 * t64 + t134 * t63 - t15 * t205 - t16 * t204 + t303 * t97 + t305 * t48) * t496 + (-t293 * t506 + t295 * t319) * t494 + t30 * t394 + t31 * t393 + ((-t305 * t66 - t13) * t282 + (t305 * t65 - t14) * t281) * t501 + (-t281 * t66 - t282 * t65) * t505 + (-t281 * t64 - t282 * t63) * t504 + ((-t305 * t64 - t15) * t282 + (t305 * t63 - t16) * t281) * t498 + (-t281 * t84 - t282 * t83) * t495 + ((-t305 * t84 - t40) * t282 + (t305 * t83 - t41) * t281) * t493;
t124 = -t254 + t284 + (pkin(3) + t488) * t282 + t492 * t281;
t142 = t404 - t530;
t143 = t291 + t266 + t196;
t318 = (-t142 * t272 + t143 * t546) * t305;
t117 = t281 * t351 - t183;
t116 = t117 * t305;
t36 = qJD(4) * t366 + t116;
t37 = qJD(4) * t365 - t467;
t51 = qJD(4) * t357 + t107 * t310 + t109 * t308;
t52 = qJD(4) * t355 + t106 * t310 + t108 * t308;
t57 = t281 * t327 + t282 * t513;
t58 = -t281 * t513 + t282 * t327;
t317 = -t356 * t503 + t84 * t505 + (t116 + ((t71 + t137 - t138 + (t156 - t460) * t281) * t281 + (-t136 - t72 + (t156 - t355) * t282 + (t459 + t461) * t281) * t282) * qJD(4)) * t391 - t133 * t98 / 0.2e1 - t190 * t118 / 0.2e1 + (t96 - (t281 * t381 + t125 + t66) * t205 + (t126 - t127 + t65 + (t144 - t463) * t281) * t204 + (t204 * t380 - t531) * t282) * t500 + (t83 + t97) * t504 + (t117 + t99) * t502 + (t31 + t482 - (-t127 + t64) * t205 + (-t125 + t63) * t204 + (-t204 * t345 - t205 * t380) * t282 + (-t204 * t381 + t531) * t281) * t499 + (t40 + t48) * t498 + (t51 + t58) * t390 + (t467 + ((t138 - t70 + (t157 - t459) * t282) * t282 + (-t136 + t69 + (t157 + t461) * t281) * t281) * qJD(4) + t37) * t389 + (t41 + t47 + t30) * t501 + (t52 + t57 + t36) * t392 + (qJD(4) * t351 + t226 * t310 + t227 * t308 - t293 * t377 + t295 * t378) * t305 + (t217 * t295 + t293 * t526 + Icges(3,3) + Icges(4,3) + t352) * t303;
t316 = (t79 * (-t253 + t527) - t78 * t254) * t305 + t542 * qJD(4);
t315 = (t59 * (-t221 * t304 - t403) + (t60 * (-t285 - t279) - t59 * t312) * t305) * t281 + (t60 * (-rSges(6,2) * t441 + t346 - t407) - t59 * t405) * t282;
t255 = pkin(4) * t446;
t135 = t281 * t151;
t87 = qJD(4) * t354 + qJD(3);
t44 = -t163 * t190 - t164 * t191 + qJDD(3) + (-t110 * t282 + t111 * t281) * qJD(4);
t22 = t281 * t512 + t282 * t330;
t21 = -t281 * t518 + t282 * t329;
t20 = t281 * t330 - t282 * t512;
t19 = t281 * t329 + t282 * t518;
t10 = -t133 * t151 - t134 * t152 - t139 * t190 + t140 * t191 + t204 * t95 - t205 * t94 + qJDD(3) + (-t112 * t282 + t113 * t281) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t317 + (t534 * (t224 + t301) + t533 * (t299 + t525) + (-t178 + t199 + t291) * t177) * m(3) + (-g(2) * t270 - g(3) * t268 + (t268 ^ 2 + t270 ^ 2) * qJDD(1)) * m(2) + (t60 * (t350 - t404) + t315 + t540 * (t301 + t115) + t537 * (t114 + t299) + (t60 + t522) * t59) * m(6) + (t79 * (t400 - t404) + t316 + t539 * (t301 + t124) + t536 * (t299 + t123) + (t79 + t521) * t78) * m(5) + (-t143 * t404 + t318 + t538 * (t181 + t301) + t535 * (t299 - t546) + (t143 + t528) * t142) * m(4); t317 + (t349 * t484 + t315 + (-t348 + t350) * t60 + t522 * t59 + t540 * t115 + t537 * t114) * m(6) + (t372 * t483 + t316 + (t399 + t400) * t79 + t521 * t78 + t539 * t124 + t536 * t123) * m(5) + (t528 * t142 - t143 * t530 + t538 * t181 - t535 * t546 + t318) * m(4) + (t177 * t199 - t178 * t198 + (-t177 * t305 + t534) * t224 + (t178 * t305 + t533) * t525) * m(3); m(4) * qJDD(3) + m(5) * t44 + m(6) * t10 + (-m(4) - m(5) - m(6)) * g(1); ((t305 * t356 - t51) * t282 + (t305 * t99 - t52) * t281) * t493 + ((t183 * t415 - t340) * t281 + (t544 + (-t326 * t282 + (-t455 + t325) * t281) * qJD(4)) * t282) * t391 + t36 * t394 + t37 * t393 + ((-t308 * t419 + t310 * t418) * t305 + ((t281 * t429 - t282 * t430) * t310 + (-t281 * t431 + t282 * t432) * t308) * qJD(4)) * t494 + ((-t305 * t72 - t19) * t282 + (t305 * t71 - t20) * t281) * t392 + ((-t305 * t70 - t21) * t282 + (t305 * t69 - t22) * t281) * t390 + ((-t414 * t455 - t340) * t282 + (-t544 + (-t325 * t281 + (t183 + t326) * t282) * qJD(4)) * t281) * t389 + t323 + (-t118 * t303 + t190 * t72 + t191 * t71 + t305 * t57 + (-t19 * t282 - t20 * t281) * qJD(4)) * t497 + (t117 * t303 + t190 * t70 + t191 * t69 + t305 * t58 + (-t21 * t282 - t22 * t281) * qJD(4)) * t496 + t365 * t503 + t366 * t502 + (t281 * t356 - t282 * t99) * t495 + (-(-t60 * t397 + ((-t281 * t59 - t282 * t60) * t305 + t53 * (-t281 ^ 2 - t282 ^ 2) * qJD(4)) * t308) * pkin(4) - g(1) * (t223 + t300) - g(3) * (t255 + t175) - g(2) * t386 * t281 + t10 * t135 + t53 * t402 + t17 * t255 + (t10 * t437 + t53 * t481 + t17 * t221 + t59 * t197 + (t139 * t53 + t386 * t60) * t305) * t282 + (t10 * t139 + t53 * t113 + t18 * t386 + t60 * (-pkin(4) * t411 - t197) + (t140 * t53 + t386 * t59) * t305) * t281 + t519) * m(6) + (-g(1) * t269 + g(2) * t188 - g(3) * t189 - t542 * t305 - (t87 * (-t188 * t281 - t189 * t282) + t364 * t269) * qJD(4) + t44 * t354 + t87 * ((t163 * t305 - t110) * t282 + (t111 - t150) * t281) + t364 * t237 + ((t45 - t483) * t282 + (-t305 * t78 - t46) * t281) * t267) * m(5); t323 + (-g(1) * t223 + g(2) * t174 - g(3) * t175 + t10 * (t152 * t282 + t135) + t53 * (-t282 * t94 + t402) + (-t281 * t60 + t282 * t59) * t197 + ((t17 - t484) * t282 + (-t305 * t59 - t18) * t281) * t221 + t519) * m(6);];
tau = t1;
