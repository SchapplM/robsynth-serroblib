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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:31:49
% EndTime: 2019-12-05 18:32:04
% DurationCPUTime: 10.50s
% Computational Cost: add. (19139->627), mult. (11926->765), div. (0->0), fcn. (9239->10), ass. (0->375)
t308 = qJ(1) + qJ(2);
t295 = pkin(9) + t308;
t290 = cos(t295);
t309 = sin(qJ(4));
t461 = t290 * t309;
t267 = rSges(5,2) * t461;
t311 = cos(qJ(4));
t460 = t290 * t311;
t289 = sin(t295);
t507 = rSges(5,3) * t289;
t352 = -rSges(5,1) * t460 - t507;
t167 = -t267 - t352;
t306 = qJD(1) + qJD(2);
t156 = t306 * t167;
t277 = rSges(5,1) * t309 + rSges(5,2) * t311;
t431 = qJD(4) * t289;
t205 = t277 * t431;
t516 = pkin(7) * t289;
t517 = pkin(3) * t290;
t212 = t516 + t517;
t299 = cos(t308);
t457 = t299 * t306;
t425 = pkin(2) * t457;
t379 = -t306 * t212 - t425;
t428 = qJD(4) * t309;
t411 = t289 * t428;
t427 = qJD(4) * t311;
t412 = rSges(5,1) * t411 + (t289 * t427 + t306 * t461) * rSges(5,2);
t591 = -t156 + t205 + t379 - t412;
t310 = sin(qJ(1));
t312 = cos(qJ(1));
t315 = qJD(1) ^ 2;
t343 = (-qJDD(1) * t310 - t312 * t315) * pkin(1);
t590 = g(3) - t343;
t283 = t289 * rSges(4,2);
t512 = rSges(4,1) * t290;
t211 = -t283 + t512;
t469 = t289 * t306;
t258 = rSges(4,2) * t469;
t518 = pkin(2) * t299;
t374 = -t512 - t518;
t589 = -t425 - t258 + (-t211 - t374) * t306;
t282 = qJDD(4) * t290;
t305 = qJD(4) + qJD(5);
t133 = qJDD(5) * t290 - t305 * t469 + t282;
t429 = qJD(4) * t306;
t193 = -t289 * t429 + t282;
t307 = qJ(4) + qJ(5);
t298 = cos(t307);
t288 = t298 * rSges(6,1);
t296 = sin(t307);
t232 = -rSges(6,2) * t296 + t288;
t200 = t232 * t305;
t210 = t290 * t305;
t230 = rSges(6,1) * t296 + rSges(6,2) * t298;
t304 = qJDD(1) + qJDD(2);
t297 = sin(t308);
t519 = pkin(2) * t297;
t471 = t289 * t296;
t556 = rSges(6,2) * t471 + t290 * rSges(6,3);
t377 = t556 - t519;
t455 = t311 * qJD(4) ^ 2;
t303 = t306 ^ 2;
t458 = t299 * t303;
t255 = pkin(4) * t411;
t313 = -pkin(8) - pkin(7);
t456 = t306 * t313;
t437 = t289 * t456 + t255;
t302 = t311 * pkin(4);
t291 = t302 + pkin(3);
t514 = pkin(3) - t291;
t113 = (t290 * t514 + t516) * t306 + t437;
t463 = t290 * t298;
t419 = rSges(6,1) * t463;
t351 = rSges(6,3) * t289 + t419;
t470 = t289 * t298;
t177 = rSges(6,1) * t471 + rSges(6,2) * t470;
t464 = t290 * t296;
t247 = rSges(6,2) * t464;
t415 = t177 * t305 + t306 * t247;
t95 = -t306 * t351 + t415;
t501 = t113 + t95;
t515 = t290 * pkin(7);
t375 = -pkin(3) * t289 + t515;
t559 = -t303 * t212 + t304 * t375;
t513 = pkin(7) + t313;
t560 = t513 * t290;
t395 = -t291 - t288;
t381 = -pkin(3) - t395;
t561 = t381 * t289;
t588 = pkin(2) * t458 + t133 * t230 + t210 * t200 - t501 * t306 - (-t193 * t309 - t290 * t455) * pkin(4) - (t377 - t560 - t561) * t304 - t559 + t590;
t372 = rSges(3,1) * t297 + rSges(3,2) * t299;
t201 = t372 * t306;
t506 = pkin(1) * qJD(1);
t422 = t310 * t506;
t179 = t201 + t422;
t241 = Icges(6,4) * t471;
t153 = -Icges(6,1) * t470 + Icges(6,5) * t290 + t241;
t242 = Icges(6,4) * t464;
t154 = Icges(6,1) * t463 + Icges(6,5) * t289 - t242;
t209 = t289 * t305;
t286 = Icges(6,4) * t298;
t364 = -Icges(6,2) * t296 + t286;
t555 = Icges(6,1) * t296 + t286;
t571 = t555 + t364;
t322 = t209 * (-Icges(6,2) * t463 + t154 - t242) + t210 * (Icges(6,2) * t470 + t153 + t241) + t306 * t571;
t346 = t364 * t289;
t151 = Icges(6,6) * t290 - t346;
t152 = Icges(6,4) * t463 - Icges(6,2) * t464 + Icges(6,6) * t289;
t497 = Icges(6,4) * t296;
t226 = Icges(6,2) * t298 + t497;
t229 = Icges(6,1) * t298 - t497;
t536 = t209 * (t290 * t555 + t152) + t210 * (-t289 * t555 + t151) + t306 * (t226 - t229);
t587 = t322 * t296 + t298 * t536;
t261 = pkin(3) * t469;
t420 = pkin(4) * t428;
t441 = -t290 * t456 - t291 * t469;
t112 = t261 + (-pkin(7) * t306 - t420) * t290 + t441;
t192 = qJDD(4) * t289 + t290 * t429;
t462 = t290 * t306;
t132 = qJD(5) * t462 + qJDD(5) * t289 + t192;
t181 = pkin(7) * t462 - t261;
t520 = pkin(1) * t312;
t521 = pkin(1) * t310;
t376 = -qJDD(1) * t520 + t315 * t521;
t354 = t303 * t519 + t376;
t397 = -t212 - t518;
t243 = t290 * t291;
t141 = t289 * t513 - t243 + t517;
t157 = -t247 + t351;
t450 = t141 - t157;
t357 = t397 + t450;
t423 = rSges(6,1) * t470;
t577 = t210 * t230;
t414 = -t306 * t423 - t577;
t94 = t306 * t556 + t414;
t18 = t132 * t230 + t200 * t209 + (t192 * t309 + t289 * t455) * pkin(4) + (-t112 - t181 - t94) * t306 + t357 * t304 + t354;
t586 = -g(2) + t18;
t408 = t290 * t427;
t409 = t290 * t428;
t467 = t289 * t311;
t424 = rSges(5,1) * t467;
t413 = -rSges(5,1) * t409 - rSges(5,2) * t408 - t306 * t424;
t468 = t289 * t309;
t266 = rSges(5,2) * t468;
t434 = t290 * rSges(5,3) + t266;
t110 = t306 * t434 + t413;
t511 = rSges(5,1) * t311;
t279 = -rSges(5,2) * t309 + t511;
t248 = t279 * qJD(4);
t378 = -t167 + t397;
t46 = t248 * t431 + t192 * t277 + (-t110 - t181) * t306 + t378 * t304 + t354;
t585 = -g(2) + t46;
t398 = -t211 - t518;
t435 = -rSges(4,1) * t469 - rSges(4,2) * t462;
t584 = t304 * t398 - t306 * t435 - g(2) + t354;
t111 = t306 * t352 + t412;
t320 = t343 + (-t297 * t304 - t458) * pkin(2);
t356 = t424 - t434;
t430 = qJD(4) * t290;
t45 = t306 * t111 - t193 * t277 - t248 * t430 - t304 * t356 + t320 + t559;
t582 = -g(3) + t45;
t371 = -rSges(4,1) * t289 - rSges(4,2) * t290;
t581 = -g(3) + t304 * t371 + t306 * (-rSges(4,1) * t462 + t258) + t320;
t233 = rSges(3,1) * t299 - t297 * rSges(3,2);
t580 = t201 * t306 - t233 * t304 - g(2) + t376;
t459 = t297 * t306;
t202 = -rSges(3,1) * t457 + rSges(3,2) * t459;
t579 = t202 * t306 - t304 * t372 - t590;
t475 = t230 * t290;
t578 = t177 * t209 + t210 * t475 + t290 * t94;
t300 = Icges(5,4) * t311;
t365 = -Icges(5,2) * t309 + t300;
t554 = Icges(5,1) * t309 + t300;
t432 = t554 + t365;
t498 = Icges(5,4) * t309;
t271 = Icges(5,2) * t311 + t498;
t274 = Icges(5,1) * t311 - t498;
t433 = t271 - t274;
t576 = (t309 * t432 + t311 * t433) * t306;
t384 = t209 * t230 + t255;
t574 = t306 * t450 + t379 + t384 - t415 - t437;
t573 = t306 * t177 - t210 * t232 - t230 * t469;
t569 = -t289 * t514 + t560;
t421 = t312 * t506;
t60 = t306 * t357 + t384 - t421;
t568 = (t200 * t289 - t209 * t232 + t230 * t462 - t306 * t475) * t60;
t150 = Icges(6,5) * t463 - Icges(6,6) * t464 + Icges(6,3) * t289;
t64 = t290 * t150 + t152 * t471 - t154 * t470;
t360 = t226 * t296 - t298 * t555;
t224 = Icges(6,5) * t296 + Icges(6,6) * t298;
t480 = t224 * t290;
t97 = t289 * t360 + t480;
t567 = t209 * t64 + t97 * t306;
t163 = Icges(5,4) * t460 - Icges(5,2) * t461 + Icges(5,6) * t289;
t265 = Icges(5,4) * t461;
t165 = Icges(5,1) * t460 + Icges(5,5) * t289 - t265;
t361 = t163 * t309 - t165 * t311;
t565 = t290 * t361;
t363 = t152 * t296 - t154 * t298;
t562 = t363 * t290;
t557 = -t243 - t518;
t477 = t226 * t305;
t552 = -Icges(6,6) * t306 + t477;
t155 = t306 * t356;
t426 = pkin(2) * t459;
t380 = -t306 * t375 + t426;
t328 = t277 * t430 + t155 + t380;
t78 = t328 + t422;
t79 = t306 * t378 + t205 - t421;
t551 = t289 * t79 + t290 * t78;
t545 = -Icges(5,6) * t306 + qJD(4) * t271;
t107 = t289 * t545 - t365 * t462;
t349 = t306 * t274;
t543 = -Icges(5,5) * t306 + qJD(4) * t554;
t109 = t289 * t543 - t290 * t349;
t270 = Icges(5,5) * t311 - Icges(5,6) * t309;
t345 = t270 * t289;
t160 = Icges(5,3) * t290 - t345;
t347 = t365 * t289;
t162 = Icges(5,6) * t290 - t347;
t264 = Icges(5,4) * t468;
t164 = -Icges(5,1) * t467 + Icges(5,5) * t290 + t264;
t99 = t162 * t311 + t164 * t309;
t550 = qJD(4) * t99 + t107 * t309 - t109 * t311 - t160 * t306;
t549 = -Icges(6,3) * t306 + t224 * t305;
t548 = -Icges(6,5) * t306 + t305 * t555;
t100 = t163 * t311 + t165 * t309;
t106 = -t290 * t545 - t306 * t347;
t108 = -t289 * t349 - t290 * t543;
t161 = Icges(5,5) * t460 - Icges(5,6) * t461 + Icges(5,3) * t289;
t547 = qJD(4) * t100 + t106 * t309 - t108 * t311 - t161 * t306;
t269 = Icges(5,5) * t309 + Icges(5,6) * t311;
t546 = -Icges(5,3) * t306 + qJD(4) * t269;
t237 = t365 * qJD(4);
t238 = t274 * qJD(4);
t359 = t271 * t311 + t309 * t554;
t544 = qJD(4) * t359 + t237 * t309 - t238 * t311 - t269 * t306;
t442 = -Icges(5,2) * t460 + t165 - t265;
t444 = t290 * t554 + t163;
t541 = t309 * t442 + t311 * t444;
t443 = Icges(5,2) * t467 + t164 + t264;
t445 = -t289 * t554 + t162;
t540 = -t309 * t443 - t311 * t445;
t385 = -t229 * t305 + t477;
t386 = t571 * t305;
t539 = -t224 * t306 + t296 * t386 + t298 * t385;
t391 = t154 * t305 - t290 * t552 - t306 * t346;
t348 = t229 * t306;
t393 = t152 * t305 + t289 * t348 + t290 * t548;
t538 = -t150 * t306 + t296 * t391 + t298 * t393;
t225 = Icges(6,5) * t298 - Icges(6,6) * t296;
t149 = Icges(6,3) * t290 - t225 * t289;
t392 = t153 * t305 + t289 * t552 - t364 * t462;
t394 = t151 * t305 - t289 * t548 + t290 * t348;
t537 = -t149 * t306 + t296 * t392 + t298 * t394;
t535 = t132 / 0.2e1;
t534 = t133 / 0.2e1;
t533 = t192 / 0.2e1;
t532 = t193 / 0.2e1;
t531 = -t209 / 0.2e1;
t530 = t209 / 0.2e1;
t529 = -t210 / 0.2e1;
t528 = t210 / 0.2e1;
t527 = t289 / 0.2e1;
t526 = t290 / 0.2e1;
t525 = t304 / 0.2e1;
t524 = -t306 / 0.2e1;
t523 = t306 / 0.2e1;
t522 = -rSges(5,3) - pkin(7);
t171 = t224 * t289;
t98 = -t290 * t360 + t171;
t502 = t98 * t306;
t184 = t269 * t289;
t358 = t271 * t309 - t311 * t554;
t118 = -t290 * t358 + t184;
t488 = t118 * t306;
t485 = t151 * t296;
t484 = t153 * t298;
t483 = t162 * t309;
t482 = t164 * t311;
t478 = t225 * t306;
t476 = t230 * t289;
t474 = t269 * t290;
t473 = t270 * t306;
t190 = t277 * t289;
t472 = t277 * t290;
t466 = t290 * t112;
t465 = t290 * t141;
t454 = t290 * t149 + t151 * t471;
t453 = -t289 * t149 - t153 * t463;
t452 = t290 * t160 + t162 * t468;
t451 = t289 * t160 + t164 * t460;
t268 = pkin(4) * t468;
t407 = -t469 / 0.2e1;
t406 = t462 / 0.2e1;
t405 = -pkin(3) - t511;
t402 = -t431 / 0.2e1;
t401 = t431 / 0.2e1;
t400 = -t430 / 0.2e1;
t399 = t430 / 0.2e1;
t388 = -t161 - t482;
t383 = pkin(4) * t409;
t382 = t261 - t413;
t280 = rSges(2,1) * t312 - t310 * rSges(2,2);
t373 = rSges(2,1) * t310 + rSges(2,2) * t312;
t69 = -t164 * t467 + t452;
t70 = t290 * t161 + t163 * t468 - t165 * t467;
t370 = t289 * t70 + t290 * t69;
t71 = -t162 * t461 + t451;
t72 = t161 * t289 - t565;
t369 = t289 * t72 + t290 * t71;
t84 = t152 * t298 + t154 * t296;
t362 = -t482 + t483;
t355 = t423 - t556;
t65 = -t151 * t464 - t453;
t183 = t283 + t374;
t353 = t560 - t556;
t180 = -t233 * t306 - t421;
t342 = t422 + t426;
t182 = t371 - t519;
t337 = t171 * t210 - t209 * t480 + t478;
t334 = -t289 * t478 - t290 * t549 + t306 * t363;
t333 = -t290 * t478 + t549 * t289 + (-t484 + t485) * t306;
t332 = -t290 * t546 + (-t345 + t361) * t306;
t331 = -t270 * t462 + t289 * t546 + t306 * t362;
t330 = t225 * t305 + t306 * t360;
t329 = t270 * qJD(4) + t306 * t358;
t327 = t383 - t414 - t441;
t325 = t290 * t167 + t289 * t356;
t13 = t333 * t289 - t290 * t537;
t14 = t334 * t289 - t290 * t538;
t15 = t289 * t537 + t333 * t290;
t16 = t289 * t538 + t334 * t290;
t63 = -t153 * t470 + t454;
t30 = t210 * t63 + t567;
t66 = t150 * t289 - t562;
t31 = t209 * t66 + t210 * t65 + t502;
t40 = -t296 * t394 + t298 * t392;
t41 = -t296 * t393 + t298 * t391;
t47 = t330 * t289 - t290 * t539;
t48 = t289 * t539 + t330 * t290;
t83 = t151 * t298 + t153 * t296;
t324 = (t13 * t210 + t132 * t66 + t133 * t65 + t14 * t209 + t304 * t98 + t306 * t47) * t527 + (t337 * t289 - t290 * t587) * t531 + (t587 * t289 + t337 * t290) * t529 + (t132 * t64 + t133 * t63 + t15 * t210 + t16 * t209 + t304 * t97 + t306 * t48) * t526 + (-t296 * t536 + t298 * t322) * t524 + t30 * t407 + t31 * t406 + ((t306 * t66 + t13) * t290 + (-t306 * t65 + t14) * t289) * t530 + (t289 * t66 + t290 * t65) * t535 + (t289 * t64 + t290 * t63) * t534 + ((t306 * t64 + t15) * t290 + (-t306 * t63 + t16) * t289) * t528 + (t289 * t84 + t290 * t83) * t525 + ((t306 * t84 + t40) * t290 + (-t306 * t83 + t41) * t289) * t523;
t115 = -t419 + t247 + (-rSges(6,3) + t313) * t289 + t557;
t123 = t289 * t405 + t434 + t515 - t519;
t124 = t289 * t522 + t290 * t405 + t267 - t518;
t114 = t289 * t395 - t290 * t313 + t377;
t321 = t577 + t380 + t383 + (t569 + t355) * t306;
t319 = t562 + (-t150 - t484) * t289 + t454;
t59 = t321 + t422;
t318 = (-t60 * t377 - t59 * (-t351 + t557)) * t306;
t117 = t289 * t358 + t474;
t116 = t117 * t306;
t36 = qJD(4) * t370 + t116;
t37 = qJD(4) * t369 + t488;
t51 = -qJD(4) * t362 + t107 * t311 + t109 * t309;
t52 = -qJD(4) * t361 + t106 * t311 + t108 * t309;
t57 = t329 * t289 - t290 * t544;
t58 = t289 * t544 + t329 * t290;
t317 = (t116 + ((t452 + t72 + t565) * t290 + (-t71 + (t388 - t483) * t290 + t70 + t451) * t289) * qJD(4)) * t402 + ((t66 + t319) * t210 + t567) * t531 + (t84 + t98) * t535 + (t83 + t97) * t534 + (t118 + t100) * t533 + (t117 + t99) * t532 + (-t502 + (t64 + (-t150 + t485) * t290 - t363 * t289 + t453) * t210 + (-t63 + t319) * t209 + t31) * t529 + (t40 + t48) * t528 + (-t488 + ((t70 + (-t161 + t483) * t290 - t451) * t290 + (t289 * t388 + t452 - t69) * t289) * qJD(4) + t37) * t400 + (t51 + t58) * t399 + (-qJD(4) * t358 + t237 * t311 + t238 * t309 - t296 * t385 + t298 * t386) * t306 + (t41 + t47 + t30) * t530 + (t52 + t57 + t36) * t401 + (t226 * t298 + t296 * t555 + Icges(3,3) + Icges(4,3) + t359) * t304;
t316 = (t79 * (-t266 + t519) - t78 * (-t507 - t516 - t518) + (-t405 * t78 + t522 * t79) * t290) * t306;
t198 = t306 * t371;
t148 = t306 * t398 - t421;
t147 = -t198 + t342;
t137 = t290 * t157;
t87 = qJD(4) * t325 + qJD(3);
t53 = qJD(3) + t210 * t157 + t209 * t355 + (t289 * t569 - t465) * qJD(4);
t44 = qJDD(3) + t193 * t167 + t192 * t356 + (t290 * t110 - t289 * t111) * qJD(4);
t22 = t289 * t547 + t332 * t290;
t21 = t289 * t550 + t331 * t290;
t20 = t332 * t289 - t290 * t547;
t19 = t331 * t289 - t290 * t550;
t10 = qJDD(3) - t193 * t141 + t192 * t569 + t133 * t157 + t210 * t94 + t132 * t355 - t209 * t95 + (-t289 * t113 + t466) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t317 + (t580 * (-t233 - t520) + t579 * (-t372 - t521) + (t180 - t202 + t421) * t179) * m(3) + ((qJDD(1) * t373 + g(3)) * t373 + (qJDD(1) * t280 + g(2)) * t280) * m(2) + (t60 * (t327 + t422) + t318 + t586 * (t115 - t520) - t588 * (t114 - t521) + (-t60 + t574) * t59) * m(6) + (t79 * (t382 + t422) + t316 + (-t79 + t591) * t78 + t585 * (t124 - t520) + t582 * (t123 - t521)) * m(5) + (t148 * (t342 - t435) + t584 * (t183 - t520) + t581 * (t182 - t521) + (-t148 + t589) * t147) * m(4); t317 + (t318 + (-t321 + t327) * t60 + t574 * t59 + t586 * t115 - t588 * t114) * m(6) + (t316 + (-t328 + t382) * t79 + t591 * t78 + t585 * t124 + t582 * t123) * m(5) + (t584 * t183 + t581 * t182 + t589 * t147 + (t198 - t435) * t148) * m(4) + (-t179 * t202 + t180 * t201 + (-t179 * t306 - t580) * t233 - (t180 * t306 + t579) * t372) * m(3); m(4) * qJDD(3) + m(5) * t44 + m(6) * t10 + (-m(4) - m(5) - m(6)) * g(1); ((t184 * t430 + t473) * t290 + (t576 + (t541 * t289 + (-t474 - t540) * t290) * qJD(4)) * t289) * t400 + t324 + t36 * t407 + t37 * t406 + ((-t309 * t433 + t311 * t432) * t306 + ((t289 * t442 + t290 * t443) * t311 + (-t289 * t444 - t290 * t445) * t309) * qJD(4)) * t524 + ((t100 * t306 + t51) * t290 + (-t306 * t99 + t52) * t289) * t523 + ((-t431 * t474 + t473) * t289 + (-t576 + (t540 * t290 + (t184 - t541) * t289) * qJD(4)) * t290) * t402 + (t118 * t304 + t192 * t72 + t193 * t71 + t306 * t57 + (t19 * t290 + t20 * t289) * qJD(4)) * t527 + (t117 * t304 + t192 * t70 + t193 * t69 + t306 * t58 + (t21 * t290 + t22 * t289) * qJD(4)) * t526 + ((t306 * t72 + t19) * t290 + (-t306 * t71 + t20) * t289) * t401 + ((t306 * t70 + t21) * t290 + (-t306 * t69 + t22) * t289) * t399 + t369 * t533 + t370 * t532 + (t100 * t289 + t290 * t99) * t525 + (-g(1) * (t232 + t302) - g(2) * (t268 + t177) + t10 * (-t465 + t137 + (t353 + t561) * t289) + t18 * (t268 + t476) - t588 * t290 * (-pkin(4) * t309 - t230) + (t466 - t501 * t289 + (t353 * t290 + (t290 * t381 + t450) * t289) * t306 - (-t289 ^ 2 - t290 ^ 2) * t420 + t578) * t53 + t568 + (-(-pkin(4) * t427 - t200) * t290 - pkin(4) * t408 + t573) * t59) * m(6) + (t44 * t325 + t87 * ((-t111 - t156) * t289 + (t110 + t155) * t290) + t46 * t190 - t45 * t472 + (t79 * t462 - t78 * t469) * t277 + t551 * t248 - (-t190 * t78 + t472 * t79) * t306 - (t87 * (-t190 * t289 - t290 * t472) + t551 * t279) * qJD(4) - g(1) * t279 - g(2) * t190 + g(3) * t472) * m(5); t324 + (t10 * (t289 * t355 + t137) + t18 * t476 - g(1) * t232 - g(2) * t177 + (-t289 * t95 + (-t289 * t157 + t290 * t355) * t306 + t578) * t53 + t568 + (t200 * t290 + t573) * t59 + t588 * t475) * m(6);];
tau = t1;
