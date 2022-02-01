% Calculate vector of inverse dynamics joint torques for
% S5RRPRR6
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
% m [6x1]
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:51
% EndTime: 2022-01-20 11:17:16
% DurationCPUTime: 17.07s
% Computational Cost: add. (26150->759), mult. (25223->1031), div. (0->0), fcn. (24420->10), ass. (0->370)
t371 = qJ(4) + qJ(5);
t364 = cos(t371);
t372 = qJ(1) + qJ(2);
t363 = sin(t372);
t369 = qJD(4) + qJD(5);
t370 = qJD(1) + qJD(2);
t374 = cos(pkin(9));
t582 = -t369 * t374 + t370;
t431 = t582 * t363;
t455 = t370 * t374 - t369;
t362 = sin(t371);
t365 = cos(t372);
t516 = t365 * t362;
t175 = t364 * t431 - t455 * t516;
t517 = t364 * t365;
t176 = t362 * t431 + t455 * t517;
t435 = rSges(6,1) * t176 + rSges(6,2) * t175;
t373 = sin(pkin(9));
t510 = t370 * t373;
t473 = t365 * t510;
t104 = rSges(6,3) * t473 + t435;
t375 = sin(qJ(4));
t519 = t363 * t375;
t346 = pkin(4) * t519;
t508 = t374 * t375;
t550 = pkin(4) * qJD(4);
t450 = t550 * t508;
t377 = cos(qJ(4));
t379 = -pkin(8) - pkin(7);
t509 = t373 * t379;
t597 = t370 * t509 + t377 * t550;
t469 = -t363 * t450 - t365 * t597;
t359 = pkin(4) * t377 + pkin(3);
t553 = pkin(3) - t359;
t554 = pkin(7) * t373;
t573 = t374 * t553 + t554;
t125 = (-t365 * t573 + t346) * t370 + t469;
t520 = t363 * t374;
t264 = t362 * t520 + t517;
t265 = t364 * t520 - t516;
t434 = t265 * rSges(6,1) - t264 * rSges(6,2);
t521 = t363 * t373;
t171 = -rSges(6,3) * t521 - t434;
t512 = t365 * t375;
t484 = pkin(4) * t512;
t442 = -t359 * t520 + t484;
t552 = pkin(7) + t379;
t556 = pkin(3) * t374;
t188 = (t373 * t552 + t556) * t363 + t442;
t486 = qJD(4) * t373;
t465 = t370 * t486;
t485 = qJDD(4) * t373;
t261 = t363 * t485 + t365 * t465;
t515 = t365 * t370;
t207 = (qJD(5) * t515 + qJDD(5) * t363) * t373 + t261;
t433 = -rSges(6,1) * t362 - rSges(6,2) * t364;
t293 = t433 * t373;
t239 = t369 * t293;
t243 = -t373 * t553 + t374 * t552;
t258 = -rSges(6,3) * t374 + (rSges(6,1) * t364 - rSges(6,2) * t362) * t373;
t452 = t369 * t373;
t284 = t363 * t452;
t367 = qJDD(1) + qJDD(2);
t308 = (-qJDD(4) - qJDD(5)) * t374 + t367;
t343 = -qJDD(4) * t374 + t367;
t348 = -qJD(4) * t374 + t370;
t316 = t365 * pkin(2) + t363 * qJ(3);
t352 = qJD(3) * t365;
t235 = t316 * t370 - t352;
t376 = sin(qJ(1));
t378 = cos(qJ(1));
t381 = qJD(1) ^ 2;
t410 = (-qJDD(1) * t376 - t378 * t381) * pkin(1);
t487 = qJD(3) * t370;
t393 = qJDD(3) * t363 + t365 * t487 + t410;
t439 = t554 + t556;
t294 = t439 * t363;
t354 = t365 * qJ(3);
t314 = pkin(2) * t363 - t354;
t492 = -t314 - t294;
t386 = (-t439 * t515 - t235) * t370 + t492 * t367 + t393;
t555 = pkin(4) * t375;
t453 = t373 ^ 2 * qJD(4) ^ 2 * t555;
t20 = -t104 * t582 - t348 * t125 + t171 * t308 + t343 * t188 + t207 * t258 + t284 * t239 + t261 * t243 - t363 * t453 + t386;
t616 = t20 - g(1);
t507 = t374 * t377;
t290 = t363 * t507 - t512;
t291 = t363 * t377 - t365 * t508;
t205 = -qJD(4) * t290 + t291 * t370;
t289 = t363 * t508 + t365 * t377;
t292 = t365 * t507 + t519;
t206 = -qJD(4) * t289 + t292 * t370;
t437 = rSges(5,1) * t206 + rSges(5,2) * t205;
t117 = rSges(5,3) * t473 + t437;
t436 = t290 * rSges(5,1) - t289 * rSges(5,2);
t191 = -rSges(5,3) * t521 - t436;
t283 = -rSges(5,3) * t374 + (rSges(5,1) * t377 - rSges(5,2) * t375) * t373;
t307 = (-rSges(5,1) * t375 - rSges(5,2) * t377) * t373;
t299 = qJD(4) * t307;
t467 = t363 * t486;
t48 = -t117 * t348 + t191 * t343 + t261 * t283 + t299 * t467 + t386;
t615 = t48 - g(1);
t614 = t171 * t582 + t348 * t188 + t243 * t467 + t258 * t284;
t613 = t348 * t191 + t283 * t467;
t179 = Icges(5,5) * t290 - Icges(5,6) * t289 + Icges(5,3) * t521;
t270 = Icges(5,4) * t290;
t182 = -Icges(5,2) * t289 + Icges(5,6) * t521 + t270;
t269 = Icges(5,4) * t289;
t186 = -Icges(5,1) * t290 - Icges(5,5) * t521 + t269;
t82 = -t179 * t374 - (t182 * t375 + t186 * t377) * t373;
t159 = Icges(6,5) * t265 - Icges(6,6) * t264 + Icges(6,3) * t521;
t245 = Icges(6,4) * t265;
t162 = -Icges(6,2) * t264 + Icges(6,6) * t521 + t245;
t244 = Icges(6,4) * t264;
t166 = -Icges(6,1) * t265 - Icges(6,5) * t521 + t244;
t75 = -t159 * t374 - (t162 * t362 + t166 * t364) * t373;
t610 = t354 - t436;
t609 = t354 + t442 - t434;
t303 = t370 * t314;
t351 = qJD(3) * t363;
t331 = qJ(3) * t515;
t490 = t331 + t351;
t608 = t490 - t351 + t303;
t255 = -Icges(6,3) * t374 + (Icges(6,5) * t364 - Icges(6,6) * t362) * t373;
t532 = Icges(6,4) * t364;
t256 = -Icges(6,6) * t374 + (-Icges(6,2) * t362 + t532) * t373;
t533 = Icges(6,4) * t362;
t257 = -Icges(6,5) * t374 + (Icges(6,1) * t364 - t533) * t373;
t513 = t365 * t374;
t523 = t363 * t364;
t266 = -t362 * t513 + t523;
t524 = t362 * t363;
t267 = t364 * t513 + t524;
t514 = t365 * t373;
t106 = t255 * t514 + t256 * t266 + t257 * t267;
t285 = t365 * t452;
t69 = t159 * t514 + t266 * t162 - t166 * t267;
t161 = Icges(6,5) * t267 + Icges(6,6) * t266 + Icges(6,3) * t514;
t534 = Icges(6,4) * t267;
t164 = Icges(6,2) * t266 + Icges(6,6) * t514 + t534;
t246 = Icges(6,4) * t266;
t167 = Icges(6,1) * t267 + Icges(6,5) * t514 + t246;
t70 = t161 * t514 + t266 * t164 + t267 * t167;
t29 = t106 * t582 + t284 * t69 + t285 * t70;
t606 = -t182 * t289 - t186 * t290;
t605 = t291 * t182 - t186 * t292;
t67 = t159 * t521 - t162 * t264 - t166 * t265;
t599 = t179 * t514;
t200 = t266 * rSges(6,1) - t267 * rSges(6,2);
t229 = t258 * t521;
t598 = -t200 * t582 + t370 * t229 + t285 * t293;
t199 = -t264 * rSges(6,1) - t265 * rSges(6,2);
t596 = t374 * t104 + t199 * t582 + t239 * t521 + t258 * t473 - t284 * t293;
t551 = pkin(1) * qJD(1);
t480 = t376 * t551;
t315 = rSges(3,1) * t363 + rSges(3,2) * t365;
t525 = t315 * t370;
t259 = -t480 - t525;
t268 = t370 * t294;
t595 = t268 + t608;
t323 = rSges(4,2) * t473;
t483 = rSges(4,1) * t520;
t489 = rSges(4,2) * t521 + t365 * rSges(4,3);
t241 = t483 - t489;
t493 = -t314 - t241;
t583 = -rSges(4,1) * t513 - t363 * rSges(4,3);
t87 = t493 * t367 + (t370 * t583 - t235 + t323) * t370 + t393;
t594 = -g(1) + t87;
t474 = t363 * t510;
t430 = t582 * t365;
t173 = t364 * t430 + t455 * t524;
t174 = t362 * t430 - t455 * t523;
t506 = t174 * rSges(6,1) + t173 * rSges(6,2);
t103 = -rSges(6,3) * t474 + t506;
t401 = t363 * t597 - t365 * t450 + t370 * t484;
t522 = t363 * t370;
t124 = t522 * t573 + t401;
t172 = t267 * rSges(6,1) + t266 * rSges(6,2) + rSges(6,3) * t514;
t295 = pkin(3) * t513 + pkin(7) * t514;
t422 = t359 * t513 - t365 * t509 + t346;
t189 = t422 - t295;
t335 = t365 * t485;
t208 = t335 + (qJDD(5) * t365 - t369 * t522) * t373;
t262 = -t363 * t465 + t335;
t366 = t378 * pkin(1);
t557 = pkin(1) * t376;
t443 = qJDD(1) * t366 - t381 * t557;
t407 = t370 * (-pkin(2) * t522 + t490) + t367 * t316 + t363 * t487 + t443;
t392 = -t294 * t370 ^ 2 + t367 * t295 + t407;
t21 = t103 * t582 + t124 * t348 + t172 * t308 + t189 * t343 - t208 * t258 - t239 * t285 - t243 * t262 + (-qJDD(3) + t453) * t365 + t392;
t593 = -g(2) + t21;
t203 = -qJD(4) * t292 + t289 * t370;
t204 = qJD(4) * t291 - t290 * t370;
t500 = t204 * rSges(5,1) + t203 * rSges(5,2);
t116 = -rSges(5,3) * t474 + t500;
t192 = t292 * rSges(5,1) + t291 * rSges(5,2) + rSges(5,3) * t514;
t49 = t116 * t348 + t192 * t343 - t262 * t283 + (-t299 * t486 - qJDD(3)) * t365 + t392;
t592 = -g(2) + t49;
t242 = -rSges(4,2) * t514 - t583;
t491 = rSges(4,2) * t474 + rSges(4,3) * t515;
t88 = -qJDD(3) * t365 + t367 * t242 + t370 * (-t370 * t483 + t491) + t407;
t591 = -g(2) + t88;
t282 = rSges(3,1) * t515 - rSges(3,2) * t522;
t590 = -t282 * t370 - t315 * t367 - g(1) + t410;
t317 = t365 * rSges(3,1) - rSges(3,2) * t363;
t589 = t317 * t367 - t370 * t525 - g(2) + t443;
t73 = t599 + t605;
t587 = t73 - t599;
t222 = t242 + t316;
t586 = t222 * t370;
t494 = t295 + t316;
t585 = t370 * t494;
t584 = t370 * t241 + t491;
t581 = t401 + t506 + t595 - t614;
t580 = t331 + t303;
t466 = t365 * t486;
t579 = t192 * t348 - t283 * t466;
t409 = -t556 - pkin(2) + (-rSges(5,3) - pkin(7)) * t373;
t441 = t351 - t480;
t400 = t370 * t492 + t441;
t89 = t400 + t613;
t479 = t378 * t551;
t440 = -t352 + t479;
t399 = t440 + t585;
t90 = t399 + t579;
t577 = (t89 * t409 * t365 + (-t89 * qJ(3) + t90 * (-rSges(5,3) * t373 - pkin(2) - t439)) * t363) * t370 + t615 * t363 * t409;
t421 = -rSges(6,3) * t373 - t359 * t374 - pkin(2);
t62 = t400 + t614;
t572 = t172 * t582 + t189 * t348 - t243 * t466 - t285 * t258;
t63 = t399 + t572;
t576 = (t62 * t421 * t365 + (t62 * (-qJ(3) - t555) + t63 * t421) * t363) * t370 + t616 * t363 * (-pkin(2) + (-rSges(6,3) + t379) * t373);
t68 = t161 * t521 - t264 * t164 + t265 * t167;
t181 = Icges(5,5) * t292 + Icges(5,6) * t291 + Icges(5,3) * t514;
t537 = Icges(5,4) * t292;
t184 = Icges(5,2) * t291 + Icges(5,6) * t514 + t537;
t271 = Icges(5,4) * t291;
t187 = Icges(5,1) * t292 + Icges(5,5) * t514 + t271;
t72 = t181 * t521 - t289 * t184 + t290 * t187;
t574 = t500 - t613;
t571 = t363 * (-Icges(5,2) * t290 - t186 - t269) + t365 * (-Icges(5,2) * t292 + t187 + t271);
t287 = (-Icges(6,2) * t364 - t533) * t373;
t570 = t284 * (-Icges(6,2) * t265 - t166 - t244) + t285 * (-Icges(6,2) * t267 + t167 + t246) + t582 * (t257 + t287);
t569 = t207 / 0.2e1;
t568 = t208 / 0.2e1;
t567 = t261 / 0.2e1;
t566 = t262 / 0.2e1;
t565 = -t284 / 0.2e1;
t564 = t284 / 0.2e1;
t562 = t285 / 0.2e1;
t560 = t363 / 0.2e1;
t559 = -t365 / 0.2e1;
t558 = -t374 / 0.2e1;
t100 = Icges(6,1) * t176 + Icges(6,4) * t175 + Icges(6,5) * t473;
t96 = Icges(6,5) * t176 + Icges(6,6) * t175 + Icges(6,3) * t473;
t98 = Icges(6,4) * t176 + Icges(6,2) * t175 + Icges(6,6) * t473;
t36 = -t374 * t96 + ((-t162 * t369 + t100) * t364 + (t166 * t369 - t98) * t362) * t373;
t547 = t36 * t284;
t71 = t179 * t521 + t606;
t546 = t363 * t71;
t95 = Icges(6,5) * t174 + Icges(6,6) * t173 - Icges(6,3) * t474;
t97 = Icges(6,4) * t174 + Icges(6,2) * t173 - Icges(6,6) * t474;
t99 = Icges(6,1) * t174 + Icges(6,4) * t173 - Icges(6,5) * t474;
t37 = -t374 * t95 + ((-t164 * t369 + t99) * t364 + (-t167 * t369 - t97) * t362) * t373;
t545 = t37 * t285;
t544 = t63 * t239;
t543 = t75 * t207;
t76 = -t161 * t374 + (-t164 * t362 + t167 * t364) * t373;
t542 = t76 * t208;
t541 = t82 * t261;
t83 = -t181 * t374 + (-t184 * t375 + t187 * t377) * t373;
t540 = t83 * t262;
t123 = -t255 * t374 + (-t256 * t362 + t257 * t364) * t373;
t286 = (-Icges(6,5) * t362 - Icges(6,6) * t364) * t373;
t236 = t369 * t286;
t237 = t369 * t287;
t288 = (-Icges(6,1) * t362 - t532) * t373;
t238 = t369 * t288;
t84 = -t236 * t374 + ((-t256 * t369 + t238) * t364 + (-t257 * t369 - t237) * t362) * t373;
t539 = t123 * t308 + t582 * t84;
t535 = Icges(5,4) * t377;
t279 = -Icges(5,6) * t374 + (-Icges(5,2) * t375 + t535) * t373;
t536 = Icges(5,4) * t375;
t280 = -Icges(5,5) * t374 + (Icges(5,1) * t377 - t536) * t373;
t304 = (-Icges(5,5) * t375 - Icges(5,6) * t377) * t373;
t296 = qJD(4) * t304;
t305 = (-Icges(5,2) * t377 - t536) * t373;
t297 = qJD(4) * t305;
t306 = (-Icges(5,1) * t375 - t535) * t373;
t298 = qJD(4) * t306;
t107 = -t296 * t374 + (-t297 * t375 + t298 * t377 + (-t279 * t377 - t280 * t375) * qJD(4)) * t373;
t278 = -Icges(5,3) * t374 + (Icges(5,5) * t377 - Icges(5,6) * t375) * t373;
t135 = -t278 * t374 + (-t279 * t375 + t280 * t377) * t373;
t538 = t107 * t348 + t135 * t343;
t118 = t278 * t521 - t279 * t289 + t280 * t290;
t530 = t118 * t348;
t503 = -t172 - t189;
t496 = -t279 + t306;
t495 = t280 + t305;
t74 = t181 * t514 + t291 * t184 + t292 * t187;
t464 = t521 / 0.2e1;
t463 = t514 / 0.2e1;
t462 = -rSges(4,1) * t374 - pkin(2);
t461 = -t486 / 0.2e1;
t460 = t486 / 0.2e1;
t459 = t171 * t370 - t103;
t458 = t285 * t199 - t200 * t284;
t449 = -t474 / 0.2e1;
t448 = t370 * t463;
t447 = t363 * t461;
t446 = t363 * t460;
t445 = t365 * t461;
t444 = t365 * t460;
t338 = rSges(2,1) * t378 - rSges(2,2) * t376;
t337 = rSges(2,1) * t376 + rSges(2,2) * t378;
t432 = t363 * t89 - t365 * t90;
t429 = t352 - t585;
t427 = -t191 * t365 - t192 * t363;
t426 = (-Icges(5,5) * t289 - Icges(5,6) * t290) * t363 + (Icges(5,5) * t291 - Icges(5,6) * t292) * t365;
t273 = t291 * pkin(4);
t420 = t352 - t437;
t413 = (t365 * t72 + t546) * t373;
t412 = (t363 * t73 + t365 * t74) * t373;
t411 = t373 * t426;
t134 = t192 + t494;
t272 = t289 * pkin(4);
t405 = (-Icges(6,5) * t264 - Icges(6,6) * t265) * t284 + (Icges(6,5) * t266 - Icges(6,6) * t267) * t285 + t286 * t582;
t108 = t427 * t486;
t403 = (Icges(5,1) * t291 - t184 - t537) * t365 + (-Icges(5,1) * t289 - t182 - t270) * t363;
t221 = t363 * t462 + t354 + t489;
t14 = -t103 * t284 + t104 * t285 - t171 * t208 - t172 * t207 - t188 * t262 - t189 * t261 + (-t124 * t363 + t125 * t365) * t486;
t57 = -t171 * t285 - t172 * t284 + (-t188 * t365 - t189 * t363) * t486;
t397 = t20 * (-t171 * t374 + t229) + (t104 * t57 - t14 * t171) * t514;
t105 = t255 * t521 - t256 * t264 + t257 * t265;
t396 = t352 - t435 - t469;
t394 = t373 * t405;
t122 = t422 + t172 + t316;
t16 = t100 * t267 + t162 * t173 - t166 * t174 + t266 * t98 + (-t159 * t522 + t365 * t96) * t373;
t17 = t164 * t173 + t167 * t174 + t266 * t97 + t267 * t99 + (-t161 * t522 + t365 * t95) * t373;
t18 = t100 * t265 + t162 * t175 - t166 * t176 - t264 * t98 + (t159 * t515 + t363 * t96) * t373;
t19 = t164 * t175 + t167 * t176 - t264 * t97 + t265 * t99 + (t161 * t515 + t363 * t95) * t373;
t389 = (Icges(6,1) * t266 - t164 - t534) * t285 + (-Icges(6,1) * t264 - t162 - t245) * t284 + (-t256 + t288) * t582;
t53 = t173 * t256 + t174 * t257 + t237 * t266 + t238 * t267 + (t236 * t365 - t255 * t522) * t373;
t54 = t175 * t256 + t176 * t257 - t237 * t264 + t238 * t265 + (t236 * t363 + t255 * t515) * t373;
t390 = (t105 * t308 + t18 * t284 + t19 * t285 + t207 * t67 + t208 * t68 + t54 * t582) * t464 + (-t264 * t570 + t265 * t389 + t363 * t394) * t565 - (t266 * t570 + t389 * t267 + t365 * t394) * t285 / 0.2e1 - (-t405 * t374 + (-t362 * t570 + t364 * t389) * t373) * t582 / 0.2e1 + (t106 * t308 + t16 * t284 + t17 * t285 + t207 * t69 + t208 * t70 + t53 * t582) * t463 + t29 * t449 + (t105 * t582 + t284 * t67 + t285 * t68) * t448 + (-t105 * t374 + (t363 * t67 + t365 * t68) * t373) * t569 + (-t106 * t374 + (t363 * t69 + t365 * t70) * t373) * t568 + (-t374 * t54 + ((t370 * t67 + t19) * t365 + (-t370 * t68 + t18) * t363) * t373) * t564 + (-t374 * t53 + ((t370 * t69 + t17) * t365 + (-t370 * t70 + t16) * t363) * t373) * t562 + t308 * (-t123 * t374 + (t363 * t75 + t365 * t76) * t373) / 0.2e1 + (t539 + t542 + t543 + t545 + t547) * t558 + t582 * (-t374 * t84 + ((t370 * t75 + t37) * t365 + (-t370 * t76 + t36) * t363) * t373) / 0.2e1;
t168 = t370 * t493 + t441;
t169 = t440 + t586;
t385 = (t168 * t462 * t365 + (t168 * (-rSges(4,3) - qJ(3)) + t169 * t462) * t363) * t370;
t119 = t278 * t514 + t279 * t291 + t280 * t292;
t109 = t119 * t348;
t41 = qJD(4) * t413 + t530;
t42 = qJD(4) * t412 + t109;
t111 = Icges(5,5) * t206 + Icges(5,6) * t205 + Icges(5,3) * t473;
t113 = Icges(5,4) * t206 + Icges(5,2) * t205 + Icges(5,6) * t473;
t115 = Icges(5,1) * t206 + Icges(5,4) * t205 + Icges(5,5) * t473;
t46 = -t111 * t374 + (-t113 * t375 + t115 * t377 + (-t182 * t377 + t186 * t375) * qJD(4)) * t373;
t110 = Icges(5,5) * t204 + Icges(5,6) * t203 - Icges(5,3) * t474;
t112 = Icges(5,4) * t204 + Icges(5,2) * t203 - Icges(5,6) * t474;
t114 = Icges(5,1) * t204 + Icges(5,4) * t203 - Icges(5,5) * t474;
t47 = -t110 * t374 + (-t112 * t375 + t114 * t377 + (-t184 * t377 - t187 * t375) * qJD(4)) * t373;
t60 = t203 * t279 + t204 * t280 + t291 * t297 + t292 * t298 + (-t278 * t522 + t296 * t365) * t373;
t61 = t205 * t279 + t206 * t280 - t289 * t297 + t290 * t298 + (t278 * t515 + t296 * t363) * t373;
t384 = t119 * t566 + t118 * t567 + t106 * t568 + t105 * t569 + t547 / 0.2e1 + t53 * t562 + t545 / 0.2e1 + t541 / 0.2e1 + t542 / 0.2e1 + t543 / 0.2e1 + t540 / 0.2e1 + t538 + t539 + (t109 + ((-t606 + t71 + t74) * t365 + t587 * t363) * t486) * t447 + t29 * t565 + (t54 + t29) * t564 + (t47 + t60) * t444 + (Icges(4,2) * t374 ^ 2 + (Icges(4,1) * t373 + 0.2e1 * Icges(4,4) * t374) * t373 + Icges(3,3)) * t367 + (t46 + t61 + t42) * t446 + (t41 - t530 + ((t587 - t605 - t72) * t365 - t546) * t486) * t445;
t260 = t317 * t370 + t479;
t217 = rSges(5,1) * t291 - rSges(5,2) * t292;
t216 = -rSges(5,1) * t289 - rSges(5,2) * t290;
t33 = -t112 * t289 + t114 * t290 + t184 * t205 + t187 * t206 + (t110 * t363 + t181 * t515) * t373;
t32 = -t113 * t289 + t115 * t290 + t182 * t205 - t186 * t206 + (t111 * t363 + t179 * t515) * t373;
t31 = t112 * t291 + t114 * t292 + t184 * t203 + t187 * t204 + (t110 * t365 - t181 * t522) * t373;
t30 = t113 * t291 + t115 * t292 + t182 * t203 - t186 * t204 + (t111 * t365 - t179 * t522) * t373;
t1 = [Icges(2,3) * qJDD(1) + t384 + (t589 * (t317 + t366) + t590 * (-t315 - t557) + (-t282 - t479 + t260) * t259) * m(3) + ((t337 ^ 2 + t338 ^ 2) * qJDD(1) + g(1) * t337 - g(2) * t338) * m(2) + (t62 * (t396 - t479) + t593 * (t122 + t366) + (t62 + t581) * t63 + t576 + t616 * (t609 - t557)) * m(6) + (t89 * (t420 - t479) + (t268 + t89 + t574 + t580) * t90 + t592 * (t366 + t134) + t577 + t615 * (t610 - t557)) * m(5) + (t168 * (t323 - t440) + t385 + t591 * (t222 + t366) + t594 * (t221 - t557) + (t168 + t580 + t584) * t169) * m(4); t384 + (t581 * t63 + (t396 - t429 + t572) * t62 + t593 * t122 + t576 + t616 * t609) * m(6) + ((t574 + t595) * t90 + (t420 - t429 + t579) * t89 + t592 * t134 + t577 + t615 * t610) * m(5) + (t385 + t591 * t222 + t594 * t221 + (t584 + t608) * t169 + (t323 + t586) * t168) * m(4) + (-t259 * t282 - t260 * t525 + (t259 * t370 + t589) * t317 + (t260 * t370 - t590) * t315) * m(3); (-m(4) - m(5) - m(6)) * (g(1) * t363 - g(2) * t365) + 0.2e1 * (t20 * t560 + t21 * t559) * m(6) + 0.2e1 * (t48 * t560 + t49 * t559) * m(5) + 0.2e1 * (t559 * t88 + t560 * t87) * m(4); ((-t289 * t495 + t290 * t496 + t304 * t521) * t348 + (-t289 * t571 + t290 * t403 + t363 * t411) * t486) * t447 + (-t374 * t60 + ((t370 * t73 + t31) * t365 + (-t370 * t74 + t30) * t363) * t373) * t444 + (-t374 * t61 + ((t370 * t71 + t33) * t365 + (-t370 * t72 + t32) * t363) * t373) * t446 + (t119 * t343 + t261 * t73 + t262 * t74 + t348 * t60 + (t30 * t363 + t31 * t365) * t486) * t463 + t348 * (-t107 * t374 + ((t370 * t82 + t47) * t365 + (-t370 * t83 + t46) * t363) * t373) / 0.2e1 - t348 * (-t374 * t304 * t348 + ((-t375 * t495 + t377 * t496) * t348 + ((-t375 * t571 + t377 * t403) * t373 - t426 * t374) * qJD(4)) * t373) / 0.2e1 + t41 * t448 + t42 * t449 + ((t291 * t495 + t292 * t496 + t304 * t514) * t348 + (t291 * t571 + t403 * t292 + t365 * t411) * t486) * t445 + (-t119 * t374 + t412) * t566 + t343 * (-t135 * t374 + (t363 * t82 + t365 * t83) * t373) / 0.2e1 + (-t118 * t374 + t413) * t567 + (t118 * t343 + t261 * t71 + t262 * t72 + t348 * t61 + (t32 * t363 + t33 * t365) * t486) * t464 + (t541 + t540 + (t363 * t46 + t365 * t47) * t486 + t538) * t558 + t390 + (-g(1) * (t273 + t200) - g(2) * (-t272 + t199) + (-t20 * t188 + t21 * t503) * t374 + t397 - t57 * ((-t272 * t365 - t273 * t363) * t486 + t458) + ((-t103 - t124) * t374 - t273 * t348 + t598) * t63 + (t125 * t374 - t272 * t348 + t596) * t62 + (-g(3) * (t433 - t555) + (t21 * (-t243 - t258) - t544 - t14 * t188 + t57 * t125 + (t62 * t243 + t503 * t57) * t370) * t365 + (t14 * t503 + t57 * (t188 * t370 - t124 + t459) + (t370 * t63 + t20) * t243) * t363) * t373) * m(6) + ((-t116 * t90 + t117 * t89 - t191 * t48 - t192 * t49) * t374 + ((-t191 * t262 - t192 * t261) * t427 + t432 * t299 + ((t370 * t89 - t49) * t365 + (t370 * t90 + t48) * t363) * t283 + (-0.2e1 * t116 * t363 + 0.2e1 * t117 * t365 + t191 * t522 - t192 * t515) * t108) * t373 - (-t216 * t89 + t217 * t90) * t348 - (t108 * (t216 * t365 - t217 * t363) + t432 * t307) * t486 - g(1) * t217 - g(2) * t216 - g(3) * t307) * m(5); t390 + (-t21 * t374 * t172 + ((-t172 * t370 * t57 - t21 * t258 - t544) * t365 + (-t14 * t172 + t459 * t57) * t363) * t373 + t397 - t57 * t458 - g(1) * t200 - g(2) * t199 - g(3) * t293 + (-t103 * t374 + t598) * t63 + t596 * t62) * m(6);];
tau = t1;
