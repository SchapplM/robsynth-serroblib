% Calculate vector of inverse dynamics joint torques for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:09
% EndTime: 2019-12-05 18:58:24
% DurationCPUTime: 9.55s
% Computational Cost: add. (23835->640), mult. (14122->781), div. (0->0), fcn. (10979->10), ass. (0->381)
t315 = qJ(1) + qJ(2);
t307 = qJ(3) + t315;
t295 = sin(t307);
t312 = qJD(4) + qJD(5);
t217 = t295 * t312;
t296 = cos(t307);
t287 = qJDD(4) * t296;
t313 = qJD(1) + qJD(2);
t301 = qJD(3) + t313;
t133 = qJDD(5) * t296 - t217 * t301 + t287;
t439 = qJD(4) * t301;
t190 = -t295 * t439 + t287;
t314 = qJ(4) + qJ(5);
t304 = cos(t314);
t294 = t304 * rSges(6,1);
t302 = sin(t314);
t239 = -rSges(6,2) * t302 + t294;
t207 = t239 * t312;
t218 = t296 * t312;
t237 = rSges(6,1) * t302 + rSges(6,2) * t304;
t311 = qJDD(1) + qJDD(2);
t300 = qJDD(3) + t311;
t316 = sin(qJ(4));
t523 = pkin(8) * t295;
t524 = pkin(3) * t296;
t224 = t523 + t524;
t303 = sin(t315);
t305 = cos(t315);
t310 = t313 ^ 2;
t317 = sin(qJ(1));
t319 = cos(qJ(1));
t322 = qJD(1) ^ 2;
t359 = (-qJDD(1) * t317 - t319 * t322) * pkin(1);
t329 = t359 + (-t303 * t311 - t305 * t310) * pkin(2);
t522 = t296 * pkin(8);
t389 = -pkin(3) * t295 + t522;
t327 = -t301 ^ 2 * t224 + t300 * t389 + t329;
t318 = cos(qJ(4));
t465 = t318 * qJD(4) ^ 2;
t471 = t296 * t304;
t429 = rSges(6,1) * t471;
t367 = rSges(6,3) * t295 + t429;
t477 = t295 * t304;
t478 = t295 * t302;
t183 = rSges(6,1) * t478 + rSges(6,2) * t477;
t472 = t296 * t302;
t260 = rSges(6,2) * t472;
t422 = t183 * t312 + t301 * t260;
t91 = -t301 * t367 + t422;
t438 = qJD(4) * t316;
t419 = t295 * t438;
t265 = pkin(4) * t419;
t320 = -pkin(9) - pkin(8);
t468 = t301 * t320;
t448 = t295 * t468 + t265;
t309 = t318 * pkin(4);
t297 = t309 + pkin(3);
t521 = pkin(3) - t297;
t99 = (t296 * t521 + t523) * t301 + t448;
t519 = t99 + t91;
t565 = rSges(6,2) * t478 + t296 * rSges(6,3);
t520 = pkin(8) + t320;
t567 = t520 * t296;
t369 = t567 - t565;
t405 = -t297 - t294;
t391 = -pkin(3) - t405;
t586 = t295 * t391 + t369;
t578 = -t133 * t237 - t218 * t207 + t519 * t301 - t586 * t300 + (-t190 * t316 - t296 * t465) * pkin(4) + t327 - g(3);
t256 = t296 * t297;
t148 = t295 * t520 - t256 + t524;
t161 = -t260 + t367;
t460 = t148 - t161;
t387 = rSges(3,1) * t303 + rSges(3,2) * t305;
t208 = t387 * t313;
t513 = pkin(1) * qJD(1);
t432 = t317 * t513;
t185 = t208 + t432;
t203 = t301 * t224;
t594 = t460 * t301 - t203 - t422 - t448;
t254 = Icges(6,4) * t478;
t159 = -Icges(6,1) * t477 + Icges(6,5) * t296 + t254;
t255 = Icges(6,4) * t472;
t160 = Icges(6,1) * t471 + Icges(6,5) * t295 - t255;
t292 = Icges(6,4) * t304;
t379 = -Icges(6,2) * t302 + t292;
t564 = Icges(6,1) * t302 + t292;
t584 = t564 + t379;
t332 = t217 * (-Icges(6,2) * t471 + t160 - t255) + t218 * (Icges(6,2) * t477 + t159 + t254) + t301 * t584;
t157 = Icges(6,6) * t296 - t295 * t379;
t158 = Icges(6,4) * t471 - Icges(6,2) * t472 + Icges(6,6) * t295;
t504 = Icges(6,4) * t302;
t233 = Icges(6,2) * t304 + t504;
t236 = Icges(6,1) * t304 - t504;
t543 = t217 * (t296 * t564 + t158) + t218 * (-t295 * t564 + t157) + t301 * (t233 - t236);
t593 = t332 * t302 + t304 * t543;
t437 = qJD(4) * t318;
t416 = t296 * t437;
t417 = t296 * t438;
t475 = t295 * t318;
t434 = rSges(5,1) * t475;
t421 = -rSges(5,1) * t417 - rSges(5,2) * t416 - t301 * t434;
t476 = t295 * t316;
t444 = rSges(5,2) * t476 + t296 * rSges(5,3);
t107 = t301 * t444 + t421;
t479 = t295 * t301;
t251 = pkin(3) * t479;
t473 = t296 * t301;
t176 = pkin(8) * t473 - t251;
t189 = qJDD(4) * t295 + t296 * t439;
t516 = rSges(5,2) * t316;
t518 = rSges(5,1) * t318;
t284 = -t516 + t518;
t248 = t284 * qJD(4);
t282 = rSges(5,1) * t316 + rSges(5,2) * t318;
t527 = pkin(1) * t319;
t528 = pkin(1) * t317;
t390 = -qJDD(1) * t527 + t322 * t528;
t525 = pkin(2) * t305;
t526 = pkin(2) * t303;
t340 = t310 * t526 - t311 * t525 + t390;
t441 = qJD(4) * t295;
t470 = t296 * t316;
t272 = rSges(5,2) * t470;
t469 = t296 * t318;
t368 = -rSges(5,1) * t469 - rSges(5,3) * t295;
t171 = -t272 - t368;
t451 = -t171 - t224;
t45 = t248 * t441 + t189 * t282 + (-t107 - t176) * t301 + t451 * t300 + t340;
t580 = -g(2) + t45;
t420 = rSges(5,1) * t419 + (t295 * t437 + t301 * t470) * rSges(5,2);
t108 = t301 * t368 + t420;
t372 = t434 - t444;
t440 = qJD(4) * t296;
t44 = t301 * t108 - t190 * t282 - t248 * t440 - t300 * t372 + t327;
t579 = -g(3) + t44;
t132 = qJD(5) * t473 + qJDD(5) * t295 + t189;
t425 = -t224 + t460;
t433 = rSges(6,1) * t477;
t589 = t218 * t237;
t423 = -t301 * t433 - t589;
t90 = t301 * t565 + t423;
t430 = pkin(4) * t438;
t450 = -t296 * t468 - t297 * t479;
t98 = t251 + (-pkin(8) * t301 - t430) * t296 + t450;
t18 = t132 * t237 + t207 * t217 + (t189 * t316 + t295 * t465) * pkin(4) + (-t176 - t90 - t98) * t301 + t425 * t300 + t340;
t577 = t18 - g(2);
t175 = -rSges(4,1) * t473 + rSges(4,2) * t479;
t386 = -rSges(4,1) * t295 - rSges(4,2) * t296;
t576 = t301 * t175 + t300 * t386 - g(3) + t329;
t174 = rSges(4,1) * t479 + rSges(4,2) * t473;
t219 = rSges(4,1) * t296 - t295 * rSges(4,2);
t575 = t174 * t301 - t219 * t300 - g(2) + t340;
t466 = t305 * t313;
t467 = t303 * t313;
t209 = -rSges(3,1) * t466 + rSges(3,2) * t467;
t592 = t209 * t313 - t311 * t387 - g(3) + t359;
t240 = rSges(3,1) * t305 - t303 * rSges(3,2);
t591 = t208 * t313 - t240 * t311 - g(2) + t390;
t483 = t237 * t296;
t590 = t183 * t217 + t218 * t483 + t296 * t90;
t306 = Icges(5,4) * t318;
t380 = -Icges(5,2) * t316 + t306;
t563 = Icges(5,1) * t316 + t306;
t442 = t563 + t380;
t505 = Icges(5,4) * t316;
t276 = Icges(5,2) * t318 + t505;
t279 = Icges(5,1) * t318 - t505;
t443 = t276 - t279;
t588 = (t316 * t442 + t318 * t443) * t301;
t435 = pkin(2) * t466;
t587 = t435 + t594;
t585 = t301 * t183 - t218 * t239 - t237 * t479;
t153 = t301 * t171;
t583 = -t153 - t203 - t420;
t582 = -t295 * t521 + t567;
t394 = t217 * t237 + t265;
t356 = t394 - t435;
t431 = t319 * t513;
t335 = t356 - t431;
t59 = t301 * t425 + t335;
t581 = (t207 * t295 - t217 * t239 + t237 * t473 - t301 * t483) * t59;
t156 = Icges(6,5) * t471 - Icges(6,6) * t472 + Icges(6,3) * t295;
t63 = t296 * t156 + t158 * t478 - t160 * t477;
t375 = t233 * t302 - t304 * t564;
t231 = Icges(6,5) * t302 + Icges(6,6) * t304;
t487 = t231 * t296;
t96 = t295 * t375 + t487;
t574 = t217 * t63 + t96 * t301;
t168 = Icges(5,4) * t469 - Icges(5,2) * t470 + Icges(5,6) * t295;
t270 = Icges(5,4) * t470;
t170 = Icges(5,1) * t469 + Icges(5,5) * t295 - t270;
t376 = t168 * t316 - t170 * t318;
t571 = t296 * t376;
t378 = t158 * t302 - t160 * t304;
t568 = t378 * t296;
t199 = t301 * t219;
t357 = t431 + t435;
t147 = -t357 - t199;
t485 = t233 * t312;
t561 = -Icges(6,6) * t301 + t485;
t152 = t301 * t372;
t202 = t301 * t389;
t370 = t282 * t440 + t152 - t202;
t436 = pkin(2) * t467;
t339 = t370 + t436;
t77 = t339 + t432;
t213 = t282 * t441;
t348 = t213 - t357;
t78 = t301 * t451 + t348;
t560 = t295 * t78 + t296 * t77;
t558 = -Icges(6,3) * t301 + t231 * t312;
t557 = -Icges(6,5) * t301 + t312 * t564;
t363 = t380 * t301;
t553 = -Icges(5,6) * t301 + qJD(4) * t276;
t103 = t295 * t553 - t296 * t363;
t365 = t279 * t301;
t551 = -Icges(5,5) * t301 + qJD(4) * t563;
t105 = t295 * t551 - t296 * t365;
t167 = Icges(5,6) * t296 - t295 * t380;
t269 = Icges(5,4) * t476;
t169 = -Icges(5,1) * t475 + Icges(5,5) * t296 + t269;
t109 = t167 * t318 + t169 * t316;
t275 = Icges(5,5) * t318 - Icges(5,6) * t316;
t165 = Icges(5,3) * t296 - t275 * t295;
t556 = qJD(4) * t109 + t103 * t316 - t105 * t318 - t165 * t301;
t102 = -t295 * t363 - t296 * t553;
t104 = -t295 * t365 - t296 * t551;
t110 = t168 * t318 + t170 * t316;
t166 = Icges(5,5) * t469 - Icges(5,6) * t470 + Icges(5,3) * t295;
t555 = qJD(4) * t110 + t102 * t316 - t104 * t318 - t166 * t301;
t274 = Icges(5,5) * t316 + Icges(5,6) * t318;
t554 = -Icges(5,3) * t301 + qJD(4) * t274;
t243 = t380 * qJD(4);
t244 = t279 * qJD(4);
t374 = t276 * t318 + t316 * t563;
t552 = qJD(4) * t374 + t243 * t316 - t244 * t318 - t274 * t301;
t413 = -pkin(3) - t518;
t529 = -rSges(5,3) - pkin(8);
t325 = ((-t516 * t78 - t529 * t77) * t295 + (-t413 * t77 + t529 * t78) * t296) * t301;
t550 = t325 + (t213 + t583) * t77;
t452 = -Icges(5,2) * t469 + t170 - t270;
t454 = t296 * t563 + t168;
t548 = t316 * t452 + t318 * t454;
t453 = Icges(5,2) * t475 + t169 + t269;
t455 = -t295 * t563 + t167;
t547 = -t316 * t453 - t318 * t455;
t395 = -t236 * t312 + t485;
t396 = t584 * t312;
t546 = -t231 * t301 + t302 * t396 + t304 * t395;
t362 = t379 * t301;
t401 = t160 * t312 - t295 * t362 - t296 * t561;
t364 = t236 * t301;
t403 = t158 * t312 + t295 * t364 + t296 * t557;
t545 = -t156 * t301 + t302 * t401 + t304 * t403;
t232 = Icges(6,5) * t304 - Icges(6,6) * t302;
t360 = t232 * t295;
t155 = Icges(6,3) * t296 - t360;
t402 = t159 * t312 + t295 * t561 - t296 * t362;
t404 = t157 * t312 - t295 * t557 + t296 * t364;
t544 = -t155 * t301 + t302 * t402 + t304 * t404;
t542 = t132 / 0.2e1;
t541 = t133 / 0.2e1;
t540 = t189 / 0.2e1;
t539 = t190 / 0.2e1;
t538 = -t217 / 0.2e1;
t537 = t217 / 0.2e1;
t536 = -t218 / 0.2e1;
t535 = t218 / 0.2e1;
t534 = t295 / 0.2e1;
t533 = t296 / 0.2e1;
t532 = t300 / 0.2e1;
t531 = -t301 / 0.2e1;
t530 = t301 / 0.2e1;
t509 = t296 * t98;
t177 = t231 * t295;
t97 = -t296 * t375 + t177;
t508 = t97 * t301;
t192 = t274 * t295;
t373 = t276 * t316 - t318 * t563;
t117 = -t296 * t373 + t192;
t495 = t117 * t301;
t492 = t157 * t302;
t491 = t159 * t304;
t490 = t167 * t316;
t489 = t169 * t318;
t484 = t237 * t295;
t482 = t274 * t296;
t481 = t275 * t301;
t200 = t282 * t295;
t480 = t282 * t296;
t474 = t296 * t148;
t464 = t296 * t155 + t157 * t478;
t463 = -t295 * t155 - t159 * t471;
t462 = t296 * t165 + t167 * t476;
t461 = t295 * t165 + t169 * t469;
t273 = pkin(4) * t476;
t415 = -t479 / 0.2e1;
t414 = t473 / 0.2e1;
t410 = -t441 / 0.2e1;
t409 = t441 / 0.2e1;
t408 = -t440 / 0.2e1;
t407 = t440 / 0.2e1;
t397 = -t166 - t489;
t393 = pkin(4) * t417;
t392 = t251 - t421;
t285 = rSges(2,1) * t319 - t317 * rSges(2,2);
t388 = rSges(2,1) * t317 + rSges(2,2) * t319;
t70 = -t169 * t475 + t462;
t71 = t296 * t166 + t168 * t476 - t170 * t475;
t385 = t295 * t71 + t296 * t70;
t72 = -t167 * t470 + t461;
t73 = t166 * t295 - t571;
t384 = t295 * t73 + t296 * t72;
t93 = t158 * t304 + t160 * t302;
t377 = -t489 + t490;
t371 = t433 - t565;
t64 = -t157 * t472 - t463;
t188 = -t219 - t525;
t186 = -t240 * t313 - t431;
t358 = t432 + t436;
t187 = t386 - t526;
t355 = t392 + t436;
t351 = t177 * t218 - t217 * t487 + t232 * t301;
t349 = t175 - t435;
t347 = -t296 * t558 + (-t360 + t378) * t301;
t346 = -t232 * t473 + t558 * t295 + (-t491 + t492) * t301;
t345 = -t481 * t295 - t296 * t554 + t301 * t376;
t344 = t295 * t554 - t481 * t296 + t301 * t377;
t343 = t232 * t312 + t301 * t375;
t342 = t275 * qJD(4) + t301 * t373;
t123 = -t429 - t256 + t260 + (-rSges(6,3) + t320) * t295;
t338 = t393 - t423 - t450;
t130 = t295 * t413 + t444 + t522;
t131 = t295 * t529 + t296 * t413 + t272;
t337 = t589 - t202 + t393 + (t582 + t371) * t301;
t336 = t296 * t171 + t295 * t372;
t122 = t295 * t405 - t296 * t320 + t565;
t13 = t346 * t295 - t296 * t544;
t14 = t347 * t295 - t296 * t545;
t15 = t295 * t544 + t346 * t296;
t16 = t295 * t545 + t347 * t296;
t62 = -t159 * t477 + t464;
t30 = t218 * t62 + t574;
t65 = t156 * t295 - t568;
t31 = t217 * t65 + t218 * t64 + t508;
t40 = -t302 * t404 + t304 * t402;
t41 = -t302 * t403 + t304 * t401;
t46 = t343 * t295 - t296 * t546;
t47 = t295 * t546 + t343 * t296;
t92 = t157 * t304 + t159 * t302;
t334 = (t13 * t218 + t132 * t65 + t133 * t64 + t14 * t217 + t300 * t97 + t301 * t46) * t534 + (t351 * t295 - t593 * t296) * t538 + (t593 * t295 + t351 * t296) * t536 + (t132 * t63 + t133 * t62 + t15 * t218 + t16 * t217 + t300 * t96 + t301 * t47) * t533 + (-t302 * t543 + t304 * t332) * t531 + t30 * t415 + t31 * t414 + ((t301 * t65 + t13) * t296 + (-t301 * t64 + t14) * t295) * t537 + (t295 * t65 + t296 * t64) * t542 + (t295 * t63 + t296 * t62) * t541 + ((t301 * t63 + t15) * t296 + (-t301 * t62 + t16) * t295) * t535 + (t295 * t93 + t296 * t92) * t532 + ((t301 * t93 + t40) * t296 + (-t301 * t92 + t41) * t295) * t530;
t115 = t123 - t525;
t124 = t130 - t526;
t125 = t131 - t525;
t114 = t122 - t526;
t331 = t338 + t436;
t330 = t337 + t436;
t328 = t568 + (-t156 - t491) * t295 + t464;
t58 = t330 + t432;
t326 = (-t59 * t565 - t58 * (-t367 - t256)) * t301;
t116 = t295 * t373 + t482;
t113 = t116 * t301;
t36 = qJD(4) * t385 + t113;
t37 = qJD(4) * t384 + t495;
t50 = -qJD(4) * t377 + t103 * t318 + t105 * t316;
t51 = -qJD(4) * t376 + t102 * t318 + t104 * t316;
t56 = t342 * t295 - t296 * t552;
t57 = t295 * t552 + t342 * t296;
t324 = ((t65 + t328) * t218 + t574) * t538 + (t113 + ((t462 + t73 + t571) * t296 + (-t72 + (t397 - t490) * t296 + t71 + t461) * t295) * qJD(4)) * t410 + (t93 + t97) * t542 + (t92 + t96) * t541 + (t117 + t110) * t540 + (t116 + t109) * t539 + (-t508 + (t63 + (-t156 + t492) * t296 - t378 * t295 + t463) * t218 + (-t62 + t328) * t217 + t31) * t536 + (t40 + t47) * t535 + (t37 - t495 + ((t71 + (-t166 + t490) * t296 - t461) * t296 + (t295 * t397 + t462 - t70) * t295) * qJD(4)) * t408 + (t50 + t57) * t407 + (-qJD(4) * t373 + t243 * t318 + t244 * t316 - t302 * t395 + t304 * t396) * t301 + (t41 + t46 + t30) * t537 + (t51 + t56 + t36) * t409 + (t233 * t304 + t302 * t564 + Icges(4,3) + t374) * t300;
t323 = Icges(3,3) * t311 + t324;
t198 = t301 * t386;
t146 = -t198 + t358;
t140 = t296 * t161;
t95 = t336 * qJD(4);
t52 = t218 * t161 + t217 * t371 + (t295 * t582 - t474) * qJD(4);
t22 = t295 * t555 + t345 * t296;
t21 = t295 * t556 + t344 * t296;
t20 = t345 * t295 - t296 * t555;
t19 = t344 * t295 - t296 * t556;
t10 = -t190 * t148 + t189 * t582 + t133 * t161 + t218 * t90 + t132 * t371 - t217 * t91 + (-t295 * t99 + t509) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t323 + (t147 * (t358 + t174) - t146 * (t349 - t431) + t575 * (t188 - t527) + t576 * (t187 - t528)) * m(4) + (t591 * (-t240 - t527) + t592 * (-t387 - t528) + (t186 - t209 + t431) * t185) * m(3) + ((qJDD(1) * t388 + g(3)) * t388 + (qJDD(1) * t285 + g(2)) * t285) * m(2) + (t59 * (t331 + t432) + t326 + (t431 + t335 - t59 + t587) * t58 + t577 * (t115 - t527) + t578 * (t114 - t528)) * m(6) + (t78 * (t355 + t432) + t325 + (t357 + t348 - t78 + t583) * t77 + t580 * (t125 - t527) + t579 * (t124 - t528)) * m(5); t323 + (t326 + (t331 - t330) * t59 + (t356 + t587) * t58 + t577 * t115 + t578 * t114) * m(6) + ((t355 - t339) * t78 + t580 * t125 + t579 * t124 + t550) * m(5) + (t575 * t188 + t576 * t187 + (-t199 - t435 - t349) * t146 + (t198 + t174) * t147) * m(4) + (-t185 * t209 + t186 * t208 + (-t185 * t313 - t591) * t240 - (t186 * t313 + t592) * t387) * m(3); t324 + (t326 + (t338 - t337) * t59 + (t394 + t594) * t58 + t577 * t123 + t578 * t122) * m(6) + ((-t370 + t392) * t78 + t580 * t131 + t579 * t130 + t550) * m(5) + (-t146 * t175 + t147 * t174 + (t147 * t301 + t576) * t386 + (-t146 * t301 - t575) * t219) * m(4); (t117 * t300 + t189 * t73 + t190 * t72 + t301 * t56 + (t19 * t296 + t20 * t295) * qJD(4)) * t534 + (t116 * t300 + t189 * t71 + t190 * t70 + t301 * t57 + (t21 * t296 + t22 * t295) * qJD(4)) * t533 + ((t192 * t440 + t481) * t296 + (t588 + (t548 * t295 + (-t482 - t547) * t296) * qJD(4)) * t295) * t408 + ((t110 * t301 + t50) * t296 + (-t109 * t301 + t51) * t295) * t530 + ((-t316 * t443 + t318 * t442) * t301 + ((t295 * t452 + t296 * t453) * t318 + (-t295 * t454 - t296 * t455) * t316) * qJD(4)) * t531 + ((-t441 * t482 + t481) * t295 + (-t588 + (t547 * t296 + (t192 - t548) * t295) * qJD(4)) * t296) * t410 + ((t301 * t73 + t19) * t296 + (-t301 * t72 + t20) * t295) * t409 + ((t301 * t71 + t21) * t296 + (-t301 * t70 + t22) * t295) * t407 + t334 + t36 * t415 + t37 * t414 + t384 * t540 + t385 * t539 + (t109 * t296 + t110 * t295) * t532 + (t10 * (t295 * t586 + t140 - t474) + t18 * (t273 + t484) - g(1) * (t239 + t309) - g(2) * (t273 + t183) + t578 * t296 * (-pkin(4) * t316 - t237) + (-(-t295 ^ 2 - t296 ^ 2) * t430 + t509 - t519 * t295 + (t369 * t296 + (t296 * t391 + t460) * t295) * t301 + t590) * t52 + t581 + (-pkin(4) * t416 - (-pkin(4) * t437 - t207) * t296 + t585) * t58) * m(6) + (-(-t200 * t77 + t480 * t78) * t301 - (t95 * (-t200 * t295 - t296 * t480) + t560 * t284) * qJD(4) + (t190 * t171 + t189 * t372 + (t296 * t107 - t295 * t108) * qJD(4)) * t336 + t95 * ((-t108 - t153) * t295 + (t107 + t152) * t296) + t45 * t200 - t44 * t480 + (t78 * t473 - t77 * t479) * t282 + t560 * t248 - g(1) * t284 - g(2) * t200 + g(3) * t480) * m(5); t334 + (t10 * (t295 * t371 + t140) + t18 * t484 - g(1) * t239 - g(2) * t183 + (-t295 * t91 + (-t295 * t161 + t296 * t371) * t301 + t590) * t52 + t581 + (t207 * t296 + t585) * t58 - t578 * t483) * m(6);];
tau = t1;
