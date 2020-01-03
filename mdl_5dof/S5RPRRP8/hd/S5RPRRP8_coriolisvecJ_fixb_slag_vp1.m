% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP8_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:08
% EndTime: 2019-12-31 18:47:27
% DurationCPUTime: 16.34s
% Computational Cost: add. (9467->517), mult. (21229->618), div. (0->0), fcn. (22089->6), ass. (0->274)
t288 = sin(qJ(4));
t289 = cos(qJ(4));
t462 = Icges(6,5) * t289;
t251 = -Icges(6,1) * t288 + t462;
t468 = Icges(5,4) * t289;
t357 = Icges(5,1) * t288 + t468;
t590 = t251 - t357;
t463 = Icges(6,5) * t288;
t348 = -Icges(6,3) * t289 + t463;
t589 = -t288 * t348 + t590 * t289;
t481 = sin(qJ(3));
t482 = sin(qJ(1));
t483 = cos(qJ(3));
t484 = cos(qJ(1));
t236 = -t482 * t481 - t484 * t483;
t237 = t484 * t481 - t482 * t483;
t356 = Icges(6,1) * t289 + t463;
t129 = Icges(6,4) * t236 + t237 * t356;
t469 = Icges(5,4) * t288;
t358 = Icges(5,1) * t289 - t469;
t133 = Icges(5,5) * t236 + t237 * t358;
t571 = t129 + t133;
t588 = t288 * t571;
t587 = t571 * t289;
t351 = Icges(5,5) * t289 - Icges(5,6) * t288;
t116 = -Icges(5,3) * t237 + t236 * t351;
t353 = Icges(6,4) * t289 + Icges(6,6) * t288;
t120 = -Icges(6,2) * t237 + t236 * t353;
t532 = t116 + t120;
t117 = Icges(5,3) * t236 + t237 * t351;
t121 = Icges(6,2) * t236 + t237 * t353;
t533 = t117 + t121;
t128 = -Icges(6,4) * t237 + t236 * t356;
t132 = -Icges(5,5) * t237 + t236 * t358;
t531 = t128 + t132;
t354 = Icges(5,2) * t289 + t469;
t586 = -t288 * t354 - t589;
t585 = t531 * t288;
t440 = t237 * t289;
t208 = Icges(6,5) * t440;
t441 = t237 * t288;
t457 = Icges(6,6) * t236;
t114 = -Icges(6,3) * t441 - t208 - t457;
t355 = -Icges(5,2) * t288 + t468;
t125 = Icges(5,6) * t236 + t237 * t355;
t584 = t114 + t125;
t444 = t236 * t289;
t209 = Icges(6,5) * t444;
t445 = t236 * t288;
t456 = Icges(6,6) * t237;
t115 = -Icges(6,3) * t445 - t209 + t456;
t124 = -Icges(5,6) * t237 + t236 * t355;
t583 = t115 + t124;
t564 = Icges(5,6) - Icges(6,6);
t566 = Icges(6,4) + Icges(5,5);
t582 = t566 * t288 + t564 * t289;
t573 = Icges(5,4) - Icges(6,5);
t581 = t573 * t288 + (Icges(5,2) + Icges(6,3)) * t289;
t580 = (Icges(5,1) + Icges(6,1)) * t288 + t573 * t289;
t579 = t115 * t441 - t532 * t236 - t531 * t440;
t578 = t114 * t445 + t533 * t237 - t571 * t444;
t350 = Icges(5,5) * t288 + Icges(5,6) * t289;
t149 = t350 * t236;
t352 = Icges(6,4) * t288 - Icges(6,6) * t289;
t151 = t352 * t236;
t577 = t586 * t237 + t149 + t151;
t148 = t350 * t237;
t150 = t352 * t237;
t576 = t586 * t236 - t148 - t150;
t575 = t289 * t531;
t349 = Icges(6,3) * t288 + t462;
t570 = (t349 - t355) * qJD(4);
t569 = (t356 + t358) * qJD(4);
t552 = rSges(6,1) + pkin(4);
t550 = -t124 * t441 - t579;
t549 = -t125 * t445 - t578;
t101 = t120 * t237;
t38 = -t115 * t445 + t128 * t444 - t101;
t475 = -t116 * t237 + t132 * t444;
t40 = -t124 * t445 + t475;
t568 = t38 + t40;
t562 = rSges(6,3) + qJ(5);
t471 = -t117 * t236 - t133 * t440;
t35 = -t125 * t441 - t471;
t399 = t114 * t441 - t121 * t236 - t129 * t440;
t567 = t399 - t35;
t508 = t289 * t584 + t588;
t541 = t289 * t583 + t585;
t565 = Icges(6,2) + Icges(5,3);
t281 = t484 * qJ(2);
t398 = t482 * pkin(1);
t258 = t398 - t281;
t278 = qJD(2) * t482;
t389 = qJD(1) * t484;
t410 = qJ(2) * t389 + t278;
t563 = qJD(1) * t258 - t278 + t410;
t287 = qJD(1) - qJD(3);
t504 = t237 * pkin(3) + t236 * pkin(7);
t181 = t287 * t504;
t190 = t287 * t236;
t404 = qJD(4) * t288;
t393 = t236 * t404;
t191 = t287 * t237;
t448 = t191 * t289;
t336 = t393 + t448;
t403 = qJD(4) * t289;
t337 = -t191 * t288 + t236 * t403;
t560 = (t565 * t190 + t566 * t336 + t564 * t337) * t237;
t559 = t236 * t349 - t124 - t456;
t536 = t237 * t349 - t125 + t457;
t558 = (t351 + t353) * qJD(4);
t557 = t569 * t289 + t570 * t288 + ((t348 - t354) * t289 + t590 * t288) * qJD(4);
t452 = t124 * t288;
t555 = -t115 * t288 - t452 + t575;
t453 = t125 * t288;
t524 = -t114 * t288 - t453 + t587;
t554 = t350 + t352;
t397 = t482 * pkin(2);
t327 = -t398 - t397;
t388 = qJD(1) * t482;
t553 = pkin(2) * t388 + t327 * qJD(1) + t563;
t546 = t577 * t287;
t545 = t576 * t287;
t279 = qJD(2) * t484;
t285 = t484 * pkin(2);
t408 = t484 * pkin(1) + t482 * qJ(2);
t394 = t285 + t408;
t318 = qJD(1) * t394 - t279;
t551 = -(-t124 * t237 - t125 * t236) * t288 + t578 + t579;
t391 = t237 * t404;
t450 = t190 * t289;
t338 = t391 - t450;
t339 = t190 * t288 + t237 * t403;
t548 = t565 * t191 + t566 * t338 + t564 * t339;
t366 = rSges(5,1) * t289 - rSges(5,2) * t288;
t240 = t366 * qJD(4);
t290 = qJD(1) ^ 2;
t422 = (-pkin(1) * t388 + t278 + t410) * qJD(1);
t335 = -t290 * t397 + t422;
t423 = t191 * pkin(3) + t190 * pkin(7);
t321 = t287 * t423 + t335;
t365 = rSges(5,1) * t288 + rSges(5,2) * t289;
t328 = t336 * rSges(5,1) + t190 * rSges(5,3);
t56 = rSges(5,2) * t337 + t328;
t22 = t287 * t56 + (t190 * t365 + t237 * t240) * qJD(4) + t321;
t496 = qJD(1) * t408 - t279;
t379 = (t279 - t496) * qJD(1);
t316 = -t285 * t290 + t379;
t329 = t338 * rSges(5,1) + t191 * rSges(5,3);
t54 = rSges(5,2) * t339 + t329;
t97 = -t190 * pkin(3) + t191 * pkin(7);
t23 = (-t54 - t97) * t287 + (-t191 * t365 + t236 * t240) * qJD(4) + t316;
t319 = t278 + (-t397 - t258) * qJD(1);
t311 = t181 + t319;
t405 = qJD(4) * t365;
t224 = t236 * rSges(5,3);
t420 = rSges(5,1) * t440 + t224;
t141 = rSges(5,2) * t441 - t420;
t540 = t287 * t141;
t59 = t236 * t405 + t311 - t540;
t222 = t237 * rSges(5,3);
t419 = -rSges(5,1) * t444 + t222;
t143 = rSges(5,2) * t445 + t419;
t503 = t236 * pkin(3) - t237 * pkin(7);
t60 = t237 * t405 + (t143 - t503) * t287 + t318;
t544 = (t288 * (-t190 * t59 - t191 * t60 + t22 * t236 - t23 * t237) + (t236 * t60 - t237 * t59) * t403) * rSges(5,2);
t542 = t399 - t38;
t368 = t237 * rSges(4,1) - t236 * rSges(4,2);
t539 = t287 * t368;
t523 = t552 * t288 - t562 * t289;
t380 = qJD(4) * t523;
t402 = qJD(5) * t288;
t205 = t236 * t402;
t225 = t236 * rSges(6,2);
t427 = -t552 * t440 - t441 * t562 - t225;
t520 = -t287 * t427 - t205;
t31 = t236 * t380 + t311 + t520;
t204 = t237 * t402;
t223 = t237 * rSges(6,2);
t426 = -t552 * t444 - t445 * t562 + t223;
t325 = t503 - t426;
t521 = -t287 * t325 - t204;
t32 = t237 * t380 + t318 + t521;
t499 = -t236 * t32 + t237 * t31;
t538 = t499 * t404 * t552;
t527 = (-t549 * t236 + t568 * t237) * qJD(4);
t526 = (t567 * t236 + t550 * t237) * qJD(4);
t411 = t562 * t288 + t552 * t289;
t406 = qJD(4) * t237;
t407 = qJD(4) * t236;
t522 = t508 * t406 - t541 * t407;
t516 = t553 - t181;
t515 = 0.2e1 * qJD(4);
t514 = t526 + t546;
t513 = t527 + t545;
t512 = -t524 * qJD(4) + t191 * t582 + t338 * t580 + t339 * t581;
t511 = t555 * qJD(4) - t190 * t582 - t336 * t580 - t337 * t581;
t510 = t190 * t586 - t191 * t554 + t236 * t558 + t237 * t557;
t509 = -t190 * t554 - t191 * t586 + t236 * t557 - t237 * t558;
t382 = t31 * t523;
t381 = t32 * t523;
t277 = qJD(5) * t289;
t28 = t277 + (t426 * t236 + t427 * t237) * qJD(4);
t506 = qJD(4) * t28;
t502 = t191 * rSges(6,2) - t562 * t339 - t450 * t552 - t204;
t501 = -t190 * rSges(6,2) + t562 * t337 - t448 * t552 + t205;
t324 = t484 * rSges(3,1) + t482 * rSges(3,3);
t500 = t408 + t324;
t431 = t354 * t237 - t133;
t435 = -t357 * t237 - t125;
t495 = t288 * t431 + t289 * t435;
t433 = -t348 * t237 - t129;
t437 = Icges(6,1) * t441 + t114 - t208;
t494 = t288 * t433 - t289 * t437;
t492 = t236 ^ 2;
t487 = -t287 / 0.2e1;
t480 = t552 * t391 + t502;
t479 = t552 * t393 - t501;
t472 = t60 * t365;
t439 = t351 * t287;
t438 = t353 * t287;
t436 = -Icges(6,1) * t445 - t115 + t209;
t434 = -t357 * t236 - t124;
t432 = -t348 * t236 - t128;
t430 = t354 * t236 - t132;
t429 = -t411 * t236 + t223;
t428 = t411 * t237 + t225;
t425 = t523 * t237;
t424 = t523 * t236;
t421 = t411 * qJD(4) - t277;
t416 = t348 + t356;
t415 = -t349 - t251;
t414 = -t354 + t358;
t413 = -t355 - t357;
t401 = -t484 / 0.2e1;
t400 = t482 / 0.2e1;
t396 = t482 * rSges(3,1);
t385 = t407 / 0.2e1;
t384 = -t406 / 0.2e1;
t374 = t503 - t419;
t369 = -t277 + t421;
t96 = t191 * rSges(4,1) - t190 * rSges(4,2);
t95 = -t190 * rSges(4,1) - t191 * rSges(4,2);
t367 = rSges(4,1) * t236 + rSges(4,2) * t237;
t359 = -t236 * t59 - t237 * t60;
t61 = (t141 * t237 + t143 * t236) * qJD(4);
t326 = -t396 - t398;
t317 = t281 + t327;
t315 = t288 * t432 + t289 * t436;
t314 = t288 * t430 + t289 * t434;
t313 = t329 + t97;
t312 = -t328 - t423;
t308 = (-t288 * t415 + t289 * t416) * t287;
t307 = (t288 * t413 + t289 * t414) * t287;
t306 = t317 + t504;
t305 = t97 + t502;
t304 = -t423 + t501;
t291 = t514 * t406 / 0.2e1 + (-t569 * t288 + t570 * t289 + t354 * t404) * t287 + (t509 + t511) * t384 + (t510 - t512 + t513) * t385 + (t589 * t287 - (t541 + t576) * t190 / 0.2e1 - (t508 + t577) * t191 / 0.2e1) * qJD(4);
t283 = t484 * rSges(3,3);
t276 = rSges(3,3) * t389;
t166 = t365 * t236;
t162 = t365 * t237;
t145 = -t290 * t324 + t379;
t144 = qJD(1) * (-rSges(3,1) * t388 + t276) + t422;
t139 = t237 * t366 + t224;
t137 = t236 * t366 - t222;
t111 = -t367 * t287 + t318;
t110 = t319 + t539;
t73 = -t287 * t95 + t316;
t72 = t287 * t96 + t335;
t11 = t191 * t402 + (-t97 - t480) * t287 + (-t191 * t523 + t236 * t369) * qJD(4) + t316;
t10 = -t190 * t402 + t479 * t287 + (t190 * t523 + t237 * t369) * qJD(4) + t321;
t1 = (t427 * t190 - t426 * t191 + t479 * t236 + t480 * t237 - t402) * qJD(4);
t2 = [-t291 + (((t536 * t289 - t588) * t237 + t541 * t236) * qJD(4) + t522) * t487 + (((t550 + t551) * t236 + (t35 + t475 + (-t288 * t536 - t587) * t237 + (-t452 - t533) * t236 - t542) * t237) * qJD(4) + t545) * t385 + ((t492 * t121 + (t40 + (t117 + t452) * t236 - t475) * t236 + ((-t555 + t533) * t237 + ((-t536 - t584) * t288 - t532) * t236 + t549) * t237) * qJD(4) - t546) * t384 + (-(-t237 * t382 + (t381 + t28 * (t427 + t428)) * t236) * qJD(4) + t11 * (t306 - t427) + t10 * (-t325 + t394) - t538 + (-t428 * t287 + t205 - t304 + t516) * t32 + (pkin(2) * t389 - t305 - t318 + t496 + t521) * t31) * m(6) + (-(t472 + t61 * (t139 + t141)) * t407 + t23 * (t306 + t420) + t22 * (-t374 + t394) + t544 + (-t313 - t318) * t59 + (-t287 * t139 - t312 + t516 + t59) * t60) * m(5) + (t73 * (t317 + t368) + t72 * (-t367 + t394) + (-t318 - t95) * t110 + (t110 + t96 - t539 + t553) * t111) * m(4) + (t145 * (t281 + t283 + t326) + t144 * t500 + (t276 + (t396 - t283 + t326) * qJD(1) + t563) * (t500 * qJD(1) - t279)) * m(3); 0.2e1 * (t10 * t401 + t11 * t400) * m(6) + 0.2e1 * (t22 * t401 + t23 * t400) * m(5) + 0.2e1 * (t400 * t73 + t401 * t72) * m(4) + 0.2e1 * (t144 * t401 + t145 * t400) * m(3); t291 + ((-t508 * t237 + (-t289 * t559 + t585) * t236) * qJD(4) + t522) * t487 + (((t101 - t35 - t471 + (t116 - t453) * t237) * t237 + ((t524 + t532) * t236 + ((t559 + t583) * t288 - t533) * t237 - t550) * t236) * qJD(4) - t545) * t385 + (((-t549 - t551) * t237 + (-t40 + t471 + (t288 * t559 + t575) * t236 + (t453 - t532) * t237 + t542) * t236) * qJD(4) + t546) * t384 + (-(-t236 * t381 + (t382 + t28 * (t426 - t429)) * t237) * qJD(4) + t11 * (-t504 + t427) + t10 * t325 + t538 + (t181 + t304 + t520) * t32 + (t204 - (-t503 + t429) * t287 + t305) * t31) * m(6) + (t23 * (-t504 - t420) + t59 * t313 + t22 * t374 - t544 - t59 * (-t137 - t503) * t287 - (-t236 * t472 + (t59 * t365 + t61 * (t137 + t143)) * t237) * qJD(4) + (t181 + t312 - t540) * t60) * m(5) + (-(-t110 * t367 - t111 * t368) * t287 + t110 * t95 - t111 * t96 - t368 * t73 + t72 * t367) * m(4); (((-t413 + t415) * t289 + (t414 + t416) * t288) * t287 + (((-t430 - t432) * t237 + (t431 + t433) * t236) * t289 + ((t434 + t436) * t237 + (-t435 + t437) * t236) * t288) * qJD(4)) * t487 + (t190 * t541 + t508 * t191 + t512 * t236 + t511 * t237) * t287 / 0.2e1 + ((t150 * t407 + t438) * t236 + (t308 + (t315 * t237 + (-t151 - t494) * t236) * qJD(4)) * t237 + (t148 * t407 + t439) * t236 + (t307 + (t314 * t237 + (-t149 - t495) * t236) * qJD(4)) * t237) * t385 + ((t151 * t406 - t438) * t237 + (t308 + (-t494 * t236 + (-t150 + t315) * t237) * qJD(4)) * t236 + (t149 * t406 - t439) * t237 + (t307 + (-t495 * t236 + (-t148 + t314) * t237) * qJD(4)) * t236) * t384 + ((-t28 * t426 - t382) * t191 - (-t28 * t427 - t381) * t190 + (t1 * t427 + t10 * t523 + t28 * t480 + t32 * t421) * t237 + (t1 * t426 + t11 * t523 + t28 * t479 + t31 * t421) * t236 - (-t28 * t288 + (-t236 * t31 - t237 * t32) * t289) * qJD(5) - (-t31 * t425 + t32 * t424) * t287 - ((t28 * t425 + t32 * t411) * t237 + (t28 * t424 + t31 * t411) * t236) * qJD(4)) * m(6) + (-(-t162 * t59 + t166 * t60) * t287 - (t61 * (t162 * t237 + t166 * t236) - t359 * t366) * qJD(4) + 0.2e1 * t61 * (t141 * t190 - t143 * t191 + t236 * t56 + t237 * t54) - t359 * t240 - (-t190 * t60 + t191 * t59 - t22 * t237 - t23 * t236) * t365) * m(5) - (t510 * t287 + (-t567 * t191 + t550 * t190 + (t190 * t555 - t532 * t191) * t237 + (-t524 * t190 + t533 * t191 + t548 * t236 - t560) * t236) * t515) * t236 / 0.2e1 + (t509 * t287 + (t549 * t191 + t568 * t190 + (-t532 * t190 - t191 * t555 + t560) * t237 + (t533 * t190 + t524 * t191 - t548 * t237) * t236) * t515) * t237 / 0.2e1 + (t513 + t527) * t190 / 0.2e1 + (t514 + t526) * t191 / 0.2e1; t1 * t289 * m(6) + 0.2e1 * (m(6) * (-t10 * t237 - t11 * t236 - t190 * t32 + t191 * t31 - t506) / 0.2e1 - m(6) * (t499 * t287 + (-t237 ^ 2 - t492) * t506) / 0.2e1) * t288;];
tauc = t2(:);
