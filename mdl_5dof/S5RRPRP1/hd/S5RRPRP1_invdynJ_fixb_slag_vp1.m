% Calculate vector of inverse dynamics joint torques for
% S5RRPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:00
% EndTime: 2019-12-05 18:22:21
% DurationCPUTime: 15.06s
% Computational Cost: add. (13377->512), mult. (9752->591), div. (0->0), fcn. (7457->8), ass. (0->303)
t591 = Icges(5,3) + Icges(6,3);
t283 = qJ(1) + qJ(2);
t272 = pkin(8) + t283;
t266 = sin(t272);
t267 = cos(t272);
t287 = cos(qJ(4));
t423 = t267 * t287;
t285 = sin(qJ(4));
t424 = t267 * t285;
t140 = Icges(6,4) * t423 - Icges(6,2) * t424 + Icges(6,6) * t266;
t142 = Icges(5,4) * t423 - Icges(5,2) * t424 + Icges(5,6) * t266;
t583 = t140 + t142;
t230 = Icges(6,4) * t424;
t144 = Icges(6,1) * t423 + Icges(6,5) * t266 - t230;
t232 = Icges(5,4) * t424;
t146 = Icges(5,1) * t423 + Icges(5,5) * t266 - t232;
t581 = t144 + t146;
t241 = Icges(6,5) * t287 - Icges(6,6) * t285;
t243 = Icges(5,5) * t287 - Icges(5,6) * t285;
t594 = t241 + t243;
t593 = Icges(5,5) + Icges(6,5);
t592 = -Icges(5,6) - Icges(6,6);
t448 = Icges(6,4) * t285;
t244 = Icges(6,2) * t287 + t448;
t449 = Icges(5,4) * t285;
t246 = Icges(5,2) * t287 + t449;
t587 = t244 + t246;
t249 = Icges(6,1) * t287 - t448;
t251 = Icges(5,1) * t287 - t449;
t590 = t249 + t251;
t275 = Icges(6,4) * t287;
t341 = -Icges(6,2) * t285 + t275;
t276 = Icges(5,4) * t287;
t342 = -Icges(5,2) * t285 + t276;
t589 = t341 + t342;
t499 = Icges(5,1) * t285 + t276;
t500 = Icges(6,1) * t285 + t275;
t586 = t499 + t500;
t588 = -t594 * t266 + t591 * t267;
t139 = Icges(6,6) * t267 - t266 * t341;
t141 = Icges(5,6) * t267 - t266 * t342;
t584 = t139 + t141;
t427 = t266 * t285;
t229 = Icges(6,4) * t427;
t426 = t266 * t287;
t143 = -Icges(6,1) * t426 + Icges(6,5) * t267 + t229;
t231 = Icges(5,4) * t427;
t145 = -Icges(5,1) * t426 + Icges(5,5) * t267 + t231;
t582 = t143 + t145;
t565 = t583 * t285 - t581 * t287;
t585 = t591 * t266 + t593 * t423 + t592 * t424;
t240 = Icges(6,5) * t285 + Icges(6,6) * t287;
t242 = Icges(5,5) * t285 + Icges(5,6) * t287;
t580 = t240 + t242;
t282 = qJD(1) + qJD(2);
t579 = t589 * t282;
t578 = t590 * t282;
t577 = t586 * qJD(4) - t593 * t282;
t576 = t587 * qJD(4) + t592 * t282;
t567 = t587 * t285 - t586 * t287;
t514 = -t588 * t266 - t582 * t423;
t513 = t588 * t267 + t584 * t427;
t575 = t582 * t287;
t512 = t584 * t285;
t574 = t565 * t267;
t531 = -t582 * t426 + t513;
t530 = t585 * t267 - t581 * t426 + t583 * t427;
t529 = -t584 * t424 - t514;
t528 = t585 * t266 - t574;
t527 = t582 * t285 + t584 * t287;
t526 = t581 * t285 + t583 * t287;
t573 = -t579 * t266 - t576 * t267;
t572 = t576 * t266 - t579 * t267;
t571 = -t578 * t266 - t577 * t267;
t570 = t577 * t266 - t578 * t267;
t569 = t589 * qJD(4);
t568 = t590 * qJD(4);
t566 = t586 * t285 + t587 * t287;
t564 = t512 - t575;
t430 = t243 * t282;
t432 = t241 * t282;
t563 = -t430 - t432;
t562 = t580 * qJD(4) - t591 * t282;
t431 = t242 * t267;
t433 = t240 * t267;
t522 = t567 * t266 + t431 + t433;
t157 = t240 * t266;
t159 = t242 * t266;
t521 = -t567 * t267 + t157 + t159;
t238 = rSges(5,2) * t424;
t456 = rSges(5,3) * t266;
t329 = -rSges(5,1) * t423 - t456;
t148 = -t238 - t329;
t132 = t282 * t148;
t255 = rSges(5,1) * t285 + rSges(5,2) * t287;
t396 = qJD(4) * t266;
t183 = t255 * t396;
t465 = pkin(7) * t266;
t466 = pkin(3) * t267;
t188 = t465 + t466;
t274 = cos(t283);
t421 = t274 * t282;
t390 = pkin(2) * t421;
t356 = t282 * t188 + t390;
t393 = qJD(4) * t285;
t377 = t266 * t393;
t392 = qJD(4) * t287;
t516 = t266 * t392 + t282 * t424;
t380 = rSges(5,1) * t377 + t516 * rSges(5,2);
t561 = -t356 + t183 - t132 - t380;
t560 = t594 * qJD(4) + t567 * t282;
t559 = t562 * t266 + t563 * t267 + t564 * t282;
t558 = t563 * t266 - t562 * t267 + t565 * t282;
t557 = t528 * t266 + t529 * t267;
t556 = t530 * t266 + t531 * t267;
t262 = t266 * rSges(4,2);
t459 = rSges(4,1) * t267;
t187 = -t262 + t459;
t428 = t266 * t282;
t221 = rSges(4,2) * t428;
t467 = pkin(2) * t274;
t352 = -t459 - t467;
t555 = -t390 - t221 + (-t187 - t352) * t282;
t554 = t566 * qJD(4) - t580 * t282 + t569 * t285 - t568 * t287;
t553 = t527 * qJD(4) - t588 * t282 + t572 * t285 - t570 * t287;
t552 = t526 * qJD(4) - t282 * t585 + t573 * t285 - t571 * t287;
t551 = t522 * t282;
t273 = sin(t283);
t350 = rSges(3,1) * t273 + rSges(3,2) * t274;
t178 = t350 * t282;
t286 = sin(qJ(1));
t455 = pkin(1) * qJD(1);
t388 = t286 * t455;
t152 = t178 + t388;
t284 = -qJ(5) - pkin(7);
t462 = pkin(7) + t284;
t503 = rSges(6,2) * t427 + t267 * rSges(6,3);
t549 = t462 * t267 - t503;
t548 = t521 * t282;
t547 = t559 * t266 - t553 * t267;
t546 = t558 * t266 - t552 * t267;
t224 = pkin(3) * t428;
t425 = t267 * t282;
t154 = pkin(7) * t425 - t224;
t394 = qJD(4) * t282;
t173 = qJDD(4) * t266 + t267 * t394;
t278 = t287 * rSges(6,1);
t257 = -rSges(6,2) * t285 + t278;
t207 = t257 * qJD(4);
t452 = t287 * rSges(6,2);
t254 = rSges(6,1) * t285 + t452;
t260 = qJD(5) * t266;
t281 = qJDD(1) + qJDD(2);
t280 = t282 ^ 2;
t290 = qJD(1) ^ 2;
t288 = cos(qJ(1));
t469 = pkin(1) * t288;
t470 = pkin(1) * t286;
t354 = -qJDD(1) * t469 + t290 * t470;
t468 = pkin(2) * t273;
t330 = t280 * t468 + t354;
t364 = -t188 - t467;
t279 = t287 * pkin(4);
t268 = t279 + pkin(3);
t206 = t267 * t268;
t237 = rSges(6,2) * t424;
t389 = rSges(6,1) * t423;
t328 = rSges(6,3) * t266 + t389;
t414 = t266 * t462 - t206 + t237 - t328 + t466;
t332 = t364 + t414;
t419 = t287 * qJD(4) ^ 2;
t374 = t267 * t392;
t420 = t282 * t284;
t375 = t267 * t393;
t517 = -t282 * t426 - t375;
t498 = t517 * rSges(6,1) - rSges(6,2) * t374 - t267 * t420 - t268 * t428 + t260;
t461 = t224 + (-pkin(4) * t393 - pkin(7) * t282) * t267 + t282 * t503 + t498;
t15 = t207 * t396 + qJDD(5) * t267 + t173 * t254 + (t173 * t285 + t266 * t419) * pkin(4) + (-t154 - t260 - t461) * t282 + t332 * t281 + t330;
t509 = -g(2) + t15;
t174 = qJDD(4) * t267 - t266 * t394;
t319 = (-qJDD(1) * t286 - t288 * t290) * pkin(1);
t296 = t319 + (-t273 * t281 - t274 * t280) * pkin(2);
t464 = t267 * pkin(7);
t353 = -pkin(3) * t266 + t464;
t295 = -t280 * t188 + t281 * t353 + t296;
t362 = -t268 - t278;
t358 = -pkin(3) - t362;
t363 = pkin(4) * t285 + t254;
t404 = pkin(4) * t377 + qJD(5) * t267;
t326 = rSges(6,1) * t377 + t516 * rSges(6,2) + t266 * t420 + t404;
t463 = pkin(3) - t268;
t460 = t326 + (t267 * t463 - t328 + t465) * t282;
t14 = t281 * t503 + t460 * t282 - t363 * t174 + (-t281 * t358 + qJDD(5)) * t266 + (-pkin(4) * t419 - qJD(4) * t207 + qJD(5) * t282 - t281 * t462) * t267 + t295;
t508 = -g(3) + t14;
t458 = rSges(5,1) * t287;
t258 = -rSges(5,2) * t285 + t458;
t208 = t258 * qJD(4);
t234 = rSges(5,2) * t427;
t401 = t267 * rSges(5,3) + t234;
t331 = rSges(5,1) * t426 - t401;
t395 = qJD(4) * t267;
t95 = t282 * t329 + t380;
t26 = -t174 * t255 - t208 * t395 - t281 * t331 + t282 * t95 + t295;
t545 = t26 - g(3);
t355 = -t148 + t364;
t382 = t517 * rSges(5,1) - rSges(5,2) * t374;
t93 = t282 * t401 + t382;
t27 = t208 * t396 + t173 * t255 + (-t154 - t93) * t282 + t355 * t281 + t330;
t544 = t27 - g(2);
t349 = -rSges(4,1) * t266 - rSges(4,2) * t267;
t543 = t281 * t349 + t282 * (-rSges(4,1) * t425 + t221) + t296 - g(3);
t365 = -t187 - t467;
t403 = -rSges(4,1) * t428 - rSges(4,2) * t425;
t542 = t281 * t365 - t282 * t403 - g(2) + t330;
t541 = t553 * t266 + t559 * t267;
t422 = t273 * t282;
t179 = -rSges(3,1) * t421 + rSges(3,2) * t422;
t540 = t179 * t282 - t281 * t350 - g(3) + t319;
t197 = rSges(3,1) * t274 - t273 * rSges(3,2);
t539 = t178 * t282 - t197 * t281 - g(2) + t354;
t538 = t552 * t266 + t558 * t267;
t537 = t556 * qJD(4) + t551;
t536 = t557 * qJD(4) + t548;
t535 = -t564 * qJD(4) + t570 * t285 + t572 * t287;
t534 = -t565 * qJD(4) + t571 * t285 + t573 * t287;
t533 = t560 * t266 - t554 * t267;
t532 = t554 * t266 + t560 * t267;
t525 = t282 * t414;
t399 = t500 + t341;
t400 = t244 - t249;
t524 = (t285 * t399 + t287 * t400) * t282;
t397 = t499 + t342;
t398 = t246 - t251;
t523 = (t285 * t397 + t287 * t398) * t282;
t384 = t254 * t396 + t404;
t519 = -t356 + t384 + t525 - t326;
t518 = -t585 - t575;
t298 = t266 * t358 + t549;
t511 = -t503 + t468;
t504 = -t206 - t467;
t169 = rSges(6,1) * t427 + rSges(6,2) * t426;
t239 = pkin(4) * t427;
t502 = t239 + t169;
t501 = t257 + t279;
t130 = t282 * t331;
t391 = pkin(2) * t422;
t357 = -t282 * t353 + t391;
t302 = t255 * t395 + t130 + t357;
t63 = t302 + t388;
t387 = t288 * t455;
t64 = t282 * t355 + t183 - t387;
t497 = t266 * t64 + t267 * t63;
t406 = -Icges(5,2) * t423 + t146 - t232;
t410 = t267 * t499 + t142;
t482 = t285 * t406 + t287 * t410;
t407 = Icges(5,2) * t426 + t145 + t231;
t411 = -t266 * t499 + t141;
t481 = -t285 * t407 - t287 * t411;
t408 = -Icges(6,2) * t423 + t144 - t230;
t412 = t267 * t500 + t140;
t480 = t285 * t408 + t287 * t412;
t409 = Icges(6,2) * t426 + t143 + t229;
t413 = -t266 * t500 + t139;
t479 = -t285 * t409 - t287 * t413;
t478 = t173 / 0.2e1;
t477 = t174 / 0.2e1;
t471 = -rSges(5,3) - pkin(7);
t170 = t255 * t266;
t429 = t255 * t267;
t371 = -pkin(3) - t458;
t369 = -t396 / 0.2e1;
t368 = t396 / 0.2e1;
t367 = -t395 / 0.2e1;
t366 = t395 / 0.2e1;
t359 = t224 - t382;
t259 = rSges(2,1) * t288 - t286 * rSges(2,2);
t351 = rSges(2,1) * t286 + rSges(2,2) * t288;
t156 = t262 + t352;
t153 = -t197 * t282 - t387;
t318 = t388 + t391;
t155 = t349 - t468;
t300 = t267 * t148 + t266 * t331;
t299 = pkin(4) * t375 - t498;
t101 = -t389 + t237 + (-rSges(6,3) + t284) * t266 + t504;
t110 = t266 * t371 + t401 + t464 - t468;
t111 = t266 * t471 + t267 * t371 + t238 - t467;
t100 = t266 * t362 - t267 * t284 - t511;
t297 = t363 * t395 - t260 + t357 + (rSges(6,1) * t426 - t266 * t463 + t549) * t282;
t294 = t266 * t298 - t267 * t414;
t49 = t297 + t388;
t50 = t282 * t332 + t384 - t387;
t293 = (t50 * t511 - t49 * (-t328 + t504)) * t282;
t292 = (((t513 + t528 + t574) * t267 + ((-t512 + t518) * t267 - t514 - t529 + t530) * t266) * qJD(4) + t551) * t369 + (-t567 * qJD(4) + t568 * t285 + t569 * t287) * t282 + (t521 + t526) * t478 + (t522 + t527) * t477 + ((((t512 - t585) * t267 + t514 + t530) * t267 + (t518 * t266 + t513 - t531) * t266) * qJD(4) + t536 - t548) * t367 + (t532 + t535) * t366 + (Icges(3,3) + Icges(4,3) + t566) * t281 + (t533 + t534 + t537) * t368;
t291 = (t64 * (-t234 + t468) - t63 * (-t456 - t465 - t467) + (-t371 * t63 + t471 * t64) * t267) * t282;
t176 = t282 * t349;
t171 = t254 * t267;
t128 = t282 * t365 - t387;
t127 = -t176 + t318;
t71 = qJD(4) * t300 + qJD(3);
t44 = qJD(4) * t294 + qJD(3);
t25 = qJDD(3) + t174 * t148 + t173 * t331 + (-t266 * t95 + t267 * t93) * qJD(4);
t5 = qJDD(3) - t414 * t174 + t298 * t173 + (-t460 * t266 + t461 * t267) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t292 + (t539 * (-t197 - t469) + t540 * (-t350 - t470) + (t153 - t179 + t387) * t152) * m(3) + ((qJDD(1) * t351 + g(3)) * t351 + (qJDD(1) * t259 + g(2)) * t259) * m(2) + (t50 * (t299 + t388) + t293 + t509 * (t101 - t469) + t508 * (t100 - t470) + (-t50 + t519) * t49) * m(6) + (t64 * (t359 + t388) + t291 + (-t64 + t561) * t63 + t544 * (t111 - t469) + t545 * (t110 - t470)) * m(5) + (t128 * (t318 - t403) + t542 * (t156 - t469) + t543 * (t155 - t470) + (-t128 + t555) * t127) * m(4); t292 + (t293 + (-t297 + t299) * t50 + t519 * t49 + t509 * t101 + t508 * t100) * m(6) + (t291 + (-t302 + t359) * t64 + t561 * t63 + t544 * t111 + t545 * t110) * m(5) + (t542 * t156 + t543 * t155 + t555 * t127 + (t176 - t403) * t128) * m(4) + (-t152 * t179 + t153 * t178 + (-t152 * t282 - t539) * t197 - (t153 * t282 + t540) * t350) * m(3); m(4) * qJDD(3) + m(5) * t25 + m(6) * t5 + (-m(4) - m(5) - m(6)) * g(1); t557 * t478 + t556 * t477 + (t533 * t282 + t521 * t281 + t529 * t174 + t528 * t173 + (t546 * t266 + t547 * t267) * qJD(4)) * t266 / 0.2e1 + (t532 * t282 + t522 * t281 + t531 * t174 + t530 * t173 + (t538 * t266 + t541 * t267) * qJD(4)) * t267 / 0.2e1 + (t526 * t266 + t527 * t267) * t281 / 0.2e1 - (((t397 + t399) * t287 + (-t398 - t400) * t285) * t282 + (((t407 + t409) * t267 + (t406 + t408) * t266) * t287 + ((-t411 - t413) * t267 + (-t410 - t412) * t266) * t285) * qJD(4)) * t282 / 0.2e1 + ((t526 * t282 + t535) * t267 + (-t527 * t282 + t534) * t266) * t282 / 0.2e1 - t537 * t428 / 0.2e1 + t536 * t425 / 0.2e1 + ((-t396 * t433 + t432) * t266 + (-t524 + (t479 * t267 + (t157 - t480) * t266) * qJD(4)) * t267 + (-t396 * t431 + t430) * t266 + (-t523 + (t481 * t267 + (t159 - t482) * t266) * qJD(4)) * t267) * t369 + ((t528 * t282 + t547) * t267 + (-t529 * t282 + t546) * t266) * t368 + ((t157 * t395 + t432) * t267 + (t524 + (t480 * t266 + (-t433 - t479) * t267) * qJD(4)) * t266 + (t159 * t395 + t430) * t267 + (t523 + (t482 * t266 + (-t431 - t481) * t267) * qJD(4)) * t266) * t367 + ((t530 * t282 + t541) * t267 + (-t531 * t282 + t538) * t266) * t366 + (-g(1) * t501 - g(2) * t502 + t15 * t239 + t5 * t294 + (t169 * t282 - t254 * t428) * t49 + (-t171 * t282 + t254 * t425) * t50 + (-(t50 * t257 - t44 * t502) * qJD(4) + t44 * (-t460 + t525) + t15 * t254 + t50 * t207) * t266 + (-g(3) * (-t452 + (-rSges(6,1) - pkin(4)) * t285) - (t49 * t501 + (-pkin(4) * t424 - t171) * t44) * qJD(4) + t44 * (t298 * t282 + t461) - t14 * t363 - t49 * (-pkin(4) * t392 - t207)) * t267) * m(6) + (t25 * t300 + t71 * ((-t95 - t132) * t266 + (t93 + t130) * t267) + t27 * t170 - t26 * t429 + (t64 * t425 - t63 * t428) * t255 + t497 * t208 - (-t170 * t63 + t429 * t64) * t282 - (t71 * (-t170 * t266 - t267 * t429) + t497 * t258) * qJD(4) - g(1) * t258 - g(2) * t170 + g(3) * t429) * m(5); (t266 * t508 + t267 * t509) * m(6);];
tau = t1;
