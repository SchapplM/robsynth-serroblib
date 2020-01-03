% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:55
% EndTime: 2019-12-31 18:19:16
% DurationCPUTime: 14.11s
% Computational Cost: add. (46044->566), mult. (39707->826), div. (0->0), fcn. (42961->10), ass. (0->334)
t362 = qJ(3) + pkin(9);
t357 = sin(t362);
t363 = qJ(1) + pkin(8);
t358 = sin(t363);
t492 = t357 * t358;
t616 = Icges(4,3) + Icges(5,3);
t359 = cos(t362);
t364 = sin(qJ(5));
t367 = cos(qJ(5));
t404 = Icges(6,5) * t367 - Icges(6,6) * t364;
t264 = -Icges(6,3) * t359 + t357 * t404;
t506 = Icges(6,4) * t367;
t407 = -Icges(6,2) * t364 + t506;
t268 = -Icges(6,6) * t359 + t357 * t407;
t507 = Icges(6,4) * t364;
t408 = Icges(6,1) * t367 - t507;
t272 = -Icges(6,5) * t359 + t357 * t408;
t360 = cos(t363);
t479 = t360 * t367;
t489 = t358 * t364;
t291 = t359 * t489 + t479;
t481 = t360 * t364;
t487 = t358 * t367;
t292 = t359 * t487 - t481;
t138 = t264 * t492 - t268 * t291 + t272 * t292;
t199 = Icges(6,5) * t292 - Icges(6,6) * t291 + Icges(6,3) * t492;
t286 = Icges(6,4) * t292;
t202 = -Icges(6,2) * t291 + Icges(6,6) * t492 + t286;
t285 = Icges(6,4) * t291;
t206 = -Icges(6,1) * t292 - Icges(6,5) * t492 + t285;
t102 = t199 * t492 - t202 * t291 - t206 * t292;
t293 = -t359 * t481 + t487;
t294 = t359 * t479 + t489;
t491 = t357 * t360;
t201 = Icges(6,5) * t294 + Icges(6,6) * t293 + Icges(6,3) * t491;
t508 = Icges(6,4) * t294;
t204 = Icges(6,2) * t293 + Icges(6,6) * t491 + t508;
t287 = Icges(6,4) * t293;
t207 = Icges(6,1) * t294 + Icges(6,5) * t491 + t287;
t103 = t201 * t492 - t291 * t204 + t292 * t207;
t403 = t102 * t358 + t103 * t360;
t12 = t359 * t138 - t357 * t403;
t521 = rSges(6,1) * t367;
t411 = -rSges(6,2) * t364 + t521;
t276 = -rSges(6,3) * t359 + t357 * t411;
t317 = pkin(4) * t357 - pkin(7) * t359;
t365 = sin(qJ(3));
t529 = pkin(3) * t365;
t416 = t276 + t317 + t529;
t195 = t416 * t358;
t197 = t416 * t360;
t577 = t358 * t195 + t197 * t360;
t315 = rSges(5,1) * t357 + rSges(5,2) * t359;
t385 = t315 + t529;
t580 = t385 * t360;
t581 = t385 * t358;
t586 = t358 * t581 + t360 * t580;
t467 = -t577 * m(6) / 0.2e1 - m(5) * t586 / 0.2e1;
t544 = rSges(6,3) + pkin(7);
t183 = (t529 - t544 * t359 + (pkin(4) + t411) * t357) * t358;
t483 = t359 * t360;
t325 = pkin(7) * t483;
t448 = t357 * rSges(6,2) * t481 + rSges(6,3) * t483;
t184 = t325 + (-t529 + (-pkin(4) - t521) * t357) * t360 + t448;
t564 = m(6) / 0.2e1;
t565 = m(5) / 0.2e1;
t469 = (t183 * t358 - t184 * t360) * t564 + t586 * t565;
t38 = t469 - t467;
t524 = m(6) * qJD(5);
t574 = -t292 * rSges(6,1) + t291 * rSges(6,2);
t209 = rSges(6,3) * t492 - t574;
t278 = rSges(6,3) * t357 + t359 * t411;
t122 = (t278 * t358 - t209) * t357;
t412 = t294 * rSges(6,1) + t293 * rSges(6,2);
t211 = rSges(6,3) * t491 + t412;
t236 = -rSges(6,1) * t357 * t479 + t448;
t123 = (-t276 * t360 - t236) * t359 + (-t278 * t360 + t211) * t357;
t73 = t122 * t358 - t123 * t360;
t615 = -t38 * qJD(1) - t73 * t524 / 0.2e1;
t368 = cos(qJ(3));
t528 = pkin(3) * t368;
t354 = pkin(2) + t528;
t526 = -qJ(4) - pkin(6);
t341 = t358 * t526;
t530 = cos(qJ(1)) * pkin(1);
t427 = -t341 + t530;
t527 = pkin(4) * t359;
t571 = -t544 * t357 - t527;
t161 = (t354 - t571) * t360 + t412 + t427;
t222 = -rSges(6,1) * t291 - rSges(6,2) * t292;
t223 = rSges(6,1) * t293 - rSges(6,2) * t294;
t447 = -t358 * t354 - t360 * t526;
t531 = sin(qJ(1)) * pkin(1);
t417 = t447 - t531;
t590 = t571 * t358 + t417 + t574;
t614 = m(6) * (t161 * t223 - t222 * t590);
t309 = Icges(5,5) * t359 - Icges(5,6) * t357;
t327 = Icges(4,5) * t368 - Icges(4,6) * t365;
t612 = (t309 + t327) * t360 + t616 * t358;
t486 = t358 * t368;
t488 = t358 * t365;
t490 = t358 * t359;
t605 = -Icges(4,5) * t486 - Icges(5,5) * t490 + Icges(4,6) * t488 + Icges(5,6) * t492 + t616 * t360;
t140 = t264 * t491 + t268 * t293 + t272 * t294;
t104 = t199 * t491 + t293 * t202 - t206 * t294;
t105 = t201 * t491 + t293 * t204 + t294 * t207;
t402 = t104 * t358 + t105 * t360;
t611 = -t359 * t140 + t357 * t402;
t53 = -t102 * t360 + t103 * t358;
t164 = t211 * t359 + t276 * t491;
t606 = t209 * t359 + t276 * t492;
t578 = t164 * t358 - t360 * t606;
t499 = t199 * t359;
t600 = t202 * t364 + t206 * t367;
t119 = t600 * t357 + t499;
t509 = Icges(5,4) * t357;
t313 = Icges(5,1) * t359 - t509;
t258 = Icges(5,5) * t358 + t313 * t360;
t510 = Icges(4,4) * t365;
t331 = Icges(4,1) * t368 - t510;
t275 = Icges(4,5) * t358 + t331 * t360;
t608 = -t258 * t490 - t275 * t486;
t607 = -t104 * t360 + t105 * t358;
t604 = -Icges(4,5) * t365 - Icges(5,5) * t357 - Icges(4,6) * t368 - Icges(5,6) * t359;
t603 = t612 * t360 + t608;
t321 = Icges(5,4) * t492;
t257 = Icges(5,1) * t490 - Icges(5,5) * t360 - t321;
t338 = Icges(4,4) * t488;
t274 = Icges(4,1) * t486 - Icges(4,5) * t360 - t338;
t478 = t360 * t368;
t602 = -t257 * t483 - t274 * t478 + t605 * t358;
t588 = t258 * t483 + t275 * t478 + t612 * t358;
t255 = Icges(5,4) * t490 - Icges(5,2) * t492 - Icges(5,6) * t360;
t270 = Icges(4,4) * t486 - Icges(4,2) * t488 - Icges(4,6) * t360;
t601 = t255 * t357 + t270 * t365;
t579 = -t161 * t358 - t360 * t590;
t328 = Icges(4,2) * t368 + t510;
t592 = (t331 / 0.2e1 - t328 / 0.2e1) * t365;
t350 = Icges(5,4) * t359;
t503 = Icges(5,2) * t357;
t256 = Icges(5,6) * t358 + (t350 - t503) * t360;
t361 = Icges(4,4) * t368;
t329 = -Icges(4,2) * t365 + t361;
t271 = Icges(4,6) * t358 + t329 * t360;
t480 = t360 * t365;
t591 = -t256 * t491 - t271 * t480 + t588;
t589 = t256 * t357 + t271 * t365 + t605;
t587 = -t255 * t491 - t256 * t492 - t270 * t480 - t271 * t488 - t602 - t603;
t551 = t358 / 0.2e1;
t548 = -t360 / 0.2e1;
t546 = t360 / 0.2e1;
t582 = t276 * t358;
t576 = t604 * t358;
t575 = t604 * t360;
t355 = t358 ^ 2;
t356 = t360 ^ 2;
t445 = t355 + t356;
t330 = Icges(4,1) * t365 + t361;
t310 = Icges(5,2) * t359 + t509;
t282 = t310 * t360;
t511 = Icges(5,1) * t357;
t409 = -t350 - t511;
t284 = t409 * t360;
t570 = (-t258 + t282) * t492 + (-t256 + t284) * t490;
t298 = (-Icges(6,2) * t367 - t507) * t357;
t301 = (-Icges(6,1) * t364 - t506) * t357;
t569 = -(t272 / 0.2e1 + t298 / 0.2e1) * t364 + (t301 / 0.2e1 - t268 / 0.2e1) * t367;
t281 = -Icges(5,2) * t490 - t321;
t283 = t409 * t358;
t299 = -Icges(4,2) * t486 - t338;
t302 = t330 * t358;
t568 = (t257 + t281) * t357 + (t255 - t283) * t359 + (t270 + t302) * t368 + (t274 + t299) * t365;
t567 = 0.4e1 * qJD(1);
t566 = 2 * qJD(3);
t563 = -t611 / 0.2e1;
t562 = t53 / 0.2e1;
t561 = t607 / 0.2e1;
t216 = -Icges(6,5) * t291 - Icges(6,6) * t292;
t462 = -Icges(6,2) * t292 - t206 - t285;
t464 = -Icges(6,1) * t291 - t202 - t286;
t89 = -t216 * t359 + (-t462 * t364 + t464 * t367) * t357;
t560 = t89 / 0.2e1;
t557 = m(6) * (t122 * t590 + t123 * t161 - t164 * t184 + t183 * t606);
t304 = (-rSges(6,1) * t364 - rSges(6,2) * t367) * t357;
t555 = m(6) * (-t195 * t223 + t197 * t222 + t579 * t304);
t553 = m(6) * (t161 * t184 + t183 * t590);
t552 = t357 / 0.2e1;
t550 = t358 / 0.4e1;
t549 = -t359 / 0.2e1;
t547 = -t360 / 0.4e1;
t353 = t360 * pkin(6);
t523 = rSges(4,1) * t368;
t428 = pkin(2) + t523;
t446 = rSges(4,2) * t488 + t360 * rSges(4,3);
t224 = -t358 * t428 + t353 + t446 - t531;
t340 = rSges(4,2) * t480;
t225 = t530 - t340 + t428 * t360 + (rSges(4,3) + pkin(6)) * t358;
t332 = rSges(4,1) * t365 + rSges(4,2) * t368;
t305 = t332 * t358;
t306 = t332 * t360;
t543 = m(4) * (t224 * t305 - t225 * t306);
t394 = rSges(5,1) * t490 - rSges(5,2) * t492 - t360 * rSges(5,3);
t212 = -t394 + t417;
t425 = -rSges(5,2) * t491 + t358 * rSges(5,3);
t522 = rSges(5,1) * t359;
t213 = (t354 + t522) * t360 + t425 + t427;
t542 = m(5) * (t212 * t581 - t213 * t580);
t541 = m(5) * (t212 * t360 + t213 * t358);
t537 = m(6) * t579;
t536 = m(6) * t578;
t155 = t222 * t358 + t223 * t360;
t532 = m(6) * t155;
t519 = t358 * t12;
t517 = t360 * t611;
t498 = t201 * t359;
t495 = t264 * t359;
t295 = (-Icges(6,5) * t364 - Icges(6,6) * t367) * t357;
t484 = t359 * t295;
t477 = t364 * t268;
t269 = Icges(6,6) * t357 + t359 * t407;
t476 = t364 * t269;
t475 = t367 * t272;
t273 = Icges(6,5) * t357 + t359 * t408;
t474 = t367 * t273;
t472 = t73 * qJD(4);
t398 = t209 * t360 - t211 * t358;
t101 = t398 * t359 + (-t236 * t358 - t360 * t582) * t357;
t78 = 0.2e1 * (t101 / 0.4e1 - t155 / 0.4e1) * m(6);
t471 = t78 * qJD(2);
t463 = Icges(6,1) * t293 - t204 - t508;
t461 = -Icges(6,2) * t294 + t207 + t287;
t456 = -t358 * (pkin(2) * t358 - t353 + t447) + t360 * (-t358 * pkin(6) - t341 + (-pkin(2) + t354) * t360);
t451 = -t268 + t301;
t450 = t272 + t298;
t444 = qJD(1) * t357;
t443 = qJD(1) * t359;
t442 = qJD(5) * t357;
t441 = qJD(5) * t359;
t437 = t563 + t611 / 0.2e1;
t435 = t492 / 0.4e1;
t429 = t309 / 0.2e1 + t327 / 0.2e1;
t426 = rSges(5,2) * t357 - t522 - t528;
t418 = t445 * t529;
t318 = pkin(7) * t357 + t527;
t415 = -t278 - t318 - t528;
t115 = -t450 * t291 + t451 * t292 + t295 * t492;
t116 = t450 * t293 + t451 * t294 + t295 * t491;
t217 = Icges(6,5) * t293 - Icges(6,6) * t294;
t90 = -t217 * t359 + (-t461 * t364 + t463 * t367) * t357;
t414 = t555 / 0.2e1 + (t116 + t90) * t550 + (t115 + t89) * t547;
t399 = -t204 * t364 + t207 * t367;
t120 = t357 * t399 - t498;
t401 = -t119 * t358 + t120 * t360;
t396 = t475 - t477;
t74 = t216 * t492 - t462 * t291 + t464 * t292;
t75 = t217 * t492 - t461 * t291 + t463 * t292;
t33 = t358 * t75 - t360 * t74;
t76 = t216 * t491 + t462 * t293 + t464 * t294;
t77 = t217 * t491 + t461 * t293 + t463 * t294;
t34 = t358 * t77 - t360 * t76;
t393 = t33 * t548 + t34 * t551;
t390 = -t264 * t358 + t600;
t389 = -t264 * t360 - t399;
t265 = Icges(6,3) * t357 + t359 * t404;
t388 = t265 - t396;
t300 = t328 * t360;
t303 = t330 * t360;
t383 = (-t271 - t303) * t368 + (-t275 + t300) * t365;
t382 = t12 * t550 + t611 * t547 - t519 / 0.4e1 + t517 / 0.4e1 + (t435 - t492 / 0.4e1) * t607;
t381 = -t607 / 0.2e1 + t561 - t588 * t358 / 0.2e1 + t591 * t551 + ((t601 + t612) * t360 + t587 + t602 + t608) * t548;
t380 = t562 - t53 / 0.2e1 + (t605 * t360 + (t257 * t359 + t274 * t368 - t601) * t358) * t548 + (t589 * t360 - t588 + t591) * t546 + (t589 * t358 + t587 + t603) * t551;
t121 = -t388 * t359 + (t264 + t474 - t476) * t357;
t165 = t357 * t396 - t495;
t231 = t268 * t358;
t233 = t272 * t358;
t84 = -t390 * t359 + (t231 * t364 - t233 * t367 + t199) * t357;
t232 = t268 * t360;
t234 = t272 * t360;
t85 = -t389 * t359 + (t232 * t364 - t234 * t367 + t201) * t357;
t373 = t357 * t388 + t495;
t98 = -t269 * t291 + t273 * t292 + t358 * t373;
t99 = t269 * t293 + t273 * t294 + t360 * t373;
t376 = t121 * t549 + t165 * t552 + t557 / 0.2e1 + (t84 + t98) * t435 + (t85 + t99) * t491 / 0.4e1 + (-t119 + t138) * t490 / 0.4e1 + (t120 + t140) * t483 / 0.4e1;
t375 = t357 * t390 + t499;
t374 = t357 * t389 + t498;
t372 = t475 / 0.2e1 - t477 / 0.2e1 - t265 / 0.2e1 + t350 + t511 / 0.2e1 - t503 / 0.2e1;
t371 = t474 / 0.2e1 - t476 / 0.2e1 + t264 / 0.2e1 + t313 / 0.2e1 - t310 / 0.2e1;
t334 = -rSges(4,2) * t365 + t523;
t263 = t426 * t360;
t261 = t426 * t358;
t215 = -t305 * t358 - t306 * t360;
t198 = t415 * t360;
t196 = t415 * t358;
t179 = -t445 * t315 - t418;
t171 = -t223 * t359 - t304 * t491;
t170 = t222 * t359 + t304 * t492;
t151 = -t532 / 0.2e1;
t149 = (t222 * t360 - t223 * t358) * t357;
t135 = t398 * t357;
t129 = -t484 + (-t450 * t364 + t451 * t367) * t357;
t125 = (-pkin(4) * t491 + t236 + t325) * t360 - t418 + (-t317 * t358 - t582) * t358;
t107 = -t536 / 0.2e1;
t92 = (t318 * t358 + t209) * t358 + (t318 * t360 + t211) * t360 + t456;
t83 = -t537 + t541;
t79 = (t101 + t155) * t564;
t72 = t73 * qJD(3) * t564;
t71 = -t232 * t293 - t234 * t294 + t360 * t374;
t70 = -t231 * t293 - t233 * t294 + t360 * t375;
t69 = t232 * t291 - t234 * t292 + t358 * t374;
t68 = t231 * t291 - t233 * t292 + t358 * t375;
t61 = -t484 / 0.2e1 + t614 + t569 * t357;
t48 = -t165 * t359 + t357 * t401;
t47 = t155 * t92 + t577 * t304;
t39 = t467 + t469;
t32 = t107 + t532 / 0.2e1;
t31 = t151 + t107;
t30 = t151 + t536 / 0.2e1;
t28 = t358 * t71 - t360 * t70;
t27 = t358 * t69 - t360 * t68;
t24 = t101 * t135 + t122 * t606 - t123 * t164;
t23 = -t116 * t359 + (t358 * t76 + t360 * t77) * t357;
t22 = -t115 * t359 + (t358 * t74 + t360 * t75) * t357;
t21 = (t330 / 0.2e1 + t329 / 0.2e1) * t368 + t592 + t543 + t542 + t553 + t372 * t359 + t371 * t357;
t16 = (-t121 + t401) * t359 + (t84 * t358 + t85 * t360 + t165) * t357;
t9 = (t402 - t99) * t359 + (t358 * t70 + t360 * t71 + t140) * t357;
t8 = (t403 - t98) * t359 + (t358 * t68 + t360 * t69 + t138) * t357;
t7 = m(6) * t47 + t393;
t6 = t437 * t492;
t5 = m(6) * t24 + (t517 / 0.2e1 - t519 / 0.2e1 - t16 / 0.2e1) * t359 + (t9 * t546 + t8 * t551 + t48 / 0.2e1) * t357;
t4 = t358 * t380 + t360 * t381;
t3 = t376 + t414;
t2 = -t557 / 0.2e1 + t382 + (-t165 / 0.2e1 + (-t99 / 0.4e1 - t85 / 0.4e1) * t360 + (-t98 / 0.4e1 - t84 / 0.4e1) * t358) * t357 + (t121 / 0.2e1 + (-t140 / 0.4e1 - t120 / 0.4e1) * t360 + (-t138 / 0.4e1 + t119 / 0.4e1) * t358) * t359 + t414;
t1 = -t555 / 0.2e1 + (-t116 / 0.4e1 - t90 / 0.4e1) * t358 + (t115 / 0.4e1 + t89 / 0.4e1) * t360 + t382 + t376;
t10 = [t21 * qJD(3) + t83 * qJD(4) + t61 * qJD(5), 0, t21 * qJD(1) + t39 * qJD(4) + t3 * qJD(5) + ((t212 * t263 + t213 * t261) * t565 + (t161 * t196 - t183 * t197 - t184 * t195 + t198 * t590) * t564) * t566 + ((m(4) * (-t224 * t334 - t305 * t332) - t84 / 0.2e1 - t98 / 0.2e1 + t429 * t360 + (-t274 / 0.2e1 - t299 / 0.2e1) * t368 + (t270 / 0.2e1 + t302 / 0.2e1) * t365 + (-t257 / 0.2e1 - t281 / 0.2e1) * t359 + (t255 / 0.2e1 - t283 / 0.2e1) * t357 - t381) * t360 + (m(4) * (-t225 * t334 + t306 * t332) + t85 / 0.2e1 + t99 / 0.2e1 + t429 * t358 + (t275 / 0.2e1 - t300 / 0.2e1) * t368 + (-t271 / 0.2e1 - t303 / 0.2e1) * t365 + (t258 / 0.2e1 - t282 / 0.2e1) * t359 + (-t256 / 0.2e1 + t284 / 0.2e1) * t357 - t380) * t358) * qJD(3), qJD(1) * t83 + qJD(3) * t39 + qJD(5) * t31, t61 * qJD(1) + t3 * qJD(3) + t31 * qJD(4) - t129 * t441 + (t161 * t171 - t164 * t223 + t170 * t590 - t222 * t606) * t524 + ((t90 / 0.2e1 + t116 / 0.2e1) * t360 + (t560 + t115 / 0.2e1 - t437) * t358) * t442; 0, 0, t79 * qJD(5) + (m(4) * t215 / 0.2e1 + t179 * t565 + t125 * t564) * t566, 0, t79 * qJD(3) + t149 * t524; t4 * qJD(3) - t38 * qJD(4) + t2 * qJD(5) + (-t553 / 0.4e1 - t542 / 0.4e1 - t543 / 0.4e1) * t567 - t372 * t443 - t371 * t444 + (-(t330 + t329) * t368 / 0.2e1 - t592) * qJD(1), -qJD(5) * t78, t4 * qJD(1) + t7 * qJD(5) + (m(6) * (t125 * t92 - t195 * t196 - t197 * t198) + m(5) * ((t358 * t394 + t360 * (rSges(5,1) * t483 + t425) + t456) * t179 - t581 * t261 - t580 * t263) + m(4) * (t332 * t334 * t445 + (t358 * (rSges(4,1) * t486 - t446) + t360 * (rSges(4,1) * t478 + t358 * rSges(4,3) - t340)) * t215) + (t28 + (t568 * t360 + (t383 - t576) * t358 + t570) * t360 + t575 * t355) * t551 + (t27 + (t383 * t358 + (t568 - t575) * t360 + t570) * t358 + t576 * t356) * t548) * qJD(3), t615, t2 * qJD(1) - t471 + t7 * qJD(3) + (t22 * t548 + t23 * t551) * qJD(5) + (t16 / 0.2e1 + (t560 + t563) * t360 + (-t90 / 0.2e1 + t12 / 0.2e1) * t358) * t441 + (-t472 / 0.2e1 + (t135 * t155 + t149 * t92 - t170 * t197 - t171 * t195 + t578 * t304 - t24) * qJD(5)) * m(6) + (-t48 / 0.2e1 + (t34 / 0.2e1 - t9 / 0.2e1) * t360 + (t33 / 0.2e1 - t8 / 0.2e1) * t358) * t442; t38 * qJD(3) + t30 * qJD(5) + (t537 / 0.4e1 - t541 / 0.4e1) * t567, 0, ((-t196 * t360 + t198 * t358) * t564 + (-t261 * t360 + t263 * t358) * t565) * t566 - t615, 0, t30 * qJD(1) + t72 + (t170 * t358 - t171 * t360) * t524; t295 * t443 / 0.2e1 + t1 * qJD(3) + t32 * qJD(4) + t6 * qJD(5) - qJD(1) * t614 - t569 * t444, qJD(3) * t78, t1 * qJD(1) + t471 + (t8 * t548 + t9 * t551 + t490 * t562 + t27 * t492 / 0.2e1 + (t119 * t360 + t120 * t358) * t552 + (t85 * t358 - t84 * t360) * t549 + t483 * t561 + t28 * t491 / 0.2e1 - t393) * qJD(3) + t5 * qJD(5) + ((t101 * t92 - t122 * t197 - t123 * t195 + t125 * t135 - t164 * t196 + t198 * t606 - t47) * qJD(3) + t472 / 0.2e1) * m(6), qJD(1) * t32 + t72, t6 * qJD(1) + t5 * qJD(3) + (m(6) * (t135 * t149 - t164 * t171 + t170 * t606) + t359 ^ 2 * t129 / 0.2e1 + (t23 * t546 + t22 * t551 + (t358 * t89 + t360 * t90) * t549) * t357) * qJD(5);];
Cq = t10;
