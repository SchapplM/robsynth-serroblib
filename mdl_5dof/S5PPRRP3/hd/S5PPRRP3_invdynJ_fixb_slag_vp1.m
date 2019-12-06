% Calculate vector of inverse dynamics joint torques for
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:34
% EndTime: 2019-12-05 15:11:15
% DurationCPUTime: 36.77s
% Computational Cost: add. (15373->783), mult. (42689->1134), div. (0->0), fcn. (46982->8), ass. (0->293)
t546 = Icges(5,4) - Icges(6,5);
t526 = Icges(5,1) + Icges(6,1);
t525 = Icges(6,4) + Icges(5,5);
t524 = Icges(5,2) + Icges(6,3);
t536 = Icges(6,2) + Icges(5,3);
t545 = Icges(5,6) - Icges(6,6);
t331 = cos(pkin(8));
t330 = sin(pkin(7));
t335 = cos(qJ(3));
t440 = t330 * t335;
t399 = t331 * t440;
t332 = cos(pkin(7));
t334 = sin(qJ(3));
t439 = t332 * t334;
t303 = t399 - t439;
t329 = sin(pkin(8));
t333 = sin(qJ(4));
t444 = t329 * t333;
t463 = cos(qJ(4));
t243 = t303 * t463 + t330 * t444;
t558 = t546 * t243;
t438 = t335 * t332;
t441 = t330 * t334;
t305 = t331 * t438 + t441;
t245 = t305 * t463 + t332 * t444;
t557 = t546 * t245;
t394 = t329 * t463;
t379 = t335 * t394;
t307 = -t331 * t333 + t379;
t556 = t546 * t307;
t442 = t329 * t335;
t306 = t331 * t463 + t333 * t442;
t555 = t546 * t306;
t353 = -t305 * t333 + t332 * t394;
t554 = t546 * t353;
t354 = -t303 * t333 + t330 * t394;
t553 = t546 * t354;
t302 = t331 * t441 + t438;
t280 = t302 * qJD(3);
t164 = qJD(4) * t243 - t280 * t333;
t165 = qJD(4) * t354 - t280 * t463;
t413 = qJD(3) * t334;
t281 = qJD(3) * t399 - t332 * t413;
t552 = -t545 * t164 + t525 * t165 + t536 * t281;
t304 = t331 * t439 - t440;
t282 = t304 * qJD(3);
t166 = qJD(4) * t245 - t282 * t333;
t167 = qJD(4) * t353 - t282 * t463;
t283 = t305 * qJD(3);
t551 = -t545 * t166 + t525 * t167 + t536 * t283;
t550 = t524 * t164 - t546 * t165 - t545 * t281;
t549 = t524 * t166 - t546 * t167 - t545 * t283;
t548 = -t546 * t164 + t526 * t165 + t525 * t281;
t547 = -t546 * t166 + t526 * t167 + t525 * t283;
t544 = -t545 * t302 - t524 * t354 - t558;
t543 = -t545 * t304 - t524 * t353 - t557;
t502 = t525 * t243 + t536 * t302 + t545 * t354;
t501 = t525 * t245 + t536 * t304 + t545 * t353;
t542 = t526 * t243 + t525 * t302 + t553;
t541 = t526 * t245 + t525 * t304 + t554;
t246 = -qJD(4) * t379 + (qJD(4) * t331 + t329 * t413) * t333;
t380 = t334 * t394;
t247 = -qJD(3) * t380 - qJD(4) * t306;
t412 = qJD(3) * t335;
t391 = t329 * t412;
t540 = -t545 * t246 - t525 * t247 - t536 * t391;
t539 = t524 * t246 + t546 * t247 + t545 * t391;
t538 = -t546 * t246 - t526 * t247 - t525 * t391;
t443 = t329 * t334;
t537 = t524 * t306 - t545 * t443 - t556;
t496 = -t545 * t306 + t525 * t307 + t536 * t443;
t530 = -t526 * t307 - t525 * t443 + t555;
t519 = rSges(6,1) + pkin(4);
t516 = t544 * t164 + t542 * t165 + t548 * t243 + t502 * t281 + t552 * t302 - t550 * t354;
t515 = t543 * t164 + t541 * t165 + t547 * t243 + t501 * t281 + t551 * t302 - t549 * t354;
t514 = t544 * t166 + t542 * t167 + t548 * t245 + t502 * t283 + t552 * t304 - t550 * t353;
t513 = t543 * t166 + t541 * t167 + t547 * t245 + t501 * t283 + t551 * t304 - t549 * t353;
t511 = (t552 * t334 + t502 * t412) * t329 + t548 * t307 + t550 * t306 + t542 * t247 - t544 * t246;
t510 = (t551 * t334 + t501 * t412) * t329 + t547 * t307 + t549 * t306 + t541 * t247 - t543 * t246;
t509 = -t537 * t164 + t530 * t165 + t538 * t243 - t496 * t281 + t540 * t302 - t539 * t354;
t508 = -t537 * t166 + t530 * t167 + t538 * t245 - t496 * t283 + t540 * t304 - t539 * t353;
t507 = (t540 * t334 - t496 * t412) * t329 + t538 * t307 + t539 * t306 + t530 * t247 + t537 * t246;
t491 = t542 * t243 + t502 * t302 - t544 * t354;
t490 = t541 * t243 + t501 * t302 - t543 * t354;
t489 = t542 * t245 + t502 * t304 - t544 * t353;
t488 = t541 * t245 + t501 * t304 - t543 * t353;
t487 = t544 * t306 + t542 * t307 + t502 * t443;
t486 = t543 * t306 + t541 * t307 + t501 * t443;
t506 = t530 * t243 - t496 * t302 + t537 * t354;
t505 = t530 * t245 - t496 * t304 + t537 * t353;
t504 = -t537 * t306 + t530 * t307 - t496 * t443;
t529 = -t524 * t333 + t546 * t463;
t528 = t545 * t333 - t525 * t463;
t527 = t546 * t333 - t526 * t463;
t406 = qJDD(3) * t329;
t387 = t330 * t406;
t169 = qJD(4) * t281 + qJDD(4) * t302 + t387;
t386 = t332 * t406;
t170 = qJD(4) * t283 + qJDD(4) * t304 + t386;
t415 = qJD(3) * t329;
t393 = t330 * t415;
t256 = qJD(4) * t302 + t393;
t392 = t332 * t415;
t257 = qJD(4) * t304 + t392;
t405 = qJDD(3) * t331;
t411 = qJD(4) * t335;
t263 = -t405 + (qJD(3) * t411 + qJDD(4) * t334) * t329;
t414 = qJD(3) * t331;
t314 = qJD(4) * t443 - t414;
t522 = t491 * t169 + t490 * t170 + t516 * t256 + t515 * t257 - t506 * t263 - t509 * t314;
t521 = t489 * t169 + t488 * t170 + t514 * t256 + t513 * t257 - t505 * t263 - t508 * t314;
t520 = t487 * t169 + t486 * t170 + t511 * t256 + t510 * t257 - t504 * t263 - t507 * t314;
t518 = t491 * t256 + t490 * t257 - t506 * t314;
t517 = t489 * t256 + t488 * t257 - t505 * t314;
t512 = t487 * t256 + t486 * t257 - t504 * t314;
t503 = rSges(6,3) + qJ(5);
t446 = t329 * t330;
t445 = t329 * t332;
t500 = t529 * t302 - t303 * t545;
t499 = t529 * t304 - t305 * t545;
t498 = t527 * t302 + t525 * t303;
t497 = t527 * t304 + t525 * t305;
t495 = (t529 * t334 - t335 * t545) * t329;
t494 = (t527 * t334 + t525 * t335) * t329;
t493 = (t530 * t463 - t537 * t333 + (t528 * t334 + t536 * t335) * t329) * t314 + (t528 * t304 + t536 * t305 - t333 * t543 - t463 * t541) * t257 + (t528 * t302 + t536 * t303 - t333 * t544 - t463 * t542) * t256;
t492 = -t333 * t503 - t519 * t463;
t482 = (t307 * t524 + t530 + t555) * t314 + (t245 * t524 - t541 - t554) * t257 + (t243 * t524 - t542 - t553) * t256;
t481 = (-t306 * t526 + t537 - t556) * t314 + (t353 * t526 + t543 - t557) * t257 + (t354 * t526 + t544 - t558) * t256;
t480 = (-t306 * t525 - t307 * t545) * t314 + (-t245 * t545 + t353 * t525) * t257 + (-t243 * t545 + t354 * t525) * t256;
t308 = (-Icges(4,5) * t334 - Icges(4,6) * t335) * t329;
t284 = qJD(3) * t308;
t219 = pkin(3) * t303 + pkin(6) * t302;
t313 = (pkin(3) * t335 + pkin(6) * t334) * t329;
t328 = qJD(2) * t330;
t397 = t219 * t414 + t313 * t393 + t328;
t409 = qJD(5) * t353;
t424 = rSges(6,2) * t443 + t306 * t503 + t307 * t519;
t437 = rSges(6,2) * t302 + t243 * t519 - t354 * t503;
t38 = t256 * t424 - t314 * t437 + t397 - t409;
t410 = qJD(5) * t354;
t461 = rSges(6,2) * t281 + t164 * t503 + t165 * t519 - t410;
t204 = -pkin(3) * t280 + pkin(6) * t281;
t294 = (-pkin(3) * t334 + pkin(6) * t335) * t415;
t327 = qJDD(2) * t330;
t368 = t204 * t414 + t219 * t405 + t294 * t393 + t313 * t387 + t327;
t408 = qJD(5) * t306;
t434 = rSges(6,2) * t391 - t246 * t503 + t247 * t519 + t408;
t8 = qJD(5) * t166 - qJDD(5) * t353 + t169 * t424 + t256 * t434 - t263 * t437 - t314 * t461 + t368;
t479 = t38 * t461 + t437 * t8;
t478 = t169 / 0.2e1;
t477 = t170 / 0.2e1;
t476 = -t256 / 0.2e1;
t475 = t256 / 0.2e1;
t474 = -t257 / 0.2e1;
t473 = t257 / 0.2e1;
t472 = t263 / 0.2e1;
t467 = -t314 / 0.2e1;
t466 = t314 / 0.2e1;
t460 = rSges(6,2) * t283 + t166 * t503 + t167 * t519 - t409;
t205 = -pkin(3) * t282 + pkin(6) * t283;
t95 = rSges(5,1) * t167 - rSges(5,2) * t166 + rSges(5,3) * t283;
t459 = -t205 - t95;
t458 = Icges(4,4) * t303;
t457 = Icges(4,4) * t305;
t456 = Icges(4,4) * t334;
t455 = Icges(4,4) * t335;
t448 = t302 * t333;
t447 = t304 * t333;
t436 = rSges(6,2) * t304 + t245 * t519 - t353 * t503;
t117 = rSges(5,1) * t245 + rSges(5,2) * t353 + rSges(5,3) * t304;
t221 = pkin(3) * t305 + pkin(6) * t304;
t435 = -t117 - t221;
t433 = t243 * t503 + t354 * t519;
t432 = t245 * t503 + t353 * t519;
t396 = t302 * t463;
t431 = t303 * rSges(6,2) - t396 * t519 - t448 * t503;
t395 = t304 * t463;
t430 = t305 * rSges(6,2) - t395 * t519 - t447 * t503;
t429 = t331 * t204 + t294 * t446;
t176 = -Icges(4,2) * t302 + Icges(4,6) * t446 + t458;
t428 = -Icges(4,1) * t302 - t176 - t458;
t177 = -Icges(4,2) * t304 + Icges(4,6) * t445 + t457;
t427 = -Icges(4,1) * t304 - t177 - t457;
t288 = Icges(4,4) * t302;
t178 = Icges(4,1) * t303 + Icges(4,5) * t446 - t288;
t426 = Icges(4,2) * t303 - t178 + t288;
t289 = Icges(4,4) * t304;
t179 = Icges(4,1) * t305 + Icges(4,5) * t445 - t289;
t425 = Icges(4,2) * t305 - t179 + t289;
t218 = -t302 * pkin(3) + pkin(6) * t303;
t323 = pkin(6) * t442;
t312 = -pkin(3) * t443 + t323;
t423 = t218 * t414 + t312 * t393;
t422 = t331 * t219 + t313 * t446;
t421 = -t306 * t519 + t307 * t503;
t322 = rSges(6,2) * t442;
t420 = t443 * t492 + t322;
t268 = -Icges(4,6) * t331 + (-Icges(4,2) * t334 + t455) * t329;
t310 = (-Icges(4,1) * t334 - t455) * t329;
t419 = t268 - t310;
t269 = -Icges(4,5) * t331 + (Icges(4,1) * t335 - t456) * t329;
t309 = (-Icges(4,2) * t335 - t456) * t329;
t418 = t269 + t309;
t400 = t333 * t443;
t417 = rSges(5,2) * t400 + rSges(5,3) * t442;
t416 = qJD(2) * t332;
t407 = qJD(5) * t333;
t404 = -t205 - t460;
t402 = -m(3) - m(4) - m(5) - m(6);
t398 = -t221 - t436;
t385 = t415 / 0.2e1;
t220 = -t304 * pkin(3) + pkin(6) * t305;
t115 = rSges(5,1) * t243 + rSges(5,2) * t354 + rSges(5,3) * t302;
t127 = t247 * rSges(5,1) + t246 * rSges(5,2) + rSges(5,3) * t391;
t192 = rSges(5,1) * t307 - rSges(5,2) * t306 + rSges(5,3) * t443;
t93 = rSges(5,1) * t165 - rSges(5,2) * t164 + rSges(5,3) * t281;
t29 = -t115 * t263 + t127 * t256 + t169 * t192 - t314 * t93 + t368;
t59 = -t115 * t314 + t192 * t256 + t397;
t376 = t29 * t115 + t59 * t93;
t194 = -Icges(4,5) * t280 - Icges(4,6) * t281;
t195 = -Icges(4,5) * t282 - Icges(4,6) * t283;
t196 = -Icges(4,4) * t280 - Icges(4,2) * t281;
t197 = -Icges(4,4) * t282 - Icges(4,2) * t283;
t198 = -Icges(4,1) * t280 - Icges(4,4) * t281;
t199 = -Icges(4,1) * t282 - Icges(4,4) * t283;
t375 = t330 * (-t176 * t281 - t178 * t280 + t194 * t446 - t196 * t302 + t198 * t303) + t332 * (-t177 * t281 - t179 * t280 + t195 * t446 - t197 * t302 + t199 * t303);
t374 = t330 * (-t176 * t283 - t178 * t282 + t194 * t445 - t196 * t304 + t198 * t305) + t332 * (-t177 * t283 - t179 * t282 + t195 * t445 - t197 * t304 + t199 * t305);
t373 = t330 * (-t331 * t194 + (-t196 * t334 + t198 * t335 + (-t176 * t335 - t178 * t334) * qJD(3)) * t329) + t332 * (-t331 * t195 + (-t197 * t334 + t199 * t335 + (-t177 * t335 - t179 * t334) * qJD(3)) * t329);
t174 = Icges(4,5) * t303 - Icges(4,6) * t302 + Icges(4,3) * t446;
t175 = Icges(4,5) * t305 - Icges(4,6) * t304 + Icges(4,3) * t445;
t372 = t330 * (t174 * t446 - t176 * t302 + t178 * t303) + t332 * (t175 * t446 - t177 * t302 + t179 * t303);
t371 = t330 * (t174 * t445 - t176 * t304 + t178 * t305) + t332 * (t175 * t445 - t177 * t304 + t179 * t305);
t200 = -rSges(4,1) * t280 - rSges(4,2) * t281;
t311 = (-rSges(4,1) * t334 - rSges(4,2) * t335) * t329;
t287 = qJD(3) * t311;
t181 = rSges(4,1) * t303 - rSges(4,2) * t302 + rSges(4,3) * t446;
t279 = -t331 * rSges(4,3) + (rSges(4,1) * t335 - rSges(4,2) * t334) * t329;
t362 = t181 * t331 + t279 * t446;
t74 = t327 + t362 * qJDD(3) + (t200 * t331 + t287 * t446) * qJD(3);
t182 = rSges(4,1) * t305 - rSges(4,2) * t304 + rSges(4,3) * t445;
t201 = -rSges(4,1) * t282 - rSges(4,2) * t283;
t75 = (-qJD(3) * t201 - qJDD(3) * t182) * t331 + (-qJDD(2) + (-qJD(3) * t287 - qJDD(3) * t279) * t329) * t332;
t370 = t330 * t74 - t332 * t75;
t369 = t330 * (-t331 * t174 + (-t176 * t334 + t178 * t335) * t329) + t332 * (-t331 * t175 + (-t177 * t334 + t179 * t335) * t329);
t100 = qJD(3) * t362 + t328;
t101 = -t416 + (-t182 * t331 - t279 * t445) * qJD(3);
t367 = t100 * t330 - t101 * t332;
t366 = t181 * t332 - t182 * t330;
t365 = t200 * t332 - t201 * t330;
t161 = -rSges(5,1) * t396 + rSges(5,2) * t448 + t303 * rSges(5,3);
t163 = -rSges(5,1) * t395 + rSges(5,2) * t447 + t305 * rSges(5,3);
t364 = -qJD(3) * t205 - qJDD(3) * t221;
t363 = t218 * t392 - t220 * t393;
t355 = t219 * t392 - t221 * t393 + qJD(1);
t35 = -t256 * t436 + t257 * t437 + t355 + t408;
t343 = t204 * t392 + t219 * t386 + t364 * t446 + qJDD(1);
t7 = -qJD(5) * t246 + qJDD(5) * t306 - t169 * t436 + t170 * t437 - t256 * t460 + t257 * t461 + t343;
t352 = t35 * t461 + t437 * t7;
t349 = t38 * t434 + t424 * t8;
t348 = -t35 * t436 + t38 * t424;
t344 = (-t221 * t331 - t313 * t445) * qJD(3) - t416;
t39 = -t257 * t424 + t314 * t436 + t344 - t410;
t347 = t35 * t437 - t39 * t424;
t346 = -t220 * t414 - t312 * t392;
t345 = t335 * (-t38 * t437 + t39 * t436);
t338 = (-qJDD(2) + (-qJD(3) * t294 - qJDD(3) * t313) * t329) * t332 + t364 * t331;
t286 = qJD(3) * t310;
t285 = qJD(3) * t309;
t267 = -Icges(4,3) * t331 + (Icges(4,5) * t335 - Icges(4,6) * t334) * t329;
t255 = -rSges(5,1) * t380 + t417;
t230 = -rSges(5,1) * t306 - rSges(5,2) * t307;
t217 = -rSges(4,1) * t304 - rSges(4,2) * t305;
t216 = -rSges(4,1) * t302 - rSges(4,2) * t303;
t211 = -Icges(4,5) * t304 - Icges(4,6) * t305;
t210 = -Icges(4,5) * t302 - Icges(4,6) * t303;
t193 = t219 * t445;
t171 = t204 * t445;
t146 = rSges(5,1) * t353 - rSges(5,2) * t245;
t142 = rSges(5,1) * t354 - rSges(5,2) * t243;
t119 = -t331 * t267 + (-t268 * t334 + t269 * t335) * t329;
t99 = t267 * t445 - t268 * t304 + t269 * t305;
t98 = t267 * t446 - t268 * t302 + t269 * t303;
t97 = t366 * t415 + qJD(1);
t96 = -t331 * t284 + (-t285 * t334 + t286 * t335 + (-t268 * t335 - t269 * t334) * qJD(3)) * t329;
t67 = -t268 * t283 - t269 * t282 + t284 * t445 - t285 * t304 + t286 * t305;
t66 = -t268 * t281 - t269 * t280 + t284 * t446 - t285 * t302 + t286 * t303;
t65 = qJDD(1) + (qJD(3) * t365 + qJDD(3) * t366) * t329;
t60 = t117 * t314 - t192 * t257 + t344;
t48 = t115 * t257 - t117 * t256 + t355;
t30 = t117 * t263 - t127 * t257 - t170 * t192 + t314 * t95 + t338;
t24 = t115 * t170 - t117 * t169 - t256 * t95 + t257 * t93 + t343;
t9 = qJD(5) * t164 - qJDD(5) * t354 - t170 * t424 - t257 * t434 + t263 * t436 + t314 * t460 + t338;
t1 = [(m(2) + m(3)) * qJDD(1) + m(4) * t65 + m(5) * t24 + m(6) * t7 + (-m(2) + t402) * g(3); t402 * (g(1) * t330 - g(2) * t332) + m(4) * t370 + m(5) * (t29 * t330 - t30 * t332) + m(6) * (t330 * t8 - t332 * t9) + m(3) * (t330 ^ 2 + t332 ^ 2) * qJDD(2); (t331 ^ 2 * t284 + (((t330 * t428 + t332 * t427) * t335 + (t330 * t426 + t332 * t425) * t334) * t329 + (-t210 * t330 - t211 * t332 + t334 * t418 + t335 * t419) * t331) * t415) * t414 / 0.2e1 - (-t119 * t331 + t369 * t329) * t405 / 0.2e1 - (t329 * t373 - t331 * t96) * t414 / 0.2e1 + (t506 * t331 + (t330 * t491 + t332 * t490) * t329) * t478 + (t505 * t331 + (t330 * t489 + t332 * t488) * t329) * t477 + ((t243 * t494 + t303 * t496 - t354 * t495) * t314 + t493 * t302 + (t243 * t497 + t303 * t501 - t354 * t499) * t257 + (t243 * t498 + t303 * t502 - t354 * t500) * t256 + (t303 * t491 + t305 * t490 - t442 * t506) * qJD(4)) * t476 + (t509 * t331 + (t330 * t516 + t332 * t515) * t329) * t475 + ((t245 * t494 + t305 * t496 - t353 * t495) * t314 + t493 * t304 + (t245 * t497 + t305 * t501 - t353 * t499) * t257 + (t245 * t498 + t305 * t502 - t353 * t500) * t256 + (t303 * t489 + t305 * t488 - t442 * t505) * qJD(4)) * t474 + (t508 * t331 + (t330 * t514 + t332 * t513) * t329) * t473 + (t504 * t331 + (t330 * t487 + t332 * t486) * t329) * t472 + ((t306 * t495 + t307 * t494) * t314 + (t306 * t499 + t307 * t497) * t257 + (t306 * t500 + t307 * t498) * t256 + (t303 * t487 + t305 * t486) * qJD(4) + ((-qJD(4) * t504 + t256 * t502 + t257 * t501 + t314 * t496) * t335 + t493 * t334) * t329) * t467 + (t507 * t331 + (t330 * t511 + t332 * t510) * t329) * t466 + (t332 * (t329 * t374 - t331 * t67) + t330 * (t329 * t375 - t331 * t66)) * t385 + (t332 * (t329 * t371 - t331 * t99) + t330 * (t329 * t372 - t331 * t98)) * t406 / 0.2e1 - ((-qJD(3) * t96 - qJDD(3) * t119) * t331 + (qJD(3) * t373 + qJDD(3) * t369) * t329 + t520) * t331 / 0.2e1 + ((-qJD(3) * t66 - qJDD(3) * t98) * t331 + (qJD(3) * t375 + qJDD(3) * t372) * t329 + t522) * t446 / 0.2e1 + ((-qJD(3) * t67 - qJDD(3) * t99) * t331 + (qJD(3) * t374 + qJDD(3) * t371) * t329 + t521) * t445 / 0.2e1 + (t8 * t422 + t7 * t193 + (t398 * t9 + t479) * t331 + ((t9 * (-t313 - t424) + t352) * t332 + (t398 * t7 + t349) * t330) * t329 - (t303 * t348 + t305 * t347 + t329 * t345) * qJD(4) - g(1) * (t220 + t430) - g(2) * (t218 + t431) - (t322 + t323 + (-pkin(3) + t492) * t443) * g(3) + (t404 * t331 + (-t294 - t434) * t445 + t302 * t407 - t346 - t430 * t314 + t420 * t257) * t39 + (qJD(5) * t400 + t430 * t256 - t431 * t257 + t404 * t446 + t171 - t363) * t35 + (-t420 * t256 + t304 * t407 + t431 * t314 - t423 + t429) * t38) * m(6) + (-g(1) * (t163 + t220) - g(2) * (t161 + t218) - g(3) * (t323 + (-rSges(5,1) * t463 - pkin(3)) * t443 + t417) - t59 * (-t314 * t161 + t256 * t255 + t423) - t60 * (t314 * t163 - t257 * t255 + t346) - t48 * (t161 * t257 - t163 * t256 + t363) - (t59 * (-t115 * t442 + t192 * t303) + t60 * (t117 * t442 - t192 * t305) + t48 * (t115 * t305 - t117 * t303)) * qJD(4) + t29 * t422 + t59 * t429 + t24 * t193 + t48 * t171 + (t30 * t435 + t459 * t60 + t376) * t331 + ((t30 * (-t192 - t313) + t60 * (-t127 - t294) + t24 * t115 + t48 * t93) * t332 + (t59 * t127 + t29 * t192 + t24 * t435 + t459 * t48) * t330) * t329) * m(5) + ((t100 * t200 - t101 * t201 + t181 * t74 - t182 * t75) * t331 + (t279 * t370 + t287 * t367 + t365 * t97 + t366 * t65) * t329 - ((t100 * t216 - t101 * t217) * t331 + (t97 * (t216 * t332 - t217 * t330) + t367 * t311) * t329) * qJD(3) - g(1) * t217 - g(2) * t216 - g(3) * t311) * m(4) - ((t332 * ((t211 * t445 + t304 * t425 + t305 * t427) * t445 + (t210 * t445 + t304 * t426 + t305 * t428) * t446 - (-t304 * t418 - t305 * t419 + t308 * t445) * t331) + t330 * ((t211 * t446 + t302 * t425 + t303 * t427) * t445 + (t210 * t446 + t302 * t426 + t303 * t428) * t446 - (-t302 * t418 - t303 * t419 + t308 * t446) * t331)) * qJD(3) ^ 2 + t512 * t411) * t329 / 0.2e1 - (t303 * t518 + t305 * t517) * qJD(4) / 0.2e1; (t302 * t491 + t490 * t304 - t443 * t506) * t478 + (t302 * t489 + t304 * t488 - t443 * t505) * t477 + (t243 * t481 + t302 * t480 - t354 * t482) * t476 + ((-t334 * t509 - t412 * t506) * t329 + t515 * t304 + t516 * t302 + t490 * t283 + t491 * t281) * t475 + (t245 * t481 + t304 * t480 - t353 * t482) * t474 + ((-t334 * t508 - t412 * t505) * t329 + t513 * t304 + t514 * t302 + t488 * t283 + t489 * t281) * t473 + (t302 * t487 + t304 * t486 - t443 * t504) * t472 + t518 * t281 / 0.2e1 + t517 * t283 / 0.2e1 + t522 * t302 / 0.2e1 + t521 * t304 / 0.2e1 + (t306 * t482 + t307 * t481 + t443 * t480) * t467 + ((-t334 * t507 - t412 * t504) * t329 + t510 * t304 + t511 * t302 + t486 * t283 + t487 * t281) * t466 + t520 * t443 / 0.2e1 + t512 * t335 * t385 + (t347 * t283 + t348 * t281 + (-t39 * t434 - t424 * t9 + t352) * t304 + (-t35 * t460 - t436 * t7 + t349) * t302 + (qJD(3) * t345 + (t39 * t460 + t436 * t9 - t479) * t334) * t329 - g(1) * t432 - g(2) * t433 - g(3) * t421 - (t243 * t39 + t245 * t38 + t307 * t35) * qJD(5) - (-t38 * t433 + t39 * t432) * t314 - (t35 * t433 - t39 * t421) * t257 - (-t35 * t432 + t38 * t421) * t256) * m(6) + (-g(1) * t146 - g(2) * t142 - g(3) * t230 + t24 * (t115 * t304 - t117 * t302) + (t302 * t59 - t304 * t60) * t127 + (t281 * t59 - t283 * t60 + t29 * t302 - t30 * t304) * t192 + ((-t115 * t59 + t117 * t60) * t412 + (t117 * t30 + t60 * t95 - t376) * t334) * t329 - t59 * (-t142 * t314 + t230 * t256) - t60 * (t146 * t314 - t230 * t257) + (t115 * t283 - t117 * t281 - t142 * t257 + t146 * t256 - t302 * t95 + t304 * t93) * t48) * m(5); (t164 * t39 + t166 * t38 - t246 * t35 + (-t38 * t256 + t39 * t257 - g(3) + t7) * t306 - (t35 * t256 - t39 * t314 - g(1) + t8) * t353 - (-t35 * t257 + t38 * t314 - g(2) + t9) * t354) * m(6);];
tau = t1;
