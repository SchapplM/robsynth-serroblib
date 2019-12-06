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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:35:42
% EndTime: 2019-12-05 18:36:09
% DurationCPUTime: 20.47s
% Computational Cost: add. (26150->795), mult. (25223->1054), div. (0->0), fcn. (24420->10), ass. (0->391)
t360 = qJ(1) + qJ(2);
t351 = sin(t360);
t362 = cos(pkin(9));
t363 = sin(qJ(4));
t520 = t362 * t363;
t479 = t351 * t520;
t353 = cos(t360);
t365 = cos(qJ(4));
t526 = t353 * t365;
t283 = t479 + t526;
t519 = t362 * t365;
t532 = t351 * t363;
t286 = t353 * t519 + t532;
t358 = qJD(1) + qJD(2);
t197 = -qJD(4) * t286 + t283 * t358;
t527 = t353 * t363;
t284 = t351 * t519 - t527;
t531 = t351 * t365;
t285 = t353 * t520 - t531;
t198 = -qJD(4) * t285 - t284 * t358;
t361 = sin(pkin(9));
t523 = t358 * t361;
t481 = t351 * t523;
t113 = rSges(5,1) * t198 + rSges(5,2) * t197 - rSges(5,3) * t481;
t432 = t286 * rSges(5,1) - t285 * rSges(5,2);
t529 = t353 * t361;
t186 = -rSges(5,3) * t529 - t432;
t327 = qJDD(4) * t529;
t488 = qJD(4) * t361;
t471 = t351 * t488;
t255 = -t358 * t471 + t327;
t277 = -rSges(5,3) * t362 + (rSges(5,1) * t365 - rSges(5,2) * t363) * t361;
t303 = (-rSges(5,1) * t363 - rSges(5,2) * t365) * t361;
t292 = qJD(4) * t303;
t355 = qJDD(1) + qJDD(2);
t334 = -qJDD(4) * t362 + t355;
t593 = qJD(4) * t362;
t339 = t358 - t593;
t522 = t358 * t362;
t480 = t351 * t522;
t235 = -pkin(3) * t480 - pkin(7) * t481;
t369 = qJD(1) ^ 2;
t366 = cos(qJ(1));
t570 = pkin(1) * t366;
t364 = sin(qJ(1));
t571 = pkin(1) * t364;
t442 = -qJDD(1) * t570 + t369 * t571;
t416 = qJDD(3) * t353 + t442;
t344 = qJD(3) * t351;
t535 = t351 * t358;
t330 = pkin(2) * t535;
t490 = t344 - t330;
t530 = t353 * t358;
t454 = -qJ(3) * t530 - t344 - t490;
t567 = pkin(7) * t361;
t569 = pkin(3) * t362;
t441 = t567 + t569;
t288 = t441 * t353;
t542 = qJ(3) * t351;
t309 = pkin(2) * t353 + t542;
t492 = -t309 - t288;
t374 = t492 * t355 + (-t235 + t454) * t358 + t416;
t469 = t353 * t488;
t49 = -t113 * t339 + t186 * t334 + t255 * t277 + t292 * t469 + t374;
t620 = -g(2) + t49;
t359 = qJ(4) + qJ(5);
t352 = cos(t359);
t357 = qJD(4) + qJD(5);
t307 = -t357 * t362 + t358;
t425 = t307 * t353;
t350 = sin(t359);
t457 = -t357 + t522;
t427 = t350 * t457;
t168 = t351 * t427 + t352 * t425;
t426 = t457 * t352;
t169 = t350 * t425 - t351 * t426;
t100 = rSges(6,1) * t169 + rSges(6,2) * t168 - rSges(6,3) * t481;
t470 = qJD(4) * t531;
t367 = -pkin(8) - pkin(7);
t521 = t361 * t367;
t476 = t358 * t521;
t347 = pkin(4) * t365 + pkin(3);
t568 = pkin(4) * t363;
t451 = t568 * t593;
t496 = -t347 * t480 - t353 * t451;
t121 = t351 * t476 + (t358 * t527 + t470) * pkin(4) - t235 + t496;
t528 = t353 * t362;
t259 = t350 * t528 - t351 * t352;
t260 = t350 * t351 + t352 * t528;
t431 = t260 * rSges(6,1) - t259 * rSges(6,2);
t167 = -rSges(6,3) * t529 - t431;
t335 = pkin(4) * t532;
t494 = -t347 * t528 - t335;
t564 = pkin(7) + t367;
t183 = (t361 * t564 + t569) * t353 + t494;
t453 = t357 * t358;
t202 = t327 + (qJDD(5) * t353 - t351 * t453) * t361;
t430 = -rSges(6,1) * t350 - rSges(6,2) * t352;
t287 = t430 * t361;
t234 = t357 * t287;
t251 = -rSges(6,3) * t362 + (rSges(6,1) * t352 - rSges(6,2) * t350) * t361;
t525 = t357 * t361;
t279 = t353 * t525;
t487 = -qJDD(4) - qJDD(5);
t304 = t362 * t487 + t355;
t565 = pkin(3) - t347;
t464 = t565 * t361;
t395 = t362 * t564 - t464;
t455 = qJD(4) ^ 2 * t361 ^ 2 * t568;
t21 = -t307 * t100 - t339 * t121 + t167 * t304 + t334 * t183 + t202 * t251 + t279 * t234 + t255 * t395 - t353 * t455 + t374;
t619 = t21 - g(2);
t387 = t361 * t395;
t381 = qJD(4) * t387;
t618 = t307 * t167 + t339 * t183 + t251 * t279 + t353 * t381;
t617 = t339 * t186 + t277 * t469;
t533 = t351 * t362;
t257 = t350 * t533 + t352 * t353;
t258 = -t353 * t350 + t352 * t533;
t534 = t351 * t361;
t154 = -Icges(6,5) * t258 + Icges(6,6) * t257 - Icges(6,3) * t534;
t155 = Icges(6,5) * t260 - Icges(6,6) * t259 + Icges(6,3) * t529;
t239 = Icges(6,4) * t260;
t158 = -Icges(6,2) * t259 + Icges(6,6) * t529 + t239;
t238 = Icges(6,4) * t259;
t162 = -Icges(6,1) * t260 - Icges(6,5) * t529 + t238;
t422 = t259 * t158 + t260 * t162;
t545 = Icges(6,4) * t258;
t157 = Icges(6,2) * t257 - Icges(6,6) * t534 - t545;
t237 = Icges(6,4) * t257;
t160 = -Icges(6,1) * t258 - Icges(6,5) * t534 + t237;
t518 = t257 * t157 - t258 * t160;
t616 = t422 + t518 + (-t154 * t351 - t155 * t353) * t361;
t248 = -Icges(6,3) * t362 + (Icges(6,5) * t352 - Icges(6,6) * t350) * t361;
t543 = Icges(6,4) * t352;
t249 = -Icges(6,6) * t362 + (-Icges(6,2) * t350 + t543) * t361;
t544 = Icges(6,4) * t350;
t250 = -Icges(6,5) * t362 + (Icges(6,1) * t352 - t544) * t361;
t103 = t248 * t529 - t249 * t259 + t250 * t260;
t278 = t351 * t525;
t68 = t154 * t529 - t259 * t157 + t260 * t160;
t613 = -t103 * t307 + t278 * t68;
t75 = -t155 * t362 - (t158 * t350 + t162 * t352) * t361;
t434 = rSges(3,1) * t351 + rSges(3,2) * t353;
t275 = t434 * t358;
t560 = pkin(1) * qJD(1);
t484 = t364 * t560;
t252 = t275 + t484;
t611 = -t431 + t494 - t542;
t610 = -t432 - t542;
t265 = Icges(5,4) * t286;
t178 = -Icges(5,2) * t285 + Icges(5,6) * t529 + t265;
t264 = Icges(5,4) * t285;
t182 = -Icges(5,1) * t286 - Icges(5,5) * t529 + t264;
t517 = t283 * t178 + t182 * t284;
t421 = -t178 * t285 - t286 * t182;
t67 = -t155 * t534 + t257 * t158 + t162 * t258;
t175 = Icges(5,5) * t286 - Icges(5,6) * t285 + Icges(5,3) * t529;
t82 = -t175 * t362 + (-t178 * t363 - t182 * t365) * t361;
t316 = rSges(4,1) * t480;
t333 = rSges(4,2) * t529;
t413 = -rSges(4,1) * t528 - rSges(4,3) * t351;
t236 = -t333 - t413;
t493 = -t309 - t236;
t561 = rSges(4,3) * t353;
t86 = t493 * t355 + (t316 + (-rSges(4,2) * t534 - t561) * t358 + t454) * t358 + t416;
t608 = -g(2) + t86;
t478 = t353 * t523;
t424 = t307 * t351;
t170 = -t352 * t424 + t353 * t427;
t171 = -t350 * t424 - t353 * t426;
t514 = t171 * rSges(6,1) + t170 * rSges(6,2);
t101 = -rSges(6,3) * t478 + t514;
t338 = pkin(4) * t526;
t474 = qJD(4) * t338 + t351 * t451 + t353 * t476;
t410 = t362 * t565 + t567;
t587 = t353 * t410;
t122 = (-t335 + t587) * t358 + t474;
t504 = -t258 * rSges(6,1) + t257 * rSges(6,2);
t165 = -rSges(6,3) * t534 + t504;
t201 = (t351 * t487 - t353 * t453) * t361;
t254 = (-qJD(4) * t530 - qJDD(4) * t351) * t361;
t345 = qJD(3) * t353;
t400 = (-qJDD(1) * t364 - t366 * t369) * pkin(1);
t541 = qJ(3) * t353;
t429 = -pkin(2) * t351 + t541;
t299 = t358 * t309;
t489 = t345 - t299;
t379 = qJDD(3) * t351 + t355 * t429 + t400 + (t345 + t489) * t358;
t411 = t351 * t441;
t376 = -t288 * t358 ^ 2 - t355 * t411 + t379;
t536 = t347 * t362;
t397 = t441 - t536;
t491 = pkin(4) * t527 + t351 * t521;
t20 = t334 * t491 + t339 * t122 - t254 * t395 + t304 * t165 + t307 * t101 - t201 * t251 + t278 * t234 + (t334 * t397 - t455) * t351 + t376;
t607 = -g(3) + t20;
t199 = qJD(4) * t284 + t285 * t358;
t200 = qJD(4) * t283 - t286 * t358;
t505 = t200 * rSges(5,1) + t199 * rSges(5,2);
t114 = -rSges(5,3) * t478 + t505;
t500 = -t284 * rSges(5,1) + t283 * rSges(5,2);
t184 = -rSges(5,3) * t534 + t500;
t48 = t114 * t339 + t184 * t334 - t254 * t277 + t292 * t471 + t376;
t606 = -g(3) + t48;
t317 = rSges(4,2) * t478;
t562 = rSges(4,2) * t361;
t563 = rSges(4,1) * t362;
t433 = t562 - t563;
t385 = t351 * t433 + t561;
t85 = t355 * t385 + (t358 * t413 + t317) * t358 + t379;
t605 = -g(3) + t85;
t310 = rSges(3,1) * t353 - t351 * rSges(3,2);
t604 = t275 * t358 - t310 * t355 - g(2) + t442;
t276 = -rSges(3,1) * t530 + rSges(3,2) * t535;
t603 = t276 * t358 - t355 * t434 - g(3) + t400;
t599 = -t358 * t236 - t317;
t262 = t358 * t288;
t472 = -t262 + t489;
t598 = -t345 + t472 - t474 - t514 + t618;
t194 = -t259 * rSges(6,1) - t260 * rSges(6,2);
t597 = t362 * t100 + t194 * t307 + t234 * t529 - t279 * t287;
t193 = t257 * rSges(6,1) + t258 * rSges(6,2);
t223 = t251 * t529;
t596 = t362 * t101 + t307 * t193 - t358 * t223 - t234 * t534 + t278 * t287;
t595 = t165 * t481 + t193 * t279 + t278 * t194;
t594 = qJ(3) * t535 - t345 - t505 + t617;
t414 = -rSges(6,3) * t361 - pkin(2) - t536;
t461 = -qJ(3) - t568;
t261 = t358 * t411;
t297 = t358 * t429;
t495 = -t297 - t344;
t475 = t261 + t495;
t390 = -t307 * t165 - t251 * t278 - t339 * (t351 * t410 + t491) - t351 * t381 + t475;
t61 = t390 + t484;
t483 = t366 * t560;
t443 = -t345 + t483;
t391 = t358 * t492 - t443;
t62 = t391 + t618;
t592 = ((-t461 * t61 - t521 * t62) * t351 + (-t414 * t61 + t461 * t62) * t353) * t358 + t619 * t353 * (-pkin(2) + (-rSges(6,3) + t367) * t361);
t399 = -t569 - pkin(2) + (-rSges(5,3) - pkin(7)) * t361;
t566 = g(2) * t353;
t444 = -t344 + t484;
t418 = -t297 + t444;
t585 = -t339 * t184 - t277 * t471;
t87 = t261 + t418 + t585;
t88 = t391 + t617;
t591 = (t49 * t399 + (-t88 * qJ(3) - t87 * (-rSges(5,3) * t361 - pkin(2) - t441)) * t358) * t353 - t399 * t566;
t589 = (-t254 * t183 - t255 * t491 - t122 * t469 - t201 * t167 - t278 * t100 - t202 * t165 - t279 * t101 + (-t121 * t488 - t255 * t397) * t351) * t361;
t102 = -t248 * t534 + t249 * t257 - t250 * t258;
t588 = t102 * t307 + t279 * t67;
t174 = -Icges(5,5) * t284 + Icges(5,6) * t283 - Icges(5,3) * t534;
t548 = Icges(5,4) * t284;
t177 = Icges(5,2) * t283 - Icges(5,6) * t534 - t548;
t263 = Icges(5,4) * t283;
t180 = -Icges(5,1) * t284 - Icges(5,5) * t534 + t263;
t72 = t174 * t529 - t285 * t177 + t286 * t180;
t392 = t351 * (Icges(5,2) * t284 + t180 + t263) - t353 * (-Icges(5,2) * t286 - t182 - t264);
t393 = t351 * (-Icges(5,1) * t283 + t177 - t548) - t353 * (Icges(5,1) * t285 + t178 + t265);
t281 = (-Icges(6,2) * t352 - t544) * t361;
t377 = t278 * (Icges(6,2) * t258 + t160 + t237) - t279 * (-Icges(6,2) * t260 - t162 - t238) - t307 * (t250 + t281);
t282 = (-Icges(6,1) * t350 - t543) * t361;
t378 = t278 * (-Icges(6,1) * t257 + t157 - t545) - t279 * (Icges(6,1) * t259 + t158 + t239) - t307 * (t249 - t282);
t582 = t362 ^ 2;
t581 = t201 / 0.2e1;
t580 = t202 / 0.2e1;
t579 = t254 / 0.2e1;
t578 = t255 / 0.2e1;
t577 = t278 / 0.2e1;
t576 = -t278 / 0.2e1;
t575 = -t279 / 0.2e1;
t574 = t279 / 0.2e1;
t572 = -t362 / 0.2e1;
t93 = Icges(6,5) * t171 + Icges(6,6) * t170 - Icges(6,3) * t478;
t95 = Icges(6,4) * t171 + Icges(6,2) * t170 - Icges(6,6) * t478;
t97 = Icges(6,1) * t171 + Icges(6,4) * t170 - Icges(6,5) * t478;
t36 = -t362 * t93 + ((-t157 * t357 + t97) * t352 + (-t160 * t357 - t95) * t350) * t361;
t558 = t36 * t278;
t92 = Icges(6,5) * t169 + Icges(6,6) * t168 - Icges(6,3) * t481;
t94 = Icges(6,4) * t169 + Icges(6,2) * t168 - Icges(6,6) * t481;
t96 = Icges(6,1) * t169 + Icges(6,4) * t168 - Icges(6,5) * t481;
t37 = -t362 * t92 + ((-t158 * t357 + t96) * t352 + (t162 * t357 - t94) * t350) * t361;
t557 = t37 * t279;
t74 = -t154 * t362 + (-t157 * t350 + t160 * t352) * t361;
t556 = t74 * t201;
t555 = t75 * t202;
t81 = -t174 * t362 + (-t177 * t363 + t180 * t365) * t361;
t554 = t81 * t254;
t553 = t82 * t255;
t552 = rSges(4,3) + qJ(3);
t120 = -t248 * t362 + (-t249 * t350 + t250 * t352) * t361;
t280 = (-Icges(6,5) * t350 - Icges(6,6) * t352) * t361;
t231 = t357 * t280;
t232 = t357 * t281;
t233 = t357 * t282;
t83 = -t231 * t362 + ((-t249 * t357 + t233) * t352 + (-t250 * t357 - t232) * t350) * t361;
t551 = t120 * t304 + t83 * t307;
t546 = Icges(5,4) * t365;
t273 = -Icges(5,6) * t362 + (-Icges(5,2) * t363 + t546) * t361;
t547 = Icges(5,4) * t363;
t274 = -Icges(5,5) * t362 + (Icges(5,1) * t365 - t547) * t361;
t300 = (-Icges(5,5) * t363 - Icges(5,6) * t365) * t361;
t289 = qJD(4) * t300;
t301 = (-Icges(5,2) * t365 - t547) * t361;
t290 = qJD(4) * t301;
t302 = (-Icges(5,1) * t363 - t546) * t361;
t291 = qJD(4) * t302;
t104 = -t289 * t362 + (-t290 * t363 + t291 * t365 + (-t273 * t365 - t274 * t363) * qJD(4)) * t361;
t272 = -Icges(5,3) * t362 + (Icges(5,5) * t365 - Icges(5,6) * t363) * t361;
t132 = -t272 * t362 + (-t273 * t363 + t274 * t365) * t361;
t550 = t104 * t339 + t132 * t334;
t116 = t272 * t529 - t273 * t285 + t274 * t286;
t539 = t116 * t339;
t516 = -t283 * t177 + t284 * t180;
t515 = -t167 * t362 + t223;
t498 = t273 - t302;
t497 = t274 + t301;
t266 = pkin(4) * t479 + t338;
t473 = -t165 - t491;
t468 = -t534 / 0.2e1;
t467 = t529 / 0.2e1;
t466 = -t523 / 0.2e1;
t465 = -pkin(2) - t563;
t463 = -t488 / 0.2e1;
t462 = t488 / 0.2e1;
t450 = t351 * t466;
t449 = t353 * t466;
t448 = t351 * t463;
t447 = t351 * t462;
t446 = t353 * t463;
t445 = t353 * t462;
t332 = rSges(2,1) * t366 - t364 * rSges(2,2);
t435 = rSges(2,1) * t364 + rSges(2,2) * t366;
t428 = -t87 * t351 + t88 * t353;
t420 = -t184 * t353 + t186 * t351;
t419 = (Icges(5,5) * t283 + Icges(5,6) * t284) * t351 - (-Icges(5,5) * t285 - Icges(5,6) * t286) * t353;
t267 = t285 * pkin(4);
t417 = t299 + t443;
t415 = t174 * t534 + t516;
t253 = -t310 * t358 - t483;
t71 = -t175 * t534 + t517;
t404 = (t351 * t415 + t353 * t71) * t361;
t73 = t175 * t529 + t421;
t403 = (-t351 * t72 + t353 * t73) * t361;
t394 = -(Icges(6,5) * t257 + Icges(6,6) * t258) * t278 + (-Icges(6,5) * t259 - Icges(6,6) * t260) * t279 + t280 * t307;
t105 = t420 * t488;
t386 = t183 - t587;
t383 = -t113 - t490 - t235;
t382 = t353 * t387;
t215 = -t351 * t552 + t353 * t465 + t333;
t16 = t157 * t168 + t160 * t169 - t259 * t95 + t260 * t97 + (-t154 * t535 + t353 * t93) * t361;
t17 = t158 * t168 - t162 * t169 - t259 * t94 + t260 * t96 + (-t155 * t535 + t353 * t92) * t361;
t18 = t157 * t170 + t160 * t171 + t257 * t95 - t258 * t97 + (-t154 * t530 - t351 * t93) * t361;
t19 = t158 * t170 - t162 * t171 + t257 * t94 - t258 * t96 + (-t155 * t530 - t351 * t92) * t361;
t66 = -t154 * t534 + t518;
t28 = -t278 * t66 + t588;
t69 = t155 * t529 - t422;
t29 = t279 * t69 - t613;
t53 = t168 * t249 + t169 * t250 - t232 * t259 + t233 * t260 + (t231 * t353 - t248 * t535) * t361;
t54 = t170 * t249 + t171 * t250 + t232 * t257 - t233 * t258 + (-t231 * t351 - t248 * t530) * t361;
t380 = (t102 * t304 - t18 * t278 + t19 * t279 + t201 * t66 + t202 * t67 + t307 * t54) * t468 + (-t257 * t377 - t258 * t378 - t394 * t534) * t577 + (t259 * t377 + t260 * t378 + t394 * t529) * t575 - (-t394 * t362 + (t350 * t377 + t352 * t378) * t361) * t307 / 0.2e1 + (t103 * t304 - t16 * t278 + t17 * t279 + t201 * t68 + t202 * t69 + t307 * t53) * t467 + t29 * t450 + t28 * t449 + (-t102 * t362 + (-t351 * t66 + t353 * t67) * t361) * t581 + (-t103 * t362 + (-t351 * t68 + t353 * t69) * t361) * t580 + (-t362 * t54 + ((-t358 * t66 + t19) * t353 + (-t358 * t67 - t18) * t351) * t361) * t576 + (-t362 * t53 + ((-t358 * t68 + t17) * t353 + (-t358 * t69 - t16) * t351) * t361) * t574 + t304 * (-t120 * t362 + (-t351 * t74 + t353 * t75) * t361) / 0.2e1 + (t551 + t555 + t556 + t557 - t558) * t572 + t307 * (-t362 * t83 + ((-t358 * t74 + t37) * t353 + (-t358 * t75 - t36) * t351) * t361) / 0.2e1;
t214 = t552 * t353 + (-pkin(2) + t433) * t351;
t131 = t351 * t399 + t500 + t541;
t375 = -pkin(4) * t470 - t100 - t490 - t496;
t119 = t351 * t414 + t491 + t504 + t541;
t228 = t358 * t385;
t163 = -t228 + t418;
t164 = t358 * t493 - t443;
t373 = ((t163 * t552 - t164 * t562) * t351 + (-t163 * t465 - t164 * t552) * t353) * t358;
t115 = -t272 * t534 + t273 * t283 - t274 * t284;
t106 = t115 * t339;
t41 = qJD(4) * t404 + t106;
t42 = qJD(4) * t403 + t539;
t108 = Icges(5,5) * t200 + Icges(5,6) * t199 - Icges(5,3) * t478;
t110 = Icges(5,4) * t200 + Icges(5,2) * t199 - Icges(5,6) * t478;
t112 = Icges(5,1) * t200 + Icges(5,4) * t199 - Icges(5,5) * t478;
t46 = -t108 * t362 + (-t110 * t363 + t112 * t365 + (-t177 * t365 - t180 * t363) * qJD(4)) * t361;
t107 = Icges(5,5) * t198 + Icges(5,6) * t197 - Icges(5,3) * t481;
t109 = Icges(5,4) * t198 + Icges(5,2) * t197 - Icges(5,6) * t481;
t111 = Icges(5,1) * t198 + Icges(5,4) * t197 - Icges(5,5) * t481;
t47 = -t107 * t362 + (-t109 * t363 + t111 * t365 + (-t178 * t365 + t182 * t363) * qJD(4)) * t361;
t59 = t197 * t273 + t198 * t274 - t285 * t290 + t286 * t291 + (-t272 * t535 + t289 * t353) * t361;
t60 = t199 * t273 + t200 * t274 + t283 * t290 - t284 * t291 + (-t272 * t530 - t289 * t351) * t361;
t371 = (t106 + (t517 * t353 + (t415 + t421 - t73) * t351) * t488) * t446 - t558 / 0.2e1 + t557 / 0.2e1 + t54 * t576 + t116 * t578 + t115 * t579 + t103 * t580 + t102 * t581 + t555 / 0.2e1 + t556 / 0.2e1 + t553 / 0.2e1 + t554 / 0.2e1 + t550 + t551 + (-(t69 + t616) * t278 + t588) * t575 + (t29 + (-t66 + t616) * t279 + t613) * t577 + (t53 + t28) * t574 + (t42 - t539 + (-(-t517 - t72) * t351 + t415 * t353 + (-t421 - t516) * t353 - t71 * t351 + (-t175 * t351 ^ 2 + (-t174 * t351 - t175 * t353) * t353) * t361) * t488) * t447 + (Icges(4,2) * t582 + (Icges(4,1) * t361 + 0.2e1 * Icges(4,4) * t362) * t361 + Icges(3,3)) * t355 + (t46 + t60) * t448 + (t47 + t59 + t41) * t445;
t222 = t251 * t534;
t210 = -rSges(5,1) * t285 - rSges(5,2) * t286;
t209 = rSges(5,1) * t283 + rSges(5,2) * t284;
t57 = -t279 * t165 + t278 * t167 + (t386 * t534 - t491 * t529) * qJD(4);
t33 = t109 * t283 - t111 * t284 + t178 * t199 - t182 * t200 + (-t107 * t351 - t175 * t530) * t361;
t32 = t110 * t283 - t112 * t284 + t177 * t199 + t180 * t200 + (-t108 * t351 - t174 * t530) * t361;
t31 = -t109 * t285 + t111 * t286 + t178 * t197 - t182 * t198 + (t107 * t353 - t175 * t535) * t361;
t30 = -t110 * t285 + t112 * t286 + t177 * t197 + t180 * t198 + (t108 * t353 - t174 * t535) * t361;
t1 = [Icges(2,3) * qJDD(1) + t371 + (t604 * (-t310 - t570) + t603 * (-t434 - t571) + (t253 - t276 + t483) * t252) * m(3) + ((qJDD(1) * t435 + g(3)) * t435 + (qJDD(1) * t332 + g(2)) * t332) * m(2) + (t62 * (t375 + t484) + t607 * (t119 - t571) + (-t62 + t598) * t61 + t592 + t619 * (t611 - t570)) * m(6) + (t88 * (t383 + t484) + (t483 - t262 - t417 - t88 + t594) * t87 + t606 * (t131 - t571) + t591 + t620 * (t610 - t570)) * m(5) + (t164 * (t316 + t330 + t444) + t373 + t608 * (t215 - t570) + t605 * (t214 - t571) + (-t164 - t417 + t443 + t599) * t163) * m(4); t371 + ((t375 - t390) * t62 + t598 * t61 + t607 * t119 + t592 + t619 * t611) * m(6) + ((t383 - t475 - t585) * t88 + (t472 + t594) * t87 + t606 * t131 + t591 + t620 * t610) * m(5) + (t373 + t608 * t215 + t605 * t214 + (t316 - t490 + t228 - t495) * t164 + (-t345 + t489 + t599) * t163) * m(4) + (-t252 * t276 + t253 * t275 + (-t252 * t358 - t604) * t310 - (t253 * t358 + t603) * t434) * m(3); (-m(4) - m(5) - m(6)) * (g(3) * t351 + t566) + m(4) * (t351 * t85 + t353 * t86) + m(5) * (t48 * t351 + t49 * t353) + m(6) * (t20 * t351 + t21 * t353); ((t283 * t497 + t284 * t498 - t300 * t534) * t339 + (-t283 * t392 - t284 * t393 + t419 * t534) * t488) * t447 + (-t115 * t362 + t404) * t579 - t339 * (-t362 * t300 * t339 + ((-t363 * t497 - t365 * t498) * t339 + ((t363 * t392 + t365 * t393) * t361 + t419 * t362) * qJD(4)) * t361) / 0.2e1 + (-t362 * t60 + ((t358 * t415 + t33) * t353 + (-t358 * t71 - t32) * t351) * t361) * t448 + t380 + t339 * (-t104 * t362 + ((-t358 * t81 + t47) * t353 + (-t358 * t82 - t46) * t351) * t361) / 0.2e1 + t41 * t449 + t42 * t450 + (-t116 * t362 + t403) * t578 + (t554 + t553 + (-t351 * t46 + t353 * t47) * t488 + t550) * t572 + ((-t285 * t497 - t286 * t498 + t300 * t529) * t339 + (t285 * t392 + t286 * t393 - t419 * t529) * t488) * t446 + t334 * (-t132 * t362 + (-t351 * t81 + t353 * t82) * t361) / 0.2e1 + (t116 * t334 + t254 * t72 + t255 * t73 + t339 * t59 + (-t30 * t351 + t31 * t353) * t488) * t467 + (-t362 * t59 + ((-t358 * t72 + t31) * t353 + (-t358 * t73 - t30) * t351) * t361) * t445 + (t115 * t334 - t254 * t415 + t255 * t71 + t339 * t60 + (-t32 * t351 + t33 * t353) * t488) * t468 + (-g(2) * (t266 + t193) - g(3) * (-t267 + t194) - g(1) * (t430 - t568) * t361 + (t473 * t353 + (t167 + t386) * t351) * t589 + t21 * (-t362 * t183 + t382 + t515) + t20 * (t222 + t473 * t362 + (-t582 * t565 + (t362 * t367 - t464) * t361) * t351) + (((-t101 - t122 + (t183 + t167) * t358) * t353 + (-t121 - t100 + (t351 * t397 + t491) * t358) * t351) * t361 - (-t266 * t353 + t267 * t351) * t488 + t595) * t57 + (t362 * t121 + (-t251 - t395) * t481 - t267 * t339 + t597) * t62 + (t362 * t122 + t266 * t339 - t358 * t382 + t596) * t61) * m(6) + (-g(1) * t303 - g(2) * t209 - g(3) * t210 + (t113 * t88 + t114 * t87 - t184 * t48 - t186 * t49) * t362 + ((-t184 * t255 - t186 * t254) * t420 + t428 * t292 + ((-t358 * t87 + t49) * t353 + (-t358 * t88 + t48) * t351) * t277 + (-0.2e1 * t113 * t351 - 0.2e1 * t114 * t353 + t184 * t535 + t186 * t530) * t105) * t361 - (-t209 * t87 - t210 * t88) * t339 - (t105 * (-t209 * t353 - t210 * t351) + t428 * t303) * t488) * m(5); t380 + (-g(1) * t287 - g(2) * t193 - g(3) * t194 + t21 * t515 + t20 * (-t362 * t165 + t222) + (-t353 * t165 + t167 * t351) * t589 + (-t251 * t361 * t535 + t597) * t62 + t596 * t61 + ((-t351 * t100 - t353 * t101 + t167 * t530) * t361 + t595) * t57) * m(6);];
tau = t1;
