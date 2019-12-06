% Calculate vector of inverse dynamics joint torques for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:37
% EndTime: 2019-12-05 17:36:23
% DurationCPUTime: 27.03s
% Computational Cost: add. (12208->573), mult. (16242->743), div. (0->0), fcn. (15498->8), ass. (0->276)
t559 = Icges(5,4) + Icges(6,4);
t531 = -Icges(5,1) - Icges(6,1);
t530 = Icges(5,5) + Icges(6,5);
t529 = Icges(5,2) + Icges(6,2);
t528 = Icges(5,6) + Icges(6,6);
t275 = cos(qJ(4));
t565 = t559 * t275;
t273 = sin(qJ(4));
t564 = t559 * t273;
t269 = qJ(1) + pkin(7);
t265 = sin(t269);
t266 = cos(t269);
t271 = cos(pkin(8));
t421 = t271 * t273;
t203 = -t265 * t275 + t266 * t421;
t420 = t271 * t275;
t428 = t265 * t273;
t204 = t266 * t420 + t428;
t172 = Icges(6,4) * t204;
t270 = sin(pkin(8));
t426 = t266 * t270;
t110 = -Icges(6,2) * t203 + Icges(6,6) * t426 + t172;
t175 = Icges(5,4) * t204;
t113 = -Icges(5,2) * t203 + Icges(5,6) * t426 + t175;
t498 = t113 + t110;
t171 = Icges(6,4) * t203;
t117 = -Icges(6,1) * t204 - Icges(6,5) * t426 + t171;
t174 = Icges(5,4) * t203;
t120 = -Icges(5,1) * t204 - Icges(5,5) * t426 + t174;
t548 = t120 + t117;
t553 = t203 * t498 + t204 * t548;
t201 = t265 * t421 + t266 * t275;
t424 = t266 * t273;
t202 = t265 * t420 - t424;
t551 = t201 * t498 + t202 * t548;
t558 = -Icges(5,3) - Icges(6,3);
t429 = t265 * t270;
t436 = Icges(6,4) * t202;
t109 = Icges(6,2) * t201 - Icges(6,6) * t429 - t436;
t439 = Icges(5,4) * t202;
t112 = Icges(5,2) * t201 - Icges(5,6) * t429 - t439;
t550 = t112 + t109;
t170 = Icges(6,4) * t201;
t115 = -Icges(6,1) * t202 - Icges(6,5) * t429 + t170;
t173 = Icges(5,4) * t201;
t118 = -Icges(5,1) * t202 - Icges(5,5) * t429 + t173;
t549 = t118 + t115;
t103 = -Icges(6,5) * t202 + Icges(6,6) * t201 - Icges(6,3) * t429;
t106 = -Icges(5,5) * t202 + Icges(5,6) * t201 - Icges(5,3) * t429;
t527 = t103 + t106;
t104 = Icges(6,5) * t204 - Icges(6,6) * t203 + Icges(6,3) * t426;
t107 = Icges(5,5) * t204 - Icges(5,6) * t203 + Icges(5,3) * t426;
t557 = t107 + t104;
t556 = (-t530 * t273 - t528 * t275) * t270;
t555 = (t529 * t275 + t564) * t270;
t554 = (t531 * t273 - t565) * t270;
t524 = t558 * t271 + (-t273 * t528 + t275 * t530) * t270;
t523 = t528 * t271 + (t273 * t529 - t565) * t270;
t522 = t530 * t271 + (t275 * t531 + t564) * t270;
t552 = t201 * t550 - t202 * t549;
t540 = -t557 * t429 + t551;
t505 = -t550 * t203 + t549 * t204 + t527 * t426;
t504 = t557 * t426 - t553;
t131 = qJD(1) * t201 - qJD(4) * t204;
t288 = t203 * qJD(4);
t132 = -qJD(1) * t202 - t288;
t388 = qJD(1) * t270;
t373 = t265 * t388;
t538 = t528 * t131 + t530 * t132 + t558 * t373;
t133 = qJD(1) * t203 + qJD(4) * t202;
t134 = -qJD(1) * t204 + qJD(4) * t201;
t372 = t266 * t388;
t537 = t528 * t133 + t530 * t134 + t558 * t372;
t536 = t529 * t131 + t559 * t132 - t528 * t373;
t535 = t529 * t133 + t559 * t134 - t528 * t372;
t534 = -t559 * t131 + t531 * t132 + t530 * t373;
t533 = -t559 * t133 + t531 * t134 + t530 * t372;
t521 = t556 * qJD(4);
t520 = t555 * qJD(4);
t519 = t554 * qJD(4);
t516 = t527 * t429 - t552;
t247 = pkin(4) * t428;
t263 = pkin(4) * t275 + pkin(3);
t425 = t266 * t271;
t547 = -t204 * rSges(6,1) + t203 * rSges(6,2) - t263 * t425 - t247;
t257 = -qJD(4) * t271 + qJD(1);
t385 = qJD(5) * t270;
t369 = t265 * t385;
t387 = qJD(4) * t270;
t370 = t266 * t387;
t272 = -qJ(5) - pkin(6);
t463 = pkin(6) + t272;
t464 = pkin(3) - t263;
t405 = (t463 - rSges(6,3)) * t271 + (rSges(6,1) * t275 - rSges(6,2) * t273 - t464) * t270;
t468 = pkin(3) * t271;
t410 = (t270 * t463 + t468) * t266 - rSges(6,3) * t426 + t547;
t546 = t257 * t410 + t405 * t370 - t369;
t501 = -t523 * t201 + t522 * t202 - t524 * t429;
t500 = t523 * t203 - t522 * t204 + t524 * t426;
t333 = t204 * rSges(5,1) - t203 * rSges(5,2);
t126 = -rSges(5,3) * t426 - t333;
t200 = -rSges(5,3) * t271 + (rSges(5,1) * t275 - rSges(5,2) * t273) * t270;
t544 = t126 * t257 + t200 * t370;
t382 = qJD(1) * qJD(4);
t169 = (qJDD(4) * t266 - t265 * t382) * t270;
t243 = t266 * t385;
t256 = -qJDD(4) * t271 + qJDD(1);
t466 = pkin(6) * t270;
t344 = t466 + t468;
t205 = t344 * t266;
t433 = qJ(3) * t265;
t230 = pkin(2) * t266 + t433;
t276 = cos(qJ(1));
t469 = pkin(1) * t276;
t360 = -t230 - t469;
t345 = -t205 + t360;
t277 = qJD(1) ^ 2;
t274 = sin(qJ(1));
t470 = pkin(1) * t274;
t264 = t277 * t470;
t392 = qJDD(3) * t266 + t264;
t282 = qJDD(1) * t345 + t392;
t389 = qJD(1) * t266;
t260 = qJD(3) * t265;
t390 = qJD(1) * t265;
t391 = -pkin(2) * t390 + t260;
t164 = qJ(3) * t389 + t391;
t165 = -pkin(6) * t373 - t390 * t468;
t338 = -t164 - t165 - t260;
t452 = rSges(6,2) * t275;
t301 = t270 * (-rSges(6,1) * t273 - t452);
t449 = pkin(4) * qJD(4);
t377 = t273 * t449;
t384 = qJD(5) * t271;
t396 = qJD(4) * t301 - t270 * t377 - t384;
t356 = t396 * qJD(4);
t422 = t270 * t272;
t430 = t263 * t271;
t314 = t422 - t430;
t483 = -rSges(6,1) * t132 - rSges(6,2) * t131 - t243;
t460 = rSges(6,3) * t373 - t314 * t390 - (t273 * t389 - t288) * pkin(4) + t165 + t483;
t11 = t460 * t257 + t410 * t256 + t405 * t169 + (-qJDD(5) * t265 + t266 * t356) * t270 + (t338 - t243) * qJD(1) + t282;
t543 = t11 - g(2);
t542 = t537 * t271 + (t533 * t275 + t535 * t273 + (t549 * t273 + t550 * t275) * qJD(4)) * t270;
t541 = -t538 * t271 + (-t534 * t275 - t536 * t273 + (t548 * t273 - t498 * t275) * qJD(4)) * t270;
t539 = -t521 * t271 + (t519 * t275 + t520 * t273 + (t522 * t273 + t523 * t275) * qJD(4)) * t270;
t532 = -t524 * t271 + (t523 * t273 - t522 * t275) * t270;
t518 = (-t505 * t265 + t504 * t266) * t270;
t517 = (t516 * t265 + t540 * t266) * t270;
t515 = rSges(6,1) + pkin(4);
t514 = t501 * t257;
t432 = qJ(3) * t266;
t329 = t432 - t470;
t511 = t500 * t257;
t168 = (-qJDD(4) * t265 - t266 * t382) * t270;
t261 = qJD(3) * t266;
t328 = -pkin(2) * t265 + t432;
t395 = qJD(1) * t230 - t261;
t419 = t276 * t277;
t280 = qJDD(1) * t328 + qJDD(3) * t265 + (-qJDD(1) * t274 - t419) * pkin(1) + (-t395 + t261) * qJD(1);
t307 = t265 * t344;
t279 = -qJDD(1) * t307 - t205 * t277 + t280;
t353 = -t202 * rSges(6,1) + t201 * rSges(6,2) + pkin(4) * t424 + t265 * t422;
t364 = t464 * t271;
t306 = t364 + t466;
t354 = t271 * t377;
t376 = t275 * t449;
t481 = t134 * rSges(6,1) + t133 * rSges(6,2) + t265 * t354 + t266 * t376 + t272 * t372 - t369;
t459 = -rSges(6,3) * t372 + (t266 * t306 - t247) * qJD(1) + t481;
t471 = rSges(6,3) - pkin(6);
t10 = qJDD(5) * t426 + t459 * t257 + t353 * t256 - t405 * t168 + (t256 * t364 + (-qJD(5) * qJD(1) - t256 * t471 + t356) * t270) * t265 + t279;
t510 = -g(3) + t10;
t509 = t517 * qJD(4) + t514;
t508 = t518 * qJD(4) + t511;
t507 = (t521 * t266 - t524 * t390) * t270 + t519 * t204 + t520 * t203 - t522 * t132 - t523 * t131;
t506 = (t521 * t265 + t524 * t389) * t270 + t519 * t202 + t520 * t201 + t522 * t134 + t523 * t133;
t39 = -t103 * t271 + (-t109 * t273 + t115 * t275) * t270;
t41 = -t106 * t271 + (-t112 * t273 + t118 * t275) * t270;
t503 = -t39 - t41;
t40 = -t104 * t271 + (-t110 * t273 - t117 * t275) * t270;
t42 = -t107 * t271 + (-t113 * t273 - t120 * t275) * t270;
t502 = t40 + t42;
t495 = (t529 * t204 + t171 + t174 + t548) * t266 + (t529 * t202 + t170 + t173 + t549) * t265;
t494 = (t531 * t203 - t172 - t175 - t498) * t266 + (t531 * t201 - t436 - t439 + t550) * t265;
t492 = (t530 * t203 + t528 * t204) * t266 + (t530 * t201 + t528 * t202) * t265;
t491 = -t540 * t265 + t516 * t266;
t490 = t542 * t265 + t541 * t266;
t489 = (-t533 * t202 - t535 * t201 - t549 * t134 - t550 * t133 + (t537 * t265 + t527 * t389) * t270) * t265 + ((-t538 * t265 - t389 * t557) * t270 + t534 * t202 + t536 * t201 - t548 * t134 + t498 * t133) * t266;
t488 = (-t534 * t204 - t536 * t203 - t548 * t132 + t498 * t131 + (t538 * t266 - t390 * t557) * t270) * t266 + ((-t537 * t266 + t527 * t390) * t270 + t533 * t204 + t535 * t203 - t549 * t132 - t550 * t131) * t265;
t487 = -t522 - t555;
t486 = -t523 - t554;
t485 = t532 * t256 + t539 * t257;
t484 = t270 * t471;
t454 = rSges(3,1) * t266;
t339 = -t454 - t469;
t213 = t265 * rSges(3,2) + t339;
t478 = t271 ^ 2;
t477 = t168 / 0.2e1;
t476 = t169 / 0.2e1;
t467 = pkin(4) * t273;
t465 = g(2) * t266;
t453 = rSges(4,1) * t271;
t451 = rSges(6,3) * t270;
t450 = pkin(1) * qJD(1);
t444 = t39 * t168;
t443 = t40 * t169;
t442 = t41 * t168;
t441 = t42 * t169;
t440 = rSges(4,3) + qJ(3);
t408 = t134 * rSges(5,1) + t133 * rSges(5,2);
t407 = t202 * rSges(6,2) + t515 * t201;
t406 = t204 * rSges(6,2) + t515 * t203;
t404 = -t202 * rSges(5,1) + t201 * rSges(5,2);
t383 = -m(4) - m(5) - m(6);
t381 = rSges(5,3) * t429;
t378 = t274 * t450;
t371 = t265 * t387;
t366 = -pkin(2) - t453;
t365 = -pkin(2) - t451;
t362 = -t387 / 0.2e1;
t361 = t387 / 0.2e1;
t358 = -qJ(3) - t467;
t357 = t405 * t266;
t350 = t265 * t362;
t349 = t265 * t361;
t348 = t266 * t362;
t347 = t266 * t361;
t244 = rSges(4,2) * t426;
t308 = -rSges(4,1) * t425 - rSges(4,3) * t265;
t160 = -t244 - t308;
t346 = -t160 + t360;
t246 = rSges(2,1) * t276 - t274 * rSges(2,2);
t337 = rSges(2,1) * t274 + rSges(2,2) * t276;
t336 = rSges(3,1) * t265 + rSges(3,2) * t266;
t335 = -rSges(4,2) * t270 + t453;
t334 = rSges(5,1) * t132 + rSges(5,2) * t131;
t330 = -t433 - t469;
t122 = -t381 + t404;
t227 = (-rSges(5,1) * t273 - rSges(5,2) * t275) * t270;
t215 = qJD(4) * t227;
t75 = -rSges(5,3) * t372 + t408;
t21 = t122 * t256 - t168 * t200 + t215 * t371 + t257 * t75 + t279;
t73 = -rSges(5,3) * t373 + t334;
t22 = qJD(1) * t338 + t126 * t256 + t169 * t200 + t215 * t370 - t257 * t73 + t282;
t325 = t21 * t265 + t22 * t266;
t310 = -qJD(1) * t328 - t260 + t378;
t296 = qJD(1) * t307 + t310;
t44 = -t257 * t122 - t200 * t371 + t296;
t290 = qJD(1) * t345 + t261;
t45 = t290 + t544;
t321 = -t265 * t44 + t266 * t45;
t320 = -t265 * t73 - t266 * t75;
t317 = -t122 * t266 + t126 * t265;
t309 = t365 - t430;
t295 = qJD(1) * t205 + t276 * t450 + t395;
t294 = -t468 - pkin(2) + (-rSges(5,3) - pkin(6)) * t270;
t291 = t336 + t470;
t289 = -t364 + t484;
t281 = rSges(4,3) * t266 - t265 * t335;
t278 = t270 * (-t353 * t266 + (t266 * t289 + t410) * t265);
t254 = rSges(3,2) * t390;
t239 = rSges(4,2) * t372;
t166 = qJD(1) * t336 + t378;
t150 = -rSges(5,1) * t203 - rSges(5,2) * t204;
t148 = rSges(5,1) * t201 + rSges(5,2) * t202;
t102 = qJD(1) * t346 + t261;
t101 = -qJD(1) * t281 + t310;
t50 = t317 * t387 + qJD(2);
t49 = t346 * qJDD(1) + (-rSges(4,3) * t389 - t164 + (qJD(1) * t335 - qJD(3)) * t265) * qJD(1) + t392;
t48 = qJDD(1) * t281 + (qJD(1) * t308 + t239) * qJD(1) + t280;
t30 = t290 + t546;
t29 = -t243 + t296 - t405 * t371 + (rSges(6,3) * t429 - t265 * t306 - t353) * t257;
t27 = qJD(4) * t278 + qJD(2) - t384;
t20 = -t122 * t169 - t126 * t168 + t320 * t387 + qJDD(2);
t1 = -qJDD(5) * t271 + qJDD(2) + (t265 * t289 - t353) * t169 - t410 * t168 + (t265 * t460 - t266 * t459) * t387;
t2 = [t444 / 0.2e1 + t441 / 0.2e1 + t442 / 0.2e1 + t443 / 0.2e1 + t501 * t477 + t500 * t476 + ((t551 * t266 + (-t504 + t516 - t553) * t265) * t387 + t514) * t348 + (-t166 * t254 + (qJD(1) * t213 * t291 - t166 * t339) * qJD(1) + (t277 * t336 - g(2) + t264) * t213 + (pkin(1) * t419 + (t454 * qJD(1) - t254) * qJD(1) + g(3)) * t291) * m(3) + (g(2) * t246 + g(3) * t337) * m(2) + (m(2) * (t246 ^ 2 + t337 ^ 2) + Icges(4,2) * t478 + (Icges(4,1) * t270 + 0.2e1 * Icges(4,4) * t271) * t270 + m(3) * (t213 ^ 2 + t291 ^ 2) + Icges(2,3) + Icges(3,3)) * qJDD(1) + (t30 * (-t265 * t376 - t391 + t483) - t29 * (t261 + t481) + (t11 * (t365 + t422) + t30 * t354) * t266 + ((t274 * t30 + t276 * t29) * pkin(1) + (-t29 * t309 + t30 * t358) * t266 + (t30 * (-t314 + t451) - t29 * t358) * t265) * qJD(1) - (-pkin(2) + (-rSges(6,3) + t272) * t270) * t465 - (t30 + t295 - t546) * t29 + t510 * (t265 * t309 + t329 + t353) + t543 * (t330 + t547)) * m(6) + ((-g(3) + t21) * (t265 * t294 + t329 + t404) + (-t334 - t391 - t165 + (t381 - t329) * qJD(1)) * t45 + (-t261 - t408 + (-t330 + (rSges(5,3) * t270 + pkin(2) + t344) * t266) * qJD(1) - t295 - t45 + t544) * t44 + (-g(2) + t22) * (t266 * t294 + t330 - t333)) * m(5) + (-(t102 + (t160 + t469) * qJD(1) + t395) * t101 - t102 * t391 - t101 * (t239 + t261) + ((t101 * t276 + t102 * t274) * pkin(1) + (-t101 * t366 - t102 * t440) * t266 + (t101 * t440 + t102 * t335) * t265) * qJD(1) + (t49 - g(2)) * (-t265 * t440 + t266 * t366 + t244 - t469) + (t48 - g(3)) * (-t470 + t440 * t266 + (-pkin(2) - t335) * t265)) * m(4) + (((t552 + t553) * t266 + (t505 + t551) * t265 + (-t557 * t265 ^ 2 + (-t527 * t265 - t266 * t557) * t266) * t270 + t491) * t387 + t508 - t511) * t349 + (t507 + t509 + t541) * t347 + t485 + ((-t557 * t271 + (-t273 * t498 - t275 * t548) * t270 - t502) * t257 - t506 - t542) * t350; (m(3) + m(4)) * qJDD(2) + m(5) * t20 + m(6) * t1 + (-m(3) + t383) * g(1); t383 * (g(3) * t265 + t465) + m(4) * (t265 * t48 + t266 * t49) + m(5) * t325 + m(6) * (t10 * t265 + t11 * t266); (-t271 * t501 + t517) * t477 + (-t271 * t500 + t518) * t476 + (-t532 * t271 + (t265 * t503 + t266 * t502) * t270) * t256 / 0.2e1 - (((-t273 * t487 - t275 * t486) * t257 + ((t273 * t495 + t275 * t494) * t270 + t492 * t271) * qJD(4)) * t270 - t556 * t257 * t271) * t257 / 0.2e1 + (-t539 * t271 + ((-t265 * t502 + t266 * t503) * qJD(1) + t490) * t270) * t257 / 0.2e1 - (t387 * t490 + t441 + t442 + t443 + t444 + t485) * t271 / 0.2e1 - (-t168 * t516 + t540 * t169 + t501 * t256 - t506 * t257 + t489 * t387) * t429 / 0.2e1 + (t168 * t505 + t169 * t504 + t256 * t500 + t257 * t507 + t387 * t488) * t426 / 0.2e1 + (t506 * t271 + (qJD(1) * t491 + t489) * t270) * t350 + ((-t201 * t495 - t202 * t494 + t429 * t492) * t387 + (t201 * t487 + t202 * t486 - t429 * t556) * t257) * t349 + ((t203 * t495 + t204 * t494 - t426 * t492) * t387 + (-t203 * t487 - t204 * t486 + t426 * t556) * t257) * t348 + (-t507 * t271 + ((-t265 * t504 - t266 * t505) * qJD(1) + t488) * t270) * t347 + (t1 * t278 + t27 * ((qJD(1) * t410 - t459) * t266 + (((-t430 + t468 - t484) * t265 + t353) * qJD(1) + t460) * t265) * t270 + t11 * (t270 * t357 - t271 * t410) + t30 * (-t460 * t271 + (t266 * t396 - t390 * t405) * t270) + t10 * (-t353 * t271 + (-t478 * t464 + (t271 * t471 + t405) * t270) * t265) - t29 * (-t459 * t271 + (qJD(1) * t357 + t265 * t396) * t270) - (-t29 * t407 + t30 * t406) * t257 - (t27 * (t265 * t406 - t266 * t407) + (t270 * t467 - t301) * (t265 * t29 - t266 * t30)) * t387 - g(2) * t407 + g(3) * t406 - g(1) * (-t515 * t273 - t452) * t270) * m(6) + ((-t122 * t21 - t126 * t22 + t44 * t75 + t45 * t73) * t271 + (t20 * t317 + t50 * (t122 * t390 + t126 * t389 + t320) + t321 * t215 + ((-t265 * t45 - t266 * t44) * qJD(1) + t325) * t200) * t270 - (-t148 * t44 - t150 * t45) * t257 - (t50 * (-t148 * t266 - t150 * t265) + t321 * t227) * t387 - g(1) * t227 - g(2) * t148 - g(3) * t150) * m(5) - (t265 * t508 + t266 * t509) * t388 / 0.2e1; ((-t1 + g(1)) * t271 + (-t543 * t265 + t510 * t266) * t270) * m(6);];
tau = t2;
