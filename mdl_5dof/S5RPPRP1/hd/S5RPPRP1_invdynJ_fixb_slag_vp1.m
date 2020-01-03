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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:25:09
% EndTime: 2020-01-03 11:25:47
% DurationCPUTime: 24.60s
% Computational Cost: add. (12208->565), mult. (16242->714), div. (0->0), fcn. (15498->8), ass. (0->265)
t549 = Icges(5,4) + Icges(6,4);
t278 = qJ(1) + pkin(7);
t272 = sin(t278);
t273 = cos(t278);
t284 = cos(qJ(4));
t280 = cos(pkin(8));
t282 = sin(qJ(4));
t410 = t280 * t282;
t197 = -t272 * t284 + t273 * t410;
t409 = t280 * t284;
t416 = t272 * t282;
t198 = t273 * t409 + t416;
t166 = Icges(6,4) * t197;
t279 = sin(pkin(8));
t414 = t273 * t279;
t115 = Icges(6,1) * t198 + Icges(6,5) * t414 - t166;
t169 = Icges(5,4) * t197;
t118 = Icges(5,1) * t198 + Icges(5,5) * t414 - t169;
t498 = t115 + t118;
t167 = Icges(6,4) * t198;
t110 = Icges(6,2) * t197 - Icges(6,6) * t414 - t167;
t170 = Icges(5,4) * t198;
t113 = Icges(5,2) * t197 - Icges(5,6) * t414 - t170;
t499 = t110 + t113;
t558 = -t197 * t499 - t198 * t498;
t537 = Icges(5,1) + Icges(6,1);
t517 = Icges(5,5) + Icges(6,5);
t548 = Icges(5,2) + Icges(6,2);
t515 = Icges(5,6) + Icges(6,6);
t547 = Icges(5,3) + Icges(6,3);
t195 = -t272 * t410 - t273 * t284;
t418 = t272 * t279;
t412 = t273 * t282;
t196 = t272 * t409 - t412;
t428 = Icges(6,4) * t196;
t108 = Icges(6,2) * t195 + Icges(6,6) * t418 + t428;
t431 = Icges(5,4) * t196;
t111 = Icges(5,2) * t195 + Icges(5,6) * t418 + t431;
t512 = t111 + t108;
t165 = Icges(6,4) * t195;
t114 = Icges(6,1) * t196 + Icges(6,5) * t418 + t165;
t168 = Icges(5,4) * t195;
t117 = Icges(5,1) * t196 + Icges(5,5) * t418 + t168;
t511 = t117 + t114;
t555 = t549 * t284;
t554 = t549 * t282;
t103 = Icges(6,5) * t198 - Icges(6,6) * t197 + Icges(6,3) * t414;
t106 = Icges(5,5) * t198 - Icges(5,6) * t197 + Icges(5,3) * t414;
t513 = t103 + t106;
t488 = t414 * t513 - t558;
t553 = -t499 * t195 + t498 * t196;
t102 = Icges(6,5) * t196 + Icges(6,6) * t195 + Icges(6,3) * t418;
t105 = Icges(5,5) * t196 + Icges(5,6) * t195 + Icges(5,3) * t418;
t514 = t102 + t105;
t542 = -t547 * t280 + (-t282 * t515 + t284 * t517) * t279;
t541 = -t515 * t280 + (-t282 * t548 + t555) * t279;
t540 = -t517 * t280 + (t284 * t537 - t554) * t279;
t480 = (-t282 * t517 - t284 * t515) * t279;
t546 = (-t284 * t548 - t554) * t279;
t545 = (-t282 * t537 - t555) * t279;
t493 = -t197 * t512 + t198 * t511;
t503 = t512 * t195 + t511 * t196 + t514 * t418;
t502 = -t513 * t418 - t553;
t489 = -t514 * t414 - t493;
t485 = t541 * t195 + t540 * t196 + t542 * t418;
t132 = qJD(1) * t195 + qJD(4) * t198;
t297 = t197 * qJD(4);
t133 = qJD(1) * t196 + t297;
t378 = qJD(1) * t279;
t363 = t272 * t378;
t524 = t515 * t132 + t517 * t133 + t547 * t363;
t134 = -qJD(1) * t197 - qJD(4) * t196;
t296 = t195 * qJD(4);
t135 = qJD(1) * t198 + t296;
t362 = t273 * t378;
t523 = t515 * t134 + t517 * t135 + t547 * t362;
t543 = t548 * t132 + t549 * t133 + t515 * t363;
t521 = t548 * t134 + t549 * t135 + t515 * t362;
t520 = -t549 * t132 - t537 * t133 - t517 * t363;
t519 = -t549 * t134 - t537 * t135 - t517 * t362;
t533 = t480 * qJD(4);
t539 = t546 * qJD(4);
t531 = t545 * qJD(4);
t538 = -t541 * t197 + t540 * t198 + t542 * t414;
t530 = (t489 * t272 - t488 * t273) * t279;
t529 = (t503 * t272 - t502 * t273) * t279;
t528 = -t523 * t280 + (-t519 * t284 - t521 * t282 + (-t511 * t282 - t512 * t284) * qJD(4)) * t279;
t527 = t524 * t280 + (t520 * t284 + t543 * t282 + (-t498 * t282 + t499 * t284) * qJD(4)) * t279;
t526 = t533 * t280 + (-t531 * t284 + t539 * t282 + (t540 * t282 + t541 * t284) * qJD(4)) * t279;
t261 = -qJD(4) * t280 + qJD(1);
t525 = t485 * t261;
t518 = t542 * t280 + (t541 * t282 - t540 * t284) * t279;
t374 = qJD(5) * t279;
t239 = t272 * t374;
t451 = pkin(3) * t280;
t281 = -qJ(5) - pkin(6);
t447 = pkin(6) + t281;
t467 = t279 * t447;
t251 = pkin(4) * t416;
t267 = pkin(4) * t284 + pkin(3);
t413 = t273 * t280;
t509 = -t198 * rSges(6,1) + t197 * rSges(6,2) - rSges(6,3) * t414 - t267 * t413 - t251;
t398 = (t451 + t467) * t273 + t509;
t510 = -t398 * t261 + t239;
t508 = t538 * t261;
t492 = rSges(6,1) + pkin(4);
t370 = qJD(1) * qJD(4);
t164 = (-qJDD(4) * t273 + t272 * t370) * t279;
t260 = -qJDD(4) * t280 + qJDD(1);
t380 = qJD(1) * t272;
t262 = qJD(3) * t272;
t379 = qJD(1) * t273;
t383 = qJ(3) * t379 + t262;
t159 = pkin(2) * t380 - t383;
t283 = sin(qJ(1));
t275 = t283 * pkin(1);
t285 = cos(qJ(1));
t286 = qJD(1) ^ 2;
t425 = pkin(1) * qJDD(1);
t341 = -t275 * t286 + t285 * t425;
t371 = qJD(1) * qJD(3);
t312 = t272 * t371 + t341;
t339 = pkin(6) * t279 + t451;
t200 = t339 * t273;
t424 = qJ(3) * t272;
t227 = pkin(2) * t273 + t424;
t390 = t200 + t227;
t287 = (-t339 * t380 - t159) * qJD(1) + t390 * qJDD(1) + t312;
t440 = rSges(6,2) * t284;
t303 = (-rSges(6,1) * t282 - t440) * t279;
t373 = qJD(5) * t280;
t376 = qJD(4) * t279;
t450 = pkin(4) * t282;
t389 = -qJD(4) * t303 + t376 * t450 + t373;
t289 = -qJDD(3) + (qJD(1) * qJD(5) + qJD(4) * t389) * t279;
t369 = qJDD(5) * t279;
t448 = -pkin(3) + t267;
t395 = (t447 - rSges(6,3)) * t280 + (rSges(6,1) * t284 - rSges(6,2) * t282 + t448) * t279;
t368 = pkin(4) * t412;
t240 = t273 * t374;
t463 = -rSges(6,1) * t133 - rSges(6,2) * t132 + t240;
t444 = rSges(6,3) * t363 + pkin(4) * t297 + (-t368 + (t280 * t448 - t467) * t272) * qJD(1) - t463;
t11 = t164 * t395 - t260 * t398 - t261 * t444 + t272 * t369 + t273 * t289 + t287;
t507 = -g(2) + t11;
t163 = (qJDD(4) * t272 + t273 * t370) * t279;
t417 = t272 * t280;
t199 = pkin(3) * t417 + pkin(6) * t418;
t266 = t272 * pkin(2);
t225 = -qJ(3) * t273 + t266;
t377 = qJD(3) * t273;
t384 = pkin(2) * t379 + qJ(3) * t380;
t332 = -t377 + t384;
t276 = t285 * pkin(1);
t381 = t286 * t276 + t283 * t425;
t348 = qJD(1) * t332 + qJDD(1) * t225 + t381;
t361 = t280 * t379;
t385 = pkin(3) * t361 + pkin(6) * t362;
t308 = qJD(1) * t385 + qJDD(1) * t199 + t348;
t411 = t279 * t281;
t466 = t196 * rSges(6,1) + t195 * rSges(6,2) + rSges(6,3) * t418 + t267 * t417 - t272 * t411;
t399 = -t368 - t199 + t466;
t468 = t135 * rSges(6,1) + t134 * rSges(6,2) + rSges(6,3) * t362 + qJD(1) * t251 + t267 * t361 + t239;
t443 = pkin(4) * t296 - t281 * t362 - t385 + t468;
t10 = (-t369 - t371) * t273 + t443 * t261 + t399 * t260 - t395 * t163 + t289 * t272 + t308;
t506 = g(3) - t10;
t505 = t529 * qJD(4) + t525;
t504 = t530 * qJD(4) - t508;
t491 = (t533 * t273 - t380 * t542) * t279 + t531 * t198 - t539 * t197 - t540 * t133 - t541 * t132;
t490 = (-t533 * t272 - t379 * t542) * t279 - t531 * t196 - t539 * t195 - t540 * t135 - t541 * t134;
t484 = t280 * t103;
t40 = t484 + (-t110 * t282 - t115 * t284) * t279;
t483 = t280 * t106;
t42 = t483 + (-t113 * t282 - t118 * t284) * t279;
t486 = t40 + t42;
t125 = t198 * rSges(5,1) - t197 * rSges(5,2) + rSges(5,3) * t414;
t501 = t125 * t261;
t194 = -rSges(5,3) * t280 + (rSges(5,1) * t284 - rSges(5,2) * t282) * t279;
t360 = t194 * t376;
t226 = t272 * rSges(3,1) + t273 * rSges(3,2);
t207 = t275 + t226;
t497 = t502 * t272 + t503 * t273;
t495 = t513 * t272;
t494 = t514 * t273;
t477 = (-t537 * t197 - t167 - t170 + t499) * t273 + (t537 * t195 - t428 - t431 - t512) * t272;
t39 = -t102 * t280 + (-t108 * t282 + t114 * t284) * t279;
t41 = -t105 * t280 + (-t111 * t282 + t117 * t284) * t279;
t487 = t39 + t41;
t478 = (-t198 * t548 - t166 - t169 + t498) * t273 + (-t196 * t548 + t165 + t168 + t511) * t272;
t476 = (t197 * t517 + t198 * t515) * t273 + (-t195 * t517 + t196 * t515) * t272;
t475 = t272 * t528 + t527 * t273;
t474 = (-t519 * t196 + t521 * t195 + t511 * t135 + t512 * t134 + (t272 * t523 + t379 * t514) * t279) * t272 + ((-t272 * t524 + t379 * t513) * t279 + t520 * t196 - t543 * t195 + t498 * t135 - t499 * t134) * t273;
t473 = (-t520 * t198 - t543 * t197 + t498 * t133 - t499 * t132 + (t273 * t524 + t380 * t513) * t279) * t273 + ((-t273 * t523 + t380 * t514) * t279 + t519 * t198 + t521 * t197 + t511 * t133 + t512 * t132) * t272;
t472 = -rSges(4,1) * t417 + rSges(4,2) * t418;
t471 = t540 + t546;
t470 = t541 - t545;
t469 = -t518 * t260 - t526 * t261;
t353 = t225 + t275;
t464 = t199 + t353;
t439 = rSges(4,3) * t273;
t462 = t439 + t472;
t457 = t163 / 0.2e1;
t456 = t164 / 0.2e1;
t449 = g(2) * t273;
t436 = t39 * t163;
t435 = t40 * t164;
t434 = t41 * t163;
t433 = t42 * t164;
t432 = rSges(4,3) + qJ(3);
t397 = -t196 * rSges(6,2) + t195 * t492;
t396 = t198 * rSges(6,2) + t197 * t492;
t246 = rSges(4,2) * t414;
t155 = rSges(4,1) * t413 + rSges(4,3) * t272 - t246;
t387 = t227 + t155;
t386 = rSges(4,1) * t361 + rSges(4,3) * t380;
t382 = t266 + t275;
t372 = -m(4) - m(5) - m(6);
t75 = t135 * rSges(5,1) + t134 * rSges(5,2) + rSges(5,3) * t362;
t123 = t196 * rSges(5,1) + t195 * rSges(5,2) + rSges(5,3) * t418;
t212 = qJD(1) * t227;
t271 = qJD(1) * t276;
t365 = qJD(1) * t200 + t212 + t271;
t357 = rSges(4,1) * t280 + pkin(2);
t355 = -t376 / 0.2e1;
t354 = t376 / 0.2e1;
t298 = qJD(1) * t464 - t262;
t333 = t395 * t376;
t29 = t261 * t399 - t272 * t333 - t240 + t298;
t352 = t29 * t395;
t228 = rSges(3,1) * t273 - t272 * rSges(3,2);
t351 = t276 + t424;
t350 = t271 - t377;
t345 = t272 * t355;
t344 = t272 * t354;
t343 = t273 * t355;
t342 = t273 * t354;
t248 = rSges(2,1) * t285 - t283 * rSges(2,2);
t247 = rSges(2,1) * t283 + rSges(2,2) * t285;
t330 = rSges(5,1) * t133 + rSges(5,2) * t132;
t221 = (-rSges(5,1) * t282 - rSges(5,2) * t284) * t279;
t210 = qJD(4) * t221;
t320 = -t210 * t376 - qJDD(3);
t21 = t123 * t260 - t163 * t194 + t261 * t75 + t272 * t320 - t273 * t371 + t308;
t73 = rSges(5,3) * t363 + t330;
t22 = t125 * t260 + t164 * t194 - t261 * t73 + t273 * t320 + t287;
t326 = -t21 * t272 - t22 * t273;
t44 = t123 * t261 - t272 * t360 + t298;
t313 = qJD(1) * t390 + t271;
t45 = t501 + (-qJD(3) - t360) * t273 + t313;
t322 = -t272 * t44 - t273 * t45;
t321 = t272 * t73 + t273 * t75;
t317 = t123 * t273 - t125 * t272;
t311 = pkin(2) + t339;
t310 = t271 + t332;
t208 = t228 + t276;
t151 = rSges(5,1) * t197 + rSges(5,2) * t198;
t149 = rSges(5,1) * t195 - rSges(5,2) * t196;
t101 = qJD(1) * t387 + t350;
t50 = t317 * t376 + qJD(2);
t49 = -qJDD(3) * t273 + t387 * qJDD(1) + (qJD(1) * t462 - t159) * qJD(1) + t312;
t48 = -qJDD(1) * t462 - qJDD(3) * t272 + ((-rSges(4,2) * t378 - qJD(3)) * t273 + t386) * qJD(1) + t348;
t30 = (-qJD(3) - t333) * t273 + t313 + t510;
t27 = -t373 + qJD(2) + (t272 * t398 + t273 * t399) * t376;
t20 = -t123 * t164 - t125 * t163 + t321 * t376 + qJDD(2);
t1 = -qJDD(5) * t280 + qJDD(2) - t399 * t164 + t398 * t163 + (t272 * t444 + t273 * t443) * t376;
t2 = [t434 / 0.2e1 + t435 / 0.2e1 + t436 / 0.2e1 + t433 / 0.2e1 + t469 - m(2) * (g(2) * t248 + g(3) * t247) + t485 * t457 - t538 * t456 + ((((t494 - t495) * t279 + t493 + t489 - t502) * t273 + (t488 + t503 + t558) * t272) * t376 + t525) * t342 + ((-t226 * t286 - g(2) + t341) * t208 + (t381 - g(3) + (0.2e1 * rSges(3,1) * t379 - 0.2e1 * rSges(3,2) * t380 - qJD(1) * t228) * qJD(1)) * t207) * m(3) + (m(3) * (t226 * t207 + t228 * t208) + Icges(4,2) * t280 ^ 2 + (Icges(4,1) * t279 + 0.2e1 * Icges(4,4) * t280) * t279 + m(2) * (t247 ^ 2 + t248 ^ 2) + Icges(2,3) + Icges(3,3)) * qJDD(1) + (-(-t30 + t365 + t510) * t29 - (-t29 * qJD(3) - t352 * t376) * t273 + t30 * (t383 + t463) + t29 * (t310 + t468) + (t195 * t29 - t197 * t30) * qJD(4) * pkin(4) + (-t29 * t273 * t411 + (t368 - t275 + (-t267 * t280 - pkin(2) + (-rSges(6,3) + t281) * t279) * t272) * t30) * qJD(1) - t506 * ((-qJ(3) - t450) * t273 + t382 + t466) + t507 * ((pkin(2) - t411) * t273 + t351 - t509)) * m(6) + ((-t330 + t383 + (-t275 + (-t451 - pkin(2) + (-rSges(5,3) - pkin(6)) * t279) * t272) * qJD(1)) * t45 - t311 * t449 + (-g(3) + t21) * (t123 + t464) + (t271 - t365 + t384 + t385 + t45 + t75 - t501) * t44 + (t22 * t311 + t360 * t44) * t273 + (-g(2) + t22) * (t351 + t125)) * m(5) + ((t49 - g(2)) * (t272 * t432 + t273 * t357 - t246 + t276) + (t48 - g(3)) * (-t273 * t432 + t382 - t472) + (t383 + (t439 - t275 + (rSges(4,2) * t279 - t357) * t272) * qJD(1)) * t101 + (t101 - t212 + t310 - t350 + t386 + (-t155 - t246) * qJD(1)) * (-t262 + (-t462 + t353) * qJD(1))) * m(4) + (-t490 + t528) * t344 + (-t491 + t505 - t527) * t343 + ((-t483 - t484 + (t282 * t499 + t284 * t498) * t279 + t486) * t261 + (((t494 + t495) * t279 + t493 + t553) * t272 + t497 + (-t503 + t488) * t273) * t376 + t504 + t508) * t345; (m(3) + m(4)) * qJDD(2) + m(5) * t20 + m(6) * t1 + (-m(3) + t372) * g(1); t372 * (-g(3) * t272 - t449) + m(4) * (-t272 * t48 - t273 * t49) + m(5) * t326 + m(6) * (-t10 * t272 - t11 * t273); (-t280 * t485 + t529) * t457 + (t280 * t538 + t530) * t456 + (t518 * t280 + (t272 * t487 - t273 * t486) * t279) * t260 / 0.2e1 - (((-t282 * t471 - t284 * t470) * t261 + ((-t282 * t478 + t284 * t477) * t279 + t476 * t280) * qJD(4)) * t279 - t480 * t261 * t280) * t261 / 0.2e1 + (t526 * t280 + ((t272 * t486 + t273 * t487) * qJD(1) + t475) * t279) * t261 / 0.2e1 - (t376 * t475 + t433 + t434 + t435 + t436 + t469) * t280 / 0.2e1 + (t503 * t163 + t502 * t164 + t485 * t260 - t490 * t261 + t474 * t376) * t418 / 0.2e1 - (t163 * t489 + t164 * t488 - t260 * t538 - t491 * t261 + t473 * t376) * t414 / 0.2e1 + ((t195 * t478 + t196 * t477 - t418 * t476) * t376 + (t195 * t471 - t196 * t470 + t418 * t480) * t261) * t345 + (t490 * t280 + (qJD(1) * t497 + t474) * t279) * t344 + (t491 * t280 + ((t272 * t488 + t273 * t489) * qJD(1) + t473) * t279) * t343 + ((t478 * t197 - t198 * t477 + t476 * t414) * t376 + (t197 * t471 + t198 * t470 - t414 * t480) * t261) * t342 + (-(t29 * t397 - t30 * t396) * t261 - (t27 * (t272 * t396 + t273 * t397) + (t279 * t450 - t303) * (t272 * t29 + t273 * t30)) * t376 - g(2) * t397 - g(3) * t396 - g(1) * (-t282 * t492 - t440) * t279 + (-t10 * t399 + t11 * t398 - t29 * t443 + t30 * t444) * t280 + ((t1 * t399 + t27 * t443 - t11 * t395 + t30 * t389 + (t27 * t398 - t352) * qJD(1)) * t273 + (t1 * t398 + t27 * t444 - t10 * t395 + t29 * t389 + (-t27 * t399 + t30 * t395) * qJD(1)) * t272) * t279) * m(6) + (-(t149 * t44 - t151 * t45) * t261 - (t50 * (t149 * t273 + t151 * t272) + t322 * t221) * t376 + (-t123 * t21 - t125 * t22 - t44 * t75 + t45 * t73) * t280 + (t20 * t317 + t50 * (-t123 * t380 - t125 * t379 + t321) + t322 * t210 + ((t272 * t45 - t273 * t44) * qJD(1) + t326) * t194) * t279 - g(1) * t221 - g(2) * t149 - g(3) * t151) * m(5) + (t272 * t504 + t273 * t505) * t378 / 0.2e1; ((-t1 + g(1)) * t280 + (t272 * t507 + t273 * t506) * t279) * m(6);];
tau = t2;
