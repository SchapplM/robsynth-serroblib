% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:36
% EndTime: 2019-12-05 15:38:04
% DurationCPUTime: 16.26s
% Computational Cost: add. (24733->655), mult. (35335->964), div. (0->0), fcn. (37951->8), ass. (0->383)
t597 = rSges(6,1) + pkin(4);
t594 = rSges(6,3) + qJ(5);
t342 = pkin(8) + qJ(4);
t339 = cos(t342);
t345 = sin(pkin(7));
t347 = cos(pkin(7));
t338 = sin(t342);
t350 = cos(qJ(2));
t472 = t350 * t338;
t288 = t347 * t339 + t345 * t472;
t484 = t347 * t338;
t485 = t345 * t350;
t289 = t339 * t485 - t484;
t349 = sin(qJ(2));
t486 = t345 * t349;
t178 = Icges(5,5) * t289 - Icges(5,6) * t288 + Icges(5,3) * t486;
t180 = Icges(6,4) * t289 + Icges(6,2) * t486 + Icges(6,6) * t288;
t627 = t178 + t180;
t290 = -t345 * t339 + t347 * t472;
t481 = t347 * t350;
t488 = t345 * t338;
t291 = t339 * t481 + t488;
t482 = t347 * t349;
t179 = Icges(5,5) * t291 - Icges(5,6) * t290 + Icges(5,3) * t482;
t181 = Icges(6,4) * t291 + Icges(6,2) * t482 + Icges(6,6) * t290;
t626 = t179 + t181;
t513 = Icges(5,4) * t339;
t401 = -Icges(5,2) * t338 + t513;
t270 = -Icges(5,6) * t350 + t401 * t349;
t474 = t349 * t339;
t329 = Icges(6,5) * t474;
t490 = t338 * t349;
t505 = Icges(6,6) * t350;
t612 = Icges(6,3) * t490 - t270 + t329 - t505;
t398 = Icges(5,5) * t339 - Icges(5,6) * t338;
t266 = -Icges(5,3) * t350 + t398 * t349;
t400 = Icges(6,4) * t339 + Icges(6,6) * t338;
t268 = -Icges(6,2) * t350 + t400 * t349;
t625 = t266 + t268;
t508 = Icges(6,5) * t338;
t404 = Icges(6,1) * t339 + t508;
t272 = -Icges(6,4) * t350 + t404 * t349;
t514 = Icges(5,4) * t338;
t405 = Icges(5,1) * t339 - t514;
t274 = -Icges(5,5) * t350 + t405 * t349;
t611 = t272 + t274;
t343 = t349 ^ 2;
t563 = t350 ^ 2;
t473 = t349 * t350;
t199 = -Icges(5,5) * t288 - Icges(5,6) * t289;
t201 = -Icges(6,4) * t288 + Icges(6,6) * t289;
t624 = t199 + t201;
t200 = -Icges(5,5) * t290 - Icges(5,6) * t291;
t202 = -Icges(6,4) * t290 + Icges(6,6) * t291;
t623 = t200 + t202;
t283 = Icges(5,4) * t290;
t187 = Icges(5,1) * t291 + Icges(5,5) * t482 - t283;
t459 = Icges(5,2) * t291 - t187 + t283;
t509 = Icges(6,5) * t290;
t185 = Icges(6,1) * t291 + Icges(6,4) * t482 + t509;
t461 = Icges(6,3) * t291 - t185 - t509;
t621 = t459 + t461;
t282 = Icges(5,4) * t288;
t186 = Icges(5,1) * t289 + Icges(5,5) * t486 - t282;
t460 = Icges(5,2) * t289 - t186 + t282;
t510 = Icges(6,5) * t288;
t184 = Icges(6,1) * t289 + Icges(6,4) * t486 + t510;
t462 = Icges(6,3) * t289 - t184 - t510;
t620 = t460 + t462;
t515 = Icges(5,4) * t291;
t183 = -Icges(5,2) * t290 + Icges(5,6) * t482 + t515;
t463 = -Icges(5,1) * t290 - t183 - t515;
t281 = Icges(6,5) * t291;
t177 = Icges(6,6) * t482 + Icges(6,3) * t290 + t281;
t465 = -Icges(6,1) * t290 + t177 + t281;
t619 = t463 + t465;
t516 = Icges(5,4) * t289;
t182 = -Icges(5,2) * t288 + Icges(5,6) * t486 + t516;
t464 = -Icges(5,1) * t288 - t182 - t516;
t280 = Icges(6,5) * t289;
t176 = Icges(6,6) * t486 + Icges(6,3) * t288 + t280;
t466 = -Icges(6,1) * t288 + t176 + t280;
t618 = t464 + t466;
t616 = t176 - t182;
t615 = t177 - t183;
t614 = t184 + t186;
t613 = t185 + t187;
t610 = t612 * t338 + t611 * t339;
t609 = t594 * t338 + t597 * t339;
t608 = t625 * t349;
t607 = t626 * t349;
t606 = t627 * t349;
t605 = t620 * t288 + t618 * t289 + t624 * t486;
t604 = t621 * t288 + t619 * t289 + t623 * t486;
t603 = t620 * t290 + t618 * t291 + t624 * t482;
t602 = t621 * t290 + t619 * t291 + t623 * t482;
t601 = t611 + (t508 - t514 + (-Icges(5,2) - Icges(6,3)) * t339) * t349;
t600 = -(-Icges(5,1) * t338 - t513) * t349 + Icges(6,1) * t490 - t329 - t612;
t295 = (-Icges(5,5) * t338 - Icges(5,6) * t339) * t349;
t296 = (-Icges(6,4) * t338 + Icges(6,6) * t339) * t349;
t599 = (-t295 - t296) * t350;
t579 = t625 * t350;
t598 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t473 + (0.2e1 * t563 - 0.2e1 * t343) * Icges(3,4);
t557 = m(5) / 0.2e1;
t555 = m(6) / 0.2e1;
t397 = Icges(6,5) * t339 + Icges(6,3) * t338;
t357 = -t397 * t349 + t505;
t235 = t357 * t345;
t241 = t270 * t345;
t243 = t272 * t345;
t245 = t274 * t345;
t390 = -t182 * t338 + t186 * t339;
t378 = -t266 * t345 - t390;
t392 = t176 * t338 + t184 * t339;
t380 = t268 * t345 + t392;
t596 = (t380 - t378) * t350 + ((-t243 - t245) * t339 + (t235 + t241) * t338 + t627) * t349;
t236 = t357 * t347;
t242 = t270 * t347;
t244 = t272 * t347;
t246 = t274 * t347;
t389 = -t183 * t338 + t187 * t339;
t377 = -t266 * t347 - t389;
t391 = t177 * t338 + t185 * t339;
t379 = t268 * t347 + t391;
t595 = (-t377 + t379) * t350 + ((-t244 - t246) * t339 + (t236 + t242) * t338 + t626) * t349;
t593 = t616 * t288 + t614 * t289 + t606 * t345;
t592 = t615 * t288 + t613 * t289 + t607 * t345;
t399 = -Icges(3,5) * t349 - Icges(3,6) * t350;
t316 = t399 * t345;
t317 = t399 * t347;
t591 = t616 * t290 + t614 * t291 + t606 * t347;
t590 = t615 * t290 + t613 * t291 + t607 * t347;
t589 = t612 * t288 + t611 * t289 + t608 * t345;
t588 = t612 * t290 + t611 * t291 + t608 * t347;
t587 = t610 * t349 - t579;
t586 = (t397 - t401) * t350 + (-Icges(5,6) + Icges(6,6)) * t349;
t585 = (t404 + t405) * t350 + (Icges(6,4) + Icges(5,5)) * t349;
t340 = t345 ^ 2;
t341 = t347 ^ 2;
t436 = t340 + t341;
t584 = rSges(6,2) * t350 - t609 * t349;
t583 = ((t398 + t400) * t350 + (Icges(6,2) + Icges(5,3)) * t349 - t610) * t350;
t500 = t180 * t350;
t112 = t392 * t349 - t500;
t499 = t181 * t350;
t113 = t391 * t349 - t499;
t502 = t178 * t350;
t114 = t390 * t349 - t502;
t501 = t179 * t350;
t115 = t389 * t349 - t501;
t582 = (t113 + t115) * t347 + (t112 + t114) * t345;
t189 = rSges(5,1) * t289 - rSges(5,2) * t288 + rSges(5,3) * t486;
t191 = rSges(5,1) * t291 - rSges(5,2) * t290 + rSges(5,3) * t482;
t132 = (t189 * t347 - t191 * t345) * t349;
t322 = t436 * t349;
t410 = rSges(5,1) * t339 - rSges(5,2) * t338;
t277 = -rSges(5,3) * t350 + t410 * t349;
t498 = t189 * t350;
t151 = t277 * t486 + t498;
t497 = t191 * t350;
t152 = -t277 * t482 - t497;
t468 = t151 * t481 + t152 * t485;
t458 = rSges(6,2) * t486 + t594 * t288 + t597 * t289;
t419 = t458 * t350;
t118 = -t486 * t584 + t419;
t457 = rSges(6,2) * t482 + t594 * t290 + t597 * t291;
t418 = t457 * t350;
t119 = t482 * t584 - t418;
t470 = t118 * t481 + t119 * t485;
t95 = (-t457 * t345 + t458 * t347) * t349;
t526 = (t322 * t95 + t470) * t555 + (t132 * t322 + t468) * t557;
t248 = t277 * t345;
t250 = t277 * t347;
t109 = (-t248 * t349 + t498) * t347 + (t250 * t349 - t497) * t345;
t279 = t349 * rSges(5,3) + t410 * t350;
t385 = t277 * t350 + t279 * t349;
t122 = -t349 * t189 - t248 * t350 + t385 * t345;
t123 = t349 * t191 + t250 * t350 - t385 * t347;
t452 = t584 * t347;
t453 = t584 * t345;
t57 = (t453 * t349 + t419) * t347 + (-t452 * t349 - t418) * t345;
t443 = t349 * rSges(6,2) + t609 * t350;
t566 = t443 * t349 - t350 * t584;
t83 = t566 * t345 - t458 * t349 + t453 * t350;
t84 = -t566 * t347 + t457 * t349 - t452 * t350;
t527 = (-t350 * t57 + (t345 * t84 + t347 * t83 + t95) * t349 + t470) * t555 + (-t109 * t350 + (t122 * t347 + t123 * t345 + t132) * t349 + t468) * t557;
t578 = t526 - t527;
t577 = -t345 / 0.2e1;
t537 = t345 / 0.2e1;
t536 = -t347 / 0.2e1;
t535 = t347 / 0.2e1;
t576 = -t349 / 0.2e1;
t534 = t349 / 0.2e1;
t575 = -t350 / 0.2e1;
t573 = (t601 * t288 + t600 * t289) * t350 + (t604 * t347 + (t599 + t605) * t345) * t349;
t572 = (t601 * t290 + t600 * t291) * t350 + ((t599 + t602) * t347 + t603 * t345) * t349;
t571 = t339 * t350;
t554 = m(6) / 0.4e1;
t556 = m(5) / 0.4e1;
t417 = m(4) / 0.4e1 + t556 + t554;
t437 = t436 * t473;
t570 = t417 * (t437 - t473);
t192 = t288 * t345 + t290 * t347;
t159 = (-t192 + t472) * t490;
t420 = t472 / 0.2e1;
t171 = (t420 - t192 / 0.2e1) * m(6);
t523 = m(6) * qJD(5);
t569 = t171 * qJD(1) + t159 * t523;
t455 = -t597 * t290 + t594 * t291;
t456 = -t597 * t288 + t594 * t289;
t110 = t456 * t345 + t455 * t347;
t217 = -rSges(5,1) * t288 - rSges(5,2) * t289;
t221 = -rSges(5,1) * t290 - rSges(5,2) * t291;
t138 = t217 * t345 + t221 * t347;
t440 = (-t597 * t338 + t594 * t339) * t349;
t195 = t440 * t345;
t196 = t440 * t347;
t302 = (-rSges(5,1) * t338 - rSges(5,2) * t339) * t349;
t519 = (-t138 * t350 - t302 * t322) * t557 + (-t110 * t350 + (-t195 * t345 - t196 * t347) * t349) * t555;
t344 = sin(pkin(8));
t346 = cos(pkin(8));
t411 = rSges(4,1) * t346 - rSges(4,2) * t344;
t567 = rSges(4,3) * t350 - t411 * t349;
t528 = pkin(3) * t346;
t370 = pkin(6) * t350 - t528 * t349;
t312 = -t344 * t485 - t347 * t346;
t483 = t347 * t344;
t313 = t346 * t485 - t483;
t314 = -t344 * t481 + t345 * t346;
t487 = t345 * t344;
t315 = t346 * t481 + t487;
t565 = (-(Icges(4,5) * t313 + Icges(4,6) * t312 + Icges(4,3) * t486) * t347 + (Icges(4,5) * t315 + Icges(4,6) * t314 + Icges(4,3) * t482) * t345) * t350 + (-((Icges(4,4) * t313 + Icges(4,2) * t312 + Icges(4,6) * t486) * t344 - (Icges(4,1) * t313 + Icges(4,4) * t312 + Icges(4,5) * t486) * t346) * t347 + ((Icges(4,4) * t315 + Icges(4,2) * t314 + Icges(4,6) * t482) * t344 - (Icges(4,1) * t315 + Icges(4,4) * t314 + Icges(4,5) * t482) * t346) * t345) * t349;
t562 = 2 * qJD(2);
t561 = 4 * qJD(2);
t560 = 2 * qJD(4);
t559 = 4 * qJD(4);
t558 = m(4) / 0.2e1;
t327 = pkin(2) * t350 + t349 * qJ(3);
t438 = t436 * t327;
t127 = t345 * (rSges(4,1) * t313 + rSges(4,2) * t312 + rSges(4,3) * t486) + t347 * (rSges(4,1) * t315 + rSges(4,2) * t314 + rSges(4,3) * t482) + t438;
t325 = t349 * pkin(2) - qJ(3) * t350;
t442 = -t325 + t567;
t228 = t442 * t345;
t230 = t442 * t347;
t454 = t228 * t485 + t230 * t481;
t553 = m(4) * (t127 * t322 + t454);
t551 = m(5) * (t109 * t132 + t122 * t151 + t123 * t152);
t450 = -t325 + t370;
t423 = -t277 + t450;
t163 = t423 * t345;
t165 = t423 * t347;
t262 = pkin(6) * t349 + t528 * t350;
t415 = t345 * (-pkin(3) * t483 + t262 * t345) + t347 * (pkin(3) * t487 + t262 * t347) + t438;
t93 = t189 * t345 + t191 * t347 + t415;
t550 = m(5) * (t138 * t93 + (-t163 * t345 - t165 * t347) * t302);
t467 = t163 * t485 + t165 * t481;
t549 = m(5) * (t322 * t93 + t467);
t394 = -t118 * t347 - t119 * t345;
t547 = m(6) * (t288 * t84 + t290 * t83 + (t350 * t95 + (t394 + t57) * t349) * t338);
t545 = m(6) * (t118 * t83 + t119 * t84 + t57 * t95);
t413 = t584 + t450;
t144 = t413 * t345;
t146 = t413 * t347;
t174 = (t288 * t347 - t290 * t345) * t349;
t232 = t288 * t350 + t343 * t488;
t233 = -t290 * t350 - t343 * t484;
t58 = t458 * t345 + t457 * t347 + t415;
t544 = m(6) * (t144 * t233 + t146 * t232 + t174 * t58 + t192 * t95 + t394 * t490);
t543 = m(6) * (t144 * t289 + t146 * t291 - t195 * t288 - t196 * t290 + (t110 * t338 + t339 * t58) * t349);
t542 = m(6) * (t110 * t58 - t144 * t195 - t146 * t196);
t469 = t144 * t485 + t146 * t481;
t540 = m(6) * (t322 * t58 + t469);
t451 = t288 * t485 + t290 * t481;
t531 = m(6) * ((-t563 + (0.1e1 - t436) * t343) * t338 + t451);
t530 = m(6) * (t322 * t490 + t451);
t489 = t343 * t338;
t529 = m(6) * (-t192 * t350 - t436 * t489);
t525 = m(6) * qJD(2);
t524 = m(6) * qJD(4);
t503 = t174 * t350;
t29 = 0.2e1 * (t57 / 0.4e1 - t110 / 0.4e1) * m(6) + 0.2e1 * (t109 / 0.4e1 - t138 / 0.4e1) * m(5);
t494 = t29 * qJD(1);
t449 = -t262 - t327;
t441 = -t349 * rSges(4,3) - t411 * t350 - t327;
t439 = t436 * t325;
t435 = qJD(4) * t349;
t381 = 0.2e1 * t417 * t322;
t416 = t555 + t557 + t558;
t170 = -t416 * t349 + t381;
t434 = t170 * qJD(1);
t352 = -t380 * t349 + t500;
t59 = t288 * t235 - t289 * t243 + t352 * t345;
t351 = -t379 * t349 + t499;
t60 = t288 * t236 - t289 * t244 + t351 * t345;
t354 = t378 * t349 + t502;
t61 = t288 * t241 - t289 * t245 + t354 * t345;
t353 = t377 * t349 + t501;
t62 = t288 * t242 - t289 * t246 + t353 * t345;
t431 = (t586 * t288 + t585 * t289) * t575 + ((-t579 + t593) * t350 + (t61 + t59 - t583) * t349) * t537 + (t592 * t350 + (t60 + t62) * t349) * t535 + t589 * t534;
t63 = t290 * t235 - t291 * t243 + t352 * t347;
t64 = t290 * t236 - t291 * t244 + t351 * t347;
t65 = t290 * t241 - t291 * t245 + t354 * t347;
t66 = t290 * t242 - t291 * t246 + t353 * t347;
t430 = (t586 * t290 + t585 * t291) * t575 + (t591 * t350 + (t63 + t65) * t349) * t537 + ((-t579 + t590) * t350 + (t64 + t66 - t583) * t349) * t535 + t588 * t534;
t429 = ((-t586 * t338 - t585 * t339 - t625) * t350 + t595 * t347 + t596 * t345 + t587) * t576 + (t582 + t583) * t575;
t428 = t605 * t535 + t604 * t577;
t427 = t603 * t536 + t602 * t537;
t426 = t589 * t575 + (t593 * t345 + t592 * t347) * t534;
t425 = t588 * t575 + (t591 * t345 + t590 * t347) * t534;
t424 = t582 * t576 + t587 * t350 / 0.2e1;
t422 = -t279 + t449;
t421 = t474 / 0.2e1;
t414 = t436 * t370 - t439;
t412 = -t443 + t449;
t326 = t349 * rSges(3,1) + rSges(3,2) * t350;
t393 = -t144 * t345 - t146 * t347;
t388 = t232 * t347 + t233 * t345;
t167 = (t421 - t174 / 0.2e1) * m(6);
t367 = m(6) * (t289 * t345 + t291 * t347 - t571);
t97 = -m(6) * t503 / 0.2e1 + 0.2e1 * (t388 * t554 - t367 / 0.4e1) * t349;
t384 = -t167 * qJD(1) + t97 * qJD(3);
t383 = -t428 - t431;
t382 = t427 - t430;
t369 = t345 * t598 + t317;
t368 = -t347 * t598 + t316;
t362 = Icges(4,5) * t350 + (-Icges(4,1) * t346 + Icges(4,4) * t344) * t349;
t359 = Icges(4,6) * t350 + (-Icges(4,4) * t346 + Icges(4,2) * t344) * t349;
t256 = t362 * t347;
t255 = t362 * t345;
t254 = t359 * t347;
t253 = t359 * t345;
t231 = t441 * t347;
t229 = t441 * t345;
t227 = t436 * t326;
t172 = m(6) * t420 + t192 * t555;
t169 = t381 + (m(4) + m(5) + m(6)) * t534;
t168 = m(6) * t421 + t174 * t555;
t166 = t422 * t347;
t164 = t422 * t345;
t161 = -t221 * t350 - t302 * t482;
t160 = t217 * t350 + t302 * t486;
t154 = t529 / 0.2e1;
t153 = t530 / 0.2e1;
t150 = 0.4e1 * t570;
t147 = t412 * t347;
t145 = t412 * t345;
t140 = t288 * t289 + t290 * t291 + t339 * t489;
t139 = t531 / 0.2e1;
t134 = (t217 * t347 - t221 * t345) * t349;
t133 = t436 * t567 - t439;
t125 = -t349 * t196 - t455 * t350;
t124 = t456 * t350 + t440 * t486;
t111 = -t248 * t345 - t250 * t347 + t414;
t106 = (-t455 * t345 + t456 * t347) * t349;
t96 = (t388 * t349 - t503) * t555 + t367 * t534;
t92 = t453 * t345 + t452 * t347 + t414;
t91 = t153 + t139 - t529 / 0.2e1;
t90 = t154 + t153 - t531 / 0.2e1;
t89 = t154 + t139 - t530 / 0.2e1;
t88 = -t200 * t350 + (t459 * t338 + t463 * t339) * t349;
t87 = -t199 * t350 + (t460 * t338 + t464 * t339) * t349;
t86 = -t202 * t350 + (t461 * t338 + t465 * t339) * t349;
t85 = -t201 * t350 + (t462 * t338 + t466 * t339) * t349;
t44 = t192 * t58 + t393 * t490;
t43 = t118 * t232 + t119 * t233 + t174 * t95;
t38 = t345 * t66 - t347 * t65;
t37 = t345 * t64 - t347 * t63;
t36 = t345 * t62 - t347 * t61;
t35 = t345 * t60 - t347 * t59;
t30 = (t109 + t138) * t557 + (t110 + t57) * t555;
t23 = t543 / 0.2e1;
t22 = t540 + t549 + t553;
t18 = t544 / 0.2e1;
t9 = t547 / 0.2e1;
t8 = t526 + t527 - t519;
t7 = t519 - t578;
t6 = t519 + t578;
t5 = t23 + t9 - t544 / 0.2e1;
t4 = t18 + t23 - t547 / 0.2e1;
t3 = t18 + t9 - t543 / 0.2e1;
t2 = t427 * t345 + t428 * t347 + t542 + t550;
t1 = t551 + t545 + (t426 * t345 + t425 * t347 + t429) * t350 + (t431 * t345 + t430 * t347 - t424) * t349;
t10 = [0, t169 * qJD(3) + t30 * qJD(4) + t172 * qJD(5) + (-m(3) * t227 / 0.2e1 + t133 * t558 + t111 * t557 + t92 * t555) * t562, t169 * qJD(2), t30 * qJD(2) + (t106 * t555 + t134 * t557) * t560 + t168 * qJD(5), qJD(2) * t172 + qJD(4) * t168; qJD(3) * t170 - qJD(4) * t29 - qJD(5) * t171, t22 * qJD(3) + t2 * qJD(4) + t44 * t523 + (m(6) * (t144 * t145 + t146 * t147 + t58 * t92) + m(5) * (t111 * t93 + t163 * t164 + t165 * t166) + m(4) * (t127 * t133 + t228 * t229 + t230 * t231) + m(3) * (-t227 + t326) * t436 * (rSges(3,1) * t350 - t349 * rSges(3,2)) + (t340 * t317 + (t314 * t254 + t315 * t256) * t345 + t38 + t37 + ((-t316 + t368) * t345 - t314 * t253 - t315 * t255 + t565 + t369 * t347) * t347) * t537 + (t341 * t316 - (t312 * t253 + t313 * t255) * t347 + t35 + t36 + ((-t317 + t369) * t347 + t312 * t254 + t313 * t256 + t565 + t368 * t345) * t345) * t536) * qJD(2), t22 * qJD(2) + t6 * qJD(4) + t90 * qJD(5) + t434 + (-0.4e1 * t570 + 0.2e1 * t416 * (-t322 * t350 + t437)) * qJD(3), -t494 + t2 * qJD(2) + t6 * qJD(3) + t4 * qJD(5) + (-t545 / 0.4e1 - t551 / 0.4e1) * t559 + ((t132 * t138 + t134 * t93 + t160 * t165 + t161 * t163 + (-t151 * t347 - t152 * t345) * t302) * t557 + (t106 * t58 + t110 * t95 - t118 * t196 - t119 * t195 + t124 * t146 + t125 * t144) * t555) * t560 + (t345 * t383 + t347 * t382 + t424) * t435 + (((t85 / 0.2e1 + t87 / 0.2e1 - t425) * t347 + (-t86 / 0.2e1 - t88 / 0.2e1 - t426) * t345 - t429) * t350 + t572 * t537 + t573 * t536) * qJD(4), t90 * qJD(3) + t4 * qJD(4) + t44 * t525 - t569; -t170 * qJD(2), -t434 + t150 * qJD(3) + t7 * qJD(4) + t89 * qJD(5) + (-t540 / 0.4e1 - t549 / 0.4e1 - t553 / 0.4e1) * t561 + ((-t350 * t92 + t469) * t555 + (-t111 * t350 + t467) * t557 + (-t133 * t350 + t454) * t558 + ((t145 * t345 + t147 * t347 + t58) * t555 + (t164 * t345 + t166 * t347 + t93) * t557 + (t229 * t345 + t231 * t347 + t127) * t558) * t349) * t562, t150 * qJD(2), t7 * qJD(2) + ((-t134 * t350 + (t160 * t347 + t161 * t345) * t349) * t557 + (-t106 * t350 + (t124 * t347 + t125 * t345) * t349) * t555) * t560 + t96 * qJD(5), qJD(2) * t89 + qJD(4) * t96; qJD(2) * t29 - qJD(5) * t167, t494 + t8 * qJD(3) + t1 * qJD(4) + t3 * qJD(5) + (-t542 / 0.4e1 - t550 / 0.4e1) * t561 + ((t109 * t93 + t111 * t132 + t122 * t165 + t123 * t163 + t151 * t166 + t152 * t164) * t557 + (t118 * t147 + t119 * t145 + t144 * t84 + t146 * t83 + t57 * t58 + t92 * t95) * t555) * t562 + (t383 * t347 - t382 * t345 + (t595 * t577 + (t592 * t345 - t593 * t347) * t537 + (t590 * t345 - t591 * t347 + t596) * t535) * t350 + ((-t114 / 0.2e1 - t112 / 0.2e1 + t38 / 0.2e1 + t37 / 0.2e1) * t347 + (t115 / 0.2e1 + t113 / 0.2e1 + t36 / 0.2e1 + t35 / 0.2e1) * t345) * t349) * qJD(2), qJD(2) * t8 + qJD(5) * t97, t1 * qJD(2) + ((t106 * t95 + t118 * t124 + t119 * t125) * t554 + (t132 * t134 + t151 * t160 + t152 * t161) * t556) * t559 + t43 * t523 + (-t295 / 0.2e1 - t296 / 0.2e1) * qJD(4) * t350 * t563 + (t573 * t537 + t572 * t535 + ((t88 + t86) * t347 + (t87 + t85) * t345 + t600 * t571 + t601 * t472) * t575) * t435, t3 * qJD(2) + t43 * t524 + (t174 * t490 + t232 * t290 + t233 * t288 - t140) * t523 + t384; qJD(2) * t171 + qJD(4) * t167, (t288 * t145 + t290 * t147 - t44 + (t350 * t58 + (t393 + t92) * t349) * t338) * t525 + t91 * qJD(3) + t5 * qJD(4) + t569, qJD(2) * t91 - qJD(4) * t97, t5 * qJD(2) + (t118 * t291 + t119 * t289 + t124 * t290 + t125 * t288 + (t106 * t338 + t339 * t95) * t349 - t43) * t524 + t140 * t523 - t384, 0.4e1 * (t159 * qJD(2) / 0.4e1 + t140 * qJD(4) / 0.4e1) * m(6);];
Cq = t10;
