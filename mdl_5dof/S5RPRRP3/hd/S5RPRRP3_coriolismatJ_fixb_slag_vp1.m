% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:11
% EndTime: 2019-12-05 18:03:31
% DurationCPUTime: 8.56s
% Computational Cost: add. (29351->389), mult. (22022->503), div. (0->0), fcn. (20324->8), ass. (0->247)
t595 = Icges(5,4) + Icges(6,4);
t592 = Icges(5,5) + Icges(6,5);
t588 = Icges(5,2) + Icges(6,2);
t591 = Icges(5,6) + Icges(6,6);
t358 = qJ(3) + qJ(4);
t353 = cos(t358);
t596 = t595 * t353;
t593 = Icges(5,1) + Icges(6,1);
t352 = sin(t358);
t594 = -t591 * t352 + t592 * t353;
t570 = t588 * t352 - t596;
t590 = Icges(5,3) + Icges(6,3);
t357 = qJ(1) + pkin(8);
t351 = cos(t357);
t443 = t351 * t352;
t589 = t595 * t443;
t350 = sin(t357);
t584 = t570 * t350 + t591 * t351;
t442 = t351 * t353;
t583 = t591 * t350 + t595 * t442 - t588 * t443;
t581 = t592 * t350 + t593 * t442 - t589;
t449 = t350 * t352;
t587 = t595 * t449;
t586 = t594 * t350;
t569 = t593 * t352 + t596;
t585 = t590 * t350 + t592 * t442 - t591 * t443;
t448 = t350 * t353;
t582 = -t592 * t351 + t593 * t448 - t587;
t580 = t595 * t352;
t578 = t588 * t353 + t580;
t577 = t593 * t353 - t580;
t576 = t592 * t352 + t591 * t353;
t575 = -t588 * t442 + t581 - t589;
t574 = t588 * t448 - t582 + t587;
t573 = t569 * t351 + t583;
t572 = t569 * t350 - t584;
t571 = t584 * t449 + (t590 * t351 - t586) * t351;
t568 = (t583 * t352 - t581 * t353) * t351;
t567 = t582 * t448 + t571;
t566 = t585 * t351 - t581 * t448 + t583 * t449;
t565 = t582 * t353 - t585;
t564 = t584 * t352;
t398 = rSges(6,1) * t443 + rSges(6,2) * t442;
t233 = t351 * t398;
t255 = rSges(6,1) * t449 + rSges(6,2) * t448;
t348 = t350 ^ 2;
t349 = t351 ^ 2;
t396 = t348 + t349;
t486 = pkin(4) * t352;
t126 = -t255 * t350 - t396 * t486 - t233;
t297 = rSges(6,1) * t352 + rSges(6,2) * t353;
t321 = pkin(4) * t449;
t209 = t350 * t297 + t321;
t359 = sin(qJ(3));
t447 = t350 * t359;
t338 = pkin(3) * t447;
t178 = t338 + t209;
t487 = pkin(3) * t359;
t320 = -t486 - t487;
t180 = (t297 - t320) * t351;
t476 = rSges(6,1) * t353;
t299 = -rSges(6,2) * t352 + t476;
t210 = pkin(4) * t448 + t350 * t299;
t485 = pkin(4) * t353;
t393 = -t299 - t485;
t212 = t393 * t351;
t426 = t178 * t210 - t180 * t212;
t363 = -pkin(7) - pkin(6);
t340 = t350 * t363;
t361 = cos(qJ(3));
t355 = t361 * pkin(3);
t347 = t355 + pkin(2);
t483 = pkin(2) - t347;
t184 = t351 * (pkin(6) * t350 + t483 * t351 + t340);
t344 = t351 * pkin(6);
t439 = t351 * t363;
t199 = t483 * t350 - t344 - t439;
t301 = t347 + t485;
t356 = -qJ(5) + t363;
t562 = -rSges(6,1) * t442 + rSges(6,2) * t443 - t351 * t301 - (rSges(6,3) - t356) * t350;
t419 = (-t347 * t351 + t340 - t562) * t351;
t400 = -rSges(6,2) * t449 - t351 * rSges(6,3);
t424 = -(-t356 + t363) * t351 - (-t301 + t347) * t350 + rSges(6,1) * t448 + t400;
t55 = -t184 + (-t199 + t424) * t350 + t419;
t32 = t55 * t126 + t426;
t256 = rSges(5,1) * t449 + rSges(5,2) * t448;
t298 = rSges(5,1) * t352 + rSges(5,2) * t353;
t257 = t298 * t351;
t160 = -t256 * t350 - t351 * t257;
t229 = t298 * t350 + t338;
t382 = (t298 + t487) * t351;
t477 = rSges(5,1) * t353;
t300 = -rSges(5,2) * t352 + t477;
t451 = t300 * t351;
t452 = t300 * t350;
t390 = rSges(5,2) * t443 - rSges(5,3) * t350;
t193 = t351 * (rSges(5,1) * t442 - t390);
t399 = -rSges(5,2) * t449 - t351 * rSges(5,3);
t226 = -rSges(5,1) * t448 - t399;
t89 = -t184 + t193 + (-t199 - t226) * t350;
t44 = t89 * t160 + t229 * t452 + t382 * t451;
t563 = -m(5) * t44 - m(6) * t32;
t561 = t584 * t443;
t560 = t585 * t350 - t568 + t571;
t518 = m(5) / 0.2e1;
t517 = m(6) / 0.2e1;
t484 = sin(qJ(1)) * pkin(1);
t158 = t484 + t351 * t356 + (t301 + t476) * t350 + t400;
t488 = cos(qJ(1)) * pkin(1);
t159 = -t488 + t562;
t559 = m(6) * (-t351 * t158 - t159 * t350);
t558 = -t350 / 0.2e1;
t499 = t350 / 0.2e1;
t497 = t351 / 0.2e1;
t553 = -t361 / 0.2e1;
t551 = (Icges(4,1) * t361 / 0.2e1 - Icges(4,4) * t359 + Icges(4,2) * t553) * t359;
t550 = t576 * t350;
t549 = t576 * t351;
t182 = -t320 * t350 + t255;
t183 = -t320 * t351 + t398;
t201 = t338 + t256;
t211 = (t297 + t486) * t351;
t391 = t347 + t477;
t161 = t391 * t350 + t399 + t439 + t484;
t162 = -t391 * t351 + t340 + t390 - t488;
t383 = t161 * t451 + t162 * t452;
t427 = -t212 * t158 + t210 * t159;
t481 = (-t182 * t211 + t183 * t209 + t427) * t517 + ((-t201 * t351 + t350 * t382) * t298 + t383) * t518;
t197 = t321 + t255;
t198 = pkin(4) * t443 + t398;
t482 = (t178 * t198 - t180 * t197 + t427) * t517 + (t229 * t257 - t256 * t382 + t383) * t518;
t547 = t481 - t482;
t545 = -t574 * t352 + t572 * t353;
t544 = t575 * t352 + t573 * t353;
t543 = (-t577 + t578) * t353 + (t569 - t570) * t352;
t542 = qJD(1) * t559;
t440 = t351 * t361;
t441 = t351 * t359;
t238 = Icges(4,4) * t440 - Icges(4,2) * t441 + Icges(4,6) * t350;
t335 = Icges(4,4) * t441;
t240 = Icges(4,1) * t440 + Icges(4,5) * t350 - t335;
t541 = (t238 * t359 - t240 * t361) * t351;
t354 = Icges(4,4) * t361;
t536 = Icges(4,1) * t359 + t354;
t533 = Icges(4,2) * t359 - t354;
t405 = -Icges(4,2) * t440 + t240 - t335;
t407 = t351 * t536 + t238;
t525 = t359 * t405 + t361 * t407;
t334 = Icges(4,4) * t447;
t446 = t350 * t361;
t239 = -Icges(4,1) * t446 + Icges(4,5) * t351 + t334;
t406 = Icges(4,2) * t446 + t239 + t334;
t237 = Icges(4,6) * t351 + t533 * t350;
t408 = -t350 * t536 + t237;
t524 = -t359 * t406 - t361 * t408;
t367 = (-t570 / 0.2e1 + t569 / 0.2e1) * t353 + (-t578 / 0.2e1 + t577 / 0.2e1) * t352;
t368 = (t566 * t350 + t567 * t351) * t558 + ((t560 + t568) * t351 + ((-t564 + t565) * t351 + t561 + t566) * t350) * t499 + (((t564 - t585) * t351 - t561 + t566) * t351 + (t565 * t350 + t560 - t567) * t350) * t497;
t523 = 0.4e1 * qJD(1);
t522 = 2 * qJD(3);
t520 = 2 * qJD(4);
t519 = m(4) / 0.2e1;
t478 = rSges(4,1) * t361;
t394 = pkin(2) + t478;
t397 = -rSges(4,2) * t447 - t351 * rSges(4,3);
t170 = t394 * t350 - t344 + t397 + t484;
t337 = rSges(4,2) * t441;
t171 = -t488 + t337 - t394 * t351 + (-rSges(4,3) - pkin(6)) * t350;
t328 = rSges(4,1) * t359 + rSges(4,2) * t361;
t280 = t328 * t350;
t281 = t328 * t351;
t515 = m(4) * (-t170 * t280 + t171 * t281);
t146 = -t226 * t350 + t193;
t80 = t396 * t298 * t300 + t146 * t160;
t79 = m(5) * t80;
t509 = m(5) * (-t161 * t201 + t162 * t382);
t508 = m(5) * (-t161 * t256 + t162 * t257);
t425 = t209 * t210 - t211 * t212;
t84 = t424 * t350 + t419;
t99 = t349 * (t320 + t487) - t233 + (t338 - t182) * t350;
t504 = m(6) * (t84 * t99 + t425);
t501 = m(6) * (-t158 * t182 + t159 * t183);
t500 = m(6) * (-t158 * t197 + t159 * t198);
t498 = -t351 / 0.2e1;
t494 = m(6) * (t182 * t350 + t183 * t351);
t493 = m(6) * (-t178 * t350 - t351 * t180);
t492 = m(6) * t126;
t491 = m(6) * (t197 * t350 + t198 * t351);
t489 = m(6) * (-t209 * t350 - t351 * t211);
t479 = m(6) * qJD(3);
t430 = t359 * t237;
t76 = 0.2e1 * (t126 / 0.4e1 - t99 / 0.4e1) * m(6);
t428 = t76 * qJD(2);
t378 = Icges(4,5) * t361 - Icges(4,6) * t359;
t235 = Icges(4,3) * t351 - t378 * t350;
t418 = t351 * t235 + t350 * t430;
t417 = t350 * t235 + t239 * t440;
t395 = t84 * t126 + t425;
t392 = (-t549 * t348 + (t545 * t351 + (-t544 + t550) * t350) * t351) * t499 + (t550 * t349 + (t544 * t350 + (-t545 - t549) * t351) * t350) * t497;
t236 = Icges(4,5) * t440 - Icges(4,6) * t441 + Icges(4,3) * t350;
t387 = -t239 * t361 - t236;
t386 = t396 * t487;
t385 = t79 + t392;
t377 = Icges(4,5) * t359 + Icges(4,6) * t361;
t113 = t351 * t236 + t238 * t447 - t240 * t446;
t365 = -t368 + (-t543 * t351 - t573 * t352 + t575 * t353 + t586) * t499 + (t543 * t350 + t594 * t351 + t572 * t352 + t574 * t353) * t497;
t364 = -t367 + (t581 * t352 + t583 * t353) * (t497 + t498);
t339 = pkin(3) * t446;
t330 = -rSges(4,2) * t359 + t478;
t275 = t377 * t351;
t274 = t350 * t377;
t232 = (-t300 - t355) * t351;
t230 = t339 + t452;
t181 = (t393 - t355) * t351;
t179 = t339 + t210;
t165 = -t280 * t350 - t281 * t351;
t155 = m(5) * t160;
t141 = t210 * t351 + t212 * t350;
t135 = t160 - t386;
t134 = t489 / 0.2e1;
t128 = m(6) * t141 * qJD(4);
t125 = t491 / 0.2e1;
t118 = t493 / 0.2e1;
t117 = t494 / 0.2e1;
t115 = t350 * t236 - t541;
t114 = -t351 * t430 + t417;
t112 = -t239 * t446 + t418;
t90 = -t386 + t99;
t75 = t114 * t351 + t115 * t350;
t74 = t112 * t351 + t113 * t350;
t73 = t134 - t491 / 0.2e1;
t72 = t134 + t125;
t71 = t125 - t489 / 0.2e1;
t48 = t155 + t492 / 0.2e1 + t99 * t517;
t47 = t118 + t117;
t46 = t118 - t494 / 0.2e1;
t45 = t117 - t493 / 0.2e1;
t28 = t367 + t500 + t508;
t25 = (t113 + (-t236 + t430) * t351 - t417) * t351 + (t387 * t350 - t112 + t418) * t350;
t24 = (t115 + t418 + t541) * t351 + (-t114 + (t387 - t430) * t351 + t113 + t417) * t350;
t18 = (t536 / 0.2e1 - t533 / 0.2e1) * t361 + t551 + t515 + t509 + t501 + t367;
t7 = t385 + t504;
t6 = t392 - t563;
t4 = t368 + t547;
t3 = t368 - t547;
t2 = (t75 / 0.2e1 + t25 / 0.2e1) * t351 + (t24 / 0.2e1 - t74 / 0.2e1) * t350 + t368;
t1 = t365 + t481 + t482;
t5 = [t18 * qJD(3) + t28 * qJD(4) + qJD(5) * t559, 0, t18 * qJD(1) + t1 * qJD(4) + t47 * qJD(5) + (((t170 * t351 + t171 * t350) * t330 + (-t280 * t351 + t281 * t350) * t328) * t519 + (-t161 * t232 + t162 * t230 + (-t201 + t229) * t382) * t518 + (-t158 * t181 + t159 * t179 + t178 * t183 - t180 * t182) * t517) * t522 + ((-t359 * t408 + t361 * t406) * t497 + t365 + t24 * t558 + (-t359 * t407 + t361 * t405 + t74) * t499 + (t75 + t25) * t498 + (t349 / 0.2e1 + t348 / 0.2e1) * t378) * qJD(3), t28 * qJD(1) + t1 * qJD(3) + t365 * qJD(4) + t72 * qJD(5) + (((-t256 * t351 + t257 * t350) * t298 + t383) * t518 + (-t197 * t211 + t198 * t209 + t427) * t517) * t520, t47 * qJD(3) + t72 * qJD(4) + t542; 0, 0, t48 * qJD(4) + (t135 * t518 + t165 * t519 + t90 * t517) * t522, t48 * qJD(3) + (t155 + t492) * qJD(4), 0; (t364 + (-t533 + t536) * t553 - t551) * qJD(1) + t2 * qJD(3) + t3 * qJD(4) + t46 * qJD(5) + (-t515 / 0.4e1 - t501 / 0.4e1 - t509 / 0.4e1) * t523, qJD(4) * t76, t2 * qJD(1) + (m(6) * (t178 * t179 - t180 * t181 + t55 * t90) + m(5) * (t135 * t89 + t229 * t230 - t232 * t382) + m(4) * (t328 * t330 * t396 + (t351 * (rSges(4,1) * t440 + rSges(4,3) * t350 - t337) - t350 * (-rSges(4,1) * t446 - t397)) * t165) + (-t348 * t275 + (t524 * t351 + (t274 - t525) * t350) * t351) * t499 + (t349 * t274 + (t525 * t350 + (-t275 - t524) * t351) * t350) * t497 + t392) * qJD(3) + t6 * qJD(4), t3 * qJD(1) + t428 + t6 * qJD(3) + ((t395 + t32) * t517 + (t44 + t80) * t518) * t520 + (t392 - t79 - t504) * qJD(4), t46 * qJD(1); t364 * qJD(1) + t4 * qJD(3) + t368 * qJD(4) + t73 * qJD(5) + (-t508 / 0.4e1 - t500 / 0.4e1) * t523, -qJD(3) * t76, t4 * qJD(1) - t428 + t7 * qJD(4) + ((t179 * t209 - t181 * t211 + t55 * t99 + t84 * t90 + t426) * t517 + (t135 * t146 + (t230 * t350 - t232 * t351) * t298 + t44) * t518) * t522 + (t392 + t563) * qJD(3), t368 * qJD(1) + t7 * qJD(3) + (m(6) * t395 + t385) * qJD(4), t73 * qJD(1); t45 * qJD(3) + t71 * qJD(4) - t542, 0, t45 * qJD(1) + (t179 * t351 + t181 * t350) * t479 + t128, t71 * qJD(1) + t141 * t479 + t128, 0;];
Cq = t5;
