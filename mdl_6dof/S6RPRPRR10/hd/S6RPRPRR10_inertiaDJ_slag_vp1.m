% Calculate time derivative of joint inertia matrix for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR10_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:40
% EndTime: 2019-03-09 04:08:12
% DurationCPUTime: 19.95s
% Computational Cost: add. (38853->1048), mult. (45031->1455), div. (0->0), fcn. (42890->10), ass. (0->493)
t642 = Icges(5,5) / 0.2e1;
t641 = Icges(5,6) / 0.2e1;
t640 = Icges(5,3) / 0.2e1;
t382 = sin(qJ(1));
t604 = -t382 / 0.2e1;
t639 = -qJD(1) / 0.2e1;
t638 = rSges(5,3) + qJ(4);
t374 = pkin(10) + qJ(5);
t363 = cos(t374);
t381 = sin(qJ(3));
t384 = cos(qJ(1));
t362 = sin(t374);
t560 = t382 * t362;
t289 = t363 * t384 - t381 * t560;
t559 = t382 * t363;
t290 = t362 * t384 + t381 * t559;
t383 = cos(qJ(3));
t556 = t382 * t383;
t210 = t290 * rSges(6,1) + t289 * rSges(6,2) - rSges(6,3) * t556;
t378 = sin(pkin(10));
t567 = t378 * t384;
t353 = pkin(4) * t567;
t379 = cos(pkin(10));
t360 = t379 * pkin(4) + pkin(3);
t380 = -pkin(8) - qJ(4);
t565 = t381 * t382;
t499 = t360 * t565 + t380 * t556 + t353;
t503 = -t210 - t499;
t524 = qJD(3) * t381;
t481 = -t524 / 0.2e1;
t526 = qJD(1) * t383;
t482 = -t526 / 0.2e1;
t637 = t382 * t482 + t384 * t481;
t480 = t524 / 0.2e1;
t636 = t382 * t480 + t384 * t482;
t522 = qJD(3) * t383;
t488 = t382 * t522;
t525 = qJD(1) * t384;
t635 = t381 * t525 + t488;
t521 = qJD(3) * t384;
t489 = t381 * t521;
t492 = t382 * t526;
t634 = t489 + t492;
t523 = qJD(3) * t382;
t490 = t381 * t523;
t491 = t383 * t525;
t392 = t490 - t491;
t585 = Icges(6,4) * t363;
t428 = -Icges(6,2) * t362 + t585;
t270 = Icges(6,6) * t381 + t383 * t428;
t586 = Icges(6,4) * t362;
t433 = Icges(6,1) * t363 - t586;
t271 = Icges(6,5) * t381 + t383 * t433;
t633 = -t270 * t362 + t271 * t363;
t364 = qJ(6) + t374;
t358 = sin(t364);
t359 = cos(t364);
t583 = Icges(7,4) * t359;
t427 = -Icges(7,2) * t358 + t583;
t262 = Icges(7,6) * t381 + t383 * t427;
t584 = Icges(7,4) * t358;
t432 = Icges(7,1) * t359 - t584;
t263 = Icges(7,5) * t381 + t383 * t432;
t632 = -t262 * t358 + t263 * t359;
t373 = -pkin(9) + t380;
t531 = t373 - t380;
t631 = t381 * t531;
t597 = pkin(3) - t360;
t630 = t381 * t597;
t458 = rSges(4,1) * t381 + rSges(4,2) * t383;
t400 = t384 * t458;
t588 = Icges(4,4) * t381;
t430 = Icges(4,2) * t383 + t588;
t295 = Icges(4,6) * t384 + t382 * t430;
t587 = Icges(4,4) * t383;
t435 = Icges(4,1) * t381 + t587;
t297 = Icges(4,5) * t384 + t382 * t435;
t411 = t295 * t383 + t297 * t381;
t398 = t411 * t384;
t322 = pkin(5) * t363 + t360;
t535 = t322 - t360;
t235 = t383 * t535 - t631;
t450 = rSges(7,1) * t359 - rSges(7,2) * t358;
t264 = rSges(7,3) * t381 + t383 * t450;
t542 = t235 + t264;
t475 = t542 * t382;
t558 = t382 * t378;
t308 = t379 * t384 - t381 * t558;
t557 = t382 * t379;
t309 = t381 * t557 + t567;
t356 = pkin(3) * t565;
t629 = -t309 * rSges(5,1) - t308 * rSges(5,2) - t356;
t496 = rSges(4,1) * t635 + rSges(4,2) * t491;
t610 = -pkin(1) - pkin(7);
t518 = -rSges(4,3) + t610;
t533 = qJ(2) * t525 + qJD(2) * t382;
t179 = (-rSges(4,2) * t524 + qJD(1) * t518) * t382 + t496 + t533;
t596 = rSges(4,2) * t381;
t338 = rSges(4,1) * t383 - t596;
t366 = qJD(2) * t384;
t180 = t366 + t338 * t521 + (t518 * t384 + (-qJ(2) - t458) * t382) * qJD(1);
t628 = -t179 * t384 + t382 * t180;
t528 = qJD(1) * t381;
t468 = qJD(5) + t528;
t487 = t383 * t521;
t627 = t382 * t468 - t487;
t375 = qJD(5) + qJD(6);
t472 = t375 + t528;
t626 = t382 * t472 - t487;
t571 = t322 * t381;
t591 = rSges(7,3) - t373;
t625 = -t383 * t591 + t571;
t425 = Icges(4,5) * t381 + Icges(4,6) * t383;
t624 = -Icges(4,3) * t382 + t384 * t425;
t623 = -Icges(4,6) * t382 + t384 * t430;
t622 = -Icges(4,5) * t382 + t384 * t435;
t621 = t384 * t468 + t488;
t620 = t384 * t472 + t488;
t562 = t382 * t358;
t276 = t359 * t384 - t381 * t562;
t561 = t382 * t359;
t277 = t358 * t384 + t381 * t561;
t200 = t277 * rSges(7,1) + t276 * rSges(7,2) - rSges(7,3) * t556;
t143 = t381 * t200 + t264 * t556;
t564 = t381 * t384;
t278 = t358 * t564 + t561;
t279 = -t359 * t564 + t562;
t451 = -t279 * rSges(7,1) - t278 * rSges(7,2);
t555 = t383 * t384;
t201 = rSges(7,3) * t555 - t451;
t245 = t264 * t555;
t144 = -t381 * t201 + t245;
t473 = t375 * t381 + qJD(1);
t408 = t384 * t473;
t169 = -t358 * t626 + t359 * t408;
t170 = t358 * t408 + t359 * t626;
t452 = t170 * rSges(7,1) + t169 * rSges(7,2);
t109 = -rSges(7,3) * t634 + t452;
t569 = t375 * t383;
t188 = (-rSges(7,1) * t358 - rSges(7,2) * t359) * t569 + (rSges(7,3) * t383 - t381 * t450) * qJD(3);
t174 = t188 * t555;
t63 = -t264 * t492 - t381 * t109 + t174 + (-t201 * t383 - t264 * t564) * qJD(3);
t409 = t382 * t473;
t171 = -t358 * t620 - t359 * t409;
t172 = -t358 * t409 + t359 * t620;
t110 = t172 * rSges(7,1) + t171 * rSges(7,2) + rSges(7,3) * t392;
t467 = t381 * t110 + t188 * t556 + t200 * t522 + t264 * t491;
t64 = -t264 * t490 + t467;
t619 = qJD(1) * (t143 * t382 + t144 * t384) + t63 * t382 - t384 * t64;
t618 = t381 * t535 + t383 * t531;
t617 = 2 * m(4);
t616 = 2 * m(5);
t615 = 2 * m(6);
t614 = 2 * m(7);
t376 = t382 ^ 2;
t377 = t384 ^ 2;
t613 = m(5) / 0.2e1;
t612 = m(6) / 0.2e1;
t611 = m(7) / 0.2e1;
t429 = Icges(5,4) * t379 - Icges(5,2) * t378;
t283 = Icges(5,6) * t381 + t383 * t429;
t609 = t283 / 0.2e1;
t434 = Icges(5,1) * t379 - Icges(5,4) * t378;
t284 = Icges(5,5) * t381 + t383 * t434;
t608 = t284 / 0.2e1;
t310 = t378 * t564 + t557;
t607 = t310 / 0.2e1;
t404 = t379 * t564 - t558;
t606 = -t404 / 0.2e1;
t605 = t381 / 0.2e1;
t603 = t382 / 0.2e1;
t602 = t384 / 0.2e1;
t601 = rSges(3,2) - pkin(1);
t600 = m(4) * t338;
t599 = pkin(3) * t381;
t598 = pkin(4) * t378;
t595 = rSges(5,3) * t383;
t594 = pkin(5) * qJD(5);
t593 = t382 * rSges(4,3);
t370 = t384 * rSges(4,3);
t592 = rSges(6,3) - t380;
t422 = Icges(7,5) * t359 - Icges(7,6) * t358;
t261 = Icges(7,3) * t381 + t383 * t422;
t137 = t261 * t381 + t383 * t632;
t183 = (-Icges(7,2) * t359 - t584) * t569 + (Icges(7,6) * t383 - t381 * t427) * qJD(3);
t182 = (-Icges(7,5) * t358 - Icges(7,6) * t359) * t569 + (Icges(7,3) * t383 - t381 * t422) * qJD(3);
t184 = (-Icges(7,1) * t358 - t583) * t569 + (Icges(7,5) * t383 - t381 * t432) * qJD(3);
t389 = t383 * t359 * t184 + t381 * t182 + t261 * t522 - t524 * t632;
t507 = t262 * t359 * t375;
t590 = t137 * t522 + ((-t507 + (-t263 * t375 - t183) * t358) * t383 + t389) * t381;
t423 = Icges(6,5) * t363 - Icges(6,6) * t362;
t269 = Icges(6,3) * t381 + t383 * t423;
t140 = t269 * t381 + t383 * t633;
t519 = qJD(5) * t383;
t212 = (-Icges(6,5) * t362 - Icges(6,6) * t363) * t519 + (Icges(6,3) * t383 - t381 * t423) * qJD(3);
t214 = (-Icges(6,1) * t362 - t585) * t519 + (Icges(6,5) * t383 - t381 * t433) * qJD(3);
t388 = t383 * t363 * t214 + t381 * t212 + t269 * t522 - t524 * t633;
t573 = t270 * t363;
t213 = (-Icges(6,2) * t363 - t586) * t519 + (Icges(6,6) * t383 - t381 * t428) * qJD(3);
t577 = t213 * t362;
t589 = t140 * t522 + ((-t577 + (-t271 * t362 - t573) * qJD(5)) * t383 + t388) * t381;
t291 = t362 * t564 + t559;
t292 = -t363 * t564 + t560;
t211 = t292 * rSges(6,1) + t291 * rSges(6,2) + rSges(6,3) * t555;
t578 = t211 * t384;
t570 = t322 * t383;
t568 = t378 * t381;
t566 = t379 * t381;
t554 = -qJ(4) - t380;
t319 = t360 * t487;
t516 = t362 * t594;
t470 = t381 * t516;
t329 = pkin(5) * t362 + t598;
t479 = -t329 + t598;
t515 = t363 * t594;
t553 = t382 * t515 + t319 + (t470 + (-t570 + t631) * qJD(3)) * t384 + (t382 * t618 - t479 * t384) * qJD(1) + t109;
t464 = t322 * t635 + t373 * t491 + t384 * t515;
t500 = t360 * t635 + t380 * t491;
t101 = (t479 * qJD(1) + (-qJD(3) * t531 - t516) * t381) * t382 + t464 - t500;
t552 = -t101 - t110;
t497 = pkin(3) * t487 + qJ(4) * t634;
t520 = qJD(4) * t383;
t217 = t384 * (qJD(1) * t356 + t384 * t520 - t497);
t551 = t384 * (t380 * t489 - t319 + (t353 + (t380 * t383 - t630) * t382) * qJD(1) + t497) + t217;
t498 = pkin(3) * t635 + qJ(4) * t490;
t228 = (-qJ(4) * t525 - qJD(4) * t382) * t383 + t498;
t512 = qJ(4) * t555;
t394 = -pkin(4) * t558 + t512;
t550 = -qJD(1) * t394 + t380 * t490 - t228 + t498 - t500;
t501 = t322 * t565 + t384 * t329 + t373 * t556;
t160 = -t499 + t501;
t549 = -t160 - t200;
t534 = t360 * t564 + t380 * t555;
t161 = (-t373 * t383 - t571) * t384 - t479 * t382 + t534;
t548 = t161 + t201;
t486 = t362 * t519;
t202 = -pkin(5) * t486 - qJD(3) * t618;
t547 = t188 + t202;
t357 = pkin(3) * t564;
t227 = t357 - t394 - t534;
t315 = -t357 + t512;
t304 = t384 * t315;
t546 = t384 * t227 + t304;
t545 = t556 * t638 + t629;
t543 = -t227 - t315;
t249 = (t383 * t554 + t630) * qJD(3);
t307 = qJD(4) * t381 + (qJ(4) * t383 - t599) * qJD(3);
t541 = -t249 - t307;
t265 = t381 * t554 - t383 * t597;
t337 = pkin(3) * t383 + qJ(4) * t381;
t527 = qJD(1) * t382;
t316 = t337 * t527;
t540 = t265 * t527 + t316;
t321 = t382 * t337;
t539 = t382 * t265 + t321;
t538 = -t265 - t337;
t537 = t382 * t307 + t337 * t525;
t532 = t384 * pkin(1) + t382 * qJ(2);
t530 = t376 + t377;
t293 = Icges(4,3) * t384 + t382 * t425;
t529 = qJD(1) * t293;
t517 = -t329 + t610;
t469 = qJD(5) * t381 + qJD(1);
t407 = t382 * t469;
t192 = -t362 * t621 - t363 * t407;
t193 = -t362 * t407 + t363 * t621;
t118 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t392;
t120 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t392;
t122 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t392;
t204 = Icges(6,5) * t290 + Icges(6,6) * t289 - Icges(6,3) * t556;
t206 = Icges(6,4) * t290 + Icges(6,2) * t289 - Icges(6,6) * t556;
t208 = Icges(6,1) * t290 + Icges(6,4) * t289 - Icges(6,5) * t556;
t416 = t206 * t362 - t208 * t363;
t33 = (qJD(3) * t416 + t118) * t381 + (qJD(3) * t204 - t120 * t362 + t122 * t363 + (-t206 * t363 - t208 * t362) * qJD(5)) * t383;
t56 = t192 * t270 + t193 * t271 - t212 * t556 + t289 * t213 + t290 * t214 + t269 * t392;
t514 = -t33 / 0.2e1 - t56 / 0.2e1;
t406 = t384 * t469;
t190 = -t362 * t627 + t363 * t406;
t191 = t362 * t406 + t363 * t627;
t117 = Icges(6,5) * t191 + Icges(6,6) * t190 - Icges(6,3) * t634;
t119 = Icges(6,4) * t191 + Icges(6,2) * t190 - Icges(6,6) * t634;
t121 = Icges(6,1) * t191 + Icges(6,4) * t190 - Icges(6,5) * t634;
t205 = Icges(6,5) * t292 + Icges(6,6) * t291 + Icges(6,3) * t555;
t207 = Icges(6,4) * t292 + Icges(6,2) * t291 + Icges(6,6) * t555;
t209 = Icges(6,1) * t292 + Icges(6,4) * t291 + Icges(6,5) * t555;
t415 = t207 * t362 - t209 * t363;
t34 = (qJD(3) * t415 + t117) * t381 + (qJD(3) * t205 - t119 * t362 + t121 * t363 + (-t207 * t363 - t209 * t362) * qJD(5)) * t383;
t55 = t190 * t270 + t191 * t271 + t212 * t555 + t291 * t213 + t292 * t214 - t269 * t634;
t513 = t34 / 0.2e1 + t55 / 0.2e1;
t218 = Icges(5,5) * t309 + Icges(5,6) * t308 - Icges(5,3) * t556;
t511 = t218 * t556;
t510 = t218 * t555;
t219 = -Icges(5,5) * t404 + Icges(5,6) * t310 + Icges(5,3) * t555;
t509 = t219 * t556;
t508 = t219 * t555;
t129 = -t269 * t556 + t270 * t289 + t271 * t290;
t94 = t204 * t381 - t383 * t416;
t506 = t94 / 0.2e1 + t129 / 0.2e1;
t130 = t269 * t555 + t291 * t270 + t292 * t271;
t95 = t205 * t381 - t383 * t415;
t505 = -t95 / 0.2e1 - t130 / 0.2e1;
t504 = t200 * t634 + t201 * t490;
t254 = -qJD(1) * t310 - t378 * t488;
t255 = qJD(1) * t404 + t379 * t488;
t502 = -t255 * rSges(5,1) - t254 * rSges(5,2) - rSges(5,3) * t490;
t299 = rSges(4,1) * t565 + rSges(4,2) * t556 + t370;
t495 = t384 * pkin(7) + t532;
t453 = rSges(6,1) * t363 - rSges(6,2) * t362;
t275 = rSges(6,3) * t381 + t383 * t453;
t493 = t275 * t527;
t485 = -t556 / 0.2e1;
t484 = t555 / 0.2e1;
t104 = Icges(7,5) * t172 + Icges(7,6) * t171 + Icges(7,3) * t392;
t106 = Icges(7,4) * t172 + Icges(7,2) * t171 + Icges(7,6) * t392;
t108 = Icges(7,1) * t172 + Icges(7,4) * t171 + Icges(7,5) * t392;
t194 = Icges(7,5) * t277 + Icges(7,6) * t276 - Icges(7,3) * t556;
t196 = Icges(7,4) * t277 + Icges(7,2) * t276 - Icges(7,6) * t556;
t198 = Icges(7,1) * t277 + Icges(7,4) * t276 - Icges(7,5) * t556;
t418 = t196 * t358 - t198 * t359;
t24 = (qJD(3) * t418 + t104) * t381 + (qJD(3) * t194 + (-t196 * t375 + t108) * t359 + (-t198 * t375 - t106) * t358) * t383;
t103 = Icges(7,5) * t170 + Icges(7,6) * t169 - Icges(7,3) * t634;
t105 = Icges(7,4) * t170 + Icges(7,2) * t169 - Icges(7,6) * t634;
t107 = Icges(7,1) * t170 + Icges(7,4) * t169 - Icges(7,5) * t634;
t195 = Icges(7,5) * t279 + Icges(7,6) * t278 + Icges(7,3) * t555;
t197 = Icges(7,4) * t279 + Icges(7,2) * t278 + Icges(7,6) * t555;
t199 = Icges(7,1) * t279 + Icges(7,4) * t278 + Icges(7,5) * t555;
t417 = t197 * t358 - t199 * t359;
t25 = (qJD(3) * t417 + t103) * t381 + (qJD(3) * t195 + (-t197 * t375 + t107) * t359 + (-t199 * t375 - t105) * t358) * t383;
t125 = -t261 * t556 + t262 * t276 + t263 * t277;
t74 = -t194 * t556 + t196 * t276 + t198 * t277;
t75 = -t195 * t556 + t197 * t276 + t199 * t277;
t444 = t382 * t74 - t384 * t75;
t38 = t125 * t381 - t383 * t444;
t126 = t261 * t555 + t278 * t262 + t279 * t263;
t18 = t104 * t555 + t278 * t106 + t279 * t108 + t169 * t196 + t170 * t198 - t194 * t634;
t19 = t103 * t555 + t278 * t105 + t279 * t107 + t169 * t197 + t170 * t199 - t195 * t634;
t76 = t194 * t555 + t278 * t196 + t279 * t198;
t77 = t195 * t555 + t278 * t197 + t279 * t199;
t443 = t382 * t76 - t384 * t77;
t47 = t169 * t262 + t170 * t263 + t182 * t555 + t278 * t183 + t279 * t184 - t261 * t634;
t54 = t77 * t382 + t384 * t76;
t4 = (qJD(3) * t443 + t47) * t381 + (-qJD(1) * t54 + qJD(3) * t126 - t18 * t382 + t19 * t384) * t383;
t90 = t194 * t381 - t383 * t418;
t91 = t195 * t381 - t383 * t417;
t438 = t90 * t382 - t91 * t384;
t439 = t91 * t382 + t384 * t90;
t483 = t4 * t555 + t38 * t490 + (t137 * t381 - t383 * t438) * t522 + t381 * (t438 * t524 + (-qJD(1) * t439 - t24 * t382 + t25 * t384) * t383 + t590);
t82 = -t204 * t556 + t206 * t289 + t208 * t290;
t83 = -t205 * t556 + t207 * t289 + t209 * t290;
t441 = t382 * t82 - t384 * t83;
t41 = t129 * t381 - t383 * t441;
t478 = t381 * t94 + t41;
t477 = t383 * t638;
t326 = t458 * qJD(3);
t476 = t326 * t530;
t474 = qJD(1) * t542;
t471 = -t598 + t610;
t466 = -t499 + t549;
t465 = t382 * t249 + t265 * t525 + t537;
t39 = t126 * t381 - t383 * t443;
t84 = t204 * t555 + t291 * t206 + t292 * t208;
t85 = t205 * t555 + t291 * t207 + t292 * t209;
t440 = t382 * t84 - t384 * t85;
t42 = t130 * t381 - t383 * t440;
t463 = -t381 * t95 - t39 - t42;
t252 = qJD(1) * t308 + t378 * t487;
t253 = qJD(1) * t309 - t379 * t487;
t457 = -t253 * rSges(5,1) - t252 * rSges(5,2);
t456 = rSges(5,1) * t404 - t310 * rSges(5,2);
t455 = rSges(5,1) * t379 - rSges(5,2) * t378;
t454 = t191 * rSges(6,1) + t190 * rSges(6,2);
t31 = t174 + (-t521 * t542 - t553) * t381 + (-qJD(1) * t475 - qJD(3) * t548 + t202 * t384) * t383;
t32 = t381 * t101 + (t382 * t202 + t235 * t525) * t383 + (t160 * t383 - t381 * t475) * qJD(3) + t467;
t449 = t31 * t382 - t32 * t384;
t123 = -rSges(6,3) * t634 + t454;
t215 = (-rSges(6,1) * t362 - rSges(6,2) * t363) * t519 + (rSges(6,3) * t383 - t381 * t453) * qJD(3);
t66 = (-t275 * t521 - t123) * t381 + (-qJD(3) * t211 + t215 * t384 - t493) * t383;
t124 = t193 * rSges(6,1) + t192 * rSges(6,2) + rSges(6,3) * t392;
t403 = t382 * t215 + t275 * t525;
t67 = (-t275 * t523 + t124) * t381 + (qJD(3) * t210 + t403) * t383;
t447 = t66 * t382 - t384 * t67;
t68 = t382 * t474 + (t541 - t547) * t384 + t540;
t69 = t382 * t547 + t384 * t474 + t465;
t446 = t69 * t382 - t384 * t68;
t393 = qJD(1) * t517 - t520;
t70 = ((-qJD(3) * t373 - t516) * t381 + t393) * t382 + t110 + t464 + t533;
t71 = t366 + (-t515 + (-qJ(2) - t625) * qJD(1)) * t382 + (-t470 + (t381 * t591 + t570) * qJD(3) + t393) * t384 - t452;
t445 = t382 * t71 - t384 * t70;
t53 = t75 * t382 + t384 * t74;
t78 = (qJD(1) * t471 - t380 * t524 - t520) * t382 + t124 + t500 + t533;
t79 = t319 + t366 + (t524 * t592 - t520) * t384 + (t471 * t384 + (-t360 * t381 + t383 * t592 - qJ(2)) * t382) * qJD(1) - t454;
t442 = t382 * t79 - t384 * t78;
t57 = t83 * t382 + t384 * t82;
t58 = t85 * t382 + t384 * t84;
t92 = t493 + (-t215 + t541) * t384 + t540;
t93 = t403 + t465;
t437 = t93 * t382 - t384 * t92;
t436 = Icges(4,1) * t383 - t588;
t431 = -Icges(4,2) * t381 + t587;
t426 = Icges(4,5) * t383 - Icges(4,6) * t381;
t424 = Icges(5,5) * t379 - Icges(5,6) * t378;
t390 = t382 * t610 - t384 * t477;
t111 = qJD(1) * t390 - t382 * t520 + t498 - t502 + t533;
t112 = t366 + (rSges(5,3) * t524 - t520) * t384 + (t610 * t384 + (-qJ(2) + t595 - t599) * t382) * qJD(1) + t457 + t497;
t421 = -t111 * t384 + t382 * t112;
t272 = (-t381 * t455 + t595) * qJD(3);
t288 = rSges(5,3) * t381 + t383 * t455;
t146 = t288 * t527 + t316 + (-t272 - t307) * t384;
t147 = t382 * t272 + t288 * t525 + t537;
t419 = -t146 * t384 + t147 * t382;
t414 = t210 * t384 + t211 * t382;
t410 = -t381 * t622 - t383 * t623;
t405 = rSges(3,3) * t384 + t382 * t601;
t399 = (t613 + t612 + t611) * t524;
t397 = t410 * t382;
t396 = qJD(3) * t436;
t395 = qJD(3) * t431;
t11 = -qJD(1) * t443 + t18 * t384 + t19 * t382;
t20 = -t104 * t556 + t276 * t106 + t277 * t108 + t171 * t196 + t172 * t198 + t194 * t392;
t21 = -t103 * t556 + t276 * t105 + t277 * t107 + t171 * t197 + t172 * t199 + t195 * t392;
t12 = -qJD(1) * t444 + t20 * t384 + t21 * t382;
t48 = t171 * t262 + t172 * t263 - t182 * t556 + t276 * t183 + t277 * t184 + t261 * t392;
t5 = (qJD(3) * t444 + t48) * t381 + (-qJD(1) * t53 + qJD(3) * t125 - t20 * t382 + t21 * t384) * t383;
t387 = t11 * t484 + t12 * t485 + t4 * t603 + t5 * t602 + (-qJD(1) * t438 + t24 * t384 + t25 * t382) * t605 - t38 * t527 / 0.2e1 + t39 * t525 / 0.2e1 + t439 * t522 / 0.2e1 + t637 * t54 + t636 * t53;
t386 = t590 + (t24 + t48) * t485 + (t25 + t47) * t484 + (t126 + t91) * t637 + (t125 + t90) * t636;
t385 = (-t382 * t5 + (-t38 * t384 - t382 * t39) * qJD(1)) * t383 - t39 * t489 + t483;
t369 = t384 * qJ(2);
t323 = t425 * qJD(3);
t303 = -rSges(3,2) * t384 + t382 * rSges(3,3) + t532;
t302 = t369 + t405;
t300 = t593 - t400;
t268 = (Icges(5,5) * t383 - t381 * t434) * qJD(3);
t267 = (Icges(5,6) * t383 - t381 * t429) * qJD(3);
t258 = t366 + (t601 * t384 + (-rSges(3,3) - qJ(2)) * t382) * qJD(1);
t257 = qJD(1) * t405 + t533;
t243 = t495 + t299;
t242 = t382 * t518 + t369 + t400;
t238 = (-t288 - t337) * t384;
t237 = t288 * t382 + t321;
t230 = qJD(1) * t624 + t426 * t523;
t229 = -t426 * t521 + t529;
t225 = rSges(5,3) * t555 - t456;
t223 = -Icges(5,1) * t404 + Icges(5,4) * t310 + Icges(5,5) * t555;
t222 = Icges(5,1) * t309 + Icges(5,4) * t308 - Icges(5,5) * t556;
t221 = -Icges(5,4) * t404 + Icges(5,2) * t310 + Icges(5,6) * t555;
t220 = Icges(5,4) * t309 + Icges(5,2) * t308 - Icges(5,6) * t556;
t165 = (-t275 + t538) * t384;
t164 = t275 * t382 + t539;
t163 = -t382 * t477 + t495 - t629;
t162 = t357 + t369 + t390 + t456;
t159 = -t382 * t624 - t384 * t410;
t158 = t382 * t293 - t398;
t157 = -t384 * t624 + t397;
t156 = t293 * t384 + t382 * t411;
t155 = Icges(5,1) * t255 + Icges(5,4) * t254 + Icges(5,5) * t392;
t154 = Icges(5,1) * t253 + Icges(5,4) * t252 - Icges(5,5) * t634;
t153 = Icges(5,4) * t255 + Icges(5,2) * t254 + Icges(5,6) * t392;
t152 = Icges(5,4) * t253 + Icges(5,2) * t252 - Icges(5,6) * t634;
t149 = -t381 * t211 + t275 * t555;
t148 = t210 * t381 + t275 * t556;
t142 = t495 - t503;
t141 = t382 * t471 - t211 + t369 + t534;
t135 = t200 + t495 + t501;
t134 = t517 * t382 + t384 * t625 + t369 + t451;
t133 = (t538 - t542) * t384;
t132 = t539 + t475;
t131 = t414 * t383;
t128 = (-t200 * t384 - t201 * t382) * t383;
t127 = t225 * t384 + t382 * t545 + t304;
t99 = t310 * t221 - t223 * t404 + t508;
t98 = t310 * t220 - t222 * t404 + t510;
t97 = t221 * t308 + t223 * t309 - t509;
t96 = t220 * t308 + t222 * t309 - t511;
t81 = t235 * t555 - t381 * t548 + t245;
t80 = t160 * t381 + t235 * t556 + t143;
t73 = t382 * t503 + t546 + t578;
t72 = (-t382 * t548 + t384 * t549) * t383;
t61 = t217 + t384 * (-rSges(5,3) * t489 - t457) + (-t228 + t502) * t382 + (t545 * t384 + (-t225 - t315) * t382) * qJD(1);
t59 = t382 * t466 + t384 * t548 + t546;
t43 = t414 * t524 + (-t123 * t382 - t124 * t384 + (t382 * t210 - t578) * qJD(1)) * t383;
t40 = (-t109 * t382 + (-qJD(1) * t201 - t110) * t384) * t383 + t504;
t30 = -t117 * t556 + t289 * t119 + t290 * t121 + t192 * t207 + t193 * t209 + t205 * t392;
t29 = -t118 * t556 + t289 * t120 + t290 * t122 + t192 * t206 + t193 * t208 + t204 * t392;
t28 = t117 * t555 + t291 * t119 + t292 * t121 + t190 * t207 + t191 * t209 - t205 * t634;
t27 = t118 * t555 + t291 * t120 + t292 * t122 + t190 * t206 + t191 * t208 - t204 * t634;
t26 = t123 * t384 + (-t124 + t550) * t382 + (t503 * t384 + (-t211 + t543) * t382) * qJD(1) + t551;
t17 = (t160 * t384 + t161 * t382) * t524 + ((qJD(1) * t160 - t553) * t382 + (-qJD(1) * t548 + t552) * t384) * t383 + t504;
t16 = t553 * t384 + (t550 + t552) * t382 + (t466 * t384 + (t543 - t548) * t382) * qJD(1) + t551;
t15 = -qJD(1) * t441 + t29 * t384 + t30 * t382;
t14 = -qJD(1) * t440 + t27 * t384 + t28 * t382;
t7 = (qJD(3) * t441 + t56) * t381 + (-t57 * qJD(1) + qJD(3) * t129 - t29 * t382 + t30 * t384) * t383;
t6 = (qJD(3) * t440 + t55) * t381 + (-qJD(1) * t58 + qJD(3) * t130 - t27 * t382 + t28 * t384) * t383;
t1 = [-t435 * t522 + 0.2e1 * m(3) * (t257 * t303 + t258 * t302) + (t179 * t243 + t180 * t242) * t617 + (t111 * t163 + t112 * t162) * t616 + (t134 * t71 + t135 * t70) * t614 + (t141 * t79 + t142 * t78) * t615 - t271 * t486 - t519 * t573 - t358 * t263 * t569 + t388 + t389 + (Icges(5,3) * t522 - t396) * t381 + (t283 * t378 - t284 * t379 - t381 * t424 + t430) * t524 + (Icges(5,3) * t524 - t183 * t358 - t267 * t378 + t268 * t379 + t424 * t522 - t395 - t507 - t577) * t383; m(7) * ((t134 * t384 + t135 * t382) * qJD(1) + t445) + m(6) * ((t141 * t384 + t142 * t382) * qJD(1) + t442) + m(5) * ((t162 * t384 + t163 * t382) * qJD(1) + t421) + m(4) * ((t242 * t384 + t243 * t382) * qJD(1) + t628) + m(3) * (-t257 * t384 + t382 * t258 + (t302 * t384 + t303 * t382) * qJD(1)); 0; (t254 * t609 + t255 * t608 + t308 * t267 / 0.2e1 + t309 * t268 / 0.2e1 + t48 / 0.2e1 + t24 / 0.2e1 - t323 * t602 - t514) * t384 + (t252 * t609 + t253 * t608 + t267 * t607 + t268 * t606 + t47 / 0.2e1 + t25 / 0.2e1 - t323 * t603 + t513) * t382 + m(4) * (t628 * t338 - (t242 * t382 - t243 * t384) * t326) + m(5) * (t111 * t238 + t112 * t237 + t146 * t163 + t147 * t162) + m(6) * (t141 * t93 + t142 * t92 + t164 * t79 + t165 * t78) + m(7) * (t132 * t71 + t133 * t70 + t134 * t69 + t135 * t68) + ((t254 * t641 + t255 * t642 + t392 * t640 + t395 * t604 + t623 * t639) * t384 + (t253 * t642 + t252 * t641 - t634 * t640 + t295 * t639 + t431 * t521 / 0.2e1) * t382) * t381 + ((qJD(1) * t297 - t152 * t378 + t154 * t379 - t436 * t521) * t603 + (qJD(1) * t622 - t153 * t378 + t155 * t379 + t382 * t396) * t602) * t383 + ((t218 * t383 + t220 * t568 - t222 * t566) * t602 - t398 / 0.2e1 + (t219 * t383 + t221 * t568 - t223 * t566) * t603 + t410 * t604) * qJD(3) + ((t243 * t600 - t283 * t308 / 0.2e1 - t284 * t309 / 0.2e1 - t125 / 0.2e1 - t90 / 0.2e1 + (t295 / 0.2e1 - t218 / 0.2e1) * t381 + (-t297 / 0.2e1 + t220 * t378 / 0.2e1 - t222 * t379 / 0.2e1) * t383 - t506) * t382 + (t283 * t607 + t284 * t606 + t242 * t600 + t126 / 0.2e1 + t91 / 0.2e1 + (t219 / 0.2e1 + t623 / 0.2e1) * t381 + (-t221 * t378 / 0.2e1 + t223 * t379 / 0.2e1 - t622 / 0.2e1) * t383 - t505) * t384) * qJD(1); m(5) * ((t237 * t384 + t238 * t382) * qJD(1) + t419) + m(6) * ((t164 * t384 + t165 * t382) * qJD(1) + t437) + m(7) * ((t132 * t384 + t133 * t382) * qJD(1) + t446) - m(4) * t476; (t132 * t69 + t133 * t68 + t59 * t16) * t614 + t384 * t12 + t382 * t11 + (t164 * t93 + t165 * t92 + t26 * t73) * t615 + t384 * t15 + t382 * t14 + (t127 * t61 + t146 * t238 + t147 * t237) * t616 + ((-t382 * t299 + t300 * t384) * (-t382 * t496 + (-t338 * t377 + t376 * t596) * qJD(3) + ((-t299 + t370) * t384 + (-t300 + t400 + t593) * t382) * qJD(1)) - t338 * t476) * t617 + t382 * ((t382 * t229 + (-t158 + t397) * qJD(1)) * t382 + (t159 * qJD(1) + (t295 * t524 - t297 * t522 + t529) * t384 + (t230 + (-t381 * t623 + t383 * t622) * qJD(3) + t411 * qJD(1)) * t382) * t384) + t382 * ((t310 * t152 - t404 * t154 + t252 * t221 + t253 * t223 + (-t98 - t509) * qJD(1)) * t382 + (t310 * t153 - t404 * t155 + t252 * t220 + t253 * t222 + (t99 - t511) * qJD(1)) * t384) + t384 * ((t308 * t153 + t309 * t155 + t254 * t220 + t255 * t222 + (t97 - t510) * qJD(1)) * t384 + (t308 * t152 + t309 * t154 + t254 * t221 + t255 * t223 + (-t96 - t508) * qJD(1)) * t382) + t384 * ((t230 * t384 + (t157 + t398) * qJD(1)) * t384 + (-t156 * qJD(1) + (-t522 * t622 + t524 * t623) * t382 + (t229 + (-t295 * t381 + t297 * t383) * qJD(3) + (-t293 + t410) * qJD(1)) * t384) * t382) + (-t53 - t57 + (-t156 - t96) * t384 + (-t157 - t97) * t382) * t527 + (t54 + t58 + (t158 + t98) * t384 + (t159 + t99) * t382) * t525; 0.2e1 * ((t134 * t382 - t135 * t384) * t611 + (t141 * t382 - t142 * t384) * t612 + (t162 * t382 - t163 * t384) * t613) * t524 + 0.2e1 * ((-t134 * t525 - t135 * t527 - t445) * t611 + (-t141 * t525 - t142 * t527 - t442) * t612 + (-t162 * t525 - t163 * t527 - t421) * t613) * t383; 0.2e1 * t530 * t399; 0.2e1 * ((t132 * t523 - t133 * t521 + t16) * t611 + (t164 * t523 - t165 * t521 + t26) * t612 + (t237 * t523 - t238 * t521 + t61) * t613) * t381 + 0.2e1 * ((qJD(3) * t59 - t132 * t525 - t133 * t527 - t446) * t611 + (qJD(3) * t73 - t164 * t525 - t165 * t527 - t437) * t612 + (qJD(3) * t127 - t237 * t525 - t238 * t527 - t419) * t613) * t383; 0.4e1 * (0.1e1 - t530) * t383 * t399; t386 + (t513 * t384 + t514 * t382 + (t382 * t505 - t384 * t506) * qJD(1)) * t383 + m(6) * (t141 * t66 + t142 * t67 + t148 * t78 + t149 * t79) + m(7) * (t134 * t31 + t135 * t32 + t70 * t80 + t71 * t81) + (t382 * t506 + t384 * t505) * t524 + t589; m(6) * ((t148 * t382 + t149 * t384) * qJD(1) + t447) + m(7) * ((t382 * t80 + t384 * t81) * qJD(1) + t449); (qJD(3) * (t95 * t382 + t384 * t94) / 0.2e1 + t14 * t602 + t15 * t604 + (-t384 * t57 / 0.2e1 + t58 * t604) * qJD(1)) * t383 + t387 + ((qJD(1) * t95 + t33) * t605 + t7 / 0.2e1 + qJD(1) * t42 / 0.2e1 + t58 * t481) * t384 + ((-qJD(1) * t94 + t34) * t605 + t6 / 0.2e1 + t41 * t639 + t57 * t480) * t382 + m(6) * (-t131 * t26 + t148 * t92 + t149 * t93 + t164 * t66 + t165 * t67 + t43 * t73) + m(7) * (t132 * t31 + t133 * t32 + t16 * t72 + t17 * t59 + t68 * t80 + t69 * t81); 0.2e1 * ((-t148 * t521 + t149 * t523 + t43) * t612 + (-t521 * t80 + t523 * t81 + t17) * t611) * t381 + 0.2e1 * ((-qJD(3) * t131 - t148 * t527 - t149 * t525 - t447) * t612 + (qJD(3) * t72 - t525 * t81 - t527 * t80 - t449) * t611) * t383; (t17 * t72 + t31 * t81 + t32 * t80) * t614 + (-t131 * t43 + t148 * t67 + t149 * t66) * t615 + ((t382 * t478 + t384 * t463) * qJD(3) + t589) * t381 + ((t34 * t381 + t6) * t384 + (-t33 * t381 - t5 - t7) * t382 + (t140 * t381 + (-t382 * t94 + t384 * t95) * t383) * qJD(3) + ((-t38 - t478) * t384 + t463 * t382) * qJD(1)) * t383 + t483; t386 + m(7) * (t134 * t63 + t135 * t64 + t143 * t70 + t144 * t71); m(7) * t619; t387 + m(7) * (t128 * t16 + t132 * t63 + t133 * t64 + t143 * t68 + t144 * t69 + t40 * t59); m(7) * ((t40 + (-t143 * t384 + t144 * t382) * qJD(3)) * t381 + (qJD(3) * t128 - t619) * t383); m(7) * (t128 * t17 + t143 * t32 + t144 * t31 + t40 * t72 + t63 * t81 + t64 * t80) + t385; (t128 * t40 + t143 * t64 + t144 * t63) * t614 + t385;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;