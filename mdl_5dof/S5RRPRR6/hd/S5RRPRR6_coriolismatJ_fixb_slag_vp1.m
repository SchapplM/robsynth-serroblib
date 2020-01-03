% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:31
% EndTime: 2020-01-03 12:05:49
% DurationCPUTime: 13.77s
% Computational Cost: add. (80557->534), mult. (74468->711), div. (0->0), fcn. (80559->10), ass. (0->326)
t481 = qJ(1) + qJ(2);
t475 = sin(t481);
t478 = cos(qJ(1)) * pkin(1);
t480 = qJ(4) + qJ(5);
t474 = sin(t480);
t476 = cos(t480);
t477 = cos(t481);
t483 = cos(pkin(9));
t580 = t477 * t483;
t420 = t474 * t580 - t475 * t476;
t421 = t474 * t475 + t476 * t580;
t482 = sin(pkin(9));
t581 = t477 * t482;
t342 = t421 * rSges(6,1) - t420 * rSges(6,2) + rSges(6,3) * t581;
t486 = cos(qJ(4));
t473 = pkin(4) * t486 + pkin(3);
t531 = t482 * (-pkin(8) - pkin(7));
t484 = sin(qJ(4));
t583 = t475 * t484;
t497 = -pkin(4) * t583 - t473 * t580 + t477 * t531;
t534 = t477 * pkin(2) + t475 * qJ(3);
t676 = -t497 + t534 + t342;
t689 = t478 + t676;
t266 = t689 * t475;
t468 = t477 * qJ(3);
t584 = t475 * t483;
t418 = -t474 * t584 - t476 * t477;
t419 = -t477 * t474 + t476 * t584;
t510 = t419 * rSges(6,1) + t418 * rSges(6,2);
t579 = t477 * t484;
t536 = -pkin(4) * t579 - t475 * t531;
t276 = -t468 + (rSges(6,3) * t482 + t473 * t483 + pkin(2)) * t475 + t510 + t536;
t611 = sin(qJ(1)) * pkin(1);
t273 = t276 + t611;
t572 = t483 * t484;
t440 = -t475 * t486 + t477 * t572;
t571 = t483 * t486;
t441 = t477 * t571 + t583;
t360 = t441 * rSges(5,1) - t440 * rSges(5,2) + rSges(5,3) * t581;
t535 = pkin(3) * t580 + pkin(7) * t581;
t708 = t360 + t534 + t535;
t713 = t478 + t708;
t285 = t713 * t475;
t438 = -t475 * t572 - t477 * t486;
t439 = t475 * t571 - t579;
t511 = t439 * rSges(5,1) + t438 * rSges(5,2);
t304 = -t468 + (pkin(3) * t483 + pkin(2) + (rSges(5,3) + pkin(7)) * t482) * t475 + t511;
t292 = t304 + t611;
t383 = rSges(4,1) * t580 - rSges(4,2) * t581 + t475 * rSges(4,3) + t534;
t381 = t383 + t478;
t370 = t381 * t475;
t377 = t383 * t475;
t585 = t475 * t482;
t382 = -rSges(4,2) * t585 - t468 - t477 * rSges(4,3) + (rSges(4,1) * t483 + pkin(2)) * t475;
t380 = t382 + t611;
t664 = m(6) / 0.2e1;
t665 = m(5) / 0.2e1;
t695 = t475 * t676;
t702 = m(4) / 0.2e1;
t717 = t475 * t708;
t527 = (t266 + t695 + (-t273 - t276) * t477) * t664 + (t285 + t717 + (-t292 - t304) * t477) * t665 + (t370 + t377 + (-t380 - t382) * t477) * t702;
t723 = t292 - t304;
t724 = t273 - t276;
t528 = (-t724 * t477 + t266 - t695) * t664 + (-t723 * t477 + t285 - t717) * t665 + (t370 - t377 + (-t380 + t382) * t477) * t702;
t18 = t528 - t527;
t725 = t18 * qJD(1);
t368 = t418 * rSges(6,1) - t419 * rSges(6,2);
t369 = t420 * rSges(6,1) + t421 * rSges(6,2);
t141 = t276 * t368 - t369 * t676;
t722 = m(6) * t141;
t138 = t273 * t368 - t369 * t689;
t721 = m(6) * t138;
t329 = t483 * t342;
t358 = t497 + t535;
t413 = -t483 * rSges(6,3) + (rSges(6,1) * t476 - rSges(6,2) * t474) * t482;
t610 = -pkin(3) + t473;
t541 = -t483 * pkin(8) + t610 * t482 + t413;
t512 = t541 * t581;
t201 = t483 * t358 - t329 - t512;
t202 = (-t342 + t358) * t483 - t512;
t720 = m(6) * (t201 - t202);
t136 = t292 * t708 - t304 * t713;
t332 = Icges(6,5) * t419 + Icges(6,6) * t418 + Icges(6,3) * t585;
t715 = t332 * t581;
t714 = -t413 * t581 - t329;
t333 = Icges(6,5) * t421 - Icges(6,6) * t420 + Icges(6,3) * t581;
t712 = t333 * t581;
t341 = rSges(6,3) * t585 + t510;
t357 = (-pkin(7) * t482 + t610 * t483) * t475 + t536;
t200 = (t341 + t357) * t483 + t541 * t585;
t711 = t200 * t689;
t600 = Icges(6,4) * t474;
t412 = -Icges(6,5) * t483 + (Icges(6,1) * t476 - t600) * t482;
t539 = t412 + (-Icges(6,2) * t476 - t600) * t482;
t709 = t539 * t474;
t114 = t273 * t676 - t276 * t689;
t434 = -t483 * rSges(5,3) + (rSges(5,1) * t486 - rSges(5,2) * t484) * t482;
t707 = -t360 * t483 - t434 * t581;
t403 = Icges(6,4) * t421;
t337 = Icges(6,2) * t420 - Icges(6,6) * t581 - t403;
t402 = Icges(6,4) * t420;
t339 = Icges(6,1) * t421 + Icges(6,5) * t581 - t402;
t706 = -t418 * t337 + t419 * t339;
t424 = Icges(5,4) * t441;
t353 = Icges(5,2) * t440 - Icges(5,6) * t581 - t424;
t423 = Icges(5,4) * t440;
t355 = Icges(5,1) * t441 + Icges(5,5) * t581 - t423;
t704 = -t438 * t353 + t439 * t355;
t507 = t353 * t440 + t355 * t441;
t645 = m(3) * (-t478 * (rSges(3,1) * t475 + rSges(3,2) * t477) + t611 * (t477 * rSges(3,1) - rSges(3,2) * t475));
t641 = m(4) * (t380 * t383 - t382 * t381);
t162 = -t333 * t585 - t706;
t593 = t333 * t475;
t594 = t332 * t477;
t694 = t482 * t477 ^ 2;
t701 = (t333 * t694 - t712 * t477 + (t162 + t706 + (t593 + t594) * t482 - t715) * t475) * t482;
t700 = t200 * t676;
t431 = -Icges(5,3) * t483 + (Icges(5,5) * t486 - Icges(5,6) * t484) * t482;
t602 = Icges(5,4) * t486;
t432 = -Icges(5,6) * t483 + (-Icges(5,2) * t484 + t602) * t482;
t603 = Icges(5,4) * t484;
t433 = -Icges(5,5) * t483 + (Icges(5,1) * t486 - t603) * t482;
t698 = (t431 * t581 - t432 * t440 + t433 * t441) * t483;
t578 = t482 * t484;
t537 = t433 + (-Icges(5,2) * t486 - t603) * t482;
t193 = t201 * t475;
t359 = rSges(5,3) * t585 + t511;
t294 = t359 * t483 + t434 * t585;
t570 = (t200 * t477 + t193) * t664 + (t294 * t477 + t707 * t475) * t665;
t605 = (-t202 * t475 + t193) * t664;
t687 = t570 - t605;
t568 = -t201 * t276 - t700;
t505 = (t202 * t276 + t568 + t700) * t664;
t569 = t202 * t273 + t711;
t606 = (t568 + t569) * t664 + (t723 * t707 + (-t708 + t713) * t294) * t665;
t686 = t606 - t505;
t445 = (-Icges(5,5) * t484 - Icges(5,6) * t486) * t482;
t573 = t483 * t445;
t685 = -t573 / 0.2e1 - t537 * t578 / 0.2e1;
t435 = (-Icges(6,5) * t474 - Icges(6,6) * t476) * t482;
t574 = t483 * t435;
t684 = -t574 / 0.2e1 - t482 * t709 / 0.2e1;
t426 = t440 * pkin(4);
t247 = t341 * t581 - t342 * t585;
t155 = (t357 * t477 + t358 * t475) * t482 + t247;
t263 = t368 * t581 + t369 * t585;
t442 = (-rSges(6,1) * t474 - rSges(6,2) * t476) * t482;
t306 = -t483 * t368 - t442 * t585;
t347 = t483 * t369;
t307 = -t442 * t581 + t347;
t362 = Icges(6,5) * t418 - Icges(6,6) * t419;
t401 = Icges(6,4) * t418;
t338 = Icges(6,1) * t419 + Icges(6,5) * t585 + t401;
t548 = -Icges(6,2) * t419 + t338 + t401;
t601 = Icges(6,4) * t419;
t335 = Icges(6,2) * t418 + Icges(6,6) * t585 + t601;
t550 = -Icges(6,1) * t418 + t335 + t601;
t142 = -t483 * t362 + (-t548 * t474 - t550 * t476) * t482;
t363 = Icges(6,5) * t420 + Icges(6,6) * t421;
t547 = Icges(6,2) * t421 - t339 + t402;
t549 = -Icges(6,1) * t420 + t337 - t403;
t143 = -t483 * t363 + (-t547 * t474 - t549 * t476) * t482;
t599 = Icges(6,4) * t476;
t411 = -Icges(6,6) * t483 + (-Icges(6,2) * t474 + t599) * t482;
t437 = (-Icges(6,1) * t474 - t599) * t482;
t540 = t411 - t437;
t175 = t539 * t418 - t419 * t540 + t435 * t585;
t176 = t539 * t420 + t421 * t540 - t435 * t581;
t519 = -t581 / 0.2e1;
t522 = t585 / 0.2e1;
t598 = (-t574 + (-t476 * t540 - t709) * t482) * t483;
t647 = -t483 / 0.2e1;
t532 = (-t598 + (t142 * t475 - t143 * t477) * t482) * t647 + (-t175 * t483 + (t362 * t585 + t548 * t418 - t550 * t419) * t585 - (t363 * t585 + t547 * t418 - t549 * t419) * t581) * t522 + (-t176 * t483 + (-t362 * t581 + t548 * t420 + t550 * t421) * t585 - (-t363 * t581 + t547 * t420 + t549 * t421) * t581) * t519;
t14 = t532 + m(6) * (t155 * t263 - t200 * t306 + t201 * t307);
t683 = t14 * qJD(5);
t604 = Icges(5,4) * t439;
t351 = Icges(5,2) * t438 + Icges(5,6) * t585 + t604;
t422 = Icges(5,4) * t438;
t354 = Icges(5,1) * t439 + Icges(5,5) * t585 + t422;
t552 = -t440 * t351 + t441 * t354;
t533 = qJD(1) + qJD(2);
t319 = -t369 - t426;
t378 = rSges(5,1) * t438 - rSges(5,2) * t439;
t379 = rSges(5,1) * t440 + rSges(5,2) * t441;
t425 = t438 * pkin(4);
t542 = -t368 - t425;
t557 = (-t319 * t477 + t475 * t542) * t664 + (-t378 * t475 + t379 * t477) * t665;
t129 = -t273 * t542 + t319 * t689;
t134 = -t276 * t542 + t319 * t676;
t154 = t292 * t378 - t379 * t713;
t160 = t304 * t378 - t379 * t708;
t675 = (t160 + t154) * t665 + (t134 + t129) * t664;
t674 = t200 * t720;
t671 = 0.4e1 * qJD(1);
t670 = 0.2e1 * qJD(2);
t669 = 0.4e1 * qJD(2);
t668 = 2 * qJD(4);
t281 = t483 * t341 + t413 * t585;
t661 = t281 * t720;
t563 = t306 * t273 + t307 * t689;
t566 = -t200 * t368 - t201 * t369;
t657 = m(6) * (t563 + t566);
t562 = t306 * t276 + t307 * t676;
t656 = m(6) * (t562 + t566);
t655 = m(6) * (t724 * t714 + (-t676 + t689) * t281);
t558 = t281 * t542 + t319 * t714;
t652 = m(6) * (t558 + t563);
t650 = m(6) * (t558 + t562);
t648 = m(6) * (t141 + t138);
t348 = Icges(5,5) * t439 + Icges(5,6) * t438 + Icges(5,3) * t585;
t180 = t348 * t585 + t438 * t351 + t439 * t354;
t349 = Icges(5,5) * t441 - Icges(5,6) * t440 + Icges(5,3) * t581;
t181 = -t349 * t585 - t704;
t595 = (t431 * t585 + t432 * t438 + t433 * t439) * t483;
t100 = -t595 + (t180 * t475 - t181 * t477) * t482;
t182 = -t348 * t581 - t552;
t183 = t349 * t581 + t507;
t101 = t698 + (t182 * t475 - t183 * t477) * t482;
t591 = t349 * t475;
t592 = t348 * t477;
t41 = -t595 + ((t180 + t183 - t507) * t475 + (t182 + (-t591 + t592) * t482 - t181 + t552) * t477) * t482;
t42 = -t698 + (t349 * t694 + t507 * t477 + (t552 + t181 + (t591 + t592) * t482 + t704) * t475) * t482;
t161 = t332 * t585 + t418 * t335 + t419 * t338;
t597 = ((-Icges(6,3) * t483 + (Icges(6,5) * t476 - Icges(6,6) * t474) * t482) * t585 + t411 * t418 + t412 * t419) * t483;
t37 = -t597 + ((t161 + t712) * t475 + ((-t593 + t594) * t482 - t162 - t715) * t477) * t482;
t518 = t581 / 0.2e1;
t93 = -t597 + (t161 * t475 - t162 * t477) * t482;
t513 = t37 * t519 + t93 * t518 + t701 * t522;
t2 = ((t100 / 0.2e1 - t41 / 0.2e1) * t477 + (t42 / 0.2e1 + t101 / 0.2e1) * t475) * t482 + t513 + t674;
t29 = t570 + t605 - t557;
t646 = t29 * qJD(3) + t2 * qJD(4);
t639 = m(4) * (-t380 * t477 + t370);
t638 = m(4) * (-t382 * t477 + t377);
t633 = m(5) * t136;
t631 = m(5) * t154;
t630 = m(5) * t160;
t629 = m(5) * (-t292 * t477 + t285);
t627 = m(5) * (-t304 * t477 + t717);
t623 = m(6) * t114;
t619 = m(6) * t129;
t618 = m(6) * t134;
t617 = m(6) * (-t273 * t477 + t266);
t616 = m(6) * (-t276 * t477 + t695);
t615 = m(6) * (t281 * t477 + t714 * t475);
t614 = m(6) * (-t306 * t475 - t307 * t477);
t612 = m(6) * (-t368 * t475 + t369 * t477);
t179 = t615 / 0.2e1;
t74 = t179 - t612 / 0.2e1;
t609 = t74 * qJD(3) + qJD(5) * t513;
t27 = t557 - t687;
t265 = t612 / 0.2e1;
t72 = t265 - t615 / 0.2e1;
t608 = t27 * qJD(4) + t72 * qJD(5);
t28 = t557 + t687;
t73 = t265 + t179;
t607 = t28 * qJD(4) + t73 * qJD(5);
t447 = (-Icges(5,1) * t484 - t602) * t482;
t538 = t432 - t447;
t596 = (-t573 + (-t484 * t537 - t486 * t538) * t482) * t483;
t582 = t476 * t482;
t577 = t482 * t486;
t567 = t200 * t542 + t201 * t319;
t556 = -t281 * t368 - t369 * t714;
t555 = -t294 * t378 - t379 * t707;
t546 = -Icges(5,1) * t438 + t351 + t604;
t545 = -Icges(5,1) * t440 + t353 - t424;
t544 = -Icges(5,2) * t439 + t354 + t422;
t543 = Icges(5,2) * t441 - t355 + t423;
t523 = -t585 / 0.2e1;
t521 = -t582 / 0.2e1;
t520 = t582 / 0.2e1;
t515 = -t577 / 0.2e1;
t514 = t577 / 0.2e1;
t509 = t411 * t521 + t437 * t520 + t684;
t506 = t661 / 0.2e1 + t513;
t503 = t648 / 0.2e1 + t509;
t502 = (pkin(4) * t578 - t442) * t482;
t496 = t411 * t520 + t437 * t521 - t684;
t494 = t432 * t515 + t447 * t514 + t509 + t685;
t493 = t37 * t518 - t598 + t701 * t523 + (t142 + t175) * t522 + (t93 + t143 + t176) * t519;
t492 = t494 + t675;
t491 = -t661 / 0.2e1 + t493;
t489 = t432 * t514 + t447 * t515 + t496 - t685;
t371 = Icges(5,5) * t438 - Icges(5,6) * t439;
t152 = -t483 * t371 + (-t544 * t484 - t546 * t486) * t482;
t372 = Icges(5,5) * t440 + Icges(5,6) * t441;
t153 = -t483 * t372 + (-t543 * t484 - t545 * t486) * t482;
t215 = t438 * t537 - t439 * t538 + t445 * t585;
t216 = t440 * t537 + t441 * t538 - t445 * t581;
t488 = t28 * qJD(3) + (t41 * t518 + t493 - t596 - t674 + (t100 + t153 + t216) * t519 + (t101 + t42) * t523 + (t152 + t215) * t522) * qJD(4);
t448 = (-rSges(5,1) * t484 - rSges(5,2) * t486) * t482;
t321 = t379 * t483 - t448 * t581;
t320 = -t378 * t483 - t448 * t585;
t262 = t483 * t426 + t477 * t502 + t347;
t261 = t502 * t475 + t542 * t483;
t219 = (t425 * t477 + t426 * t475) * t482 + t263;
t209 = qJD(5) * t614;
t107 = t509 + t722;
t106 = t509 + t721;
t105 = t616 + t627 + t638;
t102 = t617 + t629 + t639;
t78 = t650 / 0.2e1;
t75 = t652 / 0.2e1;
t70 = t73 * qJD(3);
t61 = t655 / 0.2e1;
t59 = t656 / 0.2e1;
t57 = t657 / 0.2e1;
t56 = t494 + t618 + t630;
t55 = t494 + t619 + t631;
t49 = t623 + t633 + t641 + t645;
t22 = -t655 / 0.2e1 + t503;
t21 = t61 + t503;
t20 = t527 + t528;
t17 = t61 - t648 / 0.2e1 + t496;
t16 = m(6) * (t247 * t263 - t281 * t306 + t307 * t714) + t532;
t15 = t16 * qJD(5);
t13 = t492 + t686;
t12 = t492 - t686;
t9 = t489 + t505 + t606 - t675;
t8 = t78 - t656 / 0.2e1 + t506;
t7 = t59 - t650 / 0.2e1 + t506;
t6 = t75 - t657 / 0.2e1 + t506;
t5 = t57 - t652 / 0.2e1 + t506;
t4 = t491 + t59 + t78;
t3 = t491 + t57 + t75;
t1 = [t49 * qJD(2) + t102 * qJD(3) + t55 * qJD(4) + t106 * qJD(5), t49 * qJD(1) + t20 * qJD(3) + t13 * qJD(4) + t21 * qJD(5) + (t645 / 0.2e1 + t641 / 0.2e1 + t136 * t665 + t114 * t664) * t670, qJD(1) * t102 + qJD(2) * t20 + t607, t55 * qJD(1) + t13 * qJD(2) + t3 * qJD(5) + ((t292 * t320 + t321 * t713 + t555) * t665 + (t261 * t273 + t262 * t689 + t567) * t664) * t668 + t488, t106 * qJD(1) + t21 * qJD(2) + t70 + t3 * qJD(4) + (m(6) * (t556 + t563) + t493) * qJD(5); -t18 * qJD(3) + t12 * qJD(4) + t22 * qJD(5) + (-t645 / 0.4e1 - t641 / 0.4e1 - t633 / 0.4e1 - t623 / 0.4e1) * t671, qJD(3) * t105 + qJD(4) * t56 + qJD(5) * t107, qJD(2) * t105 + t607 - t725, t12 * qJD(1) + t56 * qJD(2) + t4 * qJD(5) + ((t261 * t276 + t262 * t676 + t567) * t664 + (t304 * t320 + t321 * t708 + t555) * t665) * t668 + t488, t22 * qJD(1) + t107 * qJD(2) + t70 + t4 * qJD(4) + (m(6) * (t556 + t562) + t493) * qJD(5); t18 * qJD(2) + (-t639 / 0.4e1 - t629 / 0.4e1 - t617 / 0.4e1) * t671 + t608, t725 + (-t616 / 0.4e1 - t627 / 0.4e1 - t638 / 0.4e1) * t669 + t608, 0, ((-t320 * t475 - t321 * t477) * t665 + (-t261 * t475 - t262 * t477) * t664) * t668 + t209 + t533 * t27, qJD(4) * t614 + t533 * t72 + t209; t9 * qJD(2) + t5 * qJD(5) + (-t631 / 0.4e1 - t619 / 0.4e1) * t671 + t646 + (t489 + 0.2e1 * (-t201 * t273 + t569 - t711) * t664) * qJD(1), t9 * qJD(1) + t489 * qJD(2) + t7 * qJD(5) + (-t618 / 0.4e1 - t630 / 0.4e1) * t669 + t505 * t670 + t646, t533 * t29, (m(5) * (-t294 * t320 + t707 * t321 + (t359 * t477 - t360 * t475) * t482 ^ 2 * (t378 * t477 + t379 * t475)) + (-t596 + (t152 * t475 - t153 * t477) * t482) * t647 + (-t215 * t483 + (t371 * t585 + t438 * t544 - t439 * t546) * t585 - (t372 * t585 + t438 * t543 - t439 * t545) * t581) * t522 + (-t216 * t483 + (-t371 * t581 + t440 * t544 + t441 * t546) * t585 - (-t372 * t581 + t440 * t543 + t441 * t545) * t581) * t519 + m(6) * (t155 * t219 - t200 * t261 + t201 * t262) + t532) * qJD(4) + t683 + t533 * t2, t5 * qJD(1) + t7 * qJD(2) + t14 * qJD(4) + t683; (t496 - t721) * qJD(1) + t17 * qJD(2) + t6 * qJD(4) + t609, t17 * qJD(1) + (t496 - t722) * qJD(2) + t8 * qJD(4) + t609, t533 * t74, t6 * qJD(1) + t8 * qJD(2) + ((t219 * t247 - t261 * t281 + t262 * t714) * m(6) + t532) * qJD(4) + t15, qJD(4) * t16 + t513 * t533 + t15;];
Cq = t1;
