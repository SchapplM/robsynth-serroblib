% Calculate time derivative of joint inertia matrix for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:08:13
% EndTime: 2019-03-09 07:08:48
% DurationCPUTime: 18.99s
% Computational Cost: add. (77376->1026), mult. (61997->1396), div. (0->0), fcn. (59476->12), ass. (0->521)
t429 = pkin(11) + qJ(3);
t420 = qJ(4) + t429;
t413 = cos(t420);
t434 = qJ(5) + qJ(6);
t424 = cos(t434);
t441 = cos(qJ(1));
t627 = t424 * t441;
t423 = sin(t434);
t439 = sin(qJ(1));
t630 = t423 * t439;
t363 = -t413 * t630 - t627;
t628 = t424 * t439;
t629 = t423 * t441;
t364 = t413 * t628 - t629;
t412 = sin(t420);
t634 = t412 * t439;
t250 = Icges(7,5) * t364 + Icges(7,6) * t363 + Icges(7,3) * t634;
t252 = Icges(7,4) * t364 + Icges(7,2) * t363 + Icges(7,6) * t634;
t254 = Icges(7,1) * t364 + Icges(7,4) * t363 + Icges(7,5) * t634;
t495 = -t252 * t423 + t254 * t424;
t128 = -t250 * t413 + t412 * t495;
t507 = Icges(7,5) * t424 - Icges(7,6) * t423;
t322 = -Icges(7,3) * t413 + t412 * t507;
t657 = Icges(7,4) * t424;
t511 = -Icges(7,2) * t423 + t657;
t323 = -Icges(7,6) * t413 + t412 * t511;
t658 = Icges(7,4) * t423;
t516 = Icges(7,1) * t424 - t658;
t324 = -Icges(7,5) * t413 + t412 * t516;
t162 = t322 * t634 + t323 * t363 + t324 * t364;
t709 = -t162 - t128;
t365 = -t413 * t629 + t628;
t366 = t413 * t627 + t630;
t633 = t412 * t441;
t251 = Icges(7,5) * t366 + Icges(7,6) * t365 + Icges(7,3) * t633;
t253 = Icges(7,4) * t366 + Icges(7,2) * t365 + Icges(7,6) * t633;
t255 = Icges(7,1) * t366 + Icges(7,4) * t365 + Icges(7,5) * t633;
t494 = -t253 * t423 + t255 * t424;
t129 = -t251 * t413 + t412 * t494;
t163 = t322 * t633 + t323 * t365 + t324 * t366;
t708 = -t163 - t129;
t596 = qJD(1) * t439;
t562 = t412 * t596;
t430 = qJD(5) + qJD(6);
t542 = -t413 * t430 + qJD(1);
t481 = t441 * t542;
t597 = qJD(1) * t413;
t541 = -t430 + t597;
t431 = qJD(3) + qJD(4);
t624 = t431 * t441;
t579 = t412 * t624;
t690 = t439 * t541 + t579;
t229 = t423 * t690 + t424 * t481;
t230 = t423 * t481 - t424 * t690;
t577 = t413 * t624;
t572 = t230 * rSges(7,1) + t229 * rSges(7,2) + rSges(7,3) * t577;
t151 = -rSges(7,3) * t562 + t572;
t394 = pkin(9) * t577;
t442 = -pkin(10) - pkin(9);
t438 = sin(qJ(5));
t591 = qJD(5) * t438;
t585 = pkin(5) * t591;
t471 = -t431 * t442 - t585;
t440 = cos(qJ(5));
t590 = qJD(5) * t440;
t584 = pkin(5) * t590;
t621 = t438 * t441;
t589 = pkin(5) * t621;
t564 = qJD(1) * t589 + t439 * t584 + t442 * t562;
t415 = pkin(5) * t440 + pkin(4);
t675 = pkin(4) - t415;
t707 = t151 - t394 + (pkin(9) * t596 + t624 * t675) * t412 + (t441 * t471 + t596 * t675) * t413 + t564;
t522 = -t364 * rSges(7,1) - t363 * rSges(7,2);
t257 = rSges(7,3) * t634 - t522;
t674 = pkin(9) + t442;
t454 = -t412 * t674 - t413 * t675;
t279 = t439 * t454 - t589;
t706 = t257 + t279;
t632 = t413 * t431;
t556 = t632 / 0.2e1;
t705 = -t441 * t556 + t562 / 0.2e1;
t595 = qJD(1) * t441;
t552 = t595 / 0.2e1;
t704 = -t412 * t552 - t439 * t556;
t626 = t431 * t439;
t578 = t413 * t626;
t457 = t412 * t595 + t578;
t218 = t511 * t632 + (Icges(7,6) * t431 + (-Icges(7,2) * t424 - t658) * t430) * t412;
t703 = -t430 * t424 * t323 + (-t324 * t430 - t218) * t423;
t418 = sin(t429);
t419 = cos(t429);
t663 = Icges(4,4) * t419;
t515 = -Icges(4,2) * t418 + t663;
t353 = Icges(4,6) * t439 + t441 * t515;
t664 = Icges(4,4) * t418;
t520 = Icges(4,1) * t419 - t664;
t355 = Icges(4,5) * t439 + t441 * t520;
t484 = t353 * t418 - t355 * t419;
t702 = t439 * t484;
t661 = Icges(5,4) * t413;
t513 = -Icges(5,2) * t412 + t661;
t340 = Icges(5,6) * t439 + t441 * t513;
t662 = Icges(5,4) * t412;
t518 = Icges(5,1) * t413 - t662;
t342 = Icges(5,5) * t439 + t441 * t518;
t486 = t340 * t412 - t342 * t413;
t701 = t439 * t486;
t352 = -Icges(4,6) * t441 + t439 * t515;
t354 = -Icges(4,5) * t441 + t439 * t520;
t485 = t352 * t418 - t354 * t419;
t700 = t441 * t485;
t339 = -Icges(5,6) * t441 + t439 * t513;
t341 = -Icges(5,5) * t441 + t439 * t518;
t487 = t339 * t412 - t341 * t413;
t699 = t441 * t487;
t217 = t507 * t632 + (Icges(7,3) * t431 + (-Icges(7,5) * t423 - Icges(7,6) * t424) * t430) * t412;
t648 = t323 * t423;
t698 = -t431 * t648 - t217;
t508 = Icges(6,5) * t440 - Icges(6,6) * t438;
t238 = t508 * t632 + (Icges(6,3) * t431 + (-Icges(6,5) * t438 - Icges(6,6) * t440) * qJD(5)) * t412;
t659 = Icges(6,4) * t440;
t512 = -Icges(6,2) * t438 + t659;
t328 = -Icges(6,6) * t413 + t412 * t512;
t645 = t328 * t438;
t697 = -t431 * t645 - t238;
t258 = t366 * rSges(7,1) + t365 * rSges(7,2) + rSges(7,3) * t633;
t696 = -t439 * t257 - t441 * t258;
t619 = t440 * t441;
t622 = t438 * t439;
t372 = -t413 * t622 - t619;
t620 = t439 * t440;
t373 = t413 * t620 - t621;
t525 = -rSges(6,1) * t373 - rSges(6,2) * t372;
t291 = rSges(6,3) * t634 - t525;
t374 = -t413 * t621 + t620;
t375 = t413 * t619 + t622;
t292 = t375 * rSges(6,1) + t374 * rSges(6,2) + rSges(6,3) * t633;
t695 = -t439 * t291 - t441 * t292;
t437 = -pkin(7) - qJ(2);
t436 = cos(pkin(11));
t414 = t436 * pkin(2) + pkin(1);
t671 = rSges(4,2) * t418;
t673 = rSges(4,1) * t419;
t528 = -t671 + t673;
t474 = -t414 - t528;
t315 = (rSges(4,3) - t437) * t441 + t474 * t439;
t427 = t439 * rSges(4,3);
t602 = t441 * t673 + t427;
t316 = -t437 * t439 + (t414 - t671) * t441 + t602;
t694 = t315 * t441 + t316 * t439;
t509 = Icges(5,5) * t413 - Icges(5,6) * t412;
t337 = -Icges(5,3) * t441 + t439 * t509;
t693 = qJD(1) * t337;
t510 = Icges(4,5) * t419 - Icges(4,6) * t418;
t350 = -Icges(4,3) * t441 + t439 * t510;
t538 = -qJD(5) + t597;
t691 = t439 * t538 + t579;
t580 = t412 * t626;
t689 = t441 * t541 - t580;
t381 = Icges(5,2) * t413 + t662;
t382 = Icges(5,1) * t412 + t661;
t483 = t381 * t412 - t382 * t413;
t688 = qJD(1) * t483 + t509 * t431;
t476 = rSges(3,1) * t436 - rSges(3,2) * sin(pkin(11)) + pkin(1);
t666 = rSges(3,3) + qJ(2);
t336 = t439 * t666 + t441 * t476;
t428 = -pkin(8) + t437;
t601 = t428 - t437;
t396 = pkin(3) * t419 + t414;
t605 = t396 - t414;
t317 = t439 * t605 + t441 * t601;
t687 = 2 * m(4);
t686 = 2 * m(5);
t685 = 2 * m(6);
t684 = 2 * m(7);
t432 = t439 ^ 2;
t433 = t441 ^ 2;
t683 = -t413 / 0.2e1;
t682 = t439 / 0.2e1;
t681 = -t441 / 0.2e1;
t680 = -rSges(6,3) - pkin(9);
t390 = rSges(4,1) * t418 + rSges(4,2) * t419;
t679 = m(4) * t390;
t383 = rSges(5,1) * t412 + rSges(5,2) * t413;
t678 = m(5) * t383;
t677 = pkin(3) * t418;
t676 = pkin(4) * t413;
t672 = rSges(5,1) * t413;
t482 = t439 * t542;
t231 = -t423 * t689 + t424 * t482;
t232 = t423 * t482 + t424 * t689;
t146 = Icges(7,5) * t232 + Icges(7,6) * t231 + Icges(7,3) * t457;
t148 = Icges(7,4) * t232 + Icges(7,2) * t231 + Icges(7,6) * t457;
t150 = Icges(7,1) * t232 + Icges(7,4) * t231 + Icges(7,5) * t457;
t37 = (t431 * t495 - t146) * t413 + (t250 * t431 + (-t252 * t430 + t150) * t424 + (-t254 * t430 - t148) * t423) * t412;
t670 = t37 * t441;
t456 = -t562 + t577;
t145 = Icges(7,5) * t230 + Icges(7,6) * t229 + Icges(7,3) * t456;
t147 = Icges(7,4) * t230 + Icges(7,2) * t229 + Icges(7,6) * t456;
t149 = Icges(7,1) * t230 + Icges(7,4) * t229 + Icges(7,5) * t456;
t38 = (t431 * t494 - t145) * t413 + (t251 * t431 + (-t253 * t430 + t149) * t424 + (-t255 * t430 - t147) * t423) * t412;
t669 = t38 * t439;
t426 = t439 * rSges(5,3);
t539 = -qJD(5) * t413 + qJD(1);
t478 = t539 * t440;
t274 = t439 * t478 + (-t441 * t538 + t580) * t438;
t479 = t539 * t438;
t625 = t431 * t440;
t275 = t538 * t619 + (-t412 * t625 + t479) * t439;
t173 = Icges(6,5) * t275 + Icges(6,6) * t274 + Icges(6,3) * t457;
t175 = Icges(6,4) * t275 + Icges(6,2) * t274 + Icges(6,6) * t457;
t177 = Icges(6,1) * t275 + Icges(6,4) * t274 + Icges(6,5) * t457;
t285 = Icges(6,5) * t373 + Icges(6,6) * t372 + Icges(6,3) * t634;
t287 = Icges(6,4) * t373 + Icges(6,2) * t372 + Icges(6,6) * t634;
t289 = Icges(6,1) * t373 + Icges(6,4) * t372 + Icges(6,5) * t634;
t493 = -t287 * t438 + t289 * t440;
t47 = (t431 * t493 - t173) * t413 + (-t175 * t438 + t177 * t440 + t285 * t431 + (-t287 * t440 - t289 * t438) * qJD(5)) * t412;
t668 = t47 * t441;
t272 = t438 * t691 + t441 * t478;
t273 = -t440 * t691 + t441 * t479;
t172 = Icges(6,5) * t273 + Icges(6,6) * t272 + Icges(6,3) * t456;
t174 = Icges(6,4) * t273 + Icges(6,2) * t272 + Icges(6,6) * t456;
t176 = Icges(6,1) * t273 + Icges(6,4) * t272 + Icges(6,5) * t456;
t286 = Icges(6,5) * t375 + Icges(6,6) * t374 + Icges(6,3) * t633;
t288 = Icges(6,4) * t375 + Icges(6,2) * t374 + Icges(6,6) * t633;
t290 = Icges(6,1) * t375 + Icges(6,4) * t374 + Icges(6,5) * t633;
t492 = -t288 * t438 + t290 * t440;
t48 = (t431 * t492 - t172) * t413 + (-t174 * t438 + t176 * t440 + t286 * t431 + (-t288 * t440 - t290 * t438) * qJD(5)) * t412;
t667 = t48 * t439;
t665 = -rSges(7,3) + t442;
t660 = Icges(6,4) * t438;
t647 = t324 * t424;
t551 = -t383 - t677;
t332 = t551 * t441;
t644 = t332 * t441;
t643 = t352 * t419;
t642 = t353 * t419;
t641 = t354 * t418;
t640 = t355 * t418;
t527 = -rSges(5,2) * t412 + t672;
t361 = t527 * t431;
t639 = t361 * t439;
t638 = t381 * t431;
t637 = t382 * t431;
t635 = t412 * t431;
t631 = t413 * t441;
t239 = t512 * t632 + (Icges(6,6) * t431 + (-Icges(6,2) * t440 - t660) * qJD(5)) * t412;
t623 = t438 * t239;
t523 = t232 * rSges(7,1) + t231 * rSges(7,2);
t152 = rSges(7,3) * t457 + t523;
t618 = t152 * t633 + t257 * t577;
t521 = rSges(7,1) * t424 - rSges(7,2) * t423;
t220 = t521 * t632 + (rSges(7,3) * t431 + (-rSges(7,1) * t423 - rSges(7,2) * t424) * t430) * t412;
t559 = t412 * t591;
t243 = -pkin(5) * t559 + t431 * t454;
t616 = -t220 - t243;
t325 = -rSges(7,3) * t413 + t412 * t521;
t319 = t325 * t596;
t615 = t258 * t635 + t412 * t319;
t524 = rSges(6,1) * t440 - rSges(6,2) * t438;
t242 = t524 * t632 + (rSges(6,3) * t431 + (-rSges(6,1) * t438 - rSges(6,2) * t440) * qJD(5)) * t412;
t530 = pkin(9) * t412 + t676;
t362 = t530 * t431;
t614 = -t242 - t362;
t405 = pkin(4) * t631;
t368 = pkin(9) * t633 + t405;
t411 = pkin(5) * t622;
t475 = t415 * t631 - t442 * t633 + t411;
t280 = t475 - t368;
t612 = -t258 - t280;
t611 = -t292 - t368;
t193 = t413 * t257 + t325 * t634;
t385 = t441 * t396;
t318 = -t414 * t441 - t439 * t601 + t385;
t610 = t439 * t317 + t441 * t318;
t320 = -t412 * t675 + t413 * t674;
t609 = -t320 - t325;
t330 = -rSges(6,3) * t413 + t412 * t524;
t321 = t330 * t596;
t384 = t412 * pkin(4) - t413 * pkin(9);
t369 = t384 * t596;
t608 = t321 + t369;
t607 = -t330 - t384;
t343 = -rSges(5,3) * t441 + t439 * t527;
t344 = rSges(5,1) * t631 - rSges(5,2) * t633 + t426;
t261 = t439 * t343 + t441 * t344;
t367 = t530 * t439;
t606 = t439 * t367 + t441 * t368;
t604 = rSges(5,2) * t562 + rSges(5,3) * t595;
t594 = qJD(3) * t418;
t587 = pkin(3) * t594;
t603 = -t428 * t596 - t439 * t587;
t600 = t432 + t433;
t338 = Icges(5,3) * t439 + t441 * t509;
t599 = qJD(1) * t338;
t351 = Icges(4,3) * t439 + t441 * t510;
t598 = qJD(1) * t351;
t593 = qJD(3) * t419;
t592 = qJD(3) * t439;
t588 = t441 * t671;
t586 = pkin(3) * t593;
t581 = t415 * t635;
t219 = t516 * t632 + (Icges(7,5) * t431 + (-Icges(7,1) * t423 - t657) * t430) * t412;
t575 = t412 * t424 * t219 + t322 * t635 + t632 * t647;
t574 = -t362 + t616;
t517 = Icges(6,1) * t440 - t660;
t240 = t517 * t632 + (Icges(6,5) * t431 + (-Icges(6,1) * t438 - t659) * qJD(5)) * t412;
t327 = -Icges(6,3) * t413 + t412 * t508;
t329 = -Icges(6,5) * t413 + t412 * t517;
t573 = t412 * t440 * t240 + t413 * t329 * t625 + t327 * t635;
t458 = -t413 * t596 - t579;
t469 = t383 * t431;
t571 = t439 * (-t439 * t469 + (t441 * t527 + t426) * qJD(1)) + t441 * (rSges(5,1) * t458 - rSges(5,2) * t577 + t604) + t343 * t595;
t570 = -t368 + t612;
t569 = t273 * rSges(6,1) + t272 * rSges(6,2) + rSges(6,3) * t577;
t410 = t437 * t596;
t568 = t439 * (t595 * t605 + t410 + t603) + t441 * (-qJD(1) * t317 - t441 * t587) + t317 * t595;
t393 = pkin(4) * t580;
t567 = t439 * (pkin(9) * t457 + qJD(1) * t405 - t393) + t441 * (pkin(4) * t458 - pkin(9) * t562 + t394) + t367 * t595;
t309 = t320 * t596;
t566 = t309 + t319 + t369;
t565 = -t384 + t609;
t422 = qJD(2) * t441;
t563 = t422 - t603;
t560 = t418 * t596;
t558 = t634 / 0.2e1;
t557 = t633 / 0.2e1;
t139 = -t285 * t413 + t412 * t493;
t184 = t327 * t634 + t328 * t372 + t329 * t373;
t555 = t139 / 0.2e1 + t184 / 0.2e1;
t140 = -t286 * t413 + t412 * t492;
t185 = t327 * t633 + t328 * t374 + t329 * t375;
t554 = t185 / 0.2e1 + t140 / 0.2e1;
t553 = t596 / 0.2e1;
t550 = -t384 - t677;
t549 = t441 * t609;
t284 = t607 * t441;
t264 = -qJD(1) * t339 - t441 * t638;
t548 = t342 * t431 + t264;
t265 = qJD(1) * t340 - t439 * t638;
t547 = t341 * t431 + t265;
t266 = -qJD(1) * t341 - t441 * t637;
t546 = -t340 * t431 + t266;
t267 = qJD(1) * t342 - t439 * t637;
t545 = t339 * t431 - t267;
t544 = -t439 * t428 + t385;
t543 = -t413 * t415 - t396;
t540 = t441 * t584;
t537 = t413 * t152 + t220 * t634 + t325 * t457;
t161 = t606 - t695;
t532 = -t330 + t550;
t531 = -t362 - t586;
t206 = t565 * t441;
t187 = -t322 * t413 + (t647 - t648) * t412;
t29 = t146 * t633 + t148 * t365 + t150 * t366 + t229 * t252 + t230 * t254 + t250 * t456;
t30 = t145 * t633 + t147 * t365 + t149 * t366 + t229 * t253 + t230 * t255 + t251 * t456;
t114 = t250 * t633 + t252 * t365 + t254 * t366;
t115 = t251 * t633 + t253 * t365 + t255 * t366;
t503 = t114 * t439 + t115 * t441;
t504 = t114 * t441 - t115 * t439;
t63 = t217 * t633 + t218 * t365 + t219 * t366 + t229 * t323 + t230 * t324 + t322 * t456;
t5 = (t431 * t503 - t63) * t413 + (qJD(1) * t504 + t163 * t431 + t29 * t439 + t30 * t441) * t412;
t501 = t128 * t439 + t129 * t441;
t112 = t250 * t634 + t252 * t363 + t254 * t364;
t113 = t251 * t634 + t253 * t363 + t255 * t364;
t505 = t112 * t439 + t113 * t441;
t55 = -t162 * t413 + t412 * t505;
t56 = -t163 * t413 + t412 * t503;
t31 = t146 * t634 + t148 * t363 + t150 * t364 + t231 * t252 + t232 * t254 + t250 * t457;
t32 = t145 * t634 + t147 * t363 + t149 * t364 + t231 * t253 + t232 * t255 + t251 * t457;
t506 = t112 * t441 - t113 * t439;
t64 = t217 * t634 + t218 * t363 + t219 * t364 + t231 * t323 + t232 * t324 + t322 * t457;
t6 = (t431 * t505 - t64) * t413 + (qJD(1) * t506 + t162 * t431 + t31 * t439 + t32 * t441) * t412;
t529 = t5 * t633 + t56 * t577 + t6 * t634 + (-t187 * t413 + t412 * t501) * t635 + t457 * t55;
t526 = rSges(6,1) * t275 + rSges(6,2) * t274;
t519 = Icges(4,1) * t418 + t663;
t514 = Icges(4,2) * t419 + t664;
t380 = Icges(5,5) * t412 + Icges(5,6) * t413;
t502 = t128 * t441 - t129 * t439;
t132 = t285 * t634 + t287 * t372 + t289 * t373;
t133 = t286 * t634 + t288 * t372 + t290 * t373;
t500 = t132 * t441 - t133 * t439;
t499 = t132 * t439 + t133 * t441;
t134 = t285 * t633 + t287 * t374 + t289 * t375;
t135 = t286 * t633 + t288 * t374 + t290 * t375;
t498 = t134 * t441 - t135 * t439;
t497 = t134 * t439 + t135 * t441;
t496 = t139 * t439 + t140 * t441;
t491 = t291 * t441 - t292 * t439;
t473 = -t396 - t527;
t299 = (rSges(5,3) - t428) * t441 + t473 * t439;
t300 = t344 + t544;
t488 = t299 * t441 + t300 * t439;
t480 = t550 + t609;
t477 = -t242 + t531;
t247 = t532 * t441;
t178 = -rSges(6,3) * t562 + t569;
t179 = rSges(6,3) * t457 + t526;
t472 = t441 * t178 + t439 * t179 + t291 * t595 + t567;
t105 = t439 * t279 + t441 * t280 + t606 - t696;
t470 = -qJD(1) * t428 - t587;
t468 = t531 + t616;
t467 = qJD(3) * t390;
t464 = t431 * t380;
t463 = t279 * t441 + t439 * t612;
t462 = qJD(3) * t519;
t461 = qJD(3) * t514;
t460 = qJD(3) * (-Icges(4,5) * t418 - Icges(4,6) * t419);
t196 = t480 * t441;
t459 = t412 * t680 - t396 - t676;
t455 = t412 * t665 + t543;
t15 = qJD(1) * t503 - t29 * t441 + t30 * t439;
t197 = -t337 * t441 - t439 * t487;
t198 = -t338 * t441 - t701;
t199 = t337 * t439 - t699;
t200 = t338 * t439 - t441 * t486;
t39 = t173 * t633 + t175 * t374 + t177 * t375 + t272 * t287 + t273 * t289 + t285 * t456;
t40 = t172 * t633 + t174 * t374 + t176 * t375 + t272 * t288 + t273 * t290 + t286 * t456;
t21 = qJD(1) * t497 - t39 * t441 + t40 * t439;
t262 = -t441 * t464 - t693;
t263 = -t439 * t464 + t599;
t453 = (-t197 * t441 - t500 - t506) * t596 + (-t199 * t441 - t498 - t504) * t595 + (t15 + t21 + t198 * t596 + t200 * t595 + (t200 * qJD(1) + (t265 * t412 - t267 * t413 + t339 * t632 + t341 * t635 - t693) * t441) * t441 + ((t199 + t701) * qJD(1) + (-t263 + t546 * t413 - t548 * t412 + (t338 - t487) * qJD(1)) * t441 + t439 * t262) * t439) * t439;
t171 = -t540 + t393 + (-t581 + (-t431 * t674 - t585) * t413) * t439 + (t441 * t454 + t411) * qJD(1);
t452 = t567 + t706 * t595 + t707 * t441 + (t152 + t171) * t439;
t16 = qJD(1) * t505 - t31 * t441 + t32 * t439;
t451 = t15 * t557 + t16 * t558 + t5 * t682 + t6 * t681 + (qJD(1) * t501 + t669 - t670) * t683 + t55 * t553 + t56 * t552 - t502 * t635 / 0.2e1 + t705 * t504 + t704 * t506;
t450 = rSges(4,2) * t560 + rSges(4,3) * t595 - t441 * t467;
t449 = -t441 * t428 + t439 * t459;
t335 = -t439 * t476 + t441 * t666;
t183 = t187 * t635;
t75 = t412 * t703 + t698 * t413 + t575;
t7 = t183 + (t431 * t501 - t75) * t413 + (qJD(1) * t502 + t37 * t439 + t38 * t441) * t412;
t448 = -t413 * t7 - t56 * t562 + t529;
t447 = t183 + (t37 + t64) * t558 + (t38 + t63) * t557 + t708 * t705 + t709 * t704;
t41 = t173 * t634 + t175 * t372 + t177 * t373 + t274 * t287 + t275 * t289 + t285 * t457;
t42 = t172 * t634 + t174 * t372 + t176 * t373 + t274 * t288 + t275 * t290 + t286 * t457;
t22 = qJD(1) * t499 - t41 * t441 + t42 * t439;
t28 = (t441 * t263 + (t198 + t699) * qJD(1)) * t441 + (t197 * qJD(1) + (-t264 * t412 + t266 * t413 - t340 * t632 - t342 * t635 + t599) * t439 + (-t262 + t545 * t413 + t547 * t412 + (-t337 - t486) * qJD(1)) * t441) * t439;
t446 = (-t16 - t22 - t28) * t441 + t453;
t357 = t513 * t431;
t358 = t518 * t431;
t445 = qJD(1) * t380 + (t358 - t638) * t413 + (-t357 - t637) * t412;
t73 = t238 * t633 + t239 * t374 + t240 * t375 + t272 * t328 + t273 * t329 + t327 * t456;
t10 = (t431 * t497 - t73) * t413 + (qJD(1) * t498 + t185 * t431 + t39 * t439 + t40 * t441) * t412;
t74 = t238 * t634 + t239 * t372 + t240 * t373 + t274 * t328 + t275 * t329 + t327 * t457;
t11 = (t431 * t499 - t74) * t413 + (qJD(1) * t500 + t431 * t184 + t41 * t439 + t42 * t441) * t412;
t67 = -t184 * t413 + t412 * t499;
t68 = -t185 * t413 + t412 * t497;
t444 = t10 * t682 + t11 * t681 + t21 * t557 + t22 * t558 + (qJD(1) * t496 + t667 - t668) * t683 + t451 + t67 * t553 + t68 * t552 + (-t139 * t441 + t140 * t439) * t635 / 0.2e1 + t705 * t498 + t704 * t500;
t443 = -t670 / 0.2e1 + t669 / 0.2e1 - t668 / 0.2e1 + t667 / 0.2e1 + (t412 * t546 + t413 * t548 + t439 * t688 + t445 * t441 + t63 + t73) * t682 + (-t412 * t545 + t413 * t547 + t445 * t439 - t441 * t688 + t64 + t74) * t681 + (t339 * t413 + t341 * t412 - t380 * t441 - t439 * t483 + t139 + t184 - t709) * t553 + (t340 * t413 + t342 * t412 + t380 * t439 - t441 * t483 + t140 + t185 - t708) * t552;
t421 = qJD(2) * t439;
t400 = pkin(3) * t560;
t379 = t528 * qJD(3);
t360 = -t588 + t602;
t359 = -rSges(4,3) * t441 + t439 * t528;
t331 = t551 * t439;
t314 = -qJD(1) * t336 + t422;
t313 = qJD(1) * t335 + t421;
t294 = t439 * t460 + t598;
t293 = -qJD(1) * t350 + t441 * t460;
t283 = t607 * t439;
t246 = t532 * t439;
t236 = t257 * t633;
t235 = -t383 * t595 - t639 + (-t418 * t595 - t419 * t592) * pkin(3);
t234 = t383 * t596 + t400 + (-t361 - t586) * t441;
t228 = t410 + t422 + t390 * t592 + (t441 * t474 - t427) * qJD(1);
t227 = t421 + (-t437 * t441 + (-t414 - t673) * t439) * qJD(1) + t450;
t212 = t351 * t439 - t441 * t484;
t211 = t350 * t439 - t700;
t210 = -t351 * t441 - t702;
t209 = -t350 * t441 - t439 * t485;
t208 = t544 - t611;
t207 = t449 + t525;
t205 = t565 * t439;
t204 = -t292 * t413 - t330 * t633;
t203 = t291 * t413 + t330 * t634;
t202 = t383 * t626 + (t441 * t473 - t426) * qJD(1) + t563;
t201 = t421 + (-t396 - t672) * t596 + (-t469 + t470) * t441 + t604;
t195 = t480 * t439;
t194 = -t258 * t413 - t325 * t633;
t192 = -t327 * t413 + (t329 * t440 - t645) * t412;
t191 = t475 + t544 + t258;
t190 = (pkin(5) * t438 - t428) * t441 + t455 * t439 + t522;
t189 = t491 * t412;
t188 = t192 * t635;
t186 = -t258 * t634 + t236;
t182 = t261 + t610;
t167 = qJD(1) * t284 + t439 * t614;
t166 = t441 * t614 + t608;
t158 = qJD(1) * t247 + t439 * t477;
t157 = t441 * t477 + t400 + t608;
t138 = -t344 * t596 + t571;
t123 = t412 * t549 + t413 * t612;
t122 = t279 * t413 + t320 * t634 + t193;
t111 = t459 * t595 + t578 * t680 + t393 - t526 + t563;
t110 = t394 + t421 + (-pkin(4) * t635 - t587) * t441 + t449 * qJD(1) + t569;
t109 = t161 + t610;
t108 = t412 * t463 + t236;
t107 = qJD(1) * t206 + t439 * t574;
t106 = t441 * t574 + t566;
t104 = qJD(1) * t196 + t439 * t468;
t103 = t441 * t468 + t400 + t566;
t102 = t540 + (t581 + (t431 * t665 + t585) * t413) * t439 + (t441 * t455 - t411) * qJD(1) - t523 + t563;
t101 = t421 + (-rSges(7,3) * t412 + t543) * t596 + (t413 * t471 + t470 - t581) * t441 + t564 + t572;
t100 = (t330 * t626 + t179) * t413 + (t242 * t439 - t291 * t431 + t330 * t595) * t412;
t99 = (-t330 * t624 - t178) * t413 + (-t242 * t441 + t292 * t431 + t321) * t412;
t98 = t105 + t610;
t96 = (-t318 - t344) * t596 + t568 + t571;
t91 = -t257 * t635 + t537;
t90 = -t220 * t633 + (-t325 * t624 - t151) * t413 + t615;
t89 = t697 * t413 + (-t623 + (-t328 * t440 - t329 * t438) * qJD(5)) * t412 + t573;
t70 = t491 * t632 + (qJD(1) * t695 - t178 * t439 + t179 * t441) * t412;
t69 = t596 * t611 + t472;
t57 = -t258 * t578 + (qJD(1) * t696 - t151 * t439) * t412 + t618;
t49 = (-t318 + t611) * t596 + t472 + t568;
t44 = (t320 * t626 + t171) * t413 + (t243 * t439 + t320 * t595 - t431 * t706) * t412 + t537;
t43 = (t431 * t549 - t707) * t413 + (t280 * t431 + t441 * t616 + t309) * t412 + t615;
t26 = t570 * t596 + t452;
t25 = t463 * t632 + (t171 * t441 - t707 * t439 + (-t439 * t706 + t441 * t612) * qJD(1)) * t412 + t618;
t24 = (-t318 + t570) * t596 + t452 + t568;
t1 = [(t101 * t191 + t102 * t190) * t684 + (t110 * t208 + t111 * t207) * t685 + (t201 * t300 + t202 * t299) * t686 + (t227 * t316 + t228 * t315) * t687 - t329 * t559 + t382 * t632 + t573 + t575 - t381 * t635 + 0.2e1 * m(3) * (t313 * t336 + t314 * t335) + (-t514 + t520) * t594 + (t515 + t519) * t593 + (t357 + t697 + t698) * t413 + (-t328 * t590 + t358 - t623 + t703) * t412; m(7) * (-t101 * t441 + t102 * t439 + (t190 * t441 + t191 * t439) * qJD(1)) + m(6) * (-t110 * t441 + t111 * t439 + (t207 * t441 + t208 * t439) * qJD(1)) + m(5) * (qJD(1) * t488 - t201 * t441 + t202 * t439) + m(4) * (qJD(1) * t694 - t227 * t441 + t228 * t439) + m(3) * (-t313 * t441 + t314 * t439 + (t335 * t441 + t336 * t439) * qJD(1)); 0; (t432 / 0.2e1 + t433 / 0.2e1) * t510 * qJD(3) + m(4) * ((-t227 * t439 - t228 * t441) * t390 - t694 * t379) + m(5) * (t201 * t331 + t202 * t332 + t234 * t299 + t235 * t300) + m(6) * (t110 * t246 + t111 * t247 + t157 * t207 + t158 * t208) + m(7) * (t101 * t195 + t102 * t196 + t103 * t190 + t104 * t191) + t443 + (-qJD(3) * t485 + (qJD(1) * t353 - t439 * t461) * t419 + (qJD(1) * t355 - t439 * t462) * t418) * t681 + ((-t316 * t679 + t642 / 0.2e1 + t640 / 0.2e1) * t441 + (t315 * t679 + t643 / 0.2e1 + t641 / 0.2e1) * t439) * qJD(1) + (-qJD(3) * t484 + (-qJD(1) * t352 - t441 * t461) * t419 + (-qJD(1) * t354 - t441 * t462) * t418) * t682; m(5) * (t234 * t439 - t235 * t441 + (t331 * t439 + t644) * qJD(1)) + m(6) * (t157 * t439 - t158 * t441 + (t246 * t439 + t247 * t441) * qJD(1)) + m(7) * (t103 * t439 - t104 * t441 + (t195 * t439 + t196 * t441) * qJD(1)); ((t359 * t439 + t360 * t441) * ((qJD(1) * t359 + t450) * t441 + (-t439 * t467 + (-t360 - t588 + t427) * qJD(1)) * t439) + t600 * t390 * t379) * t687 + (t103 * t196 + t104 * t195 + t24 * t98) * t684 + (t109 * t49 + t157 * t247 + t158 * t246) * t685 + (t182 * t96 + t234 * t332 + t235 * t331) * t686 - t441 * t16 - t441 * t22 - t441 * t28 + (-t209 * t441 + t210 * t439) * t596 + (-t211 * t441 + t212 * t439) * t595 + t453 + t439 * ((t439 * t293 + (t211 + t702) * qJD(1)) * t439 + (t212 * qJD(1) + (t352 * t593 + t354 * t594) * t441 + (-t294 + (-t640 - t642) * qJD(3) + (t351 - t485) * qJD(1)) * t439) * t441) - t441 * ((t441 * t294 + (t210 + t700) * qJD(1)) * t441 + (t209 * qJD(1) + (-t353 * t593 - t355 * t594 + t598) * t439 + (-t293 + (t641 + t643) * qJD(3) - t484 * qJD(1)) * t441) * t439); t443 + m(6) * (t110 * t283 + t111 * t284 + t166 * t207 + t167 * t208) + m(7) * (t101 * t205 + t102 * t206 + t106 * t190 + t107 * t191) - m(5) * t488 * t361 + (-t201 * t439 - t202 * t441 + (t299 * t439 - t300 * t441) * qJD(1)) * t678; m(6) * (t166 * t439 - t167 * t441 + (t283 * t439 + t284 * t441) * qJD(1)) + m(7) * (t106 * t439 - t107 * t441 + (t205 * t439 + t206 * t441) * qJD(1)); t446 + (-t234 * t441 - t235 * t439 + (-t331 * t441 + t332 * t439) * qJD(1)) * t678 + m(7) * (t103 * t206 + t104 * t205 + t105 * t24 + t106 * t196 + t107 * t195 + t26 * t98) + m(6) * (t109 * t69 + t157 * t284 + t158 * t283 + t161 * t49 + t166 * t247 + t167 * t246) + m(5) * (t138 * t182 + t261 * t96 - t331 * t639 - t361 * t644); (t105 * t26 + t106 * t206 + t107 * t205) * t684 + (t161 * t69 + t166 * t284 + t167 * t283) * t685 + (t361 * t383 * t600 + t138 * t261) * t686 + t446; t188 + ((t48 / 0.2e1 + t73 / 0.2e1) * t441 + (t47 / 0.2e1 + t74 / 0.2e1) * t439 + (-t439 * t554 + t441 * t555) * qJD(1)) * t412 + (-t75 - t89 + (t439 * t555 + t441 * t554) * t431) * t413 + m(6) * (t100 * t207 + t110 * t204 + t111 * t203 + t208 * t99) + m(7) * (t101 * t123 + t102 * t122 + t190 * t44 + t191 * t43) + t447; m(6) * (t100 * t439 - t441 * t99 + (t203 * t441 + t204 * t439) * qJD(1)) + m(7) * (-t43 * t441 + t439 * t44 + (t122 * t441 + t123 * t439) * qJD(1)); t444 + m(6) * (t100 * t247 + t109 * t70 + t157 * t203 + t158 * t204 + t189 * t49 + t246 * t99) + m(7) * (t103 * t122 + t104 * t123 + t108 * t24 + t195 * t43 + t196 * t44 + t25 * t98); t444 + m(6) * (t100 * t284 + t161 * t70 + t166 * t203 + t167 * t204 + t189 * t69 + t283 * t99) + m(7) * (t105 * t25 + t106 * t122 + t107 * t123 + t108 * t26 + t205 * t43 + t206 * t44); (t108 * t25 + t122 * t44 + t123 * t43) * t684 + (t100 * t203 + t189 * t70 + t204 * t99) * t685 + (t89 * t413 - t188 - t7 + (-t413 * t496 + t439 * t67 + t441 * t68) * t431) * t413 + (t441 * t10 + t439 * t11 + t496 * t635 + (-t192 * t431 - t47 * t439 - t48 * t441) * t413 + ((-t139 * t413 + t67) * t441 + (t140 * t413 - t56 - t68) * t439) * qJD(1)) * t412 + t529; -t75 * t413 + t447 + m(7) * (t101 * t194 + t102 * t193 + t190 * t91 + t191 * t90); m(7) * (t439 * t91 - t441 * t90 + (t193 * t441 + t194 * t439) * qJD(1)); t451 + m(7) * (t103 * t193 + t104 * t194 + t186 * t24 + t195 * t90 + t196 * t91 + t57 * t98); t451 + m(7) * (t105 * t57 + t106 * t193 + t107 * t194 + t186 * t26 + t205 * t90 + t206 * t91); m(7) * (t108 * t57 + t122 * t91 + t123 * t90 + t186 * t25 + t193 * t44 + t194 * t43) + t448; (t186 * t57 + t193 * t91 + t194 * t90) * t684 + t448;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
