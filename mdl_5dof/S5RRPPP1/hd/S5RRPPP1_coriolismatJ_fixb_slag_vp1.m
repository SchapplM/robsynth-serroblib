% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:23
% EndTime: 2019-12-31 19:24:00
% DurationCPUTime: 22.66s
% Computational Cost: add. (22700->637), mult. (62492->880), div. (0->0), fcn. (72626->8), ass. (0->339)
t756 = Icges(5,5) + Icges(6,4);
t770 = -Icges(4,6) + t756;
t487 = sin(qJ(1));
t486 = sin(qJ(2));
t488 = cos(qJ(2));
t623 = cos(pkin(8));
t624 = cos(pkin(5));
t513 = t624 * t623;
t501 = t488 * t513;
t622 = sin(pkin(8));
t433 = t486 * t622 - t501;
t512 = t624 * t622;
t434 = t486 * t623 + t488 * t512;
t485 = sin(pkin(5));
t599 = t485 * t488;
t447 = pkin(2) * t486 - qJ(3) * t599;
t566 = pkin(3) * t434 + qJ(4) * t433 + t447;
t651 = rSges(6,1) + pkin(4);
t740 = rSges(6,3) + qJ(5);
t521 = rSges(6,2) * t433 + t740 * t434 - t651 * t599 + t566;
t192 = t521 * t487;
t489 = cos(qJ(1));
t194 = t521 * t489;
t543 = -rSges(5,1) * t599 - rSges(5,2) * t434 + rSges(5,3) * t433 + t566;
t233 = t543 * t487;
t235 = t543 * t489;
t570 = rSges(4,1) * t434 - rSges(4,2) * t433 - rSges(4,3) * t599 + t447;
t306 = t570 * t487;
t308 = t570 * t489;
t598 = t486 * t487;
t437 = t485 * t598 - t489 * t624;
t597 = t486 * t489;
t438 = t485 * t597 + t487 * t624;
t600 = t485 * t486;
t448 = pkin(2) * t488 + qJ(3) * t600;
t500 = pkin(1) + t448;
t533 = t624 * qJ(3);
t491 = (t533 + pkin(7)) * t487 + t500 * t489;
t435 = t486 * t513 + t488 * t622;
t390 = -t485 * t487 * t623 + t435 * t489;
t502 = t486 * t512;
t535 = t489 * t623;
t537 = t487 * t622;
t391 = t485 * t537 + t488 * t535 - t489 * t502;
t496 = t391 * rSges(4,1) - t390 * rSges(4,2) + t438 * rSges(4,3);
t212 = t491 + t496;
t595 = t488 * t489;
t545 = t485 * t595;
t596 = t487 * t488;
t546 = t485 * t596;
t388 = t435 * t487 + t485 * t535;
t436 = t488 * t623 - t502;
t534 = t489 * t622;
t389 = t436 * t487 - t485 * t534;
t285 = t389 * rSges(4,1) - t388 * rSges(4,2) + t437 * rSges(4,3);
t471 = t489 * t533;
t484 = t489 * pkin(7);
t494 = -t487 * t500 + t471 + t484;
t755 = -t285 + t494;
t579 = t212 * t546 + t545 * t755;
t378 = t388 * qJ(4);
t493 = -t378 + t494;
t510 = -t437 * rSges(5,1) - t388 * rSges(5,3);
t650 = rSges(5,2) - pkin(3);
t162 = t389 * t650 + t493 + t510;
t379 = t390 * qJ(4);
t490 = t379 + t491;
t509 = t438 * rSges(5,1) + t390 * rSges(5,3);
t163 = -t391 * t650 + t490 + t509;
t588 = t162 * t545 + t163 * t546;
t551 = pkin(3) + t740;
t763 = t388 * rSges(6,2) + t437 * t651;
t134 = -t389 * t551 + t493 - t763;
t497 = t390 * rSges(6,2) + t391 * qJ(5) + t438 * t651;
t135 = (rSges(6,3) + pkin(3)) * t391 + t490 + t497;
t592 = t134 * t545 + t135 * t546;
t691 = m(6) / 0.2e1;
t692 = m(5) / 0.2e1;
t694 = m(4) / 0.2e1;
t548 = (-t438 * t192 + t437 * t194 + t592) * t691 + (-t438 * t233 + t437 * t235 + t588) * t692 + (-t438 * t306 + t437 * t308 + t579) * t694;
t520 = -pkin(2) * t597 + qJ(3) * t545;
t418 = t486 * t534 - t489 * t501;
t419 = t434 * t489;
t529 = -t419 * pkin(3) - t418 * qJ(4);
t495 = t520 + t529;
t498 = -t418 * rSges(6,2) - t740 * t419 + t651 * t545;
t172 = t495 + t498;
t417 = t434 * t487;
t416 = t486 * t537 - t487 * t501;
t406 = t416 * qJ(4);
t477 = pkin(2) * t598;
t559 = t406 + t477;
t627 = rSges(6,2) * t416;
t173 = t627 + (-qJ(3) - t651) * t546 + t551 * t417 + t559;
t626 = rSges(5,3) * t416;
t206 = t626 + (-rSges(5,1) - qJ(3)) * t546 - t650 * t417 + t559;
t539 = rSges(5,1) * t545 + t419 * rSges(5,2) - t418 * rSges(5,3);
t207 = t495 + t539;
t511 = rSges(4,1) * t417 - rSges(4,2) * t416;
t294 = t477 + (-rSges(4,3) - qJ(3)) * t546 + t511;
t540 = -t419 * rSges(4,1) + t418 * rSges(4,2) + rSges(4,3) * t545;
t295 = t520 + t540;
t549 = (t172 * t437 + t173 * t438 + t592) * t691 + (t206 * t438 + t207 * t437 + t588) * t692 + (t294 * t438 + t295 * t437 + t579) * t694;
t4 = t549 - t548;
t769 = t4 * qJD(1);
t505 = -t134 * t489 - t135 * t487;
t632 = (-t390 * t192 + t388 * t194 + t433 * t505) * t691 + ((-t162 * t489 - t163 * t487) * t433 - t390 * t233 + t388 * t235) * t692;
t633 = (-t134 * t418 - t135 * t416 + t172 * t388 + t173 * t390) * t691 + (-t162 * t418 - t163 * t416 + t206 * t390 + t207 * t388) * t692;
t7 = t633 - t632;
t768 = t7 * qJD(1);
t757 = Icges(4,5) + Icges(6,5);
t751 = Icges(5,4) - t757;
t767 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t753 = Icges(5,1) + Icges(4,3) + Icges(6,1);
t752 = Icges(4,4) + Icges(5,6) - Icges(6,6);
t764 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t719 = -t752 * t388 + t767 * t389 - t437 * t751;
t718 = t764 * t388 - t752 * t389 + t437 * t770;
t717 = t388 * t770 - t389 * t751 + t753 * t437;
t716 = t390 * t770 - t391 * t751 + t753 * t438;
t715 = t764 * t390 - t752 * t391 + t438 * t770;
t714 = -t752 * t390 + t767 * t391 - t438 * t751;
t760 = m(6) * (-t389 * t438 + t437 * t391);
t247 = -t760 / 0.2e1;
t249 = t760 / 0.2e1;
t761 = m(6) * (-t388 * t391 + t389 * t390);
t198 = -t761 / 0.2e1;
t199 = t761 / 0.2e1;
t302 = t487 * t388 + t390 * t489;
t303 = t487 * t389 + t391 * t489;
t741 = m(6) * (t434 * t302 - t433 * t303);
t152 = -t741 / 0.2e1;
t153 = t741 / 0.2e1;
t747 = m(5) + m(6);
t762 = t747 / 0.2e1;
t255 = -t388 * t438 + t437 * t390;
t552 = m(5) / 0.4e1 + m(6) / 0.4e1;
t733 = -0.2e1 * t552;
t759 = t255 * t733;
t758 = t747 * t255;
t748 = (-t718 * t390 - t719 * t391 - t717 * t438) * t489 + (t715 * t390 + t714 * t391 + t716 * t438) * t487;
t541 = t389 * rSges(5,2) + t510;
t347 = Icges(6,5) * t600 + Icges(6,6) * t435 + Icges(6,3) * t436;
t348 = -Icges(6,5) * t599 + Icges(6,6) * t433 + Icges(6,3) * t434;
t349 = Icges(5,5) * t600 - Icges(5,6) * t436 + Icges(5,3) * t435;
t350 = -Icges(5,5) * t599 - Icges(5,6) * t434 + Icges(5,3) * t433;
t351 = Icges(6,4) * t600 + Icges(6,2) * t435 + Icges(6,6) * t436;
t352 = -Icges(6,4) * t599 + Icges(6,2) * t433 + Icges(6,6) * t434;
t353 = Icges(5,4) * t600 - Icges(5,2) * t436 + Icges(5,6) * t435;
t354 = -Icges(5,4) * t599 - Icges(5,2) * t434 + Icges(5,6) * t433;
t355 = Icges(6,1) * t600 + Icges(6,4) * t435 + Icges(6,5) * t436;
t357 = Icges(5,1) * t600 - Icges(5,4) * t436 + Icges(5,5) * t435;
t360 = Icges(4,5) * t436 - Icges(4,6) * t435 + Icges(4,3) * t600;
t361 = Icges(4,4) * t434 - Icges(4,2) * t433 - Icges(4,6) * t599;
t362 = Icges(4,4) * t436 - Icges(4,2) * t435 + Icges(4,6) * t600;
t363 = Icges(4,1) * t434 - Icges(4,4) * t433 - Icges(4,5) * t599;
t364 = Icges(4,1) * t436 - Icges(4,4) * t435 + Icges(4,5) * t600;
t620 = Icges(3,4) * t486;
t458 = Icges(3,2) * t488 + t620;
t461 = Icges(3,1) * t488 - t620;
t734 = -(t486 * (-t753 * t599 / 0.2e1 + (-Icges(5,4) / 0.2e1 + t757 / 0.2e1) * t434 + (-Icges(4,6) / 0.2e1 + t756 / 0.2e1) * t433) - t488 * (t357 / 0.2e1 + t360 / 0.2e1 + t355 / 0.2e1)) * t485 + (t458 / 0.2e1 - t461 / 0.2e1) * t486 + (t353 / 0.2e1 - t364 / 0.2e1 - t347 / 0.2e1) * t434 + (t361 / 0.2e1 - t350 / 0.2e1 - t352 / 0.2e1) * t435 + (t362 / 0.2e1 - t349 / 0.2e1 - t351 / 0.2e1) * t433 - (t363 / 0.2e1 + t348 / 0.2e1 - t354 / 0.2e1) * t436;
t560 = t437 * t546 + t438 * t545;
t725 = (m(4) / 0.4e1 + t552) * (-t485 ^ 2 * t486 * t488 + t560);
t724 = t552 * (-t388 * t416 - t390 * t418 + t433 * t435);
t530 = -t389 * pkin(3) - t378;
t482 = Icges(3,4) * t488;
t459 = -Icges(3,2) * t486 + t482;
t460 = Icges(3,1) * t486 + t482;
t713 = t752 * t418 - t419 * t767 - t751 * t545;
t712 = -t418 * t764 + t752 * t419 + t545 * t770;
t711 = t416 * t770 - t751 * t417 - t753 * t546;
t710 = -t418 * t770 + t751 * t419 + t753 * t545;
t709 = -t416 * t764 + t752 * t417 + t546 * t770;
t708 = t752 * t416 - t417 * t767 - t751 * t546;
t707 = t347 + t364 - t353;
t706 = t349 - t362 + t351;
t705 = -t350 + t361 - t352;
t704 = t360 + t355 + t357;
t703 = -t363 - t348 + t354;
t571 = t740 * t389 + t763;
t700 = t487 ^ 2;
t699 = t489 ^ 2;
t698 = 0.2e1 * qJD(1);
t697 = 0.4e1 * qJD(1);
t696 = 2 * qJD(2);
t564 = t487 * (t448 * t487 - t471) + t489 * (pkin(2) * t595 + qJ(3) * t438);
t137 = t487 * t285 + t489 * t496 + t564;
t375 = t487 * t437 + t438 * t489;
t572 = -t306 * t546 - t308 * t545;
t687 = m(4) * (t137 * t375 + t572);
t525 = -t487 * t530 + t489 * (t391 * pkin(3) + t379) + t564;
t103 = -t487 * t541 + t489 * (-t391 * rSges(5,2) + t509) + t525;
t576 = -t233 * t546 - t235 * t545;
t679 = m(5) * (t103 * t375 + t576);
t677 = m(5) * (t162 * t206 + t163 * t207);
t148 = t163 * t438;
t98 = -t162 * t437 + t148;
t676 = m(5) * t98;
t670 = m(6) * (-t134 * t419 - t135 * t417 + t172 * t389 + t173 * t391);
t668 = m(6) * (-t391 * t192 + t389 * t194 + t434 * t505);
t582 = -t192 * t546 - t194 * t545;
t88 = (t391 * rSges(6,3) + t497) * t489 + t571 * t487 + t525;
t666 = m(6) * (t375 * t88 + t582);
t664 = m(6) * (t134 * t173 + t135 * t172);
t131 = t135 * t438;
t87 = -t134 * t437 + t131;
t663 = m(6) * t87;
t655 = t487 / 0.2e1;
t653 = -t489 / 0.2e1;
t628 = rSges(3,1) * t488;
t538 = pkin(1) + t628;
t553 = rSges(3,2) * t598 + t489 * rSges(3,3);
t396 = -t487 * t538 + t484 + t553;
t476 = rSges(3,2) * t597;
t397 = -t476 + t538 * t489 + (rSges(3,3) + pkin(7)) * t487;
t463 = rSges(3,1) * t486 + rSges(3,2) * t488;
t445 = t463 * t487;
t446 = t463 * t489;
t649 = m(3) * (t396 * t445 - t397 * t446);
t648 = m(4) * (t212 * t295 + t294 * t755);
t647 = m(4) * (t212 * t438 - t437 * t755);
t646 = m(6) * (-t388 * t417 - t389 * t416 - t390 * t419 - t391 * t418 + t433 * t436 + t434 * t435);
t567 = t389 * t546 + t391 * t545;
t644 = m(6) * (-t417 * t437 - t419 * t438 + (t434 * t486 - t436 * t488) * t485 + t567);
t610 = t375 * t434;
t640 = m(6) * (-t303 * t599 - t610);
t639 = m(6) * (t567 + t610);
t631 = m(6) * qJD(1);
t630 = m(6) * qJD(2);
t629 = m(6) * qJD(5);
t125 = t135 * t390;
t126 = t135 * t391;
t140 = t163 * t390;
t611 = t375 * t433;
t424 = Icges(3,4) * t596 - Icges(3,2) * t598 - Icges(3,6) * t489;
t605 = t424 * t486;
t79 = -t134 * t388 + t125;
t80 = -t134 * t389 + t126;
t492 = t494 + t530;
t136 = t492 - t571;
t591 = t134 - t136;
t94 = -t162 * t388 + t140;
t568 = t388 * t546 + t390 * t545;
t145 = -t416 * t437 - t418 * t438 + (t433 * t486 - t435 * t488) * t485 + t568;
t590 = t145 * t762;
t164 = t492 + t541;
t587 = t162 - t164;
t179 = -t302 * t599 - t611;
t584 = t179 * t762;
t188 = t568 + t611;
t581 = t188 * t762;
t575 = -t758 / 0.2e1;
t574 = t758 / 0.2e1;
t569 = -rSges(4,1) * t436 + rSges(4,2) * t435 - rSges(4,3) * t600 - t448;
t565 = -pkin(3) * t436 - qJ(4) * t435 - t448;
t422 = Icges(3,5) * t596 - Icges(3,6) * t598 - Icges(3,3) * t489;
t474 = Icges(3,4) * t598;
t426 = Icges(3,1) * t596 - Icges(3,5) * t489 - t474;
t562 = -t487 * t422 - t426 * t595;
t507 = Icges(3,5) * t488 - Icges(3,6) * t486;
t423 = Icges(3,3) * t487 + t489 * t507;
t427 = Icges(3,5) * t487 + t461 * t489;
t561 = t487 * t423 + t427 * t595;
t558 = t487 * (qJ(3) * t546 - t477) + t489 * t520;
t557 = t460 * t487 + t424;
t425 = Icges(3,6) * t487 + t459 * t489;
t556 = -t460 * t489 - t425;
t555 = -Icges(3,2) * t596 + t426 - t474;
t554 = -t458 * t489 + t427;
t187 = -t389 * t417 - t391 * t419 + t434 * t436;
t547 = t187 * t629;
t544 = -rSges(5,1) * t600 + rSges(5,2) * t436 - rSges(5,3) * t435 + t565;
t532 = t556 * t487;
t531 = t554 * t487;
t400 = t427 * t596;
t528 = t489 * t423 - t400;
t527 = t425 * t486 - t422;
t523 = t487 * (-pkin(3) * t417 - t406) + t489 * t529 + t558;
t522 = -rSges(6,2) * t435 - t740 * t436 - t651 * t600 + t565;
t504 = t192 * t487 + t194 * t489;
t19 = m(5) * (t103 * t302 + (t233 * t487 + t235 * t489) * t433) + m(6) * (t88 * t302 + t433 * t504);
t32 = m(5) * t94 + m(6) * t79;
t506 = -Icges(3,5) * t486 - Icges(3,6) * t488;
t465 = -rSges(3,2) * t486 + t628;
t440 = t506 * t489;
t439 = t506 * t487;
t309 = t569 * t489;
t307 = t569 * t487;
t293 = -t425 * t597 + t561;
t292 = -t424 * t597 - t562;
t291 = -t425 * t598 - t528;
t236 = t544 * t489;
t234 = t544 * t487;
t195 = t522 * t489;
t193 = t522 * t487;
t184 = t639 / 0.2e1;
t175 = t640 / 0.2e1;
t171 = -t292 * t489 + t293 * t487;
t170 = -(-t487 * (-t426 * t488 + t605) - t489 * t422) * t489 + t291 * t487;
t169 = t487 * (rSges(4,3) * t546 - t511) + t489 * t540 + t558;
t154 = 0.4e1 * t725;
t143 = t644 / 0.2e1;
t122 = t487 * (rSges(5,1) * t546 + rSges(5,2) * t417 - t626) + t489 * t539 + t523;
t119 = t646 / 0.2e1;
t118 = 0.4e1 * t724;
t105 = t498 * t489 + (-t417 * t740 + t546 * t651 - t627) * t487 + t523;
t101 = t249 + t247;
t100 = 0.2e1 * t249;
t99 = 0.2e1 * t247;
t91 = (t291 - t400 + (t423 + t605) * t489 + t562) * t489 + t561 * t487;
t90 = (t489 * t527 + t293 - t561) * t489 + (t487 * t527 + t292 + t528) * t487;
t86 = t175 + t184 - t644 / 0.2e1;
t85 = t175 + t143 - t639 / 0.2e1;
t84 = t184 + t143 - t640 / 0.2e1;
t83 = 0.2e1 * t199;
t82 = t198 + t199;
t81 = 0.2e1 * t198;
t70 = 0.2e1 * t153 + t119;
t69 = t152 + t153 - t646 / 0.2e1;
t68 = 0.2e1 * t152 + t119;
t52 = t574 - t759;
t51 = t574 + t575;
t50 = t575 + t759;
t44 = t668 / 0.2e1;
t42 = t88 * t303 + t434 * t504;
t36 = t670 / 0.2e1;
t31 = t188 * t733 + t584 + t590;
t30 = t145 * t733 + t581 + t584;
t29 = t179 * t733 + t581 + t590;
t21 = t647 + t663 + t676;
t18 = t666 + t679 + t687;
t11 = t36 - t668 / 0.2e1;
t10 = t44 + t36;
t9 = t44 - t670 / 0.2e1;
t8 = (t459 / 0.2e1 + t460 / 0.2e1) * t488 + t649 + t648 + t677 + t664 - t734;
t6 = t632 + t633;
t3 = t548 + t549;
t1 = (t171 / 0.2e1 - t91 / 0.2e1) * t489 + (t90 / 0.2e1 + t170 / 0.2e1) * t487;
t2 = [(-m(5) * t587 * t163 / 0.4e1 - m(6) * t591 * t135 / 0.4e1) * t697 + t8 * qJD(2) + t21 * qJD(3) + t32 * qJD(4) + t80 * t629, t8 * qJD(1) + t3 * qJD(3) + t6 * qJD(4) + t10 * qJD(5) + (m(3) * ((-t396 * t489 - t397 * t487) * t465 + (-t445 * t489 + t446 * t487) * t463) / 0.2e1 + (t212 * t307 - t294 * t308 - t295 * t306 + t309 * t755) * t694 + (t162 * t236 + t163 * t234 - t206 * t235 - t207 * t233) * t692 + (t134 * t195 + t135 * t193 - t172 * t192 - t173 * t194) * t691) * t696 + ((t91 + t748) * t489 / 0.2e1 + (t390 * t706 + t391 * t707 + t418 * t705 + t419 * t703 + t433 * t712 + t434 * t713 + t435 * t715 + t436 * t714 + t438 * t704 + t486 * t556 + t488 * t554) * t655 - (t170 + t90) * t487 / 0.2e1 + (t388 * t706 + t389 * t707 + t416 * t705 + t417 * t703 + t433 * t709 + t434 * t708 + t435 * t718 + t436 * t719 + t437 * t704 - t486 * t557 + t488 * t555 + t171 + t748) * t653 + (t700 / 0.2e1 + t699 / 0.2e1) * t507 + ((t486 * t716 - t488 * t710) * t655 + (t486 * t717 + t488 * t711) * t653) * t485) * qJD(2), t21 * qJD(1) + t3 * qJD(2) + t51 * qJD(4) + t101 * qJD(5), t32 * qJD(1) + t6 * qJD(2) + t51 * qJD(3) + t82 * qJD(5), t10 * qJD(2) + t101 * qJD(3) + t82 * qJD(4) + t631 * t80; t1 * qJD(2) - t4 * qJD(3) - t7 * qJD(4) + t9 * qJD(5) + (t591 * t192 * t691 + t587 * t233 * t692) * t698 + (-t649 / 0.4e1 - t648 / 0.4e1 - t677 / 0.4e1 - t664 / 0.4e1) * t697 + (-(t459 + t460) * t488 / 0.2e1 + t734) * qJD(1), t1 * qJD(1) + t18 * qJD(3) + t19 * qJD(4) + t42 * t629 + (m(6) * (t105 * t88 - t192 * t193 - t194 * t195) + m(5) * (t103 * t122 - t233 * t234 - t235 * t236) + m(4) * (t137 * t169 - t306 * t307 - t308 * t309) + m(3) * ((t487 * (rSges(3,1) * t596 - t553) + t489 * (rSges(3,1) * t595 + t487 * rSges(3,3) - t476)) * (-t487 * t445 - t446 * t489) + (t699 + t700) * t465 * t463) + (t700 * t440 + (-t531 * t486 + t532 * t488 + (t555 * t486 + t557 * t488) * t489 - t717 * t545 + t711 * t438 + t719 * t419 + t718 * t418 - t708 * t391 - t709 * t390) * t489 + (t390 * t712 + t391 * t713 - t418 * t715 - t419 * t714 + t438 * t710 - t439 * t489 + t545 * t716) * t487) * t655 + (t699 * t439 + (-t388 * t709 - t389 * t708 + t416 * t718 + t417 * t719 + t437 * t711 - t546 * t717) * t489 + (-t489 * t440 + (t489 * t557 + t532) * t488 + (t489 * t555 - t531) * t486 + t716 * t546 + t710 * t437 - t714 * t417 - t715 * t416 + t713 * t389 + t712 * t388) * t487) * t653) * qJD(2), -t769 + t18 * qJD(2) + t30 * qJD(4) + t86 * qJD(5) + (-0.4e1 * t725 + 0.2e1 * (t691 + t692 + t694) * (-t375 * t599 + t560)) * qJD(3), t19 * qJD(2) + t30 * qJD(3) - 0.4e1 * qJD(4) * t724 + t69 * qJD(5) - t768, t9 * qJD(1) + t86 * qJD(3) + t69 * qJD(4) + t42 * t630 - t547; t4 * qJD(2) + t52 * qJD(4) + t100 * qJD(5) + (-t663 / 0.4e1 - t647 / 0.4e1 - t676 / 0.4e1) * t697 + ((t136 * t437 - t131 + t87) * t691 + (t164 * t437 - t148 + t98) * t692) * t698, t769 + t154 * qJD(3) + t31 * qJD(4) + t85 * qJD(5) + 0.4e1 * (-t666 / 0.4e1 - t679 / 0.4e1 - t687 / 0.4e1) * qJD(2) + ((t193 * t437 + t195 * t438 + (-t105 * t488 + t486 * t88) * t485 + t582) * t691 + (t234 * t437 + t236 * t438 + (t103 * t486 - t122 * t488) * t485 + t576) * t692 + (t307 * t437 + t309 * t438 + (t137 * t486 - t169 * t488) * t485 + t572) * t694) * t696, qJD(2) * t154, qJD(1) * t52 + qJD(2) * t31, qJD(1) * t100 + qJD(2) * t85; (m(6) * (t136 * t388 - t125 + t79) + m(5) * (t164 * t388 - t140 + t94) - t32) * qJD(1) + t7 * qJD(2) + t50 * qJD(3) + t81 * qJD(5), t768 + (m(6) * (t105 * t433 + t192 * t416 + t193 * t388 + t194 * t418 + t195 * t390 + t435 * t88) + m(5) * (t103 * t435 + t122 * t433 + t233 * t416 + t234 * t388 + t235 * t418 + t236 * t390) - t19) * qJD(2) + t29 * qJD(3) + t118 * qJD(4) + t68 * qJD(5), qJD(1) * t50 + qJD(2) * t29, qJD(2) * t118, qJD(1) * t81 + qJD(2) * t68; (t136 * t389 - t126) * t631 + t11 * qJD(2) + t99 * qJD(3) + t83 * qJD(4), t11 * qJD(1) + (t105 * t434 + t192 * t417 + t193 * t389 + t194 * t419 + t195 * t391 + t436 * t88 - t42) * t630 + t84 * qJD(3) + t70 * qJD(4) + t547, qJD(1) * t99 + qJD(2) * t84, qJD(1) * t83 + qJD(2) * t70, t187 * t630;];
Cq = t2;
