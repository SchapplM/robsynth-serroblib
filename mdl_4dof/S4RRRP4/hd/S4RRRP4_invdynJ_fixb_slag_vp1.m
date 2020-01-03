% Calculate vector of inverse dynamics joint torques for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:10
% EndTime: 2019-12-31 17:15:40
% DurationCPUTime: 26.54s
% Computational Cost: add. (11504->661), mult. (15913->846), div. (0->0), fcn. (12614->6), ass. (0->367)
t734 = Icges(4,4) + Icges(5,4);
t735 = Icges(4,1) + Icges(5,1);
t731 = Icges(4,2) + Icges(5,2);
t728 = Icges(4,5) + Icges(5,5);
t727 = Icges(4,6) + Icges(5,6);
t371 = qJ(2) + qJ(3);
t346 = cos(t371);
t733 = t734 * t346;
t345 = sin(t371);
t732 = t734 * t345;
t730 = t735 * t346 - t732;
t729 = -t731 * t345 + t733;
t726 = Icges(4,3) + Icges(5,3);
t718 = t731 * t346 + t732;
t722 = t735 * t345 + t733;
t373 = sin(qJ(1));
t559 = t345 * t373;
t725 = t734 * t559;
t375 = cos(qJ(1));
t724 = t727 * t375;
t723 = t728 * t375;
t556 = t346 * t373;
t639 = t734 * t556 - t731 * t559 - t724;
t638 = t727 * t373 + t729 * t375;
t707 = t728 * t373 + t730 * t375;
t689 = t735 * t556 - t723 - t725;
t706 = -t727 * t345 + t728 * t346;
t719 = t728 * t345 + t727 * t346;
t717 = t726 * t375;
t716 = t729 + t722;
t715 = t730 - t718;
t714 = -t718 * t375 + t707;
t713 = -t722 * t375 - t638;
t712 = t722 * t373 + t639;
t708 = -t728 * t556 + t727 * t559 + t717;
t711 = t726 * t373 + t706 * t375;
t710 = t718 * t345 - t722 * t346;
t705 = t639 * t345 - t689 * t346;
t709 = t707 * t556;
t368 = qJD(2) + qJD(3);
t704 = t715 * t368;
t703 = t716 * t368;
t702 = -qJD(1) * t707 + t368 * t712;
t701 = t713 * t368 + (-t373 * t730 + t723) * qJD(1);
t292 = t368 * t373;
t700 = qJD(1) * t638 - t718 * t292 + t689 * t368;
t699 = t714 * t368 + (-t373 * t729 + t724) * qJD(1);
t698 = t719 * t375;
t697 = t719 * t373;
t696 = t638 * t345;
t661 = -t710 * t373 - t698;
t660 = -t710 * t375 + t697;
t376 = -pkin(6) - pkin(5);
t367 = -qJ(4) + t376;
t695 = t367 - rSges(5,3);
t694 = t705 * t373;
t693 = t711 * t375 - t709;
t692 = t711 * qJD(1);
t555 = t346 * t375;
t634 = -t711 * t373 - t707 * t555;
t691 = t708 * t373 - t689 * t555;
t690 = t708 * t375;
t665 = t690 - t694;
t664 = -t638 * t559 - t693;
t558 = t345 * t375;
t663 = -t639 * t558 - t691;
t662 = -t638 * t558 - t634;
t688 = qJD(1) * t719 - t703 * t345 + t704 * t346;
t687 = -t699 * t345 + t701 * t346 + t692;
t686 = t708 * qJD(1) + t700 * t345 + t702 * t346;
t293 = t368 * t375;
t685 = qJD(1) * t715 + t292 * t713 + t293 * t712;
t684 = (t731 * t556 - t689 + t725) * t293 + t714 * t292 + t716 * qJD(1);
t683 = t710 * qJD(1) + t706 * t368;
t682 = t705 * qJD(1) - t697 * t368 + t692;
t681 = -t698 * t368 + (-t707 * t346 - t706 * t373 + t696 + t717) * qJD(1);
t680 = t660 * qJD(1);
t335 = pkin(3) * t346;
t374 = cos(qJ(2));
t364 = t374 * pkin(2);
t336 = t364 + pkin(1);
t495 = t335 + t336;
t679 = -rSges(5,1) * t556 + rSges(5,2) * t559 - t373 * t495 - t695 * t375;
t357 = t373 * rSges(5,3);
t678 = rSges(5,1) * t555 - rSges(5,2) * t558 + t375 * t495 + t357;
t372 = sin(qJ(2));
t599 = pkin(2) * qJD(2);
t494 = t372 * t599;
t609 = pkin(3) * t345;
t247 = -t368 * t609 - t494;
t343 = qJD(4) * t373;
t502 = qJD(1) * t373;
t487 = t345 * t502;
t501 = qJD(1) * t375;
t677 = rSges(5,2) * t487 + rSges(5,3) * t501 + t375 * t247 + t343;
t676 = t661 * qJD(1);
t675 = -t373 * t682 + t375 * t686;
t674 = t373 * t681 + t375 * t687;
t673 = t373 * t686 + t375 * t682;
t672 = t373 * t687 - t375 * t681;
t671 = t664 * t292 - t293 * t665 + t676;
t670 = t292 * t662 - t293 * t663 + t680;
t669 = t373 * t683 + t375 * t688;
t668 = t373 * t688 - t375 * t683;
t667 = t702 * t345 - t700 * t346;
t666 = t701 * t345 + t699 * t346;
t659 = t345 * t689 + t346 * t639;
t658 = t707 * t345 + t638 * t346;
t657 = rSges(3,2) * t372;
t402 = -t293 * t345 - t346 * t502;
t499 = qJD(2) * t375;
t483 = t372 * t499;
t448 = pkin(2) * t483;
t450 = t336 - t495;
t492 = t293 * t346;
t507 = t376 - t367;
t595 = rSges(5,1) * t402 - rSges(5,2) * t492 + t448 + (t373 * t450 + t375 * t507) * qJD(1) + t677;
t598 = t346 * rSges(5,2);
t274 = rSges(5,1) * t345 + t598;
t241 = t274 * t373;
t344 = qJD(4) * t375;
t452 = -t373 * t247 + t344;
t320 = t373 * t494;
t509 = t376 * t502 + t320;
t552 = t373 * t367;
t331 = t346 * rSges(5,1);
t600 = rSges(5,2) * t345;
t637 = t331 - t600;
t594 = -t368 * t241 - t452 + t509 + (t357 - t552 + (-t450 + t637) * t375) * qJD(1);
t340 = t375 * t376;
t511 = t373 * t336 + t340;
t542 = t511 + t679;
t314 = t375 * t336;
t541 = t373 * t507 - t314 + t678;
t646 = -t345 * t684 + t346 * t685;
t645 = t706 * qJD(1) - t698 * t292 + t697 * t293;
t275 = rSges(4,1) * t345 + rSges(4,2) * t346;
t242 = t275 * t373;
t244 = t275 * t375;
t332 = t346 * rSges(4,1);
t636 = -rSges(4,2) * t345 + t332;
t201 = rSges(4,1) * t556 - rSges(4,2) * t559 - t375 * rSges(4,3);
t358 = t373 * rSges(4,3);
t203 = rSges(4,1) * t555 - rSges(4,2) * t558 + t358;
t365 = t375 * pkin(5);
t305 = pkin(1) * t373 - t365;
t186 = t305 - t511;
t363 = t373 * pkin(5);
t306 = t375 * pkin(1) + t363;
t451 = -t373 * t376 + t314;
t187 = t451 - t306;
t500 = qJD(2) * t373;
t535 = -t186 * t500 + t187 * t499;
t70 = t201 * t292 + t203 * t293 + t535;
t404 = -t275 * t293 - t448;
t530 = t186 - t305;
t488 = -t201 + t530;
t71 = qJD(1) * t488 + t404;
t529 = -t187 - t203;
t72 = -t275 * t292 - t320 + (t306 - t529) * qJD(1);
t644 = -t71 * (qJD(1) * t242 - t293 * t636) - t70 * (-t292 * t242 - t244 * t293) - t72 * (-qJD(1) * t244 - t292 * t636);
t497 = qJD(1) * qJD(2);
t333 = t373 * t497;
t496 = qJD(1) * qJD(3);
t207 = t373 * t496 + t333 + (-qJDD(2) - qJDD(3)) * t375;
t218 = t637 * t368;
t286 = -qJDD(2) * t375 + t333;
t549 = t374 * qJD(2) ^ 2;
t610 = pkin(2) * t372;
t414 = -pkin(2) * t375 * t549 + t286 * t610;
t446 = t530 + t542;
t604 = pkin(1) - t336;
t146 = (-t375 * t604 - t363) * qJD(1) - t509;
t284 = t306 * qJD(1);
t543 = -t146 - t284;
t557 = t346 * t368;
t12 = qJDD(4) * t373 + t207 * t274 - t218 * t293 + (t207 * t345 - t293 * t557) * pkin(3) + t446 * qJDD(1) + (t344 + t543 - t594) * qJD(1) + t414;
t643 = t12 - g(1);
t285 = qJDD(2) * t373 + t375 * t497;
t206 = qJDD(3) * t373 + t375 * t496 + t285;
t341 = pkin(5) * t501;
t145 = -t448 - t341 + (t373 * t604 - t340) * qJD(1);
t518 = qJD(1) * (-pkin(1) * t502 + t341) + qJDD(1) * t306;
t394 = qJD(1) * t145 + qJDD(1) * t187 + (-t285 * t372 - t373 * t549) * pkin(2) + t518;
t13 = -qJDD(4) * t375 - t206 * t274 - t218 * t292 + t541 * qJDD(1) + (-t206 * t345 - t292 * t557) * pkin(3) + (t343 + t595) * qJD(1) + t394;
t642 = t13 - g(2);
t474 = -t274 - t609;
t489 = -t187 - t541;
t61 = -t320 - t344 + t474 * t292 + (t306 - t489) * qJD(1);
t641 = qJD(1) * t61 + t12;
t288 = qJD(1) * t305;
t640 = qJD(1) * t186 - t288;
t356 = Icges(3,4) * t374;
t431 = -Icges(3,2) * t372 + t356;
t299 = Icges(3,1) * t372 + t356;
t635 = t696 + t708;
t633 = g(1) * t375 + g(2) * t373;
t551 = t373 * t374;
t554 = t372 * t373;
t578 = Icges(3,3) * t375;
t220 = Icges(3,5) * t551 - Icges(3,6) * t554 - t578;
t325 = Icges(3,4) * t554;
t587 = Icges(3,5) * t375;
t224 = Icges(3,1) * t551 - t325 - t587;
t581 = Icges(3,6) * t375;
t222 = Icges(3,4) * t551 - Icges(3,2) * t554 - t581;
t568 = t222 * t372;
t422 = -t224 * t374 + t568;
t85 = -t375 * t220 - t373 * t422;
t296 = Icges(3,5) * t374 - Icges(3,6) * t372;
t295 = Icges(3,5) * t372 + Icges(3,6) * t374;
t405 = qJD(2) * t295;
t590 = Icges(3,4) * t372;
t300 = Icges(3,1) * t374 - t590;
t225 = Icges(3,5) * t373 + t300 * t375;
t223 = Icges(3,6) * t373 + t375 * t431;
t567 = t223 * t372;
t421 = -t225 * t374 + t567;
t628 = -t375 * t405 + (-t296 * t373 + t421 + t578) * qJD(1);
t221 = Icges(3,3) * t373 + t296 * t375;
t504 = qJD(1) * t221;
t627 = qJD(1) * t422 - t373 * t405 + t504;
t297 = Icges(3,2) * t374 + t590;
t416 = t297 * t372 - t299 * t374;
t624 = t416 * qJD(1) + t296 * qJD(2);
t623 = t373 * (-t297 * t375 + t225) - t375 * (-Icges(3,2) * t551 + t224 - t325);
t620 = t206 / 0.2e1;
t619 = t207 / 0.2e1;
t618 = t285 / 0.2e1;
t617 = t286 / 0.2e1;
t616 = -t292 / 0.2e1;
t615 = t292 / 0.2e1;
t614 = -t293 / 0.2e1;
t613 = t293 / 0.2e1;
t612 = t373 / 0.2e1;
t611 = -t375 / 0.2e1;
t606 = -qJD(1) / 0.2e1;
t605 = qJD(1) / 0.2e1;
t603 = rSges(3,1) * t374;
t602 = rSges(3,2) * t374;
t359 = t373 * rSges(3,3);
t597 = t71 * t275;
t596 = qJDD(1) / 0.2e1;
t43 = -t292 * t542 + t293 * t541 + t535;
t575 = qJD(1) * t43;
t301 = rSges(3,1) * t372 + t602;
t484 = t301 * t499;
t508 = rSges(3,2) * t554 + t375 * rSges(3,3);
t227 = rSges(3,1) * t551 - t508;
t521 = -t227 - t305;
t121 = qJD(1) * t521 - t484;
t573 = t121 * t373;
t572 = t121 * t375;
t550 = t374 * t375;
t553 = t372 * t375;
t228 = rSges(3,1) * t550 - rSges(3,2) * t553 + t359;
t168 = t228 + t306;
t122 = qJD(1) * t168 - t301 * t500;
t258 = t301 * t375;
t571 = t122 * t258;
t562 = t292 * t346;
t561 = t295 * t373;
t560 = t295 * t375;
t534 = -t373 * t186 + t375 * t187;
t533 = t373 * t201 + t375 * t203;
t532 = -t373 * t220 - t224 * t550;
t531 = t373 * t221 + t225 * t550;
t294 = pkin(3) * t487;
t520 = t274 * t502 + t294;
t514 = rSges(4,2) * t487 + rSges(4,3) * t501;
t513 = -t297 + t300;
t512 = t299 + t431;
t510 = rSges(3,3) * t501 + t502 * t657;
t503 = qJD(1) * t296;
t125 = -t373 * t416 - t560;
t498 = t125 * qJD(1);
t493 = t374 * t599;
t116 = rSges(4,1) * t402 - rSges(4,2) * t492 + t514;
t118 = -t368 * t242 + (t375 * t636 + t358) * qJD(1);
t491 = t375 * t116 + t373 * t118 + t201 * t501;
t490 = t375 * t145 + t373 * t146 - t186 * t501;
t485 = t372 * t501;
t481 = -pkin(1) - t603;
t480 = t502 / 0.2e1;
t479 = t501 / 0.2e1;
t478 = -t500 / 0.2e1;
t477 = t500 / 0.2e1;
t476 = -t499 / 0.2e1;
t475 = t499 / 0.2e1;
t403 = -t275 - t610;
t473 = t372 * (-t373 ^ 2 - t375 ^ 2);
t243 = t274 * t375;
t462 = -t292 * t241 - t243 * t293;
t181 = t225 * t551;
t461 = t375 * t221 - t181;
t454 = -qJD(1) * t243 - t292 * t637;
t453 = -t220 + t567;
t447 = -t542 * t373 + t541 * t375;
t445 = -pkin(3) * t557 - t218;
t219 = t636 * t368;
t444 = -t219 - t493;
t443 = t637 + t335;
t289 = -t609 - t610;
t442 = -pkin(3) * t492 + qJD(1) * t241 - t293 * t637;
t304 = rSges(2,1) * t375 - rSges(2,2) * t373;
t302 = rSges(2,1) * t373 + rSges(2,2) * t375;
t303 = t603 - t657;
t124 = t223 * t374 + t225 * t372;
t406 = qJD(2) * t297;
t138 = -t375 * t406 + (-t373 * t431 + t581) * qJD(1);
t407 = qJD(2) * t299;
t140 = -t375 * t407 + (-t300 * t373 + t587) * qJD(1);
t384 = -qJD(2) * t124 - t138 * t372 + t140 * t374 + t504;
t123 = t222 * t374 + t224 * t372;
t139 = qJD(1) * t223 - t373 * t406;
t141 = qJD(1) * t225 - t373 * t407;
t385 = qJD(1) * t220 - qJD(2) * t123 - t139 * t372 + t141 * t374;
t439 = -(t373 * t627 + t385 * t375) * t375 + (t373 * t628 + t384 * t375) * t373;
t438 = -(t385 * t373 - t375 * t627) * t375 + (t384 * t373 - t375 * t628) * t373;
t437 = -t373 * t72 - t375 * t71;
t86 = -t223 * t554 - t461;
t436 = t373 * t86 - t375 * t85;
t87 = -t222 * t553 - t532;
t88 = -t223 * t553 + t531;
t435 = t373 * t88 - t375 * t87;
t428 = -t122 * t373 - t572;
t142 = -t499 * t602 + (-t374 * t502 - t483) * rSges(3,1) + t510;
t257 = t301 * t373;
t143 = -qJD(2) * t257 + (t303 * t375 + t359) * qJD(1);
t427 = t142 * t375 + t143 * t373;
t420 = t227 * t373 + t228 * t375;
t417 = t297 * t374 + t299 * t372;
t415 = -t495 - t331;
t413 = t594 * t373 + t595 * t375 - t542 * t501;
t412 = -t274 + t289;
t411 = t145 * t499 + t146 * t500 - t285 * t186 - t187 * t286;
t401 = t445 - t493;
t398 = t222 * t375 - t223 * t373;
t397 = (-t372 * t512 + t374 * t513) * qJD(1);
t396 = t293 * t474 + t343 - t448;
t279 = t431 * qJD(2);
t280 = t300 * qJD(2);
t383 = qJD(1) * t295 - qJD(2) * t417 - t279 * t372 + t280 * t374;
t382 = -t372 * t623 + t398 * t374;
t381 = (t662 * t373 - t663 * t375) * t620 + (t664 * t373 - t665 * t375) * t619 + (t645 * t373 + t646 * t375) * t616 + (t675 * t375 + t674 * t373 + (t663 * t373 + t662 * t375) * qJD(1)) * t615 + (t673 * t375 + t672 * t373 + (t665 * t373 + t664 * t375) * qJD(1)) * t614 + (t646 * t373 - t645 * t375) * t613 + (t669 * qJD(1) + t660 * qJDD(1) + t662 * t206 + t663 * t207 + t674 * t292 + t675 * t293) * t612 + (t668 * qJD(1) + t661 * qJDD(1) + t664 * t206 + t665 * t207 + t672 * t292 + t673 * t293) * t611 + (t345 * t685 + t346 * t684) * t606 + (t667 * t375 + t666 * t373 + (t659 * t373 + t658 * t375) * qJD(1)) * t605 + (t658 * t373 - t659 * t375) * t596 + t671 * t480 + t670 * t479;
t281 = t303 * qJD(2);
t260 = t375 * t289;
t259 = t373 * t289;
t209 = pkin(2) * t553 + t260;
t208 = pkin(2) * t554 + t259;
t126 = -t375 * t416 + t561;
t120 = t126 * qJD(1);
t119 = t420 * qJD(2);
t67 = qJD(1) * t142 + qJDD(1) * t228 - t281 * t500 - t285 * t301 + t518;
t66 = -t281 * t499 + t286 * t301 + t521 * qJDD(1) + (-t143 - t284) * qJD(1);
t65 = t383 * t373 - t375 * t624;
t64 = t373 * t624 + t383 * t375;
t63 = -qJD(2) * t421 + t138 * t374 + t140 * t372;
t62 = -t422 * qJD(2) + t139 * t374 + t141 * t372;
t60 = qJD(1) * t446 + t396;
t49 = qJD(2) * t435 + t120;
t48 = qJD(2) * t436 + t498;
t42 = qJD(1) * t116 + qJDD(1) * t203 - t206 * t275 - t219 * t292 + t394;
t41 = t207 * t275 - t219 * t293 + t488 * qJDD(1) + (-t118 + t543) * qJD(1) + t414;
t18 = t116 * t293 + t118 * t292 + t201 * t206 - t203 * t207 + t411;
t9 = -t206 * t542 - t207 * t541 + t292 * t594 + t293 * t595 + t411;
t1 = [(t120 + ((t86 - t181 + (t221 + t568) * t375 + t532) * t375 + t531 * t373) * qJD(2)) * t475 - m(2) * (-g(1) * t302 + g(2) * t304) + (t126 + t124) * t618 + (t125 + t123) * t617 + (((t373 * t638 + t375 * t639) * t345 + t664 + t691 + t693) * t293 + (-t689 * t556 + (t373 * t639 - t375 * t638) * t345 - t634 + t665 - t690) * t292 + t680) * t613 + (-t498 + ((t375 * t453 - t531 + t88) * t375 + (t373 * t453 + t461 + t87) * t373) * qJD(2) + t48) * t478 + (t63 + t64) * t477 + (t62 + t65 + t49) * t476 + (-qJD(2) * t416 + t279 * t374 + t280 * t372 + t704 * t345 + t703 * t346) * qJD(1) + (t60 * t452 + t61 * t677 + (t241 * t60 - t243 * t61) * t368 + ((t60 * (t415 + t600) - t61 * t367) * t375 + (t61 * t415 + t60 * t695) * t373) * qJD(1) - (qJD(1) * t542 + t396 - t60 + t640) * t61 + t642 * (-t552 + t678) + t643 * t679) * m(5) + (t71 * t509 + t72 * (-t448 + t514) + (-t244 * t72 + t373 * t597) * t368 + ((-t71 * rSges(4,3) + t72 * (-t336 - t332)) * t373 + (t71 * (-t336 - t636) - t72 * t376) * t375) * qJD(1) - (-qJD(1) * t201 + t404 + t640 - t71) * t72 + (t42 - g(2)) * (t203 + t451) + (t41 - g(1)) * (-t201 - t511)) * m(4) + (t122 * (t341 + t510) + (t301 * t573 - t571) * qJD(2) + ((-pkin(1) - t303) * t572 + (t121 * (-rSges(3,3) - pkin(5)) + t122 * t481) * t373) * qJD(1) - (-qJD(1) * t227 - t121 - t288 - t484) * t122 + (t67 - g(2)) * t168 + (t66 - g(1)) * (t481 * t373 + t365 + t508)) * m(3) + (t658 + t660) * t620 + (t659 + t661) * t619 + ((t375 * t635 + t634 + t662 - t694) * t293 + ((t705 + t711) * t375 + t635 * t373 + t663 - t709) * t292 + t671 - t676) * t616 + (t666 + t669) * t615 + (m(2) * (t302 ^ 2 + t304 ^ 2) + t417 + Icges(2,3) + t718 * t346 + t722 * t345) * qJDD(1) + (-t667 + t668 + t670) * t614; t49 * t479 + ((-t499 * t561 - t503) * t375 + (t397 + (t375 * t560 + t382) * qJD(2)) * t373) * t475 + ((t372 * t513 + t374 * t512) * qJD(1) + (t398 * t372 + t374 * t623) * qJD(2)) * t606 + (t373 * t63 - t375 * t62 + (t123 * t373 + t124 * t375) * qJD(1)) * t605 + (qJD(1) * t64 + qJD(2) * t439 + qJDD(1) * t126 + t285 * t88 + t286 * t87) * t612 + (qJD(1) * t65 + qJD(2) * t438 + qJDD(1) * t125 + t285 * t86 + t286 * t85) * t611 + ((t85 * t373 + t86 * t375) * qJD(1) + t438) * t476 + ((-t500 * t560 + t503) * t373 + (t397 + (t373 * t561 + t382) * qJD(2)) * t375) * t478 + ((t87 * t373 + t88 * t375) * qJD(1) + t439) * t477 + t381 + t48 * t480 + (-t123 * t375 + t124 * t373) * t596 + t435 * t618 + t436 * t617 + (-g(1) * (t260 - t243) - g(2) * (t259 - t241) - g(3) * (t364 + t443) + t60 * t520 + t9 * (t447 + t534) + t43 * (t413 + t490) + (t401 * t60 + t412 * t641) * t375 + (t13 * t412 + t401 * t61 + t489 * t575) * t373 - t60 * (-qJD(1) * t208 + t442) - t61 * (-pkin(3) * t562 + qJD(1) * t209 + t454) - t43 * (t208 * t292 + t209 * t293 + t462) - (-t61 * t485 + ((-t373 * t61 - t375 * t60) * t374 + t43 * t473) * qJD(2)) * pkin(2)) * m(5) + (t18 * (t533 + t534) + t70 * (t490 + t491) + (t444 * t71 + (qJD(1) * t72 + t41) * t403) * t375 + (t42 * t403 + t72 * t444 + (t529 * t70 + t597) * qJD(1)) * t373 - (-t72 * t485 + (t374 * t437 + t473 * t70) * qJD(2)) * pkin(2) - g(3) * (t636 + t364) - t633 * t403 + t644) * m(4) + ((qJD(2) * t427 + t227 * t285 - t228 * t286) * t420 + t119 * ((t227 * t375 - t228 * t373) * qJD(1) + t427) + t428 * t281 + (-t67 * t373 - t66 * t375 + (-t122 * t375 + t573) * qJD(1)) * t301 - (t121 * t257 - t571) * qJD(1) - (t119 * (-t257 * t373 - t258 * t375) + t428 * t303) * qJD(2) + g(1) * t258 + g(2) * t257 - g(3) * t303) * m(3); t381 + (t9 * t447 + t43 * t413 + (t445 * t61 - t541 * t575) * t373 - t61 * t454 - t43 * t462 - (-t61 * t562 + (-t61 * t501 + t43 * (-t292 * t373 - t293 * t375)) * t345) * pkin(3) - g(3) * t443 + (t13 * t373 + t375 * t641) * t474 + (t375 * t445 - t294 - t442 + t520) * t60 - t633 * (-t598 + (-rSges(5,1) - pkin(3)) * t345)) * m(5) + (t18 * t533 + t70 * (-t203 * t502 + t491) + t437 * t219 + (-t42 * t373 - t41 * t375 + (t373 * t71 - t375 * t72) * qJD(1)) * t275 + g(1) * t244 + g(2) * t242 - g(3) * t636 + t644) * m(4); ((t292 * t43 - t642) * t375 + (-t293 * t43 + t643) * t373) * m(5);];
tau = t1;
