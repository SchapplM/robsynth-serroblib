% Calculate vector of inverse dynamics joint torques for
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:35
% EndTime: 2022-01-23 09:30:15
% DurationCPUTime: 34.39s
% Computational Cost: add. (19387->718), mult. (16589->885), div. (0->0), fcn. (12988->8), ass. (0->393)
t752 = Icges(5,6) + Icges(6,6);
t761 = Icges(5,1) + Icges(6,1);
t762 = Icges(5,4) + Icges(6,4);
t758 = Icges(5,5) + Icges(6,5);
t760 = Icges(5,2) + Icges(6,2);
t383 = qJ(1) + pkin(8);
t375 = cos(t383);
t759 = t752 * t375;
t757 = Icges(5,3) + Icges(6,3);
t374 = sin(t383);
t384 = qJ(3) + qJ(4);
t377 = cos(t384);
t584 = t374 * t377;
t376 = sin(t384);
t585 = t374 * t376;
t732 = t762 * t584 - t760 * t585 - t759;
t618 = Icges(6,4) * t376;
t619 = Icges(5,4) * t376;
t756 = t761 * t377 - t618 - t619;
t755 = t762 * t585;
t754 = t762 * t377;
t753 = t758 * t375;
t730 = t761 * t584 - t753 - t755;
t728 = -t752 * t376 + t758 * t377;
t282 = Icges(6,2) * t377 + t618;
t284 = Icges(5,2) * t377 + t619;
t751 = t282 + t284;
t727 = t761 * t376 + t754;
t749 = -t760 * t376 + t754;
t729 = t758 * t374 + t756 * t375;
t748 = t758 * t376 + t752 * t377;
t747 = -t284 + t756;
t746 = t732 * t376;
t745 = t757 * t375;
t744 = t749 + t727;
t733 = -t758 * t584 + t752 * t585 + t745;
t700 = t757 * t374 + t728 * t375;
t731 = t752 * t374 + t749 * t375;
t743 = t376 * t751 - t727 * t377;
t724 = -t730 * t377 + t746;
t742 = t729 * t584;
t382 = qJD(3) + qJD(4);
t741 = t727 * t382;
t588 = t282 * t382;
t740 = t747 * t382 - t588;
t739 = t744 * t382;
t719 = t748 * t375;
t718 = t748 * t374;
t738 = -t284 * t375 + t729;
t702 = -t743 * t374 - t719;
t701 = -t743 * t375 + t718;
t737 = t724 * t374;
t736 = t731 * t585 - t742;
t579 = t375 * t377;
t735 = t700 * t374 + t729 * t579;
t734 = t733 * t374 - t730 * t579;
t684 = t733 * t375 - t737;
t683 = -t700 * t375 - t736;
t580 = t375 * t376;
t682 = -t732 * t580 - t734;
t681 = -t731 * t580 + t735;
t726 = t748 * qJD(1) - t739 * t376 + t740 * t377;
t520 = rSges(5,1) * t584;
t389 = -pkin(7) - pkin(6);
t341 = t375 * t389;
t387 = cos(qJ(3));
t379 = t387 * pkin(3);
t370 = t379 + pkin(2);
t546 = t374 * t370 + t341;
t386 = sin(qJ(1));
t640 = pkin(1) * t386;
t725 = -t520 - t640 - t546;
t723 = -t729 * qJD(1) + t741 * t374 + t732 * t382;
t722 = -t731 * t382 - t741 * t375 + (-t374 * t756 + t753) * qJD(1);
t272 = t374 * t382;
t721 = t731 * qJD(1) - t272 * t284 - t374 * t588 + t730 * t382;
t720 = -t375 * t588 + t738 * t382 + (-t374 * t749 + t759) * qJD(1);
t717 = t731 * t376;
t716 = t743 * qJD(1) + t728 * t382;
t388 = cos(qJ(1));
t380 = t388 * pkin(1);
t715 = t701 * qJD(1);
t381 = -qJ(5) + t389;
t714 = t381 - rSges(6,3);
t274 = rSges(3,1) * t374 + rSges(3,2) * t375;
t260 = -t274 - t640;
t713 = t702 * qJD(1);
t712 = t700 * qJD(1);
t696 = t717 + t733;
t711 = t696 * t375 - t735 - t737;
t525 = qJD(1) * qJD(3);
t336 = t374 * t525;
t524 = qJD(1) * qJD(4);
t190 = t374 * t524 + t336 + (-qJDD(3) - qJDD(4)) * t375;
t367 = t377 * rSges(6,1);
t629 = rSges(6,2) * t376;
t666 = t367 - t629;
t248 = t666 * t382;
t273 = t375 * t382;
t627 = t377 * rSges(6,2);
t290 = rSges(6,1) * t376 + t627;
t345 = qJD(5) * t375;
t268 = -qJDD(3) * t375 + t336;
t391 = qJD(1) ^ 2;
t523 = t391 * t380;
t574 = t387 * qJD(3) ^ 2;
t385 = sin(qJ(3));
t639 = pkin(3) * t385;
t410 = -pkin(3) * t375 * t574 + t268 * t639 - t523;
t363 = t375 * pkin(6);
t276 = pkin(2) * t374 - t363;
t185 = t276 - t546;
t500 = -t276 - t640;
t471 = t185 + t500;
t369 = pkin(4) * t377;
t521 = t369 + t370;
t689 = -rSges(6,1) * t584 + rSges(6,2) * t585 - t374 * t521 - t714 * t375;
t571 = -t546 - t689;
t429 = t471 - t571;
t362 = t374 * pkin(6);
t530 = qJD(3) * t385;
t519 = pkin(3) * t530;
t317 = t374 * t519;
t535 = qJD(1) * t374;
t544 = t389 * t535 + t317;
t633 = pkin(2) - t370;
t150 = (-t375 * t633 - t362) * qJD(1) - t544;
t277 = t375 * pkin(2) + t362;
t266 = t277 * qJD(1);
t572 = -t150 - t266;
t575 = t377 * t382;
t233 = t290 * t374;
t576 = t376 * t382;
t269 = -pkin(4) * t576 - t519;
t355 = t374 * rSges(6,3);
t477 = t370 - t521;
t583 = t374 * t381;
t624 = -t382 * t233 + t374 * t269 - t345 + t544 + (t355 - t583 + (-t477 + t666) * t375) * qJD(1);
t12 = qJDD(5) * t374 + t190 * t290 - t248 * t273 + (t190 * t376 - t273 * t575) * pkin(4) + t429 * qJDD(1) + (t345 + t572 - t624) * qJD(1) + t410;
t710 = t12 - g(1);
t267 = qJDD(3) * t374 + t375 * t525;
t189 = qJDD(4) * t374 + t375 * t524 + t267;
t344 = qJD(5) * t374;
t534 = qJD(1) * t375;
t342 = pkin(6) * t534;
t509 = t375 * t530;
t475 = pkin(3) * t509;
t149 = -t475 - t342 + (t374 * t633 - t341) * qJD(1);
t310 = t375 * t370;
t478 = -t374 * t389 + t310;
t186 = t478 - t277;
t473 = qJDD(1) * t380 - t391 * t640;
t430 = qJD(1) * (-pkin(2) * t535 + t342) + qJDD(1) * t277 + t473;
t395 = qJD(1) * t149 + qJDD(1) * t186 + (-t267 * t385 - t374 * t574) * pkin(3) + t430;
t540 = t389 - t381;
t691 = rSges(6,1) * t579 - rSges(6,2) * t580 + t375 * t521 + t355;
t570 = t374 * t540 - t310 + t691;
t517 = t375 * t576;
t415 = -t377 * t535 - t517;
t514 = t376 * t535;
t516 = t375 * t575;
t690 = rSges(6,3) * t534 + (t514 - t516) * rSges(6,2) + t375 * t269 + t344;
t625 = rSges(6,1) * t415 + t475 + (t374 * t477 + t375 * t540) * qJD(1) + t690;
t13 = -qJDD(5) * t375 - t189 * t290 - t272 * t248 + t570 * qJDD(1) + (-t189 * t376 - t272 * t575) * pkin(4) + (t344 + t625) * qJD(1) + t395;
t709 = t13 - g(2);
t708 = t683 * t272 - t273 * t684 + t713;
t707 = t272 * t681 - t273 * t682 + t715;
t706 = t376 * t723 - t377 * t721;
t705 = t376 * t722 + t377 * t720;
t704 = t374 * t716 + t375 * t726;
t703 = t374 * t726 - t375 * t716;
t680 = t376 * t730 + t377 * t732;
t679 = t376 * t729 + t377 * t731;
t699 = -t376 * t720 + t377 * t722 + t712;
t698 = qJD(1) * t733 + t376 * t721 + t377 * t723;
t697 = (t374 * t727 + t732) * t273 + (-t375 * t727 - t731) * t272 + (-t282 + t747) * qJD(1);
t695 = (t760 * t584 - t730 + t755) * t273 + (-t282 * t375 + t738) * t272 + t744 * qJD(1);
t694 = qJD(1) * t724 - t382 * t718 + t712;
t693 = -t719 * t382 + (-t374 * t728 - t377 * t729 + t717 + t745) * qJD(1);
t688 = -t694 * t374 + t698 * t375;
t687 = t693 * t374 + t699 * t375;
t686 = t698 * t374 + t694 * t375;
t685 = t699 * t374 - t693 * t375;
t674 = -t695 * t376 + t697 * t377;
t673 = qJD(1) * t728 - t272 * t719 + t273 * t718;
t467 = t666 + t369;
t291 = rSges(5,1) * t376 + rSges(5,2) * t377;
t234 = t291 * t374;
t236 = t291 * t375;
t665 = t377 * rSges(5,1) - rSges(5,2) * t376;
t667 = -rSges(5,2) * t585 - t375 * rSges(5,3);
t204 = t520 + t667;
t356 = t374 * rSges(5,3);
t206 = rSges(5,1) * t579 - rSges(5,2) * t580 + t356;
t531 = qJD(3) * t375;
t532 = qJD(3) * t374;
t508 = -t185 * t532 + t186 * t531 + qJD(2);
t67 = t204 * t272 + t206 * t273 + t508;
t417 = -t273 * t291 - t475;
t435 = -t204 + t471;
t72 = qJD(1) * t435 + t417;
t499 = t277 + t380;
t470 = t186 + t499;
t73 = -t272 * t291 - t317 + (t206 + t470) * qJD(1);
t670 = -(qJD(1) * t234 - t273 * t665) * t72 - t67 * (-t272 * t234 - t236 * t273) - t73 * (-qJD(1) * t236 - t272 * t665);
t638 = pkin(4) * t376;
t498 = -t290 - t638;
t61 = -t317 - t345 + t498 * t272 + (t470 + t570) * qJD(1);
t604 = qJD(1) * t61;
t669 = t12 + t604;
t271 = qJD(1) * t276;
t668 = qJD(1) * t185 - t271;
t357 = t374 * rSges(4,3);
t577 = t375 * t387;
t578 = t375 * t385;
t218 = rSges(4,1) * t577 - rSges(4,2) * t578 + t357;
t154 = t218 + t499;
t275 = t375 * rSges(3,1) - rSges(3,2) * t374;
t261 = t275 + t380;
t378 = Icges(4,4) * t387;
t452 = -Icges(4,2) * t385 + t378;
t323 = Icges(4,1) * t385 + t378;
t664 = g(1) * t375 + g(2) * t374;
t581 = t374 * t387;
t582 = t374 * t385;
t608 = Icges(4,3) * t375;
t209 = Icges(4,5) * t581 - Icges(4,6) * t582 - t608;
t332 = Icges(4,4) * t582;
t617 = Icges(4,5) * t375;
t213 = Icges(4,1) * t581 - t332 - t617;
t611 = Icges(4,6) * t375;
t211 = Icges(4,4) * t581 - Icges(4,2) * t582 - t611;
t595 = t211 * t385;
t443 = -t213 * t387 + t595;
t82 = -t209 * t375 - t374 * t443;
t320 = Icges(4,5) * t387 - Icges(4,6) * t385;
t319 = Icges(4,5) * t385 + Icges(4,6) * t387;
t419 = qJD(3) * t319;
t620 = Icges(4,4) * t385;
t324 = Icges(4,1) * t387 - t620;
t214 = Icges(4,5) * t374 + t324 * t375;
t212 = Icges(4,6) * t374 + t375 * t452;
t594 = t212 * t385;
t442 = -t214 * t387 + t594;
t658 = -t375 * t419 + (-t320 * t374 + t442 + t608) * qJD(1);
t210 = Icges(4,3) * t374 + t320 * t375;
t537 = qJD(1) * t210;
t657 = qJD(1) * t443 - t374 * t419 + t537;
t321 = Icges(4,2) * t387 + t620;
t437 = t321 * t385 - t323 * t387;
t654 = qJD(1) * t437 + t320 * qJD(3);
t554 = -Icges(4,2) * t581 + t213 - t332;
t556 = t323 * t374 + t211;
t653 = -t385 * t554 - t387 * t556;
t650 = t189 / 0.2e1;
t649 = t190 / 0.2e1;
t648 = t267 / 0.2e1;
t647 = t268 / 0.2e1;
t646 = -t272 / 0.2e1;
t645 = t272 / 0.2e1;
t644 = -t273 / 0.2e1;
t643 = t273 / 0.2e1;
t642 = t374 / 0.2e1;
t641 = -t375 / 0.2e1;
t635 = -qJD(1) / 0.2e1;
t634 = qJD(1) / 0.2e1;
t632 = rSges(4,1) * t387;
t628 = t375 * t73;
t626 = qJDD(1) / 0.2e1;
t43 = t272 * t571 + t273 * t570 + t508;
t605 = qJD(1) * t43;
t541 = rSges(4,2) * t582 + t375 * rSges(4,3);
t217 = rSges(4,1) * t581 - t541;
t469 = -t217 + t500;
t325 = rSges(4,1) * t385 + rSges(4,2) * t387;
t510 = t325 * t531;
t105 = qJD(1) * t469 - t510;
t603 = t105 * t374;
t602 = t105 * t375;
t106 = qJD(1) * t154 - t325 * t532;
t258 = t325 * t375;
t601 = t106 * t258;
t593 = t272 * t377;
t587 = t319 * t374;
t586 = t319 * t375;
t564 = -t374 * t185 + t375 * t186;
t563 = t374 * t204 + t375 * t206;
t562 = -t374 * t209 - t213 * t577;
t561 = t374 * t210 + t214 * t577;
t555 = -t323 * t375 - t212;
t553 = -t321 * t375 + t214;
t299 = pkin(4) * t514;
t551 = t290 * t535 + t299;
t547 = rSges(5,2) * t514 + rSges(5,3) * t534;
t533 = qJD(1) * t385;
t545 = rSges(4,2) * t374 * t533 + rSges(4,3) * t534;
t543 = -t321 + t324;
t542 = t323 + t452;
t536 = qJD(1) * t320;
t529 = qJD(3) * t387;
t143 = -t374 * t437 - t586;
t526 = t143 * qJD(1);
t120 = rSges(5,1) * t415 - rSges(5,2) * t516 + t547;
t428 = t291 * t382;
t122 = -t374 * t428 + (t375 * t665 + t356) * qJD(1);
t522 = t375 * t120 + t374 * t122 + t204 * t534;
t518 = pkin(3) * t529;
t515 = t375 * t149 + t374 * t150 - t185 * t534;
t512 = t375 * t533;
t507 = -pkin(2) - t632;
t506 = t535 / 0.2e1;
t505 = t534 / 0.2e1;
t504 = -t532 / 0.2e1;
t503 = t532 / 0.2e1;
t502 = -t531 / 0.2e1;
t501 = t531 / 0.2e1;
t416 = -t291 - t639;
t497 = t385 * (-t374 ^ 2 - t375 ^ 2);
t235 = t290 * t375;
t488 = -t272 * t233 - t235 * t273;
t182 = t214 * t581;
t487 = t210 * t375 - t182;
t484 = -t209 + t594;
t483 = -qJD(1) * t235 - t272 * t666;
t474 = t374 * t571 + t375 * t570;
t472 = -pkin(4) * t575 - t248;
t249 = t665 * t382;
t468 = -t249 - t518;
t315 = -t638 - t639;
t464 = qJD(1) * t233 - t273 * t467;
t328 = rSges(2,1) * t388 - rSges(2,2) * t386;
t326 = rSges(2,1) * t386 + rSges(2,2) * t388;
t327 = -rSges(4,2) * t385 + t632;
t126 = t212 * t387 + t214 * t385;
t420 = qJD(3) * t321;
t139 = -t375 * t420 + (-t374 * t452 + t611) * qJD(1);
t421 = qJD(3) * t323;
t141 = -t375 * t421 + (-t324 * t374 + t617) * qJD(1);
t397 = -qJD(3) * t126 - t139 * t385 + t141 * t387 + t537;
t125 = t211 * t387 + t213 * t385;
t140 = qJD(1) * t212 - t374 * t420;
t142 = qJD(1) * t214 - t374 * t421;
t398 = qJD(1) * t209 - qJD(3) * t125 - t140 * t385 + t142 * t387;
t460 = -(t374 * t657 + t398 * t375) * t375 + (t374 * t658 + t397 * t375) * t374;
t459 = -(t398 * t374 - t375 * t657) * t375 + (t397 * t374 - t375 * t658) * t374;
t458 = -t374 * t73 - t375 * t72;
t83 = -t212 * t582 - t487;
t457 = t374 * t83 - t375 * t82;
t84 = -t211 * t578 - t562;
t85 = -t212 * t578 + t561;
t456 = t374 * t85 - t375 * t84;
t449 = -t106 * t374 - t602;
t145 = -rSges(4,2) * t375 * t529 + (-t387 * t535 - t509) * rSges(4,1) + t545;
t257 = t325 * t374;
t146 = -qJD(3) * t257 + (t327 * t375 + t357) * qJD(1);
t448 = t145 * t375 + t146 * t374;
t441 = t217 * t374 + t218 * t375;
t438 = t321 * t387 + t323 * t385;
t436 = -t521 - t367;
t434 = t374 * t624 + t375 * t625 + t534 * t571;
t433 = -t290 + t315;
t418 = t149 * t531 + t150 * t532 - t267 * t185 - t186 * t268 + qJDD(2);
t414 = t472 - t518;
t411 = -t385 * t553 + t387 * t555;
t409 = (-t385 * t542 + t387 * t543) * qJD(1);
t408 = t273 * t498 + t344 - t475;
t295 = t452 * qJD(3);
t296 = t324 * qJD(3);
t396 = qJD(1) * t319 - qJD(3) * t438 - t295 * t385 + t296 * t387;
t394 = (t374 * t681 - t375 * t682) * t650 + (t374 * t683 - t375 * t684) * t649 + (t374 * t673 + t375 * t674) * t646 + (t688 * t375 + t687 * t374 + (t374 * t682 + t375 * t681) * qJD(1)) * t645 + (t686 * t375 + t685 * t374 + (t374 * t684 + t375 * t683) * qJD(1)) * t644 + (t374 * t674 - t375 * t673) * t643 + (t704 * qJD(1) + t701 * qJDD(1) + t681 * t189 + t682 * t190 + t687 * t272 + t688 * t273) * t642 + (t703 * qJD(1) + t702 * qJDD(1) + t683 * t189 + t684 * t190 + t685 * t272 + t686 * t273) * t641 + (t697 * t376 + t695 * t377) * t635 + (t706 * t375 + t705 * t374 + (t374 * t680 + t375 * t679) * qJD(1)) * t634 + (t374 * t679 - t375 * t680) * t626 + t708 * t506 + t707 * t505;
t300 = t327 * qJD(3);
t265 = t375 * t315;
t264 = t374 * t315;
t216 = pkin(3) * t578 + t265;
t215 = pkin(3) * t582 + t264;
t144 = -t375 * t437 + t587;
t131 = t144 * qJD(1);
t100 = qJD(3) * t441 + qJD(2);
t71 = t396 * t374 - t375 * t654;
t70 = t374 * t654 + t396 * t375;
t66 = qJD(1) * t145 + qJDD(1) * t218 - t267 * t325 - t300 * t532 + t430;
t65 = -t523 - t300 * t531 + t268 * t325 + (-t146 - t266) * qJD(1) + t469 * qJDD(1);
t64 = -qJD(3) * t442 + t139 * t387 + t141 * t385;
t63 = -t443 * qJD(3) + t140 * t387 + t142 * t385;
t62 = qJD(3) * t448 + t217 * t267 - t218 * t268 + qJDD(2);
t60 = qJD(1) * t429 + t408;
t49 = qJD(3) * t456 + t131;
t48 = qJD(3) * t457 + t526;
t38 = qJD(1) * t120 + qJDD(1) * t206 - t189 * t291 - t249 * t272 + t395;
t37 = t190 * t291 - t249 * t273 + (-t122 + t572) * qJD(1) + t435 * qJDD(1) + t410;
t18 = t120 * t273 + t122 * t272 + t189 * t204 - t190 * t206 + t418;
t9 = t189 * t571 - t190 * t570 + t272 * t624 + t273 * t625 + t418;
t1 = [(t131 + ((t83 - t182 + (t210 + t595) * t375 + t562) * t375 + t561 * t374) * qJD(3)) * t501 - m(2) * (-g(1) * t326 + g(2) * t328) + (t144 + t126) * t648 + (t143 + t125) * t647 + (((t700 + t746) * t375 + t683 + t734 + t736) * t273 + (t684 - t711) * t272 + t715) * t643 + (-t526 + ((t375 * t484 - t561 + t85) * t375 + (t374 * t484 + t487 + t84) * t374) * qJD(3) + t48) * t504 + (t64 + t70) * t503 + ((-t274 * t391 - g(2) + t473) * t261 + (-g(1) - t523 + (-0.2e1 * t275 - t380 + t261) * t391) * t260) * m(3) + (t63 + t71 + t49) * t502 + (-qJD(3) * t437 + t295 * t387 + t296 * t385 + t740 * t376 + t739 * t377) * qJD(1) + (-(-t60 + t408 + t668) * t61 - (-t571 - t640) * t604 + t60 * t345 + t61 * (-rSges(6,1) * t517 + t690) + t60 * (t290 * t382 - t269) * t374 + ((-t386 * t61 - t388 * t60) * pkin(1) + (t60 * (t436 + t629) - t61 * t381) * t375 + (t61 * t436 + t60 * t714) * t374) * qJD(1) + t709 * (t380 - t583 + t691) + t710 * (-t640 + t689)) * m(6) + ((-t428 - t519) * t628 + (t544 + (rSges(5,1) * t576 + rSges(5,2) * t575) * t374 + (-t356 - t380 + (-t370 - t665) * t375) * qJD(1)) * t72 + (t547 + t72 - t417 - t668 + (t204 + t640 + t725) * qJD(1)) * t73 + (t38 - g(2)) * (t206 + t380 + t478) + (t37 - g(1)) * (-t667 + t725)) * m(5) + (t106 * (t342 + t545) + (t325 * t603 - t601) * qJD(3) + ((-t105 * t388 - t106 * t386) * pkin(1) + (-pkin(2) - t327) * t602 + (t105 * (-rSges(4,3) - pkin(6)) + t106 * t507) * t374) * qJD(1) - (-t510 - t105 - t271 + (-t217 - t640) * qJD(1)) * t106 + (t66 - g(2)) * t154 + (t65 - g(1)) * (t507 * t374 + t363 + t541 - t640)) * m(4) + (t701 + t679) * t650 + (t702 + t680) * t649 + ((t681 + t711) * t273 + ((t700 + t724) * t375 + t696 * t374 + t682 - t742) * t272 + t708 - t713) * t646 + (t704 + t705) * t645 + (m(3) * (t260 ^ 2 + t275 * t261) + t438 + m(2) * (t326 ^ 2 + t328 ^ 2) + Icges(3,3) + Icges(2,3) + t751 * t377 + t727 * t376) * qJDD(1) + (t703 - t706 + t707) * t644; m(3) * qJDD(2) + m(4) * t62 + m(5) * t18 + m(6) * t9 + (-m(3) - m(4) - m(5) - m(6)) * g(3); ((t385 * t543 + t387 * t542) * qJD(1) + ((t374 * t553 - t375 * t554) * t387 + (t374 * t555 + t375 * t556) * t385) * qJD(3)) * t635 + t456 * t648 + t457 * t647 + (-t125 * t375 + t126 * t374) * t626 + (qJD(1) * t70 + qJD(3) * t460 + qJDD(1) * t144 + t267 * t85 + t268 * t84) * t642 + (t374 * t64 - t375 * t63 + (t125 * t374 + t126 * t375) * qJD(1)) * t634 + (qJD(1) * t71 + qJD(3) * t459 + qJDD(1) * t143 + t267 * t83 + t268 * t82) * t641 + ((t84 * t374 + t85 * t375) * qJD(1) + t460) * t503 + t394 + ((t82 * t374 + t83 * t375) * qJD(1) + t459) * t502 + t48 * t506 + ((-t532 * t586 + t536) * t374 + (t409 + (-t653 * t375 + (t587 + t411) * t374) * qJD(3)) * t375) * t504 + ((-t531 * t587 - t536) * t375 + (t409 + (t411 * t374 + (t586 - t653) * t375) * qJD(3)) * t374) * t501 + t49 * t505 + (t60 * t551 + t9 * (t474 + t564) + t43 * (t434 + t515) + (t414 * t60 + t433 * t669) * t375 + (t13 * t433 + t61 * t414 + (-t186 - t570) * t605) * t374 - t60 * (-qJD(1) * t215 + t464) - t61 * (-pkin(4) * t593 + qJD(1) * t216 + t483) - t43 * (t215 * t272 + t216 * t273 + t488) - (-t61 * t512 + ((-t374 * t61 - t375 * t60) * t387 + t43 * t497) * qJD(3)) * pkin(3) - g(1) * (t265 - t235) - g(2) * (t264 - t233) - g(3) * (t379 + t467)) * m(6) + (-g(3) * (t665 + t379) - t664 * t416 + t18 * (t563 + t564) + t67 * (t515 + t522) + (t468 * t72 + (qJD(1) * t73 + t37) * t416) * t375 + (t38 * t416 + t73 * t468 + (t72 * t291 + t67 * (-t186 - t206)) * qJD(1)) * t374 - (-t73 * t512 + (t387 * t458 + t497 * t67) * qJD(3)) * pkin(3) + t670) * m(5) + (t62 * t441 + t100 * ((t217 * t375 - t218 * t374) * qJD(1) + t448) + t449 * t300 + (-t66 * t374 - t65 * t375 + (-t106 * t375 + t603) * qJD(1)) * t325 - (t105 * t257 - t601) * qJD(1) - (t100 * (-t257 * t374 - t258 * t375) + t449 * t327) * qJD(3) + g(1) * t258 + g(2) * t257 - g(3) * t327) * m(4); t394 + (-g(3) * t467 + t9 * t474 + t43 * t434 + (t472 * t61 - t570 * t605) * t374 - t61 * t483 - t43 * t488 - (-t61 * t593 + (-t61 * t534 + t43 * (-t272 * t374 - t273 * t375)) * t376) * pkin(4) + (t13 * t374 + t375 * t669) * t498 + (t375 * t472 - t299 - t464 + t551) * t60 - t664 * (-t627 + (-rSges(6,1) - pkin(4)) * t376)) * m(6) + (t18 * t563 + t67 * (-t206 * t535 + t522) + t458 * t249 + (-t37 * t375 - t38 * t374 + (t374 * t72 - t628) * qJD(1)) * t291 + g(1) * t236 + g(2) * t234 - g(3) * t665 + t670) * m(5); ((t272 * t43 - t709) * t375 + (-t273 * t43 + t710) * t374) * m(6);];
tau = t1;
