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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:46:55
% EndTime: 2020-01-03 11:47:41
% DurationCPUTime: 36.11s
% Computational Cost: add. (19387->737), mult. (16589->889), div. (0->0), fcn. (12988->8), ass. (0->390)
t768 = Icges(5,4) + Icges(6,4);
t769 = Icges(5,1) + Icges(6,1);
t767 = Icges(5,2) + Icges(6,2);
t763 = Icges(5,6) + Icges(6,6);
t764 = Icges(5,5) + Icges(6,5);
t372 = qJ(3) + qJ(4);
t363 = cos(t372);
t766 = t768 * t363;
t362 = sin(t372);
t765 = t768 * t362;
t762 = t769 * t363 - t765;
t757 = -t767 * t362 + t766;
t371 = qJ(1) + pkin(8);
t360 = sin(t371);
t761 = t763 * t360;
t760 = Icges(5,3) + Icges(6,3);
t361 = cos(t371);
t563 = t361 * t363;
t564 = t361 * t362;
t740 = t768 * t563 - t767 * t564 + t761;
t751 = t767 * t363 + t765;
t758 = t768 * t564;
t756 = t764 * t360;
t747 = t769 * t362 + t766;
t741 = t757 * t360 - t763 * t361;
t739 = t762 * t360 - t764 * t361;
t752 = t769 * t563 + t756 - t758;
t755 = t764 * t362 + t763 * t363;
t730 = -t763 * t362 + t764 * t363;
t754 = t747 + t757;
t753 = -t762 + t751;
t750 = -t751 * t362 + t747 * t363;
t749 = t740 * t362;
t748 = t760 * t360;
t746 = -t751 * t360 + t739;
t745 = t747 * t361 + t740;
t744 = t747 * t360 + t741;
t743 = t730 * t360 - t760 * t361;
t742 = t764 * t563 - t763 * t564 + t748;
t731 = t755 * t361;
t729 = t752 * t363 - t749;
t370 = qJD(3) + qJD(4);
t738 = t753 * t370;
t737 = t754 * t370;
t724 = t755 * t360;
t698 = t750 * t360 - t731;
t697 = t747 * t563 - t751 * t564 + t724;
t567 = t360 * t363;
t736 = t739 * t567;
t735 = t752 * t567;
t734 = t741 * t564;
t728 = t744 * t370 + (-t361 * t762 - t756) * qJD(1);
t727 = t739 * qJD(1) + t745 * t370;
t726 = t746 * t370 + (t757 * t361 + t761) * qJD(1);
t562 = t361 * t370;
t725 = -t741 * qJD(1) + t752 * t370 - t751 * t562;
t723 = t739 * t363;
t722 = t741 * t362;
t568 = t360 * t362;
t683 = -t743 * t361 - t741 * t568 + t736;
t682 = t742 * t361 + t740 * t568 - t735;
t681 = -t743 * t360 - t739 * t563 + t734;
t680 = t742 * t360 + t729 * t361;
t721 = t750 * qJD(1) - t730 * t370;
t720 = -t729 - t743;
t719 = t743 * qJD(1);
t718 = -t755 * qJD(1) + t737 * t362 + t738 * t363;
t717 = t698 * qJD(1);
t263 = t360 * t370;
t715 = -t746 * t562 + (-t767 * t563 + t752 - t758) * t263 + t754 * qJD(1);
t714 = -t724 * t370 + (t730 * t361 + t722 - t723 + t748) * qJD(1);
t713 = t729 * qJD(1) + t731 * t370 + t719;
t712 = t697 * qJD(1);
t711 = -qJD(1) * t742 + t725 * t362 + t727 * t363;
t710 = t726 * t362 + t728 * t363 - t719;
t666 = -t753 * qJD(1) - t745 * t263 + t744 * t562;
t709 = -t715 * t362 + t363 * t666;
t708 = (-t742 + t723) * t361 - t734;
t503 = qJD(1) * qJD(3);
t334 = t360 * t503;
t502 = -qJDD(3) - qJDD(4);
t512 = qJD(1) * t360;
t177 = qJD(4) * t512 + t361 * t502 + t334;
t616 = rSges(6,2) * t362;
t653 = t363 * rSges(6,1) - t616;
t238 = t653 * t370;
t615 = rSges(6,2) * t363;
t282 = rSges(6,1) * t362 + t615;
t375 = cos(qJ(3));
t367 = t375 * pkin(3);
t353 = t367 + pkin(2);
t511 = qJD(1) * t361;
t286 = t353 * t511;
t377 = -pkin(7) - pkin(6);
t373 = sin(qJ(3));
t508 = qJD(3) * t373;
t497 = pkin(3) * t508;
t523 = pkin(2) * t511 + pkin(6) * t512;
t144 = t286 + (-qJD(1) * t377 - t497) * t360 - t523;
t347 = t360 * pkin(2);
t268 = -pkin(6) * t361 + t347;
t339 = t361 * t377;
t528 = t360 * t353 + t339;
t174 = -t268 + t528;
t259 = -qJDD(3) * t361 + t334;
t376 = cos(qJ(1));
t368 = t376 * pkin(1);
t374 = sin(qJ(1));
t379 = qJD(1) ^ 2;
t586 = pkin(1) * qJDD(1);
t522 = t379 * t368 + t374 * t586;
t454 = qJD(1) * t523 + qJDD(1) * t268 + t522;
t500 = qJD(3) ^ 2 * t367;
t623 = pkin(3) * t373;
t388 = qJD(1) * t144 + qJDD(1) * t174 - t259 * t623 + t361 * t500 + t454;
t506 = qJD(5) * t361;
t352 = pkin(4) * t363;
t295 = t352 + t353;
t369 = -qJ(5) + t377;
t612 = rSges(6,3) * t361;
t694 = rSges(6,1) * t567 - rSges(6,2) * t568 + t360 * t295 + t361 * t369 - t612;
t554 = -t528 + t694;
t558 = t363 * t370;
t492 = t360 * t558;
t408 = -t362 * t511 - t492;
t521 = t369 - t377;
t559 = t362 * t370;
t260 = -pkin(4) * t559 - t497;
t488 = t363 * t511;
t493 = t360 * t559;
t692 = rSges(6,3) * t512 + (t488 - t493) * rSges(6,1) + t360 * t260 + t295 * t511;
t605 = rSges(6,2) * t408 - t506 - t286 + (-qJD(1) * t521 + t497) * t360 + t692;
t12 = -qJDD(5) * t360 - t177 * t282 + t238 * t562 + t554 * qJDD(1) + (-t177 * t362 + t558 * t562) * pkin(4) + (-t506 + t605) * qJD(1) + t388;
t707 = t12 - g(3);
t176 = -qJD(1) * t562 + t360 * t502;
t341 = pkin(6) * t511;
t257 = pkin(2) * t512 - t341;
t343 = qJD(5) * t360;
t258 = -qJDD(3) * t360 - t361 * t503;
t366 = t374 * pkin(1);
t453 = -t366 * t379 + t376 * t586;
t390 = t258 * t623 - t360 * t500 + t453;
t624 = pkin(2) * t361;
t269 = pkin(6) * t360 + t624;
t303 = t361 * t353;
t175 = t624 - t303 + (pkin(6) + t377) * t360;
t693 = rSges(6,1) * t563 - rSges(6,2) * t564 + t361 * t295;
t553 = t303 - t693 + (-rSges(6,3) + t521) * t360;
t490 = -t175 - t553;
t455 = t269 + t490;
t456 = t361 * t497;
t143 = t456 + t341 + (t339 + (-pkin(2) + t353) * t360) * qJD(1);
t667 = -t295 - t653;
t606 = t282 * t562 - t343 + (-t260 - t497) * t361 + (t521 * t361 - t612 + (-t353 - t667) * t360) * qJD(1);
t495 = -t143 - t606;
t13 = -qJDD(5) * t361 + t176 * t282 - t238 * t263 + (t176 * t362 - t263 * t558) * pkin(4) + t455 * qJDD(1) + (-t257 + t495 + t343) * qJD(1) + t390;
t706 = t13 - g(2);
t705 = -t682 * t263 - t683 * t562 + t717;
t704 = -t680 * t263 - t681 * t562 - t712;
t703 = t728 * t362 - t726 * t363;
t702 = -t727 * t362 + t725 * t363;
t701 = t721 * t360 + t718 * t361;
t700 = -t718 * t360 + t361 * t721;
t679 = t739 * t362 + t741 * t363;
t699 = -t362 * t752 - t740 * t363;
t444 = -t360 * rSges(3,1) - t361 * rSges(3,2);
t691 = t444 - t366;
t689 = t720 * t562;
t688 = (t742 + t722) * t360 - t736;
t687 = t714 * t360 - t361 * t710;
t686 = t713 * t360 + t361 * t711;
t685 = t710 * t360 + t361 * t714;
t684 = -t711 * t360 + t361 * t713;
t283 = rSges(5,1) * t362 + rSges(5,2) * t363;
t677 = t283 * t562;
t670 = pkin(4) * t558 + t238;
t665 = -t730 * qJD(1) + t731 * t263 - t724 * t562;
t447 = t653 + t352;
t364 = Icges(4,4) * t375;
t432 = -Icges(4,2) * t373 + t364;
t651 = Icges(4,1) * t373 + t364;
t524 = t651 + t432;
t601 = Icges(4,4) * t373;
t319 = Icges(4,2) * t375 + t601;
t322 = Icges(4,1) * t375 - t601;
t525 = t319 - t322;
t661 = (t373 * t524 + t375 * t525) * qJD(1);
t262 = qJD(1) * t269;
t656 = -qJD(1) * t175 + t262;
t617 = rSges(5,2) * t362;
t652 = t363 * rSges(5,1) - t617;
t613 = rSges(5,3) * t361;
t191 = rSges(5,1) * t567 - rSges(5,2) * t568 - t613;
t478 = t268 + t366;
t451 = t174 + t478;
t71 = t456 + t677 + (t191 + t451) * qJD(1);
t357 = qJD(1) * t368;
t486 = t360 * t508;
t417 = -pkin(3) * t486 + t357;
t402 = -t263 * t283 + t417;
t450 = rSges(5,1) * t563 - rSges(5,2) * t564;
t193 = rSges(5,3) * t360 + t450;
t551 = -t175 + t193;
t489 = t269 + t551;
t72 = qJD(1) * t489 + t402;
t650 = t360 * t71 + t361 * t72;
t224 = t283 * t360;
t226 = rSges(5,1) * t564 + rSges(5,2) * t563;
t509 = qJD(3) * t361;
t510 = qJD(3) * t360;
t414 = t174 * t510 - t175 * t509 + qJD(2);
t66 = t191 * t263 + t193 * t562 + t414;
t649 = -t71 * (-qJD(1) * t224 + t562 * t652) - t66 * (-t224 * t263 - t226 * t562) - t72 * (-qJD(1) * t226 - t263 * t652);
t198 = -Icges(4,6) * t361 + t360 * t432;
t200 = -Icges(4,5) * t361 + t322 * t360;
t121 = t198 * t375 + t200 * t373;
t242 = t319 * t360;
t592 = Icges(4,6) * t360;
t134 = -qJD(3) * t242 + (t361 * t432 + t592) * qJD(1);
t244 = t651 * t360;
t598 = Icges(4,5) * t360;
t136 = -qJD(3) * t244 + (t322 * t361 + t598) * qJD(1);
t318 = Icges(4,5) * t375 - Icges(4,6) * t373;
t196 = -Icges(4,3) * t361 + t318 * t360;
t516 = qJD(1) * t196;
t648 = qJD(3) * t121 + t134 * t373 - t136 * t375 - t516;
t288 = t432 * qJD(3);
t289 = t322 * qJD(3);
t317 = Icges(4,5) * t373 + Icges(4,6) * t375;
t419 = t319 * t375 + t373 * t651;
t647 = -qJD(1) * t317 + qJD(3) * t419 + t288 * t373 - t289 * t375;
t133 = qJD(1) * t198 + t319 * t509;
t245 = t651 * t361;
t135 = qJD(1) * t200 + qJD(3) * t245;
t560 = t361 * t375;
t561 = t361 * t373;
t589 = Icges(4,3) * t360;
t197 = Icges(4,5) * t560 - Icges(4,6) * t561 + t589;
t199 = Icges(4,4) * t560 - Icges(4,2) * t561 + t592;
t329 = Icges(4,4) * t561;
t201 = Icges(4,1) * t560 - t329 + t598;
t424 = t199 * t375 + t201 * t373;
t646 = -qJD(1) * t197 + qJD(3) * t424 - t133 * t373 + t135 * t375;
t635 = t176 / 0.2e1;
t634 = t177 / 0.2e1;
t633 = t258 / 0.2e1;
t632 = t259 / 0.2e1;
t631 = -t263 / 0.2e1;
t630 = t263 / 0.2e1;
t629 = t562 / 0.2e1;
t628 = -t562 / 0.2e1;
t627 = -t360 / 0.2e1;
t626 = -t361 / 0.2e1;
t625 = rSges(4,3) + pkin(6);
t622 = pkin(4) * t362;
t621 = g(2) * t360;
t620 = -qJD(1) / 0.2e1;
t619 = qJD(1) / 0.2e1;
t618 = rSges(4,1) * t375;
t614 = rSges(4,3) * t361;
t609 = qJDD(1) / 0.2e1;
t608 = rSges(5,3) - t377;
t607 = rSges(6,3) - t369;
t565 = t360 * t375;
t566 = t360 * t373;
t448 = rSges(4,1) * t565 - rSges(4,2) * t566;
t204 = t448 - t614;
t323 = rSges(4,1) * t373 + rSges(4,2) * t375;
t101 = t323 * t509 + (t204 + t478) * qJD(1);
t246 = t323 * t360;
t585 = t101 * t246;
t445 = -t323 * t510 + t357;
t331 = rSges(4,2) * t561;
t498 = rSges(4,1) * t560;
t205 = rSges(4,3) * t360 - t331 + t498;
t537 = t205 + t269;
t102 = qJD(1) * t537 + t445;
t584 = t102 * t360;
t583 = t102 * t361;
t576 = t198 * t373;
t575 = t199 * t373;
t574 = t200 * t375;
t573 = t562 * t363;
t312 = -t622 - t623;
t570 = t312 * t361;
t569 = t317 * t360;
t241 = t317 * t361;
t116 = t677 + (t360 * t652 - t613) * qJD(1);
t556 = -t116 - t143;
t555 = t554 * t360;
t542 = t447 * t263;
t541 = t198 + t244;
t540 = t199 + t245;
t539 = t200 - t242;
t538 = -Icges(4,2) * t560 + t201 - t329;
t536 = t670 * t361;
t529 = rSges(5,1) * t488 + rSges(5,3) * t512;
t225 = rSges(6,1) * t564 + rSges(6,2) * t563;
t313 = pkin(4) * t564;
t527 = t361 * t282 + t313;
t526 = rSges(4,3) * t512 + qJD(1) * t498;
t513 = qJD(1) * t318;
t507 = qJD(3) * t375;
t138 = -t319 * t561 + t560 * t651 + t569;
t504 = t138 * qJD(1);
t501 = pkin(3) * t566;
t333 = pkin(3) * t561;
t496 = pkin(3) * t507;
t118 = -rSges(5,1) * t493 + rSges(5,2) * t408 + t529;
t494 = t360 * t118 + t191 * t511 - t193 * t512;
t491 = t360 * t144 + t174 * t511 + t175 * t512;
t485 = t360 * t507;
t484 = t512 / 0.2e1;
t483 = -t511 / 0.2e1;
t482 = -t510 / 0.2e1;
t481 = t510 / 0.2e1;
t480 = -t509 / 0.2e1;
t479 = t509 / 0.2e1;
t477 = -t283 - t623;
t476 = -t282 - t622;
t267 = rSges(3,1) * t361 - t360 * rSges(3,2);
t223 = t282 * t360;
t463 = -qJD(1) * t223 + t562 * t653;
t457 = qJD(3) * (-t360 ^ 2 - t361 ^ 2);
t326 = rSges(2,1) * t376 - t374 * rSges(2,2);
t324 = rSges(2,1) * t374 + rSges(2,2) * t376;
t325 = -rSges(4,2) * t373 + t618;
t425 = t574 - t576;
t394 = qJD(3) * t569 + (-t318 * t361 + t425 - t589) * qJD(1);
t423 = -t201 * t375 + t575;
t395 = qJD(1) * t423 - qJD(3) * t241 - t516;
t441 = -(t394 * t360 + t361 * t648) * t361 - (t395 * t360 - t361 * t646) * t360;
t440 = -(-t360 * t648 + t394 * t361) * t361 - (t360 * t646 + t395 * t361) * t360;
t171 = t200 * t565;
t81 = -t196 * t361 - t198 * t566 + t171;
t172 = t201 * t565;
t82 = t197 * t361 + t199 * t566 - t172;
t439 = -t360 * t82 - t361 * t81;
t173 = t198 * t561;
t83 = -t196 * t360 - t200 * t560 + t173;
t84 = t197 * t360 - t423 * t361;
t438 = -t360 * t84 - t361 * t83;
t437 = qJD(1) * t477;
t436 = qJD(1) * t476;
t429 = t101 * t361 - t584;
t247 = t323 * t361;
t139 = qJD(3) * t247 + (t325 * t360 - t614) * qJD(1);
t140 = -rSges(4,1) * t486 + (-t373 * t511 - t485) * rSges(4,2) + t526;
t428 = -t139 * t361 + t140 * t360;
t422 = t204 * t360 + t205 * t361;
t418 = -t319 * t373 + t375 * t651;
t416 = -t282 + t312;
t415 = t360 * t605 + t511 * t554 + t512 * t553;
t410 = t650 * qJD(1);
t409 = qJD(1) * t416;
t401 = t373 * t539 + t375 * t541;
t400 = t373 * t538 + t375 * t540;
t391 = t418 * qJD(1) - t318 * qJD(3);
t389 = -t143 * t509 + t144 * t510 - t174 * t258 + t259 * t175 + qJDD(2);
t381 = t263 * t476 + t417 - t506;
t380 = (-t360 * t680 - t361 * t681) * t635 + (-t360 * t682 - t361 * t683) * t634 + (t687 * t361 + t686 * t360 + (t360 * t681 - t361 * t680) * qJD(1)) * t631 + (t665 * t360 - t361 * t709) * t630 + (t709 * t360 + t665 * t361) * t629 + (t685 * t361 + t684 * t360 + (t360 * t683 - t361 * t682) * qJD(1)) * t628 + (t701 * qJD(1) - t697 * qJDD(1) + t680 * t176 + t681 * t177 + t686 * t263 + t687 * t562) * t627 + (t700 * qJD(1) + t698 * qJDD(1) + t682 * t176 + t683 * t177 + t684 * t263 + t685 * t562) * t626 + (t666 * t362 + t363 * t715) * t620 + (t703 * t361 + t702 * t360 + (t360 * t679 - t361 * t699) * qJD(1)) * t619 + (-t360 * t699 - t361 * t679) * t609 + t705 * t484 + t704 * t483;
t292 = t325 * qJD(3);
t256 = t360 * t312;
t252 = t267 + t368;
t239 = t652 * t370;
t203 = -t333 - t570;
t202 = t256 + t501;
t168 = t360 * t191;
t161 = t562 * t225;
t160 = t360 * t174;
t137 = t360 * t418 - t241;
t125 = t137 * qJD(1);
t96 = qJD(3) * t422 + qJD(2);
t70 = -t360 * t647 + t391 * t361;
t69 = t391 * t360 + t361 * t647;
t65 = -t292 * t510 + t258 * t323 + t537 * qJDD(1) + (-t139 - t257) * qJD(1) + t453;
t64 = qJD(1) * t140 + qJDD(1) * t204 - t259 * t323 + t292 * t509 + t454;
t63 = qJD(3) * t423 + t133 * t375 + t135 * t373;
t62 = qJD(3) * t425 + t134 * t375 + t136 * t373;
t61 = qJD(3) * t428 - t204 * t258 - t205 * t259 + qJDD(2);
t60 = qJD(1) * t455 + t381;
t59 = t456 - t343 - t476 * t562 + (t451 + t554) * qJD(1);
t49 = qJD(3) * t438 - t504;
t48 = qJD(3) * t439 + t125;
t43 = t263 * t554 - t553 * t562 + t414;
t38 = t176 * t283 - t239 * t263 + t489 * qJDD(1) + (-t257 + t556) * qJD(1) + t390;
t37 = qJD(1) * t118 + qJDD(1) * t191 - t177 * t283 + t239 * t562 + t388;
t18 = -t116 * t562 + t118 * t263 - t176 * t191 - t177 * t193 + t389;
t9 = -t176 * t554 + t177 * t553 + t263 * t605 - t562 * t606 + t389;
t1 = [-t424 * t633 - m(2) * (g(2) * t326 + g(3) * t324) - t258 * t138 / 0.2e1 + (t125 + ((t83 + t172 - t173 + (t196 - t575) * t360) * t360 + (-t171 - t84 + (t196 - t423) * t361 + (t574 + t576) * t360) * t361) * qJD(3)) * t481 - t697 * t176 / 0.2e1 + t699 * t635 + (t137 + t121) * t632 + ((t688 - t680) * t562 - t689 * t361 + ((t743 - t749) * t360 + t681 + t708 + t735) * t263 + t717) * t630 + (t62 + t70) * t480 + (t49 + t504 + ((t173 - t82 + (t197 - t574) * t361) * t361 + (-t171 + t81 + (t197 + t576) * t360) * t360) * qJD(3)) * t479 + (qJD(3) * t418 + t288 * t375 + t289 * t373 - t738 * t362 + t737 * t363) * qJD(1) + ((t379 * t444 - g(2) + t453) * t252 + (g(3) - t522 + (-0.2e1 * rSges(3,1) * t511 + 0.2e1 * rSges(3,2) * t512 + qJD(1) * t267) * qJD(1)) * t691) * m(3) + (t63 + t69 + t48) * t482 + (-(-qJD(1) * t553 + t381 - t60 + t656) * t59 + t60 * t343 + t59 * (-rSges(6,2) * t492 + t357 + t692) + (t60 * (-rSges(6,1) * t559 - rSges(6,2) * t558 + t260) - t59 * qJD(5)) * t361 + (-t60 * t366 + (-t59 * t616 + t60 * t607) * t361 + (-t59 * t369 + t60 * t667) * t360) * qJD(1) + t706 * (t360 * t607 + t368 + t693) + t707 * (t366 + t694)) * m(6) + (-(qJD(1) * t193 + t402 + t656 - t72) * t71 + t71 * (t286 + t357 + t529) + (-t72 * t366 + (t608 * t72 - t617 * t71) * t361 + (t72 * (-t353 - t652) - t71 * t377) * t360) * qJD(1) + t650 * (-t283 * t370 - t497) + (t38 - g(2)) * (t360 * t608 + t303 + t368 + t450) + (t37 - g(3)) * (t191 + t366 + t528)) * m(5) + (t102 * t341 + (-t323 * t583 - t585) * qJD(3) + (t102 * (t614 - t366) + (-pkin(2) - t325) * t584) * qJD(1) + (t65 - g(2)) * (-t331 + t368 + (pkin(2) + t618) * t361 + t625 * t360) + (t64 - g(3)) * (-t361 * t625 + t347 + t366 + t448) + (t102 - t262 + t357 - t445 + t523 + t526 + (-t205 - t331) * qJD(1)) * t101) * m(4) + (t698 + t679) * t634 + (t689 * t360 + (-t682 - t708) * t562 + (-t361 * t720 + t683 + t688) * t263 + t704 + t712) * t629 + (t700 - t703) * t628 + (m(3) * (t252 * t267 + t444 * t691) + t419 + m(2) * (t324 ^ 2 + t326 ^ 2) + Icges(2,3) + Icges(3,3) + t751 * t363 + t747 * t362) * qJDD(1) + (t701 - t702 + t705) * t631; m(3) * qJDD(2) + m(4) * t61 + m(5) * t18 + m(6) * t9 + (-m(3) - m(4) - m(5) - m(6)) * g(1); ((t360 * t81 - t82 * t361) * qJD(1) + t440) * t480 + ((-t373 * t525 + t375 * t524) * qJD(1) + ((t360 * t538 - t361 * t539) * t375 + (-t360 * t540 + t361 * t541) * t373) * qJD(3)) * t620 + ((-t509 * t569 - t513) * t361 + (-t661 + (-t400 * t360 + (t241 + t401) * t361) * qJD(3)) * t360) * t479 + t438 * t633 + t439 * t632 + (-t121 * t361 + t360 * t424) * t609 + (-t360 * t63 - t361 * t62 + (t121 * t360 + t361 * t424) * qJD(1)) * t619 + ((t83 * t360 - t361 * t84) * qJD(1) + t441) * t482 + (qJD(1) * t70 + qJD(3) * t440 + qJDD(1) * t137 + t258 * t82 + t259 * t81) * t626 + t380 + (qJD(1) * t69 + qJD(3) * t441 - qJDD(1) * t138 + t258 * t84 + t259 * t83) * t627 + t48 * t484 + t49 * t483 + ((t241 * t510 - t513) * t360 + (t661 + (-t401 * t361 + (-t569 + t400) * t360) * qJD(3)) * t361) * t481 + (-g(1) * (t367 + t447) - g(2) * (t256 - t223) - g(3) * (t225 - t570) + t12 * (t333 + t527) + t13 * t416 * t360 + (t361 * t490 + t160 + t555) * t9 + (t409 * t361 + (-t670 - t496) * t360 + pkin(3) * t485 + t542 - (-t203 - t225 - t333) * qJD(1)) * t60 + (t536 + t409 * t360 - pkin(4) * t573 - t463 - (t202 - t501) * qJD(1)) * t59 + (t415 + t491 + t495 * t361 + t203 * t562 + t161 - (t202 - t223) * t263 - t457 * t623) * t43) * m(6) + (-g(1) * (t652 + t367) - g(3) * (t333 + t226) - t477 * t621 + t18 * (t160 + t168) + t66 * (t491 + t494) + t37 * t333 + (t18 * t551 + t71 * t239 + t37 * t283 + t437 * t72 + t556 * t66) * t361 + (t38 * t477 + t72 * (-t239 - t496) + t71 * t437) * t360 - (-t72 * t485 + (t457 * t66 - t410) * t373) * pkin(3) + t649) * m(5) + (-g(1) * t325 + g(2) * t246 - g(3) * t247 + t61 * t422 + t96 * ((t204 * t361 - t205 * t360) * qJD(1) + t428) + t429 * t292 + (-t65 * t360 + t64 * t361 + (-t101 * t360 - t583) * qJD(1)) * t323 - (-t102 * t247 - t585) * qJD(1) - (t96 * (-t246 * t360 - t247 * t361) + t429 * t325) * qJD(3)) * m(4); t380 + (-g(1) * t447 - g(3) * (t313 + t225) - (-t615 + (-rSges(6,1) - pkin(4)) * t362) * t621 + t9 * t555 + t43 * t415 + t12 * t527 + t59 * t536 + (-t43 * t606 + t436 * t60 - t553 * t9) * t361 + (t13 * t476 + t436 * t59 - t60 * t670) * t360 - t43 * (-t223 * t263 - t161) - t60 * (-qJD(1) * t225 - t542) - t59 * t463 - (t59 * t573 + (t43 * (-t263 * t360 - t361 * t562) + (-t360 * t59 - t361 * t60) * qJD(1)) * t362) * pkin(4)) * m(6) + (-g(1) * t652 + g(2) * t224 - g(3) * t226 + t18 * (t193 * t361 + t168) + t66 * (-t116 * t361 + t494) + (-t360 * t72 + t361 * t71) * t239 + (-t38 * t360 + t37 * t361 - t410) * t283 + t649) * m(5); ((t263 * t43 - t706) * t361 + (-t43 * t562 - t707) * t360) * m(6);];
tau = t1;
