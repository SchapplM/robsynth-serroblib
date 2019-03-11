% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:20
% EndTime: 2019-03-09 06:14:48
% DurationCPUTime: 15.21s
% Computational Cost: add. (32759->735), mult. (65917->960), div. (0->0), fcn. (74996->8), ass. (0->382)
t399 = sin(pkin(10));
t400 = cos(pkin(10));
t403 = sin(qJ(3));
t636 = cos(qJ(3));
t381 = t399 * t403 - t400 * t636;
t383 = t399 * t636 + t403 * t400;
t401 = sin(qJ(5));
t402 = sin(qJ(4));
t404 = cos(qJ(5));
t405 = cos(qJ(4));
t444 = t401 * t402 - t404 * t405;
t294 = t444 * t383;
t279 = t294 * qJ(6);
t497 = -pkin(2) * t400 - pkin(1);
t626 = pkin(8) * t383;
t311 = pkin(3) * t381 + t497 - t626;
t611 = pkin(7) + qJ(2);
t389 = t611 * t400;
t489 = t611 * t399;
t328 = t389 * t636 - t403 * t489;
t182 = t405 * t311 - t402 * t328;
t544 = t383 * t405;
t153 = -pkin(9) * t544 + t182;
t123 = t381 * pkin(4) + t153;
t183 = t402 * t311 + t328 * t405;
t545 = t383 * t402;
t154 = -pkin(9) * t545 + t183;
t147 = t401 * t154;
t78 = t404 * t123 - t147;
t56 = t279 + t78;
t48 = pkin(5) * t381 + t56;
t88 = t404 * t153 - t147;
t66 = t279 + t88;
t759 = t48 - t66;
t612 = Ifges(7,4) + Ifges(6,4);
t737 = t444 * t612;
t445 = t401 * t405 + t404 * t402;
t738 = t445 * t612;
t758 = t612 * t294;
t296 = t445 * t383;
t757 = t612 * t296;
t675 = -pkin(9) - pkin(8);
t390 = t675 * t402;
t391 = t675 * t405;
t336 = t390 * t401 - t404 * t391;
t291 = -qJ(6) * t444 + t336;
t713 = Ifges(7,6) + Ifges(6,6);
t714 = Ifges(7,5) + Ifges(6,5);
t479 = -t714 * t444 - t713 * t445;
t709 = t404 * t390 + t401 * t391;
t720 = -t445 * qJ(6) + t709;
t756 = -t336 * mrSges(6,1) - t291 * mrSges(7,1) - t709 * mrSges(6,2) - t720 * mrSges(7,2) + t479;
t218 = mrSges(7,1) * t381 + t294 * mrSges(7,3);
t219 = mrSges(6,1) * t381 + t294 * mrSges(6,3);
t651 = t336 / 0.2e1;
t754 = t291 / 0.2e1;
t755 = t218 * t754 + t219 * t651;
t643 = -t445 / 0.2e1;
t752 = mrSges(7,3) + mrSges(6,3);
t751 = Ifges(6,1) + Ifges(7,1);
t731 = Ifges(6,2) + Ifges(7,2);
t565 = qJ(6) * t296;
t149 = t404 * t154;
t79 = t123 * t401 + t149;
t57 = t79 - t565;
t87 = -t153 * t401 - t149;
t65 = t87 + t565;
t750 = t57 + t65;
t749 = t78 - t88;
t608 = mrSges(7,3) * t444;
t745 = t296 * t444 / 0.2e1 + t294 * t643;
t322 = mrSges(7,1) * t444 + mrSges(7,2) * t445;
t628 = pkin(4) * t405;
t396 = -pkin(3) - t628;
t342 = pkin(5) * t444 + t396;
t744 = m(7) * t342 + t322;
t743 = t79 + t87;
t607 = mrSges(7,3) * t445;
t574 = t444 * mrSges(6,3);
t663 = -t720 / 0.2e1;
t734 = t663 + t720 / 0.2e1;
t629 = pkin(4) * t404;
t295 = t445 * t381;
t297 = t444 * t381;
t729 = t731 * t295 + t612 * t297 + t713 * t383;
t728 = t612 * t295 + t751 * t297 + t714 * t383;
t395 = pkin(5) + t629;
t727 = t395 - t629;
t726 = -t731 * t444 + t738;
t725 = -t731 * t296 - t758;
t724 = -t751 * t444 - t738;
t723 = t751 * t445 - t737;
t722 = -t751 * t294 - t757;
t627 = pkin(8) * t381;
t631 = pkin(3) * t383;
t321 = t627 + t631;
t327 = t389 * t403 + t636 * t489;
t196 = t405 * t321 + t402 * t327;
t197 = t402 * t321 - t405 * t327;
t721 = -t196 * t402 + t197 * t405;
t684 = mrSges(7,3) / 0.2e1;
t685 = mrSges(6,3) / 0.2e1;
t516 = t684 + t685;
t428 = -t218 / 0.2e1 - t219 / 0.2e1 + t516 * t294;
t610 = t48 - t56;
t715 = m(7) * t610;
t719 = -t715 / 0.2e1 + t428;
t572 = t445 * mrSges(6,3);
t718 = t291 * t607 + t336 * t572;
t717 = -t731 * t445 + t723 - t737;
t245 = pkin(4) * t545 + t327;
t157 = pkin(5) * t296 + t245;
t716 = -t157 / 0.2e1;
t630 = pkin(4) * t401;
t707 = -t197 * mrSges(5,2) / 0.2e1 + t196 * mrSges(5,1) / 0.2e1;
t480 = t294 * t713 - t296 * t714;
t706 = -mrSges(5,1) * t405 + t402 * mrSges(5,2);
t510 = t88 / 0.2e1 - t78 / 0.2e1;
t698 = t405 ^ 2;
t699 = t402 ^ 2;
t705 = t699 + t698;
t704 = -mrSges(5,3) * t705 / 0.2e1;
t606 = Ifges(5,4) * t402;
t459 = Ifges(5,2) * t405 + t606;
t605 = Ifges(5,4) * t405;
t467 = Ifges(5,1) * t402 + t605;
t638 = -t405 / 0.2e1;
t639 = t402 / 0.2e1;
t703 = t459 * t639 + t467 * t638;
t280 = t296 * mrSges(7,2);
t487 = -t294 * mrSges(7,1) - t280;
t702 = -t396 * (-mrSges(6,1) * t294 - mrSges(6,2) * t296) / 0.2e1 - t342 * t487 / 0.2e1 - t479 * t381 / 0.4e1;
t584 = t296 * mrSges(6,3);
t215 = -mrSges(6,2) * t381 - t584;
t583 = t296 * mrSges(7,3);
t214 = -mrSges(7,2) * t381 - t583;
t647 = -t444 / 0.2e1;
t532 = t214 * t647 + t218 * t643;
t701 = t215 * t647 + t219 * t643 - t752 * t745 + t532;
t700 = 0.2e1 * m(7);
t697 = m(5) / 0.2e1;
t696 = -m(6) / 0.2e1;
t695 = m(6) / 0.2e1;
t694 = -m(7) / 0.2e1;
t693 = m(7) / 0.2e1;
t692 = pkin(5) / 0.2e1;
t691 = pkin(8) / 0.2e1;
t690 = m(7) * pkin(5);
t689 = mrSges(6,1) / 0.2e1;
t688 = -mrSges(6,2) / 0.2e1;
t687 = -mrSges(7,2) / 0.2e1;
t683 = t48 / 0.2e1;
t546 = t381 * t405;
t133 = t383 * pkin(4) + pkin(9) * t546 + t196;
t547 = t381 * t402;
t163 = pkin(9) * t547 + t197;
t82 = t404 * t133 - t163 * t401;
t51 = pkin(5) * t383 - qJ(6) * t297 + t82;
t682 = t51 / 0.2e1;
t681 = -t56 / 0.2e1;
t680 = -t57 / 0.2e1;
t679 = -t65 / 0.2e1;
t677 = -t79 / 0.2e1;
t674 = m(7) * t57;
t673 = pkin(3) * mrSges(5,1);
t672 = pkin(3) * mrSges(5,2);
t166 = mrSges(7,1) * t296 - mrSges(7,2) * t294;
t671 = -t166 / 0.2e1;
t670 = t166 / 0.2e1;
t669 = -t214 / 0.2e1;
t668 = -t215 / 0.2e1;
t216 = mrSges(7,1) * t383 - t297 * mrSges(7,3);
t667 = t216 / 0.2e1;
t576 = t383 * mrSges(6,1);
t217 = -t297 * mrSges(6,3) + t576;
t666 = t217 / 0.2e1;
t665 = t245 / 0.2e1;
t662 = -t291 / 0.2e1;
t657 = t295 / 0.2e1;
t656 = -t296 / 0.2e1;
t654 = t297 / 0.2e1;
t653 = -t322 / 0.2e1;
t652 = -t327 / 0.2e1;
t650 = -t709 / 0.2e1;
t648 = t383 / 0.2e1;
t645 = t445 / 0.2e1;
t641 = -t395 / 0.2e1;
t640 = -t402 / 0.2e1;
t637 = t405 / 0.2e1;
t635 = m(7) * t157;
t634 = m(7) * t445;
t633 = m(7) * t395;
t632 = m(7) * t401;
t625 = pkin(8) * t402;
t624 = pkin(8) * t405;
t623 = t445 * pkin(5);
t622 = t402 * pkin(4);
t621 = t56 * mrSges(7,2);
t620 = t57 * mrSges(7,1);
t619 = t65 * mrSges(7,1);
t618 = t66 * mrSges(7,2);
t617 = t78 * mrSges(6,2);
t616 = t79 * mrSges(6,1);
t615 = t87 * mrSges(6,1);
t614 = t88 * mrSges(6,2);
t613 = mrSges(6,2) + mrSges(7,2);
t596 = Ifges(5,5) * t381;
t595 = Ifges(5,5) * t405;
t594 = Ifges(5,6) * t381;
t593 = Ifges(5,6) * t402;
t592 = pkin(4) * qJD(4);
t591 = t157 * mrSges(7,1);
t590 = t245 * mrSges(6,2);
t585 = t295 * mrSges(7,1);
t582 = t297 * mrSges(7,2);
t246 = -pkin(4) * t547 + t328;
t158 = -t295 * pkin(5) + t246;
t164 = t582 - t585;
t165 = -t295 * mrSges(6,1) + t297 * mrSges(6,2);
t167 = mrSges(6,1) * t296 - mrSges(6,2) * t294;
t212 = -mrSges(7,2) * t383 + t295 * mrSges(7,3);
t575 = t383 * mrSges(6,2);
t213 = t295 * mrSges(6,3) - t575;
t460 = -Ifges(5,2) * t402 + t605;
t234 = Ifges(5,6) * t383 - t381 * t460;
t468 = Ifges(5,1) * t405 - t606;
t236 = Ifges(5,5) * t383 - t381 * t468;
t568 = t405 * mrSges(5,2);
t570 = t402 * mrSges(5,1);
t470 = t568 + t570;
t309 = t470 * t381;
t310 = t470 * t383;
t316 = -mrSges(5,2) * t383 + mrSges(5,3) * t547;
t317 = -mrSges(5,2) * t381 - mrSges(5,3) * t545;
t318 = t383 * mrSges(5,1) + mrSges(5,3) * t546;
t319 = t381 * mrSges(5,1) - mrSges(5,3) * t544;
t367 = t383 * mrSges(4,1);
t439 = -t595 / 0.2e1 + t593 / 0.2e1;
t433 = t439 * t381;
t514 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t515 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t436 = t468 * t383;
t237 = t436 + t596;
t534 = t405 * t237;
t435 = t460 * t383;
t235 = t435 + t594;
t537 = t402 * t235;
t83 = t401 * t133 + t404 * t163;
t61 = qJ(6) * t295 + t83;
t3 = (t236 * t637 + t234 * t640 + (-Ifges(4,4) - t439) * t383 - t514 * t296 - t515 * t294) * t383 - t728 * t294 / 0.2e1 + t729 * t656 + t722 * t654 + t725 * t657 + (-t497 * mrSges(4,2) + Ifges(4,4) * t381 - t534 / 0.2e1 + t537 / 0.2e1 + t433 + t714 * t297 + t713 * t295 + (Ifges(5,3) - Ifges(4,1) + Ifges(4,2) + Ifges(7,3) + Ifges(6,3)) * t383) * t381 - t327 * t309 + t328 * t310 + t197 * t317 + t182 * t318 + t196 * t319 + t183 * t316 + t245 * t165 + t246 * t167 + t51 * t218 + t82 * t219 + t57 * t212 + t79 * t213 + t61 * t214 + t83 * t215 + t48 * t216 + t78 * t217 + t157 * t164 + t158 * t166 + t497 * t367 + m(5) * (t182 * t196 + t183 * t197 + t327 * t328) + m(6) * (t245 * t246 + t78 * t82 + t79 * t83) + m(7) * (t157 * t158 + t48 * t51 + t57 * t61);
t581 = t3 * qJD(1);
t573 = t445 * mrSges(7,1);
t571 = t401 * t61;
t210 = pkin(4) * t544 - t294 * pkin(5);
t411 = -t157 * t280 - (t381 * t515 + t590 - t757) * t296 + t48 * t583 + t78 * t584 + t480 * t381 / 0.2e1;
t488 = t731 - t751;
t415 = -t245 * mrSges(6,1) + t79 * mrSges(6,3) + t57 * mrSges(7,3) - t296 * t488 + t381 * t514 - t591 - t758;
t446 = t182 * t402 - t183 * t405;
t451 = Ifges(5,5) * t402 + Ifges(5,6) * t405;
t6 = t182 * t317 - t183 * t319 + t65 * t218 + t87 * t219 + t210 * t166 + t66 * t214 + t88 * t215 + m(6) * (t78 * t87 + t79 * t88) + m(7) * (t157 * t210 + t48 * t65 + t57 * t66) + (-t381 * t451 / 0.2e1 - t327 * t706 + t235 * t638 + t237 * t640 + t703 * t383 + (m(6) * t245 + t167) * t628 + t446 * mrSges(5,3)) * t383 + t415 * t294 + t411;
t567 = t6 * qJD(1);
t7 = t56 * t214 + t78 * t215 - t79 * t219 + (-t218 - t715) * t57 + ((-t166 - t635) * pkin(5) + t415) * t294 + t411;
t566 = t7 * qJD(1);
t408 = -m(5) * (-t627 * t705 - t631) / 0.2e1 + (t295 * t709 + t297 * t336 + t383 * t396) * t696 + (t291 * t297 + t295 * t720 + t342 * t383) * t694 + t383 * t653 - t706 * t648;
t410 = (t196 * t405 + t402 * t197) * t697 + (-t444 * t82 + t445 * t83) * t695 + (-t444 * t51 + t445 * t61) * t693 + t316 * t639 + t318 * t637;
t490 = t667 + t666;
t492 = t213 / 0.2e1 + t212 / 0.2e1;
t15 = t367 - (t575 / 0.2e1 - t516 * t295 - t492) * t445 - (t576 / 0.2e1 - t516 * t297 + t490) * t444 + (-mrSges(4,2) + (t699 / 0.2e1 + t698 / 0.2e1) * mrSges(5,3)) * t381 + t408 + t410;
t563 = qJD(1) * t15;
t30 = -t296 * t214 + t294 * t218 + m(7) * (t294 * t48 - t296 * t57);
t562 = qJD(1) * t30;
t533 = t405 * t317;
t536 = t402 * t319;
t550 = t327 * t383;
t14 = (t214 + t215) * t297 + (t218 + t219) * t295 + (mrSges(4,3) * t381 - t533 + t536) * t381 + (mrSges(4,3) * t383 + t166 + t167 + t310) * t383 + m(6) * (t245 * t383 + t295 * t78 + t297 * t79) + m(5) * (t381 * t446 + t550) + m(4) * (-t328 * t381 + t550) + m(7) * (t157 * t383 + t295 * t48 + t297 * t57) + (m(3) * qJ(2) + mrSges(3,3)) * (t399 ^ 2 + t400 ^ 2);
t561 = t14 * qJD(1);
t558 = t245 * t402;
t557 = t291 * t294;
t555 = t291 * t395;
t554 = t294 * t322;
t552 = t294 * t401;
t549 = t336 * t294;
t548 = t342 * t294;
t543 = t444 * t294;
t542 = t445 * t296;
t541 = t445 * t401;
t539 = t395 * t294;
t102 = (t383 / 0.4e1 + t543 / 0.4e1 + t542 / 0.4e1) * t700;
t525 = t102 * qJD(1);
t523 = t690 / 0.2e1;
t522 = t296 * t630;
t521 = mrSges(5,3) * t691;
t361 = t444 * t630;
t308 = -t395 * t445 - t361;
t519 = t308 * t693;
t517 = mrSges(7,1) / 0.2e1 + t689;
t513 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t512 = t680 + t679;
t511 = -t87 / 0.2e1 + t677;
t509 = mrSges(7,1) + t690;
t505 = t608 / 0.2e1;
t504 = -t608 / 0.2e1;
t503 = -t607 / 0.2e1;
t500 = t584 / 0.2e1;
t499 = t583 / 0.2e1;
t498 = -t572 / 0.2e1;
t491 = t669 + t668;
t370 = t444 * mrSges(7,2);
t486 = t370 - t573;
t371 = t445 * mrSges(6,1);
t485 = -mrSges(6,2) * t444 + t371;
t481 = t342 * t486 - t396 * t485 + t718;
t478 = m(7) * t682 + t667;
t477 = Ifges(5,3) / 0.2e1 + t513;
t476 = (t688 + t687) * t297;
t475 = t513 * t383;
t474 = t641 + t629 / 0.2e1;
t473 = t370 - t485;
t472 = t622 + t623;
t469 = mrSges(6,1) * t444 + mrSges(6,2) * t445;
t452 = -t593 + t595;
t262 = t297 * t630;
t409 = t517 * t295 + t476 + (t570 / 0.2e1 + t568 / 0.2e1) * t381 + (t295 * t629 + t262) * t695 + (t295 * t395 + t262) * t693;
t413 = (-t743 * t444 - t749 * t445) * t696 + (-t750 * t444 - t445 * t759) * t694 + t536 / 0.2e1 - t533 / 0.2e1;
t425 = (-t296 * t516 + t491) * t444;
t10 = -t428 * t445 + t409 + t413 - t425;
t449 = t10 * qJD(1);
t421 = t476 + (t523 + t517) * t295;
t12 = -t719 * t445 + t421 - t425;
t448 = t12 * qJD(1);
t416 = t745 * mrSges(7,3) + (-t291 * t296 + t294 * t720 - t444 * t57 - t445 * t48) * t693 + t532;
t426 = t158 * t694 + t585 / 0.2e1 - t582 / 0.2e1;
t21 = t416 + t426;
t95 = m(7) * (-t291 * t444 - t445 * t720) + (t444 ^ 2 + t445 ^ 2) * mrSges(7,3);
t447 = t21 * qJD(1) + t95 * qJD(3);
t145 = (t308 / 0.4e1 - t622 / 0.4e1 - t623 / 0.4e1) * t700 + t486;
t59 = (-t522 / 0.4e1 + t539 / 0.4e1 - t210 / 0.4e1) * t700 - t487;
t443 = qJD(1) * t59 + qJD(3) * t145;
t104 = t294 * t509 + t280;
t251 = -t445 * t509 + t370;
t442 = qJD(1) * t104 + qJD(3) * t251;
t441 = t692 + t474;
t440 = Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1 - Ifges(6,2) / 0.2e1 - Ifges(7,2) / 0.2e1;
t438 = pkin(4) * t469;
t434 = (t401 * t83 + t404 * t82) * t695;
t432 = t469 * t648;
t417 = mrSges(7,1) * t682 + t295 * t514 + t297 * t515 + t61 * t687 + t688 * t83 + t689 * t82;
t406 = (t720 * t656 - t557 / 0.2e1) * mrSges(7,3) + (t709 * t656 - t549 / 0.2e1) * mrSges(6,3) - t370 * t716 - t371 * t665 + t417 + t702;
t414 = t591 / 0.2e1 + (-Ifges(6,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t381 + t758 + (t635 / 0.2e1 + t670) * pkin(5) - t440 * t296;
t418 = t757 + (-Ifges(6,5) / 0.4e1 - Ifges(7,5) / 0.4e1) * t381 + t440 * t294 - t590 / 0.2e1;
t419 = t210 * t653 + t668 * t709 + t669 * t720 + t755;
t422 = t342 * t210 - t291 * t759 + t750 * t720;
t424 = -t749 * t336 + t743 * t709;
t427 = t401 * t492 + t404 * t666;
t1 = -(mrSges(6,3) * t511 + mrSges(7,3) * t512 + t414) * t445 + t406 + t419 + (-0.3e1 / 0.4e1 * t596 + mrSges(5,2) * t652 - t237 / 0.4e1 + t319 * t691 + (-t438 / 0.2e1 + t673 / 0.2e1 + 0.3e1 / 0.2e1 * t606 + pkin(4) * t396 * t696 + (t521 + Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.4e1) * t405) * t383) * t405 + t395 * t667 + (pkin(4) * t558 + t424) * t696 + pkin(4) * t434 - ((-t66 / 0.2e1 + t683) * mrSges(7,3) - t510 * mrSges(6,3) + t418) * t444 + ((-t672 / 0.2e1 + (t521 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.4e1) * t402) * t402 + t477) * t383 + (t235 / 0.4e1 + mrSges(5,1) * t652 + 0.3e1 / 0.4e1 * t594 + t317 * t691 + (-t167 / 0.2e1 + t671) * pkin(4)) * t402 + t427 * pkin(4) + (pkin(4) * t571 + t395 * t51) * t693 + (t157 * t622 + t422) * t694 + t707;
t19 = -(t336 * mrSges(6,3) + mrSges(7,3) * t291 - t738) * t445 - m(6) * t396 * t622 + (t672 - t605) * t405 - (t445 * t488 + t737) * t444 + t481 + (-t438 + t606 + t673 + (-Ifges(5,1) + Ifges(5,2)) * t405) * t402 - t744 * t472;
t431 = -t1 * qJD(1) - t19 * qJD(3);
t22 = t744 * t623 + t726 * t643 + t724 * t645 + t717 * t647 - t481 + t718;
t420 = t214 * t663 + t215 * t650 + t755;
t5 = t406 + t475 + (pkin(5) * t682 + t48 * t754 + t548 * t692 + t56 * t662 + t57 * t734) * m(7) - t414 * t445 + (t667 + t554 / 0.2e1) * pkin(5) - ((t681 + t683) * mrSges(7,3) + t418) * t444 + t420;
t430 = -t5 * qJD(1) + t22 * qJD(3);
t429 = t291 * t629;
t142 = t441 * t634;
t261 = t613 * t629 + (t727 * m(7) + mrSges(6,1) + mrSges(7,1)) * t630;
t32 = t734 * mrSges(7,2) + (t650 + t709 / 0.2e1) * mrSges(6,2) + (t662 + t754) * mrSges(7,1) + (-t336 / 0.2e1 + t651) * mrSges(6,1) + (pkin(5) * t754 - t555 / 0.2e1 + t429 / 0.2e1) * m(7) - t441 * t608;
t412 = ((t674 / 0.2e1 + t500 - t491) * t404 + t719 * t401) * pkin(4);
t9 = -(t692 + t641) * t583 + (t681 + t66 / 0.2e1) * mrSges(7,2) + t510 * mrSges(6,2) + t512 * mrSges(7,1) + t511 * mrSges(6,1) + (pkin(5) * t679 + t57 * t641) * m(7) + t412;
t423 = t9 * qJD(1) + t142 * qJD(2) + t32 * qJD(3) - t261 * qJD(4);
t407 = t720 * t499 + t709 * t500 + t48 * t505 + t485 * t665 + t486 * t716 + t79 * t498 + t57 * t503 + t549 * t685 + t557 * t684 + t417 - t702 - t724 * t294 / 0.4e1 + t726 * t294 / 0.4e1 + (-t296 * t751 + t758) * t445 / 0.4e1 - (t713 * t381 + t725) * t445 / 0.4e1 - t717 * t296 / 0.4e1 - (t731 * t294 + t714 * t381 + t722 - t757) * t444 / 0.4e1;
t199 = t472 * t693 + t519;
t103 = (-t542 - t543) * t693 + m(7) * t648;
t96 = -t727 * t634 / 0.2e1 - (mrSges(7,1) + t523) * t445 + t473;
t91 = (-t522 + t539 + t210) * t693;
t31 = pkin(5) * t505 - t291 * t523 + (t429 - t555) * t693 - t474 * t608 + t756;
t20 = t416 - t426;
t16 = t381 * t704 - t490 * t444 + t492 * t445 - t408 + t410 + t432 + t752 * (t295 * t643 + t297 * t647);
t13 = t643 * t715 + t421 + t701;
t11 = t409 - t413 + t701;
t8 = -t614 / 0.2e1 + t615 / 0.2e1 - t621 / 0.2e1 - t620 / 0.2e1 - t617 / 0.2e1 - t616 / 0.2e1 - t618 / 0.2e1 + t619 / 0.2e1 + t65 * t523 + pkin(5) * t499 + (-t674 / 0.2e1 + t499) * t395 + t412 + t480;
t4 = t407 + t475 - t610 * t291 * t693 + t56 * t504 - t623 * t671 - t572 * t677 - t607 * t680 - t420 + ((t157 * t445 - t548) * t693 - t554 / 0.2e1 + t478) * pkin(5);
t2 = t707 + t407 - t419 + (t571 * t693 + t427 + t434) * pkin(4) - t319 * t624 / 0.2e1 - t317 * t625 / 0.2e1 + t167 * t622 / 0.2e1 + t65 * t503 + t66 * t504 + t87 * t498 - t467 * t545 / 0.2e1 - t459 * t544 / 0.2e1 + ((t396 * t544 + t558) * pkin(4) + t424) * t695 + (t157 * t472 + t422) * t693 + t472 * t670 + t432 * t628 + t626 * t704 + t534 / 0.4e1 + t706 * t631 / 0.2e1 - t510 * t574 + t477 * t383 + t478 * t395 + t327 * t470 / 0.2e1 + t381 * t452 / 0.4e1 - t537 / 0.4e1 + t433 - t402 * t435 / 0.4e1 + t405 * t436 / 0.4e1;
t17 = [qJD(2) * t14 + qJD(3) * t3 + qJD(4) * t6 + qJD(5) * t7 + qJD(6) * t30, t561 + 0.2e1 * (t693 + t695) * (-t295 * t444 + t297 * t445) * qJD(2) + t16 * qJD(3) + t11 * qJD(4) + t13 * qJD(5) + t103 * qJD(6), t16 * qJD(2) + t2 * qJD(4) + t4 * qJD(5) + t20 * qJD(6) + t581 + (t726 * t657 + t728 * t645 + t729 * t647 + t721 * mrSges(5,3) + 0.2e1 * (-pkin(3) * t328 + t721 * pkin(8)) * t697 + t723 * t654 + (-t444 * t713 + t445 * t714 + t451) * t648 + t709 * t217 + 0.2e1 * (t246 * t396 + t336 * t83 + t709 * t82) * t695 + t396 * t165 - Ifges(4,6) * t383 + t336 * t213 + t342 * t164 + t327 * mrSges(4,2) - t328 * mrSges(4,1) + t158 * t322 + pkin(3) * t309 + t291 * t212 - t318 * t625 + t316 * t624 + t234 * t637 + t236 * t639 - t51 * t607 - t61 * t608 - t82 * t572 - t83 * t574 + t328 * t706 + t720 * t216 + 0.2e1 * (t158 * t342 + t291 * t61 + t51 * t720) * t693 + t246 * t469 + (-Ifges(4,5) + t703) * t381) * qJD(3), t567 + t11 * qJD(2) + t2 * qJD(3) + (-t183 * mrSges(5,1) - t182 * mrSges(5,2) - Ifges(5,5) * t545 - Ifges(5,6) * t544 + t395 * t583 + t633 * t65 + t480 - t614 + t615 - t618 + t619) * qJD(4) + t8 * qJD(5) + t91 * qJD(6) + (mrSges(7,3) * t552 + (t296 * t404 + t552) * mrSges(6,3) + t66 * t632 + m(6) * (t401 * t88 + t404 * t87)) * t592, t566 + t13 * qJD(2) + t4 * qJD(3) + t8 * qJD(4) + (-t616 - t620 - t617 - t621 + (t583 - t674) * pkin(5) + t480) * qJD(5), qJD(2) * t103 + qJD(3) * t20 + qJD(4) * t91 + t562; qJD(3) * t15 - qJD(4) * t10 - qJD(5) * t12 - qJD(6) * t102 - t561, 0, t563 (-t470 + t473 - t573) * qJD(4) + t96 * qJD(5) + 0.2e1 * ((-t445 * t629 - t361) * t695 + t519) * qJD(4) - t449, t96 * qJD(4) + (t251 - t485) * qJD(5) - t448, -t525; -qJD(2) * t15 - qJD(4) * t1 - qJD(5) * t5 + qJD(6) * t21 - t581, -t563, -qJD(4) * t19 + qJD(5) * t22 + qJD(6) * t95 (pkin(8) * t706 - t291 * t633 + t395 * t608 + t452 + t756) * qJD(4) + t31 * qJD(5) + t199 * qJD(6) + (-mrSges(7,3) * t541 + (t404 * t444 - t541) * mrSges(6,3) + t720 * t632 + m(6) * (-t336 * t404 + t401 * t709)) * t592 + t431, t31 * qJD(4) + ((-m(7) * t291 + t608) * pkin(5) + t756) * qJD(5) + t430, t199 * qJD(4) + t447; qJD(2) * t10 + qJD(3) * t1 + qJD(5) * t9 + qJD(6) * t59 - t567, qJD(5) * t142 + t449, qJD(5) * t32 + qJD(6) * t145 - t431, -t261 * qJD(5) (-t613 * t404 + (-mrSges(6,1) - t509) * t401) * qJD(5) * pkin(4) + t423, t443; qJD(2) * t12 + qJD(3) * t5 - qJD(4) * t9 + qJD(6) * t104 - t566, -qJD(4) * t142 + t448, -qJD(4) * t32 + qJD(6) * t251 - t430, -t423, 0, t442; qJD(2) * t102 - qJD(3) * t21 - qJD(4) * t59 - qJD(5) * t104 - t562, t525, -t145 * qJD(4) - t251 * qJD(5) - t447, -t443, -t442, 0;];
Cq  = t17;
