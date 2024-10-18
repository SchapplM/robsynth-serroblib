% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR15_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:24:08
% EndTime: 2024-09-27 22:24:26
% DurationCPUTime: 10.75s
% Computational Cost: add. (55338->666), mult. (148374->906), div. (0->0), fcn. (165064->12), ass. (0->379)
t420 = cos(qJ(4));
t652 = t420 * pkin(3);
t407 = pkin(4) + t652;
t411 = sin(pkin(6));
t412 = sin(pkin(5));
t417 = sin(qJ(3));
t418 = sin(qJ(2));
t421 = cos(qJ(3));
t422 = cos(qJ(2));
t378 = (-t417 * t418 + t421 * t422) * t412;
t379 = (-t417 * t422 - t418 * t421) * t412;
t416 = sin(qJ(4));
t318 = t378 * t416 - t379 * t420;
t419 = cos(qJ(5));
t523 = t420 * t378 + t379 * t416;
t413 = cos(pkin(6));
t415 = sin(qJ(5));
t567 = t413 * t415;
t722 = -t318 * t567 + t419 * t523;
t566 = t413 * t419;
t723 = -t318 * t566 - t415 * t523;
t514 = mrSges(6,1) * t723 - mrSges(6,2) * t722;
t760 = t411 * t514;
t763 = t760 / 0.2e1;
t764 = t407 * t763;
t744 = Ifges(6,6) * t723;
t745 = Ifges(6,5) * t722;
t492 = t745 / 0.2e1 + t744 / 0.2e1;
t728 = t318 * t411;
t746 = t728 / 0.2e1;
t697 = Ifges(6,3) * t746 + t492;
t710 = Ifges(5,5) * t523;
t729 = Ifges(5,6) * t318;
t638 = Ifges(6,4) * t415;
t373 = t413 * Ifges(6,6) + (Ifges(6,2) * t419 + t638) * t411;
t739 = t723 * t373;
t570 = t411 * t419;
t403 = Ifges(6,4) * t570;
t571 = t411 * t415;
t374 = Ifges(6,1) * t571 + t413 * Ifges(6,5) + t403;
t740 = t722 * t374;
t734 = Ifges(6,1) * t722 + Ifges(6,4) * t723 + Ifges(6,5) * t728;
t752 = t734 * t571;
t733 = Ifges(6,4) * t722 + Ifges(6,2) * t723 + Ifges(6,6) * t728;
t753 = t733 * t570;
t750 = Ifges(6,3) * t728 + t744 + t745;
t759 = t413 * t750;
t762 = t710 / 0.2e1 - t729 / 0.2e1 + t739 / 0.4e1 + t740 / 0.4e1 + t752 / 0.4e1 + t753 / 0.4e1 + t759 / 0.4e1;
t372 = Ifges(6,3) * t413 + (Ifges(6,5) * t415 + Ifges(6,6) * t419) * t411;
t761 = t372 * t746 + t739 / 0.2e1 + t740 / 0.2e1 + t752 / 0.2e1 + t753 / 0.2e1 + t759 / 0.2e1;
t414 = cos(pkin(5));
t589 = t523 * t413;
t500 = t411 * t414 + t589;
t229 = -t318 * t415 + t500 * t419;
t230 = t318 * t419 + t500 * t415;
t590 = t523 * t411;
t287 = t413 * t414 - t590;
t541 = -pkin(2) * t422 - pkin(1);
t520 = t541 * t412;
t476 = t378 * pkin(3) - t520;
t639 = Ifges(6,4) * t230;
t477 = Ifges(6,2) * t229 + Ifges(6,6) * t287 + t639;
t228 = Ifges(6,4) * t229;
t493 = Ifges(6,1) * t230 + Ifges(6,5) * t287 + t228;
t659 = pkin(1) * t414;
t404 = t422 * t659;
t569 = t412 * t418;
t686 = pkin(8) + pkin(9);
t361 = -t686 * t569 + t404;
t345 = pkin(2) * t414 + t361;
t555 = t418 * t659;
t568 = t412 * t422;
t362 = t686 * t568 + t555;
t349 = t417 * t362;
t286 = t421 * t345 - t349;
t375 = t379 * pkin(10);
t252 = t286 + t375;
t245 = pkin(3) * t414 + t252;
t351 = t421 * t362;
t499 = t345 * t417 + t351;
t655 = pkin(10) * t378;
t253 = t499 + t655;
t248 = t420 * t253;
t501 = t245 * t416 + t248;
t246 = t416 * t253;
t183 = t420 * t245 - t246;
t654 = pkin(11) * t413;
t303 = t318 * t654;
t141 = t183 - t303;
t138 = pkin(4) * t414 + t141;
t209 = -pkin(4) * t523 - pkin(11) * t728 - t476;
t506 = t138 * t411 - t209 * t413;
t131 = t500 * pkin(11) + t501;
t507 = t138 * t413 + t209 * t411;
t57 = -t131 * t415 + t507 * t419;
t58 = t131 * t419 + t507 * t415;
t640 = Ifges(5,4) * t318;
t677 = -t523 / 0.2e1;
t682 = -t230 / 0.2e1;
t732 = -t318 / 0.2e1;
t736 = -mrSges(6,2) * t728 + mrSges(6,3) * t723;
t737 = mrSges(6,1) * t728 - mrSges(6,3) * t722;
t756 = -t229 * t733 / 0.2e1 + (t476 * mrSges(5,1) + mrSges(5,3) * t501 + Ifges(5,6) * t414 / 0.2e1 + t640 / 0.2e1 + (-t677 + t523 / 0.2e1) * Ifges(5,2)) * t318 - t57 * t737 - t58 * t736 + t734 * t682 - t493 * t722 / 0.2e1 - t723 * t477 / 0.2e1 + t476 * mrSges(5,2) * t523 + (Ifges(5,1) * t318 + 0.2e1 * Ifges(5,4) * t523 + Ifges(5,5) * t414) * t677 - t750 * t287 / 0.2e1 + (-t640 + t411 * (Ifges(6,5) * t230 + Ifges(6,6) * t229 + Ifges(6,3) * t287) + Ifges(5,1) * t523) * t732 - t506 * t514;
t647 = mrSges(5,3) * t318;
t394 = pkin(4) * t567 + pkin(11) * t570;
t754 = t394 * t736;
t542 = Ifges(6,3) * t732;
t742 = t372 * t728;
t658 = pkin(2) * t421;
t408 = pkin(3) + t658;
t563 = t417 * t420;
t386 = pkin(2) * t563 + t408 * t416;
t741 = t386 * t647;
t133 = Ifges(6,5) * t229 - Ifges(6,6) * t230;
t388 = (mrSges(6,1) * t415 + mrSges(6,2) * t419) * t411;
t389 = Ifges(6,5) * t570 - Ifges(6,6) * t571;
t398 = mrSges(6,1) * t413 - mrSges(6,3) * t571;
t399 = -mrSges(6,2) * t413 + mrSges(6,3) * t570;
t687 = -t58 / 0.2e1;
t442 = (t133 / 0.4e1 + t209 * t388 / 0.2e1) * t413 + t287 * t389 / 0.4e1 + t57 * t399 / 0.2e1 + t398 * t687;
t511 = -t729 + t710;
t730 = pkin(4) * t318;
t713 = mrSges(5,3) * t523;
t727 = -t742 / 0.4e1;
t726 = -t741 / 0.2e1;
t177 = -mrSges(6,2) * t287 + mrSges(6,3) * t229;
t178 = mrSges(6,1) * t287 - mrSges(6,3) * t230;
t188 = -t252 * t416 - t248;
t189 = t420 * t252 - t246;
t409 = t411 * pkin(11);
t365 = t386 + t409;
t564 = t416 * t417;
t385 = -pkin(2) * t564 + t420 * t408;
t380 = pkin(4) + t385;
t301 = -t365 * t415 + t380 * t566;
t302 = t365 * t419 + t380 * t567;
t396 = (-t416 * t421 - t563) * pkin(2);
t397 = (t420 * t421 - t564) * pkin(2);
t325 = t396 * t566 - t397 * t415;
t326 = t396 * t567 + t397 * t419;
t515 = mrSges(6,1) * t229 - mrSges(6,2) * t230;
t474 = -t411 * t515 / 0.2e1;
t482 = t411 * t506;
t483 = t499 * mrSges(4,1);
t512 = mrSges(5,2) * t414 - t713;
t618 = t378 * mrSges(4,3);
t513 = mrSges(4,2) * t414 - t618;
t649 = mrSges(4,3) * t379;
t518 = mrSges(4,1) * t414 + t649;
t609 = t419 * mrSges(6,1);
t610 = t415 * mrSges(6,2);
t387 = (-t609 + t610) * t411;
t556 = t523 * t654;
t151 = t188 - t556;
t657 = pkin(3) * t379;
t225 = -t409 * t523 - t657 + t730;
t96 = -t151 * t411 + t225 * t413;
t602 = t96 * t387;
t152 = -t303 + t189;
t503 = t151 * t413 + t225 * t411;
t72 = t152 * t419 + t503 * t415;
t606 = t72 * t399;
t71 = -t152 * t415 + t503 * t419;
t607 = t71 * t398;
t629 = t189 * mrSges(5,2);
t630 = t188 * mrSges(5,1);
t519 = t630 / 0.2e1 - t629 / 0.2e1 + t607 / 0.2e1 + t606 / 0.2e1 + t602 / 0.2e1;
t548 = t380 * t760;
t552 = t385 * t713;
t580 = t380 * t411;
t581 = t379 * t417;
t592 = t302 * t736;
t593 = t301 * t737;
t626 = t286 * mrSges(4,2);
t672 = -t326 / 0.2e1;
t516 = mrSges(5,1) * t414 - t647;
t715 = -t516 / 0.2e1;
t725 = -m(5) * (t396 * t183 + t385 * t188 + t386 * t189 + t397 * t501) / 0.2e1 - m(6) * (t301 * t71 + t302 * t72 + t325 * t57 + t326 * t58 + t396 * t482 - t96 * t580) / 0.2e1 + t626 / 0.2e1 - t593 / 0.2e1 - t592 / 0.2e1 - t325 * t178 / 0.2e1 + t177 * t672 + t396 * t715 + t397 * t512 / 0.2e1 + t483 / 0.2e1 + t552 / 0.2e1 - t726 - t548 / 0.2e1 + t396 * t474 + (t513 + t618) * t658 / 0.2e1 + (t417 * t518 / 0.2e1 - t581 * mrSges(4,3) / 0.2e1) * pkin(2) - t519;
t390 = -Ifges(6,2) * t571 + t403;
t529 = t374 / 0.4e1 + t390 / 0.4e1;
t391 = (Ifges(6,1) * t419 - t638) * t411;
t531 = -t373 / 0.4e1 + t391 / 0.4e1;
t614 = t394 * mrSges(6,3);
t392 = pkin(4) * t566 - pkin(11) * t571;
t615 = t392 * mrSges(6,3);
t663 = -t394 / 0.2e1;
t664 = t392 / 0.2e1;
t724 = (-t615 / 0.2e1 + t529) * t229 + (-t614 / 0.2e1 + t531) * t230 + t177 * t664 + t178 * t663 + t442;
t661 = -t414 / 0.2e1;
t716 = t736 / 0.2e1;
t711 = Ifges(3,5) * t568;
t709 = Ifges(3,6) * t569;
t572 = t411 * t387;
t708 = mrSges(5,1) - t572;
t549 = -t652 / 0.2e1;
t705 = t549 * t713;
t704 = -Ifges(6,2) * t230 + t228 + t493;
t293 = -t361 * t417 - t351;
t258 = t293 - t655;
t294 = t421 * t361 - t349;
t259 = t375 + t294;
t192 = t420 * t258 - t259 * t416;
t164 = t192 - t556;
t558 = pkin(2) * t569;
t217 = t225 + t558;
t97 = -t164 * t411 + t217 * t413;
t601 = t97 * t387;
t193 = t416 * t258 + t420 * t259;
t165 = -t303 + t193;
t502 = t164 * t413 + t217 * t411;
t78 = t165 * t419 + t502 * t415;
t604 = t78 * t399;
t77 = -t165 * t415 + t502 * t419;
t605 = t77 * t398;
t627 = t193 * mrSges(5,2);
t628 = t192 * mrSges(5,1);
t444 = t628 / 0.2e1 - t627 / 0.2e1 + t605 / 0.2e1 + t604 / 0.2e1 + t601 / 0.2e1;
t656 = pkin(4) * t411;
t691 = m(6) / 0.2e1;
t703 = t444 + (t392 * t77 + t394 * t78 - t97 * t656) * t691;
t702 = (t392 * t71 + t394 * t72 - t96 * t656) * t691 + t519;
t701 = t512 + t713;
t699 = Ifges(4,5) * t378 + Ifges(4,6) * t379 + t511;
t584 = t326 * t399;
t585 = t325 * t398;
t612 = t397 * mrSges(5,2);
t696 = t584 + t585 - t612 + (-mrSges(4,1) * t417 - mrSges(4,2) * t421) * pkin(2);
t695 = (pkin(8) * t568 + t555) * mrSges(3,1) + (-pkin(8) * t569 + t404) * mrSges(3,2);
t142 = -pkin(11) * t589 - t501;
t241 = -pkin(11) * t590 + t730;
t505 = t142 * t413 + t241 * t411;
t599 = t141 * t419;
t458 = t505 * t415 + t599;
t446 = t458 * t399;
t600 = t141 * t415;
t459 = t505 * t419 - t600;
t447 = t459 * t398;
t504 = t142 * t411 - t241 * t413;
t480 = t504 * t387;
t484 = t501 * mrSges(5,1);
t631 = t183 * mrSges(5,2);
t694 = t446 / 0.2e1 + t447 / 0.2e1 - t480 / 0.2e1 - t484 / 0.2e1 - t631 / 0.2e1 + t762;
t693 = t412 ^ 2;
t692 = m(5) / 0.2e1;
t690 = -pkin(4) / 0.2e1;
t689 = -mrSges(6,1) / 0.2e1;
t688 = mrSges(6,2) / 0.2e1;
t680 = t301 / 0.2e1;
t679 = -t302 / 0.2e1;
t319 = -t385 * t415 - t386 * t566;
t675 = t319 / 0.2e1;
t320 = t385 * t419 - t386 * t567;
t674 = -t320 / 0.2e1;
t673 = t325 / 0.2e1;
t653 = t416 * pkin(3);
t400 = t409 + t653;
t347 = -t400 * t415 + t407 * t566;
t671 = t347 / 0.2e1;
t348 = t400 * t419 + t407 * t567;
t670 = -t348 / 0.2e1;
t669 = -t378 / 0.2e1;
t668 = -t380 / 0.2e1;
t383 = (-t415 * t420 - t416 * t566) * pkin(3);
t667 = -t383 / 0.2e1;
t384 = (-t416 * t567 + t419 * t420) * pkin(3);
t666 = t384 / 0.2e1;
t662 = -t407 / 0.2e1;
t645 = Ifges(3,4) * t418;
t644 = Ifges(3,4) * t422;
t625 = t293 * mrSges(4,1);
t624 = t294 * mrSges(4,2);
t623 = t301 * mrSges(6,3);
t622 = t302 * mrSges(6,3);
t620 = t347 * mrSges(6,3);
t619 = t348 * mrSges(6,3);
t617 = t385 * mrSges(5,2);
t616 = t386 * mrSges(5,1);
t613 = t396 * mrSges(5,1);
t423 = -(-mrSges(4,1) * t379 + mrSges(4,2) * t378) * t520 - t499 * t649 + t286 * t618 + t183 * t713 + (0.2e1 * Ifges(4,4) * t378 + Ifges(4,5) * t414) * t669 + (Ifges(4,4) * t379 + Ifges(4,6) * t661 + (Ifges(4,1) - Ifges(4,2)) * (t378 / 0.2e1 - t669)) * t379 + t699 * t661 + t756;
t498 = t558 - t657;
t517 = -mrSges(5,1) * t523 + mrSges(5,2) * t318;
t4 = -(-mrSges(4,1) * t378 - mrSges(4,2) * t379) * t558 - m(4) * (t541 * pkin(2) * t418 * t693 + t286 * t293 + t499 * t294) - m(5) * (t183 * t192 + t501 * t193 - t476 * t498) - t498 * t517 - m(6) * (-t506 * t97 + t57 * t77 + t58 * t78) + t693 * pkin(1) * (mrSges(3,1) * t418 + mrSges(3,2) * t422) + t294 * t513 - t293 * t518 + t193 * t512 - t192 * t516 + t97 * t515 - t77 * t178 - t78 * t177 + t423 - (t418 * (Ifges(3,1) * t422 - t645) + t422 * (-Ifges(3,2) * t418 + t644)) * t693 / 0.2e1 + (-(Ifges(3,1) * t418 + t644) * t568 / 0.2e1 + (Ifges(3,2) * t422 + t645) * t569 / 0.2e1 + (Ifges(3,5) * t422 - Ifges(3,6) * t418) * t661) * t412 + (-t711 / 0.2e1 + t709 / 0.2e1 + t695) * t414;
t611 = t4 * qJD(1);
t7 = t517 * t657 + t499 * t518 - m(6) * (-t506 * t96 + t57 * t71 + t58 * t72) - m(5) * (t183 * t188 + t501 * t189 + t476 * t657) + t286 * t513 + t189 * t512 - t188 * t516 + t96 * t515 - t71 * t178 - t72 * t177 + t423;
t608 = t7 * qJD(1);
t9 = t511 * t661 + t501 * t516 - t504 * t515 - m(6) * (t58 * t458 + t57 * t459 + t506 * t504) - t458 * t177 - t459 * t178 + t701 * t183 + t756;
t603 = t9 * qJD(1);
t132 = mrSges(6,1) * t230 + mrSges(6,2) * t229;
t135 = Ifges(6,1) * t229 - t639;
t16 = -t506 * t132 + t57 * t177 - t58 * t178 + t230 * t135 / 0.2e1 + t287 * t133 / 0.2e1 + t477 * t682 + (-t229 * t57 - t230 * t58) * mrSges(6,3) + t704 * t229 / 0.2e1;
t598 = t16 * qJD(1);
t587 = t319 * t398;
t586 = t320 * t399;
t583 = t347 * t737;
t582 = t348 * t736;
t579 = t383 * t398;
t578 = t384 * t399;
t410 = t411 ^ 2;
t577 = t386 * t410;
t576 = t396 * t410;
t575 = t407 * t411;
t565 = t416 * t318;
t562 = t420 * t523;
t557 = t410 * t653;
t550 = -t653 / 0.2e1;
t547 = t386 * t572;
t546 = t396 * t572;
t535 = t413 * t389 / 0.2e1;
t534 = -t565 / 0.2e1;
t533 = t671 + t680;
t532 = t670 + t679;
t530 = t374 / 0.2e1 + t390 / 0.2e1;
t528 = t391 / 0.2e1 - t373 / 0.2e1;
t527 = t664 + t680;
t526 = t664 + t671;
t525 = t663 + t679;
t524 = t663 + t670;
t522 = t572 * t653;
t429 = t737 * t664 + t754 / 0.2e1 + pkin(4) * t763 - t727 + t762;
t424 = t429 - t694 + t727;
t481 = t411 * t504;
t427 = (t301 * t459 + t302 * t458 + t319 * t57 + t320 * t58 + t380 * t481) * t691 + t737 * t680 + t302 * t716 + t178 * t675 + t320 * t177 / 0.2e1 + t726 + t380 * t763 - t701 * t385 / 0.2e1 + (-t482 * t691 + t474 + t715) * t386;
t2 = -t427 + t424 + t703;
t445 = t547 + t586 + t587 - t616 - t617;
t80 = m(6) * (t301 * t319 + t302 * t320 - t380 * t577) + t445;
t510 = -t2 * qJD(1) + t80 * qJD(2);
t433 = (t347 * t77 + t348 * t78 - t97 * t575) * t691 + t625 / 0.2e1 - t624 / 0.2e1 + t583 / 0.2e1 + t582 / 0.2e1 + t764 + t444;
t475 = (t192 * t420 + t193 * t416) * t692;
t8 = t705 + t433 + (mrSges(5,3) * t534 + t475) * pkin(3) + t725;
t81 = t708 * t396 + m(6) * (t301 * t325 + t302 * t326 + t380 * t576) + m(5) * (t385 * t396 + t386 * t397) + t696;
t509 = -t8 * qJD(1) + t81 * qJD(2);
t439 = -t138 * t388 / 0.2e1 + (t687 * mrSges(6,3) + t135 / 0.4e1 - t477 / 0.4e1) * t415 + (-t57 * mrSges(6,3) / 0.2e1 + t704 / 0.4e1) * t419;
t435 = t132 * t668 + t439;
t437 = (-t623 / 0.2e1 + t529) * t229 + (-t622 / 0.2e1 + t531) * t230 + t177 * t680 + t178 * t679 + t442;
t496 = t78 * t688 + t77 * t689;
t10 = (t542 + t435) * t411 + t437 + t496 - t492;
t54 = t301 * t399 - t302 * t398 + t535 + (-t380 * t388 + (t530 - t623) * t419 + (t528 - t622) * t415) * t411;
t508 = t10 * qJD(1) + t54 * qJD(2);
t497 = t72 * t688 + t71 * t689;
t491 = mrSges(6,1) * t675 + mrSges(6,2) * t674;
t490 = mrSges(6,1) * t673 + mrSges(6,2) * t672;
t489 = mrSges(6,1) * t667 + mrSges(6,2) * t666;
t434 = t132 * t662 + t439;
t436 = (-t620 / 0.2e1 + t529) * t229 + (-t619 / 0.2e1 + t531) * t230 + t177 * t671 + t178 * t670 + t442;
t12 = (t542 + t434) * t411 + t436 + t497 - t492;
t453 = t528 * t415 + t530 * t419;
t430 = t533 * t399 + t532 * t398 + ((t668 + t662) * t388 + (t532 * t415 - t533 * t419) * mrSges(6,3) + t453) * t411 + t535;
t31 = t430 - t490;
t82 = t347 * t399 - t348 * t398 + t535 + (-t407 * t388 + (t530 - t620) * t419 + (t528 - t619) * t415) * t411;
t471 = t12 * qJD(1) + t31 * qJD(2) + t82 * qJD(3);
t443 = -mrSges(5,1) * t653 - mrSges(5,2) * t652 + t522 + t578 + t579;
t163 = m(6) * (t347 * t383 + t348 * t384 - t407 * t557) + t443;
t441 = m(6) * (t301 * t383 + t302 * t384 + t319 * t347 + t320 * t348 + (-t380 * t653 - t386 * t407) * t410);
t448 = (pkin(4) * t576 + t325 * t392 + t326 * t394) * t691;
t29 = t448 - t441 / 0.2e1 + (t326 / 0.2e1 - t384 / 0.2e1 + t674) * t399 + (t673 + t667 - t319 / 0.2e1) * t398 + (-t397 / 0.2e1 + t385 / 0.2e1 + t652 / 0.2e1) * mrSges(5,2) + t708 * (t396 / 0.2e1 + t386 / 0.2e1 + t653 / 0.2e1);
t426 = (t347 * t459 + t348 * t458 + t383 * t57 + t384 * t58 + t407 * t481 - t482 * t653) * t691 + t737 * t671 + t348 * t716 + t383 * t178 / 0.2e1 + t177 * t666 + t764 + t512 * t549 + t474 * t653 + t705 + (t516 + t647) * t550;
t6 = t424 - t426 + t702;
t470 = -t6 * qJD(1) - t29 * qJD(2) + t163 * qJD(3);
t457 = t699 + t761;
t104 = t392 * t399 - t394 * t398 + t535 + (-pkin(4) * t388 + (t530 - t615) * t419 + (t528 - t614) * t415) * t411;
t438 = t132 * t690 + t439;
t14 = (t142 * t567 + t599) * t688 + (t142 * t566 - t600) * t689 + (t542 + (t610 / 0.2e1 - t609 / 0.2e1) * t241 + t438) * t411 - t492 + t724;
t432 = t527 * t399 + t525 * t398 + ((t668 + t690) * t388 + (t525 * t415 - t527 * t419) * mrSges(6,3) + t453) * t411 + t535;
t33 = t432 - t491;
t431 = t526 * t399 + t524 * t398 + ((t662 + t690) * t388 + (t524 * t415 - t526 * t419) * mrSges(6,3) + t453) * t411 + t535;
t51 = t431 + t489;
t454 = t14 * qJD(1) + t33 * qJD(2) + t51 * qJD(3) + t104 * qJD(4);
t428 = t694 + t742 / 0.4e1 + t429;
t50 = t431 - t489;
t34 = t432 + t491;
t32 = t430 + t490;
t30 = t547 / 0.2e1 + t578 / 0.2e1 + t579 / 0.2e1 + t441 / 0.2e1 + t522 / 0.2e1 + t586 / 0.2e1 + t587 / 0.2e1 - t617 / 0.2e1 - t616 / 0.2e1 + mrSges(5,2) * t549 + mrSges(5,1) * t550 - t546 / 0.2e1 + t448 + t584 / 0.2e1 + t585 / 0.2e1 - t612 / 0.2e1 + t613 / 0.2e1;
t15 = -t458 * mrSges(6,2) / 0.2e1 + t459 * mrSges(6,1) / 0.2e1 + t438 * t411 + t724 + t697;
t13 = t434 * t411 + t436 - t497 + t697;
t11 = t435 * t411 + t437 - t496 + t697;
t5 = t428 + t426 + t702;
t3 = (t475 + (t534 - t562 / 0.2e1) * mrSges(5,3)) * pkin(3) + t457 + t433 - t725;
t1 = t428 + t427 + t703;
t17 = [-qJD(2) * t4 - qJD(3) * t7 - qJD(4) * t9 + qJD(5) * t16, t3 * qJD(3) + t1 * qJD(4) + t11 * qJD(5) - t611 + (t711 - t709 - t741 + t457 + t548 - t552 + t592 + t593 + t601 + t604 + t605 - t624 + t625 - t627 + t628 + 0.2e1 * (t192 * t385 + t193 * t386) * t692 + 0.2e1 * (t301 * t77 + t302 * t78 - t97 * t580) * t691 + (m(4) * (t293 * t421 + t294 * t417) + (-t378 * t421 + t581) * mrSges(4,3)) * pkin(2) - t695) * qJD(2), t3 * qJD(2) + t5 * qJD(4) + t13 * qJD(5) - t608 + (t407 * t760 - t483 + m(6) * (t347 * t71 + t348 * t72 - t96 * t575) + t607 + t606 + t602 + t583 + t582 - t626 + t630 - t629 + t457 + (m(5) * (t188 * t420 + t189 * t416) + (-t562 - t565) * mrSges(5,3)) * pkin(3)) * qJD(3), -t603 + t1 * qJD(2) + t5 * qJD(3) + (-t480 + m(6) * (pkin(4) * t481 + t392 * t459 + t394 * t458) + t754 + t392 * t737 + pkin(4) * t760 + t446 + t447 - t484 - t631 + t511 + t761) * qJD(4) + t15 * qJD(5), t598 + t11 * qJD(2) + t13 * qJD(3) + t15 * qJD(4) + (-mrSges(6,1) * t58 - mrSges(6,2) * t57 + t133) * qJD(5); -qJD(3) * t8 - qJD(4) * t2 + qJD(5) * t10 + t611, qJD(3) * t81 + qJD(4) * t80 + qJD(5) * t54, (-t546 + t613 + t696) * qJD(3) + t30 * qJD(4) + t32 * qJD(5) + 0.2e1 * ((t325 * t347 + t326 * t348 + t407 * t576) * t691 + (t396 * t420 + t397 * t416) * pkin(3) * t692) * qJD(3) + t509, t30 * qJD(3) + (m(6) * (-pkin(4) * t577 + t319 * t392 + t320 * t394) + t445) * qJD(4) + t34 * qJD(5) + t510, t32 * qJD(3) + t34 * qJD(4) + (-mrSges(6,1) * t302 - mrSges(6,2) * t301 + t389) * qJD(5) + t508; qJD(2) * t8 - qJD(4) * t6 + qJD(5) * t12 + t608, -qJD(4) * t29 + qJD(5) * t31 - t509, qJD(4) * t163 + qJD(5) * t82, (m(6) * (-pkin(4) * t557 + t383 * t392 + t384 * t394) + t443) * qJD(4) + t50 * qJD(5) + t470, t50 * qJD(4) + (-mrSges(6,1) * t348 - mrSges(6,2) * t347 + t389) * qJD(5) + t471; qJD(2) * t2 + qJD(3) * t6 + qJD(5) * t14 + t603, qJD(3) * t29 + qJD(5) * t33 - t510, qJD(5) * t51 - t470, t104 * qJD(5), (-mrSges(6,1) * t394 - mrSges(6,2) * t392 + t389) * qJD(5) + t454; -qJD(2) * t10 - qJD(3) * t12 - qJD(4) * t14 - t598, -qJD(3) * t31 - qJD(4) * t33 - t508, -qJD(4) * t51 - t471, -t454, 0;];
Cq = t17;
