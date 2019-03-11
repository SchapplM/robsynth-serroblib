% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRPP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:43:02
% EndTime: 2019-03-08 22:43:21
% DurationCPUTime: 12.11s
% Computational Cost: add. (14483->745), mult. (35428->1015), div. (0->0), fcn. (36860->10), ass. (0->377)
t455 = sin(qJ(4));
t456 = sin(qJ(3));
t459 = cos(qJ(3));
t553 = t459 * pkin(3) + pkin(2);
t458 = cos(qJ(4));
t579 = t458 * t459;
t563 = pkin(8) * t579;
t636 = -qJ(5) - pkin(9);
t475 = t455 * (t456 * t636 - t553) + t563;
t603 = sin(pkin(11));
t759 = t475 * t603;
t604 = cos(pkin(11));
t758 = t604 * t475;
t437 = -pkin(4) * t458 - pkin(3);
t526 = t603 * t455;
t481 = t604 * t458 - t526;
t528 = t604 * t455;
t482 = t458 * t603 + t528;
t241 = -pkin(5) * t481 - qJ(6) * t482 + t437;
t276 = -mrSges(7,1) * t481 - mrSges(7,3) * t482;
t757 = m(7) * t241 + t276;
t657 = t481 / 0.2e1;
t654 = t482 / 0.2e1;
t450 = t455 ^ 2;
t452 = t458 ^ 2;
t731 = t450 + t452;
t755 = pkin(9) * t731;
t742 = Ifges(7,4) + Ifges(6,5);
t754 = -Ifges(6,6) + Ifges(7,6);
t498 = -pkin(9) * t456 - t553;
t395 = t458 * t498;
t442 = t458 * qJ(5);
t510 = -t442 * t456 + t395;
t265 = (-pkin(8) * t455 - pkin(4)) * t459 + t510;
t127 = t603 * t265 + t758;
t111 = -qJ(6) * t459 + t127;
t126 = t265 * t604 - t759;
t117 = t459 * pkin(5) - t126;
t454 = sin(pkin(6));
t457 = sin(qJ(2));
t586 = t454 * t457;
t605 = cos(pkin(6));
t375 = t456 * t605 + t459 * t586;
t460 = cos(qJ(2));
t585 = t454 * t460;
t287 = -t375 * t455 - t458 * t585;
t288 = t375 * t458 - t455 * t585;
t136 = -t604 * t287 + t288 * t603;
t357 = t481 * t456;
t358 = t482 * t456;
t581 = t456 * t458;
t565 = pkin(4) * t581;
t197 = pkin(5) * t357 + qJ(6) * t358 + t565;
t374 = t456 * t586 - t459 * t605;
t583 = t455 * t456;
t401 = mrSges(5,2) * t459 - mrSges(5,3) * t583;
t483 = t287 * t603 + t288 * t604;
t582 = t455 * t459;
t564 = pkin(8) * t582;
t289 = t510 - t564;
t139 = t289 * t603 + t758;
t140 = t289 * t604 - t759;
t523 = t136 * t139 + t140 * t483;
t304 = -mrSges(6,1) * t459 - t357 * mrSges(6,3);
t622 = t357 * mrSges(7,2);
t305 = mrSges(7,1) * t459 + t622;
t529 = -t304 / 0.2e1 + t305 / 0.2e1;
t443 = t459 * mrSges(7,3);
t621 = t358 * mrSges(7,2);
t301 = -t443 - t621;
t302 = mrSges(6,2) * t459 - t358 * mrSges(6,3);
t530 = t301 / 0.2e1 + t302 / 0.2e1;
t594 = t288 * t458;
t595 = t287 * t455;
t631 = mrSges(5,3) * t456;
t403 = -mrSges(5,1) * t459 - mrSges(5,3) * t581;
t650 = -t403 / 0.2e1;
t666 = -t358 / 0.2e1;
t667 = -t357 / 0.2e1;
t676 = t287 / 0.2e1;
t702 = m(7) / 0.2e1;
t704 = m(6) / 0.2e1;
t743 = mrSges(6,3) + mrSges(7,2);
t753 = t743 * (t136 * t666 + t483 * t667) + (-t126 * t483 - t127 * t136 + t374 * t565 + t523) * t704 + (-t111 * t136 + t117 * t483 + t197 * t374 + t523) * t702 + t401 * t676 + t288 * t650 + (t595 / 0.2e1 - t594 / 0.2e1) * t631 + t529 * t483 - t530 * t136;
t577 = t459 * t460;
t333 = (-t455 * t577 + t457 * t458) * t454;
t334 = (t455 * t457 + t458 * t577) * t454;
t183 = -t333 * t604 + t334 * t603;
t184 = t333 * t603 + t334 * t604;
t448 = t458 * pkin(9);
t568 = t448 + t442;
t293 = -t636 * t528 + t568 * t603;
t716 = t636 * t526 + t568 * t604;
t499 = t293 * t183 + t184 * t716;
t555 = t456 * t585;
t592 = t334 * t458;
t593 = t333 * t455;
t703 = -m(7) / 0.2e1;
t705 = -m(6) / 0.2e1;
t752 = -t743 * (t183 * t654 + t184 * t657) + (t593 / 0.2e1 - t592 / 0.2e1) * mrSges(5,3) - m(5) * (-pkin(3) * t555 + (t592 - t593) * pkin(9)) / 0.2e1 + (t437 * t555 + t499) * t705 + (t241 * t555 + t499) * t703;
t751 = 0.2e1 * t702;
t447 = t456 * pkin(8);
t405 = pkin(4) * t583 + t447;
t171 = pkin(5) * t358 - qJ(6) * t357 + t405;
t220 = t357 * mrSges(6,1) - t358 * mrSges(6,2);
t221 = mrSges(7,1) * t358 - mrSges(7,3) * t357;
t638 = pkin(4) * t455;
t247 = pkin(5) * t482 - qJ(6) * t481 + t638;
t724 = -mrSges(5,1) * t458 + mrSges(5,2) * t455;
t385 = t724 * t456;
t522 = t139 * t293 + t140 * t716;
t640 = -t459 / 0.4e1;
t646 = t437 / 0.2e1;
t275 = mrSges(6,1) * t482 + mrSges(6,2) * t481;
t679 = t275 / 0.2e1;
t274 = mrSges(7,1) * t482 - mrSges(7,3) * t481;
t680 = t274 / 0.2e1;
t681 = t247 / 0.2e1;
t219 = t357 * mrSges(7,1) + t358 * mrSges(7,3);
t684 = t219 / 0.2e1;
t687 = t197 / 0.2e1;
t728 = t481 * t742 + t482 * t754;
t749 = (-t111 * t293 + t117 * t716 + t171 * t247 + t197 * t241 + t522) * t702 + pkin(3) * t385 / 0.2e1 + t171 * t680 + t276 * t687 + t241 * t684 + t221 * t681 + t405 * t679 + t220 * t646 + t728 * t640 + t529 * t716 - t530 * t293;
t745 = m(6) + m(7);
t744 = mrSges(6,1) + mrSges(7,1);
t444 = Ifges(5,5) * t458;
t623 = Ifges(5,6) * t455;
t741 = Ifges(4,4) - t444 / 0.2e1 + t623 / 0.2e1;
t736 = t219 + t220;
t735 = t274 + t275;
t277 = -mrSges(6,1) * t481 + mrSges(6,2) * t482;
t734 = t276 + t277;
t625 = Ifges(7,5) * t481;
t283 = Ifges(7,1) * t482 - t625;
t384 = Ifges(6,4) * t481;
t285 = Ifges(6,1) * t482 + t384;
t733 = t285 + t283;
t445 = Ifges(5,4) * t458;
t732 = -Ifges(5,2) * t455 + t445;
t417 = Ifges(5,1) * t455 + t445;
t727 = t357 * t754 - t358 * t742;
t421 = pkin(3) * t456 - pkin(9) * t459;
t347 = pkin(8) * t583 + t458 * t421;
t348 = -pkin(8) * t581 + t455 * t421;
t726 = -t347 * t455 + t348 * t458;
t609 = t459 * mrSges(4,2);
t414 = t456 * mrSges(4,1) + t609;
t389 = t417 * t456;
t651 = -t401 / 0.2e1;
t723 = pkin(9) * t651 - t389 / 0.4e1;
t629 = Ifges(5,4) * t455;
t415 = Ifges(5,2) * t458 + t629;
t388 = t456 * t415;
t722 = pkin(9) * t650 - t388 / 0.4e1;
t697 = mrSges(7,3) / 0.2e1;
t698 = -mrSges(6,2) / 0.2e1;
t721 = t697 + t698;
t720 = t704 + t702;
t719 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t551 = t603 * pkin(4);
t430 = t551 + qJ(6);
t398 = m(7) * t430 + mrSges(7,3);
t552 = t604 * pkin(4);
t436 = -t552 - pkin(5);
t701 = m(6) * pkin(4);
t566 = t701 / 0.2e1;
t712 = t436 * t702 - t604 * t566;
t711 = m(7) * t436 - t604 * t701 - t744;
t710 = t603 * t701 - mrSges(6,2) + t398;
t709 = t430 * t702 + t566 * t603 + t721;
t707 = 2 * qJD(3);
t706 = m(5) / 0.2e1;
t700 = mrSges(5,1) / 0.2e1;
t699 = -mrSges(5,2) / 0.2e1;
t696 = t111 / 0.2e1;
t695 = -t117 / 0.2e1;
t271 = pkin(4) * t456 - qJ(5) * t579 + t347;
t291 = -qJ(5) * t582 + t348;
t132 = t271 * t604 - t291 * t603;
t120 = -t456 * pkin(5) - t132;
t694 = -t120 / 0.2e1;
t693 = t126 / 0.2e1;
t692 = t127 / 0.2e1;
t691 = t483 / 0.2e1;
t690 = -t139 / 0.2e1;
t689 = -t140 / 0.2e1;
t628 = Ifges(6,4) * t357;
t206 = -Ifges(6,2) * t358 - t459 * Ifges(6,6) + t628;
t686 = t206 / 0.4e1;
t212 = t482 * t374;
t685 = t212 / 0.2e1;
t222 = mrSges(6,1) * t358 + mrSges(6,2) * t357;
t683 = t222 / 0.2e1;
t626 = Ifges(7,5) * t358;
t225 = Ifges(7,3) * t357 - t626;
t682 = t225 / 0.4e1;
t278 = Ifges(7,3) * t482 + t625;
t678 = t278 / 0.4e1;
t627 = Ifges(6,4) * t482;
t281 = Ifges(6,2) * t481 + t627;
t677 = t281 / 0.4e1;
t359 = t482 * t459;
t300 = -mrSges(7,2) * t359 + mrSges(7,3) * t456;
t675 = t300 / 0.2e1;
t303 = -mrSges(6,2) * t456 - mrSges(6,3) * t359;
t672 = t303 / 0.2e1;
t610 = t456 * mrSges(7,1);
t360 = t481 * t459;
t617 = t360 * mrSges(7,2);
t307 = -t610 + t617;
t669 = t307 / 0.2e1;
t664 = -t359 / 0.2e1;
t663 = t359 / 0.2e1;
t662 = t360 / 0.2e1;
t661 = t374 / 0.2e1;
t660 = -t385 / 0.2e1;
t655 = -t481 / 0.2e1;
t652 = -t482 / 0.2e1;
t632 = mrSges(5,2) * t458;
t413 = mrSges(5,1) * t455 + t632;
t649 = t413 / 0.2e1;
t648 = -t415 / 0.4e1;
t418 = Ifges(5,1) * t458 - t629;
t647 = t418 / 0.4e1;
t645 = -t455 / 0.2e1;
t644 = t455 / 0.2e1;
t643 = t456 / 0.2e1;
t642 = t458 / 0.2e1;
t639 = m(7) * t139;
t449 = t459 * pkin(8);
t634 = mrSges(4,2) * t456;
t620 = t359 * mrSges(6,1);
t619 = t359 * mrSges(7,1);
t618 = t360 * mrSges(6,2);
t616 = t360 * mrSges(7,3);
t615 = t481 * mrSges(7,2);
t614 = t481 * mrSges(6,3);
t613 = t482 * mrSges(7,2);
t612 = t482 * mrSges(6,3);
t608 = t459 * Ifges(5,5);
t607 = t459 * Ifges(5,6);
t606 = -mrSges(4,1) + t724;
t213 = t481 * t374;
t264 = t374 * t375;
t14 = m(5) * (t264 + (-t594 + t595) * t374) + t745 * (-t136 * t212 - t213 * t483 + t264);
t602 = t14 * qJD(1);
t308 = t374 * t555;
t15 = m(5) * (t287 * t333 + t288 * t334 + t308) + m(4) * (t374 * t456 + t375 * t459 - t586) * t585 + t745 * (t136 * t183 + t184 * t483 + t308);
t601 = t15 * qJD(1);
t589 = t357 * t374;
t588 = t374 * t455;
t587 = t405 * t455;
t352 = t456 * t732 - t607;
t584 = t455 * t352;
t354 = t418 * t456 - t608;
t580 = t458 * t354;
t578 = t459 * t483;
t133 = t603 * t271 + t604 * t291;
t574 = t302 + t301;
t573 = t304 - t305;
t406 = pkin(4) * t582 + t449;
t562 = t638 / 0.2e1;
t561 = pkin(8) * t649;
t558 = -mrSges(6,3) / 0.2e1 - mrSges(7,2) / 0.2e1;
t557 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t556 = -Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t554 = t608 / 0.2e1;
t536 = t374 * t652;
t533 = t459 * t655;
t343 = Ifges(7,5) * t357;
t204 = -t459 * Ifges(7,6) + Ifges(7,3) * t358 + t343;
t532 = t204 / 0.2e1 - t206 / 0.2e1;
t208 = Ifges(7,1) * t357 - t459 * Ifges(7,4) + t626;
t346 = Ifges(6,4) * t358;
t210 = Ifges(6,1) * t357 - t459 * Ifges(6,5) - t346;
t531 = t208 / 0.2e1 + t210 / 0.2e1;
t525 = t444 - t623;
t521 = -t212 * t293 - t213 * t716;
t519 = t293 * t357 - t358 * t716;
t516 = t565 / 0.2e1;
t515 = mrSges(6,3) * t552;
t514 = mrSges(6,3) * t551;
t509 = 0.2e1 * (m(6) / 0.4e1 + m(7) / 0.4e1) * t375;
t119 = qJ(6) * t456 + t133;
t172 = pkin(5) * t359 - qJ(6) * t360 + t406;
t223 = -t616 + t619;
t224 = t618 + t620;
t306 = mrSges(6,1) * t456 - mrSges(6,3) * t360;
t335 = t395 - t564;
t336 = t455 * t498 + t563;
t386 = t413 * t456;
t387 = t413 * t459;
t402 = -mrSges(5,2) * t456 - mrSges(5,3) * t582;
t404 = mrSges(5,1) * t456 - mrSges(5,3) * t579;
t461 = -t530 * t213 + (t672 + t675) * t483 - t529 * t212 + (-t306 / 0.2e1 + t669) * t136 + (t386 / 0.2e1 + t683 + t221 / 0.2e1) * t375 + (t387 / 0.2e1 + t224 / 0.2e1 + t223 / 0.2e1 + t403 * t644 + t458 * t651) * t374 + (t375 * t447 + t287 * t347 + t288 * t348 + (t335 * t455 - t336 * t458 + t449) * t374) * t706 + (t126 * t212 - t127 * t213 - t132 * t136 + t133 * t483 + t374 * t406 + t375 * t405) * t704 + (-t111 * t213 - t117 * t212 + t119 * t483 + t120 * t136 + t171 * t375 + t172 * t374) * t702 + t404 * t676 + t288 * t402 / 0.2e1;
t3 = t461 + (t609 / 0.2e1 - t414 / 0.2e1 + (mrSges(4,1) / 0.2e1 - t724 / 0.2e1 - t277 / 0.2e1 - t276 / 0.2e1) * t456) * t585 + t752;
t205 = Ifges(7,5) * t360 + t456 * Ifges(7,6) + Ifges(7,3) * t359;
t207 = Ifges(6,4) * t360 - Ifges(6,2) * t359 + t456 * Ifges(6,6);
t209 = Ifges(7,1) * t360 + Ifges(7,4) * t456 + Ifges(7,5) * t359;
t211 = Ifges(6,1) * t360 - Ifges(6,4) * t359 + Ifges(6,5) * t456;
t353 = Ifges(5,6) * t456 + t459 * t732;
t355 = Ifges(5,5) * t456 + t418 * t459;
t5 = t531 * t360 + t532 * t359 + t335 * t404 + t405 * t224 + t406 * t222 - pkin(2) * t414 + t348 * t401 + t336 * t402 + t347 * t403 + t133 * t302 + t127 * t303 + t132 * t304 + t120 * t305 + t126 * t306 + t117 * t307 + t111 * t300 + t119 * t301 + t172 * t221 + t171 * t223 + (pkin(8) * t387 + t353 * t645 + t355 * t642 - t741 * t456 - t556 * t358 - t557 * t357 + (m(5) * pkin(8) ^ 2 + Ifges(4,1) - Ifges(4,2) - t719) * t459) * t456 + m(5) * (t335 * t347 + t336 * t348) + (-t207 / 0.2e1 + t205 / 0.2e1) * t358 - (-t209 / 0.2e1 - t211 / 0.2e1) * t357 + m(7) * (t111 * t119 + t117 * t120 + t171 * t172) + m(6) * (t126 * t132 + t127 * t133 + t405 * t406) + (pkin(8) * t386 - t584 / 0.2e1 + t580 / 0.2e1 + t557 * t360 + t556 * t359 + t741 * t459) * t459;
t506 = t3 * qJD(1) + t5 * qJD(2);
t467 = t709 * t184 + t333 * t700 + t334 * t699 + (-t744 / 0.2e1 + t712) * t183;
t6 = t374 * t660 + t736 * t661 - t467 + t753;
t226 = -Ifges(6,2) * t357 - t346;
t227 = -Ifges(7,1) * t358 + t343;
t228 = -Ifges(6,1) * t358 - t628;
t8 = t405 * t220 + t335 * t401 - t336 * t403 + t171 * t219 + t197 * t221 + t574 * t140 - t573 * t139 + m(7) * (t111 * t140 + t117 * t139 + t171 * t197) + m(6) * (-t126 * t139 + t127 * t140) + (-t226 / 0.2e1 + t225 / 0.2e1 + t126 * mrSges(6,3) - t117 * mrSges(7,2) - t531) * t358 - (-t227 / 0.2e1 - t228 / 0.2e1 + t127 * mrSges(6,3) + t111 * mrSges(7,2) - t532) * t357 + (-pkin(8) * t385 + (t554 + t335 * mrSges(5,3) - t354 / 0.2e1 + t388 / 0.2e1) * t455 + (t607 / 0.2e1 - t336 * mrSges(5,3) - t389 / 0.2e1 - t352 / 0.2e1 + (m(6) * t405 + t222) * pkin(4)) * t458) * t456 - t727 * t459 / 0.2e1;
t505 = t6 * qJD(1) + t8 * qJD(2);
t22 = -t574 * t358 - t573 * t357 + m(6) * (-t126 * t357 - t127 * t358) + m(7) * (-t111 * t358 + t117 * t357);
t474 = t720 * (t136 * t357 - t358 * t483);
t492 = t720 * t555;
t28 = t492 - t474;
t504 = -qJD(1) * t28 + qJD(2) * t22;
t36 = t459 * t301 - m(7) * (-t111 * t459 - t171 * t357) + t357 * t221;
t49 = 0.2e1 * (t183 / 0.4e1 + t578 / 0.4e1 + t589 / 0.4e1) * m(7);
t503 = -qJD(1) * t49 - qJD(2) * t36;
t471 = (t357 * t436 - t358 * t430) * t702 + (-t357 * t604 - t358 * t603) * t566;
t487 = m(6) * t516 + m(7) * t687;
t35 = -t471 + t487 + t736;
t470 = (t430 * t481 + t436 * t482) * t702 + (t481 * t603 - t482 * t604) * t566;
t489 = m(6) * t562 + m(7) * t681;
t42 = -t470 + t489 + t735;
t502 = qJD(2) * t35 + qJD(3) * t42;
t200 = m(7) * t357;
t244 = m(7) * t482;
t497 = -qJD(2) * t200 - qJD(3) * t244;
t495 = t172 * t703 + t406 * t705;
t494 = m(7) * t694 + t610 / 0.2e1;
t465 = -t212 * t712 - t213 * t709 + t588 * t700 + t632 * t661 + t685 * t744;
t469 = pkin(4) * t588 * t704 + t247 * t374 * t702;
t10 = t374 * t649 + t735 * t661 - t465 + t469 + t743 * (t652 + t654) * t483;
t381 = Ifges(7,5) * t482;
t279 = -Ifges(7,3) * t481 + t381;
t280 = -Ifges(6,2) * t482 + t384;
t282 = Ifges(7,1) * t481 + t381;
t284 = Ifges(6,1) * t481 - t627;
t11 = -pkin(3) * t413 + t241 * t274 + (t732 / 0.2e1 + t417 / 0.2e1) * t458 + (-t415 / 0.2e1 + t418 / 0.2e1 + pkin(4) * t277) * t455 - (-t279 / 0.2e1 - t282 / 0.2e1 - t284 / 0.2e1 + t281 / 0.2e1) * t482 - (-t285 / 0.2e1 - t280 / 0.2e1 - t283 / 0.2e1 + t278 / 0.2e1) * t481 + (m(6) * t638 + t275) * t437 + t757 * t247;
t462 = (t119 * t430 + t120 * t436) * t702 + Ifges(6,6) * t664 + Ifges(7,6) * t663 + t119 * t697 + mrSges(7,1) * t694 + t132 * mrSges(6,1) / 0.2e1 + t133 * t698 + t347 * t700 + t348 * t699 + t430 * t675 + t436 * t669 + (t132 * t604 + t133 * t603) * t566 + t458 * t554 - Ifges(5,6) * t582 / 0.2e1 + t306 * t552 / 0.2e1 + t551 * t672 + t742 * t662 + t719 * t643;
t480 = -t126 * t716 - t127 * t293 + t522;
t2 = -t584 / 0.4e1 + t580 / 0.4e1 + (t284 + t282 + t279) * t357 / 0.4e1 - t631 * t755 / 0.2e1 - t462 + t277 * t516 + t749 + (t226 + t210 + t208) * t481 / 0.4e1 - t481 * t682 + (t228 + t227 + t204) * t482 / 0.4e1 - t482 * t686 - (t417 + t732) * t583 / 0.4e1 + ((t437 * t581 + t587) * pkin(4) + t480) * t704 - t614 * t693 - t615 * t695 - t613 * t696 - t612 * t692 - t357 * t677 + t358 * t678 + t581 * t647 + t581 * t648 + t456 * t561 + t222 * t562 + t525 * t640 + t722 * t458 + t723 * t455 - (t280 + t733) * t358 / 0.4e1 + t743 * (t139 * t654 + t140 * t657 + t293 * t666 + t667 * t716);
t486 = t10 * qJD(1) + t2 * qJD(2) + t11 * qJD(3);
t464 = -(-t358 * t558 - t530) * t481 - (t357 * t558 - t529) * t482 + (-t126 * t482 + t127 * t481 + t519) * t704 + (t111 * t481 + t117 * t482 + t519) * t702;
t13 = t721 * t360 + (-mrSges(7,1) / 0.2e1 - mrSges(6,1) / 0.2e1) * t359 + t464 + t495;
t473 = (t705 + t703) * (t136 * t482 + t481 * t483);
t30 = t509 + t473;
t34 = (t481 ^ 2 + t482 ^ 2) * t743 + t745 * (t293 * t482 + t481 * t716);
t485 = -qJD(1) * t30 + qJD(2) * t13 + qJD(3) * t34;
t100 = (t536 + t685) * m(7);
t472 = (-t171 * t482 - t241 * t357 - t459 * t716) * t702 + t276 * t667 + t221 * t652;
t27 = (t533 - t360 / 0.2e1) * mrSges(7,2) + t472 + t494;
t67 = t757 * t482;
t484 = qJD(1) * t100 + qJD(2) * t27 - qJD(3) * t67;
t479 = -t443 + ((-qJ(6) - t430) * t459 + t127) * t702;
t41 = -t639 / 0.2e1 + t479;
t478 = -qJD(2) * t41 - qJD(4) * t398;
t453 = t459 ^ 2;
t451 = t456 ^ 2;
t411 = t451 * pkin(8) * t585;
t99 = m(7) * t536 - t212 * t702;
t68 = t716 * t751 + t615;
t58 = t470 + t489;
t56 = m(7) * t691 + t483 * t702;
t50 = (-t578 - t589 + t183) * t702;
t46 = t471 + t487;
t40 = -t621 + t639 / 0.2e1 + t479;
t31 = t509 - t473;
t29 = t492 + t474;
t26 = mrSges(7,2) * t533 + t617 / 0.2e1 + t472 - t494;
t12 = -t616 / 0.2e1 + t619 / 0.2e1 + t618 / 0.2e1 + t620 / 0.2e1 + t464 - t495;
t9 = (t649 + t679 + t680) * t374 + t465 + t469 - t743 * (-t483 / 0.2e1 + t691) * t482;
t7 = t467 + (t660 + t220 / 0.2e1 + t684) * t374 + t753;
t4 = t461 - t414 * t585 + (t724 + t734) * t555 / 0.2e1 - t752;
t1 = t462 + t444 * t640 + (t561 + (-t417 / 0.4e1 - t732 / 0.4e1) * t455 + (-t452 / 0.2e1 - t450 / 0.2e1) * pkin(9) * mrSges(5,3) + (t647 + t648 + (t277 / 0.2e1 + m(6) * t646) * pkin(4)) * t458) * t456 + (t354 / 0.4e1 + t722) * t458 + (pkin(4) * t587 + t480) * t704 + (t607 / 0.4e1 - t352 / 0.4e1 + pkin(4) * t683 + t723) * t455 + (-t285 / 0.4e1 - t283 / 0.4e1 - t280 / 0.4e1 + t678 + t558 * t293) * t358 - (-t284 / 0.4e1 - t282 / 0.4e1 - t279 / 0.4e1 + t677 - t558 * t716) * t357 - (t686 - t228 / 0.4e1 - t227 / 0.4e1 - t204 / 0.4e1 + (t692 + t690) * mrSges(6,3) + (t690 + t696) * mrSges(7,2)) * t482 - (-t226 / 0.4e1 - t210 / 0.4e1 - t208 / 0.4e1 + t682 + (t689 + t693) * mrSges(6,3) + (t689 + t695) * mrSges(7,2)) * t481 + t749;
t16 = [qJD(2) * t15 + qJD(3) * t14, t4 * qJD(3) + t7 * qJD(4) + t29 * qJD(5) + t50 * qJD(6) + t601 + (-t183 * t304 + t183 * t305 + t184 * t301 + t184 * t302 + t333 * t403 + t334 * t401 + ((-mrSges(4,1) * t459 - mrSges(3,1) + t634) * t457 + (-mrSges(3,2) + (t451 + t453) * mrSges(4,3) + (t221 + t222 + t386) * t456) * t460) * t454 + 0.2e1 * (-t126 * t183 + t127 * t184 + t405 * t555) * t704 + (t111 * t184 + t117 * t183 + t171 * t555) * t751 + 0.2e1 * (t333 * t335 + t334 * t336 + t411) * t706 + m(4) * (t411 + (pkin(8) * t453 * t460 - pkin(2) * t457) * t454)) * qJD(2), t602 + t4 * qJD(2) + t9 * qJD(4) + t31 * qJD(5) + t99 * qJD(6) + ((t375 * t437 + t521) * t704 + (t241 * t375 + t521) * t702 + (-pkin(3) * t375 - t374 * t755) * t706) * t707 + ((-mrSges(5,3) * t731 + mrSges(4,2)) * t374 + (t606 + t734) * t375 + t743 * (-t212 * t482 - t213 * t481)) * qJD(3), t7 * qJD(2) + t9 * qJD(3) + (-t288 * mrSges(5,1) - t287 * mrSges(5,2) - t136 * t710 + t483 * t711) * qJD(4) + t56 * qJD(6), qJD(2) * t29 + qJD(3) * t31, qJD(2) * t50 + qJD(3) * t99 + qJD(4) * t56; qJD(3) * t3 + qJD(4) * t6 - qJD(5) * t28 - qJD(6) * t49 - t601, qJD(3) * t5 + qJD(4) * t8 + qJD(5) * t22 - qJD(6) * t36, t1 * qJD(4) + t12 * qJD(5) + t26 * qJD(6) + ((t119 * t716 + t120 * t293 + t172 * t241) * t702 + (-t132 * t293 + t133 * t716 + t406 * t437) * t704 + (-pkin(3) * t449 + pkin(9) * t726) * t706) * t707 + t506 + ((Ifges(5,5) * t455 + Ifges(5,6) * t458 - t481 * t754 + t742 * t482) * t643 + t120 * t613 + t133 * t614 + t119 * t615 + (-t306 + t307) * t293 - Ifges(4,6) * t456 + t437 * t224 + t406 * t277 - pkin(3) * t387 + t172 * t276 + t241 * t223 + (t303 + t300) * t716 + t279 * t663 + t281 * t664 + t205 * t655 + t207 * t657 + (pkin(8) * t606 + t415 * t645 + t417 * t642 + Ifges(4,5)) * t459 - t455 * pkin(9) * t404 - t132 * t612 + pkin(8) * t634 + t402 * t448 + t353 * t642 + t355 * t644 + (t211 + t209) * t654 + t726 * mrSges(5,3) + t733 * t662) * qJD(3), t1 * qJD(3) + (-t336 * mrSges(5,1) - t335 * mrSges(5,2) - Ifges(5,5) * t583 - Ifges(5,6) * t581 + t139 * t711 + t140 * t710 - t357 * t514 + t358 * t515 - t430 * t622 - t436 * t621 + t727) * qJD(4) + t46 * qJD(5) + t40 * qJD(6) + t505, qJD(3) * t12 + qJD(4) * t46 + t504, qJD(3) * t26 + qJD(4) * t40 + t503; -qJD(2) * t3 + qJD(4) * t10 - qJD(5) * t30 + qJD(6) * t100 - t602, qJD(4) * t2 + qJD(5) * t13 + qJD(6) * t27 - t506, qJD(4) * t11 + qJD(5) * t34 - qJD(6) * t67 (pkin(9) * t724 - t293 * t710 - t430 * t613 + t436 * t615 - t481 * t515 - t482 * t514 + t711 * t716 + t525 + t728) * qJD(4) + t58 * qJD(5) + t68 * qJD(6) + t486, qJD(4) * t58 + t485, qJD(4) * t68 + t484; -qJD(2) * t6 - qJD(3) * t10, -qJD(3) * t2 - qJD(5) * t35 + qJD(6) * t41 - t505, -qJD(5) * t42 - t486, t398 * qJD(6), -t502, -t478; qJD(2) * t28 + qJD(3) * t30, -qJD(3) * t13 + qJD(4) * t35 - qJD(6) * t200 - t504, qJD(4) * t42 - qJD(6) * t244 - t485, t502, 0, t497; qJD(2) * t49 - qJD(3) * t100, -qJD(3) * t27 - qJD(4) * t41 + qJD(5) * t200 - t503, qJD(5) * t244 - t484, t478, -t497, 0;];
Cq  = t16;
