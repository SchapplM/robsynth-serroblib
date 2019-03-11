% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRP11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:31
% EndTime: 2019-03-09 12:45:54
% DurationCPUTime: 13.90s
% Computational Cost: add. (21258->720), mult. (40424->925), div. (0->0), fcn. (39661->6), ass. (0->379)
t445 = sin(qJ(2));
t447 = cos(qJ(2));
t443 = sin(qJ(5));
t444 = sin(qJ(4));
t446 = cos(qJ(4));
t678 = cos(qJ(5));
t483 = t443 * t446 + t444 * t678;
t355 = t483 * t447;
t330 = t355 * qJ(6);
t448 = -pkin(2) - pkin(8);
t609 = qJ(3) * t445;
t375 = t447 * t448 - pkin(1) - t609;
t705 = pkin(3) + pkin(7);
t416 = t705 * t445;
t265 = -t375 * t444 + t446 * t416;
t579 = t444 * t447;
t235 = pkin(9) * t579 + t265;
t205 = pkin(4) * t445 + t235;
t266 = t375 * t446 + t416 * t444;
t572 = t446 * t447;
t236 = -pkin(9) * t572 + t266;
t209 = t443 * t236;
t96 = t678 * t205 - t209;
t76 = t330 + t96;
t70 = pkin(5) * t445 + t76;
t104 = t678 * t235 - t209;
t87 = t330 + t104;
t792 = t70 - t87;
t535 = t678 * t446;
t397 = t443 * t444 - t535;
t772 = Ifges(6,6) + Ifges(7,6);
t773 = Ifges(6,5) + Ifges(7,5);
t512 = t397 * t772 - t483 * t773;
t578 = t444 * t448;
t409 = -t444 * pkin(9) + t578;
t410 = (-pkin(9) + t448) * t446;
t732 = -t443 * t409 + t678 * t410;
t751 = t732 * mrSges(6,2);
t749 = t397 * qJ(6) + t732;
t769 = t749 * mrSges(7,2);
t286 = t678 * t409 + t443 * t410;
t770 = t286 * mrSges(6,1);
t201 = -t483 * qJ(6) + t286;
t787 = t201 * mrSges(7,1);
t791 = t512 - t751 - t769 - t770 - t787;
t790 = -t751 / 0.2e1 - t769 / 0.2e1 - t770 / 0.2e1 - t787 / 0.2e1;
t789 = Ifges(6,1) + Ifges(7,1);
t788 = Ifges(6,2) + Ifges(7,2);
t379 = t483 * mrSges(7,1);
t625 = t397 * mrSges(7,2);
t276 = t379 - t625;
t431 = t444 * pkin(4) + qJ(3);
t671 = t483 * pkin(5);
t329 = t431 + t671;
t782 = m(7) * t329;
t786 = t276 + t782;
t563 = t678 * pkin(4);
t433 = t563 + pkin(5);
t510 = t563 - t433;
t630 = t355 * mrSges(7,3);
t301 = mrSges(7,1) * t445 + t630;
t631 = t355 * mrSges(6,3);
t302 = mrSges(6,1) * t445 + t631;
t354 = t397 * t447;
t632 = t354 * mrSges(7,3);
t297 = -mrSges(7,2) * t445 + t632;
t702 = t297 / 0.2e1;
t783 = t286 / 0.2e1;
t785 = t201 * t301 / 0.2e1 + t302 * t783 - t749 * t702;
t211 = t678 * t236;
t103 = -t235 * t443 - t211;
t218 = -t354 * mrSges(7,1) - t355 * mrSges(7,2);
t417 = t705 * t447;
t370 = pkin(4) * t572 + t417;
t243 = -pkin(5) * t354 + t370;
t673 = pkin(5) * t355;
t293 = -pkin(4) * t579 - t673;
t672 = pkin(5) * t397;
t674 = pkin(4) * t446;
t351 = -t672 + t674;
t615 = t446 * mrSges(5,2);
t621 = t444 * mrSges(5,1);
t506 = -t615 - t621;
t372 = t506 * t447;
t622 = t483 * mrSges(7,3);
t539 = t622 / 0.2e1;
t541 = t630 / 0.2e1;
t543 = -t632 / 0.2e1;
t628 = t355 * Ifges(7,4);
t183 = t354 * Ifges(7,2) + t445 * Ifges(7,6) - t628;
t629 = t355 * Ifges(6,4);
t185 = t354 * Ifges(6,2) + t445 * Ifges(6,6) - t629;
t336 = Ifges(7,4) * t354;
t187 = -t355 * Ifges(7,1) + t445 * Ifges(7,5) + t336;
t337 = Ifges(6,4) * t354;
t189 = -t355 * Ifges(6,1) + t445 * Ifges(6,5) + t337;
t331 = t354 * mrSges(7,2);
t214 = -t355 * mrSges(7,1) + t331;
t215 = -mrSges(6,1) * t355 + mrSges(6,2) * t354;
t220 = Ifges(7,2) * t355 + t336;
t221 = Ifges(6,2) * t355 + t337;
t222 = Ifges(7,1) * t354 + t628;
t223 = Ifges(6,1) * t354 + t629;
t505 = -mrSges(6,1) * t397 - mrSges(6,2) * t483;
t378 = t483 * mrSges(7,2);
t523 = -t397 * mrSges(7,1) - t378;
t633 = t354 * mrSges(6,3);
t655 = Ifges(7,4) * t397;
t280 = -Ifges(7,2) * t483 - t655;
t656 = Ifges(6,4) * t397;
t281 = -Ifges(6,2) * t483 - t656;
t739 = t280 + t281;
t753 = -t732 / 0.2e1;
t755 = t445 / 0.4e1;
t767 = -t483 * t789 + t655 + t656;
t385 = Ifges(7,4) * t483;
t386 = Ifges(6,4) * t483;
t282 = -Ifges(7,1) * t397 - t385;
t283 = -Ifges(6,1) * t397 - t386;
t738 = t282 + t283;
t777 = t397 * t788 - t385 - t386 + t738;
t719 = t370 * t505 / 0.2e1 + t243 * t523 / 0.2e1 + t631 * t783 + t633 * t753 + t431 * t215 / 0.2e1 + t329 * t214 / 0.2e1 - (t223 + t222) * t397 / 0.4e1 + (t185 + t183) * t397 / 0.4e1 + t512 * t755 + t777 * t354 / 0.4e1 - (t221 + t220 + t189 + t187) * t483 / 0.4e1 + (-t767 / 0.4e1 + t739 / 0.4e1) * t355;
t624 = t397 * mrSges(6,3);
t659 = mrSges(7,3) * t397;
t605 = t354 * qJ(6);
t97 = t443 * t205 + t211;
t77 = t97 + t605;
t726 = -t97 * t624 / 0.2e1 - t77 * t659 / 0.2e1;
t456 = t201 * t541 + t70 * t539 + t543 * t749 + t719 - t726;
t689 = -t483 / 0.2e1;
t691 = t397 / 0.2e1;
t696 = t351 / 0.2e1;
t298 = -mrSges(6,2) * t445 + t633;
t701 = t298 / 0.2e1;
t703 = t293 / 0.2e1;
t710 = m(7) / 0.2e1;
t86 = t103 - t605;
t781 = t77 + t86;
t784 = (-t201 * t792 + t243 * t351 + t293 * t329 + t749 * t781) * t710 + qJ(3) * t372 / 0.2e1 + t218 * t696 + t732 * t701 + t276 * t703 + t456 + (t87 * t689 + t86 * t691) * mrSges(7,3) - t785 + (t103 * t691 - t483 * (t104 / 0.2e1 - t96 / 0.2e1)) * mrSges(6,3);
t743 = mrSges(6,3) + mrSges(7,3);
t780 = t104 - t96;
t677 = m(6) * t431;
t774 = mrSges(6,2) + mrSges(7,2);
t771 = Ifges(7,3) + Ifges(6,3);
t768 = t103 + t97;
t736 = t297 + t298;
t735 = t301 + t302;
t762 = mrSges(6,1) / 0.2e1;
t761 = mrSges(7,1) / 0.2e1;
t580 = t444 * t445;
t353 = -t443 * t580 + t445 * t535;
t747 = t353 / 0.2e1;
t356 = t483 * t445;
t745 = t356 / 0.2e1;
t680 = t447 / 0.2e1;
t752 = Ifges(6,4) + Ifges(7,4);
t520 = t445 * pkin(2) - qJ(3) * t447;
t396 = pkin(8) * t445 + t520;
t401 = t446 * t417;
t208 = pkin(4) * t447 + t401 + (-pkin(9) * t445 - t396) * t444;
t275 = t446 * t396 + t444 * t417;
t575 = t445 * t446;
t238 = pkin(9) * t575 + t275;
t100 = t678 * t208 - t238 * t443;
t101 = t443 * t208 + t678 * t238;
t274 = -t396 * t444 + t401;
t300 = mrSges(6,1) * t447 - t356 * mrSges(6,3);
t73 = pkin(5) * t447 - qJ(6) * t356 + t100;
t82 = qJ(6) * t353 + t101;
t469 = t73 * t761 - t82 * mrSges(7,2) / 0.2e1 + t100 * t762 - t101 * mrSges(6,2) / 0.2e1 + t772 * t747 + t773 * t745 + t771 * t680;
t516 = t563 / 0.2e1;
t675 = pkin(4) * t443;
t555 = t675 / 0.2e1;
t299 = mrSges(7,1) * t447 - t356 * mrSges(7,3);
t700 = t299 / 0.2e1;
t707 = -mrSges(5,2) / 0.2e1;
t709 = m(6) * pkin(4);
t295 = -mrSges(7,2) * t447 + t353 * mrSges(7,3);
t296 = -mrSges(6,2) * t447 + t353 * mrSges(6,3);
t737 = t295 + t296;
t750 = (t433 * t73 + t675 * t82) * t710 + Ifges(5,3) * t680 + t274 * mrSges(5,1) / 0.2e1 + t275 * t707 + t433 * t700 + (t100 * t678 + t101 * t443) * t709 / 0.2e1 + Ifges(5,5) * t580 / 0.2e1 + Ifges(5,6) * t575 / 0.2e1 + t300 * t516 + t737 * t555 + t469;
t744 = -t447 / 0.2e1;
t660 = t70 - t76;
t547 = m(7) * t660;
t650 = Ifges(5,6) * t446;
t654 = Ifges(5,5) * t444;
t742 = Ifges(3,4) - t654 / 0.2e1 - t650 / 0.2e1 + Ifges(4,6);
t741 = t774 * t397;
t740 = t201 * t659 + t286 * t624;
t366 = t397 * t675;
t734 = t483 * t563 + t366;
t733 = t433 * t483 + t366;
t731 = t774 * t678;
t513 = t354 * t773 + t355 * t772;
t729 = t444 ^ 2 + t446 ^ 2;
t716 = t483 ^ 2;
t717 = t397 ^ 2;
t728 = -t716 - t717;
t219 = -t354 * mrSges(6,1) - t355 * mrSges(6,2);
t658 = Ifges(5,4) * t444;
t415 = t446 * Ifges(5,1) - t658;
t374 = t447 * t415;
t560 = mrSges(5,3) * t572;
t618 = t445 * mrSges(5,2);
t405 = -t560 - t618;
t679 = t448 / 0.2e1;
t727 = t405 * t679 + pkin(4) * t219 / 0.2e1 - t374 / 0.4e1;
t725 = t743 * t355 * t675;
t724 = (-t201 * t710 + t539) * pkin(5) + t790;
t642 = t104 * mrSges(6,2);
t643 = t103 * mrSges(6,1);
t665 = t87 * mrSges(7,2);
t666 = t86 * mrSges(7,1);
t723 = t643 / 0.2e1 - t642 / 0.2e1 + t666 / 0.2e1 - t665 / 0.2e1 + (t86 * t710 + t543) * pkin(5);
t722 = t469 + (t73 * t710 + t700) * pkin(5);
t690 = t483 / 0.2e1;
t692 = -t397 / 0.2e1;
t720 = t735 * t689 + t736 * t692 + t743 * (t354 * t691 + t355 * t690);
t718 = 0.2e1 * m(7);
t713 = -m(6) / 0.2e1;
t712 = m(6) / 0.2e1;
t711 = -m(7) / 0.2e1;
t708 = m(7) * pkin(5);
t501 = Ifges(5,2) * t446 + t658;
t697 = Ifges(5,6) * t744 - t445 * t501 / 0.2e1;
t657 = Ifges(5,4) * t446;
t414 = -t444 * Ifges(5,2) + t657;
t688 = t414 / 0.4e1;
t687 = -t415 / 0.4e1;
t686 = t417 / 0.2e1;
t685 = -t444 / 0.2e1;
t683 = t444 / 0.2e1;
t682 = t445 / 0.2e1;
t681 = t446 / 0.2e1;
t676 = m(7) * t483;
t669 = t76 * mrSges(7,2);
t668 = t77 * mrSges(7,1);
t664 = t96 * mrSges(6,2);
t663 = t97 * mrSges(6,1);
t651 = Ifges(5,6) * t445;
t634 = t353 * mrSges(7,1);
t627 = t356 * mrSges(7,2);
t626 = t397 * mrSges(6,2);
t620 = t444 * mrSges(5,2);
t619 = t445 * mrSges(5,1);
t617 = t445 * Ifges(5,5);
t616 = t446 * mrSges(5,1);
t614 = t447 * mrSges(5,1);
t613 = t447 * mrSges(5,2);
t216 = t627 - t634;
t217 = -mrSges(6,1) * t353 + mrSges(6,2) * t356;
t369 = (-t674 - t705) * t445;
t242 = -pkin(5) * t353 + t369;
t504 = Ifges(5,1) * t444 + t657;
t349 = Ifges(5,5) * t447 + t445 * t504;
t507 = t616 - t620;
t371 = t507 * t445;
t402 = -mrSges(5,3) * t580 + t614;
t403 = mrSges(5,3) * t579 + t619;
t404 = mrSges(5,3) * t575 - t613;
t499 = -pkin(2) * t447 - t609;
t411 = -pkin(1) + t499;
t412 = t447 * mrSges(4,2) - t445 * mrSges(4,3);
t524 = m(4) * t411 + t412;
t525 = t187 / 0.2e1 + t189 / 0.2e1;
t526 = -t752 * t353 / 0.2e1 - t789 * t356 / 0.2e1 + t773 * t744;
t527 = t183 / 0.2e1 + t185 / 0.2e1;
t528 = t772 * t680 + t752 * t745 + t747 * t788;
t549 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t550 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t478 = t501 * t447;
t348 = -t478 + t651;
t574 = t446 * t348;
t479 = t504 * t447;
t350 = -t479 + t617;
t582 = t444 * t350;
t5 = m(5) * (t265 * t274 + t266 * t275 - t416 * t417) + m(6) * (t100 * t96 + t101 * t97 + t369 * t370) + m(7) * (t242 * t243 + t70 * t73 + t77 * t82) + (t574 / 0.2e1 + t582 / 0.2e1 - pkin(1) * mrSges(3,1) - t411 * mrSges(4,2) - t742 * t445 + t550 * t356 + t549 * t353) * t445 - t417 * t371 + t265 * t402 + t274 * t403 + t266 * t404 + t275 * t405 + t369 * t219 + t370 * t217 + t97 * t296 + t82 * t297 + t101 * t298 + t70 * t299 + t96 * t300 + t73 * t301 + t100 * t302 + t77 * t295 + t242 * t218 + t243 * t216 + (-t411 * mrSges(4,3) + t446 * t697 + t349 * t685 - t416 * t507 - pkin(1) * mrSges(3,2) - t550 * t355 + t549 * t354 + (-Ifges(3,2) - Ifges(4,3) + Ifges(3,1) + Ifges(5,3) + Ifges(4,2) + t771) * t445 + t742 * t447) * t447 + t525 * t356 + t526 * t355 + t527 * t353 + t528 * t354 + t524 * t520;
t612 = t5 * qJD(1);
t373 = t447 * t414;
t427 = Ifges(5,6) * t579;
t457 = t243 * t214 + t370 * t215 + (t220 / 0.2e1 + t221 / 0.2e1 + t525) * t354 - t70 * t632 - t96 * t633 + t513 * t682;
t465 = -t223 / 0.2e1 - t222 / 0.2e1 + t77 * mrSges(7,3) + t97 * mrSges(6,3) + t527;
t6 = m(6) * (t103 * t96 + t104 * t97) + m(7) * (t243 * t293 + t70 * t86 + t77 * t87) + t427 * t682 + t465 * t355 + t417 * t372 - t266 * t403 + t265 * t405 + t87 * t297 + t104 * t298 + t86 * t301 + t103 * t302 + t293 * t218 + ((-t617 / 0.2e1 + t373 / 0.2e1 - t350 / 0.2e1 + t265 * mrSges(5,3)) * t446 + (t374 / 0.2e1 + t348 / 0.2e1 + t266 * mrSges(5,3) + (-m(6) * t370 - t219) * pkin(4)) * t444) * t447 + t457;
t611 = t6 * qJD(1);
t7 = t76 * t297 + t96 * t298 - t97 * t302 + (-t301 - t547) * t77 + ((-m(7) * t243 - t218) * pkin(5) + t465) * t355 + t457;
t610 = t7 * qJD(1);
t573 = t446 * t405;
t581 = t444 * t403;
t19 = t735 * t356 - t736 * t353 + m(7) * (-t353 * t77 + t356 * t70) + m(6) * (-t353 * t97 + t356 * t96) + (m(5) * (t265 * t444 - t266 * t446) - t573 + t581 - t524) * t445;
t608 = t19 * qJD(1);
t606 = t353 * t483;
t602 = t354 * t483;
t599 = t355 * t397;
t598 = t356 * t397;
t597 = t370 * t446;
t586 = t433 * t355;
t585 = t433 * t397;
t380 = t483 * mrSges(6,1);
t569 = -t379 - t380;
t559 = (-t599 + t602) * t710;
t554 = t675 / 0.4e1;
t553 = t762 + t761;
t552 = mrSges(6,2) / 0.2e1 + mrSges(7,2) / 0.2e1;
t551 = mrSges(7,3) / 0.2e1 + mrSges(6,3) / 0.2e1;
t548 = mrSges(7,1) + t708;
t542 = t632 / 0.2e1;
t538 = -t622 / 0.2e1;
t536 = mrSges(5,3) * t679;
t531 = -t579 / 0.2e1;
t530 = -t448 * t403 / 0.2e1;
t517 = -t563 / 0.2e1;
t514 = -t329 * t523 - t431 * t505 - t740;
t511 = t552 * t353;
t508 = t563 * t633;
t500 = -t650 - t654;
t310 = t353 * t675;
t458 = t553 * t356 + t511 + (t356 * t563 - t310) * t712 + (t356 * t433 - t310) * t710;
t460 = (-t768 * t397 + t483 * t780) * t713 + (-t781 * t397 - t483 * t792) * t711;
t468 = (-t354 * t551 + t701 + t702) * t397;
t471 = -t301 / 0.2e1 - t302 / 0.2e1 + t551 * t355;
t10 = (t618 / 0.2e1 - t405 / 0.2e1 - t560 / 0.2e1) * t446 + (t619 / 0.2e1 + t403 / 0.2e1 + mrSges(5,3) * t531) * t444 - t471 * t483 + t468 + t458 + t460;
t495 = t10 * qJD(1);
t462 = t511 + (t708 / 0.2e1 + t553) * t356;
t13 = t468 - (-t547 / 0.2e1 + t471) * t483 + t462;
t494 = t13 * qJD(1);
t455 = -t553 * t354 - t552 * t355 + m(5) * t686 + (-t286 * t353 + t356 * t732 + t370) * t712 + (-t201 * t353 + t356 * t749 + t243) * t710;
t490 = t274 * t446 + t275 * t444;
t459 = -m(5) * t490 / 0.2e1 + (-t100 * t397 + t101 * t483) * t713 + (-t397 * t73 + t483 * t82) * t711;
t16 = (-t402 / 0.2e1 + t614 / 0.2e1) * t446 + (-t404 / 0.2e1 - t613 / 0.2e1) * t444 - (t295 / 0.2e1 + t296 / 0.2e1 - t551 * t353) * t483 + (t700 + t300 / 0.2e1 + t551 * t356) * t397 + t455 + t459;
t474 = -t506 - t569;
t74 = mrSges(4,3) - t741 + (m(5) + m(4)) * qJ(3) + t677 + t782 + t474;
t493 = qJD(1) * t16 + qJD(2) * t74;
t461 = (t201 * t354 + t355 * t749 + t397 * t70 - t483 * t77) * t711 + t301 * t692 + t297 * t690;
t470 = t242 * t710 - t634 / 0.2e1 + t627 / 0.2e1;
t23 = (t602 / 0.2e1 - t599 / 0.2e1) * mrSges(7,3) + t461 + t470;
t51 = -m(7) * (-t201 * t483 + t397 * t749) + t728 * mrSges(7,3);
t492 = -qJD(1) * t23 - qJD(2) * t51;
t59 = (t354 * t554 + t586 / 0.4e1 - t293 / 0.4e1) * t718 - t214;
t91 = (-t483 * t554 + t585 / 0.4e1 - t351 / 0.4e1) * t718 - t523;
t491 = qJD(1) * t59 + qJD(2) * t91;
t126 = t355 * t548 - t331;
t179 = t397 * t548 + t378;
t489 = qJD(1) * t126 + qJD(2) * t179;
t487 = t569 + t741;
t484 = t414 * t681 + t415 * t683;
t37 = t355 * t301 + m(7) * (t354 * t77 + t355 * t70) + t354 * t297;
t477 = t37 * qJD(1) + qJD(3) * t559;
t132 = (-t717 / 0.2e1 - t716 / 0.2e1 - 0.1e1 / 0.2e1) * m(7);
t476 = qJD(1) * t559 + t132 * qJD(2);
t277 = t380 - t626;
t467 = t286 * t780 + t768 * t732;
t1 = t729 * t447 * t536 + t444 * t530 - t574 / 0.4e1 + (t478 / 0.4e1 + t727) * t446 - (-t479 - t373) * t444 / 0.4e1 + t507 * t686 + t572 * t687 + t579 * t688 + t500 * t755 + pkin(4) * t277 * t531 + ((-t431 * t579 + t597) * pkin(4) + t467) * t712 - t582 / 0.4e1 - t750 + t784;
t12 = -qJ(3) * t507 + t504 * t681 + t501 * t685 - (-t283 / 0.2e1 - t282 / 0.2e1 + t386 / 0.2e1 + t385 / 0.2e1) * t483 + (-t280 / 0.2e1 - t281 / 0.2e1 + t286 * mrSges(6,3) + t201 * mrSges(7,3) + (Ifges(6,4) / 0.2e1 + Ifges(7,4) / 0.2e1) * t397 - (Ifges(7,1) / 0.2e1 - Ifges(6,2) / 0.2e1 - Ifges(7,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t483) * t397 + t484 + t514 + (-t277 - t677) * t674 - t786 * t351;
t473 = t1 * qJD(1) - t12 * qJD(2);
t17 = t672 * t786 + t690 * t777 + t767 * t691 + t739 * t692 + t514 + t740;
t453 = (-t660 * t201 + (-t243 * t397 - t329 * t355) * pkin(5)) * t711 + t298 * t753 + t276 * t673 / 0.2e1 + t218 * t672 / 0.2e1 + t76 * t539 + t785;
t4 = -t719 + t70 * t538 + t749 * t542 - t201 * t630 / 0.2e1 + t453 + t722;
t472 = -t4 * qJD(1) - t17 * qJD(2);
t198 = t731 * pkin(4) + (-m(7) * t510 + mrSges(6,1) + mrSges(7,1)) * t675;
t454 = t510 * t711 * t201 + t433 * t538 - t517 * t622 - t790;
t20 = t454 + t724;
t451 = (-t433 * t77 + (-t443 * t660 + t678 * t77) * pkin(4)) * t711 + t669 / 0.2e1 + t668 / 0.2e1 + t664 / 0.2e1 + t663 / 0.2e1 + t433 * t542 + t508 / 0.2e1 + t735 * t555 - t725 / 0.2e1 + t736 * t517;
t9 = t451 + t723;
t90 = (pkin(5) / 0.2e1 + t516 - t433 / 0.2e1) * t676;
t466 = t9 * qJD(1) + t20 * qJD(2) - t90 * qJD(3) + t198 * qJD(4);
t131 = (0.1e1 + t728) * t710;
t115 = (-t483 * t675 + t585) * t710 + m(7) * t696;
t113 = qJD(6) * t559;
t102 = (t354 * t675 + t586) * t710 + m(7) * t703;
t52 = t487 + (t510 / 0.2e1 - pkin(5) / 0.2e1) * t676;
t24 = t354 * t538 + t397 * t541 - t461 + t470;
t18 = -t454 + t512 + t724;
t15 = t404 * t683 + t402 * t681 + (m(4) * pkin(7) + t616 / 0.2e1 - t620 / 0.2e1 + mrSges(4,1)) * t447 + t455 - t459 + (t299 + t300) * t692 + t737 * t690 + t743 * (t606 / 0.2e1 + t598 / 0.2e1);
t14 = t547 * t689 + t462 + t720;
t11 = t573 / 0.2e1 - t581 / 0.2e1 + (t615 / 0.2e1 + t621 / 0.2e1) * t445 + t458 - t460 + t729 * mrSges(5,3) * t680 + t720;
t8 = -t451 + t513 + t723;
t3 = -t453 + t456 + t722 + t726;
t2 = (-t651 / 0.4e1 + mrSges(5,1) * t686 - t348 / 0.4e1 + t727) * t446 + (-t617 / 0.4e1 + t417 * t707 + t373 / 0.4e1 - t350 / 0.4e1 + t530) * t444 + (pkin(4) * t597 + t467) * t712 + ((t687 + (Ifges(5,2) / 0.4e1 + t536) * t446) * t446 + (t657 / 0.2e1 + t688 + (Ifges(5,1) / 0.4e1 + t536) * t444 + (-t677 / 0.2e1 - t277 / 0.2e1) * pkin(4)) * t444) * t447 + t750 + t784;
t21 = [qJD(2) * t5 + qJD(3) * t19 + qJD(4) * t6 + qJD(5) * t7 + qJD(6) * t37, t15 * qJD(3) + t2 * qJD(4) + t3 * qJD(5) + t24 * qJD(6) + t612 + (-qJ(3) * t371 + t749 * t299 + t201 * t295 + t329 * t216 + t431 * t217 + t242 * t276 + t369 * t277 + t732 * t300 + t286 * t296 + (-t416 * mrSges(5,2) + t349 / 0.2e1 + t448 * t402 - t274 * mrSges(5,3)) * t446 + (-t416 * mrSges(5,1) - t275 * mrSges(5,3) + t448 * t404 + t697) * t444 - (t101 * mrSges(6,3) + t82 * mrSges(7,3) + t528) * t483 + (t100 * mrSges(6,3) + t73 * mrSges(7,3) + t526) * t397 + (-pkin(2) * mrSges(4,1) + Ifges(5,5) * t681 + Ifges(5,6) * t685 - t550 * t397 - t483 * t549 - Ifges(4,4) + Ifges(3,5)) * t447 + m(5) * (-qJ(3) * t416 + t448 * t490) + 0.2e1 * (t100 * t732 + t101 * t286 + t369 * t431) * t712 + 0.2e1 * (t201 * t82 + t242 * t329 + t73 * t749) * t710 + (-qJ(3) * mrSges(4,1) + Ifges(4,5) - Ifges(3,6) + t484) * t445 + (m(4) * t499 - t447 * mrSges(3,1) + t445 * mrSges(3,2) + t412) * pkin(7) + t739 * t747 + t738 * t745) * qJD(2), t608 + t15 * qJD(2) + 0.2e1 * (t710 + t712) * (-t598 - t606) * qJD(3) + t11 * qJD(4) + t14 * qJD(5) + t113, t611 + t2 * qJD(2) + t11 * qJD(3) + (-t433 * t632 + m(7) * (t433 * t86 + t675 * t87) - t642 + t643 - t665 + t666 + (t103 * t678 + t104 * t443) * t709 - t265 * mrSges(5,2) - t266 * mrSges(5,1) - Ifges(5,5) * t572 + t427 - t508 + t513 + t725) * qJD(4) + t8 * qJD(5) + t102 * qJD(6), t610 + t3 * qJD(2) + t14 * qJD(3) + t8 * qJD(4) + (-t663 - t668 - t664 - t669 + (-m(7) * t77 - t632) * pkin(5) + t513) * qJD(5), t24 * qJD(2) + t102 * qJD(4) + t477; qJD(3) * t16 + qJD(4) * t1 - qJD(5) * t4 - qJD(6) * t23 - t612, qJD(3) * t74 - qJD(4) * t12 - qJD(5) * t17 - qJD(6) * t51, qJD(6) * t131 + t493 (m(7) * (-t433 * t201 + t675 * t749) + (-t286 * t678 + t443 * t732) * t709 - t448 * t615 - mrSges(5,1) * t578 + t500 + t733 * mrSges(7,3) + t734 * mrSges(6,3) + t791) * qJD(4) + t18 * qJD(5) + t115 * qJD(6) + t473, t18 * qJD(4) + ((-m(7) * t201 + t622) * pkin(5) + t791) * qJD(5) + t472, qJD(3) * t131 + qJD(4) * t115 + t492; -qJD(2) * t16 - qJD(4) * t10 - qJD(5) * t13 + t113 - t608, qJD(6) * t132 - t493, 0 (-t474 + t625 + t626) * qJD(4) + t52 * qJD(5) + 0.2e1 * (-t710 * t733 - t712 * t734) * qJD(4) - t495, t52 * qJD(4) + (-m(7) * t671 + t487) * qJD(5) - t494, t476; -qJD(2) * t1 + qJD(3) * t10 - qJD(5) * t9 + qJD(6) * t59 - t611, -qJD(5) * t20 + qJD(6) * t91 - t473, qJD(5) * t90 + t495, -t198 * qJD(5) ((-mrSges(6,1) - t548) * t443 - t731) * qJD(5) * pkin(4) - t466, t491; qJD(2) * t4 + qJD(3) * t13 + qJD(4) * t9 + qJD(6) * t126 - t610, qJD(4) * t20 + qJD(6) * t179 - t472, -qJD(4) * t90 + t494, t466, 0, t489; t23 * qJD(2) - t59 * qJD(4) - t126 * qJD(5) - t477, -qJD(3) * t132 - qJD(4) * t91 - qJD(5) * t179 - t492, -t476, -t491, -t489, 0;];
Cq  = t21;
