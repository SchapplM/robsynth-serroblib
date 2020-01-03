% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR12_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:37:12
% EndTime: 2019-12-31 21:37:48
% DurationCPUTime: 17.90s
% Computational Cost: add. (25308->841), mult. (62337->1231), div. (0->0), fcn. (67368->10), ass. (0->407)
t517 = sin(qJ(3));
t520 = cos(qJ(3));
t514 = sin(pkin(5));
t518 = sin(qJ(2));
t622 = t514 * t518;
t633 = cos(pkin(5));
t462 = t517 * t633 + t520 * t622;
t513 = sin(pkin(10));
t515 = cos(pkin(10));
t521 = cos(qJ(2));
t621 = t514 * t521;
t393 = -t462 * t515 + t513 * t621;
t516 = sin(qJ(5));
t519 = cos(qJ(5));
t551 = -t462 * t513 - t515 * t621;
t224 = t393 * t519 - t516 * t551;
t786 = t224 / 0.2e1;
t787 = mrSges(6,1) * t786;
t783 = t393 * t516 + t519 * t551;
t784 = -t783 / 0.2e1;
t785 = mrSges(6,2) * t784;
t461 = t517 * t622 - t520 * t633;
t655 = t224 * Ifges(6,4);
t90 = Ifges(6,2) * t783 + Ifges(6,6) * t461 - t655;
t782 = t90 / 0.4e1;
t216 = Ifges(6,4) * t783;
t91 = -Ifges(6,1) * t224 + Ifges(6,5) * t461 + t216;
t781 = t91 / 0.4e1;
t780 = -Ifges(5,3) - Ifges(6,3);
t779 = Ifges(3,6) * t622;
t597 = pkin(1) * t633;
t465 = pkin(7) * t621 + t518 * t597;
t442 = pkin(8) * t633 + t465;
t589 = -pkin(2) * t521 - pkin(8) * t518;
t443 = (-pkin(1) + t589) * t514;
t317 = -t442 * t517 + t443 * t520;
t318 = t442 * t520 + t443 * t517;
t776 = t318 * mrSges(4,1) + t317 * mrSges(4,2) + Ifges(4,6) * t462;
t692 = mrSges(5,2) * t515;
t693 = mrSges(5,1) * t513;
t775 = -t692 / 0.2e1 - t693 / 0.2e1;
t494 = pkin(3) * t517 - qJ(4) * t520;
t624 = t513 * t517;
t434 = pkin(8) * t624 + t494 * t515;
t620 = t515 * t517;
t435 = -pkin(8) * t620 + t494 * t513;
t774 = -t434 * t513 + t435 * t515;
t619 = t515 * t519;
t563 = t513 * t516 - t619;
t661 = Ifges(6,6) * t563;
t564 = t513 * t519 + t515 * t516;
t671 = Ifges(6,5) * t564;
t388 = -t661 + t671;
t668 = Ifges(5,6) * t515;
t679 = Ifges(5,5) * t513;
t491 = t668 + t679;
t773 = t491 / 0.2e1 + t388 / 0.2e1;
t463 = -pkin(7) * t622 + t521 * t597;
t464 = (pkin(2) * t518 - pkin(8) * t521) * t514;
t363 = t463 * t520 + t464 * t517;
t330 = qJ(4) * t622 + t363;
t369 = t494 * t621 + t465;
t188 = t330 * t515 + t369 * t513;
t615 = t520 * t521;
t426 = (-t513 * t615 + t515 * t518) * t514;
t427 = (t513 * t518 + t515 * t615) * t514;
t616 = t517 * t521;
t607 = t514 * t616;
t259 = Ifges(5,4) * t427 + Ifges(5,2) * t426 + Ifges(5,6) * t607;
t759 = mrSges(5,3) / 0.2e1;
t772 = t188 * t759 + t259 / 0.4e1;
t771 = t465 * mrSges(3,1) + t463 * mrSges(3,2);
t770 = t514 ^ 2;
t512 = t515 ^ 2;
t769 = 2 * qJD(3);
t768 = -m(5) / 0.2e1;
t767 = m(5) / 0.2e1;
t766 = -m(6) / 0.2e1;
t765 = m(6) / 0.2e1;
t764 = -pkin(2) / 0.2e1;
t763 = -pkin(8) / 0.2e1;
t761 = mrSges(6,2) / 0.2e1;
t760 = -mrSges(5,3) / 0.2e1;
t187 = -t330 * t513 + t369 * t515;
t121 = pkin(4) * t607 - pkin(9) * t427 + t187;
t138 = pkin(9) * t426 + t188;
t64 = t121 * t519 - t138 * t516;
t758 = t64 / 0.2e1;
t757 = -t90 / 0.2e1;
t756 = t91 / 0.2e1;
t755 = pkin(2) * mrSges(4,1);
t754 = pkin(2) * mrSges(4,2);
t753 = pkin(8) * mrSges(4,2);
t294 = t426 * t519 - t427 * t516;
t295 = t426 * t516 + t427 * t519;
t127 = Ifges(6,4) * t295 + Ifges(6,2) * t294 + Ifges(6,6) * t607;
t752 = t127 / 0.4e1;
t128 = Ifges(6,1) * t295 + Ifges(6,4) * t294 + Ifges(6,5) * t607;
t751 = -t128 / 0.4e1;
t487 = -pkin(3) * t520 - qJ(4) * t517 - pkin(2);
t475 = t515 * t487;
t380 = -pkin(9) * t620 + t475 + (-pkin(8) * t513 - pkin(4)) * t520;
t618 = t515 * t520;
t429 = pkin(8) * t618 + t487 * t513;
t397 = -pkin(9) * t624 + t429;
t208 = t380 * t516 + t397 * t519;
t750 = -t208 / 0.2e1;
t362 = -t463 * t517 + t464 * t520;
t332 = -pkin(3) * t622 - t362;
t217 = -pkin(4) * t426 + t332;
t749 = t217 / 0.2e1;
t381 = pkin(4) * t517 - pkin(9) * t618 + t434;
t623 = t513 * t520;
t400 = -pkin(9) * t623 + t435;
t218 = t381 * t519 - t400 * t516;
t748 = -t218 / 0.2e1;
t626 = t461 * t513;
t223 = -pkin(4) * t626 + t318;
t747 = -t223 / 0.2e1;
t635 = t515 * Ifges(5,4);
t636 = t513 * Ifges(5,2);
t574 = t635 - t636;
t254 = Ifges(5,6) * t462 - t461 * t574;
t742 = t254 / 0.2e1;
t684 = Ifges(5,4) * t513;
t691 = Ifges(5,1) * t515;
t578 = -t684 + t691;
t255 = Ifges(5,5) * t462 - t461 * t578;
t741 = t255 / 0.2e1;
t260 = Ifges(5,1) * t427 + Ifges(5,4) * t426 + Ifges(5,5) * t607;
t739 = t260 / 0.4e1;
t640 = t461 * mrSges(5,2);
t301 = mrSges(5,3) * t551 - t640;
t738 = -t301 / 0.2e1;
t641 = t461 * mrSges(5,1);
t302 = mrSges(5,3) * t393 + t641;
t737 = -t302 / 0.2e1;
t448 = t564 * t517;
t447 = t516 * t624 - t517 * t619;
t644 = t447 * Ifges(6,4);
t326 = -Ifges(6,2) * t448 - Ifges(6,6) * t520 - t644;
t736 = -t326 / 0.2e1;
t735 = -t326 / 0.4e1;
t441 = Ifges(6,4) * t448;
t327 = -Ifges(6,1) * t447 - Ifges(6,5) * t520 - t441;
t734 = -t327 / 0.2e1;
t337 = mrSges(6,1) * t447 + mrSges(6,2) * t448;
t733 = -t337 / 0.2e1;
t683 = Ifges(6,4) * t564;
t390 = -Ifges(6,2) * t563 + t683;
t731 = t390 / 0.2e1;
t730 = t390 / 0.4e1;
t469 = Ifges(6,4) * t563;
t392 = Ifges(6,1) * t564 - t469;
t729 = t392 / 0.2e1;
t728 = t392 / 0.4e1;
t694 = pkin(9) + qJ(4);
t488 = t694 * t513;
t490 = t694 * t515;
t403 = -t488 * t519 - t490 * t516;
t727 = t403 / 0.2e1;
t413 = mrSges(6,2) * t520 - mrSges(6,3) * t448;
t726 = -t413 / 0.2e1;
t725 = t413 / 0.2e1;
t645 = t447 * mrSges(6,3);
t414 = -mrSges(6,1) * t520 + t645;
t724 = -t414 / 0.2e1;
t723 = -t447 / 0.2e1;
t722 = t447 / 0.2e1;
t721 = t448 / 0.2e1;
t720 = -t461 / 0.2e1;
t718 = t461 / 0.2e1;
t717 = t461 / 0.4e1;
t715 = -t563 / 0.2e1;
t714 = t564 / 0.2e1;
t713 = -t564 / 0.2e1;
t480 = mrSges(5,2) * t520 - mrSges(5,3) * t624;
t712 = -t480 / 0.2e1;
t481 = -mrSges(5,1) * t520 - mrSges(5,3) * t620;
t711 = -t481 / 0.2e1;
t598 = pkin(4) * t513 + pkin(8);
t482 = t598 * t517;
t710 = t482 / 0.2e1;
t483 = t598 * t520;
t709 = -t483 / 0.2e1;
t509 = -pkin(4) * t515 - pkin(3);
t707 = t509 / 0.2e1;
t706 = -t513 / 0.2e1;
t705 = t513 / 0.2e1;
t704 = -t514 / 0.2e1;
t703 = t515 / 0.2e1;
t702 = -t520 / 0.4e1;
t701 = t520 / 0.2e1;
t700 = t520 / 0.4e1;
t699 = t521 / 0.2e1;
t698 = m(5) * qJ(4);
t697 = m(5) * t332;
t696 = pkin(8) * t517;
t695 = pkin(8) * t520;
t690 = Ifges(5,1) * t512;
t689 = Ifges(3,4) * t518;
t688 = Ifges(3,4) * t521;
t687 = Ifges(4,4) * t462;
t686 = Ifges(4,4) * t517;
t685 = Ifges(4,4) * t520;
t682 = Ifges(4,5) * t518;
t681 = Ifges(4,5) * t520;
t680 = Ifges(5,5) * t393;
t678 = Ifges(5,5) * t515;
t677 = Ifges(5,5) * t517;
t676 = Ifges(6,5) * t224;
t675 = Ifges(6,5) * t295;
t329 = t563 * t461;
t674 = Ifges(6,5) * t329;
t673 = Ifges(6,5) * t447;
t450 = t563 * t520;
t672 = Ifges(6,5) * t450;
t670 = Ifges(5,6) * t551;
t669 = Ifges(5,6) * t513;
t667 = Ifges(5,6) * t517;
t666 = Ifges(6,6) * t783;
t665 = Ifges(6,6) * t294;
t328 = t564 * t461;
t664 = Ifges(6,6) * t328;
t663 = Ifges(6,6) * t448;
t449 = t564 * t520;
t662 = Ifges(6,6) * t449;
t660 = Ifges(4,3) * t518;
t659 = Ifges(6,3) * t462;
t658 = Ifges(6,3) * t517;
t657 = t783 * mrSges(6,2);
t656 = t224 * mrSges(6,1);
t654 = t294 * mrSges(6,1);
t653 = t295 * mrSges(6,2);
t537 = -pkin(2) * t633 - t463;
t270 = pkin(3) * t461 - qJ(4) * t462 + t537;
t277 = -qJ(4) * t621 + t318;
t124 = t270 * t515 - t277 * t513;
t125 = t270 * t513 + t277 * t515;
t148 = t653 - t654;
t189 = -mrSges(6,2) * t461 + mrSges(6,3) * t783;
t190 = mrSges(6,1) * t461 + mrSges(6,3) * t224;
t281 = pkin(3) * t621 - t317;
t194 = -pkin(4) * t551 + t281;
t244 = -mrSges(6,2) * t607 + mrSges(6,3) * t294;
t245 = mrSges(6,1) * t607 - mrSges(6,3) * t295;
t646 = t427 * mrSges(5,2);
t647 = t426 * mrSges(5,1);
t315 = t646 - t647;
t374 = -mrSges(5,2) * t607 + mrSges(5,3) * t426;
t375 = mrSges(5,1) * t607 - mrSges(5,3) * t427;
t80 = pkin(4) * t461 + pkin(9) * t393 + t124;
t98 = pkin(9) * t551 + t125;
t46 = -t516 * t98 + t519 * t80;
t47 = t516 * t80 + t519 * t98;
t532 = t537 * (mrSges(4,1) * t517 + mrSges(4,2) * t520);
t536 = t517 * (-Ifges(4,2) * t461 - Ifges(4,6) * t621 + t687);
t570 = t666 - t676;
t539 = t517 * (Ifges(6,3) * t461 + t570);
t540 = t517 * (Ifges(5,3) * t461 + t670 - t680);
t456 = Ifges(4,4) * t461;
t541 = t520 * (Ifges(4,1) * t462 - Ifges(4,5) * t621 - t456);
t575 = -Ifges(5,4) * t393 + Ifges(5,2) * t551;
t545 = Ifges(5,6) * t461 + t575;
t579 = -Ifges(5,1) * t393 + Ifges(5,4) * t551;
t546 = Ifges(5,5) * t461 + t579;
t549 = t521 * (-Ifges(4,6) * t517 + t681);
t639 = t461 * mrSges(4,3);
t558 = mrSges(4,2) * t621 - t639;
t638 = t462 * mrSges(4,3);
t560 = -mrSges(4,1) * t621 - t638;
t572 = -Ifges(5,5) * t427 - Ifges(5,6) * t426;
t576 = -Ifges(4,2) * t517 + t685;
t580 = Ifges(4,1) * t520 - t686;
t585 = -mrSges(6,1) * t783 - mrSges(6,2) * t224;
t587 = -mrSges(5,1) * t551 - mrSges(5,2) * t393;
t590 = Ifges(6,3) * t607;
t595 = t621 / 0.2e1;
t596 = -t621 / 0.2e1;
t65 = t121 * t516 + t138 * t519;
t3 = t393 * t260 / 0.2e1 - m(4) * (t317 * t362 + t318 * t363 + t465 * t537) + t462 * (t521 * t580 + t682) * t704 + t536 * t595 + ((t549 + t660) * t699 + pkin(1) * (mrSges(3,1) * t518 + mrSges(3,2) * t521)) * t770 + (t541 + t540 + t539) * t596 + (Ifges(5,3) * t607 - t572 + t590 + t665 + t675) * t720 - t465 * (mrSges(4,1) * t461 + mrSges(4,2) * t462) - t125 * t374 - t124 * t375 - t281 * t315 - t295 * t91 / 0.2e1 - t188 * t301 - t187 * t302 - t47 * t244 - t46 * t245 - t194 * t148 - t65 * t189 - t64 * t190 + t128 * t786 - (Ifges(4,5) * t462 - Ifges(4,6) * t461 - Ifges(4,3) * t621) * t622 / 0.2e1 - t532 * t621 - (t521 * (-Ifges(3,2) * t518 + t688) + t518 * (Ifges(3,1) * t521 - t689)) * t770 / 0.2e1 - t551 * t259 / 0.2e1 - m(5) * (t124 * t187 + t125 * t188 + t281 * t332) - m(6) * (t194 * t217 + t46 * t64 + t47 * t65) + ((Ifges(4,6) * t518 + t521 * t576) * t718 - t317 * (mrSges(4,1) * t518 - mrSges(4,3) * t615) - t318 * (-mrSges(4,2) * t518 - mrSges(4,3) * t616) + (Ifges(3,1) * t518 + t688) * t596 + (Ifges(3,2) * t521 + t689) * t622 / 0.2e1) * t514 + t294 * t757 - t426 * t545 / 0.2e1 - t427 * t546 / 0.2e1 + ((Ifges(3,5) * t521 - Ifges(3,6) * t518) * t704 + Ifges(3,5) * t596 + t779 / 0.2e1 + t771) * t633 + t127 * t784 - t363 * t558 - t362 * t560 - t217 * t585 - t332 * t587;
t652 = t3 * qJD(1);
t651 = t328 * mrSges(6,1);
t650 = t329 * mrSges(6,2);
t136 = Ifges(6,4) * t329 + Ifges(6,2) * t328 + Ifges(6,6) * t462;
t137 = Ifges(6,1) * t329 + Ifges(6,4) * t328 + Ifges(6,5) * t462;
t168 = t650 - t651;
t364 = pkin(3) * t462 + qJ(4) * t461;
t180 = -t317 * t513 + t364 * t515;
t181 = t317 * t515 + t364 * t513;
t234 = -mrSges(6,2) * t462 + mrSges(6,3) * t328;
t235 = mrSges(6,1) * t462 - mrSges(6,3) * t329;
t586 = t692 + t693;
t336 = t586 * t461;
t354 = -mrSges(5,2) * t462 + mrSges(5,3) * t626;
t625 = t461 * t515;
t355 = mrSges(5,1) * t462 + mrSges(5,3) * t625;
t455 = Ifges(4,5) * t461;
t571 = -t669 + t678;
t547 = t571 * t461;
t555 = t674 / 0.2e1 + t664 / 0.2e1;
t591 = Ifges(4,2) - Ifges(4,1) - t780;
t118 = pkin(4) * t462 + pkin(9) * t625 + t180;
t139 = pkin(9) * t626 + t181;
t60 = t118 * t519 - t139 * t516;
t61 = t118 * t516 + t139 * t519;
t4 = t551 * t742 + t318 * t587 - t393 * t741 + t125 * t354 + t124 * t355 - t281 * t336 + t328 * t90 / 0.2e1 + t329 * t756 + t181 * t301 + t180 * t302 + t47 * t234 + t46 * t235 + t783 * t136 / 0.2e1 + t223 * t585 - t224 * t137 / 0.2e1 + t194 * t168 + t61 * t189 + t60 * t190 + (-t680 / 0.2e1 + t670 / 0.2e1 - t676 / 0.2e1 + t666 / 0.2e1 + t537 * mrSges(4,1) - t687) * t462 + m(5) * (t124 * t180 + t125 * t181 + t281 * t318) + m(6) * (t194 * t223 + t46 * t60 + t47 * t61) + (t455 / 0.2e1 + t776) * t621 + (-t537 * mrSges(4,2) + Ifges(4,5) * t595 + t456 - t515 * t579 / 0.2e1 + t575 * t705 - t547 + t591 * t462 + t555) * t461;
t649 = t4 * qJD(1);
t404 = -t488 * t516 + t490 * t519;
t648 = t404 * mrSges(6,3);
t643 = t449 * mrSges(6,1);
t642 = t450 * mrSges(6,2);
t637 = t564 * mrSges(6,3);
t105 = -t656 + t657;
t106 = Ifges(6,5) * t783 + Ifges(6,6) * t224;
t107 = Ifges(6,2) * t224 + t216;
t108 = Ifges(6,1) * t783 + t655;
t9 = t106 * t718 - t47 * t190 + t46 * t189 + t194 * t105 - (t757 - t47 * mrSges(6,3) + t108 / 0.2e1) * t224 + (t107 / 0.2e1 - t46 * mrSges(6,3) + t756) * t783;
t634 = t9 * qJD(1);
t17 = t783 * t189 + t224 * t190 + t551 * t301 + t393 * t302 + m(6) * (t224 * t46 + t47 * t783) + m(5) * (t124 * t393 + t125 * t551);
t632 = qJD(1) * t17;
t631 = t318 * t517;
t630 = t393 * t513;
t629 = t551 * t515;
t617 = t517 * t336;
t614 = -Ifges(6,5) * t448 + Ifges(6,6) * t447;
t387 = -Ifges(6,5) * t563 - Ifges(6,6) * t564;
t604 = t660 / 0.2e1;
t603 = t187 * t760;
t601 = -t635 / 0.2e1;
t600 = t107 / 0.4e1 + t781;
t599 = t782 - t108 / 0.4e1;
t338 = Ifges(6,2) * t447 - t441;
t594 = t327 / 0.4e1 + t338 / 0.4e1;
t339 = -Ifges(6,1) * t448 + t644;
t593 = t339 / 0.4e1 + t735;
t592 = t388 / 0.4e1 + t491 / 0.4e1;
t584 = mrSges(6,1) * t448 - mrSges(6,2) * t447;
t583 = -t642 + t643;
t582 = mrSges(6,1) * t517 + mrSges(6,3) * t450;
t581 = -mrSges(6,2) * t517 - mrSges(6,3) * t449;
t577 = -Ifges(6,1) * t450 - Ifges(6,4) * t449;
t573 = -Ifges(6,4) * t450 - Ifges(6,2) * t449;
t569 = t664 + t674;
t568 = t663 + t673;
t207 = t380 * t519 - t397 * t516;
t219 = t381 * t516 + t400 * t519;
t428 = -pkin(8) * t623 + t475;
t548 = t586 * t520;
t557 = -mrSges(5,2) * t517 - mrSges(5,3) * t623;
t559 = mrSges(5,1) * t517 - mrSges(5,3) * t618;
t522 = -t532 / 0.2e1 - t482 * t168 / 0.2e1 + t447 * t137 / 0.4e1 + t448 * t136 / 0.4e1 - t428 * t355 / 0.2e1 - t429 * t354 / 0.2e1 - t329 * t327 / 0.4e1 - t207 * t235 / 0.2e1 - t219 * t189 / 0.2e1 - t783 * (Ifges(6,6) * t517 + t573) / 0.4e1 + t450 * t781 + t449 * t782 - t586 * t631 / 0.2e1 + t254 * t624 / 0.4e1 - t255 * t620 / 0.4e1 + t585 * t709 + t180 * t711 + t181 * t712 + t60 * t724 + t61 * t726 + t328 * t735 + t434 * t737 + t435 * t738 + t584 * t747 + t190 * t748 + t234 * t750 + (t194 * t483 + t207 * t60 + t208 * t61 + t218 * t46 + t219 * t47 + t223 * t482) * t766 + (t124 * t434 + t125 * t435 + t180 * t428 + t181 * t429 + (t281 * t520 + t631) * pkin(8)) * t768 + t224 * (Ifges(6,5) * t517 + t577) / 0.4e1 - t281 * t548 / 0.2e1 - t125 * t557 / 0.2e1 - t124 * t559 / 0.2e1 - t47 * t581 / 0.2e1 - t46 * t582 / 0.2e1 - t194 * t583 / 0.2e1;
t386 = mrSges(6,1) * t563 + mrSges(6,2) * t564;
t489 = -mrSges(5,1) * t515 + mrSges(5,2) * t513;
t492 = Ifges(5,2) * t515 + t684;
t493 = Ifges(5,1) * t513 + t635;
t524 = (-pkin(3) * t332 + (-t187 * t513 + t188 * t515) * qJ(4)) * t767 + (t217 * t509 + t403 * t64 + t404 * t65) * t765 - pkin(3) * t315 / 0.2e1 + t386 * t749 + t294 * t730 + t295 * t728 + t332 * t489 / 0.2e1 + t362 * mrSges(4,1) / 0.2e1 - t363 * mrSges(4,2) / 0.2e1 + t245 * t727 + t404 * t244 / 0.2e1 + t426 * t492 / 0.4e1 + t427 * t493 / 0.4e1 + t148 * t707;
t550 = t520 * t587;
t562 = -Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t1 = t524 + t522 + t569 * t700 + t513 * t739 + t513 * t603 + (t604 + ((mrSges(4,1) * t763 + Ifges(4,5)) * t520 + (-Ifges(4,6) + t753 / 0.2e1 + t592) * t517) * t521) * t514 + (t672 / 0.4e1 + t662 / 0.4e1 - t754 / 0.2e1 + (t690 / 0.4e1 + (t601 + t636 / 0.4e1) * t513 - t562) * t517 + (Ifges(4,4) / 0.2e1 - t571) * t520) * t461 + (t673 / 0.4e1 + t663 / 0.4e1 + t755 / 0.2e1 + (Ifges(4,4) - t678 / 0.4e1 + t669 / 0.4e1) * t517 + t562 * t520) * t462 + (-t667 / 0.2e1 + (t601 + t636 / 0.2e1) * t520) * t551 - (-t677 / 0.2e1 + (-t691 / 0.2e1 + t684 / 0.2e1) * t520) * t393 - t517 * t570 / 0.4e1 - (mrSges(6,3) * t758 + t751) * t564 + t456 * t701 - (t752 + t65 * mrSges(6,3) / 0.2e1) * t563 + (-t550 / 0.2e1 + t617 / 0.2e1) * pkin(8) + (t374 * t703 + t375 * t706) * qJ(4) + t772 * t515;
t554 = -t672 / 0.2e1 - t662 / 0.2e1;
t561 = -Ifges(4,4) + t571;
t14 = -t219 * t413 - t208 * t581 - t207 * t582 + t573 * t721 + t577 * t722 - 0.2e1 * t548 * t696 - t435 * t480 - t449 * t736 - t450 * t734 - t434 * t481 - t428 * t559 - t429 * t557 - t218 * t414 - t482 * t583 - t483 * t584 - m(6) * (t207 * t218 + t208 * t219 + t482 * t483) - m(5) * (t428 * t434 + t429 * t435) + (t520 * t561 + t554 + t754) * t520 + (t755 - t561 * t517 + (-t690 - m(5) * pkin(8) ^ 2 + (0.2e1 * t635 - t636) * t513 + t591) * t520 + t568) * t517;
t567 = -t1 * qJD(1) - t14 * qJD(2);
t22 = t614 * t701 + t482 * t337 - t207 * t413 + (t339 / 0.2e1 + t736) * t447 - (t734 - t338 / 0.2e1 + t207 * mrSges(6,3)) * t448 + (t414 - t645) * t208;
t523 = -t600 * t448 + t594 * t783 - t593 * t224 + t599 * t447 + (t207 * t784 + t208 * t786 + t46 * t721 + t47 * t722) * mrSges(6,3) + t194 * t733 + t207 * t189 / 0.2e1 + t190 * t750 + t46 * t725 + t614 * t717 + t47 * t724 + t105 * t710 + t106 * t702;
t556 = -t675 / 0.2e1 - t665 / 0.2e1;
t531 = mrSges(6,1) * t758 - t65 * mrSges(6,2) / 0.2e1 + t590 / 0.2e1 - t556;
t6 = t523 - t531;
t566 = qJD(1) * t6 - qJD(2) * t22;
t525 = (t393 * t428 + t551 * t429 + (-t124 * t515 - t125 * t513) * t517) * t768 + (t207 * t224 + t208 * t783 + t447 * t46 - t448 * t47) * t766 + t224 * t724 + t783 * t726 + t393 * t711 + t551 * t712 + t190 * t723 + t189 * t721;
t530 = t697 / 0.2e1 + m(6) * t749 - t654 / 0.2e1 + t653 / 0.2e1 - t647 / 0.2e1 + t646 / 0.2e1;
t10 = (t301 * t705 + t302 * t703) * t517 + t525 + t530;
t34 = t447 * t414 - t448 * t413 + m(6) * (t207 * t447 - t208 * t448) + (-t513 * t480 - t515 * t481 + m(5) * (-t428 * t515 - t429 * t513)) * t517;
t565 = -qJD(1) * t10 + qJD(2) * t34;
t389 = -Ifges(6,2) * t564 - t469;
t553 = -t403 * mrSges(6,3) / 0.2e1 + t389 / 0.4e1 + t728;
t391 = -Ifges(6,1) * t563 - t683;
t552 = t648 / 0.2e1 - t391 / 0.4e1 + t730;
t385 = mrSges(6,1) * t564 - mrSges(6,2) * t563;
t526 = t385 * t710 + t387 * t702 + t403 * t725 + t404 * t724 + t447 * t552 - t448 * t553 + t509 * t733 - t563 * t594 + t564 * t593;
t533 = -t658 / 0.2e1 + mrSges(6,1) * t748 + t219 * t761 - t554;
t19 = t526 + t533;
t30 = t404 * t637 - t509 * t385 - (-t392 / 0.2e1 - t389 / 0.2e1) * t563 - (t648 + t391 / 0.2e1 - t390 / 0.2e1) * t564;
t527 = -t600 * t563 - t599 * t564 + t553 * t783 + t552 * t224 + t194 * t385 / 0.2e1 + t189 * t727 - t404 * t190 / 0.2e1 + t387 * t717 + t105 * t707;
t534 = -t659 / 0.2e1 - t60 * mrSges(6,1) / 0.2e1 + t61 * t761 - t555;
t8 = t527 + t534;
t544 = qJD(1) * t8 + qJD(2) * t19 - qJD(3) * t30;
t528 = (t224 * t713 + t715 * t783) * mrSges(6,3) + (-t124 * t513 + t125 * t515 + (t629 - t630) * qJ(4)) * t767 + (t224 * t403 + t404 * t783 - t46 * t564 - t47 * t563) * t765 + t189 * t715 + t190 * t713;
t535 = t318 * t768 + m(6) * t747 + t651 / 0.2e1 - t650 / 0.2e1;
t13 = (t301 / 0.2e1 + t551 * t759 + t640 / 0.2e1) * t515 + (t737 + t393 * t760 + t641 / 0.2e1) * t513 + t528 + t535;
t529 = (t447 * t713 - t448 * t715) * mrSges(6,3) + (-t428 * t513 + t429 * t515) * t767 + (-t207 * t564 - t208 * t563 + t403 * t447 - t404 * t448) * t765 + t413 * t715 + t414 * t713 + t481 * t706 + t480 * t703;
t538 = m(6) * t709 - t643 / 0.2e1 + t642 / 0.2e1;
t28 = (m(5) * t763 + t775) * t520 + t529 + t538;
t89 = (t563 ^ 2 + t564 ^ 2) * mrSges(6,3) + m(6) * (-t403 * t564 - t404 * t563) + (mrSges(5,3) + t698) * (t513 ^ 2 + t512);
t543 = qJD(1) * t13 + qJD(2) * t28 + qJD(3) * t89;
t51 = 0.2e1 * t785 + 0.2e1 * t787;
t542 = qJD(1) * t51 + qJD(2) * t337 - qJD(3) * t385;
t52 = t657 / 0.2e1 - t656 / 0.2e1 + t787 + t785;
t27 = t695 * t767 + mrSges(5,2) * t618 / 0.2e1 + mrSges(5,1) * t623 / 0.2e1 + t529 - t538;
t18 = t526 - t533;
t12 = t301 * t703 + t302 * t706 + (t629 / 0.2e1 - t630 / 0.2e1) * mrSges(5,3) + t528 - t535 + t775 * t461;
t11 = t620 * t737 + t624 * t738 - t525 + t530;
t7 = t527 - t534;
t5 = t523 + t531;
t2 = t540 / 0.4e1 + t539 / 0.4e1 - t393 * (t520 * t578 + t677) / 0.4e1 - t563 * t752 - t536 / 0.4e1 - t564 * t751 + t524 + t541 / 0.4e1 - t522 + (t739 + t603 - qJ(4) * t375 / 0.2e1) * t513 + (t64 * t713 + t65 * t715) * mrSges(6,3) - t462 * (Ifges(4,2) * t520 + t686) / 0.4e1 + t517 * (-Ifges(4,1) * t461 - t687) / 0.4e1 - (-Ifges(5,5) * t520 + t517 * t578) * t625 / 0.4e1 + (-Ifges(5,6) * t520 + t517 * t574) * t626 / 0.4e1 - t545 * t623 / 0.4e1 + t546 * t618 / 0.4e1 + (t604 + (t681 / 0.2e1 + (-Ifges(4,6) / 0.2e1 + t592) * t517) * t521 - t549 / 0.4e1) * t514 + (Ifges(5,3) * t517 + t571 * t520 + t658 - t662 - t672) * t717 + (Ifges(5,3) * t462 - t547 + t569 + t659) * t702 + (-t617 + t550) * pkin(8) / 0.2e1 - (t560 + t638) * t695 / 0.2e1 - (t558 + t639) * t696 / 0.2e1 - (Ifges(4,1) * t517 + t576 + t685) * t461 / 0.4e1 + t551 * (t520 * t574 + t667) / 0.4e1 + (-Ifges(4,2) * t462 - t456) * t700 + (mrSges(4,1) * t462 - mrSges(4,2) * t461) * t764 + (qJ(4) * t374 / 0.2e1 + t772) * t515 + (t571 * t517 + t780 * t520 - t568 + t580) * t462 / 0.4e1;
t15 = [-qJD(2) * t3 + qJD(3) * t4 + qJD(4) * t17 + qJD(5) * t9, t2 * qJD(3) + t11 * qJD(4) + t5 * qJD(5) - t652 + (t188 * t480 + t187 * t481 + t482 * t148 + Ifges(3,5) * t621 - t779 + t217 * t584 - t448 * t127 / 0.2e1 + t128 * t723 + t428 * t375 + t429 * t374 + t65 * t413 + t64 * t414 + t294 * t326 / 0.2e1 + t295 * t327 / 0.2e1 + t208 * t244 + t207 * t245 - t771 + 0.2e1 * (t187 * t428 + t188 * t429) * t767 + 0.2e1 * (t207 * t64 + t208 * t65 + t217 * t482) * t765 + 0.2e1 * m(4) * t465 * t764 + (t332 * t586 + t259 * t706 - t362 * mrSges(4,3) + t260 * t703 + t426 * t574 / 0.2e1 + t427 * t578 / 0.2e1 + t465 * mrSges(4,2) + (-m(4) * t362 + t315 + t697) * pkin(8) + (-t521 * t568 / 0.2e1 + t682 + (-Ifges(4,4) * t521 + t571 * t699) * t517 + t589 * mrSges(4,1)) * t514) * t517 + (-t465 * mrSges(4,1) + Ifges(4,6) * t622 + t363 * mrSges(4,3) + (m(4) * t363 - mrSges(4,2) * t622) * pkin(8) + t556 + t572 + (-t517 * t591 + t685 - t754) * t621) * t520) * qJD(2), t649 + t2 * qJD(2) + t12 * qJD(4) + t7 * qJD(5) + ((-pkin(3) * t318 + (-t180 * t513 + t181 * t515) * qJ(4)) * t767 + (t223 * t509 + t403 * t60 + t404 * t61) * t765) * t769 + (-t61 * t563 * mrSges(6,3) + pkin(3) * t336 + t136 * t715 + t137 * t714 + t509 * t168 + t223 * t386 + t404 * t234 + t403 * t235 + t318 * t489 + t328 * t731 + t329 * t729 - t60 * t637 - t455 + (mrSges(5,3) * t181 + qJ(4) * t354 + t493 * t720 + t742) * t515 + (-mrSges(5,3) * t180 - qJ(4) * t355 + t492 * t718 + t741) * t513 + t773 * t462 - t776) * qJD(3), qJD(2) * t11 + qJD(3) * t12 + qJD(5) * t52 + t632, t634 + t5 * qJD(2) + t7 * qJD(3) + t52 * qJD(4) + (-mrSges(6,1) * t47 - mrSges(6,2) * t46 + t106) * qJD(5); -qJD(3) * t1 - qJD(4) * t10 + qJD(5) * t6 + t652, -qJD(3) * t14 + qJD(4) * t34 - qJD(5) * t22, t27 * qJD(4) + t18 * qJD(5) + ((t218 * t403 + t219 * t404 + t483 * t509) * t765 + t774 * t698 / 0.2e1) * t769 + t567 + (t774 * mrSges(5,3) + t483 * t386 - t449 * t731 - t450 * t729 + t509 * t583 + t573 * t715 + t577 * t714 + (-t218 * t564 - t219 * t563 + t403 * t450 - t404 * t449) * mrSges(6,3) + (t668 / 0.2e1 + t679 / 0.2e1 + t671 / 0.2e1 - t661 / 0.2e1 + t403 * mrSges(6,1) - t404 * mrSges(6,2) - Ifges(4,6) + t753 - t586 * qJ(4) + t773) * t517 + (t578 * t705 - pkin(3) * t586 + Ifges(4,5) + t492 * t706 + (-m(5) * pkin(3) - mrSges(4,1) + t489) * pkin(8) + (t574 + t493) * t703) * t520) * qJD(3), qJD(3) * t27 + t565, t18 * qJD(3) + (-mrSges(6,1) * t208 - mrSges(6,2) * t207 + t614) * qJD(5) + t566; qJD(2) * t1 + qJD(4) * t13 + qJD(5) * t8 - t649, qJD(4) * t28 + qJD(5) * t19 - t567, qJD(4) * t89 - qJD(5) * t30, t543, (-mrSges(6,1) * t404 - mrSges(6,2) * t403 + t387) * qJD(5) + t544; qJD(2) * t10 - qJD(3) * t13 - qJD(5) * t51 - t632, -qJD(3) * t28 - qJD(5) * t337 - t565, qJD(5) * t385 - t543, 0, -t542; -qJD(2) * t6 - qJD(3) * t8 + qJD(4) * t51 - t634, -qJD(3) * t19 + qJD(4) * t337 - t566, -qJD(4) * t385 - t544, t542, 0;];
Cq = t15;
