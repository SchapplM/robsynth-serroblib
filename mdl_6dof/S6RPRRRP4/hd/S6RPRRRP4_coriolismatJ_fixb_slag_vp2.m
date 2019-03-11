% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRP4
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:07:04
% EndTime: 2019-03-09 06:07:29
% DurationCPUTime: 13.58s
% Computational Cost: add. (31412->624), mult. (61654->799), div. (0->0), fcn. (72749->8), ass. (0->348)
t384 = sin(qJ(5));
t386 = cos(qJ(5));
t610 = Ifges(7,4) * t386;
t466 = -Ifges(7,2) * t384 + t610;
t612 = Ifges(6,4) * t386;
t468 = -Ifges(6,2) * t384 + t612;
t777 = t468 + t466;
t611 = Ifges(7,4) * t384;
t470 = Ifges(7,1) * t386 - t611;
t613 = Ifges(6,4) * t384;
t472 = Ifges(6,1) * t386 - t613;
t776 = t470 + t472;
t619 = Ifges(7,5) + Ifges(6,5);
t704 = Ifges(7,6) + Ifges(6,6);
t718 = -t384 * t704 + t386 * t619;
t382 = sin(pkin(10));
t383 = cos(pkin(10));
t627 = sin(qJ(3));
t629 = cos(qJ(3));
t356 = -t382 * t627 + t383 * t629;
t357 = t382 * t629 + t383 * t627;
t385 = sin(qJ(4));
t628 = cos(qJ(4));
t452 = t385 * t356 + t357 * t628;
t318 = -t628 * t356 + t357 * t385;
t725 = t318 * t384;
t739 = -mrSges(7,2) * t452 + mrSges(7,3) * t725;
t741 = -mrSges(6,2) * t452 + mrSges(6,3) * t725;
t775 = (t741 / 0.2e1 + t739 / 0.2e1) * t384;
t381 = t386 ^ 2;
t550 = t384 ^ 2 + t381;
t537 = t628 * pkin(3);
t374 = -t537 - pkin(4);
t621 = t386 * pkin(5);
t360 = t374 - t621;
t622 = t385 * pkin(3);
t373 = pkin(9) + t622;
t363 = -mrSges(7,1) * t386 + t384 * mrSges(7,2);
t475 = mrSges(6,1) * t386 - t384 * mrSges(6,2);
t443 = t475 * t452;
t709 = t452 / 0.2e1;
t727 = mrSges(6,3) + mrSges(7,3);
t733 = t318 / 0.2e1;
t454 = -t727 * t550 * t733 - t443 / 0.2e1 + t363 * t709;
t491 = t550 * t318;
t660 = m(5) * pkin(3);
t544 = t660 / 0.2e1;
t663 = m(7) / 0.2e1;
t665 = m(6) / 0.2e1;
t555 = qJ(6) + t373;
t346 = t555 * t386;
t330 = t346 * t386;
t345 = t555 * t384;
t683 = t345 * t384 + t330;
t774 = (-t373 * t491 + t374 * t452) * t665 + (-t318 * t683 + t360 * t452) * t663 + (-t318 * t385 - t452 * t628) * t544 + t454;
t617 = pkin(7) + qJ(2);
t494 = t617 * t383;
t495 = t617 * t382;
t322 = -t494 * t627 - t495 * t629;
t271 = -t357 * pkin(8) + t322;
t323 = t494 * t629 - t495 * t627;
t272 = t356 * pkin(8) + t323;
t152 = t628 * t271 - t385 * t272;
t756 = t152 * mrSges(5,2);
t684 = t271 * t385 + t628 * t272;
t742 = -pkin(5) * t725 + t684;
t765 = t742 * t363;
t773 = t765 - t756;
t469 = Ifges(7,1) * t384 + t610;
t444 = t386 * t469;
t471 = Ifges(6,1) * t384 + t612;
t445 = t386 * t471;
t465 = Ifges(7,2) * t386 + t611;
t446 = t384 * t465;
t467 = Ifges(6,2) * t386 + t613;
t447 = t384 * t467;
t772 = (-t445 / 0.4e1 + t446 / 0.4e1 + t447 / 0.4e1 - t444 / 0.4e1) * t318 - t756 / 0.2e1 + t765 / 0.2e1;
t664 = -m(7) / 0.2e1;
t723 = t386 * t318;
t738 = mrSges(7,1) * t452 + mrSges(7,3) * t723;
t648 = t738 / 0.2e1;
t570 = t452 * t384;
t104 = pkin(5) * t570 - t152;
t771 = t104 * t742;
t770 = t346 * t739;
t769 = t360 * t742;
t616 = -qJ(6) - pkin(9);
t362 = t616 * t384;
t768 = t362 * t738;
t375 = -pkin(4) - t621;
t767 = t375 * t742;
t643 = -t452 / 0.2e1;
t766 = t550 * t643;
t764 = t742 * t663;
t763 = t777 * t452;
t762 = t776 * t452;
t499 = -t723 / 0.2e1;
t506 = t725 / 0.2e1;
t569 = t452 * t386;
t516 = -pkin(2) * t383 - pkin(1);
t326 = -pkin(3) * t356 + t516;
t154 = pkin(4) * t318 - pkin(9) * t452 + t326;
t77 = t386 * t154 - t384 * t684;
t437 = -qJ(6) * t569 + t77;
t51 = t318 * pkin(5) + t437;
t615 = mrSges(7,3) * t386;
t521 = t615 / 0.2e1;
t597 = t386 * mrSges(7,2);
t601 = t384 * mrSges(7,1);
t473 = t597 + t601;
t216 = t473 * t452;
t598 = t386 * mrSges(6,2);
t474 = t384 * mrSges(6,1) + t598;
t505 = t570 / 0.2e1;
t522 = -t615 / 0.2e1;
t543 = pkin(5) * t569;
t623 = pkin(5) * t384;
t728 = t569 / 0.2e1;
t731 = t384 / 0.4e1;
t746 = t469 + t471;
t747 = t465 + t467;
t748 = t318 * t619 + t762;
t749 = t318 * t704 + t763;
t751 = t718 * t318;
t669 = t437 * t522 + t152 * t474 / 0.2e1 - t104 * t473 / 0.2e1 - t363 * t543 / 0.2e1 - t216 * t623 / 0.2e1 - t751 / 0.4e1 + t747 * t728 + t746 * t505 + (t749 + t763) * t731 - (t748 + t762) * t386 / 0.4e1;
t726 = Ifges(6,3) + Ifges(7,3);
t759 = t499 * t619 + t506 * t704 + t51 * t521 + t709 * t726 + t669;
t631 = t386 / 0.2e1;
t633 = t384 / 0.2e1;
t740 = mrSges(6,1) * t452 + mrSges(6,3) * t723;
t758 = (t740 + t738) * t631 + (t741 + t739) * t633;
t217 = t474 * t452;
t580 = t684 * t386;
t64 = t580 + (-qJ(6) * t452 + t154) * t384;
t722 = t473 * t318;
t78 = t384 * t154 + t580;
t757 = -t104 * t722 + t742 * t216 + t684 * t217 + t51 * t738 + t64 * t739 + t77 * t740 + t78 * t741;
t755 = qJ(6) * t725;
t754 = t152 * t385;
t583 = t152 * t684;
t753 = t384 * t152;
t752 = t386 * t152;
t536 = mrSges(6,3) * t570;
t225 = -mrSges(6,2) * t318 - t536;
t535 = mrSges(6,3) * t569;
t231 = t318 * mrSges(6,1) - t535;
t534 = mrSges(7,3) * t570;
t224 = -mrSges(7,2) * t318 - t534;
t533 = mrSges(7,3) * t569;
t230 = t318 * mrSges(7,1) - t533;
t730 = -t386 / 0.2e1;
t453 = t224 * t730 + t230 * t633;
t634 = -t384 / 0.2e1;
t418 = t437 - t51;
t707 = m(7) * t418;
t750 = t225 * t631 + t231 * t634 + t633 * t707 - t453;
t424 = t550 * t628;
t365 = t616 * t386;
t462 = -t384 * t51 + t386 * t64;
t745 = t453 + ((-t362 * t386 + t365 * t384) * t452 + t462) * t664;
t744 = t453 + ((t345 * t386 - t346 * t384) * t452 + t462) * t664;
t712 = Ifges(6,5) * t452 - t318 * t472;
t713 = Ifges(7,5) * t452 - t318 * t470;
t497 = t712 / 0.2e1 + t713 / 0.2e1;
t714 = Ifges(6,6) * t452 - t318 * t468;
t715 = Ifges(7,6) * t452 - t318 * t466;
t498 = -t715 / 0.2e1 - t714 / 0.2e1;
t737 = pkin(5) * t452 + qJ(6) * t723;
t735 = -t318 / 0.2e1;
t732 = t375 / 0.2e1;
t602 = t318 * mrSges(5,3);
t724 = t384 * (t231 + t230);
t721 = t474 * t318;
t620 = Ifges(6,4) + Ifges(7,4);
t496 = t620 * t386;
t717 = t496 - (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t384;
t624 = pkin(4) * t452;
t234 = pkin(9) * t318 + t624;
t599 = t384 * mrSges(7,3);
t517 = t599 / 0.2e1;
t518 = -t599 / 0.2e1;
t600 = t384 * mrSges(6,3);
t519 = t600 / 0.2e1;
t520 = -t600 / 0.2e1;
t527 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t711 = t318 * (t384 * (Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1) - t386 * (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1)) + t527 * t452 + t51 * t522 - t669 + (t519 + t520) * t78 + (t517 + t518) * t64;
t708 = Ifges(5,4) * t643;
t705 = -mrSges(6,2) - mrSges(7,2);
t699 = mrSges(5,3) * t452;
t694 = -mrSges(5,1) - t475;
t657 = m(7) * pkin(5);
t526 = mrSges(7,1) + t657;
t693 = -mrSges(6,1) - t526;
t687 = t216 + t217;
t686 = t224 + t225;
t358 = t365 * t386;
t682 = -t362 * t384 - t358;
t366 = Ifges(7,5) * t384 + Ifges(7,6) * t386;
t367 = Ifges(6,5) * t384 + Ifges(6,6) * t386;
t679 = t367 / 0.2e1 + t366 / 0.2e1;
t301 = mrSges(7,2) * t570;
t122 = t526 * t569 - t301;
t344 = t384 * t526 + t597;
t677 = qJD(1) * t122 + (qJD(3) + qJD(4)) * t344;
t674 = mrSges(6,3) * t766 + t225 * t634 + t231 * t730;
t88 = t386 * t234 - t753;
t58 = t737 + t88;
t651 = t58 / 0.2e1;
t654 = -mrSges(7,2) / 0.2e1;
t655 = -mrSges(6,2) / 0.2e1;
t656 = mrSges(6,1) / 0.2e1;
t89 = t384 * t234 + t752;
t70 = t89 + t755;
t673 = mrSges(7,1) * t651 + t654 * t70 + t655 * t89 + t656 * t88 + (m(7) * t651 + t648) * pkin(5);
t625 = pkin(3) * t357;
t163 = t234 + t625;
t81 = t386 * t163 - t753;
t54 = t737 + t81;
t652 = t54 / 0.2e1;
t82 = t384 * t163 + t752;
t67 = t82 + t755;
t672 = mrSges(7,1) * t652 + t654 * t67 + t655 * t82 + t656 * t81 + (m(7) * t652 + t648) * pkin(5);
t667 = t357 ^ 2;
t666 = -m(6) / 0.2e1;
t662 = -pkin(4) / 0.2e1;
t661 = pkin(4) / 0.2e1;
t659 = m(6) * pkin(3);
t658 = m(7) * pkin(3);
t653 = -Ifges(5,4) / 0.2e1;
t650 = pkin(4) * mrSges(6,1);
t640 = -t345 / 0.2e1;
t639 = t360 / 0.2e1;
t638 = t365 / 0.2e1;
t635 = t374 / 0.2e1;
t626 = m(7) * t452;
t609 = Ifges(5,5) * t318;
t608 = Ifges(5,6) * t452;
t607 = t684 * mrSges(5,1);
t316 = Ifges(5,4) * t318;
t233 = Ifges(5,1) * t452 - t316;
t342 = t357 * mrSges(4,1);
t477 = t708 + t718 * t709 + (Ifges(5,2) + t726) * t733;
t479 = Ifges(5,2) / 0.2e1 + t527;
t315 = t452 * mrSges(5,1);
t603 = t318 * mrSges(5,2);
t493 = t315 - t603;
t531 = Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1;
t3 = t757 + t516 * t342 + (t516 * mrSges(4,2) + Ifges(4,4) * t356 + (-Ifges(4,2) + Ifges(4,1)) * t357) * t356 + (m(5) * t625 + t493) * t326 + m(6) * (t77 * t81 + t78 * t82 - t583) + (mrSges(5,1) * t625 - t152 * mrSges(5,3) - (t653 + t718) * t318) * t318 + t233 * t735 - (-t721 - t602) * t152 + m(7) * (t51 * t54 + t64 * t67 + t771) + (t477 + mrSges(5,2) * t625 + Ifges(5,1) * t735 + (-t531 * t723 + t497) * t386 + (t717 * t318 + t498) * t384 + t479 * t318 + t452 * t653) * t452 - t667 * Ifges(4,4) + t67 * t224 + t82 * t225 + t54 * t230 + t81 * t231;
t604 = t3 * qJD(1);
t4 = t152 * t721 + t70 * t224 + t89 * t225 + t58 * t230 + t88 * t231 + m(6) * (t77 * t88 + t78 * t89 - t583) + m(7) * (t51 * t58 + t64 * t70 + t771) + (t326 * mrSges(5,1) + t498 * t384 + t497 * t386 + t477 + t708) * t452 + (-t326 * mrSges(5,2) + t316 / 0.2e1 - t233 / 0.2e1 - t751 + (-Ifges(5,1) / 0.2e1 - t531 * t381 + t717 * t384 + t479) * t452) * t318 + t757;
t596 = t4 * qJD(1);
t595 = t81 * t384;
t594 = t82 * t386;
t593 = t89 * t386;
t478 = mrSges(7,1) * t569 - t301;
t9 = t152 * t443 - t216 * t543 - t437 * t224 - t51 * t534 + (t231 + t535) * t78 + (-t536 - t225) * t77 + (-m(7) * t543 - t478) * t104 + (t230 + t533 - t707) * t64 + t748 * t505 + t749 * t728 + (t679 * t318 + (t631 * t746 + t634 * t747) * t452) * t452;
t592 = t9 * qJD(1);
t581 = t152 * t452;
t10 = (t356 ^ 2 + t667) * mrSges(4,3) + (t687 + t699) * t452 - (t386 * t686 - t602 - t724) * t318 + m(6) * (-t581 - (-t384 * t77 + t386 * t78) * t318) + m(7) * (t104 * t452 - t318 * t462) + m(5) * (-t318 * t684 - t581) + m(4) * (-t322 * t357 + t323 * t356) + (m(3) * qJ(2) + mrSges(3,3)) * (t382 ^ 2 + t383 ^ 2);
t591 = qJD(1) * t10;
t397 = (t384 * t82 + t386 * t81) * t665 + (t384 * t67 + t386 * t54) * t663 + t357 * t544;
t11 = t356 * mrSges(4,2) + t342 + t397 + t493 + t758 - t774;
t590 = qJD(1) * t11;
t395 = -t315 / 0.2e1 + (-pkin(9) * t491 - t624) * t665 + (-t318 * t682 + t375 * t452) * t663 + t454;
t399 = (t384 * t89 + t386 * t88) * t666 + (t384 * t70 + t386 * t58) * t664 + mrSges(5,1) * t643;
t15 = -t775 + (-t738 / 0.2e1 - t740 / 0.2e1) * t386 + 0.2e1 * t733 * mrSges(5,2) + t395 + t399;
t589 = qJD(1) * t15;
t27 = (t386 * t230 + t384 * t224 - m(7) * (-t384 * t64 - t51 * t386)) * t452;
t588 = qJD(1) * t27;
t586 = t104 * t384;
t500 = t723 / 0.2e1;
t553 = mrSges(7,1) * t506 + mrSges(7,2) * t500;
t12 = mrSges(6,2) * t500 + t553 + (mrSges(6,1) + t657) * t506 - t750;
t584 = t12 * qJD(1);
t582 = t684 * t475;
t576 = t452 * t366;
t575 = t452 * t367;
t566 = t384 * t713;
t565 = t384 * t712;
t564 = t384 * t740;
t561 = t386 * t715;
t560 = t386 * t714;
t559 = t386 * t741;
t552 = t550 * mrSges(7,3);
t548 = qJD(5) * t384;
t547 = qJD(5) * t386;
t456 = m(7) * t766;
t108 = -t626 / 0.2e1 + t456;
t546 = t108 * qJD(1);
t541 = mrSges(6,3) * t594;
t538 = -t628 / 0.2e1;
t532 = t622 / 0.2e1;
t515 = t384 * t628;
t514 = t386 * t628;
t492 = (t345 - t362) * t384;
t487 = mrSges(7,3) * pkin(5) - t619;
t486 = (t360 + t375) * t657;
t485 = -t537 / 0.2e1;
t482 = mrSges(7,3) * t505;
t481 = t452 * t522;
t480 = mrSges(7,2) * t538;
t476 = t764 - mrSges(7,1) * t725 / 0.2e1 + mrSges(7,2) * t499;
t461 = t594 - t595;
t388 = (t367 + t366) * t452 / 0.4e1 + t687 * t532 + t58 * t518 + t88 * t520 + t70 * t521 + Ifges(5,6) * t643 + t485 * t724 + (t374 * t665 - t475 / 0.2e1 - mrSges(5,1) / 0.2e1) * t684 + Ifges(5,5) * t735 + (t714 + t715) * t386 / 0.4e1 + (t712 + t713) * t731 - t722 * t639 - t721 * t635 + (-t345 * t58 + t346 * t70 + t769) * t663 + t770 / 0.2e1 + t738 * t640 + (t741 * t631 + t740 * t634 + (-t88 * t384 + t593) * t665) * t373 + ((t104 * t385 - t51 * t515 + t514 * t64) * t663 + (t514 * t78 - t515 * t77 - t754) * t665 + t686 * t514 / 0.2e1) * pkin(3) + mrSges(6,3) * t593 / 0.2e1 + t772;
t389 = -pkin(4) * t684 * t666 - t768 / 0.2e1 + t607 / 0.2e1 - t772 - t560 / 0.4e1 - t561 / 0.4e1 - t575 / 0.4e1 - t576 / 0.4e1 + t81 * t519 + t67 * t522 + (t461 * t666 - t559 / 0.2e1 + t564 / 0.2e1) * pkin(9) - t565 / 0.4e1 - t541 / 0.2e1 + t608 / 0.2e1 + t54 * t517 - t566 / 0.4e1 + t582 / 0.2e1 + t722 * t732 - t721 * t661 + t609 / 0.2e1 + t739 * t638 + (t362 * t54 - t365 * t67 + t767) * t664;
t1 = t388 + t389;
t392 = -mrSges(5,2) * t537 + (t363 + t694) * t622 + t727 * t424 * pkin(3);
t87 = (t345 * t515 + t346 * t514 + t360 * t385) * t658 + (t373 * t424 + t374 * t385) * t659 + t392;
t458 = t1 * qJD(1) + t87 * qJD(3);
t21 = t476 + t744;
t239 = m(7) * t683 + t552;
t457 = -qJD(1) * t21 + qJD(3) * t239;
t329 = t360 * t473;
t338 = t374 * t474;
t339 = t375 * t473;
t355 = t363 * t623;
t455 = -t329 / 0.2e1 - t338 / 0.2e1 - t339 / 0.2e1 - t355;
t188 = (t532 - t330 / 0.2e1 + t358 / 0.2e1 - t492 / 0.2e1) * m(7) - t552;
t406 = (-t597 / 0.2e1 - t601 / 0.2e1) * t318 + t764;
t24 = t406 + t745;
t279 = m(7) * t682 + t552;
t419 = -qJD(1) * t24 - qJD(3) * t188 + qJD(4) * t279;
t412 = t620 * t384 + (Ifges(6,2) - Ifges(7,1) + Ifges(7,2) - Ifges(6,1)) * t386;
t46 = -t329 - t338 - t355 - t620 * t381 + (-t360 * t657 + t412) * t384;
t390 = (t360 * t569 + t586) * pkin(5) * t663 + t224 * t640 + t478 * t639 + t443 * t635 - t345 * t482 + (t418 * t663 - t230 / 0.2e1 + t481) * t346 + t674 * t373;
t6 = -t390 + t672 + t759;
t417 = t6 * qJD(1) + t46 * qJD(3);
t402 = t693 * t538;
t36 = (pkin(3) * t480 - t496 + (t485 + t661) * mrSges(6,2)) * t386 + (t650 / 0.2e1 - t486 / 0.2e1 - t402 * pkin(3) + t412) * t384 + t455;
t62 = -t339 - t355 + (pkin(4) * mrSges(6,2) - t496) * t386 + (-t375 * t657 + t412 + t650) * t384;
t391 = (-t418 * t365 + (t375 * t569 + t586) * pkin(5)) * t663 + t230 * t638 + t478 * t732 + t443 * t662 - t365 * t481 + (t224 / 0.2e1 + t482) * t362 + t674 * pkin(9);
t7 = -t391 + t673 + t759;
t403 = t7 * qJD(1) + t36 * qJD(3) + t62 * qJD(4);
t394 = t447 / 0.2e1 + t446 / 0.2e1 - t445 / 0.2e1 - t444 / 0.2e1;
t337 = t344 * qJD(6);
t336 = t344 * qJD(5);
t189 = (t330 - t358 + t492) * t663 + m(7) * t532 + t552;
t107 = t626 / 0.2e1 + t456;
t37 = t474 * t662 + ((mrSges(6,2) * t538 + t480) * t386 - t402 * t384) * pkin(3) - t394 - t455 + t777 * t631 + (t486 + t776) * t633;
t25 = t406 - t745;
t22 = t476 - t744;
t16 = t603 / 0.2e1 + mrSges(5,2) * t735 + t395 - t399 + t758;
t14 = (t648 + t740 / 0.2e1) * t386 + t775 + t397 + t774;
t13 = -(-t598 / 0.2e1 + (-mrSges(6,1) / 0.2e1 - t657 / 0.2e1) * t384) * t318 + t553 + t750;
t8 = t391 + t673 + t711;
t5 = t390 + t672 + t711;
t2 = t388 - t389;
t17 = [qJD(2) * t10 + qJD(3) * t3 + qJD(4) * t4 - qJD(5) * t9 - qJD(6) * t27, qJD(3) * t14 + qJD(4) * t16 + qJD(5) * t13 + qJD(6) * t107 + t591, t604 + t14 * qJD(2) + ((m(6) * t461 + t559 - t564) * t373 - t607 + t560 / 0.2e1 + t561 / 0.2e1 + t575 / 0.2e1 + t576 / 0.2e1 + t773 - t622 * t699 + t67 * t615 + t565 / 0.2e1 + t541 + t770 - t608 + t566 / 0.2e1 - t582 + (t445 + t444) * t735 + (t447 + t446) * t733 - t360 * t722 + (m(6) * t684 - t721) * t374 - Ifges(4,6) * t357 + Ifges(4,5) * t356 - t323 * mrSges(4,1) - t322 * mrSges(4,2) - t609 - t345 * t738 + (-t628 * t684 + t754) * t660 - mrSges(6,3) * t595 - t54 * t599 + t537 * t602 + m(7) * (-t345 * t54 + t346 * t67 + t769)) * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t22 * qJD(6), t16 * qJD(2) + t2 * qJD(3) + t8 * qJD(5) + t25 * qJD(6) + t596 + (m(7) * (t362 * t58 - t365 * t70 + t767) - t375 * t722 - t365 * t739 + t768 + pkin(4) * t721 + (-m(6) * pkin(4) + t694) * t684 + (t89 * mrSges(6,3) + t70 * mrSges(7,3) + (m(6) * t89 + t741) * pkin(9) - t498) * t386 + (-t88 * mrSges(6,3) - t58 * mrSges(7,3) + (-m(6) * t88 - t740) * pkin(9) + t497) * t384 + (-Ifges(5,6) + t679) * t452 + (-Ifges(5,5) + t394) * t318 + t773) * qJD(4), t13 * qJD(2) + t5 * qJD(3) + t8 * qJD(4) - t592 + (-t78 * mrSges(6,1) + ((mrSges(7,2) * qJ(6) - t704) * t386 + t487 * t384) * t452 + t705 * t77 - t526 * t64) * qJD(5), qJD(2) * t107 + qJD(3) * t22 + qJD(4) * t25 - t588; qJD(3) * t11 - qJD(4) * t15 - qJD(5) * t12 + qJD(6) * t108 - t591, 0, t590, -t589, t547 * t705 + t548 * t693 - t584, t546; -qJD(2) * t11 + qJD(4) * t1 - qJD(5) * t6 - qJD(6) * t21 - t604, -t590, qJD(4) * t87 - qJD(5) * t46 + qJD(6) * t239 ((-t362 * t515 - t365 * t514 + t375 * t385) * t658 + (-pkin(4) * t385 + pkin(9) * t424) * t659 + t392) * qJD(4) + t37 * qJD(5) + t189 * qJD(6) + t458, t37 * qJD(4) + (t345 * mrSges(7,2) - t346 * t526) * qJD(5) + (mrSges(6,2) * t373 - t704) * t548 + (-mrSges(6,1) * t373 - t487) * t547 - t417, qJD(4) * t189 + t457; qJD(2) * t15 - qJD(3) * t1 - qJD(5) * t7 - qJD(6) * t24 - t596, t589, -qJD(5) * t36 - qJD(6) * t188 - t458, -qJD(5) * t62 + qJD(6) * t279 (-t362 * mrSges(7,2) + t365 * t526) * qJD(5) + (mrSges(6,2) * pkin(9) - t704) * t548 + (-mrSges(6,1) * pkin(9) - t487) * t547 - t403, t419; qJD(2) * t12 + qJD(3) * t6 + qJD(4) * t7 - qJD(6) * t122 + t592, t584, t36 * qJD(4) - t337 + t417, -t337 + t403, 0, -t677; -qJD(2) * t108 + qJD(3) * t21 + qJD(4) * t24 + qJD(5) * t122 + t588, -t546, qJD(4) * t188 + t336 - t457, t336 - t419, t677, 0;];
Cq  = t17;
