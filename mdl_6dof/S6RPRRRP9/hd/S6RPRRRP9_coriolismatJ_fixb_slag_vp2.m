% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRP9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:14
% EndTime: 2019-03-09 06:26:39
% DurationCPUTime: 15.68s
% Computational Cost: add. (19735->788), mult. (40692->1037), div. (0->0), fcn. (40015->6), ass. (0->417)
t457 = sin(qJ(4));
t459 = cos(qJ(4));
t701 = sin(qJ(5));
t702 = cos(qJ(5));
t490 = t457 * t702 + t459 * t701;
t633 = t490 * mrSges(7,2);
t768 = t701 * t457 - t702 * t459;
t636 = t768 * mrSges(7,1);
t279 = t633 + t636;
t693 = pkin(4) * t459;
t445 = -pkin(3) - t693;
t339 = pkin(5) * t768 + t445;
t697 = m(7) * t339;
t816 = t279 + t697;
t740 = -pkin(9) - pkin(8);
t430 = t740 * t457;
t432 = t740 * t459;
t299 = -t701 * t430 + t702 * t432;
t646 = mrSges(6,3) * t299;
t219 = t768 * qJ(6) + t299;
t653 = mrSges(7,3) * t219;
t815 = t646 + t653;
t798 = Ifges(6,6) + Ifges(7,6);
t810 = Ifges(6,5) + Ifges(7,5);
t508 = -t490 * t798 - t768 * t810;
t578 = t702 * t430 + t701 * t432;
t778 = t578 * mrSges(6,2);
t773 = -qJ(6) * t490 + t578;
t788 = t773 * mrSges(7,2);
t796 = t299 * mrSges(6,1);
t812 = t219 * mrSges(7,1);
t814 = t508 - t778 - t788 + t796 + t812;
t813 = -t778 / 0.2e1 - t788 / 0.2e1 + t796 / 0.2e1 + t812 / 0.2e1;
t460 = cos(qJ(3));
t691 = pkin(8) * t460;
t458 = sin(qJ(3));
t695 = pkin(3) * t458;
t421 = qJ(2) - t691 + t695;
t405 = t459 * t421;
t461 = -pkin(1) - pkin(7);
t522 = -t457 * t461 + pkin(4);
t588 = t459 * t460;
t571 = pkin(9) * t588;
t252 = t458 * t522 + t405 - t571;
t451 = t458 * t461;
t320 = t457 * t421 + t451 * t459;
t594 = t457 * t460;
t276 = -pkin(9) * t594 + t320;
t261 = t702 * t276;
t107 = t252 * t701 + t261;
t373 = t490 * t460;
t608 = t373 * qJ(6);
t79 = t107 - t608;
t319 = -t451 * t457 + t405;
t275 = t319 - t571;
t117 = -t275 * t701 - t261;
t90 = t117 + t608;
t789 = t90 + t79;
t259 = t701 * t276;
t106 = t702 * t252 - t259;
t371 = t768 * t460;
t340 = t371 * qJ(6);
t78 = t340 + t106;
t72 = t458 * pkin(5) + t78;
t652 = t219 * t72;
t341 = t373 * mrSges(7,2);
t224 = -t371 * mrSges(7,1) - t341;
t225 = -mrSges(6,1) * t371 - mrSges(6,2) * t373;
t441 = pkin(4) * t594;
t587 = t460 * t461;
t407 = t441 - t587;
t690 = t373 * pkin(5);
t263 = t407 + t690;
t640 = t373 * mrSges(7,3);
t305 = -mrSges(7,2) * t458 - t640;
t306 = -mrSges(6,2) * t458 - t373 * mrSges(6,3);
t706 = t458 / 0.4e1;
t710 = t445 / 0.2e1;
t310 = mrSges(6,1) * t458 + t371 * mrSges(6,3);
t730 = t310 / 0.2e1;
t309 = mrSges(7,1) * t458 + t371 * mrSges(7,3);
t731 = t309 / 0.2e1;
t278 = mrSges(6,1) * t490 - mrSges(6,2) * t768;
t736 = t278 / 0.2e1;
t396 = t768 * mrSges(7,2);
t277 = mrSges(7,1) * t490 - t396;
t737 = t277 / 0.2e1;
t811 = t219 * t731 + t299 * t730 + t773 * t305 / 0.2e1 + t578 * t306 / 0.2e1 + t263 * t737 + t339 * t224 / 0.2e1 + t407 * t736 + t225 * t710 + t508 * t706;
t809 = Ifges(6,3) + Ifges(7,3);
t118 = t702 * t275 - t259;
t808 = t106 - t118;
t431 = pkin(3) * t460 + pkin(8) * t458;
t414 = t459 * t431;
t593 = t458 * t459;
t264 = pkin(9) * t593 + t460 * t522 + t414;
t327 = t457 * t431 + t459 * t587;
t595 = t457 * t458;
t295 = pkin(9) * t595 + t327;
t110 = t702 * t264 - t295 * t701;
t111 = t701 * t264 + t702 * t295;
t372 = t490 * t458;
t303 = -mrSges(7,2) * t460 + mrSges(7,3) * t372;
t304 = -mrSges(6,2) * t460 + mrSges(6,3) * t372;
t374 = t768 * t458;
t307 = mrSges(7,1) * t460 - t374 * mrSges(7,3);
t308 = mrSges(6,1) * t460 - t374 * mrSges(6,3);
t326 = -t457 * t587 + t414;
t569 = t702 * pkin(4);
t444 = t569 + pkin(5);
t514 = t569 / 0.2e1;
t570 = pkin(4) * t701;
t517 = t570 / 0.2e1;
t629 = t458 * Ifges(5,5);
t541 = -t629 / 0.2e1;
t671 = Ifges(5,6) * t457;
t549 = t671 / 0.2e1;
t744 = m(6) * pkin(4);
t575 = t744 / 0.2e1;
t703 = t460 / 0.2e1;
t711 = t444 / 0.2e1;
t747 = m(7) / 0.2e1;
t723 = t374 / 0.2e1;
t689 = t460 * pkin(5);
t76 = -t374 * qJ(6) + t110 + t689;
t803 = -mrSges(7,1) / 0.2e1;
t804 = -mrSges(6,1) / 0.2e1;
t86 = qJ(6) * t372 + t111;
t757 = t111 * mrSges(6,2) / 0.2e1 + t110 * t804 + t86 * mrSges(7,2) / 0.2e1 + t76 * t803;
t802 = t372 / 0.2e1;
t753 = t703 * t809 + t723 * t810 + t798 * t802 - t757;
t807 = (t444 * t76 + t570 * t86) * t747 + Ifges(5,3) * t703 + t326 * mrSges(5,1) / 0.2e1 - t327 * mrSges(5,2) / 0.2e1 + t307 * t711 + (t110 * t702 + t111 * t701) * t575 + t459 * t541 + t458 * t549 + t308 * t514 + (t303 + t304) * t517 + t753;
t749 = m(6) / 0.2e1;
t801 = -mrSges(6,1) - mrSges(7,1);
t800 = mrSges(7,2) + mrSges(6,2);
t455 = t457 ^ 2;
t456 = t459 ^ 2;
t577 = t455 + t456;
t797 = mrSges(5,3) * t577;
t538 = t701 * t490;
t540 = t702 * t768;
t627 = t459 * mrSges(5,1);
t630 = t457 * mrSges(5,2);
t634 = t490 * mrSges(6,2);
t637 = t768 * mrSges(6,1);
t759 = -t636 / 0.2e1 - t637 / 0.2e1 - t633 / 0.2e1 - t634 / 0.2e1;
t519 = pkin(4) * t538;
t771 = t444 * t768 - t519;
t794 = -t771 * t747 - t630 / 0.2e1 + t627 / 0.2e1 + (t538 - t540) * t575 + t759;
t765 = t569 - t444;
t790 = mrSges(6,3) + mrSges(7,3);
t786 = t117 + t107;
t675 = Ifges(7,4) * t490;
t282 = -Ifges(7,2) * t768 + t675;
t677 = Ifges(6,4) * t490;
t284 = -Ifges(6,2) * t768 + t677;
t777 = t284 + t282;
t401 = Ifges(7,4) * t768;
t286 = Ifges(7,1) * t490 - t401;
t402 = Ifges(6,4) * t768;
t288 = Ifges(6,1) * t490 - t402;
t785 = t288 + t286;
t500 = -t326 * t457 + t327 * t459;
t779 = -t72 + t78;
t682 = mrSges(7,3) * t768;
t285 = -Ifges(7,1) * t768 - t675;
t287 = -Ifges(6,1) * t768 - t677;
t776 = t287 + t285;
t775 = t305 + t306;
t774 = t309 + t310;
t454 = Ifges(5,4) * t459;
t772 = -Ifges(5,2) * t457 + t454;
t426 = Ifges(5,1) * t457 + t454;
t418 = mrSges(5,1) * t458 - mrSges(5,3) * t588;
t589 = t459 * t418;
t416 = -mrSges(5,2) * t458 - mrSges(5,3) * t594;
t597 = t457 * t416;
t770 = -t589 / 0.2e1 - t597 / 0.2e1;
t767 = -t653 / 0.2e1 - t646 / 0.2e1;
t393 = t426 * t460;
t642 = t373 * mrSges(6,1);
t645 = t371 * mrSges(6,2);
t229 = t642 - t645;
t738 = t229 / 0.2e1;
t766 = pkin(4) * t738 - t393 / 0.4e1;
t509 = t371 * t798 - t373 * t810;
t322 = t371 * t570;
t764 = t444 * t373 + t322;
t762 = t630 - t627;
t281 = -Ifges(7,2) * t490 - t401;
t283 = -Ifges(6,2) * t490 - t402;
t760 = t283 + t281 + t785;
t704 = -t460 / 0.2e1;
t758 = t704 * t797 + t770;
t756 = t701 * t801 - t702 * t800;
t552 = t682 / 0.2e1;
t755 = (t219 * t747 + t552) * pkin(5) + t813;
t497 = t372 * t800 - t374 * t801;
t559 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t560 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t754 = t559 * t372 - t560 * t374;
t752 = 0.2e1 * m(7);
t751 = 2 * qJD(3);
t750 = m(5) / 0.2e1;
t748 = -m(7) / 0.2e1;
t746 = -pkin(5) / 0.2e1;
t745 = -pkin(8) / 0.2e1;
t743 = m(7) * pkin(5);
t742 = -mrSges(6,3) / 0.2e1;
t741 = -t72 / 0.2e1;
t739 = t118 / 0.2e1;
t280 = t634 + t637;
t735 = t280 / 0.2e1;
t572 = pkin(4) * t588;
t297 = -pkin(5) * t371 + t572;
t734 = t297 / 0.2e1;
t694 = pkin(4) * t457;
t370 = pkin(5) * t490 + t694;
t729 = t370 / 0.2e1;
t389 = t762 * t460;
t722 = t389 / 0.2e1;
t679 = Ifges(5,4) * t457;
t424 = Ifges(5,2) * t459 + t679;
t392 = t460 * t424;
t721 = -t392 / 0.4e1;
t719 = -t768 / 0.2e1;
t717 = t490 / 0.2e1;
t715 = -t490 / 0.2e1;
t713 = -t424 / 0.4e1;
t427 = Ifges(5,1) * t459 - t679;
t712 = t427 / 0.4e1;
t709 = -t457 / 0.2e1;
t708 = t457 / 0.2e1;
t707 = t458 / 0.2e1;
t705 = t459 / 0.2e1;
t700 = m(7) * (-t371 * t372 + t373 * t374);
t698 = m(7) * t263;
t696 = m(7) * t374;
t692 = pkin(5) * t307;
t687 = t78 * mrSges(7,2);
t686 = t79 * mrSges(7,1);
t684 = t90 * mrSges(7,1);
t91 = t340 + t118;
t683 = t91 * mrSges(7,2);
t681 = mrSges(7,3) * t490;
t678 = Ifges(6,4) * t371;
t676 = Ifges(7,4) * t371;
t453 = Ifges(5,5) * t459;
t666 = t106 * mrSges(6,2);
t665 = t107 * mrSges(6,1);
t662 = t117 * mrSges(6,1);
t661 = t118 * mrSges(6,2);
t657 = t773 * mrSges(7,3);
t651 = t219 * t78;
t644 = t371 * mrSges(7,2);
t643 = t372 * mrSges(7,1);
t641 = t373 * mrSges(7,1);
t639 = t373 * t72;
t638 = t374 * mrSges(7,2);
t635 = t768 * mrSges(6,3);
t632 = t490 * mrSges(6,3);
t631 = t457 * mrSges(5,1);
t628 = t458 * Ifges(5,6);
t626 = t459 * mrSges(5,2);
t195 = Ifges(7,4) * t374 + Ifges(7,2) * t372 + Ifges(7,6) * t460;
t197 = Ifges(6,4) * t374 + Ifges(6,2) * t372 + Ifges(6,6) * t460;
t199 = Ifges(7,1) * t374 + Ifges(7,4) * t372 + Ifges(7,5) * t460;
t201 = Ifges(6,1) * t374 + Ifges(6,4) * t372 + Ifges(6,5) * t460;
t226 = t638 - t643;
t227 = -mrSges(6,1) * t372 + mrSges(6,2) * t374;
t228 = t641 - t644;
t406 = -pkin(4) * t595 + t451;
t262 = -pkin(5) * t372 + t406;
t366 = Ifges(5,6) * t460 - t458 * t772;
t368 = t460 * Ifges(5,5) - t427 * t458;
t423 = t626 + t631;
t390 = t423 * t458;
t391 = t460 * t423;
t415 = -mrSges(5,2) * t460 + mrSges(5,3) * t595;
t417 = mrSges(5,1) * t460 + mrSges(5,3) * t593;
t491 = Ifges(4,4) - t453 / 0.2e1 + t549;
t348 = Ifges(7,4) * t373;
t200 = -Ifges(7,1) * t371 + t458 * Ifges(7,5) - t348;
t349 = Ifges(6,4) * t373;
t202 = -Ifges(6,1) * t371 + t458 * Ifges(6,5) - t349;
t526 = -t200 / 0.2e1 - t202 / 0.2e1;
t196 = -Ifges(7,2) * t373 + t458 * Ifges(7,6) - t676;
t198 = -Ifges(6,2) * t373 + t458 * Ifges(6,6) - t678;
t527 = t196 / 0.2e1 + t198 / 0.2e1;
t369 = t460 * t427 + t629;
t591 = t459 * t369;
t367 = t460 * t772 + t628;
t598 = t457 * t367;
t5 = (-qJ(2) * mrSges(4,2) + t461 * t391 - t591 / 0.2e1 + t598 / 0.2e1 + t491 * t458 + (-m(5) * t461 ^ 2 - Ifges(4,1) + Ifges(4,2) + Ifges(5,3) + t809) * t460 + t754) * t458 + (qJ(2) * mrSges(4,1) + t366 * t709 + t368 * t705 + t560 * t371 - t559 * t373 + t461 * t390 - t491 * t460) * t460 - t526 * t374 + t527 * t372 - (t195 / 0.2e1 + t197 / 0.2e1) * t373 + (-t199 / 0.2e1 - t201 / 0.2e1) * t371 + m(5) * (t319 * t326 + t320 * t327) + m(6) * (t106 * t110 + t107 * t111 + t406 * t407) + m(7) * (t262 * t263 + t72 * t76 + t79 * t86) + t320 * t415 + t327 * t416 + t319 * t417 + t326 * t418 + t406 * t229 + t407 * t227 + t79 * t303 + t107 * t304 + t86 * t305 + t111 * t306 + t72 * t307 + t106 * t308 + t76 * t309 + t110 * t310 + t262 * t228 + t263 * t226;
t625 = t5 * qJD(1);
t230 = Ifges(7,2) * t371 - t348;
t231 = Ifges(6,2) * t371 - t349;
t621 = t106 * t373;
t472 = t263 * t224 + t407 * t225 - (t230 / 0.2e1 + t231 / 0.2e1 - t526) * t373 + mrSges(7,3) * t639 + mrSges(6,3) * t621 + t509 * t707;
t232 = -Ifges(7,1) * t373 + t676;
t233 = -Ifges(6,1) * t373 + t678;
t483 = t79 * mrSges(7,3) + t107 * mrSges(6,3) - t233 / 0.2e1 - t232 / 0.2e1 + t527;
t6 = t483 * t371 + (t461 * t389 + (t541 + t392 / 0.2e1 - t369 / 0.2e1 + t319 * mrSges(5,3)) * t457 + (-t628 / 0.2e1 - t393 / 0.2e1 - t367 / 0.2e1 - t320 * mrSges(5,3) + (m(6) * t407 + t229) * pkin(4)) * t459) * t460 + m(6) * (t106 * t117 + t107 * t118) + m(7) * (t263 * t297 + t72 * t90 + t79 * t91) + t319 * t416 - t320 * t418 + t91 * t305 + t118 * t306 + t90 * t309 + t117 * t310 + t297 * t228 + t472;
t624 = t6 * qJD(1);
t7 = t78 * t305 + t106 * t306 - t107 * t310 + (m(7) * t779 - t309) * t79 + ((-t228 - t698) * pkin(5) + t483) * t371 + t472;
t623 = t7 * qJD(1);
t622 = -mrSges(4,1) + t762;
t620 = t106 * t768;
t59 = t374 * t72;
t473 = (-t460 ^ 2 * t693 + t374 * t808) * t749 + (-t297 * t460 - t374 * t91 + t59) * t747 + (-t747 * t789 - t749 * t786) * t372;
t610 = t372 * t373;
t612 = t371 * t374;
t12 = t473 - t775 * t372 / 0.2e1 + t774 * t723 + (-t389 + t225 + t224) * t704 + t758 * t458 + t790 * (-t612 / 0.2e1 - t610 / 0.2e1) - t794;
t619 = t12 * qJD(1);
t523 = t730 + t731;
t524 = -t306 / 0.2e1 - t305 / 0.2e1;
t561 = mrSges(6,3) / 0.2e1 + mrSges(7,3) / 0.2e1;
t475 = (-t373 * t561 + t524) * t372 - (t371 * t561 - t523) * t374;
t525 = -t225 / 0.2e1 - t224 / 0.2e1;
t467 = t525 * t460 + (t371 * t689 - t374 * t78 + t59) * t747 + t475;
t574 = -t743 / 0.2e1;
t499 = t804 + t803 + t574;
t562 = -mrSges(6,2) / 0.2e1 - mrSges(7,2) / 0.2e1;
t15 = -t490 * t562 - t499 * t768 + t467;
t618 = t15 * qJD(1);
t617 = t219 * t490;
t24 = t458 * mrSges(4,1) + t460 * mrSges(4,2) + t597 + t589 + mrSges(3,3) + t775 * t490 - t774 * t768 + (m(4) + m(3)) * qJ(2) + m(7) * (t490 * t79 - t72 * t768) + m(6) * (t107 * t490 - t620) + m(5) * (t319 * t459 + t320 * t457);
t616 = t24 * qJD(1);
t615 = t299 * t490;
t609 = t372 * t490;
t606 = t374 * t768;
t604 = t407 * t457;
t602 = t444 * t371;
t599 = t444 * t490;
t596 = t457 * t417;
t592 = t458 * t460;
t590 = t459 * t415;
t576 = qJD(3) * t458;
t573 = t743 / 0.2e1;
t568 = t700 / 0.2e1;
t565 = -t690 / 0.2e1;
t558 = mrSges(7,1) + t743;
t551 = -t682 / 0.2e1;
t545 = t578 * t742;
t543 = t640 / 0.2e1;
t528 = -t461 * t423 / 0.2e1;
t521 = t453 - t671;
t520 = t373 * t569;
t516 = t570 / 0.4e1;
t515 = -t569 / 0.2e1;
t511 = -mrSges(6,3) * t615 - mrSges(7,3) * t617 - t339 * t277 - t445 * t278;
t510 = -t641 / 0.2e1 - t642 / 0.2e1 + t644 / 0.2e1 + t645 / 0.2e1;
t507 = mrSges(5,3) * (-t456 / 0.2e1 - t455 / 0.2e1);
t506 = mrSges(6,3) * t520;
t45 = 0.4e1 * (m(6) / 0.4e1 + m(7) / 0.4e1) * (-t592 + t610 + t612) + m(5) * (-0.1e1 + t577) * t592;
t462 = -(t304 / 0.2e1 + t303 / 0.2e1) * t374 + t524 * t371 + (-t308 / 0.2e1 - t307 / 0.2e1) * t372 - t523 * t373 + (-t227 / 0.2e1 - t226 / 0.2e1 + t390 / 0.2e1 + t416 * t705 + t418 * t709) * t460 + (t738 + t228 / 0.2e1 + t391 / 0.2e1 + t590 / 0.2e1 - t596 / 0.2e1) * t458 + ((-t319 * t457 + t320 * t459) * t460 + (t500 - 0.2e1 * t587) * t458) * t750 + (-t107 * t371 - t110 * t372 - t111 * t374 - t406 * t460 + t407 * t458 - t621) * t749 + (-t262 * t460 + t263 * t458 - t371 * t79 - t372 * t76 - t374 * t86 - t639) * t747;
t479 = (-t578 * t768 - t615) * t749 + (-t768 * t773 - t617) * t747;
t9 = t462 - t479;
t502 = t9 * qJD(1) + t45 * qJD(2);
t56 = (-t373 * t516 + t602 / 0.4e1 - t297 / 0.4e1) * t752 - t224;
t75 = (-t768 * t516 - t599 / 0.4e1 - t370 / 0.4e1) * t752 - t277;
t501 = qJD(1) * t56 + qJD(3) * t75;
t132 = t371 * t558 + t341;
t194 = -t490 * t558 + t396;
t498 = qJD(1) * t132 + qJD(3) * t194;
t496 = -t626 / 0.2e1 - t631 / 0.2e1;
t495 = (t277 + t278) * t704;
t494 = t196 / 0.4e1 + t198 / 0.4e1 - t232 / 0.4e1 - t233 / 0.4e1;
t493 = t200 / 0.4e1 + t202 / 0.4e1 + t230 / 0.4e1 + t231 / 0.4e1;
t39 = m(7) * (t371 * t72 - t373 * t79) - t373 * t305 + t371 * t309;
t492 = t39 * qJD(1) + qJD(2) * t568;
t465 = (-t219 * t91 + t263 * t370 + t297 * t339 + t773 * t789 + t652) * t747 + pkin(3) * t722 + t279 * t734 + t228 * t729 + t811;
t485 = t299 * t808 + t578 * t786;
t1 = -t807 - t789 * t681 / 0.2e1 - t786 * t632 / 0.2e1 + t591 / 0.4e1 + t758 * pkin(8) - (t198 + t196) * t490 / 0.4e1 + (t233 + t232) * t490 / 0.4e1 - (t426 + t772) * t594 / 0.4e1 + t773 * t543 + (-t760 / 0.4e1 - t545) * t373 + (-t776 / 0.4e1 + t777 / 0.4e1 + t767) * t371 + t91 * t551 + t72 * t552 + t460 * t528 - t598 / 0.4e1 + t766 * t457 - (t231 + t230 + t202 + t200) * t768 / 0.4e1 + t465 + t521 * t706 + t588 * t712 + t588 * t713 + t459 * t721 + t572 * t735 - t635 * t739 - t620 * t742 + ((t445 * t588 + t604) * pkin(4) + t485) * t749;
t16 = (t772 / 0.2e1 + t426 / 0.2e1) * t459 + (-t424 / 0.2e1 + t427 / 0.2e1 + pkin(4) * t280) * t457 - pkin(3) * t423 - t511 + m(6) * t445 * t694 - t815 * t490 + t776 * t717 + t777 * t715 + t760 * t719 + t816 * t370;
t471 = -t370 * t460 * t748 + t441 * t749;
t477 = (-t322 - t520) * t749 - t764 * t747 + t510;
t18 = (t736 + t737 + t423 / 0.2e1 + t496) * t460 + t471 + t477;
t489 = t1 * qJD(1) - t18 * qJD(2) + t16 * qJD(3);
t17 = t773 * t682 - (t287 / 0.2e1 - t284 / 0.2e1 + t285 / 0.2e1 - t282 / 0.2e1 + t816 * pkin(5) - t815) * t490 - (t657 - t288 / 0.2e1 - t283 / 0.2e1 - t286 / 0.2e1 - t281 / 0.2e1) * t768 + t511;
t478 = -t490 * t689 * t747 + t495;
t21 = t371 * t562 - t373 * t499 + t478;
t480 = (t288 / 0.4e1 + t286 / 0.4e1 + t283 / 0.4e1 + t281 / 0.4e1 + t545 - t657 / 0.2e1) * t373;
t481 = -t287 / 0.4e1 - t285 / 0.4e1 + t284 / 0.4e1 + t282 / 0.4e1 + t767;
t463 = -((t78 / 0.2e1 + t741) * mrSges(7,3) + t493) * t768 - ((-t228 / 0.2e1 - t698 / 0.2e1) * pkin(5) + t494) * t490 - t480 + ((-t697 / 0.2e1 - t279 / 0.2e1) * pkin(5) + t481) * t371 + t811;
t4 = (t652 / 0.2e1 - t651 / 0.2e1 + t76 * t746) * m(7) - t692 / 0.2e1 + t463 + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1) * t460 - t754 + t757;
t488 = t4 * qJD(1) + t21 * qJD(2) - t17 * qJD(3);
t103 = (t609 / 0.2e1 + t606 / 0.2e1 - t458 / 0.2e1) * m(7);
t469 = (t371 * t715 - t373 * t719) * mrSges(7,3) + (t219 * t373 + t371 * t773 - t490 * t72 - t768 * t79) * t747 + t305 * t719 + t309 * t715;
t486 = t262 * t747 - t643 / 0.2e1 + t638 / 0.2e1;
t28 = t469 - t486;
t53 = -m(7) * (t219 * t768 - t490 * t773) + (-t490 ^ 2 - t768 ^ 2) * mrSges(7,3);
t487 = qJD(1) * t28 + qJD(2) * t103 - qJD(3) * t53;
t464 = (-t444 * t79 + (t701 * t779 + t702 * t79) * pkin(4)) * t747 - t666 / 0.2e1 - t665 / 0.2e1 - t687 / 0.2e1 - t686 / 0.2e1 + t444 * t543 + t506 / 0.2e1 - t774 * t570 / 0.2e1 + t775 * t514 + t790 * t371 * t517;
t476 = -t662 / 0.2e1 + t661 / 0.2e1 - t684 / 0.2e1 + t683 / 0.2e1 + t90 * t574 + mrSges(7,3) * t565;
t13 = t464 + t476;
t212 = (m(7) * t701 * t765 + t756) * pkin(4);
t468 = -t765 * t748 * t219 + t444 * t551 - t515 * t682 - t813;
t23 = t468 + t755;
t57 = (t746 + t515 + t711) * t696;
t484 = t13 * qJD(1) + t57 * qJD(2) - t23 * qJD(3) + t212 * qJD(4);
t321 = t372 * t570;
t120 = (-t570 * t768 - t599) * t747 + m(7) * t729;
t102 = (t606 + t609) * t747 + m(7) * t707;
t96 = qJD(6) * t568;
t93 = (-t373 * t570 + t602) * t747 + m(7) * t734;
t43 = -t765 * t696 / 0.2e1 - t374 * t574 + t497;
t27 = t469 + t486;
t22 = m(7) * t565 + t478 + t510;
t20 = -t468 + t508 + t755;
t19 = t423 * t704 + t460 * t496 - t471 + t477 + t495;
t14 = -t573 * t768 + t467 + t759;
t11 = (t722 + t525) * t460 + (t460 * t507 + t770) * t458 + t473 + t475 + t794;
t10 = t464 - t476 + t509;
t8 = t462 + t479;
t3 = t692 / 0.2e1 + t76 * t573 + t753 + t463 + (-t651 + t652) * t747;
t2 = -t480 - ((t90 / 0.2e1 + t79 / 0.2e1) * mrSges(7,3) + (t117 / 0.2e1 + t107 / 0.2e1) * mrSges(6,3) + t494) * t490 + t465 + t453 * t706 + (pkin(4) * t604 + t485) * t749 + (t528 + (-t426 / 0.4e1 - t772 / 0.4e1) * t457 + pkin(8) * t507 + (t712 + t713 + (m(6) * t710 + t735) * pkin(4)) * t459) * t460 - ((t91 / 0.2e1 + t741) * mrSges(7,3) + (t739 - t106 / 0.2e1) * mrSges(6,3) + t493) * t768 + (t721 + t369 / 0.4e1 + t418 * t745) * t459 + (-t367 / 0.4e1 - t628 / 0.4e1 + t416 * t745 + t766) * t457 + t481 * t371 + t807;
t25 = [qJD(2) * t24 + qJD(3) * t5 + qJD(4) * t6 + qJD(5) * t7 + qJD(6) * t39, t616 + 0.2e1 * (t747 + t749) * (t372 * t768 - t374 * t490) * qJD(2) + t8 * qJD(3) + t11 * qJD(4) + t14 * qJD(5) + t96, t625 + t8 * qJD(2) + t2 * qJD(4) + t3 * qJD(5) + t27 * qJD(6) + (-Ifges(4,5) - t459 * t426 / 0.2e1 + t424 * t708 + t622 * t461) * t576 + ((-pkin(3) * t451 + pkin(8) * t500) * t750 + (t110 * t578 - t111 * t299 + t406 * t445) * t749 + (-t219 * t86 + t262 * t339 + t76 * t773) * t747) * t751 + (-t111 * t635 - t110 * t632 - mrSges(4,2) * t587 + t578 * t308 - Ifges(4,6) * t460 + t445 * t227 + t406 * t280 + pkin(3) * t390 + t339 * t226 - t219 * t303 - t299 * t304 + t773 * t307 + t262 * t279 + (Ifges(5,5) * t457 + Ifges(5,6) * t459 + t490 * t810 - t768 * t798) * t703 - t76 * t681 - t86 * t682 + t366 * t705 + t368 * t708 + t777 * t802 + t785 * t723 + (t197 + t195) * t719 + (t201 + t199) * t717 + (t590 - t596) * pkin(8) + t500 * mrSges(5,3)) * qJD(3), t624 + t11 * qJD(2) + t2 * qJD(3) + (-Ifges(5,5) * t594 - Ifges(5,6) * t588 + m(7) * (t444 * t90 + t570 * t91) - t683 + t684 - t661 + t662 + (t117 * t702 + t118 * t701) * t744 - t319 * mrSges(5,2) - t320 * mrSges(5,1) + t506 + mrSges(6,3) * t322 + t509 + t764 * mrSges(7,3)) * qJD(4) + t10 * qJD(5) + t93 * qJD(6), t623 + t14 * qJD(2) + t3 * qJD(3) + t10 * qJD(4) + (-t665 - t686 - t666 - t687 + (-m(7) * t79 + t640) * pkin(5) + t509) * qJD(5), t27 * qJD(3) + t93 * qJD(4) + t492; qJD(3) * t9 + qJD(4) * t12 + qJD(5) * t15 - t616 + t96, qJD(3) * t45, t19 * qJD(4) + t22 * qJD(5) + t102 * qJD(6) + (t279 + t280 + t622) * t576 + ((t299 * t371 - t373 * t578 + t445 * t458) * t749 + (t219 * t371 + t339 * t458 - t373 * t773) * t747 + (t577 * t691 - t695) * t750) * t751 + t502 + ((-mrSges(4,2) + t797) * t460 + t790 * (t371 * t768 + t373 * t490)) * qJD(3), t619 + t19 * qJD(3) + (-mrSges(5,1) * t593 + mrSges(5,2) * t595 + t497) * qJD(4) + t43 * qJD(5) + 0.2e1 * ((t374 * t569 - t321) * t749 + (t444 * t374 - t321) * t747) * qJD(4), t618 + t22 * qJD(3) + t43 * qJD(4) + (pkin(5) * t696 + t497) * qJD(5), qJD(1) * t568 + t102 * qJD(3); -qJD(2) * t9 + qJD(4) * t1 + qJD(5) * t4 + qJD(6) * t28 - t625, -qJD(4) * t18 + qJD(5) * t21 + qJD(6) * t103 - t502, qJD(4) * t16 - qJD(5) * t17 - qJD(6) * t53 ((t299 * t702 + t578 * t701) * t744 + m(7) * (t444 * t219 + t570 * t773) + t521 + t762 * pkin(8) + t771 * mrSges(7,3) + (pkin(4) * t540 - t519) * mrSges(6,3) + t814) * qJD(4) + t20 * qJD(5) + t120 * qJD(6) + t489, t20 * qJD(4) + ((m(7) * t219 + t682) * pkin(5) + t814) * qJD(5) + t488, qJD(4) * t120 + t487; -qJD(2) * t12 - qJD(3) * t1 + qJD(5) * t13 + qJD(6) * t56 - t624, qJD(3) * t18 + qJD(5) * t57 - t619, -qJD(5) * t23 + qJD(6) * t75 - t489, t212 * qJD(5) (-t701 * t743 + t756) * qJD(5) * pkin(4) + t484, t501; -qJD(2) * t15 - qJD(3) * t4 - qJD(4) * t13 + qJD(6) * t132 - t623, -qJD(3) * t21 - qJD(4) * t57 - t618, qJD(4) * t23 + qJD(6) * t194 - t488, -t484, 0, t498; -t28 * qJD(3) - t56 * qJD(4) - t132 * qJD(5) - t492, -qJD(1) * t700 / 0.2e1 - t103 * qJD(3), -qJD(4) * t75 - qJD(5) * t194 - t487, -t501, -t498, 0;];
Cq  = t25;
