% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:30:51
% EndTime: 2019-03-09 16:31:15
% DurationCPUTime: 14.88s
% Computational Cost: add. (33262->584), mult. (64540->751), div. (0->0), fcn. (75182->8), ass. (0->340)
t403 = sin(qJ(5));
t650 = Ifges(7,4) + Ifges(6,4);
t776 = t650 * t403;
t406 = cos(qJ(5));
t779 = t650 * t406;
t775 = Ifges(6,1) + Ifges(7,1);
t774 = Ifges(6,2) + Ifges(7,2);
t404 = sin(qJ(3));
t405 = sin(qJ(2));
t407 = cos(qJ(3));
t408 = cos(qJ(2));
t379 = -t404 * t405 + t407 * t408;
t380 = -t404 * t408 - t407 * t405;
t402 = sin(pkin(10));
t617 = cos(pkin(10));
t314 = t617 * t379 + t380 * t402;
t601 = t314 * t406;
t530 = t601 / 0.2e1;
t580 = t403 * t314;
t682 = m(7) / 0.2e1;
t672 = -pkin(8) - pkin(7);
t387 = t672 * t405;
t388 = t672 * t408;
t332 = t387 * t404 - t388 * t407;
t284 = qJ(4) * t379 + t332;
t697 = t407 * t387 + t404 * t388;
t716 = t380 * qJ(4) + t697;
t730 = t617 * t284 + t402 * t716;
t743 = pkin(5) * t580 + t730;
t778 = mrSges(7,1) * t580 / 0.2e1 + mrSges(7,2) * t530 + t743 * t682;
t463 = t402 * t379 - t380 * t617;
t703 = t463 * mrSges(5,3);
t600 = t463 * t403;
t553 = mrSges(7,3) * t600;
t229 = mrSges(7,2) * t314 - t553;
t599 = t463 * t406;
t552 = mrSges(7,3) * t599;
t233 = -mrSges(7,1) * t314 - t552;
t658 = t403 / 0.2e1;
t739 = -t406 / 0.2e1;
t474 = t229 * t739 + t233 * t658;
t773 = t474 + t778;
t745 = -t774 * t403 + t779;
t744 = t775 * t406 - t776;
t772 = -t474 + t778;
t657 = t406 / 0.2e1;
t668 = t463 / 0.2e1;
t659 = -t403 / 0.2e1;
t734 = t775 * t403 + t779;
t735 = t774 * t406 + t776;
t713 = t657 * t734 + t659 * t735;
t649 = Ifges(6,5) + Ifges(7,5);
t705 = Ifges(7,6) + Ifges(6,6);
t718 = t403 * t649 + t705 * t406;
t720 = t744 * t314 + t463 * t649;
t721 = t314 * t745 + t463 * t705;
t423 = Ifges(4,5) * t379 + Ifges(4,6) * t380 - Ifges(5,6) * t463 + t718 * t668 + t720 * t658 + t721 * t657 + (Ifges(5,5) + t713) * t314;
t722 = t697 * mrSges(4,2);
t723 = t332 * mrSges(4,1);
t499 = mrSges(6,1) * t406 - mrSges(6,2) * t403;
t750 = t730 * t499;
t758 = t730 * mrSges(5,1);
t729 = -t402 * t284 + t617 * t716;
t759 = t729 * mrSges(5,2);
t384 = -mrSges(7,1) * t406 + mrSges(7,2) * t403;
t764 = t743 * t384;
t771 = t423 - t722 - t723 - t750 - t758 - t759 + t764;
t770 = -t722 / 0.2e1 - t723 / 0.2e1 - t750 / 0.2e1 - t758 / 0.2e1 - t759 / 0.2e1 + t764 / 0.2e1;
t736 = -t403 * t705 + t406 * t649;
t395 = -pkin(2) * t408 - pkin(1);
t344 = -t379 * pkin(3) + t395;
t769 = m(5) * t344;
t400 = t403 ^ 2;
t401 = t406 ^ 2;
t568 = t400 + t401;
t706 = mrSges(7,3) + mrSges(6,3);
t768 = t568 * t706 - mrSges(5,2);
t122 = pkin(5) * t600 - t729;
t767 = t122 * t743;
t656 = pkin(2) * t404;
t389 = t402 * t656;
t655 = pkin(2) * t407;
t394 = pkin(3) + t655;
t353 = t394 * t617 - t389;
t349 = -pkin(4) - t353;
t651 = t406 * pkin(5);
t334 = t349 - t651;
t766 = t334 * t743;
t532 = t617 * pkin(3);
t391 = -t532 - pkin(4);
t383 = t391 - t651;
t765 = t383 * t743;
t622 = t406 * mrSges(7,2);
t497 = t403 * mrSges(7,1) + t622;
t222 = t497 * t314;
t623 = t406 * mrSges(6,2);
t498 = t403 * mrSges(6,1) + t623;
t223 = t498 * t314;
t224 = t497 * t463;
t225 = t498 * t463;
t227 = -mrSges(7,2) * t463 - mrSges(7,3) * t580;
t228 = -mrSges(6,2) * t463 - mrSges(6,3) * t580;
t231 = mrSges(7,1) * t463 - mrSges(7,3) * t601;
t232 = mrSges(6,1) * t463 - mrSges(6,3) * t601;
t307 = t314 * mrSges(5,2);
t235 = mrSges(5,1) * t463 + t307;
t186 = -pkin(4) * t314 - pkin(9) * t463 + t344;
t84 = t406 * t186 - t403 * t730;
t462 = -qJ(6) * t599 + t84;
t56 = -pkin(5) * t314 + t462;
t602 = t730 * t406;
t67 = t602 + (-qJ(6) * t463 + t186) * t403;
t704 = Ifges(7,3) + Ifges(6,3);
t724 = t314 * mrSges(5,3);
t85 = t186 * t403 + t602;
t762 = t743 * t224 + t122 * t222 + t67 * t227 + t85 * t228 + t56 * t231 + t84 * t232 + t344 * t235 + t395 * (-mrSges(4,1) * t380 + mrSges(4,2) * t379) - (t223 + t724) * t729 + ((Ifges(5,1) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t401 + ((Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t403 - t779) * t403) * t314 - Ifges(5,4) * t463 + t736 * t668 + t721 * t659 + t720 * t657) * t463 - (-mrSges(5,3) * t729 + (-Ifges(5,4) + t736) * t314 + (Ifges(5,2) + t704) * t463) * t314 + t730 * t225;
t757 = t349 * t730;
t756 = t391 * t730;
t755 = t403 * t729;
t754 = t406 * t729;
t749 = t730 * t729;
t557 = mrSges(6,3) * t600;
t230 = mrSges(6,2) * t314 - t557;
t578 = t406 * t230;
t556 = mrSges(6,3) * t599;
t234 = -mrSges(6,1) * t314 - t556;
t581 = t403 * t234;
t442 = t462 - t56;
t710 = m(7) * t442;
t748 = -t581 / 0.2e1 + t578 / 0.2e1 + t658 * t710 - t474;
t747 = t745 * t463;
t746 = t744 * t463;
t517 = t617 * t404;
t354 = pkin(2) * t517 + t402 * t394;
t742 = -t353 * t730 + t354 * t729;
t741 = t402 * t729 - t617 * t730;
t738 = -t705 * t314 + t747;
t737 = -t314 * t649 + t746;
t731 = t706 * (t400 / 0.2e1 + t401 / 0.2e1);
t725 = t599 / 0.2e1;
t529 = t600 / 0.2e1;
t645 = mrSges(7,3) * t406;
t538 = -t645 / 0.2e1;
t506 = (t227 + t228) * t658 + (t231 + t232) * t657;
t502 = pkin(5) * t463 - qJ(6) * t601;
t654 = pkin(3) * t380;
t198 = pkin(4) * t463 - pkin(9) * t314 - t654;
t90 = t406 * t198 - t755;
t60 = t502 + t90;
t684 = m(6) / 0.2e1;
t687 = -m(5) / 0.2e1;
t544 = qJ(6) * t580;
t91 = t403 * t198 + t754;
t71 = -t544 + t91;
t719 = (t403 * t91 + t406 * t90) * t684 + (t403 * t71 + t406 * t60) * t682 + t654 * t687 + t506;
t715 = (-mrSges(4,1) * t404 - mrSges(4,2) * t407) * pkin(2);
t712 = -t463 / 0.2e1;
t503 = mrSges(7,1) * t599 - mrSges(7,2) * t600;
t711 = t503 / 0.2e1;
t653 = pkin(3) * t402;
t390 = pkin(9) + t653;
t573 = qJ(6) + t390;
t364 = t573 * t406;
t351 = t364 * t406;
t363 = t573 * t403;
t477 = t363 * t403 + t351;
t709 = m(7) * t477;
t708 = pkin(5) * t682;
t707 = -mrSges(6,2) - mrSges(7,2);
t680 = m(7) * pkin(5);
t545 = mrSges(7,1) + t680;
t350 = pkin(9) + t354;
t574 = qJ(6) + t350;
t328 = t574 * t406;
t319 = t328 * t406;
t327 = t574 * t403;
t478 = t327 * t403 + t319;
t520 = -t499 / 0.2e1 + t384 / 0.2e1;
t695 = mrSges(7,3) * t529 + t229 / 0.2e1;
t694 = -mrSges(5,1) + t384 - t499;
t693 = t224 + t225 + t703;
t692 = mrSges(6,3) * t568 * t712 + t230 * t659;
t561 = pkin(5) * t599;
t134 = -m(7) * t561 - t503;
t365 = t403 * t545 + t622;
t691 = -qJD(1) * t134 + (qJD(2) + qJD(3)) * t365;
t690 = t442 * t682 + t463 * t538 - t233 / 0.2e1;
t689 = t462 * t538 + t729 * t498 / 0.2e1 - t122 * t497 / 0.2e1 - t384 * t561 / 0.2e1 + t736 * t314 / 0.4e1 + t735 * t725 + t734 * t529 + (t738 + t747) * t403 / 0.4e1 - (t737 + t746) * t406 / 0.4e1;
t686 = m(5) / 0.2e1;
t685 = -m(6) / 0.2e1;
t683 = -m(7) / 0.2e1;
t681 = m(5) * pkin(3);
t679 = -mrSges(6,1) / 0.2e1;
t678 = mrSges(6,1) / 0.2e1;
t677 = -mrSges(7,1) / 0.2e1;
t676 = -mrSges(6,2) / 0.2e1;
t675 = -mrSges(7,2) / 0.2e1;
t399 = t405 * pkin(2);
t187 = t198 + t399;
t86 = t406 * t187 - t755;
t57 = t502 + t86;
t674 = t57 / 0.2e1;
t673 = t60 / 0.2e1;
t670 = t231 / 0.2e1;
t362 = t617 * t655 - t389;
t665 = -t362 / 0.2e1;
t652 = pkin(5) * t403;
t644 = Ifges(4,4) * t380;
t236 = -mrSges(5,1) * t314 + mrSges(5,2) * t463;
t348 = t399 - t654;
t501 = Ifges(4,4) * t379 + (-Ifges(4,1) + Ifges(4,2)) * t380;
t87 = t403 * t187 + t754;
t69 = -t544 + t87;
t2 = (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t408 + (Ifges(3,1) - Ifges(3,2)) * t405) * t408 + (-mrSges(4,1) * t399 + t501) * t379 + (-mrSges(4,2) * t399 - t644) * t380 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t405) * t405 + m(4) * t395 * t399 + m(6) * (t84 * t86 + t85 * t87 - t749) + m(7) * (t56 * t57 + t67 * t69 + t767) + t69 * t229 + t87 * t230 + t57 * t233 + t86 * t234 + t762 + (t769 + t236) * t348;
t635 = t2 * qJD(1);
t4 = -t654 * t769 + m(6) * (t84 * t90 + t85 * t91 - t749) + m(7) * (t56 * t60 + t67 * t71 + t767) + t501 * t379 + (-pkin(3) * t236 - t644) * t380 + t71 * t229 + t91 * t230 + t60 * t233 + t90 * t234 + t762;
t626 = t4 * qJD(1);
t625 = t403 * mrSges(6,3);
t624 = t403 * mrSges(7,3);
t621 = t87 * t406;
t469 = t499 * t463;
t9 = t729 * t469 - t224 * t561 - t462 * t229 - t56 * t553 + (t556 + t234) * t85 + (-t230 - t557) * t84 + t134 * t122 + (t233 + t552 - t710) * t67 + t737 * t529 + t738 * t725 + t718 * t314 * t712 + t713 * t463 ^ 2;
t620 = t9 * qJD(1);
t619 = t90 * t403;
t618 = t91 * t406;
t483 = -t403 * t84 + t406 * t85;
t484 = -t403 * t56 + t406 * t67;
t604 = t729 * t463;
t10 = t693 * t463 + (t724 + (t229 + t230) * t406 + (-t233 - t234) * t403) * t314 + m(6) * (t314 * t483 - t604) + m(7) * (t122 * t463 + t314 * t484) + m(5) * (t314 * t730 - t604);
t616 = qJD(1) * t10;
t516 = t568 * t350;
t417 = (t314 * t354 - t353 * t463) * t686 + (t314 * t516 + t349 * t463) * t684 + (t314 * t478 + t334 * t463) * t682 + t314 * t731;
t421 = t348 * t687 + (t403 * t87 + t406 * t86) * t685 + (t403 * t69 + t406 * t57) * t683;
t11 = -t307 + (-t232 / 0.2e1 - t231 / 0.2e1) * t406 + (-t228 / 0.2e1 - t227 / 0.2e1) * t403 + (-mrSges(5,1) + t520) * t463 + t417 + t421;
t615 = qJD(1) * t11;
t515 = t568 * t390;
t562 = t681 / 0.2e1;
t414 = (t383 * t682 + t391 * t684 - t562 * t617 + t520) * t463 + (t402 * t562 + t477 * t682 + t515 * t684 + t731) * t314;
t13 = t414 - t235 - t719;
t614 = qJD(1) * t13;
t27 = (m(7) * (-t403 * t67 - t56 * t406) - t403 * t229 - t406 * t233) * t463;
t613 = qJD(1) * t27;
t610 = t122 * t403;
t523 = -t580 / 0.2e1;
t531 = -t601 / 0.2e1;
t570 = mrSges(7,1) * t523 + mrSges(7,2) * t531;
t14 = mrSges(6,2) * t531 + t570 + (mrSges(6,1) + t680) * t523 - t748;
t609 = t14 * qJD(1);
t361 = (t402 * t407 + t517) * pkin(2);
t603 = t729 * t361;
t598 = t327 * t231;
t596 = t328 * t227;
t595 = t334 * t222;
t594 = t349 * t223;
t591 = t363 * t231;
t589 = t364 * t227;
t588 = t383 * t222;
t587 = t391 * t223;
t583 = t403 * t232;
t575 = t406 * t390;
t569 = t568 * mrSges(7,3);
t566 = qJD(5) * t403;
t565 = qJD(5) * t406;
t504 = m(7) * t568 * t463;
t549 = m(7) * t668;
t124 = t504 / 0.2e1 + t549;
t564 = t124 * qJD(1);
t560 = mrSges(6,3) * t618;
t559 = m(7) * t674;
t558 = m(7) * t673;
t555 = t353 * t724;
t554 = t354 * t703;
t548 = t361 * t682;
t547 = t652 / 0.2e1;
t546 = t675 + t676;
t543 = t390 * t583;
t542 = t228 * t575;
t537 = t645 / 0.2e1;
t536 = -t625 / 0.2e1;
t535 = t625 / 0.2e1;
t534 = -t624 / 0.2e1;
t533 = t624 / 0.2e1;
t521 = t350 * t739;
t514 = t703 * t653;
t512 = mrSges(7,3) * pkin(5) - t649;
t507 = t679 - t680 / 0.2e1;
t378 = t384 * t652;
t500 = -t401 * t650 - t378;
t482 = -t86 * t403 + t621;
t481 = t618 - t619;
t480 = t532 * t724;
t410 = t578 * t665 - t594 / 0.2e1 + t555 / 0.2e1 + t554 / 0.2e1 + (t483 * t685 + t484 * t683 + t687 * t730 + t474) * t362 + t60 * t533 + t90 * t535 + t71 * t538 + t228 * t521 + (t122 * t361 - t327 * t60 + t328 * t71 + t766) * t683 - t560 / 0.2e1 + (t350 * t481 - t603 + t757) * t685 + (-t603 + t742) * t687 + t598 / 0.2e1 - t596 / 0.2e1 - t595 / 0.2e1 + (-t724 + t581) * t362 / 0.2e1 - t693 * t361 / 0.2e1 + t350 * t583 / 0.2e1 - t770;
t411 = -t514 / 0.2e1 + t542 / 0.2e1 + t589 / 0.2e1 - t591 / 0.2e1 + t57 * t534 + t86 * t536 + t69 * t537 - t480 / 0.2e1 + (-t363 * t57 + t364 * t69 + t765) * t682 + (t390 * t482 + t756) * t684 - t543 / 0.2e1 + t741 * t562 + t587 / 0.2e1 + t588 / 0.2e1 + mrSges(6,3) * t621 / 0.2e1 + t770;
t3 = t410 + t411;
t35 = t715 + t694 * t361 + t768 * t362 + m(7) * (t334 * t361 + t362 * t478) + m(6) * (t349 * t361 + t362 * t516) + m(5) * (-t353 * t361 + t354 * t362);
t479 = -t3 * qJD(1) + t35 * qJD(2);
t425 = m(7) * ((t327 * t406 - t328 * t403) * t463 + t484);
t22 = -t425 / 0.2e1 + t773;
t226 = m(7) * t478 + t569;
t476 = -qJD(1) * t22 + qJD(2) * t226;
t318 = t334 * t497;
t326 = t349 * t498;
t343 = t383 * t497;
t360 = t391 * t498;
t461 = -t318 / 0.2e1 - t326 / 0.2e1 - t343 / 0.2e1 - t360 / 0.2e1 - t378;
t444 = t469 / 0.2e1;
t433 = (t319 + t351 + (t327 + t363) * t403) * t683 - t569;
t118 = t548 + t433;
t424 = m(7) * ((t363 * t406 - t364 * t403) * t463 + t484);
t24 = -t424 / 0.2e1 + t773;
t250 = t569 + t709;
t443 = -qJD(1) * t24 - qJD(2) * t118 + qJD(3) * t250;
t435 = t776 + (t774 - t775) * t406;
t41 = -t318 - t326 + (-t334 * t680 + t435) * t403 + t500;
t409 = t704 * t668 + t649 * t530 + t705 * t523 + t56 * t537 + t689 + pkin(5) * t670 - t224 * t652 / 0.2e1;
t413 = (t334 * t599 + t610) * t708 + t334 * t711 + t349 * t444 + t234 * t521 + t690 * t328 - t695 * t327 + t692 * t350;
t428 = mrSges(7,1) * t674 + t675 * t69 + t676 * t87 + t678 * t86;
t5 = pkin(5) * t559 + t409 - t413 + t428;
t440 = t5 * qJD(1) + t41 * qJD(2);
t427 = mrSges(7,1) * t673 + t675 * t71 + t676 * t91 + t678 * t90;
t30 = (t362 * t546 - t779) * t406 + ((t677 + t679) * t362 + (t665 - t334 / 0.2e1 - t383 / 0.2e1) * t680 + t435) * t403 + t461;
t48 = -t343 - t360 + (-t383 * t680 + t435) * t403 + t500;
t412 = (t383 * t599 + t610) * t708 + t383 * t711 + t391 * t444 - t234 * t575 / 0.2e1 + t690 * t364 - t695 * t363 + t692 * t390;
t7 = pkin(5) * t558 + t409 - t412 + t427;
t426 = t7 * qJD(1) + t30 * qJD(2) + t48 * qJD(3);
t416 = (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t463 + ((Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1) * t406 + (-Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t403) * t314 + t56 * t538 - t689 + t224 * t547 + (t535 + t536) * t85 + (t533 + t534) * t67;
t356 = t365 * qJD(6);
t355 = t365 * qJD(5);
t125 = t549 - t504 / 0.2e1;
t119 = t548 - t433;
t31 = m(7) * (t334 + t383) * t547 + (t546 * t406 + (t677 + t507) * t403) * t362 - t461 + t744 * t658 + t745 * t657 + t713;
t25 = t424 / 0.2e1 + t772;
t23 = t425 / 0.2e1 + t772;
t16 = t414 + t719;
t15 = (-t623 / 0.2e1 + t507 * t403) * t314 + t570 + t748;
t12 = t463 * t520 + t417 - t421 + t506;
t8 = (t558 + t670) * pkin(5) + t412 + t416 + t427;
t6 = (t559 + t670) * pkin(5) + t413 + t416 + t428;
t1 = t423 - t410 + t411;
t17 = [qJD(2) * t2 + qJD(3) * t4 + qJD(4) * t10 - qJD(5) * t9 + qJD(6) * t27, t1 * qJD(3) + t12 * qJD(4) + t6 * qJD(5) + t23 * qJD(6) + t635 + (m(4) * (-t332 * t407 + t404 * t697) * pkin(2) + t594 - t555 - t554 + 0.2e1 * (-t327 * t57 + t328 * t69 + t766) * t682 + (-mrSges(3,1) * t408 + mrSges(3,2) * t405) * pkin(7) + 0.2e1 * (t350 * t482 + t757) * t684 + 0.2e1 * t742 * t686 - t598 + t596 + t595 + (-t379 * t655 + t380 * t656) * mrSges(4,3) + (t87 * mrSges(6,3) + t69 * mrSges(7,3) + t350 * t228) * t406 + (-t86 * mrSges(6,3) - t57 * mrSges(7,3) - t350 * t232) * t403 - Ifges(3,6) * t405 + Ifges(3,5) * t408 + t771) * qJD(2), t626 + t1 * qJD(2) + (-t514 + t542 + t71 * t645 + t589 - t591 - t480 + m(7) * (-t363 * t60 + t364 * t71 + t765) + t560 + m(6) * (t390 * t481 + t756) - t543 + t741 * t681 + t587 + t588 - mrSges(6,3) * t619 - t60 * t624 + t771) * qJD(3) + t16 * qJD(4) + t8 * qJD(5) + t25 * qJD(6), qJD(2) * t12 + qJD(3) * t16 + qJD(5) * t15 + qJD(6) * t125 + t616, t6 * qJD(2) + t8 * qJD(3) + t15 * qJD(4) - t620 + (-t85 * mrSges(6,1) + ((mrSges(7,2) * qJ(6) - t705) * t406 + t512 * t403) * t463 + t707 * t84 - t545 * t67) * qJD(5), qJD(2) * t23 + qJD(3) * t25 + qJD(4) * t125 + t613; -qJD(3) * t3 + qJD(4) * t11 - qJD(5) * t5 - qJD(6) * t22 - t635, qJD(3) * t35 - qJD(5) * t41 + qJD(6) * t226, t31 * qJD(5) + t119 * qJD(6) + t479 + ((m(6) * t391 + m(7) * t383 - t617 * t681 + t694) * t361 + t715 + (m(6) * t515 + t402 * t681 + t709 + t768) * t362) * qJD(3), t615, t31 * qJD(3) + (mrSges(7,2) * t327 - t328 * t545) * qJD(5) + (mrSges(6,2) * t350 - t705) * t566 + (-mrSges(6,1) * t350 - t512) * t565 - t440, qJD(3) * t119 + t476; qJD(2) * t3 + qJD(4) * t13 - qJD(5) * t7 - qJD(6) * t24 - t626, -qJD(5) * t30 - qJD(6) * t118 - t479, -qJD(5) * t48 + qJD(6) * t250, t614 (mrSges(7,2) * t363 - t364 * t545) * qJD(5) + (mrSges(6,2) * t390 - t705) * t566 + (-mrSges(6,1) * t390 - t512) * t565 - t426, t443; -qJD(2) * t11 - qJD(3) * t13 - qJD(5) * t14 - qJD(6) * t124 - t616, -t615, -t614, 0, -t609 + t707 * t565 + (-mrSges(6,1) - t545) * t566, -t564; qJD(2) * t5 + qJD(3) * t7 + qJD(4) * t14 + qJD(6) * t134 + t620, t30 * qJD(3) - t356 + t440, -t356 + t426, t609, 0, -t691; qJD(2) * t22 + qJD(3) * t24 + qJD(4) * t124 - qJD(5) * t134 - t613, qJD(3) * t118 + t355 - t476, t355 - t443, t564, t691, 0;];
Cq  = t17;
