% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:35
% EndTime: 2019-03-09 09:28:58
% DurationCPUTime: 12.90s
% Computational Cost: add. (17496->847), mult. (40952->1137), div. (0->0), fcn. (39834->8), ass. (0->432)
t461 = sin(qJ(6));
t453 = t461 ^ 2;
t464 = cos(qJ(6));
t454 = t464 ^ 2;
t585 = t453 + t454;
t741 = mrSges(7,3) * t585;
t455 = sin(pkin(6));
t463 = sin(qJ(2));
t462 = sin(qJ(5));
t466 = cos(qJ(2));
t602 = t462 * t466;
t316 = (-t461 * t602 - t463 * t464) * t455;
t699 = t316 / 0.2e1;
t592 = t464 * t466;
t317 = (-t461 * t463 + t462 * t592) * t455;
t698 = t317 / 0.2e1;
t540 = (t454 / 0.2e1 + t453 / 0.2e1) * mrSges(7,3);
t740 = mrSges(6,2) - t540;
t629 = t464 * mrSges(7,1);
t638 = t461 * mrSges(7,2);
t534 = t629 - t638;
t739 = t462 * t534;
t720 = m(7) * pkin(5);
t738 = -t720 - mrSges(6,1) - t534;
t456 = cos(pkin(6));
t465 = cos(qJ(5));
t615 = t455 * t463;
t374 = t456 * t465 + t462 * t615;
t261 = -t374 * t461 + t455 * t592;
t614 = t455 * t466;
t262 = t374 * t464 + t461 * t614;
t373 = t456 * t462 - t465 * t615;
t647 = t373 * mrSges(7,2);
t655 = t261 * mrSges(7,3);
t174 = -t647 + t655;
t648 = t373 * mrSges(7,1);
t653 = t262 * mrSges(7,3);
t175 = t648 - t653;
t679 = -t464 / 0.2e1;
t685 = -t461 / 0.2e1;
t504 = t174 * t685 + t175 * t679;
t683 = t461 / 0.2e1;
t737 = (t261 * t683 + t262 * t679) * mrSges(7,3) + t504;
t588 = t465 * t466;
t575 = t455 * t588;
t736 = -Ifges(7,3) * t575 / 0.2e1;
t642 = t374 * mrSges(6,3);
t134 = -mrSges(7,1) * t261 + mrSges(7,2) * t262;
t278 = mrSges(6,1) * t614 - t642;
t587 = -t134 + t278;
t451 = t456 * pkin(2);
t734 = t456 * qJ(4) + t451;
t452 = t465 * mrSges(6,1);
t411 = -t462 * mrSges(6,2) + t452;
t628 = t464 * mrSges(7,2);
t639 = t461 * mrSges(7,1);
t412 = t628 + t639;
t193 = t412 * t373;
t646 = t373 * mrSges(6,3);
t733 = t646 / 0.2e1 + t193 / 0.2e1;
t670 = pkin(5) * t465;
t423 = pkin(10) * t462 + t670;
t457 = qJ(3) - pkin(9);
t604 = t461 * t465;
t302 = t423 * t464 - t457 * t604;
t593 = t464 * t465;
t303 = t423 * t461 + t457 * t593;
t517 = -t302 * t461 + t303 * t464;
t634 = t462 * mrSges(7,1);
t399 = -mrSges(7,3) * t593 + t634;
t594 = t464 * t399;
t676 = -t465 / 0.2e1;
t732 = -t594 / 0.2e1 + t676 * t741;
t731 = Ifges(7,5) * t698 + Ifges(7,6) * t699;
t729 = t462 ^ 2;
t728 = t465 ^ 2;
t727 = m(4) / 0.2e1;
t726 = -m(5) / 0.2e1;
t725 = m(5) / 0.2e1;
t724 = -m(6) / 0.2e1;
t723 = m(6) / 0.2e1;
t722 = -m(7) / 0.2e1;
t721 = m(7) / 0.2e1;
t719 = -mrSges(6,1) / 0.2e1;
t718 = -mrSges(7,1) / 0.2e1;
t717 = -mrSges(7,2) / 0.2e1;
t716 = mrSges(7,2) / 0.2e1;
t450 = t456 * qJ(3);
t714 = pkin(3) + pkin(8);
t584 = pkin(4) + t714;
t536 = t584 * t614;
t674 = pkin(1) * t463;
t208 = t450 + (-pkin(9) + t674) * t456 + t536;
t458 = pkin(2) + qJ(4);
t539 = -t458 * t466 - pkin(1);
t226 = (-t457 * t463 + t539) * t455;
t96 = t208 * t465 - t226 * t462;
t84 = -pkin(5) * t614 - t96;
t715 = t84 / 0.2e1;
t712 = t174 / 0.2e1;
t179 = -mrSges(7,1) * t316 + mrSges(7,2) * t317;
t711 = t179 / 0.2e1;
t709 = t261 / 0.2e1;
t708 = t261 / 0.4e1;
t707 = t262 / 0.2e1;
t706 = t262 / 0.4e1;
t664 = Ifges(6,4) * t465;
t533 = Ifges(6,1) * t462 + t664;
t705 = (-Ifges(6,5) * t463 + t466 * t533) * t455 / 0.2e1;
t277 = -mrSges(6,2) * t614 - t646;
t704 = t277 / 0.2e1;
t703 = -t278 / 0.2e1;
t671 = pkin(5) * t462;
t407 = -pkin(10) * t465 + t458 + t671;
t605 = t461 * t462;
t289 = t407 * t464 - t457 * t605;
t702 = t289 / 0.2e1;
t603 = t462 * t464;
t290 = t407 * t461 + t457 * t603;
t701 = -t290 / 0.2e1;
t700 = t302 / 0.2e1;
t697 = t373 / 0.2e1;
t696 = -t374 / 0.2e1;
t382 = t465 * t412;
t695 = t382 / 0.2e1;
t580 = mrSges(7,3) * t604;
t632 = t462 * mrSges(7,2);
t397 = -t580 - t632;
t694 = t397 / 0.2e1;
t693 = -t399 / 0.2e1;
t692 = t399 / 0.2e1;
t691 = t534 / 0.2e1;
t690 = t412 / 0.2e1;
t414 = Ifges(7,5) * t461 + Ifges(7,6) * t464;
t689 = t414 / 0.2e1;
t663 = Ifges(7,4) * t461;
t415 = Ifges(7,2) * t464 + t663;
t688 = -t415 / 0.4e1;
t687 = t415 / 0.4e1;
t662 = Ifges(7,4) * t464;
t418 = Ifges(7,1) * t461 + t662;
t686 = -t418 / 0.4e1;
t684 = -t461 / 0.4e1;
t682 = t461 / 0.4e1;
t681 = t462 / 0.2e1;
t680 = t462 / 0.4e1;
t678 = t464 / 0.2e1;
t677 = t464 / 0.4e1;
t675 = t465 / 0.2e1;
t673 = pkin(5) * t374;
t380 = t534 * t465;
t672 = pkin(5) * t380;
t669 = pkin(10) * t373;
t668 = mrSges(5,2) + mrSges(4,3);
t667 = m(7) * qJD(4);
t666 = Ifges(6,4) * t374;
t665 = Ifges(6,4) * t462;
t249 = Ifges(7,4) * t261;
t660 = Ifges(7,5) * t464;
t658 = Ifges(7,6) * t461;
t657 = Ifges(7,3) * t374;
t656 = Ifges(7,3) * t465;
t654 = t261 * Ifges(7,2);
t652 = t262 * Ifges(7,1);
t651 = t262 * Ifges(7,4);
t645 = t373 * Ifges(7,5);
t100 = t249 + t645 + t652;
t434 = qJ(4) * t615;
t442 = pkin(2) * t615;
t254 = -t457 * t614 + t434 + t442;
t443 = t456 * t466 * pkin(1);
t259 = -t584 * t615 + t443;
t116 = -t254 * t462 + t259 * t465;
t103 = pkin(5) * t615 - t116;
t117 = t465 * t254 + t462 * t259;
t147 = Ifges(7,4) * t317 + Ifges(7,2) * t316 - Ifges(7,6) * t575;
t148 = Ifges(7,1) * t317 + Ifges(7,4) * t316 - Ifges(7,5) * t575;
t362 = Ifges(6,4) * t373;
t581 = Ifges(6,5) * t614;
t189 = Ifges(6,1) * t374 - t362 + t581;
t207 = t259 + t734;
t643 = t374 * mrSges(6,2);
t649 = t373 * mrSges(6,1);
t220 = t643 + t649;
t230 = mrSges(7,2) * t575 + mrSges(7,3) * t316;
t231 = -mrSges(7,1) * t575 - mrSges(7,3) * t317;
t577 = t455 * t714;
t245 = t463 * t577 - t443 - t734;
t583 = t456 * t674;
t260 = t536 + t583;
t623 = qJ(3) * t463;
t265 = (t539 - t623) * t455;
t299 = t466 * t577 + t583;
t272 = t450 + t299;
t376 = -qJ(3) * t614 + t442;
t295 = t376 + t434;
t298 = -t615 * t714 + t443;
t319 = t411 * t614;
t379 = pkin(8) * t614 + t583;
t320 = -t379 - t450;
t321 = (-pkin(2) * t466 - pkin(1) - t623) * t455;
t378 = pkin(8) * t615 - t443;
t325 = t378 - t451;
t353 = (-mrSges(6,1) * t463 - mrSges(6,3) * t602) * t455;
t354 = (mrSges(6,2) * t463 + mrSges(6,3) * t588) * t455;
t375 = (mrSges(4,2) * t466 - mrSges(4,3) * t463) * t455;
t377 = (-mrSges(5,2) * t463 - mrSges(5,3) * t466) * t455;
t391 = mrSges(5,1) * t615 - t456 * mrSges(5,3);
t582 = mrSges(4,1) * t614;
t392 = -t456 * mrSges(4,3) - t582;
t435 = mrSges(5,1) * t614;
t393 = t456 * mrSges(5,2) + t435;
t437 = Ifges(5,5) * t614;
t438 = Ifges(4,5) * t615;
t439 = Ifges(3,5) * t614;
t440 = Ifges(5,4) * t615;
t105 = pkin(5) * t373 - pkin(10) * t374 + t207;
t97 = t208 * t462 + t226 * t465;
t85 = pkin(10) * t614 + t97;
t54 = t105 * t464 - t461 * t85;
t55 = t105 * t461 + t464 * t85;
t531 = Ifges(6,2) * t465 + t665;
t553 = t736 - (-Ifges(6,6) * t463 + t466 * t531) * t455 / 0.2e1 + t731;
t188 = -Ifges(6,2) * t373 + Ifges(6,6) * t614 + t666;
t98 = Ifges(7,5) * t262 + Ifges(7,6) * t261 + Ifges(7,3) * t373;
t567 = t188 / 0.2e1 - t98 / 0.2e1;
t641 = t378 * mrSges(3,2);
t104 = -pkin(10) * t615 + t117;
t173 = -t583 + (-t423 - t584) * t614;
t69 = -t104 * t461 + t173 * t464;
t70 = t104 * t464 + t173 * t461;
t644 = t373 * Ifges(7,6);
t99 = t644 + t651 + t654;
t3 = t553 * t373 + m(5) * (t245 * t299 + t265 * t295 + t272 * t298) + m(6) * (t116 * t96 + t117 * t97 - t207 * t260) + m(7) * (t103 * t84 + t54 * t69 + t55 * t70) + m(4) * (t320 * t378 + t321 * t376 + t325 * t379) + (t439 / 0.2e1 + t438 / 0.2e1 + t641 + t440 / 0.2e1 + t437 / 0.2e1 + (-mrSges(3,1) + mrSges(4,2)) * t379) * t456 + ((-t321 * mrSges(4,2) + t265 * mrSges(5,3) + Ifges(6,5) * t696 + Ifges(6,6) * t697 + Ifges(5,6) * t615 - t272 * mrSges(5,1) + (-pkin(1) * mrSges(3,1) + (-Ifges(3,4) - Ifges(4,6)) * t463) * t455 + (t379 + t320) * mrSges(4,1) + (Ifges(5,4) / 0.2e1 - Ifges(3,6) + Ifges(4,5) / 0.2e1) * t456) * t463 + (Ifges(3,4) * t614 - t321 * mrSges(4,3) - t265 * mrSges(5,2) + t245 * mrSges(5,1) + t325 * mrSges(4,1) + t189 * t681 + t567 * t465 + (-pkin(1) * mrSges(3,2) + (Ifges(6,5) * t681 + Ifges(6,6) * t675 + Ifges(4,6) - Ifges(5,6)) * t466) * t455 + (Ifges(3,5) / 0.2e1 - Ifges(4,4) + Ifges(5,5) / 0.2e1) * t456 + (Ifges(3,1) - Ifges(3,2) + Ifges(4,2) - Ifges(5,2) - Ifges(4,3) + Ifges(5,3) - Ifges(6,3)) * t615) * t466) * t455 + t100 * t698 + t99 * t699 + t374 * t705 + t148 * t707 + t147 * t709 + t299 * t391 + t378 * t392 + t298 * t393 + t376 * t375 + t295 * t377 + t96 * t353 + t97 * t354 - t207 * t319 + t103 * t134 + t70 * t174 + t69 * t175 + t84 * t179 + t55 * t230 + t54 * t231 - t260 * t220 + t117 * t277 + t116 * t278;
t650 = t3 * qJD(1);
t529 = -t658 + t660;
t143 = -t373 * t529 + t657;
t530 = -Ifges(7,2) * t461 + t662;
t144 = Ifges(7,6) * t374 - t373 * t530;
t532 = Ifges(7,1) * t464 - t663;
t145 = Ifges(7,5) * t374 - t373 * t532;
t637 = t461 * mrSges(7,3);
t209 = -mrSges(7,2) * t374 + t373 * t637;
t627 = t464 * mrSges(7,3);
t210 = mrSges(7,1) * t374 + t373 * t627;
t219 = t374 * mrSges(6,1) - t373 * mrSges(6,2);
t221 = -Ifges(6,2) * t374 - t362;
t222 = -Ifges(6,1) * t373 - t666;
t564 = t614 / 0.2e1;
t586 = -Ifges(6,5) * t373 - Ifges(6,6) * t374;
t601 = t464 * t100;
t635 = t461 * t99;
t223 = t669 + t673;
t74 = t223 * t464 - t461 * t96;
t75 = t223 * t461 + t464 * t96;
t4 = t144 * t709 + m(7) * (t54 * t74 + t55 * t75) - t84 * t193 + t75 * t174 + t55 * t209 + t74 * t175 + t54 * t210 + t145 * t707 + t207 * t219 + t96 * t277 + t586 * t564 + (t222 / 0.2e1 - t567) * t374 + (t143 / 0.2e1 - t189 / 0.2e1 - t221 / 0.2e1 + t635 / 0.2e1 - t601 / 0.2e1 + t96 * mrSges(6,3)) * t373 + (m(7) * t84 - t587 - t642) * t97;
t640 = t4 * qJD(1);
t636 = t461 * t54;
t631 = t462 * Ifges(7,5);
t630 = t462 * Ifges(7,6);
t626 = t464 * t55;
t625 = t465 * t97;
t133 = mrSges(7,1) * t262 + mrSges(7,2) * t261;
t135 = Ifges(7,5) * t261 - Ifges(7,6) * t262;
t136 = -Ifges(7,2) * t262 + t249;
t137 = Ifges(7,1) * t261 - t651;
t7 = t84 * t133 + t54 * t174 - t55 * t175 + t135 * t697 + (-t99 / 0.2e1 - t55 * mrSges(7,3) + t137 / 0.2e1) * t262 + (t136 / 0.2e1 - t54 * mrSges(7,3) + t100 / 0.2e1) * t261;
t624 = t7 * qJD(1);
t528 = t626 - t636;
t600 = t464 * t174;
t610 = t461 * t175;
t16 = (-t392 + t393) * t456 + (-t375 - t377) * t615 + t587 * t374 + (t277 + t600 - t610) * t373 + m(7) * (t373 * t528 - t374 * t84) + m(6) * (t373 * t97 + t374 * t96) + m(5) * (-t265 * t615 + t272 * t456) + m(4) * (-t320 * t456 - t321 * t615);
t622 = qJD(1) * t16;
t351 = t456 * t464 + t461 * t575;
t352 = t456 * t461 - t464 * t575;
t576 = t455 * t602;
t591 = t465 * t277;
t17 = t352 * t174 + t351 * t175 + (t220 - t391) * t456 + (t587 * t462 - t377 - t591) * t614 + m(7) * (t351 * t54 + t352 * t55 - t576 * t84) + m(6) * (t207 * t456 + (t462 * t96 - t625) * t614) + m(5) * (-t245 * t456 - t265 * t614);
t621 = t17 * qJD(1);
t620 = t290 * t464;
t617 = t351 * t461;
t616 = t352 * t464;
t613 = t456 * t458;
t612 = t457 * t462;
t611 = t457 * t465;
t609 = t461 * t210;
t497 = t530 * t465;
t341 = t497 + t630;
t608 = t461 * t341;
t607 = t461 * t397;
t398 = mrSges(7,1) * t465 + mrSges(7,3) * t603;
t606 = t461 * t398;
t599 = t464 * t209;
t598 = t464 * t230;
t498 = t532 * t465;
t343 = t498 + t631;
t597 = t464 * t343;
t396 = -mrSges(7,2) * t465 + mrSges(7,3) * t605;
t596 = t464 * t396;
t595 = t464 * t397;
t590 = t465 * t374;
t589 = t465 * t380;
t579 = pkin(10) * t685;
t578 = Ifges(7,2) / 0.4e1 - Ifges(7,1) / 0.4e1;
t572 = t653 / 0.2e1;
t570 = -t642 / 0.2e1;
t569 = -t637 / 0.2e1;
t568 = t627 / 0.2e1;
t566 = -t614 / 0.2e1;
t565 = -t614 / 0.4e1;
t561 = -t605 / 0.2e1;
t560 = -t602 / 0.2e1;
t558 = t598 / 0.2e1;
t555 = t412 * t676;
t554 = t134 / 0.2e1 + t703;
t339 = Ifges(7,3) * t462 + t465 * t529;
t417 = -Ifges(6,2) * t462 + t664;
t552 = -t339 / 0.2e1 + t417 / 0.2e1;
t550 = t728 + t729;
t549 = t465 * t585;
t548 = t585 * t462;
t544 = t455 * t560;
t543 = t462 * t564;
t542 = t465 * t564;
t541 = m(6) * t550;
t538 = t712 - t655 / 0.2e1;
t537 = t572 + t175 / 0.2e1;
t535 = -t533 / 0.4e1 - t417 / 0.4e1 + t339 / 0.4e1;
t413 = t462 * mrSges(6,1) + t465 * mrSges(6,2);
t527 = -t69 * t461 + t70 * t464;
t526 = -t461 * t74 + t464 * t75;
t32 = t290 * t399 + (t457 * t380 + t343 * t683 + t341 * t678 + t414 * t681 + mrSges(7,3) * t620 + (t415 * t685 + t418 * t678) * t465) * t465 + (-t397 - t580) * t289;
t471 = -t373 * t414 / 0.4e1 + t262 * t686 + t261 * t688 + t137 * t677 - t457 * t133 / 0.2e1 - t464 * t99 / 0.4e1 + (-t626 / 0.2e1 + t636 / 0.2e1) * mrSges(7,3) + (t136 + t100) * t684;
t473 = (t262 * t701 - t289 * t261 / 0.2e1) * mrSges(7,3) + t343 * t708 - t262 * t341 / 0.4e1 + t174 * t702 + t175 * t701 + t135 * t680 + t54 * t694 + t55 * t693 + t380 * t715;
t485 = t69 * t718 + t70 * t716 - t731;
t5 = (Ifges(7,3) * t564 + t471) * t465 + t473 + t485;
t525 = t5 * qJD(1) - t32 * qJD(2);
t448 = 0.2e1 * t450;
t501 = t399 * t683 - t595 / 0.2e1;
t515 = t373 * t462 + t590;
t518 = -t289 * t461 + t620;
t468 = t668 * t456 - t501 * t373 + (t379 + t448) * t727 + (t448 + t299) * t725 + (t457 * t515 + t462 * t97 + t465 * t96) * t723 + ((t374 * t457 - t84) * t465 + t528 * t462 + t518 * t373) * t721 + t382 * t696;
t475 = -m(4) * t379 / 0.2e1 + t299 * t726 + t260 * t724 + (-t461 * t70 - t464 * t69) * t722 + t230 * t683 + t231 * t678;
t503 = t610 / 0.2e1 - t600 / 0.2e1;
t494 = t704 - t503;
t489 = -t646 / 0.2e1 + t494;
t505 = t570 - t554;
t11 = (mrSges(6,1) * t566 + t505) * t465 + (mrSges(6,2) * t564 + t489) * t462 + t468 + t475;
t57 = t399 * t605 + t465 * t382 - t462 * t595 + t550 * mrSges(6,3) - m(7) * (t457 * t728 + t462 * t518) - t457 * t541 - t668 + (-m(5) - m(4)) * qJ(3);
t524 = -qJD(1) * t11 + qJD(2) * t57;
t469 = (t413 / 0.2e1 + mrSges(5,3)) * t456 + (t298 + t613 + t734) * t725 + (t207 + t613) * t723 + (t289 * t351 + t290 * t352 + t461 * t55 + t464 * t54 + t575 * t612) * t721 + t351 * t692 + t352 * t694 + t649 / 0.2e1 + t643 / 0.2e1 - t504;
t519 = t116 * t465 + t117 * t462;
t474 = t298 * t726 + t519 * t724 + (-t103 * t465 + t462 * t527) * t722;
t15 = t469 + (t711 - t353 / 0.2e1) * t465 + (t382 * t566 - t598 / 0.2e1 + t231 * t683 - t354 / 0.2e1) * t462 + t474;
t76 = t607 + t594 + mrSges(5,3) + m(7) * (t289 * t464 + t290 * t461) + 0.4e1 * (m(6) / 0.4e1 + m(5) / 0.4e1) * t458 + t413;
t523 = -qJD(1) * t15 - qJD(2) * t76;
t478 = t133 * t676 + t462 * t737;
t510 = t351 * t718 + t352 * t716;
t22 = t478 + t510;
t509 = -t638 / 0.2e1 + t629 / 0.2e1;
t58 = t589 / 0.2e1 + (t607 / 0.2e1 + t594 / 0.2e1 + t465 * t540) * t462 + t509;
t522 = t22 * qJD(1) - t58 * qJD(2);
t23 = (-t647 / 0.2e1 + t538) * t464 + (-t648 / 0.2e1 - t537) * t461;
t77 = (-t632 / 0.2e1 + t694) * t464 + (-t634 / 0.2e1 + t693) * t461;
t521 = qJD(1) * t23 + qJD(2) * t77;
t482 = (-t461 * t75 - t464 * t74) * t721 + t209 * t685 + t210 * t679;
t491 = (t585 * t669 + t673) * t721;
t26 = (t691 + mrSges(6,1)) * t374 - t740 * t373 + t491 - t482;
t480 = (-t302 * t464 - t303 * t461) * t721 + t396 * t685 + t398 * t679;
t483 = (pkin(10) * t548 + t670) * t721 - t534 * t676;
t60 = t462 * t740 - t452 + t480 - t483;
t520 = qJD(1) * t26 - qJD(2) * t60;
t516 = t616 - t617;
t481 = t541 / 0.2e1 + (t585 * t729 + t728) * t721;
t513 = t585 * t722 + t724;
t180 = -m(5) - t481 + t513;
t479 = t515 * t723 + (t373 * t548 + t590) * t721;
t499 = m(7) * (-t351 * t464 - t352 * t461);
t53 = (m(5) + t723) * t456 - t499 / 0.2e1 + t479;
t514 = qJD(1) * t53 - qJD(2) * t180;
t512 = -t660 / 0.2e1 + t658 / 0.2e1;
t511 = mrSges(7,1) * t700 + t303 * t717;
t508 = -t628 / 0.2e1 - t639 / 0.2e1;
t507 = t635 / 0.4e1 - t601 / 0.4e1;
t506 = t570 + t554;
t502 = t608 / 0.4e1 - t597 / 0.4e1;
t500 = t415 * t683 + t418 * t679;
t496 = t516 * pkin(10);
t340 = Ifges(7,6) * t465 - t462 * t530;
t342 = Ifges(7,5) * t465 - t462 * t532;
t381 = t412 * t462;
t470 = (t289 * t74 + t290 * t75 + t302 * t54 + t303 * t55) * t721 + t207 * t411 / 0.2e1 + t340 * t708 + t342 * t706 + t210 * t702 + t290 * t209 / 0.2e1 + t175 * t700 + t303 * t712 + t458 * t219 / 0.2e1 + t54 * t398 / 0.2e1 + t55 * t396 / 0.2e1 + t74 * t692 + t75 * t694 - t381 * t715 + t97 * t695;
t472 = (-pkin(5) * t103 + pkin(10) * t527) * t722 + pkin(5) * t711 + t103 * t691 + t116 * t719 + t117 * mrSges(6,2) / 0.2e1 + t316 * t688 + t317 * t686 + Ifges(6,3) * t615 / 0.2e1;
t338 = -t462 * t529 + t656;
t420 = Ifges(6,1) * t465 - t665;
t486 = -t420 / 0.4e1 + t531 / 0.4e1 + t338 / 0.4e1 + t502;
t487 = -t188 / 0.4e1 + t222 / 0.4e1 + t98 / 0.4e1 + t144 * t684 + t145 * t677;
t488 = t143 / 0.4e1 - t189 / 0.4e1 - t221 / 0.4e1 + t507;
t2 = t470 + t486 * t373 + (-t147 / 0.4e1 - pkin(10) * t230 / 0.2e1 - t70 * mrSges(7,3) / 0.2e1) * t464 + (-t148 / 0.4e1 + pkin(10) * t231 / 0.2e1 + t69 * mrSges(7,3) / 0.2e1) * t461 + ((-0.3e1 / 0.4e1 * Ifges(6,6) + t414 / 0.4e1) * t614 + (t722 * t97 + t704 + t733) * t457 + t487) * t465 + (-0.3e1 / 0.4e1 * t581 + (m(7) * t715 + t506) * t457 + t488) * t462 + t535 * t374 + t472;
t28 = t303 * t397 + t290 * t396 + t302 * t399 + t289 * t398 + m(7) * (t289 * t302 + t290 * t303) + t458 * t411 + (t457 * t381 + t342 * t678 + t340 * t685 - t533 / 0.2e1 - t552) * t465 + (-t597 / 0.2e1 + t608 / 0.2e1 + t338 / 0.2e1 - t420 / 0.2e1 + t531 / 0.2e1 + (-m(7) * t611 + t382) * t457) * t462;
t41 = (t381 / 0.2e1 + t518 * t721 - t501) * t465 + (t695 + (t517 - 0.2e1 * t611) * t721 + t596 / 0.2e1 - t606 / 0.2e1) * t462;
t495 = t2 * qJD(1) + t28 * qJD(2) + t41 * qJD(4);
t476 = (t526 + t84) * t721 + t599 / 0.2e1 - t609 / 0.2e1 + t506;
t477 = (t528 - t97) * t721 + t494 + t733;
t13 = t496 * t722 + (t617 / 0.2e1 - t616 / 0.2e1) * mrSges(7,3) + ((-t720 / 0.2e1 - t534 / 0.2e1 + t719) * t614 + t476) * t462 + (mrSges(6,2) * t566 + t477) * t465;
t291 = (-0.1e1 + t585) * t465 * t462;
t493 = t13 * qJD(1) + t41 * qJD(2) + t291 * t667;
t492 = t657 / 0.2e1 + t74 * mrSges(7,1) / 0.2e1 + t75 * t717;
t126 = pkin(5) * t412 + t530 * t679 + t532 * t685 + t500;
t257 = (t690 + t508) * t465;
t30 = t672 / 0.2e1 + (Ifges(7,3) / 0.2e1 + t457 * t690 + pkin(10) * t540) * t465 + (-0.3e1 / 0.4e1 * t631 + pkin(10) * t692 - t343 / 0.4e1 + (t464 * t578 + t687) * t465) * t464 + (0.3e1 / 0.4e1 * t630 + pkin(10) * t694 + t341 / 0.4e1 + (t662 + t418 / 0.4e1 - t578 * t461) * t465) * t461 + t511;
t484 = pkin(5) * t133 / 0.2e1 + t261 * t686 + t262 * t687 - t84 * t412 / 0.2e1;
t8 = (-0.3e1 / 0.4e1 * t645 - t136 / 0.4e1 - t652 / 0.4e1 - t249 / 0.4e1 - t100 / 0.4e1 + t537 * pkin(10)) * t464 + (0.3e1 / 0.4e1 * t644 + t99 / 0.4e1 + t651 / 0.4e1 + t654 / 0.4e1 - t137 / 0.4e1 + t538 * pkin(10)) * t461 + t484 + t492;
t490 = t8 * qJD(1) + t30 * qJD(2) + t257 * qJD(4) + t126 * qJD(5);
t258 = t465 * t508 + t555;
t186 = t481 + t513;
t78 = t462 * t508 + t501;
t62 = t681 * t741 + t480 + t483;
t59 = -t589 / 0.2e1 + t397 * t561 + t509 + t732 * t462;
t56 = t499 / 0.2e1 + t456 * t724 + t479;
t40 = t41 * qJD(5);
t31 = -t672 / 0.2e1 + t529 * t680 + t397 * t579 - t418 * t604 / 0.2e1 - t415 * t593 / 0.2e1 + t457 * t555 + t497 * t684 + t498 * t677 + t656 / 0.2e1 + t512 * t462 - t502 + t511 + t732 * pkin(10);
t27 = t373 * t540 + t374 * t691 + t482 + t491;
t24 = t261 * t568 + t373 * t508 + t461 * t572 + t503;
t21 = t478 - t510;
t14 = t469 + (-t463 * mrSges(5,1) + t382 * t560) * t455 + t353 * t675 + t179 * t676 + t354 * t681 + t231 * t561 + t462 * t558 - t474;
t12 = (pkin(5) * t576 + t496) * t721 - t534 * t544 + t351 * t569 + t352 * t568 + mrSges(6,2) * t542 + mrSges(6,1) * t543 + t477 * t465 + t476 * t462;
t10 = mrSges(6,1) * t542 + mrSges(6,2) * t544 + t462 * t489 + t465 * t505 + t435 + t468 - t475 + t582;
t9 = t136 * t677 + t137 * t682 + t530 * t708 + t532 * t706 - t484 + t492 - t507 + (t529 / 0.4e1 + t512) * t373 + t737 * pkin(10);
t6 = t471 * t465 + t473 - t485 + t736;
t1 = t470 + (Ifges(6,5) * t565 + t488) * t462 + (-mrSges(6,3) * t612 / 0.2e1 + t535) * t374 + (mrSges(6,3) * t611 / 0.2e1 + t486) * t373 + (t462 * t703 + t134 * t681 + t591 / 0.2e1 - t193 * t676 + (t462 * t84 - t625) * t721) * t457 + t147 * t677 + t148 * t682 + Ifges(6,5) * t543 + Ifges(6,6) * t542 + t231 * t579 + t69 * t569 + pkin(10) * t558 + t70 * t568 - t472 + (t487 + (Ifges(6,6) + t414) * t565) * t465;
t18 = [qJD(2) * t3 + qJD(3) * t16 + qJD(4) * t17 + qJD(5) * t4 + qJD(6) * t7, t10 * qJD(3) + t14 * qJD(4) + t1 * qJD(5) + t6 * qJD(6) + t650 + (t641 + t438 + t439 + t440 + 0.2e1 * (qJ(3) * t298 - t299 * t458) * t725 + 0.2e1 * (-pkin(2) * t379 - qJ(3) * t378) * t727 + 0.2e1 * (-t103 * t611 + t289 * t69 + t290 * t70) * t721 + 0.2e1 * (-t260 * t458 + t457 * t519) * t723 + t341 * t699 + t343 * t698 + t437 + (-t117 * mrSges(6,3) + t457 * t354 + t553) * t462 + (t705 - t116 * mrSges(6,3) + t147 * t685 + t148 * t678 + (t353 - t179) * t457) * t465 + ((-Ifges(3,6) + Ifges(6,5) * t676 + Ifges(6,6) * t681 + (-mrSges(5,1) - mrSges(4,1)) * qJ(3)) * t463 + (-pkin(2) * mrSges(4,1) - t458 * mrSges(5,1) + t420 * t681 + t465 * t552 - Ifges(4,4)) * t466) * t455 - t458 * t319 - t260 * t413 + t70 * t397 + t69 * t399 - t379 * mrSges(3,1) + t103 * t382 - t378 * mrSges(4,3) + t379 * mrSges(4,2) + t289 * t231 + t290 * t230 + t298 * mrSges(5,2) - t299 * mrSges(5,3)) * qJD(2), qJD(2) * t10 + qJD(4) * t56 + qJD(5) * t27 + qJD(6) * t24 + t622, t621 + t14 * qJD(2) + t56 * qJD(3) + t12 * qJD(5) + t21 * qJD(6) + (t516 + t575) * t462 * t667, t640 + t1 * qJD(2) + t27 * qJD(3) + t12 * qJD(4) + (-t96 * mrSges(6,2) + t526 * mrSges(7,3) + pkin(5) * t193 + t144 * t678 + t145 * t683 + t500 * t373 + t374 * t689 + t586 + t738 * t97 + (m(7) * t526 + t599 - t609) * pkin(10)) * qJD(5) + t9 * qJD(6), t624 + t6 * qJD(2) + t24 * qJD(3) + t21 * qJD(4) + t9 * qJD(5) + (-mrSges(7,1) * t55 - mrSges(7,2) * t54 + t135) * qJD(6); qJD(3) * t11 + qJD(4) * t15 + qJD(5) * t2 + qJD(6) * t5 - t650, -qJD(3) * t57 + qJD(4) * t76 + qJD(5) * t28 - qJD(6) * t32, qJD(4) * t186 + qJD(5) * t62 + qJD(6) * t78 - t524, qJD(3) * t186 + qJD(6) * t59 + t40 - t523, t62 * qJD(3) + t31 * qJD(6) + t495 + (pkin(5) * t381 + t340 * t678 + t342 * t683 + (m(7) * t517 + t596 - t606) * pkin(10) + (t457 * t738 - Ifges(6,5) + t500) * t462 + (-t457 * mrSges(6,2) - Ifges(6,6) + t689) * t465 + t517 * mrSges(7,3)) * qJD(5), t78 * qJD(3) + t59 * qJD(4) + t31 * qJD(5) + (-mrSges(7,1) * t290 - mrSges(7,2) * t289 - t465 * t414) * qJD(6) + t525; -qJD(2) * t11 - qJD(4) * t53 - qJD(5) * t26 - qJD(6) * t23 - t622, qJD(4) * t180 + qJD(5) * t60 - qJD(6) * t77 + t524, 0, -t514, -t520, qJD(6) * t412 - t521; -qJD(2) * t15 + qJD(3) * t53 + qJD(5) * t13 + qJD(6) * t22 - t621, -qJD(3) * t180 - qJD(6) * t58 + t40 + t523, t514, m(7) * t291 * qJD(5) (-t739 + m(7) * (pkin(10) * t549 - t671) + mrSges(7,3) * t549 - t413) * qJD(5) + t258 * qJD(6) + t493, t258 * qJD(5) - qJD(6) * t739 + t522; -qJD(2) * t2 + qJD(3) * t26 - qJD(4) * t13 - qJD(6) * t8 - t640, -qJD(3) * t60 - qJD(6) * t30 - t495, t520, -t257 * qJD(6) - t493, -t126 * qJD(6) (-pkin(10) * t534 + t529) * qJD(6) - t490; -qJD(2) * t5 + qJD(3) * t23 - qJD(4) * t22 + qJD(5) * t8 - t624, qJD(3) * t77 + qJD(4) * t58 + qJD(5) * t30 - t525, t521, qJD(5) * t257 - t522, t490, 0;];
Cq  = t18;
