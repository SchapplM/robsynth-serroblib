% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRR7
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:16
% EndTime: 2019-03-09 09:17:37
% DurationCPUTime: 10.82s
% Computational Cost: add. (18359->830), mult. (41960->1125), div. (0->0), fcn. (40762->8), ass. (0->436)
t764 = -Ifges(7,5) / 0.2e1;
t763 = Ifges(7,6) / 0.2e1;
t493 = cos(qJ(6));
t662 = t493 * mrSges(7,1);
t490 = sin(qJ(6));
t664 = t490 * mrSges(7,2);
t424 = -t662 + t664;
t655 = -mrSges(6,1) + t424;
t745 = m(7) * pkin(5);
t762 = t655 - t745;
t487 = sin(pkin(6));
t492 = sin(qJ(2));
t641 = t487 * t492;
t495 = cos(qJ(2));
t640 = t487 * t495;
t331 = -t487 * pkin(1) - pkin(2) * t640 - qJ(3) * t641;
t282 = pkin(3) * t640 - t331;
t201 = (pkin(4) * t492 + pkin(9) * t495) * t487 + t282;
t449 = qJ(4) * t641;
t488 = cos(pkin(6));
t698 = pkin(1) * t488;
t464 = t495 * t698;
t569 = pkin(8) * t641 - t464;
t312 = -t449 + t569;
t496 = -pkin(2) - pkin(3);
t482 = -pkin(9) + t496;
t231 = t482 * t488 + t312;
t491 = sin(qJ(5));
t494 = cos(qJ(5));
t97 = t201 * t494 - t231 * t491;
t88 = -pkin(5) * t641 - t97;
t737 = m(7) * t88;
t386 = t488 * t491 + t494 * t640;
t601 = t493 * t641;
t270 = t386 * t490 + t601;
t271 = -t386 * t493 + t490 * t641;
t706 = -t493 / 0.2e1;
t710 = t490 / 0.2e1;
t761 = mrSges(7,3) * (t270 * t710 + t271 * t706);
t668 = t386 * mrSges(6,3);
t669 = t386 * mrSges(6,2);
t385 = -t488 * t494 + t491 * t640;
t671 = t385 * mrSges(6,1);
t233 = -t669 - t671;
t475 = t488 * mrSges(5,1);
t605 = mrSges(5,3) * t640;
t406 = -t475 + t605;
t760 = t406 - t233;
t480 = Ifges(7,5) * t493;
t680 = Ifges(7,6) * t490;
t759 = t480 - t680;
t481 = Ifges(7,4) * t493;
t758 = -Ifges(7,2) * t490 + t481;
t432 = Ifges(7,1) * t490 + t481;
t483 = t490 ^ 2;
t485 = t493 ^ 2;
t613 = t483 + t485;
t620 = t493 * t494;
t328 = (t490 * t495 + t492 * t620) * t487;
t630 = t491 * t492;
t602 = t487 * t630;
t247 = mrSges(7,1) * t602 - mrSges(7,3) * t328;
t628 = t492 * t494;
t600 = t487 * t628;
t327 = t490 * t600 - t493 * t640;
t246 = -mrSges(7,2) * t602 - mrSges(7,3) * t327;
t704 = t493 / 0.2e1;
t588 = t246 * t704;
t711 = -t490 / 0.2e1;
t757 = t247 * t711 + t588;
t695 = pkin(10) * t494;
t697 = pkin(5) * t491;
t437 = t695 - t697;
t633 = t490 * t491;
t315 = t437 * t493 + t482 * t633;
t629 = t491 * t493;
t316 = t437 * t490 - t482 * t629;
t756 = -t315 * t490 + t316 * t493;
t755 = t327 * t763 + t328 * t764;
t184 = -mrSges(7,1) * t385 - t271 * mrSges(7,3);
t625 = t493 * t184;
t183 = mrSges(7,2) * t385 + t270 * mrSges(7,3);
t638 = t490 * t183;
t534 = t638 / 0.2e1 + t625 / 0.2e1;
t514 = -t534 + t761;
t451 = qJ(3) * t640;
t313 = t496 * t641 + t451;
t473 = t488 * qJ(3);
t754 = 0.2e1 * t473;
t753 = m(4) / 0.2e1;
t752 = -m(5) / 0.2e1;
t751 = m(5) / 0.2e1;
t750 = -m(6) / 0.2e1;
t749 = m(6) / 0.2e1;
t748 = -m(7) / 0.2e1;
t747 = m(7) / 0.2e1;
t746 = -pkin(10) / 0.2e1;
t744 = -mrSges(7,1) / 0.2e1;
t743 = mrSges(7,1) / 0.2e1;
t742 = -mrSges(7,2) / 0.2e1;
t741 = -mrSges(7,3) / 0.2e1;
t740 = mrSges(7,3) / 0.2e1;
t739 = Ifges(6,6) / 0.2e1;
t738 = t88 / 0.2e1;
t259 = Ifges(7,4) * t270;
t104 = t271 * Ifges(7,1) - Ifges(7,5) * t385 + t259;
t735 = t104 / 0.2e1;
t677 = t271 * Ifges(7,4);
t150 = Ifges(7,1) * t270 - t677;
t734 = -t150 / 0.4e1;
t733 = t184 / 0.2e1;
t672 = t328 * mrSges(7,2);
t673 = t327 * mrSges(7,1);
t187 = t672 + t673;
t732 = t187 / 0.2e1;
t484 = t491 ^ 2;
t486 = t494 ^ 2;
t612 = t484 - t486;
t256 = (0.1e1 - t613) * t612;
t731 = -t256 / 0.2e1;
t730 = t270 / 0.2e1;
t729 = t271 / 0.2e1;
t489 = qJ(3) + pkin(4);
t696 = pkin(5) * t494;
t420 = pkin(10) * t491 + t489 + t696;
t632 = t490 * t494;
t295 = t420 * t493 - t482 * t632;
t728 = t295 / 0.2e1;
t727 = t316 / 0.2e1;
t726 = -t327 / 0.2e1;
t467 = Ifges(7,4) * t633;
t683 = Ifges(7,5) * t494;
t361 = -Ifges(7,1) * t629 + t467 + t683;
t725 = t361 / 0.2e1;
t724 = -t385 / 0.2e1;
t723 = -t385 / 0.4e1;
t722 = -t386 / 0.2e1;
t661 = t493 * mrSges(7,2);
t665 = t490 * mrSges(7,1);
t426 = t661 + t665;
t392 = t491 * t426;
t721 = -t392 / 0.2e1;
t407 = -mrSges(7,2) * t494 + mrSges(7,3) * t633;
t720 = -t407 / 0.2e1;
t719 = t407 / 0.2e1;
t409 = mrSges(7,1) * t494 + mrSges(7,3) * t629;
t718 = -t409 / 0.2e1;
t717 = t409 / 0.2e1;
t716 = -t424 / 0.2e1;
t715 = -t426 / 0.2e1;
t685 = Ifges(7,4) * t490;
t429 = Ifges(7,2) * t493 + t685;
t714 = t429 / 0.4e1;
t713 = -t432 / 0.4e1;
t712 = t482 / 0.2e1;
t709 = t490 / 0.4e1;
t708 = -t491 / 0.2e1;
t707 = t491 / 0.2e1;
t705 = -t493 / 0.4e1;
t703 = t493 / 0.4e1;
t702 = -t494 / 0.2e1;
t701 = t494 / 0.2e1;
t700 = t494 / 0.4e1;
t699 = -t495 / 0.2e1;
t390 = pkin(8) * t640 + t492 * t698;
t314 = -qJ(4) * t640 + t390;
t274 = -t314 - t473;
t254 = t488 * pkin(4) - t274;
t124 = -pkin(5) * t385 + pkin(10) * t386 + t254;
t98 = t201 * t491 + t231 * t494;
t89 = pkin(10) * t641 + t98;
t60 = t124 * t493 - t490 * t89;
t694 = t60 * mrSges(7,3);
t61 = t124 * t490 + t493 * t89;
t693 = t61 * mrSges(7,3);
t692 = -Ifges(5,5) - Ifges(3,6);
t530 = m(7) * t756;
t408 = mrSges(7,2) * t491 + mrSges(7,3) * t632;
t622 = t493 * t408;
t410 = -mrSges(7,1) * t491 + mrSges(7,3) * t620;
t634 = t490 * t410;
t513 = t530 / 0.2e1 + t622 / 0.2e1 - t634 / 0.2e1 + t721;
t393 = t494 * t426;
t525 = t393 / 0.2e1 + t409 * t711 + t407 * t704;
t296 = t420 * t490 + t482 * t620;
t551 = -t295 * t490 + t296 * t493;
t42 = m(7) * t612 * t712 + (t551 * t747 + t525) * t494 + t513 * t491;
t391 = t491 * t424;
t621 = t493 * t409;
t635 = t490 * t407;
t532 = -t635 / 0.2e1 - t621 / 0.2e1;
t568 = mrSges(7,3) * (t485 / 0.2e1 + t483 / 0.2e1);
t69 = t391 * t702 + t484 * t568 + t491 * t532;
t691 = t42 * qJD(5) + t69 * qJD(6);
t690 = mrSges(7,3) * t385;
t688 = Ifges(6,4) * t386;
t687 = Ifges(6,4) * t491;
t686 = Ifges(6,4) * t494;
t479 = Ifges(6,6) * t491;
t679 = Ifges(7,3) * t386;
t678 = Ifges(7,3) * t491;
t676 = t295 * mrSges(7,3);
t675 = t296 * mrSges(7,3);
t103 = t270 * Ifges(7,2) - Ifges(7,6) * t385 + t677;
t643 = t482 * t492;
t230 = t451 + (pkin(4) * t495 + t643) * t487;
t125 = t230 * t494 - t314 * t491;
t108 = -pkin(5) * t640 - t125;
t126 = t491 * t230 + t494 * t314;
t147 = -mrSges(7,1) * t270 + mrSges(7,2) * t271;
t158 = Ifges(7,4) * t328 - Ifges(7,2) * t327 + Ifges(7,6) * t602;
t159 = Ifges(7,1) * t328 - Ifges(7,4) * t327 + Ifges(7,5) * t602;
t376 = Ifges(6,4) * t385;
t606 = Ifges(6,5) * t641;
t198 = -Ifges(6,1) * t386 + t376 + t606;
t255 = t488 * t496 + t312;
t564 = Ifges(6,1) * t494 - t687;
t273 = (Ifges(6,5) * t495 + t492 * t564) * t487;
t670 = t385 * mrSges(6,3);
t286 = -mrSges(6,2) * t641 + t670;
t287 = mrSges(6,1) * t641 + t668;
t659 = t494 * mrSges(6,2);
t565 = t491 * mrSges(6,1) + t659;
t329 = t565 * t641;
t330 = t473 + t390;
t336 = -pkin(2) * t488 + t569;
t369 = (-mrSges(6,2) * t495 - mrSges(6,3) * t630) * t487;
t370 = (mrSges(6,1) * t495 - mrSges(6,3) * t628) * t487;
t387 = (-mrSges(4,1) * t495 - mrSges(4,3) * t492) * t487;
t388 = pkin(2) * t641 - t451;
t405 = mrSges(5,2) * t488 - mrSges(5,3) * t641;
t455 = Ifges(4,6) * t641;
t457 = Ifges(3,5) * t640;
t458 = Ifges(4,4) * t640;
t562 = -Ifges(6,2) * t491 + t686;
t573 = Ifges(7,3) * t602;
t583 = -t573 / 0.2e1 + (Ifges(6,6) * t495 + t492 * t562) * t487 / 0.2e1 + t755;
t102 = Ifges(7,5) * t271 + Ifges(7,6) * t270 - Ifges(7,3) * t385;
t197 = Ifges(6,2) * t385 + Ifges(6,6) * t641 - t688;
t585 = t102 / 0.2e1 - t197 / 0.2e1;
t453 = mrSges(4,2) * t640;
t474 = t488 * mrSges(4,3);
t614 = t474 + t453;
t615 = mrSges(5,1) * t640 + mrSges(5,2) * t641;
t667 = t569 * mrSges(3,2);
t109 = pkin(10) * t640 + t126;
t423 = pkin(5) * t602;
t200 = t423 + t449 + t464 + (-pkin(8) - t695) * t641;
t80 = -t109 * t490 + t200 * t493;
t81 = t109 * t493 + t200 * t490;
t3 = ((Ifges(3,4) * t640 - t331 * mrSges(4,3) - t313 * mrSges(5,2) + Ifges(6,5) * t722 + t385 * t739 - t255 * mrSges(5,3) + t336 * mrSges(4,2) + (-pkin(1) * mrSges(3,2) + (Ifges(5,4) - Ifges(4,5)) * t495) * t487 + (Ifges(3,5) / 0.2e1 + Ifges(5,6) + Ifges(4,4) / 0.2e1) * t488) * t495 + (Ifges(4,5) * t641 + t331 * mrSges(4,1) + t313 * mrSges(5,1) + t198 * t701 - t274 * mrSges(5,3) + t585 * t491 + (-pkin(1) * mrSges(3,1) + (Ifges(6,5) * t701 - t479 / 0.2e1 - Ifges(3,4) - Ifges(5,4)) * t492) * t487 + (t390 - t330) * mrSges(4,2) + (Ifges(4,6) / 0.2e1 + t692) * t488 + (Ifges(3,1) + Ifges(4,1) - Ifges(5,1) - Ifges(3,2) + Ifges(5,2) - Ifges(4,3) + Ifges(6,3)) * t640) * t492) * t487 + t282 * t615 + (t667 + t457 / 0.2e1 + t458 / 0.2e1 + t455 / 0.2e1 + (-mrSges(4,1) - mrSges(3,1)) * t390) * t488 + t583 * t385 - t569 * t614 + m(4) * (-t330 * t569 + t331 * t388 + t336 * t390) + t273 * t722 + t103 * t726 + t159 * t729 + t158 * t730 + t328 * t735 + t314 * t405 + t388 * t387 + t97 * t370 + t98 * t369 + t254 * t329 + t126 * t286 + t125 * t287 + t61 * t246 + t60 * t247 + t88 * t187 + t81 * t183 + t80 * t184 + t108 * t147 + m(6) * (t125 * t97 + t126 * t98 - t254 * t312) + m(5) * (t255 * t314 + t274 * t312 + t282 * t313) + m(7) * (t108 * t88 + t60 * t80 + t61 * t81) + t760 * t312;
t674 = t3 * qJD(1);
t154 = t385 * t759 - t679;
t155 = -Ifges(7,6) * t386 + t385 * t758;
t563 = Ifges(7,1) * t493 - t685;
t156 = -Ifges(7,5) * t386 + t385 * t563;
t202 = t426 * t385;
t221 = mrSges(7,2) * t386 - t490 * t690;
t222 = -mrSges(7,1) * t386 - t493 * t690;
t232 = -mrSges(6,1) * t386 + mrSges(6,2) * t385;
t234 = Ifges(6,2) * t386 + t376;
t235 = Ifges(6,1) * t385 + t688;
t594 = t641 / 0.2e1;
t616 = Ifges(6,5) * t385 + Ifges(6,6) * t386;
t627 = t493 * t104;
t639 = t490 * t103;
t236 = -pkin(5) * t386 - pkin(10) * t385;
t76 = t236 * t493 - t490 * t97;
t77 = t236 * t490 + t493 * t97;
t4 = t88 * t202 + m(7) * (t60 * t76 + t61 * t77) + t76 * t184 + t60 * t222 + t77 * t183 + t61 * t221 + t155 * t730 + t156 * t729 + t97 * t286 + t254 * t232 + t616 * t594 + (-t235 / 0.2e1 - t585) * t386 + (-t154 / 0.2e1 + t198 / 0.2e1 + t234 / 0.2e1 - t639 / 0.2e1 + t627 / 0.2e1 - t97 * mrSges(6,3)) * t385 + (t147 - t287 + t668 + t737) * t98;
t666 = t4 * qJD(1);
t663 = t490 * t80;
t660 = t493 * t81;
t658 = t494 * mrSges(6,3);
t657 = t494 * Ifges(7,6);
t146 = mrSges(7,1) * t271 + mrSges(7,2) * t270;
t148 = Ifges(7,5) * t270 - Ifges(7,6) * t271;
t149 = -Ifges(7,2) * t271 + t259;
t7 = t88 * t146 + t148 * t724 - t61 * t184 + t60 * t183 + (-t693 - t103 / 0.2e1 + t150 / 0.2e1) * t271 + (-t694 + t149 / 0.2e1 + t735) * t270;
t656 = t7 * qJD(1);
t425 = t494 * mrSges(6,1) - t491 * mrSges(6,2);
t603 = t484 * t641;
t499 = -t313 * t751 + (t295 * t327 - t296 * t328 - t482 * t603) * t747 + t327 * t717 + t328 * t720 + (t425 * t699 + (-t392 * t708 + (t486 / 0.2e1 + t484 / 0.2e1) * mrSges(6,3)) * t492 + (-t489 * t495 + (-t484 - t486) * t643) * t749) * t487;
t558 = t660 - t663;
t505 = t313 * t752 + (t125 * t494 + t126 * t491) * t750 + (-t108 * t494 + t491 * t558) * t748;
t524 = t369 / 0.2e1 + t757;
t582 = t732 - t370 / 0.2e1;
t18 = -t491 * t524 + t494 * t582 + t499 + t505 - t615;
t654 = qJD(1) * t18;
t439 = t486 * t641;
t647 = t328 * t493;
t648 = t327 * t490;
t550 = -t647 - t648;
t516 = m(7) * (t494 * t550 - t603);
t368 = t488 * t490 + t491 * t601;
t645 = t368 * t493;
t367 = t488 * t493 - t490 * t602;
t646 = t367 * t490;
t549 = t645 - t646;
t521 = (t491 * t549 + t439) * t747;
t48 = m(5) * t641 + t521 - t516 / 0.2e1 + 0.2e1 * (t484 * t594 + t439 / 0.2e1) * m(6);
t653 = qJD(1) * t48;
t619 = t494 * t147;
t16 = -t367 * t184 - t368 * t183 - t286 * t602 - t287 * t600 - t487 ^ 2 * t492 * (mrSges(5,1) * t492 - mrSges(5,2) * t495) - m(7) * (t367 * t60 + t368 * t61 - t600 * t88) + (t387 + t619 - m(6) * (t491 * t98 + t494 * t97) - m(5) * t282 + m(4) * t331) * t641 + (-m(4) * t330 + m(5) * t274 - m(6) * t254 - t614 + t760) * t488;
t652 = t16 * qJD(1);
t618 = t494 * t286;
t631 = t491 * t147;
t17 = -m(7) * (t327 * t60 - t328 * t61) - t327 * t184 + t328 * t183 - t287 * t602 - t406 * t640 + (t495 * t233 + (t405 + t618 + t631) * t492 + t630 * t737 - m(6) * (-t254 * t495 - t628 * t98 + t630 * t97) - m(5) * (-t255 * t492 + t274 * t495)) * t487;
t651 = t17 * qJD(1);
t644 = t482 * t491;
t642 = t482 * t494;
t637 = t490 * t222;
t359 = -t491 * t758 + t657;
t636 = t490 * t359;
t626 = t493 * t183;
t624 = t493 * t221;
t623 = t493 * t361;
t617 = t494 * t424;
t394 = Ifges(7,5) * t633 + Ifges(7,6) * t629;
t611 = qJD(5) * t494;
t574 = t613 * t494;
t303 = (-t494 + t574) * t491;
t610 = t303 * qJD(3);
t609 = t737 / 0.2e1;
t608 = t98 * t747;
t607 = m(7) * t303 * qJD(5);
t604 = t256 * t747;
t599 = t480 / 0.2e1;
t598 = mrSges(6,3) * t707;
t597 = t658 / 0.2e1;
t596 = t716 + mrSges(6,1) / 0.2e1;
t595 = -t641 / 0.2e1;
t591 = -t633 / 0.2e1;
t590 = -t632 / 0.2e1;
t589 = t492 * t479 / 0.4e1;
t586 = t620 / 0.2e1;
t584 = t147 / 0.2e1 - t287 / 0.2e1;
t581 = t202 / 0.2e1 - t286 / 0.2e1;
t396 = t432 * t491;
t580 = -t359 / 0.4e1 + t396 / 0.4e1;
t395 = Ifges(7,2) * t629 + t467;
t579 = t361 / 0.4e1 + t395 / 0.4e1;
t357 = Ifges(7,3) * t494 - t491 * t759;
t430 = -Ifges(6,2) * t494 - t687;
t578 = -t430 / 0.2e1 + t357 / 0.2e1;
t575 = t613 * t491;
t572 = t491 * t595;
t571 = t494 * t595;
t570 = t494 * t594;
t566 = t430 / 0.4e1 + t564 / 0.4e1 - t357 / 0.4e1;
t559 = -t490 * t76 + t493 * t77;
t537 = t668 / 0.2e1 + t584;
t547 = t559 + t88;
t506 = t547 * t748 + t637 / 0.2e1 - t624 / 0.2e1 - t537;
t510 = (-t647 / 0.2e1 - t648 / 0.2e1) * mrSges(7,3) + (pkin(10) * t550 + t423) * t747;
t536 = t670 / 0.2e1 + t581;
t512 = t184 * t711 + t626 / 0.2e1 + mrSges(6,2) * t595 - t536;
t548 = t490 * t60 - t493 * t61 + t98;
t10 = (t548 * t747 - t512) * t494 + (t596 * t641 + t506) * t491 + t510;
t557 = t10 * qJD(1) - t42 * qJD(2);
t526 = (t645 / 0.2e1 - t646 / 0.2e1) * mrSges(7,3);
t529 = t549 * pkin(10);
t12 = t529 * t747 + t526 + ((t745 / 0.2e1 + t596) * t641 + t506) * t494 + (t548 * t748 + t512) * t491;
t44 = t513 * t494 + ((-t551 + 0.2e1 * t642) * t747 - t525) * t491;
t556 = t12 * qJD(1) - t44 * qJD(2);
t497 = -t474 - t475 - m(4) * (t754 + t390) / 0.2e1 + (t314 + t754) * t752 + (t488 * t489 + t254) * t750 + (t295 * t367 + t296 * t368 + t490 * t61 + t493 * t60 - t600 * t644) * t748 + t367 * t718 + t368 * t720 + t671 / 0.2e1 + t669 / 0.2e1 - t488 * t425 / 0.2e1 - t534;
t552 = -t125 * t491 + t126 * t494;
t502 = t582 * t491 + t390 * t753 + t314 * t751 + t552 * t749 + (t108 * t491 + t494 * t558) * t747;
t15 = t497 + (-t392 * t594 + t524) * t494 + t502;
t82 = t635 + t621 + mrSges(5,1) + mrSges(4,3) + (m(4) + m(5)) * qJ(3) + m(7) * (t295 * t493 + t296 * t490) + m(6) * t489 + t425;
t555 = -qJD(1) * t15 + qJD(2) * t82;
t508 = t146 * t702 + t491 * t514;
t543 = -t673 / 0.2e1 - t672 / 0.2e1;
t24 = t508 + t543;
t554 = -qJD(1) * t24 - qJD(2) * t69;
t507 = t146 * t707 + t494 * t514;
t542 = t367 * t743 + t368 * t742;
t26 = t507 - t542;
t511 = (t491 * t568 + t532) * t494 + t391 * t707;
t541 = t664 / 0.2e1 - t662 / 0.2e1;
t63 = t511 + t541;
t553 = t26 * qJD(1) + t63 * qJD(2);
t546 = -pkin(5) * t391 / 0.2e1 + t480 * t700;
t545 = qJD(3) * t731 - t303 * qJD(4);
t544 = mrSges(7,2) * t727 + t315 * t744;
t540 = t661 / 0.2e1 + t665 / 0.2e1;
t539 = pkin(10) * t720 + t580;
t538 = pkin(10) * t718 + t579;
t535 = t639 / 0.4e1 - t627 / 0.4e1;
t531 = t429 * t710 + t432 * t706;
t528 = t715 + t540;
t34 = t394 * t701 - t296 * t409 + t295 * t407 + (t482 * t391 + (t359 / 0.2e1 - t396 / 0.2e1 + t675) * t493 + (t395 / 0.2e1 + t725 - t676) * t490) * t491;
t501 = (-t676 / 0.2e1 + t579) * t270 + (-t675 / 0.2e1 + t580) * t271 + t183 * t728 - t296 * t184 / 0.2e1 + t394 * t723 + t148 * t700 + t60 * t719 + t61 * t718 + t391 * t738;
t509 = t146 * t712 + (t734 + t103 / 0.4e1 + t693 / 0.2e1) * t493 + (t149 / 0.4e1 + t104 / 0.4e1 - t694 / 0.2e1) * t490;
t517 = t80 * t744 + t81 * mrSges(7,2) / 0.2e1 + t755;
t5 = (Ifges(7,3) * t595 + t509) * t491 + t501 + t517;
t527 = t5 * qJD(1) + t34 * qJD(2) + t69 * qJD(4);
t523 = -t679 / 0.2e1 + t76 * t743 + t77 * t742;
t428 = Ifges(7,5) * t490 + Ifges(7,6) * t493;
t360 = -Ifges(7,6) * t491 - t494 * t758;
t362 = -Ifges(7,5) * t491 - t494 * t563;
t498 = (t295 * t76 + t296 * t77 + t315 * t60 + t316 * t61) * t747 - t254 * t565 / 0.2e1 + t270 * t360 / 0.4e1 + t271 * t362 / 0.4e1 + t222 * t728 + t296 * t221 / 0.2e1 + t315 * t733 + t183 * t727 + t489 * t232 / 0.2e1 + t60 * t410 / 0.2e1 + t61 * t408 / 0.2e1 + t76 * t717 + t77 * t719 - t393 * t738 + t98 * t721;
t503 = (-pkin(5) * t108 + pkin(10) * t558) * t748 + pkin(5) * t732 + t108 * t716 - t125 * mrSges(6,1) / 0.2e1 + t126 * mrSges(6,2) / 0.2e1 + t327 * t714 + t328 * t713;
t358 = -t494 * t759 - t678;
t433 = -Ifges(6,1) * t491 - t686;
t518 = t623 / 0.4e1 - t636 / 0.4e1 - t562 / 0.4e1 + t433 / 0.4e1 - t358 / 0.4e1;
t519 = t154 / 0.4e1 - t198 / 0.4e1 - t234 / 0.4e1 + t535;
t520 = -t102 / 0.4e1 + t197 / 0.4e1 - t235 / 0.4e1 + t155 * t709 + t156 * t705;
t1 = t498 + t566 * t386 + (-0.3e1 / 0.4e1 * t606 + (t609 + t537) * t482 + t519) * t494 + ((-t428 / 0.4e1 + t739) * t641 + (t608 + t536) * t482 + t520) * t491 + (t246 * t746 + t81 * t741 - t158 / 0.4e1) * t493 + (pkin(10) * t247 / 0.2e1 + t80 * t740 - t159 / 0.4e1) * t490 + t518 * t385 + (Ifges(6,3) * t699 + t589) * t487 + t503;
t28 = t315 * t409 + t295 * t410 + m(7) * (t295 * t315 + t296 * t316) + t316 * t407 + t296 * t408 - t489 * t565 + (t636 / 0.2e1 - t482 * t392 - t623 / 0.2e1 + t358 / 0.2e1 + t562 / 0.2e1 - t433 / 0.2e1) * t494 + (t360 * t710 + t362 * t706 + t564 / 0.2e1 + (m(7) * t642 - t393) * t482 - t578) * t491;
t522 = t1 * qJD(1) + t28 * qJD(2) + t44 * qJD(3) + t42 * qJD(4);
t139 = pkin(5) * t426 + t563 * t711 + t706 * t758 + t531;
t266 = t528 * t491;
t267 = t528 * t494;
t504 = pkin(10) * t568 + t426 * t712 + t429 * t703 + t563 * t705 + (t432 + t758) * t709;
t32 = (t683 / 0.2e1 + t538) * t493 + (-0.3e1 / 0.4e1 * t657 + t539) * t490 + (Ifges(7,3) / 0.2e1 + t504) * t491 + t544 + t546;
t500 = pkin(5) * t146 / 0.2e1 + t490 * t734 + t149 * t705 + t88 * t715 + t535 + (t714 - t563 / 0.4e1) * t271 + (t713 - t758 / 0.4e1) * t270;
t8 = t500 + (-0.3e1 / 0.4e1 * t680 + t480 / 0.4e1 + t599) * t385 - t514 * pkin(10) + t523;
t515 = t8 * qJD(1) - t32 * qJD(2) + t266 * qJD(3) - t267 * qJD(4) + t139 * qJD(5);
t269 = t426 * t702 - t494 * t540;
t268 = t426 * t707 + t491 * t540;
t223 = qJD(5) * t604;
t62 = t511 - t541;
t55 = t516 / 0.2e1 + t521;
t43 = t44 * qJD(5);
t33 = -t678 / 0.2e1 + t620 * t764 + t632 * t763 + (-t657 / 0.4e1 + t539) * t490 + t538 * t493 + t504 * t491 - t544 + t546;
t27 = t507 + t542;
t25 = t508 - t543;
t19 = t187 * t702 + t247 * t591 + t369 * t707 + t370 * t701 + t491 * t588 + t499 - t505;
t14 = -t392 * t571 + t494 * t524 + t453 - t497 + t502 - t605;
t13 = t619 / 0.2e1 + t202 * t707 + t633 * t733 + t222 * t590 + t221 * t586 + t386 * t597 + t287 * t702 + t385 * t598 + t424 * t571 + mrSges(6,2) * t572 + mrSges(6,1) * t570 + t526 + (pkin(5) * t600 + t491 * t548 + t494 * t547 + t529) * t747 + (t626 + t286) * t708;
t11 = t631 / 0.2e1 + t202 * t702 + (t491 * t547 - t494 * t548) * t747 + t184 * t590 + t222 * t591 + t183 * t586 + t624 * t707 + t618 / 0.2e1 + t386 * t598 + t287 * t708 + t658 * t724 + (t659 / 0.2e1 + t596 * t491) * t641 + t510;
t9 = -t500 + t759 * t723 + (t599 - t680 / 0.2e1) * t385 + t523 + (t625 + t638) * t746 + pkin(10) * t761;
t6 = t501 + t509 * t491 + t573 / 0.2e1 - t517;
t2 = t498 + (t482 * t597 + t566) * t386 + (t482 * t598 + t518) * t385 + Ifges(6,3) * t640 / 0.2e1 + t660 * t740 + t663 * t741 + t487 * t589 + Ifges(6,5) * t570 + ((t608 + t581) * t482 + t520) * t491 + (-t606 / 0.4e1 + (t609 + t584) * t482 + t519) * t494 + t158 * t703 + t159 * t709 + t428 * t602 / 0.4e1 + Ifges(6,6) * t572 - t503 + t757 * pkin(10);
t20 = [qJD(2) * t3 - qJD(3) * t16 - qJD(4) * t17 + qJD(5) * t4 + qJD(6) * t7, t14 * qJD(3) + t19 * qJD(4) + t2 * qJD(5) + t6 * qJD(6) + t674 + (t458 + t455 + t457 + t667 + 0.2e1 * (-pkin(2) * t390 - qJ(3) * t569) * t753 - t569 * mrSges(4,3) + 0.2e1 * (-t312 * t489 + t482 * t552) * t749 + 0.2e1 * (-qJ(3) * t312 + t314 * t496) * t751 + 0.2e1 * (t108 * t644 + t295 * t80 + t296 * t81) * t747 + t328 * t725 + t359 * t726 + (-t126 * mrSges(6,3) + t482 * t369 - t583) * t494 + (t159 * t706 + t158 * t710 + t125 * mrSges(6,3) - t273 / 0.2e1 + (-t370 + t187) * t482) * t491 + ((-pkin(2) * mrSges(4,2) - t496 * mrSges(5,3) + Ifges(6,5) * t708 + Ifges(6,6) * t702 + Ifges(5,6)) * t495 + (t433 * t701 + t578 * t491 + (mrSges(5,3) - mrSges(4,2)) * qJ(3) + t692) * t492) * t487 + t489 * t329 - t312 * t425 + t81 * t407 + t80 * t409 - t108 * t392 - t390 * mrSges(4,1) - t390 * mrSges(3,1) + t314 * mrSges(5,2) - t312 * mrSges(5,1) + t295 * t247 + t296 * t246) * qJD(2), -t652 + t14 * qJD(2) + t55 * qJD(4) + t13 * qJD(5) + t27 * qJD(6) + m(7) * (t549 - t602) * qJD(3) * t494, -t651 + t19 * qJD(2) + t55 * qJD(3) + t11 * qJD(5) + t25 * qJD(6) + m(7) * (t550 + t600) * qJD(4) * t491, t666 + t2 * qJD(2) + t13 * qJD(3) + t11 * qJD(4) + (-t97 * mrSges(6,2) + t559 * mrSges(7,3) - pkin(5) * t202 + t155 * t704 + t156 * t710 - t531 * t385 + t428 * t722 + t616 + t762 * t98 + (m(7) * t559 + t624 - t637) * pkin(10)) * qJD(5) + t9 * qJD(6), t656 + t6 * qJD(2) + t27 * qJD(3) + t25 * qJD(4) + t9 * qJD(5) + (-mrSges(7,1) * t61 - mrSges(7,2) * t60 + t148) * qJD(6); -qJD(3) * t15 + qJD(4) * t18 + qJD(5) * t1 + qJD(6) * t5 - t674, qJD(3) * t82 + qJD(5) * t28 + qJD(6) * t34, qJD(6) * t62 + t43 + t555, t654 + t691, t33 * qJD(6) + (t482 * t762 - Ifges(6,5) + t531) * t611 + t522 + (mrSges(6,2) * t644 + pkin(5) * t393 + t360 * t704 + t362 * t710 + t428 * t708 + t479 + (t530 + t622 - t634) * pkin(10) + t756 * mrSges(7,3)) * qJD(5), t62 * qJD(3) + t33 * qJD(5) + (-mrSges(7,1) * t296 - mrSges(7,2) * t295 + t394) * qJD(6) + t527; qJD(2) * t15 - qJD(4) * t48 - qJD(5) * t12 + qJD(6) * t26 + t652, qJD(6) * t63 + t43 - t555, -t607, t223 - t653, -m(7) * t610 + qJD(4) * t604 + (t617 + m(7) * (-pkin(10) * t575 - t696) - mrSges(7,3) * t575 - t425) * qJD(5) + t268 * qJD(6) - t556, t268 * qJD(5) + qJD(6) * t617 + t553; -qJD(2) * t18 + qJD(3) * t48 - qJD(5) * t10 + qJD(6) * t24 + t651, -t654 + t691, t223 + t653, t607, t269 * qJD(6) + t655 * qJD(5) * t491 + (mrSges(7,3) * t613 - mrSges(6,2)) * t611 + ((pkin(10) * t574 - t697) * qJD(5) - t545) * m(7) - t557, qJD(5) * t269 + qJD(6) * t391 - t554; -qJD(2) * t1 + qJD(3) * t12 + qJD(4) * t10 - qJD(6) * t8 - t666, qJD(6) * t32 - t522, -t266 * qJD(6) + (qJD(4) * t731 + t610) * m(7) + t556, m(7) * t545 + t267 * qJD(6) + t557, -t139 * qJD(6) (pkin(10) * t424 + t759) * qJD(6) - t515; -qJD(2) * t5 - qJD(3) * t26 - qJD(4) * t24 + qJD(5) * t8 - t656, -qJD(3) * t63 - qJD(5) * t32 - t527, qJD(5) * t266 - t553, -qJD(5) * t267 + t554, t515, 0;];
Cq  = t20;
