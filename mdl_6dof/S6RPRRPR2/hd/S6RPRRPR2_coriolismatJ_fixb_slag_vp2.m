% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:31
% EndTime: 2019-03-09 05:00:53
% DurationCPUTime: 14.58s
% Computational Cost: add. (27955->719), mult. (58038->983), div. (0->0), fcn. (61612->10), ass. (0->353)
t455 = sin(pkin(11));
t458 = sin(qJ(4));
t461 = cos(qJ(4));
t584 = cos(pkin(11));
t513 = t584 * t461;
t401 = -t455 * t458 + t513;
t402 = -t455 * t461 - t458 * t584;
t457 = sin(qJ(6));
t460 = cos(qJ(6));
t316 = t401 * t457 - t402 * t460;
t508 = t460 * t401 + t402 * t457;
t557 = Ifges(7,5) * t508 - Ifges(7,6) * t316;
t626 = -qJ(5) - pkin(8);
t415 = t626 * t458;
t417 = t626 * t461;
t335 = t455 * t415 - t584 * t417;
t267 = -pkin(9) * t401 - t335;
t708 = t584 * t415 + t455 * t417;
t724 = pkin(9) * t402 + t708;
t155 = t267 * t460 - t457 * t724;
t731 = t267 * t457 + t460 * t724;
t755 = t155 * mrSges(7,1) - t731 * mrSges(7,2);
t26 = t557 + t755;
t757 = t26 * qJD(6);
t459 = sin(qJ(3));
t571 = t458 * t459;
t374 = t455 * t571 - t459 * t513;
t375 = t402 * t459;
t259 = t374 * t460 - t375 * t457;
t509 = t374 * t457 + t460 * t375;
t756 = -t259 * mrSges(7,1) + t509 * mrSges(7,2);
t749 = qJD(6) * t756;
t462 = cos(qJ(3));
t617 = Ifges(7,4) * t259;
t136 = Ifges(7,2) * t509 - Ifges(7,6) * t462 - t617;
t243 = Ifges(7,4) * t509;
t138 = -Ifges(7,1) * t259 - t462 * Ifges(7,5) + t243;
t212 = mrSges(7,2) * t462 + mrSges(7,3) * t509;
t602 = t259 * mrSges(7,3);
t214 = -mrSges(7,1) * t462 + t602;
t442 = sin(pkin(10)) * pkin(1) + pkin(7);
t423 = t459 * t442;
t384 = pkin(4) * t571 + t423;
t288 = -pkin(5) * t375 + t384;
t444 = -pkin(4) * t461 - pkin(3);
t352 = -pkin(5) * t401 + t444;
t685 = t155 / 0.2e1;
t723 = Ifges(7,2) * t259 + t243;
t725 = t316 * mrSges(7,1);
t735 = t508 * mrSges(7,2) + t725;
t736 = Ifges(7,1) * t509 + t617;
t753 = t731 / 0.2e1;
t754 = t214 * t685 + (t138 / 0.4e1 + t723 / 0.4e1) * t508 + t212 * t753 + t756 * t352 / 0.2e1 + t735 * t288 / 0.2e1 + (t736 / 0.4e1 - t136 / 0.4e1) * t316;
t531 = -cos(pkin(10)) * pkin(1) - pkin(2);
t398 = -pkin(3) * t462 - t459 * pkin(8) + t531;
t378 = t461 * t398;
t569 = t459 * t461;
t540 = qJ(5) * t569;
t265 = -t540 + t378 + (-t442 * t458 - pkin(4)) * t462;
t424 = t462 * t442;
t328 = t458 * t398 + t424 * t461;
t287 = -qJ(5) * t571 + t328;
t271 = t455 * t287;
t159 = t584 * t265 - t271;
t634 = t374 * pkin(9);
t111 = -pkin(5) * t462 + t159 + t634;
t515 = t584 * t287;
t160 = t455 * t265 + t515;
t636 = pkin(9) * t375;
t115 = t160 + t636;
t60 = t111 * t457 + t115 * t460;
t327 = -t424 * t458 + t378;
t286 = t327 - t540;
t163 = -t286 * t455 - t515;
t118 = t163 - t636;
t164 = t584 * t286 - t271;
t119 = t164 + t634;
t64 = t118 * t460 - t119 * t457;
t752 = t64 + t60;
t595 = t316 * mrSges(7,3);
t746 = t288 * t756;
t745 = t352 * t735;
t633 = t459 * pkin(3);
t637 = pkin(8) * t462;
t427 = t633 - t637;
t340 = t458 * t423 + t461 * t427;
t565 = t461 * t462;
t298 = t459 * pkin(4) - qJ(5) * t565 + t340;
t341 = t458 * t427 - t442 * t569;
t570 = t458 * t462;
t318 = -qJ(5) * t570 + t341;
t175 = t584 * t298 - t318 * t455;
t176 = t455 * t298 + t584 * t318;
t376 = t402 * t462;
t377 = t401 * t462;
t261 = t376 * t460 - t377 * t457;
t213 = -mrSges(7,2) * t459 + mrSges(7,3) * t261;
t264 = t376 * t457 + t377 * t460;
t215 = mrSges(7,1) * t459 - mrSges(7,3) * t264;
t337 = -mrSges(6,2) * t459 + mrSges(6,3) * t376;
t339 = mrSges(6,1) * t459 - mrSges(6,3) * t377;
t530 = t584 * pkin(4);
t443 = t530 + pkin(5);
t639 = pkin(4) * t455;
t382 = t443 * t460 - t457 * t639;
t383 = t443 * t457 + t460 * t639;
t487 = -Ifges(6,5) * t377 / 0.2e1 - Ifges(6,6) * t376 / 0.2e1;
t486 = -Ifges(7,5) * t264 / 0.2e1 - Ifges(7,6) * t261 / 0.2e1;
t644 = t459 / 0.2e1;
t121 = pkin(5) * t459 - pkin(9) * t377 + t175;
t140 = pkin(9) * t376 + t176;
t68 = t121 * t460 - t140 * t457;
t69 = t121 * t457 + t140 * t460;
t503 = Ifges(7,3) * t644 - t69 * mrSges(7,2) / 0.2e1 + t68 * mrSges(7,1) / 0.2e1 - t486;
t522 = t565 / 0.2e1;
t523 = -t570 / 0.2e1;
t689 = m(6) * pkin(4);
t547 = t689 / 0.2e1;
t653 = t383 / 0.2e1;
t688 = -mrSges(5,2) / 0.2e1;
t691 = m(7) / 0.2e1;
t718 = Ifges(5,3) + Ifges(6,3);
t744 = (t382 * t68 + t383 * t69) * t691 + t175 * mrSges(6,1) / 0.2e1 - t176 * mrSges(6,2) / 0.2e1 + t340 * mrSges(5,1) / 0.2e1 + t341 * t688 + t382 * t215 / 0.2e1 + t213 * t653 + (t175 * t584 + t176 * t455) * t547 + Ifges(5,5) * t522 + Ifges(5,6) * t523 + t337 * t639 / 0.2e1 + t339 * t530 / 0.2e1 - t487 + t718 * t644 + t503;
t601 = t264 * mrSges(7,2);
t605 = t261 * mrSges(7,1);
t561 = t605 / 0.2e1 - t601 / 0.2e1;
t592 = t377 * mrSges(6,2);
t593 = t376 * mrSges(6,1);
t743 = t561 + t593 / 0.2e1 - t592 / 0.2e1;
t616 = Ifges(7,4) * t316;
t187 = Ifges(7,2) * t508 + t616;
t737 = Ifges(7,1) * t508 - t616;
t742 = t187 / 0.4e1 - t737 / 0.4e1;
t538 = t725 / 0.2e1;
t712 = t508 * mrSges(7,3);
t719 = t509 / 0.2e1;
t740 = t712 * t719;
t739 = t160 + t163;
t438 = pkin(4) * t570;
t385 = t424 + t438;
t289 = -pkin(5) * t376 + t385;
t692 = -m(7) / 0.2e1;
t694 = -m(6) / 0.2e1;
t738 = -t289 * t692 - t385 * t694 - t743;
t714 = Ifges(7,5) * t509;
t728 = Ifges(7,6) * t259;
t559 = t714 + t728;
t307 = Ifges(7,4) * t508;
t190 = Ifges(7,1) * t316 + t307;
t722 = -Ifges(7,2) * t316 + t307;
t733 = t722 / 0.4e1 + t190 / 0.4e1;
t520 = t728 / 0.2e1 + t714 / 0.2e1;
t674 = -t259 / 0.2e1;
t730 = t259 / 0.2e1;
t669 = -t316 / 0.2e1;
t575 = t382 * t509;
t729 = mrSges(7,3) * t575;
t590 = t402 * mrSges(6,3);
t489 = -t340 * t458 + t341 * t461;
t720 = t508 / 0.2e1;
t450 = Ifges(5,5) * t461;
t613 = Ifges(5,6) * t458;
t710 = Ifges(4,4) - t450 / 0.2e1 + t613 / 0.2e1;
t321 = -t402 * mrSges(6,1) + t401 * mrSges(6,2);
t709 = t735 + t321;
t451 = Ifges(5,4) * t461;
t707 = -Ifges(5,2) * t458 + t451;
t421 = Ifges(5,1) * t458 + t451;
t453 = t458 ^ 2;
t454 = t461 ^ 2;
t551 = t453 + t454;
t354 = t374 * mrSges(6,1);
t512 = t375 * mrSges(6,2) - t354;
t482 = t756 + t512;
t706 = t335 * t374 / 0.2e1 - t708 * t375 / 0.2e1;
t620 = Ifges(5,4) * t458;
t419 = Ifges(5,2) * t461 + t620;
t396 = t459 * t419;
t412 = -mrSges(5,1) * t462 - mrSges(5,3) * t569;
t690 = -pkin(8) / 0.2e1;
t705 = t412 * t690 - t396 / 0.4e1;
t397 = t459 * t421;
t410 = mrSges(5,2) * t462 - mrSges(5,3) * t571;
t704 = t410 * t690 - t397 / 0.4e1;
t703 = Ifges(6,5) * t401 + Ifges(6,6) * t402 + t557;
t702 = Ifges(6,5) * t375 + Ifges(6,6) * t374 + t559;
t701 = -mrSges(5,1) * t461 + mrSges(5,2) * t458;
t700 = t164 / 0.2e1 - t159 / 0.2e1;
t699 = t595 * t674 + t740;
t146 = -mrSges(7,1) * t509 - mrSges(7,2) * t259;
t184 = -mrSges(7,1) * t508 + mrSges(7,2) * t316;
t619 = Ifges(6,4) * t374;
t253 = Ifges(6,2) * t375 - t462 * Ifges(6,6) - t619;
t278 = Ifges(6,1) * t375 + t619;
t394 = Ifges(6,4) * t401;
t323 = Ifges(6,2) * t402 + t394;
t618 = Ifges(6,4) * t402;
t324 = Ifges(6,2) * t401 - t618;
t325 = Ifges(6,1) * t401 + t618;
t326 = -Ifges(6,1) * t402 + t394;
t546 = pkin(4) * t569;
t331 = -pkin(5) * t374 + t546;
t638 = pkin(4) * t458;
t358 = -pkin(5) * t402 + t638;
t59 = t111 * t460 - t115 * t457;
t641 = -t462 / 0.4e1;
t647 = t444 / 0.2e1;
t65 = t118 * t457 + t119 * t460;
t659 = t358 / 0.2e1;
t338 = -mrSges(6,1) * t462 + t374 * mrSges(6,3);
t661 = t338 / 0.2e1;
t594 = t375 * mrSges(6,3);
t336 = t462 * mrSges(6,2) + t594;
t662 = t336 / 0.2e1;
t663 = t331 / 0.2e1;
t686 = -t731 / 0.2e1;
t697 = t703 * t641 + (t323 + t326) * t375 / 0.4e1 - t335 * t661 + t708 * t662 + (t253 / 0.4e1 - t278 / 0.4e1) * t402 + (t324 / 0.4e1 - t325 / 0.4e1) * t374 + (t288 * t358 + t331 * t352 + t752 * t731 + (t59 - t65) * t155) * t691 + (-t155 * t730 + t509 * t686 + t65 * t720 - t59 * t508 / 0.2e1 + t752 * t669) * mrSges(7,3) + t742 * t259 + t512 * t647 + t384 * t321 / 0.2e1 + t146 * t659 + t184 * t663 + t733 * t509 + t754;
t696 = 2 * qJD(3);
t695 = m(5) / 0.2e1;
t693 = m(6) / 0.2e1;
t677 = t261 / 0.2e1;
t673 = t264 / 0.2e1;
t275 = -mrSges(6,1) * t375 - mrSges(6,2) * t374;
t672 = t275 / 0.2e1;
t666 = t316 / 0.2e1;
t658 = -t374 / 0.2e1;
t657 = t375 / 0.2e1;
t655 = t376 / 0.2e1;
t654 = t377 / 0.2e1;
t649 = -t419 / 0.4e1;
t422 = Ifges(5,1) * t461 - t620;
t648 = t422 / 0.4e1;
t646 = -t458 / 0.2e1;
t645 = t458 / 0.2e1;
t643 = t461 / 0.2e1;
t642 = -t462 / 0.2e1;
t640 = t462 / 0.2e1;
t632 = t59 * mrSges(7,2);
t631 = t60 * mrSges(7,1);
t630 = t64 * mrSges(7,1);
t629 = t65 * mrSges(7,2);
t623 = mrSges(5,3) * t459;
t622 = mrSges(7,3) * t383;
t591 = t401 * mrSges(6,3);
t589 = t458 * mrSges(5,1);
t588 = t461 * mrSges(5,2);
t587 = t462 * Ifges(5,5);
t586 = t462 * Ifges(5,6);
t585 = -mrSges(4,1) + t701;
t18 = t212 * t719 + t214 * t730 + t756 * t642 + (-t259 ^ 2 / 0.2e1 - t509 ^ 2 / 0.2e1) * mrSges(7,3);
t582 = t18 * qJD(1);
t581 = t259 * t214;
t580 = t509 * t212;
t577 = t374 * t338;
t576 = t375 * t336;
t574 = t384 * t458;
t370 = t459 * t707 - t586;
t573 = t458 * t370;
t413 = t459 * mrSges(5,1) - mrSges(5,3) * t565;
t572 = t458 * t413;
t568 = t459 * t462;
t372 = t422 * t459 - t587;
t567 = t461 * t372;
t411 = -t459 * mrSges(5,2) - mrSges(5,3) * t570;
t566 = t461 * t411;
t418 = t588 + t589;
t395 = t462 * t418;
t550 = qJD(3) * t462;
t110 = -mrSges(7,1) * t383 - mrSges(7,2) * t382;
t549 = t110 * qJD(6);
t548 = mrSges(6,3) * t639;
t543 = t638 / 0.2e1;
t539 = t602 / 0.2e1;
t537 = t590 / 0.2e1;
t532 = t374 * t537 + t699;
t525 = t442 * t418 / 0.2e1;
t521 = t735 * t642;
t517 = t454 / 0.2e1 + t453 / 0.2e1;
t510 = t450 - t613;
t506 = t546 / 0.2e1;
t505 = mrSges(6,3) * t530;
t137 = Ifges(7,4) * t264 + Ifges(7,2) * t261 + t459 * Ifges(7,6);
t139 = Ifges(7,1) * t264 + Ifges(7,4) * t261 + t459 * Ifges(7,5);
t147 = t601 - t605;
t254 = Ifges(6,4) * t377 + Ifges(6,2) * t376 + t459 * Ifges(6,6);
t357 = Ifges(6,4) * t375;
t255 = -Ifges(6,1) * t374 - Ifges(6,5) * t462 + t357;
t256 = Ifges(6,1) * t377 + Ifges(6,4) * t376 + t459 * Ifges(6,5);
t276 = t592 - t593;
t371 = Ifges(5,6) * t459 + t462 * t707;
t373 = Ifges(5,5) * t459 + t422 * t462;
t4 = (t531 * mrSges(4,2) + t567 / 0.2e1 - t573 / 0.2e1 + t486 + t487 + t710 * t462) * t462 + t68 * t214 + t59 * t215 + t69 * t212 + t60 * t213 + (t373 * t643 + t371 * t646 + t442 * t395 + Ifges(7,5) * t674 + Ifges(7,6) * t719 + Ifges(6,5) * t658 + Ifges(6,6) * t657 + t531 * mrSges(4,1) - t710 * t459 + (-Ifges(4,2) + Ifges(4,1) - Ifges(7,3) + (m(5) * t442 + t418) * t442 - t718) * t462) * t459 + t137 * t719 + m(5) * (t327 * t340 + t328 * t341) + t288 * t147 + t289 * t146 + m(6) * (t159 * t175 + t160 * t176 + t384 * t385) + m(7) * (t288 * t289 + t59 * t68 + t60 * t69) + t255 * t654 + t253 * t655 + t254 * t657 + t176 * t336 + t160 * t337 + t175 * t338 + t159 * t339 + t384 * t276 + t385 * t275 + t341 * t410 + t328 * t411 + t340 * t412 + t327 * t413 + t256 * t658 + t138 * t673 + t139 * t674 + t136 * t677;
t9 = t212 * t673 + t213 * t674 + t214 * t677 + t215 * t719 + t336 * t654 + t337 * t658 + t338 * t655 + t339 * t657 + (t410 * t643 + t412 * t646 - t147 / 0.2e1 - t395 / 0.2e1 - t276 / 0.2e1) * t462 + (t566 / 0.2e1 - t572 / 0.2e1 + t146 / 0.2e1 + t672 + (t589 / 0.2e1 + t588 / 0.2e1) * t459) * t459 + (t159 * t376 + t160 * t377 + t175 * t375 - t176 * t374 + t384 * t459 - t385 * t462) * t693 + (-t259 * t69 + t59 * t261 + t60 * t264 + t288 * t459 - t289 * t462 + t509 * t68) * t691 + ((-t327 * t458 + t328 * t461 - t424) * t462 + (t489 + t423) * t459) * t695;
t500 = t4 * qJD(1) + t9 * qJD(2);
t499 = t59 * t259 + t509 * t60;
t490 = t159 * t374 + t160 * t375;
t10 = -t580 / 0.2e1 - t581 / 0.2e1 - t576 / 0.2e1 - t577 / 0.2e1 + (t259 * t730 + t509 * t719) * mrSges(7,3) + (t374 ^ 2 / 0.2e1 + t375 ^ 2 / 0.2e1) * mrSges(6,3) + (t163 * t375 - t164 * t374 + t490) * t694 + (-t259 * t65 - t331 * t462 + t509 * t64 + t499) * t692 + (t410 * t645 + t412 * t643 + t517 * t623 + t522 * t689) * t459 + (-t459 * t701 + t482) * t640;
t277 = Ifges(6,2) * t374 + t357;
t5 = t331 * t146 + t136 * t730 + t64 * t214 + t65 * t212 + t736 * t674 + t746 + t164 * t336 + t163 * t338 - t384 * t354 + t327 * t410 - t328 * t412 + (t259 * t60 - t509 * t59) * mrSges(7,3) + m(6) * (t159 * t163 + t160 * t164) + m(7) * (t288 * t331 + t59 * t64 + t60 * t65) + (-t159 * mrSges(6,3) + t255 / 0.2e1 + t277 / 0.2e1 + t384 * mrSges(6,2)) * t375 + (t160 * mrSges(6,3) + t253 / 0.2e1 - t278 / 0.2e1) * t374 + ((t327 * mrSges(5,3) + t396 / 0.2e1 - t372 / 0.2e1 + t587 / 0.2e1 - mrSges(5,2) * t423) * t458 + (-t328 * mrSges(5,3) - t397 / 0.2e1 - t370 / 0.2e1 + t586 / 0.2e1 + mrSges(5,1) * t423 + (m(6) * t384 + t275) * pkin(4)) * t461) * t459 + t702 * t642 + (t723 + t138) * t719;
t498 = t5 * qJD(1) - t10 * qJD(2);
t8 = t59 * t212 - t60 * t214 + t559 * t642 + t746 - (-t60 * mrSges(7,3) + t736 / 0.2e1 - t136 / 0.2e1) * t259 + (-t59 * mrSges(7,3) + t138 / 0.2e1 + t723 / 0.2e1) * t509;
t497 = t8 * qJD(1) + t18 * qJD(2);
t43 = m(6) * (-t374 * t377 + t375 * t376 - t568) + m(7) * (-t259 * t264 + t261 * t509 - t568) + m(5) * (-0.1e1 + t551) * t568;
t496 = t9 * qJD(1) + t43 * qJD(2);
t495 = t10 * qJD(1);
t27 = m(6) * t490 + m(7) * t499 + t576 + t577 + t580 + t581;
t494 = qJD(1) * t27;
t477 = (t374 * t584 + t375 * t455) * t689;
t484 = m(7) * (t259 * t382 + t383 * t509);
t470 = t484 / 0.2e1 + t477 / 0.2e1;
t481 = m(6) * t506 + m(7) * t663;
t42 = -t470 + t481 + t482;
t469 = (-t316 * t382 + t383 * t508) * t691 + (t401 * t455 + t402 * t584) * t547;
t483 = m(6) * t543 + m(7) * t659;
t45 = -t469 + t483 + t709;
t493 = qJD(1) * t42 + qJD(3) * t45;
t79 = 0.2e1 * t720 * mrSges(7,2) + 0.2e1 * t538;
t492 = qJD(1) * t756 + qJD(3) * t79;
t322 = -mrSges(6,1) * t401 - mrSges(6,2) * t402;
t13 = t187 * t669 + t737 * t666 + t745 - pkin(3) * t418 + (t707 / 0.2e1 + t421 / 0.2e1) * t461 + (t324 / 0.2e1 - t325 / 0.2e1) * t402 + (t323 / 0.2e1 + t326 / 0.2e1) * t401 + (pkin(4) * t322 - t419 / 0.2e1 + t422 / 0.2e1) * t458 + (t722 + t190) * t720 + (m(6) * t638 + t321) * t444 + (m(7) * t352 + t184) * t358;
t466 = -t438 * t694 - t358 * t462 * t692 + t395 / 0.2e1 + t740 - t316 * t539 + t709 * t640;
t468 = (t261 * t382 + t264 * t383) * t691 + (t376 * t584 + t377 * t455) * t547 + mrSges(5,1) * t523 + t565 * t688 + t743;
t16 = t466 + t468 - t699;
t475 = t739 * t708 + (-t159 + t164) * t335;
t2 = t551 * t623 * t690 + t697 - (t421 + t707) * t571 / 0.4e1 + (t277 + t255) * t401 / 0.4e1 + t510 * t641 + t275 * t543 + t322 * t506 + t701 * t633 / 0.2e1 + t569 * t648 + t569 * t649 + t567 / 0.4e1 - t573 / 0.4e1 + t459 * t525 + ((t444 * t569 + t574) * pkin(4) + t475) * t693 + t706 * mrSges(6,3) + t704 * t458 + t705 * t461 + t739 * t537 + t700 * t591 - t744;
t480 = t2 * qJD(1) - t16 * qJD(2) + t13 * qJD(3);
t21 = t745 + (t737 / 0.2e1 - t187 / 0.2e1) * t316 + (t190 / 0.2e1 + t722 / 0.2e1) * t508;
t29 = t521 - t561;
t465 = -(mrSges(7,3) * t685 - t742) * t259 + (mrSges(7,3) * t686 + t733) * t509 + t557 * t641 + t754;
t6 = t465 - t503;
t479 = t6 * qJD(1) + t29 * qJD(2) + t21 * qJD(3);
t467 = (t594 / 0.2e1 + t662) * t401 + (t159 * t402 + t160 * t401 + t335 * t375 + t374 * t708) * t693 + (-t155 * t509 + t259 * t731 - t316 * t59 + t508 * t60) * t691 + t214 * t669 + t212 * t720 + t402 * t661 + t532;
t15 = t467 - t738;
t34 = (t316 ^ 2 + t508 ^ 2) * mrSges(7,3) + (t401 ^ 2 + t402 ^ 2) * mrSges(6,3) + m(7) * (-t155 * t508 - t316 * t731) + m(6) * (t335 * t401 + t402 * t708);
t471 = (-t374 * t401 + t375 * t402) * t693 + (-t259 * t508 - t316 * t509) * t691;
t56 = (t694 + t692) * t459 + t471;
t478 = -qJD(1) * t15 - qJD(2) * t56 - qJD(3) * t34;
t476 = -t382 * t212 / 0.2e1 + t214 * t653 - t520;
t11 = (-t259 * t653 + t575 / 0.2e1) * mrSges(7,3) + (-t65 / 0.2e1 + t59 / 0.2e1) * mrSges(7,2) + (t64 / 0.2e1 + t60 / 0.2e1) * mrSges(7,1) + t476 + t520;
t25 = (t686 + t753) * mrSges(7,2) + (t685 - t155 / 0.2e1) * mrSges(7,1);
t474 = t11 * qJD(1) - t25 * qJD(3) - t110 * qJD(4);
t80 = t538 - t725 / 0.2e1;
t74 = t469 + t483;
t61 = t470 + t481;
t55 = t471 + (m(6) + m(7)) * t644;
t30 = t521 + t561;
t17 = t590 * t658 - t466 + t468 + t532;
t14 = t467 + t738;
t12 = -t632 / 0.2e1 - t631 / 0.2e1 + t383 * t539 - t729 / 0.2e1 - t629 / 0.2e1 + t630 / 0.2e1 - t476 + t520;
t7 = t465 + t503;
t3 = qJD(3) * t9 - qJD(4) * t10 + qJD(6) * t18;
t1 = ((t163 / 0.2e1 + t160 / 0.2e1) * t402 + t700 * t401 + t706) * mrSges(6,3) + (pkin(4) * t672 - t370 / 0.4e1 + t586 / 0.4e1 + t704) * t458 + (t525 - t517 * pkin(8) * mrSges(5,3) + (-t421 / 0.4e1 - t707 / 0.4e1 + pkin(3) * mrSges(5,2) / 0.2e1) * t458 + (t648 + t649 - pkin(3) * mrSges(5,1) / 0.2e1 + (t322 / 0.2e1 + m(6) * t647) * pkin(4)) * t461) * t459 + (t277 / 0.4e1 + t255 / 0.4e1) * t401 + (pkin(4) * t574 + t475) * t693 + (t372 / 0.4e1 + t705) * t461 + t450 * t641 + t697 + t744;
t19 = [qJD(3) * t4 + qJD(4) * t5 + qJD(5) * t27 + qJD(6) * t8, t3, t1 * qJD(4) + t14 * qJD(5) + t7 * qJD(6) + (t419 * t646 + t421 * t643 + t442 * t585 + Ifges(4,5)) * t550 + ((-pkin(3) * t424 + pkin(8) * t489) * t695 + (t175 * t708 + t176 * t335 + t385 * t444) * t693 + (-t155 * t69 + t289 * t352 + t68 * t731) * t691) * t696 + t500 + (t489 * mrSges(5,3) + t69 * t712 + t137 * t720 + t708 * t339 + t731 * t215 + t371 * t643 + mrSges(4,2) * t423 + t175 * t590 + t176 * t591 + t289 * t184 - t155 * t213 - t68 * t595 + t373 * t645 + t326 * t654 + t324 * t655 + t335 * t337 + t352 * t147 + t385 * t322 - pkin(3) * t395 + t401 * t254 / 0.2e1 - t402 * t256 / 0.2e1 + t139 * t666 + t190 * t673 + t187 * t677 + t444 * t276 - Ifges(4,6) * t459 + (t566 - t572) * pkin(8) + (Ifges(5,5) * t458 - Ifges(6,5) * t402 + Ifges(7,5) * t316 + Ifges(5,6) * t461 + Ifges(6,6) * t401 + Ifges(7,6) * t508) * t644) * qJD(3), t1 * qJD(3) + (t374 * t548 - t375 * t505 - t629 + t630 + t259 * t622 - t729 + m(7) * (t382 * t64 + t383 * t65) + t163 * mrSges(6,1) + (t163 * t584 + t164 * t455) * t689 - t164 * mrSges(6,2) - t327 * mrSges(5,2) - t328 * mrSges(5,1) - Ifges(5,5) * t571 - Ifges(5,6) * t569 + t702) * qJD(4) + t61 * qJD(5) + t12 * qJD(6) + t498, qJD(3) * t14 + qJD(4) * t61 + t494, t7 * qJD(3) + t12 * qJD(4) + (t559 - t631 - t632) * qJD(6) + t497; t3, t43 * qJD(3), t17 * qJD(4) + t55 * qJD(5) + t30 * qJD(6) + (mrSges(5,3) * t551 - mrSges(4,2)) * t550 + ((t335 * t377 + t376 * t708 + t444 * t459) * t693 + (-t155 * t264 + t261 * t731 + t352 * t459) * t691 + (t551 * t637 - t633) * t695) * t696 + t496 + (-t261 * t595 + t264 * t712 + t376 * t590 + t377 * t591 + (t184 + t322 + t585) * t459) * qJD(3), t17 * qJD(3) + (-mrSges(5,1) * t569 + mrSges(5,2) * t571 + t477 - t482 + t484) * qJD(4) - t749 - t495, qJD(3) * t55, t30 * qJD(3) - qJD(4) * t756 + t582 - t749; qJD(4) * t2 + qJD(5) * t15 + qJD(6) * t6 - t500, -qJD(4) * t16 + qJD(5) * t56 + qJD(6) * t29 - t496, qJD(4) * t13 + qJD(5) * t34 + qJD(6) * t21 (m(7) * (t155 * t382 + t383 * t731) - t708 * mrSges(6,2) - t335 * mrSges(6,1) + (-t335 * t584 + t455 * t708) * t689 + t402 * t548 - t401 * t505 - t316 * t622 - t382 * t712 + t510 + t701 * pkin(8) + t703 + t755) * qJD(4) + t74 * qJD(5) + t757 + t480, qJD(4) * t74 + qJD(6) * t80 - t478, t26 * qJD(4) + t80 * qJD(5) + t479 + t757; -qJD(3) * t2 - qJD(5) * t42 - qJD(6) * t11 - t498, qJD(3) * t16 + t495, -qJD(5) * t45 + qJD(6) * t25 - t480, t549, -t493, -t474 + t549; -qJD(3) * t15 + qJD(4) * t42 - t494 + t749, -qJD(3) * t56, qJD(4) * t45 + qJD(6) * t79 + t478, t493, 0, t492; -qJD(3) * t6 + qJD(4) * t11 - qJD(5) * t756 - t497, -t29 * qJD(3) - t582, -qJD(4) * t25 - qJD(5) * t79 - t479, t474, -t492, 0;];
Cq  = t19;
