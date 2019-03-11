% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:41
% EndTime: 2019-03-09 13:13:12
% DurationCPUTime: 22.11s
% Computational Cost: add. (61308->636), mult. (117294->810), div. (0->0), fcn. (145067->10), ass. (0->379)
t549 = sin(pkin(11));
t550 = cos(pkin(11));
t618 = sin(qJ(2));
t621 = cos(qJ(2));
t346 = -t549 * t618 + t550 * t621;
t415 = -t549 * t621 - t550 * t618;
t617 = sin(qJ(4));
t620 = cos(qJ(4));
t323 = -t617 * t346 + t415 * t620;
t368 = sin(qJ(5));
t619 = cos(qJ(5));
t405 = t617 * t415;
t702 = t620 * t346 + t405;
t273 = t368 * t323 + t619 * t702;
t367 = sin(qJ(6));
t369 = cos(qJ(6));
t460 = mrSges(7,1) * t369 - t367 * mrSges(7,2);
t516 = t618 * pkin(7);
t432 = -qJ(3) * t618 - t516;
t518 = t621 * pkin(7);
t433 = qJ(3) * t621 + t518;
t329 = t549 * t432 + t433 * t550;
t385 = t346 * pkin(8) + t329;
t667 = t550 * t432 - t433 * t549;
t699 = pkin(8) * t415 + t667;
t214 = t385 * t620 + t617 * t699;
t193 = pkin(9) * t702 + t214;
t748 = -t385 * t617 + t620 * t699;
t765 = pkin(9) * t323 + t748;
t791 = t193 * t619 + t368 * t765;
t795 = t791 * t460;
t792 = -t368 * t193 + t619 * t765;
t798 = t792 * mrSges(6,2);
t799 = t791 * mrSges(6,1);
t418 = -t798 / 0.2e1 - t795 / 0.2e1 - t799 / 0.2e1;
t599 = Ifges(7,4) * t369;
t457 = Ifges(7,1) * t367 + t599;
t434 = t369 * t457;
t600 = Ifges(7,4) * t367;
t454 = Ifges(7,2) * t369 + t600;
t435 = t367 * t454;
t674 = t435 / 0.4e1 - t434 / 0.4e1;
t349 = Ifges(7,5) * t367 + t369 * Ifges(7,6);
t718 = -t323 * t619 + t368 * t702;
t743 = t718 * t349;
t745 = Ifges(6,6) * t718;
t455 = -Ifges(7,2) * t367 + t599;
t588 = Ifges(7,6) * t718;
t752 = t273 * t455 + t588;
t772 = t369 * t752;
t458 = t369 * Ifges(7,1) - t600;
t595 = Ifges(7,5) * t718;
t751 = t273 * t458 + t595;
t773 = t367 * t751;
t807 = t674 * t273 - t418 - t743 / 0.4e1 + t745 / 0.2e1 - t772 / 0.4e1 - t773 / 0.4e1;
t719 = -t795 - t799 - t798 + t743 / 0.2e1 - t745 + t772 / 0.2e1 + t773 / 0.2e1;
t728 = Ifges(5,6) * t323;
t729 = Ifges(5,5) * t702;
t758 = Ifges(6,5) * t273;
t775 = t748 * mrSges(5,2);
t786 = t214 * mrSges(5,1);
t806 = t719 + t728 + t729 + t758 - t775 - t786;
t767 = -t758 / 0.2e1 + t807;
t361 = -pkin(2) * t621 - pkin(1);
t331 = -t346 * pkin(3) + t361;
t280 = -pkin(4) * t702 + t331;
t805 = m(6) * t280 - mrSges(6,1) * t273 + mrSges(6,2) * t718;
t710 = t368 * t792;
t803 = -t619 * t791 + t710;
t489 = t550 * pkin(2);
t358 = t489 + pkin(3);
t494 = pkin(2) * t549;
t335 = t620 * t358 - t494 * t617;
t334 = pkin(4) + t335;
t336 = t358 * t617 + t494 * t620;
t528 = t368 * t336;
t304 = t334 * t619 - t528;
t490 = t619 * t336;
t305 = t368 * t334 + t490;
t712 = t305 * t792;
t802 = -t304 * t791 + t712;
t120 = -pkin(5) * t273 - pkin(10) * t718 + t280;
t57 = t120 * t369 - t367 * t791;
t58 = t367 * t120 + t369 * t791;
t565 = t369 * mrSges(7,2);
t567 = t367 * mrSges(7,1);
t459 = t565 + t567;
t174 = t459 * t718;
t755 = t459 * t273;
t721 = t791 * t174 - t792 * t755;
t566 = t367 * mrSges(7,3);
t687 = mrSges(7,2) * t718;
t749 = -t273 * t566 - t687;
t564 = t369 * mrSges(7,3);
t688 = mrSges(7,1) * t718;
t750 = -t273 * t564 + t688;
t801 = t57 * t750 + t58 * t749 + t721;
t800 = pkin(5) * t791;
t517 = t619 * pkin(4);
t360 = -t517 - pkin(5);
t681 = t360 * t791;
t711 = t367 * t792;
t709 = t369 * t792;
t794 = t791 * t792;
t793 = m(7) * t791 + t755;
t789 = t728 / 0.2e1 + t729 / 0.2e1 - t775 / 0.2e1 - t786 / 0.2e1;
t771 = t749 * t369;
t788 = t771 / 0.2e1;
t787 = mrSges(5,3) * t748;
t530 = t367 * t750;
t782 = t771 - t530;
t322 = Ifges(5,4) * t702;
t715 = t323 / 0.2e1;
t692 = Ifges(5,2) * t715 + t322 / 0.2e1;
t760 = -t273 / 0.2e1;
t777 = mrSges(6,2) * t760;
t776 = Ifges(7,5) * t760;
t581 = t273 * mrSges(6,3);
t626 = t360 / 0.2e1;
t770 = t755 * t626;
t657 = pkin(5) / 0.2e1;
t769 = t755 * t657;
t302 = -pkin(5) - t304;
t303 = pkin(10) + t305;
t365 = t367 ^ 2;
t366 = t369 ^ 2;
t522 = t365 + t366;
t477 = t522 * t273;
t603 = mrSges(7,3) * t366;
t604 = mrSges(7,3) * t365;
t700 = t604 / 0.2e1 + t603 / 0.2e1;
t742 = t718 * t460;
t495 = -t742 / 0.2e1 + t700 * t273;
t655 = m(4) * pkin(2);
t659 = m(7) / 0.2e1;
t661 = m(6) / 0.2e1;
t663 = m(5) / 0.2e1;
t768 = (t323 * t335 + t336 * t702) * t663 + (t273 * t305 - t304 * t718) * t661 + (t302 * t718 + t303 * t477) * t659 + (t346 * t549 + t415 * t550) * t655 / 0.2e1 + t495;
t759 = t273 / 0.2e1;
t764 = (t273 / 0.4e1 + t759) * Ifges(7,6);
t362 = Ifges(7,5) * t369;
t586 = Ifges(7,6) * t367;
t675 = t362 - t586;
t746 = -t718 / 0.2e1;
t691 = Ifges(6,4) * t746;
t693 = t718 / 0.2e1;
t480 = t675 * t693 + t691 + (Ifges(6,2) + Ifges(7,3)) * t760;
t650 = Ifges(7,3) / 0.2e1;
t510 = t650 + Ifges(6,2) / 0.2e1;
t763 = Ifges(6,1) * t759 - t273 * t510 + t480;
t622 = t369 / 0.2e1;
t625 = -t367 / 0.2e1;
t761 = t751 * t622 + t752 * t625;
t757 = Ifges(7,3) * t746;
t756 = (-Ifges(6,5) / 0.2e1 + t674) * t273;
t704 = (t362 / 0.2e1 - t586 / 0.2e1) * t273;
t614 = pkin(5) * t718;
t186 = -t273 * pkin(10) + t614;
t726 = t702 * mrSges(5,2);
t727 = t323 * mrSges(5,1);
t479 = t726 - t727;
t587 = Ifges(7,6) * t273;
t112 = t455 * t718 - t587;
t594 = Ifges(7,5) * t273;
t115 = t458 * t718 - t594;
t670 = t280 * mrSges(6,2) + Ifges(6,1) * t693 + Ifges(6,4) * t273 + t112 * t625 + t115 * t622;
t747 = (t675 * t760 + t670) * t273 + t331 * t479 + (-t715 - t323 / 0.2e1) * Ifges(5,4) * t323 + (-Ifges(5,1) * t323 + t692) * t702;
t741 = t280 * mrSges(6,1) + t691;
t738 = -t726 / 0.2e1 + t727 / 0.2e1;
t736 = -t415 * mrSges(4,1) + t346 * mrSges(4,2);
t364 = t618 * pkin(2);
t733 = m(4) * t364;
t616 = t323 * pkin(4);
t634 = t302 / 0.2e1;
t722 = t755 * t634;
t713 = mrSges(6,1) * t746;
t669 = mrSges(7,3) * t522;
t708 = -mrSges(6,2) + t669;
t333 = -pkin(3) * t415 + t364;
t281 = t333 - t616;
t121 = t186 + t281;
t59 = t121 * t369 - t711;
t60 = t367 * t121 + t709;
t450 = -t59 * t367 + t60 * t369;
t676 = -t435 / 0.2e1 + t434 / 0.2e1;
t579 = t718 * mrSges(6,3);
t683 = -mrSges(6,1) - t460;
t682 = t302 * t791;
t677 = t174 + t579;
t673 = t788 - t530 / 0.2e1;
t515 = t718 * t566;
t178 = mrSges(7,2) * t273 - t515;
t523 = t369 * t178;
t182 = -mrSges(7,1) * t273 - t564 * t718;
t529 = t367 * t182;
t436 = t529 / 0.2e1 - t523 / 0.2e1;
t72 = t186 * t369 - t711;
t73 = t367 * t186 + t709;
t448 = -t72 * t367 + t73 * t369;
t624 = t367 / 0.2e1;
t666 = t455 * t622 + t458 * t624 + t676;
t538 = t273 * t369;
t180 = -mrSges(7,3) * t538 + t688;
t137 = t186 - t616;
t66 = t367 * t137 + t709;
t558 = t66 * t369;
t65 = t137 * t369 - t711;
t559 = t65 * t367;
t449 = t558 - t559;
t539 = t273 * t367;
t176 = -mrSges(7,3) * t539 - t687;
t525 = t369 * t176;
t660 = -m(7) / 0.2e1;
t665 = t449 * t660 - t525 / 0.2e1 + t180 * t624;
t662 = -m(6) / 0.2e1;
t658 = -pkin(5) / 0.2e1;
t656 = -pkin(10) / 0.2e1;
t654 = m(6) * pkin(4);
t653 = m(7) * pkin(4);
t652 = -mrSges(7,1) / 0.2e1;
t651 = mrSges(7,2) / 0.2e1;
t649 = -t72 / 0.2e1;
t648 = t73 / 0.2e1;
t647 = -t792 / 0.2e1;
t286 = t302 * t459;
t635 = -t286 / 0.2e1;
t633 = -t303 / 0.2e1;
t632 = -t304 / 0.2e1;
t307 = t335 * t619 - t528;
t630 = -t307 / 0.2e1;
t338 = t360 * t459;
t628 = -t338 / 0.2e1;
t615 = pkin(4) * t368;
t359 = pkin(10) + t615;
t627 = -t359 / 0.2e1;
t623 = -t369 / 0.2e1;
t580 = t273 * mrSges(6,2);
t267 = t718 * mrSges(6,1);
t388 = t415 ^ 2;
t3 = t405 * t692 + (t691 + t761 + t763) * t718 + m(7) * (t57 * t59 + t58 * t60 - t794) + t747 + (Ifges(3,1) - Ifges(3,2)) * t621 * t618 + t280 * t267 + t60 * t178 + t59 * t182 + (-t702 + t405) * t787 + (t733 + t736) * t361 - Ifges(4,4) * t388 - pkin(1) * (mrSges(3,1) * t618 + mrSges(3,2) * t621) + t333 * (-mrSges(5,1) * t405 - mrSges(5,2) * t323) + (Ifges(4,4) * t346 - mrSges(4,1) * t364 + (Ifges(4,2) - Ifges(4,1)) * t415 + (-t333 * mrSges(5,1) + t692 + t787) * t620) * t346 - t415 * mrSges(4,2) * t364 + (-t618 ^ 2 + t621 ^ 2) * Ifges(3,4) + m(5) * t331 * t333 + t801 + t805 * t281;
t577 = t3 * qJD(1);
t576 = t304 * mrSges(6,2);
t575 = t304 * mrSges(6,3);
t574 = t305 * mrSges(6,1);
t573 = t305 * mrSges(6,3);
t306 = t335 * t368 + t490;
t572 = t306 * mrSges(6,1);
t571 = t306 * t792;
t570 = t307 * mrSges(6,2);
t563 = t369 * t58;
t4 = t702 * (Ifges(5,2) * t323 + t322) / 0.2e1 + t58 * t176 + t66 * t178 + t57 * t180 + t65 * t182 + m(7) * (t57 * t65 + t58 * t66 - t794) + ((Ifges(7,1) * t538 + t595) * t622 + (-Ifges(7,2) * t539 + t588) * t625 - t538 * t600 + t741 + t763) * t718 - t805 * t616 + t721 + t747;
t562 = t4 * qJD(1);
t9 = t73 * t178 + t72 * t182 + m(7) * (t57 * t72 + t58 * t73 - t794) + (t480 + t741 + t761) * t718 - (t704 + (-Ifges(6,1) / 0.2e1 + t510) * t718 - t670) * t273 + t801;
t554 = t9 * qJD(1);
t552 = t792 * t718;
t451 = -t367 * t57 + t563;
t17 = t677 * t718 + (t323 ^ 2 + t702 ^ 2) * mrSges(5,3) + (t346 ^ 2 + t388) * mrSges(4,3) + (t523 - t529 + t581) * t273 + m(5) * (t214 * t702 + t323 * t748) + m(4) * (t329 * t346 + t415 * t667) + m(7) * (t273 * t451 - t552) + m(6) * (t273 * t791 - t552);
t548 = qJD(1) * t17;
t10 = t792 * t742 + t58 * t182 + (t115 * t624 + t112 * t622 + mrSges(7,3) * t563 + t349 * t760 + (t454 * t625 + t457 * t622) * t718) * t718 + (-t178 - t515) * t57;
t547 = t10 * qJD(1);
t383 = t333 * t663 + t281 * t661 + (t367 * t60 + t369 * t59) * t659 + t749 * t624 + t750 * t622 + t733 / 0.2e1;
t19 = t580 + t267 + t383 + t479 + t736 - t768;
t545 = t19 * qJD(1);
t377 = (t367 * t66 + t369 * t65) * t660 + t713 + t777 + t176 * t625 + t180 * t623 - t616 * t662 + t738;
t463 = -t267 / 0.2e1 + t495;
t474 = t522 * t359;
t498 = -t580 / 0.2e1;
t521 = t654 / 0.2e1;
t380 = (t273 * t474 + t360 * t718) * t659 + t498 + (t273 * t368 - t619 * t718) * t521 + t463 + t738;
t21 = t377 + t380;
t544 = t21 * qJD(1);
t396 = (t367 * t73 + t369 * t72) * t660 + t713 + t749 * t625 + t750 * t623;
t399 = (pkin(10) * t477 - t614) * t659 + t463;
t25 = t396 + t399 + 0.2e1 * t777;
t543 = t25 * qJD(1);
t439 = -t567 / 0.2e1 - t565 / 0.2e1;
t430 = t439 * t273;
t28 = t430 + t436;
t536 = t28 * qJD(1);
t535 = t305 * t460;
t534 = t306 * t460;
t533 = t359 * t367;
t520 = mrSges(7,3) * t559;
t519 = mrSges(7,3) * t558;
t514 = -t615 / 0.2e1;
t513 = t615 / 0.2e1;
t512 = t178 * t656;
t511 = t182 * t656;
t509 = t632 + t657;
t501 = t587 / 0.4e1;
t497 = -t566 / 0.2e1;
t496 = t564 / 0.2e1;
t493 = t367 * t619;
t492 = t369 * t619;
t486 = -t273 * t362 / 0.4e1;
t476 = t522 * t304;
t475 = t522 * t307;
t472 = -t517 / 0.2e1;
t471 = t517 / 0.2e1;
t469 = -t493 / 0.2e1;
t468 = (-Ifges(7,1) / 0.2e1 + Ifges(7,2) / 0.4e1) * t367;
t467 = (-Ifges(7,2) / 0.2e1 + Ifges(7,1) / 0.4e1) * t369;
t466 = -t112 / 0.4e1 + mrSges(7,1) * t647;
t465 = t115 / 0.4e1 + mrSges(7,2) * t647;
t464 = mrSges(7,3) * (-t365 / 0.2e1 - t366 / 0.2e1);
t461 = -t336 * mrSges(5,1) - t335 * mrSges(5,2);
t452 = t682 - t712;
t447 = t459 * t658 + t666;
t371 = (-t571 + t802) * t662 + (-t571 + t682) * t660 - t722 + t575 * t759 + t573 * t693 - t677 * t306 / 0.2e1 + t665 * t303 + (t791 * t662 + t451 * t660 - t581 / 0.2e1 + t436) * t307 - t789;
t372 = t520 / 0.2e1 - t519 / 0.2e1 + t767;
t407 = t60 * t496 + t59 * t497 - t767;
t373 = t659 * t681 + t770 + t803 * t521 + t407 + (t450 * t659 + t673) * t359 + (t273 * t472 + t514 * t718) * mrSges(6,3) + t789;
t1 = t371 + t372 + t373;
t44 = t683 * t306 + t708 * t307 + m(7) * (t302 * t306 + t303 * t475) + m(6) * (-t304 * t306 + t305 * t307) + t461;
t446 = -t1 * qJD(1) + t44 * qJD(2);
t406 = mrSges(7,3) * t476 - t535 - t574 - t576;
t45 = m(7) * (t302 * t305 + t303 * t476) + t406;
t382 = -t756 + t722 + t305 * t174 / 0.2e1;
t400 = (pkin(10) * t450 - t800) * t660 + t769;
t423 = t751 / 0.4e1 + t750 * t633 + t182 * t632;
t424 = t752 / 0.4e1 + t303 * t749 / 0.2e1 + t304 * t178 / 0.2e1;
t6 = (t746 + t693) * Ifges(6,6) + (t647 + t792 / 0.2e1) * mrSges(6,2) + (-t752 / 0.4e1 + t749 * t656 + t424) * t369 + (-t751 / 0.4e1 + pkin(10) * t750 / 0.2e1 + t423) * t367 + (t303 * t448 + t304 * t451 + t452) * t659 + t756 + ((t648 - t60 / 0.2e1) * t369 + (t649 + t59 / 0.2e1) * t367) * mrSges(7,3) + t382 + t400;
t445 = t6 * qJD(1) + t45 * qJD(2);
t444 = t472 + t657;
t398 = (-0.3e1 / 0.2e1 * t600 + t467) * t718 + t465;
t394 = t182 * t633 + t398;
t413 = t468 * t718 + t466;
t402 = t178 * t633 + t413;
t442 = t718 * t464;
t404 = t303 * t442 + t634 * t742 + t486;
t419 = t59 * t652 + t60 * t651 + t757;
t11 = (t764 + t402) * t367 + (t394 + t776) * t369 + t404 + t419;
t438 = t600 + (-Ifges(7,1) + Ifges(7,2)) * t369;
t411 = -Ifges(7,4) * t366 + t367 * t438;
t203 = -t286 + t411;
t443 = t11 * qJD(1) - t203 * qJD(2);
t441 = mrSges(7,1) * t649 + mrSges(7,2) * t648;
t440 = t658 * t742 + t486;
t431 = t522 * t619;
t389 = t517 * t708 + t683 * t615;
t251 = (t359 * t431 + t360 * t368) * t653 + t389;
t375 = (t360 * t305 + (t302 * t368 + t303 * t431) * pkin(4)) * t659 - t576 / 0.2e1 - t574 / 0.2e1 - t535 / 0.2e1 + mrSges(6,1) * t514 - t460 * t513 + mrSges(6,2) * t472 + t471 * t669 + (t474 * t659 + t700) * t304;
t381 = (-pkin(5) * t306 + pkin(10) * t475) * t660 + t572 / 0.2e1 + t534 / 0.2e1 + t570 / 0.2e1 - t700 * t307;
t36 = t375 + t381;
t370 = (t681 + t448 * t359 + (t492 * t58 - t493 * t57 - t710) * pkin(4)) * t659 + Ifges(6,5) * t759 + t770 - t750 * t533 / 0.2e1 + t174 * t513 + t359 * t788 + t72 * t497 + t73 * t496 + pkin(4) * t182 * t469 + t471 * t523 - t807;
t387 = pkin(10) * t665 - t660 * t800 + t769;
t7 = t370 + t372 + t387;
t425 = t7 * qJD(1) + t36 * qJD(2) + t251 * qJD(4);
t116 = t635 + t628 + (mrSges(7,2) * t630 - t599) * t369 + (mrSges(7,1) * t630 + t438) * t367;
t393 = t182 * t627 + t398;
t401 = t178 * t627 + t413;
t403 = t359 * t442 + t626 * t742 + t486;
t420 = t65 * t652 + t651 * t66 + t757;
t13 = (t764 + t401) * t367 + (t393 + t776) * t369 + t403 + t420;
t277 = -t338 + t411;
t422 = t13 * qJD(1) - t116 * qJD(2) - t277 * qJD(4);
t118 = t635 + (mrSges(7,2) * t509 - t599) * t369 + (mrSges(7,1) * t509 + t438) * t367;
t15 = (-Ifges(7,3) / 0.2e1 + pkin(10) * t464) * t718 + (t511 - t594 / 0.2e1 + t718 * t467 + t465) * t369 + (0.3e1 / 0.4e1 * t587 + t512 + (-0.3e1 / 0.2e1 * t599 + t468) * t718 + t466) * t367 + t440 + t441;
t204 = t628 + (mrSges(7,2) * t444 - t599) * t369 + (mrSges(7,1) * t444 + t438) * t367;
t278 = (pkin(5) * mrSges(7,2) - t599) * t369 + (pkin(5) * mrSges(7,1) + t438) * t367;
t412 = t15 * qJD(1) - t118 * qJD(2) - t204 * qJD(4) - t278 * qJD(5);
t408 = t66 * t496 + t65 * t497 - t767;
t337 = t338 / 0.2e1;
t282 = t286 / 0.2e1;
t205 = t337 + (-mrSges(7,2) * t492 / 0.2e1 + mrSges(7,1) * t469) * pkin(4) + t447;
t119 = t304 * t439 + t282 + t447;
t117 = t307 * t439 + t282 + t337 + t666;
t35 = t375 - t381;
t29 = t430 - t436;
t26 = mrSges(6,2) * t759 - t396 + t399 + t498;
t23 = t383 + t768;
t22 = -t377 + t380;
t16 = t718 * t650 + pkin(10) * t442 + (t501 + t512 + t413) * t367 + (t511 + t398) * t369 + t440 - t441 + t704;
t14 = (t501 + t401) * t367 + t393 * t369 + t403 - t420 + t704;
t12 = (t501 + t402) * t367 + t394 * t369 + t404 - t419 + t704;
t8 = t370 + t408 - t387;
t5 = ((-t303 * t72 - t304 * t57) * t659 + mrSges(7,3) * t649 + t423) * t367 + t452 * t659 + (mrSges(7,3) * t648 + (t303 * t73 + t304 * t58) * t659 + t424) * t369 + (t349 / 0.4e1 - Ifges(6,6) / 0.2e1) * t718 - t400 + t407 + t382 + t418 + t673 * pkin(10);
t2 = -t371 + t373 + t408;
t18 = [qJD(2) * t3 + qJD(3) * t17 + qJD(4) * t4 + qJD(5) * t9 - qJD(6) * t10, t23 * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t12 * qJD(6) + t577 + (-t718 * t573 + (mrSges(4,3) * t494 + Ifges(4,6)) * t415 + (-mrSges(4,3) * t489 + Ifges(4,5)) * t346 - Ifges(3,6) * t618 - mrSges(3,1) * t518 + (-t575 + t676) * t273 - t667 * mrSges(4,2) + m(5) * (-t214 * t335 + t336 * t748) + mrSges(3,2) * t516 + t450 * mrSges(7,3) + m(6) * t802 + t793 * t302 + Ifges(3,5) * t621 + (t323 * t336 - t335 * t702) * mrSges(5,3) - t329 * mrSges(4,1) + (-t329 * t550 + t549 * t667) * t655 + (m(7) * t450 + t782) * t303 + t806) * qJD(2), qJD(2) * t23 + qJD(4) * t22 + qJD(5) * t26 + qJD(6) * t29 + t548, t2 * qJD(2) + t22 * qJD(3) + t8 * qJD(5) + t14 * qJD(6) + t562 + (t803 * t654 - t520 - t180 * t533 + t360 * t755 - t615 * t579 + m(7) * (t359 * t449 + t681) + (-mrSges(6,3) * t517 + t676) * t273 + t359 * t525 + t519 + t806) * qJD(4), t5 * qJD(2) + t26 * qJD(3) + t8 * qJD(4) + t16 * qJD(6) + t554 + (-(-Ifges(6,5) - t676) * t273 - t793 * pkin(5) + (m(7) * t448 + t782) * pkin(10) + t448 * mrSges(7,3) + t719) * qJD(5), -t547 + t12 * qJD(2) + t29 * qJD(3) + t14 * qJD(4) + t16 * qJD(5) + (-t58 * mrSges(7,1) - t57 * mrSges(7,2) - t743) * qJD(6); -qJD(3) * t19 - qJD(4) * t1 + qJD(5) * t6 + qJD(6) * t11 - t577, qJD(4) * t44 + qJD(5) * t45 - qJD(6) * t203, -t545 (t461 - t534 - t570 - t572 + (m(7) * t360 - t619 * t654) * t306 + (m(7) * t474 + t368 * t654 + t603 + t604) * t307) * qJD(4) + t35 * qJD(5) + t117 * qJD(6) + t446, t35 * qJD(4) + (m(7) * (-pkin(5) * t305 + pkin(10) * t476) + t406) * qJD(5) + t119 * qJD(6) + t445, t117 * qJD(4) + t119 * qJD(5) + (-t303 * t460 + t675) * qJD(6) + t443; qJD(2) * t19 - qJD(4) * t21 - qJD(5) * t25 - qJD(6) * t28 - t548, t545, 0, -t544, -t543, -qJD(6) * t459 - t536; qJD(2) * t1 + qJD(3) * t21 + qJD(5) * t7 + qJD(6) * t13 - t562, qJD(5) * t36 - qJD(6) * t116 - t446, t544, qJD(5) * t251 - qJD(6) * t277 ((-pkin(5) * t368 + pkin(10) * t431) * t653 + t389) * qJD(5) + t205 * qJD(6) + t425, t205 * qJD(5) + (-t359 * t460 + t675) * qJD(6) + t422; -qJD(2) * t6 + qJD(3) * t25 - qJD(4) * t7 + qJD(6) * t15 - t554, -qJD(4) * t36 - qJD(6) * t118 - t445, t543, -qJD(6) * t204 - t425, -t278 * qJD(6) (-pkin(10) * t460 + t675) * qJD(6) + t412; -qJD(2) * t11 + qJD(3) * t28 - qJD(4) * t13 - qJD(5) * t15 + t547, qJD(4) * t116 + qJD(5) * t118 - t443, t536, qJD(5) * t204 - t422, -t412, 0;];
Cq  = t18;
