% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRP7
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:45
% EndTime: 2019-03-09 06:19:14
% DurationCPUTime: 15.45s
% Computational Cost: add. (32309->746), mult. (64642->961), div. (0->0), fcn. (72710->8), ass. (0->381)
t487 = sin(pkin(10));
t488 = cos(pkin(10));
t490 = sin(qJ(3));
t718 = cos(qJ(3));
t455 = t487 * t490 - t488 * t718;
t457 = t487 * t718 + t490 * t488;
t589 = -pkin(2) * t488 - pkin(1);
t360 = pkin(3) * t455 - pkin(8) * t457 + t589;
t698 = pkin(7) + qJ(2);
t463 = t698 * t487;
t464 = t698 * t488;
t403 = -t490 * t463 + t464 * t718;
t489 = sin(qJ(4));
t491 = cos(qJ(4));
t230 = t491 * t360 - t489 * t403;
t622 = t457 * t491;
t193 = -pkin(9) * t622 + t230;
t167 = t455 * pkin(4) + t193;
t716 = sin(qJ(5));
t585 = t716 * t167;
t628 = t403 * t491;
t194 = t628 + (-pkin(9) * t457 + t360) * t489;
t717 = cos(qJ(5));
t588 = t717 * t194;
t100 = t588 + t585;
t111 = t193 * t716 + t588;
t584 = t716 * t194;
t112 = t193 * t717 - t584;
t583 = t716 * t489;
t458 = -t717 * t491 + t583;
t587 = t717 * t489;
t519 = t491 * t716 + t587;
t623 = t457 * t489;
t599 = mrSges(5,3) * t623;
t367 = -t455 * mrSges(5,2) - t599;
t614 = t491 * t367;
t369 = t455 * mrSges(5,1) - mrSges(5,3) * t622;
t615 = t489 * t369;
t761 = m(7) / 0.2e1;
t763 = m(6) / 0.2e1;
t345 = t458 * t457;
t673 = t345 * mrSges(6,3);
t277 = mrSges(6,1) * t455 + t673;
t674 = t345 * mrSges(7,2);
t278 = -mrSges(7,1) * t455 - t674;
t347 = t519 * t457;
t724 = t519 / 0.2e1;
t668 = t347 * mrSges(6,3);
t274 = -mrSges(6,2) * t455 - t668;
t655 = t455 * mrSges(7,3);
t669 = t347 * mrSges(7,2);
t279 = t655 - t669;
t784 = t279 + t274;
t788 = mrSges(7,2) + mrSges(6,3);
t824 = t458 / 0.2e1;
t829 = -t519 / 0.2e1;
t768 = -t277 * t724 - t278 * t829 - t784 * t824 - t788 * (t345 * t829 + t347 * t824);
t444 = t455 * qJ(6);
t78 = t444 + t100;
t99 = t167 * t717 - t584;
t79 = -t455 * pkin(5) - t99;
t830 = -(-(-t112 + t99) * t519 + (-t100 + t111) * t458) * t763 - (-(-t112 - t79) * t519 + (t111 - t78) * t458) * t761 + t615 / 0.2e1 - t614 / 0.2e1 - t768;
t484 = t491 * pkin(8);
t608 = t491 * pkin(9) + t484;
t753 = -pkin(9) - pkin(8);
t411 = -t753 * t587 + t608 * t716;
t828 = -t411 / 0.2e1;
t773 = t753 * t583 + t608 * t717;
t733 = t773 / 0.2e1;
t707 = pkin(4) * t491;
t478 = -pkin(3) - t707;
t619 = t519 * qJ(6);
t533 = t458 * pkin(5) - t619;
t361 = t478 + t533;
t380 = mrSges(7,1) * t458 - mrSges(7,3) * t519;
t804 = m(7) * t361 + t380;
t384 = -Ifges(6,5) * t458 - Ifges(6,6) * t519;
t385 = -Ifges(7,4) * t458 + Ifges(7,6) * t519;
t546 = t385 + t384;
t806 = t773 * mrSges(7,1);
t807 = t773 * mrSges(6,1);
t820 = t411 * mrSges(7,3);
t821 = t411 * mrSges(6,2);
t827 = t546 - t806 - t807 + t821 - t820;
t825 = t806 / 0.2e1 + t807 / 0.2e1 + t820 / 0.2e1 - t821 / 0.2e1;
t700 = Ifges(7,4) + Ifges(6,5);
t699 = -Ifges(7,6) + Ifges(6,6);
t822 = Ifges(6,3) + Ifges(7,2);
t819 = t784 * t828;
t342 = Ifges(6,4) * t347;
t213 = Ifges(6,2) * t345 - t342;
t667 = t347 * Ifges(7,5);
t178 = -t345 * Ifges(7,1) + t455 * Ifges(7,4) + t667;
t179 = -t345 * Ifges(6,1) + t455 * Ifges(6,5) - t342;
t786 = t179 + t178;
t818 = t213 + t786;
t481 = Ifges(5,5) * t491;
t560 = -Ifges(5,6) * t489 + t481;
t815 = t345 * t733 + t347 * t828;
t813 = -pkin(5) * t773 - qJ(6) * t411;
t563 = t279 / 0.2e1 + t274 / 0.2e1;
t747 = t278 / 0.2e1;
t748 = -t277 / 0.2e1;
t564 = t747 + t748;
t812 = -t563 * t411 + t564 * t773;
t378 = mrSges(7,1) * t519 + t458 * mrSges(7,3);
t379 = mrSges(6,1) * t519 - t458 * mrSges(6,2);
t651 = t458 * Ifges(7,5);
t382 = Ifges(7,3) * t519 - t651;
t450 = Ifges(7,5) * t519;
t383 = t458 * Ifges(7,3) + t450;
t453 = Ifges(6,4) * t458;
t386 = -Ifges(6,2) * t519 - t453;
t649 = t519 * Ifges(6,4);
t387 = -t458 * Ifges(6,2) + t649;
t388 = -Ifges(7,1) * t458 + t450;
t389 = Ifges(7,1) * t519 + t651;
t390 = -Ifges(6,1) * t458 - t649;
t391 = Ifges(6,1) * t519 - t453;
t811 = (-t386 / 0.2e1 + t382 / 0.2e1 - t391 / 0.2e1 - t389 / 0.2e1) * t458 - (t387 / 0.2e1 - t383 / 0.2e1 - t390 / 0.2e1 - t388 / 0.2e1) * t519 + t361 * t378 + t478 * t379;
t346 = t519 * t455;
t741 = -t346 / 0.2e1;
t742 = t346 / 0.2e1;
t348 = t455 * t458;
t810 = t348 / 0.2e1;
t727 = t457 / 0.2e1;
t803 = t391 + t389;
t482 = Ifges(5,4) * t491;
t468 = Ifges(5,1) * t489 + t482;
t690 = Ifges(5,2) * t489;
t783 = t482 - t690;
t802 = t468 + t783;
t705 = pkin(8) * t455;
t710 = pkin(3) * t457;
t372 = t705 + t710;
t402 = t718 * t463 + t464 * t490;
t255 = t491 * t372 + t489 * t402;
t624 = t455 * t491;
t175 = t457 * pkin(4) + pkin(9) * t624 + t255;
t256 = t489 * t372 - t491 * t402;
t625 = t455 * t489;
t206 = pkin(9) * t625 + t256;
t104 = t175 * t717 - t206 * t716;
t105 = t716 * t175 + t717 * t206;
t273 = -mrSges(6,2) * t457 + t346 * mrSges(6,3);
t275 = mrSges(6,1) * t457 - t348 * mrSges(6,3);
t600 = t716 * pkin(4);
t473 = t600 + qJ(6);
t601 = t717 * pkin(4);
t477 = -t601 - pkin(5);
t89 = qJ(6) * t457 + t105;
t90 = -t457 * pkin(5) - t104;
t513 = t104 * mrSges(6,1) / 0.2e1 - t105 * mrSges(6,2) / 0.2e1 + Ifges(7,6) * t741 + Ifges(6,6) * t742 + t89 * mrSges(7,3) / 0.2e1 - t90 * mrSges(7,1) / 0.2e1 + t700 * t810 + t822 * t727;
t549 = t600 / 0.2e1;
t571 = t625 / 0.2e1;
t758 = m(6) * pkin(4);
t604 = t758 / 0.2e1;
t280 = t346 * mrSges(7,2) + mrSges(7,3) * t457;
t746 = t280 / 0.2e1;
t654 = t457 * mrSges(7,1);
t665 = t348 * mrSges(7,2);
t276 = -t654 + t665;
t749 = t276 / 0.2e1;
t801 = (t473 * t89 + t477 * t90) * t761 + Ifges(5,3) * t727 + t255 * mrSges(5,1) / 0.2e1 - t256 * mrSges(5,2) / 0.2e1 + t473 * t746 + t477 * t749 + (t104 * t717 + t105 * t716) * t604 - Ifges(5,5) * t624 / 0.2e1 + Ifges(5,6) * t571 + t275 * t601 / 0.2e1 + t273 * t549 + t513;
t800 = Ifges(4,4) - t560;
t480 = m(7) * qJ(6) + mrSges(7,3);
t798 = qJD(5) * t480;
t797 = t480 * qJD(6);
t339 = Ifges(7,5) * t345;
t176 = t455 * Ifges(7,6) + t347 * Ifges(7,3) - t339;
t214 = -Ifges(7,1) * t347 - t339;
t672 = t345 * Ifges(6,4);
t215 = -Ifges(6,1) * t347 + t672;
t796 = t176 + t215 + t214;
t792 = mrSges(7,2) / 0.2e1;
t177 = -t347 * Ifges(6,2) + t455 * Ifges(6,6) - t672;
t791 = -t177 / 0.4e1;
t210 = mrSges(7,1) * t347 + mrSges(7,3) * t345;
t751 = t210 / 0.2e1;
t790 = -t468 / 0.4e1;
t465 = t489 * mrSges(5,1) + t491 * mrSges(5,2);
t787 = t402 * t465;
t647 = t491 * mrSges(5,1);
t648 = t489 * mrSges(5,2);
t777 = t648 - t647;
t522 = t777 * t457;
t785 = -t277 + t278;
t782 = (Ifges(6,1) + Ifges(7,1)) * t348 + (Ifges(6,4) - Ifges(7,5)) * t346;
t580 = t111 * t724;
t726 = -t458 / 0.2e1;
t779 = t112 * t726 + t580;
t547 = t699 * t345 - t700 * t347;
t778 = -t255 * t489 + t256 * t491;
t765 = t491 ^ 2;
t766 = t489 ^ 2;
t775 = t765 + t766;
t308 = pkin(4) * t623 + t402;
t643 = qJ(6) * t345;
t706 = pkin(5) * t347;
t535 = t643 + t706;
t141 = t308 + t535;
t774 = m(7) * t141 + t210;
t772 = m(7) * t78 + t784;
t771 = m(7) * t79 + t785;
t207 = -t345 * pkin(5) + t347 * qJ(6);
t603 = pkin(4) * t622;
t169 = t207 + t603;
t376 = -pkin(5) * t519 - qJ(6) * t458;
t708 = pkin(4) * t489;
t365 = -t376 + t708;
t212 = -Ifges(7,3) * t345 - t667;
t653 = t458 * mrSges(7,2);
t443 = -t653 / 0.2e1;
t650 = t519 * mrSges(7,2);
t590 = -t650 / 0.2e1;
t640 = t100 * t519;
t728 = t455 / 0.4e1;
t755 = mrSges(6,3) / 0.2e1;
t208 = -mrSges(7,1) * t345 + mrSges(7,3) * t347;
t209 = -mrSges(6,1) * t345 - mrSges(6,2) * t347;
t769 = -t478 * t209 / 0.2e1 - t361 * t208 / 0.2e1 - t308 * t379 / 0.2e1 - t141 * t378 / 0.2e1;
t499 = -t769 + t78 * t590 + t79 * t443 - t640 * t755 + t347 * t382 / 0.4e1 + t345 * t387 / 0.4e1 + t546 * t728 - (t390 + t388 + t383) * t345 / 0.4e1 - (t386 + t803) * t347 / 0.4e1 + t788 * t815 + (t791 + t796 / 0.4e1) * t519 + (t212 / 0.4e1 - t818 / 0.4e1) * t458;
t559 = t411 * t111 + t112 * t773;
t734 = t380 / 0.2e1;
t770 = mrSges(7,2) * t779 + t787 / 0.2e1 + t169 * t734 + t365 * t751 + (t141 * t365 + t169 * t361 - t411 * t78 + t773 * t79 + t559) * t761 + t499;
t764 = m(5) / 0.2e1;
t762 = -m(7) / 0.2e1;
t760 = -pkin(3) / 0.2e1;
t759 = -pkin(8) / 0.2e1;
t757 = m(7) * pkin(4);
t756 = mrSges(5,2) / 0.2e1;
t754 = t99 / 0.2e1;
t752 = t100 / 0.2e1;
t211 = mrSges(6,1) * t347 - mrSges(6,2) * t345;
t750 = t211 / 0.2e1;
t743 = t345 / 0.2e1;
t740 = -t347 / 0.2e1;
t738 = t347 / 0.2e1;
t737 = -t348 / 0.2e1;
t735 = -t380 / 0.2e1;
t730 = -t773 / 0.2e1;
t729 = -t455 / 0.2e1;
t694 = Ifges(5,4) * t489;
t466 = Ifges(5,2) * t491 + t694;
t722 = -t466 / 0.2e1;
t469 = Ifges(5,1) * t491 - t694;
t721 = t469 / 0.4e1;
t720 = t489 / 0.2e1;
t715 = m(6) * t308;
t714 = m(7) * t111;
t713 = m(7) * t376;
t712 = m(7) * t773;
t711 = m(7) * t519;
t709 = pkin(3) * t465;
t702 = t99 * mrSges(6,2);
t701 = t99 * mrSges(7,3);
t696 = t79 + t99;
t692 = Ifges(5,5) * t455;
t688 = Ifges(5,6) * t455;
t683 = t100 * mrSges(6,1);
t682 = t100 * mrSges(7,1);
t679 = t111 * mrSges(6,1);
t678 = t111 * mrSges(7,1);
t677 = t112 * mrSges(6,2);
t676 = t112 * mrSges(7,3);
t309 = -pkin(4) * t625 + t403;
t534 = pkin(5) * t346 + t348 * qJ(6);
t142 = t309 - t534;
t231 = t489 * t360 + t628;
t366 = -t457 * mrSges(5,2) + mrSges(5,3) * t625;
t368 = t457 * mrSges(5,1) + mrSges(5,3) * t624;
t445 = t457 * mrSges(4,1);
t536 = Ifges(7,5) * t348 - Ifges(7,3) * t346;
t537 = Ifges(6,4) * t348 + Ifges(6,2) * t346;
t664 = t348 * mrSges(7,3);
t670 = t346 * mrSges(7,1);
t541 = -t664 - t670;
t666 = t348 * mrSges(6,2);
t671 = t346 * mrSges(6,1);
t542 = t666 - t671;
t3 = t786 * t737 + t782 * t743 + (t589 * mrSges(4,2) + (-Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1) * t348 - (Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t346 + t787 - t800 * t455) * t455 - m(5) * (t230 * t255 + t231 * t256 + t402 * t403) - m(6) * (t100 * t105 + t104 * t99 + t308 * t309) - m(7) * (t141 * t142 + t78 * t89 + t79 * t90) - t589 * t445 + t536 * t740 + t177 * t741 + t176 * t742 + t537 * t738 - t231 * t366 - t256 * t367 - t230 * t368 - t255 * t369 - t309 * t211 - t100 * t273 - t105 * t274 - t99 * t275 - t79 * t276 - t104 * t277 - t90 * t278 - t89 * t279 - t78 * t280 - t142 * t210 + ((Ifges(5,1) * t765 + Ifges(4,1) - Ifges(4,2) - Ifges(5,3) + (t690 - 0.2e1 * t482) * t489 - t822) * t455 + t700 * t345 + t699 * t347 - t465 * t403 + t800 * t457) * t457 - t141 * t541 - t308 * t542;
t675 = t3 * qJD(1);
t652 = t458 * mrSges(6,3);
t510 = -t100 * t673 - t141 * t208 - t308 * t209 + t547 * t729 + t79 * t669;
t6 = ((-Ifges(5,4) * t623 + t692) * t489 + (t688 - pkin(4) * t715 + t231 * mrSges(5,3) + (t482 + (Ifges(5,1) - Ifges(5,2)) * t489) * t457) * t491) * t457 + (-t212 / 0.2e1 + t213 / 0.2e1 + t179 / 0.2e1 + t178 / 0.2e1 - t99 * mrSges(6,3)) * t347 + (t176 / 0.2e1 - t177 / 0.2e1 + t215 / 0.2e1 + t214 / 0.2e1 - t78 * mrSges(7,2)) * t345 - t211 * t603 + t231 * t369 + t402 * t522 + t510 + (-t599 - t367) * t230 - t774 * t169 + (-m(6) * t100 - t772) * t112 + (m(6) * t99 - t771) * t111;
t646 = t6 * qJD(1);
t7 = t177 * t743 + t212 * t738 + t78 * t674 - t510 + t774 * t207 + t771 * t100 + (t668 + t772) * t99 - t796 * t345 / 0.2e1 + t818 * t740;
t645 = t7 * qJD(1);
t644 = t100 - t78;
t381 = mrSges(6,1) * t458 + mrSges(6,2) * t519;
t529 = -t346 * t411 + t348 * t773;
t497 = -m(5) * (-t705 * t775 - t710) / 0.2e1 - m(6) * (t457 * t478 + t529) / 0.2e1 + (t361 * t457 + t529) * t762;
t498 = (t255 * t491 + t489 * t256) * t764 + (-t104 * t458 + t105 * t519) * t763 + (t458 * t90 + t519 * t89) * t761 + t366 * t720 + t491 * t368 / 0.2e1;
t562 = -t280 / 0.2e1 - t273 / 0.2e1;
t565 = -t275 / 0.2e1 + t749;
t594 = t792 + t755;
t15 = t445 + (t735 - t381 / 0.2e1 + t647 / 0.2e1 - t648 / 0.2e1) * t457 + (-mrSges(4,2) + (t765 / 0.2e1 + t766 / 0.2e1) * mrSges(5,3)) * t455 - (-t346 * t594 + t562) * t519 + (t348 * t594 + t565) * t458 + t497 + t498;
t642 = qJD(1) * t15;
t36 = m(7) * (t141 * t345 + t455 * t78) + t345 * t210 + t455 * t279;
t641 = qJD(1) * t36;
t511 = t671 / 0.2e1 + t670 / 0.2e1 - t666 / 0.2e1 + t664 / 0.2e1;
t504 = t534 * t761 + t511;
t12 = -(t345 * t594 + t696 * t761 + t564) * t519 + (t347 * t594 + t644 * t762 + t563) * t458 + t504;
t639 = t12 * qJD(1);
t630 = t402 * t457;
t14 = t784 * t348 - t785 * t346 + (mrSges(4,3) * t455 - t614 + t615) * t455 + (t211 + t210 + (mrSges(4,3) + t465) * t457) * t457 + m(7) * (t141 * t457 - t346 * t79 + t348 * t78) + m(6) * (t100 * t348 + t308 * t457 + t346 * t99) + m(5) * (t630 + (t230 * t489 - t231 * t491) * t455) + m(4) * (-t403 * t455 + t630) + (m(3) * qJ(2) + mrSges(3,3)) * (t487 ^ 2 + t488 ^ 2);
t638 = t14 * qJD(1);
t613 = qJD(6) * t711;
t566 = t455 * t724;
t221 = (t566 + t742) * m(7);
t605 = t221 * qJD(1);
t602 = mrSges(5,3) * t759;
t598 = t473 * t674;
t597 = t713 / 0.2e1;
t596 = t367 * t759;
t595 = t369 * t759;
t593 = -Ifges(5,1) / 0.4e1 + Ifges(5,2) / 0.4e1;
t592 = t669 / 0.2e1;
t582 = t643 / 0.2e1;
t572 = t455 * t726;
t570 = t381 * t727;
t554 = mrSges(6,3) * t601;
t553 = mrSges(6,3) * t600;
t551 = -t601 / 0.2e1;
t550 = -t600 / 0.2e1;
t548 = 0.2e1 * t443;
t545 = t347 * t554;
t544 = t345 * t553;
t23 = -t804 * t376 + t811;
t502 = (-pkin(5) * t90 + qJ(6) * t89) * t761 - pkin(5) * t276 / 0.2e1 + qJ(6) * t746 + t513;
t503 = (-t141 * t376 + t207 * t361 + t411 * t644 + t696 * t773) * t762 + t207 * t735 + t376 * t751;
t4 = (t213 / 0.4e1 + t179 / 0.4e1 + t178 / 0.4e1 - t212 / 0.4e1) * t458 - (t791 + t215 / 0.4e1 + t214 / 0.4e1 + t176 / 0.4e1) * t519 + (-t382 / 0.4e1 + t386 / 0.4e1 + t389 / 0.4e1 + t391 / 0.4e1 + t411 * t755) * t347 + (t383 / 0.4e1 - t387 / 0.4e1 + t388 / 0.4e1 + t390 / 0.4e1 + mrSges(6,3) * t730) * t345 + (-t385 / 0.4e1 - t384 / 0.4e1) * t455 + t502 + (-(t752 - t78 / 0.2e1) * t519 + (t79 / 0.2e1 + t754) * t458 - t815) * mrSges(7,2) + t503 + t769 - t812;
t532 = -t4 * qJD(1) + t23 * qJD(3);
t496 = (-t346 * t477 + t348 * t473) * t761 + (t346 * t717 + t348 * t716) * t604 + mrSges(5,1) * t571 + t624 * t756 + t511;
t10 = t496 + t830;
t530 = t10 * qJD(1);
t527 = -t473 * t458 + t477 * t519;
t187 = t804 * t519;
t506 = (-t141 * t519 + t345 * t361 + t455 * t773) * t761 + t345 * t734 - t519 * t751;
t525 = t90 * t762 + t654 / 0.2e1;
t33 = (t572 + t737) * mrSges(7,2) + t506 + t525;
t526 = -qJD(1) * t33 + qJD(3) * t187;
t521 = -t378 - t379;
t520 = -t458 * t716 - t519 * t717;
t17 = -t709 + (t783 / 0.2e1 + t468 / 0.2e1) * t491 + (pkin(4) * t381 + t722 + t469 / 0.2e1) * t489 + m(6) * t478 * t708 + t811 + t804 * t365;
t514 = -t100 * t411 - t773 * t99 + t559;
t2 = t775 * t457 * t602 + t819 + t779 * mrSges(6,3) - t801 + t770 + t491 * t595 + t489 * t596 - t489 * (t457 * t783 + t688) / 0.4e1 + t773 * t748 - t522 * t760 + ((t308 * t489 + t478 * t622) * pkin(4) + t514) * t763 + t491 * (t457 * t469 + t692) / 0.4e1 + t708 * t750 + t652 * t754 + t560 * t728 + t278 * t733 + t570 * t707 + t622 * t721 + t622 * t722 + (-t802 / 0.4e1 + t790) * t623;
t518 = t2 * qJD(1) + t17 * qJD(3);
t507 = t655 + (t473 * t455 + t78) * t761;
t38 = -t714 / 0.2e1 + t507;
t462 = m(7) * t473 + mrSges(7,3);
t517 = -qJD(1) * t38 - qJD(4) * t462;
t40 = t655 + 0.2e1 * (t588 / 0.4e1 + t585 / 0.4e1 + t444 / 0.2e1 - t100 / 0.4e1) * m(7);
t516 = qJD(1) * t40 + qJD(4) * t480 + t798;
t508 = m(7) * (-pkin(4) * t520 + t527);
t122 = t597 - t508 / 0.2e1;
t509 = (-mrSges(6,1) - mrSges(7,1)) * t600 + (-mrSges(6,2) + mrSges(7,3)) * t601;
t301 = -(t473 * t717 + t477 * t716) * t757 - t509;
t494 = ((t600 - t473) * t411 + (t477 + t601) * t773) * t761 + t473 * t590 + t477 * t443 + t551 * t653 - t550 * t650 - t825;
t501 = pkin(5) * t443 + t619 * t792 + t813 * t762 + t825;
t35 = t494 + t501;
t492 = (t477 * t100 + t473 * t99 + (t716 * t79 + t717 * t78) * pkin(4)) * t762 + t683 / 0.2e1 + t682 / 0.2e1 + t702 / 0.2e1 - t701 / 0.2e1 - t598 / 0.2e1 + t477 * t592 + t277 * t549 + t278 * t550 - t545 / 0.2e1 - t544 / 0.2e1 + t784 * t551;
t505 = (-pkin(5) * t111 + qJ(6) * t112) * t761 - t679 / 0.2e1 - t678 / 0.2e1 - t677 / 0.2e1 + t676 / 0.2e1;
t9 = mrSges(7,2) * t582 + pkin(5) * t592 + t492 + t505;
t512 = -t9 * qJD(1) - t122 * qJD(2) + t35 * qJD(3) - t301 * qJD(4);
t441 = mrSges(7,3) + (0.2e1 * t549 + qJ(6)) * m(7);
t269 = t548 + t712;
t220 = (t566 + t741) * m(7);
t190 = t712 / 0.2e1 + m(7) * t733 + t548;
t85 = t508 / 0.2e1 + t597 + t521;
t39 = (0.2e1 * t444 + t100) * t761 + m(7) * t752 + t279;
t37 = -t669 + t714 / 0.2e1 + t507;
t34 = t494 - t501 + t546;
t32 = mrSges(7,2) * t572 + t665 / 0.2e1 + t506 - t525;
t16 = t380 * t727 + t570 + t522 / 0.2e1 - t562 * t519 + t565 * t458 - t497 + t498 + t775 * mrSges(5,3) * t729 + t788 * (-t346 * t724 + t348 * t726);
t13 = (t458 * t644 + t519 * t696) * t761 + t504 + t768;
t11 = t496 - t830;
t8 = -t492 + (t706 / 0.2e1 + t582) * mrSges(7,2) + t505 + t547;
t5 = t277 * t730 + t773 * t747 + t99 * t443 + t499 + t502 - t503 + t819 + t788 * t640 / 0.2e1;
t1 = (t596 - t688 / 0.2e1 + (t715 / 0.2e1 + t750) * pkin(4) + (t790 - t783 / 0.4e1 + pkin(3) * t756 - t482 + (t602 + t593) * t489) * t457) * t489 + (t595 + t692 / 0.4e1 + (t721 - t466 / 0.4e1 + mrSges(5,1) * t760 + (t602 - t593) * t491 + (t478 * t763 + t381 / 0.2e1) * pkin(4)) * t457) * t491 + (t580 + (-t112 / 0.2e1 + t754) * t458) * mrSges(6,3) + t514 * t763 + t481 * t728 + t770 + t801 + t812;
t18 = [qJD(2) * t14 - qJD(3) * t3 - qJD(4) * t6 + qJD(5) * t7 + qJD(6) * t36, t16 * qJD(3) + t11 * qJD(4) + t13 * qJD(5) + t220 * qJD(6) + t638 + 0.2e1 * (t761 + t763) * qJD(2) * (-t346 * t458 + t348 * t519) t16 * qJD(2) + t1 * qJD(4) + t5 * qJD(5) + t32 * qJD(6) - t675 + (t778 * mrSges(5,3) + 0.2e1 * (-pkin(3) * t403 + pkin(8) * t778) * t764 + t782 * t724 + (t466 * t720 - t489 * t469 / 0.2e1 + t709 - Ifges(4,5) - t802 * t491 / 0.2e1) * t455 + t90 * t650 - t105 * t652 - t89 * t653 + (Ifges(5,5) * t489 + Ifges(5,6) * t491 - t458 * t699 + t519 * t700 - Ifges(4,6)) * t457 + 0.2e1 * (t142 * t361 + t411 * t90 + t773 * t89) * t761 + 0.2e1 * (-t104 * t411 + t105 * t773 + t309 * t478) * t763 + (t273 + t280) * t773 - t104 * t519 * mrSges(6,3) + t383 * t741 + t387 * t742 + t537 * t726 + t366 * t484 - t489 * pkin(8) * t368 + t402 * mrSges(4,2) + t142 * t380 + t309 * t381 + t803 * t810 + (-mrSges(4,1) + t777) * t403 + (t276 - t275) * t411 + t361 * t541 + t478 * t542 + t536 * t824) * qJD(3), -t646 + t11 * qJD(2) + t1 * qJD(3) + (-Ifges(5,5) * t623 - Ifges(5,6) * t622 - t477 * t669 + m(7) * (t111 * t477 + t112 * t473) + t598 + t676 + (-t111 * t717 + t112 * t716) * t758 - t677 - t679 - t678 - t230 * mrSges(5,2) - t231 * mrSges(5,1) + t545 + t544 + t547) * qJD(4) + t8 * qJD(5) + t37 * qJD(6), t645 + t13 * qJD(2) + t5 * qJD(3) + t8 * qJD(4) + (m(7) * (-pkin(5) * t100 + qJ(6) * t99) + t701 - t702 - t683 - t682 + t535 * mrSges(7,2) + t547) * qJD(5) + t39 * qJD(6), qJD(2) * t220 + qJD(3) * t32 + qJD(4) * t37 + qJD(5) * t39 + t641; qJD(3) * t15 - qJD(4) * t10 - qJD(5) * t12 + qJD(6) * t221 - t638, 0, t642 (m(7) * t527 + t520 * t758 - t465 + t521) * qJD(4) + t85 * qJD(5) - t530 + t613, -t639 + t85 * qJD(4) + (t521 + t713) * qJD(5) + t613, t605 - 0.2e1 * (-qJD(4) / 0.2e1 - qJD(5) / 0.2e1) * t711; -qJD(2) * t15 + qJD(4) * t2 - qJD(5) * t4 + qJD(6) * t33 + t675, -t642, qJD(4) * t17 + qJD(5) * t23 - qJD(6) * t187 (-t519 * t553 + t458 * t554 - t477 * t653 - t473 * t650 + (-t411 * t716 - t717 * t773) * t758 + m(7) * (-t411 * t473 + t477 * t773) + t560 + t777 * pkin(8) + t827) * qJD(4) + t34 * qJD(5) + t190 * qJD(6) + t518, t34 * qJD(4) + (m(7) * t813 + t533 * mrSges(7,2) + t827) * qJD(5) + t269 * qJD(6) + t532, qJD(4) * t190 + qJD(5) * t269 - t526; qJD(2) * t10 - qJD(3) * t2 - qJD(5) * t9 + qJD(6) * t38 + t646, -qJD(5) * t122 + t530, qJD(5) * t35 - t518, -qJD(5) * t301 + qJD(6) * t462 ((-pkin(5) * t716 + qJ(6) * t717) * t757 + t509) * qJD(5) + t441 * qJD(6) + t512, qJD(5) * t441 - t517; qJD(2) * t12 + qJD(3) * t4 + qJD(4) * t9 + qJD(6) * t40 - t645, qJD(4) * t122 + t639, -qJD(4) * t35 - t532, -t512 + t797, t797, t516; -qJD(2) * t221 - qJD(3) * t33 - qJD(4) * t38 - qJD(5) * t40 - t641, -t605, t526, t517 - t798, -t516, 0;];
Cq  = t18;
