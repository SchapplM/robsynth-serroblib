% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR14_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_invdynJ_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:57:21
% EndTime: 2019-03-09 14:59:33
% DurationCPUTime: 74.63s
% Computational Cost: add. (88561->1382), mult. (258839->1970), div. (0->0), fcn. (228868->18), ass. (0->622)
t522 = sin(pkin(14));
t525 = sin(pkin(6));
t526 = cos(pkin(14));
t537 = cos(qJ(2));
t528 = cos(pkin(7));
t533 = sin(qJ(2));
t723 = t528 * t533;
t454 = (-t522 * t537 - t526 * t723) * t525;
t443 = qJD(1) * t454;
t455 = (-t522 * t723 + t526 * t537) * t525;
t446 = qJD(1) * t455;
t532 = sin(qJ(4));
t536 = cos(qJ(4));
t523 = sin(pkin(8));
t524 = sin(pkin(7));
t695 = qJD(1) * t525;
t660 = t533 * t695;
t631 = t524 * t660;
t606 = t523 * t631;
t527 = cos(pkin(8));
t724 = t527 * t536;
t320 = -t443 * t724 + t446 * t532 - t536 * t606;
t725 = t527 * t532;
t737 = t523 * t532;
t419 = t528 * t737 + (t522 * t536 + t526 * t725) * t524;
t403 = t419 * qJD(4);
t702 = t320 - t403;
t566 = t443 * t527 + t606;
t321 = t446 * t536 + t532 * t566;
t735 = t523 * t536;
t894 = t524 * (-t522 * t532 + t526 * t724) + t528 * t735;
t402 = t894 * qJD(4);
t701 = t321 - t402;
t529 = cos(pkin(6));
t801 = pkin(1) * t529;
t516 = t537 * t801;
t504 = qJD(1) * t516;
t761 = qJ(3) * t528;
t643 = pkin(10) + t761;
t624 = t643 * t533;
t589 = t525 * t624;
t425 = -qJD(1) * t589 + t504;
t515 = t533 * t801;
t729 = t525 * t537;
t868 = t643 * t729 + t515;
t426 = qJD(1) * t868;
t762 = qJ(3) * t524;
t588 = pkin(2) * t533 - t537 * t762;
t456 = t588 * t695;
t740 = t522 * t528;
t741 = t522 * t524;
t299 = t526 * t425 - t426 * t740 + t456 * t741;
t245 = pkin(11) * t566 + t299;
t727 = t526 * t528;
t734 = t524 * t526;
t298 = -t425 * t522 - t426 * t727 + t456 * t734;
t799 = pkin(11) * t527;
t247 = pkin(3) * t631 - t446 * t799 + t298;
t349 = t426 * t524 + t528 * t456;
t800 = pkin(11) * t523;
t286 = -pkin(3) * t443 - t446 * t800 + t349;
t479 = pkin(2) * t740 + qJ(3) * t734;
t733 = t524 * t527;
t408 = (t523 * t528 + t526 * t733) * pkin(11) + t479;
t509 = pkin(2) * t727;
t642 = qJ(3) + t799;
t422 = pkin(3) * t528 - t642 * t741 + t509;
t742 = t522 * t523;
t449 = (-pkin(3) * t526 - pkin(11) * t742 - pkin(2)) * t524;
t593 = t422 * t527 + t449 * t523;
t290 = -t532 * t408 + t593 * t536;
t693 = qJD(3) * t524;
t707 = -t536 * t245 - t247 * t725 - t286 * t737 + (-t522 * t725 + t526 * t536) * t693 + t290 * qJD(4);
t376 = -t443 * t523 + t527 * t631;
t910 = pkin(12) * t376 - t707;
t167 = -t247 * t523 + t527 * t286;
t657 = qJD(3) * t742;
t909 = -pkin(4) * t702 + t701 * pkin(12) + t524 * t657 - t167;
t803 = sin(qJ(1));
t661 = t803 * t537;
t538 = cos(qJ(1));
t715 = t538 * t533;
t483 = t529 * t715 + t661;
t662 = t803 * t533;
t716 = t537 * t538;
t482 = -t529 * t716 + t662;
t728 = t525 * t538;
t585 = t482 * t528 + t524 * t728;
t357 = -t483 * t526 + t522 * t585;
t430 = -t482 * t524 + t528 * t728;
t867 = t483 * t522 + t526 * t585;
t879 = t430 * t523 + t527 * t867;
t237 = t357 * t536 + t532 * t879;
t312 = t430 * t527 - t867 * t523;
t531 = sin(qJ(5));
t535 = cos(qJ(5));
t908 = t237 * t535 + t312 * t531;
t907 = t237 * t531 - t312 * t535;
t396 = t536 * t408;
t706 = -(t522 * t724 + t526 * t532) * t693 - (t532 * t593 + t396) * qJD(4) + t532 * t245 - t536 * (t247 * t527 + t286 * t523);
t262 = t321 * t531 - t535 * t376;
t474 = -t523 * t734 + t527 * t528;
t354 = t419 * t535 + t474 * t531;
t303 = qJD(5) * t354 + t402 * t531;
t906 = t262 - t303;
t263 = t321 * t535 + t376 * t531;
t353 = t419 * t531 - t535 * t474;
t302 = -qJD(5) * t353 + t402 * t535;
t703 = t263 - t302;
t659 = t537 * t695;
t510 = qJD(1) * t529 + qJD(2);
t743 = t510 * t524;
t453 = t528 * t659 + t743;
t591 = (-t524 ^ 2 - t528 ^ 2) * t660;
t359 = -t453 * t522 + t526 * t591;
t360 = t453 * t526 + t522 * t591;
t280 = t359 * t725 + t360 * t536;
t690 = qJD(4) * t536;
t653 = t523 * t690;
t905 = t280 - t653;
t338 = -t422 * t523 + t527 * t449;
t261 = -pkin(4) * t894 - pkin(12) * t419 + t338;
t291 = t422 * t725 + t449 * t737 + t396;
t269 = pkin(12) * t474 + t291;
t688 = qJD(5) * t535;
t689 = qJD(5) * t531;
t764 = t261 * t688 - t269 * t689 + t909 * t531 - t535 * t910;
t904 = t357 * t532;
t705 = pkin(4) * t376 - t706;
t853 = m(7) / 0.2e1;
t903 = 0.2e1 * t853;
t680 = qJD(1) * qJD(2);
t469 = (qJDD(1) * t533 + t537 * t680) * t525;
t468 = (qJDD(1) * t537 - t533 * t680) * t525;
t678 = qJDD(1) * t529;
t507 = qJDD(2) + t678;
t592 = t468 * t528 + t507 * t524;
t345 = -t469 * t522 + t526 * t592;
t807 = t345 / 0.2e1;
t346 = t469 * t526 + t522 * t592;
t806 = t346 / 0.2e1;
t415 = -t468 * t524 + t507 * t528;
t805 = t415 / 0.2e1;
t902 = -pkin(13) * t702 + t764;
t901 = -pkin(5) * t906 + pkin(13) * t703 + t705;
t480 = -t535 * t527 + t531 * t737;
t899 = -t359 * t523 * t531 + qJD(5) * t480 + t535 * t905;
t279 = -t359 * t724 + t360 * t532;
t691 = qJD(4) * t532;
t654 = t523 * t691;
t898 = -t279 + t654;
t484 = -t529 * t662 + t716;
t572 = t529 * t661 + t715;
t663 = t525 * t803;
t553 = t524 * t663 - t528 * t572;
t546 = t484 * t522 - t526 * t553;
t897 = t572 * t524 + t528 * t663;
t541 = t546 * t523 + t527 * t897;
t722 = t528 * t537;
t726 = t526 * t533;
t386 = t510 * t741 + (t522 * t722 + t726) * t695;
t567 = (-t522 * t533 + t526 * t722) * t525;
t385 = qJD(1) * t567 + t510 * t734;
t452 = t510 * t528 - t524 * t659;
t594 = t385 * t527 + t452 * t523;
t282 = -t386 * t532 + t536 * t594;
t283 = t386 * t536 + t532 * t594;
t534 = cos(qJ(6));
t530 = sin(qJ(6));
t720 = t530 * t535;
t178 = -t282 * t720 + t283 * t534;
t896 = t530 * t688 + t178;
t595 = t345 * t527 + t415 * t523;
t161 = -qJD(4) * t283 - t346 * t532 + t595 * t536;
t619 = -mrSges(7,1) * t534 + mrSges(7,2) * t530;
t574 = m(7) * pkin(5) - t619;
t621 = -mrSges(6,1) * t535 + mrSges(6,2) * t531;
t674 = m(7) * pkin(13) + mrSges(7,3);
t880 = t531 * t674 + t535 * t574 - t621;
t892 = -m(7) - m(6);
t895 = -pkin(4) * t892 + mrSges(5,1) + t880;
t328 = -t385 * t523 + t452 * t527 + qJD(4);
t199 = -t283 * t531 + t328 * t535;
t281 = qJD(5) - t282;
t200 = t283 * t535 + t328 * t531;
t778 = Ifges(6,4) * t200;
t111 = t199 * Ifges(6,2) + t281 * Ifges(6,6) + t778;
t144 = -t200 * t530 + t281 * t534;
t145 = t200 * t534 + t281 * t530;
t193 = qJD(6) - t199;
t73 = t145 * Ifges(7,5) + t144 * Ifges(7,6) + t193 * Ifges(7,3);
t874 = t73 / 0.2e1 - t111 / 0.2e1;
t159 = qJDD(5) - t161;
t828 = t159 / 0.2e1;
t160 = qJD(4) * t282 + t346 * t536 + t532 * t595;
t295 = -t345 * t523 + t415 * t527 + qJDD(4);
t95 = -qJD(5) * t200 - t160 * t531 + t295 * t535;
t834 = t95 / 0.2e1;
t94 = qJD(5) * t199 + t160 * t535 + t295 * t531;
t835 = t94 / 0.2e1;
t847 = Ifges(6,1) * t835 + Ifges(6,4) * t834 + Ifges(6,5) * t828;
t49 = qJD(6) * t144 + t159 * t530 + t534 * t94;
t846 = t49 / 0.2e1;
t50 = -qJD(6) * t145 + t159 * t534 - t530 * t94;
t845 = t50 / 0.2e1;
t93 = qJDD(6) - t95;
t836 = t93 / 0.2e1;
t827 = t160 / 0.2e1;
t826 = t161 / 0.2e1;
t810 = t295 / 0.2e1;
t682 = 0.2e1 * m(6);
t891 = -t682 / 0.2e1;
t383 = qJ(3) * t743 + t426;
t388 = pkin(2) * t510 + t425;
t731 = t524 * t533;
t440 = (-pkin(2) * t537 - qJ(3) * t731 - pkin(1)) * t695;
t267 = t526 * t383 + t388 * t740 + t440 * t741;
t198 = pkin(11) * t594 + t267;
t266 = -t383 * t522 + t388 * t727 + t440 * t734;
t201 = pkin(3) * t452 - t386 * t799 + t266;
t322 = -t388 * t524 + t528 * t440 + qJD(3);
t228 = -pkin(3) * t385 - t386 * t800 + t322;
t601 = t201 * t527 + t228 * t523;
t99 = -t532 * t198 + t536 * t601;
t87 = -pkin(4) * t328 - t99;
t890 = t87 * mrSges(6,2);
t889 = mrSges(5,2) - mrSges(6,3);
t17 = -mrSges(7,1) * t50 + mrSges(7,2) * t49;
t66 = mrSges(6,1) * t159 - mrSges(6,3) * t94;
t888 = t17 - t66;
t887 = mrSges(6,1) + t574;
t865 = mrSges(6,2) - t674;
t782 = mrSges(6,3) * t200;
t147 = mrSges(6,1) * t281 - t782;
t86 = -mrSges(7,1) * t144 + mrSges(7,2) * t145;
t886 = t147 - t86;
t732 = t524 * t529;
t584 = t525 * t722 + t732;
t697 = pkin(10) * t729 + t515;
t413 = qJ(3) * t584 + t697;
t424 = pkin(2) * t529 + t516 - t589;
t730 = t525 * t533;
t668 = t524 * t730;
t698 = pkin(2) * t729 + qJ(3) * t668;
t802 = pkin(1) * t525;
t450 = -t698 - t802;
t294 = t526 * t413 + t424 * t740 + t450 * t741;
t475 = -t524 * t729 + t528 * t529;
t558 = t526 * t732 + t567;
t556 = t558 * t527;
t548 = t475 * t523 + t556;
t222 = pkin(11) * t548 + t294;
t293 = -t413 * t522 + t424 * t727 + t450 * t734;
t417 = t522 * t584 + t525 * t726;
t229 = pkin(3) * t475 - t417 * t799 + t293;
t340 = -t524 * t424 + t528 * t450;
t260 = -pkin(3) * t558 - t417 * t800 + t340;
t119 = t536 * t222 + t229 * t725 + t260 * t737;
t350 = t475 * t527 - t523 * t558;
t107 = pkin(12) * t350 + t119;
t153 = -t229 * t523 + t527 * t260;
t750 = t417 * t532;
t305 = -t475 * t735 - t536 * t556 + t750;
t297 = t305 * pkin(4);
t306 = t417 * t536 + t532 * t548;
t641 = pkin(12) * t306 - t297;
t115 = t153 - t641;
t885 = t535 * t107 + t531 * t115;
t131 = -mrSges(6,1) * t199 + mrSges(6,2) * t200;
t784 = mrSges(5,3) * t283;
t203 = mrSges(5,1) * t328 - t784;
t714 = t131 - t203;
t685 = qJD(6) * t534;
t884 = t531 * t685 + t896;
t718 = t534 * t535;
t179 = t282 * t718 + t283 * t530;
t686 = qJD(6) * t531;
t883 = t530 * t686 - t534 * t688 + t179;
t882 = t531 * t261 + t535 * t269;
t881 = -t523 * t897 + t546 * t527;
t758 = t282 * t531;
t877 = t689 - t758;
t655 = qJD(3) * t722;
t676 = qJD(2) * t801;
t638 = qJD(1) * t676;
t679 = qJDD(1) * t525;
t671 = pkin(10) * t679;
t672 = pkin(1) * t678;
t664 = t533 * t672 + (t638 + t671) * t537;
t694 = qJD(2) * t533;
t304 = t510 * t693 + (-pkin(10) * t694 + t655) * t695 + t592 * qJ(3) + t664;
t467 = t697 * qJD(2);
t502 = t537 * t672;
t692 = qJD(3) * t533;
t656 = t525 * t692;
t629 = t528 * t656;
t317 = -t533 * t671 - t469 * t761 + pkin(2) * t507 + t502 + (-t467 - t629) * qJD(1);
t673 = pkin(1) * t679;
t341 = -t673 - pkin(2) * t468 + (-qJ(3) * t469 - qJD(1) * t656) * t524;
t171 = t526 * t304 + t317 * t740 + t341 * t741;
t138 = pkin(11) * t595 + t171;
t170 = -t304 * t522 + t317 * t727 + t341 * t734;
t140 = pkin(3) * t415 - t346 * t799 + t170;
t240 = -t317 * t524 + t528 * t341 + qJDD(3);
t164 = -pkin(3) * t345 - t346 * t800 + t240;
t651 = t527 * t690;
t34 = t536 * t138 + t140 * t725 + t164 * t737 - t198 * t691 + t201 * t651 + t228 * t653;
t28 = pkin(12) * t295 + t34;
t85 = -t140 * t523 + t527 * t164;
t46 = -pkin(4) * t161 - pkin(12) * t160 + t85;
t100 = t198 * t536 + t532 * t601;
t88 = pkin(12) * t328 + t100;
t139 = -t201 * t523 + t527 * t228;
t96 = -pkin(4) * t282 - pkin(12) * t283 + t139;
t7 = t535 * t28 + t531 * t46 + t96 * t688 - t689 * t88;
t44 = t531 * t96 + t535 * t88;
t8 = -qJD(5) * t44 - t28 * t531 + t46 * t535;
t876 = -t531 * t8 + t535 * t7;
t652 = t527 * t691;
t35 = t536 * (t140 * t527 + t164 * t523) - t532 * t138 - t198 * t690 - t201 * t652 - t228 * t654;
t29 = -pkin(4) * t295 - t35;
t14 = -pkin(5) * t95 - pkin(13) * t94 + t29;
t42 = pkin(13) * t281 + t44;
t59 = -pkin(5) * t199 - pkin(13) * t200 + t87;
t18 = -t42 * t530 + t534 * t59;
t5 = pkin(13) * t159 + t7;
t1 = qJD(6) * t18 + t14 * t530 + t5 * t534;
t19 = t42 * t534 + t530 * t59;
t2 = -qJD(6) * t19 + t14 * t534 - t5 * t530;
t875 = t1 * t534 - t2 * t530;
t625 = -g(1) * t484 - g(2) * t483;
t618 = mrSges(7,1) * t530 + mrSges(7,2) * t534;
t577 = -t618 + t889;
t118 = -t532 * t222 + t536 * (t229 * t527 + t260 * t523);
t505 = t537 * t676;
t362 = t529 * t693 + t505 + (-qJD(2) * t624 + t655) * t525;
t384 = -qJD(2) * t868 - t629;
t416 = (qJD(2) * t588 - t524 * t692) * t525;
t251 = t526 * t362 + t384 * t740 + t416 * t741;
t444 = qJD(2) * t454;
t658 = t525 * t694;
t630 = t524 * t658;
t605 = t523 * t630;
t565 = t444 * t527 + t605;
t207 = pkin(11) * t565 + t251;
t250 = -t362 * t522 + t384 * t727 + t416 * t734;
t445 = qJD(2) * t455;
t210 = pkin(3) * t630 - t445 * t799 + t250;
t319 = -t384 * t524 + t528 * t416;
t254 = -pkin(3) * t444 - t445 * t800 + t319;
t61 = t536 * (t210 * t527 + t254 * t523) - t532 * t207 - t222 * t690 - t229 * t652 - t260 * t654;
t377 = -t444 * t523 + t527 * t630;
t60 = t536 * t207 + t210 * t725 - t222 * t691 + t229 * t651 + t254 * t737 + t260 * t653;
t55 = pkin(12) * t377 + t60;
t142 = -t210 * t523 + t527 * t254;
t211 = qJD(4) * t306 - t444 * t724 + t445 * t532 - t536 * t605;
t212 = t445 * t536 + t565 * t532 + (t536 * t548 - t750) * qJD(4);
t84 = pkin(4) * t211 - pkin(12) * t212 + t142;
t16 = -qJD(5) * t885 - t531 * t55 + t535 * t84;
t573 = m(7) * pkin(12) + t618;
t864 = -m(6) * pkin(12) + mrSges(5,2) - t573;
t863 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t862 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t763 = -qJD(5) * t882 + t531 * t910 + t909 * t535;
t819 = Ifges(4,5) * t806 + Ifges(4,6) * t807 + Ifges(4,3) * t805;
t861 = -mrSges(6,1) * t87 - mrSges(7,1) * t18 + mrSges(7,2) * t19;
t43 = -t531 * t88 + t535 * t96;
t860 = t139 * mrSges(5,1) + t43 * mrSges(6,1) - t44 * mrSges(6,2);
t825 = -t193 / 0.2e1;
t830 = -t145 / 0.2e1;
t832 = -t144 / 0.2e1;
t859 = Ifges(7,5) * t830 + Ifges(7,6) * t832 + Ifges(7,3) * t825;
t858 = mrSges(6,3) * t44 + t861;
t857 = -mrSges(6,3) * t8 + 0.2e1 * t847;
t37 = t94 * Ifges(6,4) + t95 * Ifges(6,2) + t159 * Ifges(6,6);
t9 = Ifges(7,5) * t49 + Ifges(7,6) * t50 + Ifges(7,3) * t93;
t856 = -mrSges(6,3) * t7 - Ifges(6,4) * t835 + Ifges(7,5) * t846 - Ifges(6,2) * t834 - Ifges(6,6) * t828 + Ifges(7,6) * t845 + Ifges(7,3) * t836 - t37 / 0.2e1 + t9 / 0.2e1 + t862;
t36 = Ifges(6,5) * t94 + Ifges(6,6) * t95 + Ifges(6,3) * t159;
t855 = -mrSges(5,3) * t34 + Ifges(6,5) * t835 + Ifges(6,6) * t834 + Ifges(6,3) * t828 + t36 / 0.2e1 + t863 + (-t810 - t295 / 0.2e1) * Ifges(5,6) + (-t826 - t161 / 0.2e1) * Ifges(5,2) + (-t827 - t160 / 0.2e1) * Ifges(5,4);
t854 = t525 ^ 2;
t10 = t49 * Ifges(7,4) + t50 * Ifges(7,2) + t93 * Ifges(7,6);
t851 = t10 / 0.2e1;
t850 = Ifges(7,1) * t846 + Ifges(7,4) * t845 + Ifges(7,5) * t836;
t775 = Ifges(7,4) * t145;
t74 = t144 * Ifges(7,2) + t193 * Ifges(7,6) + t775;
t843 = -t74 / 0.2e1;
t842 = t74 / 0.2e1;
t143 = Ifges(7,4) * t144;
t75 = t145 * Ifges(7,1) + t193 * Ifges(7,5) + t143;
t841 = -t75 / 0.2e1;
t840 = t75 / 0.2e1;
t80 = Ifges(5,5) * t160 + Ifges(5,6) * t161 + Ifges(5,3) * t295;
t839 = t80 / 0.2e1;
t837 = Ifges(5,1) * t827 + Ifges(5,4) * t826 + Ifges(5,5) * t810;
t831 = t144 / 0.2e1;
t829 = t145 / 0.2e1;
t824 = t193 / 0.2e1;
t823 = -t199 / 0.2e1;
t822 = t199 / 0.2e1;
t821 = -t200 / 0.2e1;
t820 = t200 / 0.2e1;
t818 = Ifges(4,4) * t806 + Ifges(4,2) * t807 + Ifges(4,6) * t805;
t817 = Ifges(4,1) * t806 + Ifges(4,4) * t807 + Ifges(4,5) * t805;
t816 = -t281 / 0.2e1;
t815 = t281 / 0.2e1;
t814 = -t282 / 0.2e1;
t813 = t282 / 0.2e1;
t812 = -t283 / 0.2e1;
t811 = t283 / 0.2e1;
t809 = -t328 / 0.2e1;
t808 = t328 / 0.2e1;
t804 = t533 / 0.2e1;
t6 = -pkin(5) * t159 - t8;
t793 = t531 * t6;
t788 = -mrSges(7,3) + mrSges(6,2);
t156 = -pkin(13) * t894 + t882;
t268 = -pkin(4) * t474 - t290;
t181 = pkin(5) * t353 - pkin(13) * t354 + t268;
t101 = -t156 * t530 + t181 * t534;
t787 = qJD(6) * t101 + t530 * t901 + t534 * t902;
t102 = t156 * t534 + t181 * t530;
t786 = -qJD(6) * t102 - t530 * t902 + t534 * t901;
t785 = mrSges(5,3) * t282;
t783 = mrSges(6,3) * t199;
t781 = Ifges(3,4) * t533;
t780 = Ifges(3,4) * t537;
t779 = Ifges(5,4) * t283;
t777 = Ifges(6,4) * t531;
t776 = Ifges(6,4) * t535;
t774 = Ifges(7,4) * t530;
t773 = Ifges(7,4) * t534;
t772 = Ifges(3,2) * t537;
t769 = t385 * Ifges(4,6);
t768 = t386 * Ifges(4,5);
t766 = t452 * Ifges(4,3);
t765 = pkin(5) * t702 - t763;
t185 = pkin(4) * t283 - pkin(12) * t282;
t77 = t531 * t185 + t535 * t99;
t760 = t199 * t530;
t759 = t199 * t534;
t757 = t282 * t535;
t378 = t482 * t522 - t483 * t727;
t754 = t378 * t523;
t380 = -t484 * t727 + t522 * t572;
t753 = t380 * t523;
t747 = t454 * t523;
t738 = t523 * t524;
t736 = t523 * t535;
t721 = t530 * t531;
t719 = t531 * t534;
t503 = Ifges(3,4) * t659;
t717 = t537 * (Ifges(3,1) * t660 + t510 * Ifges(3,5) + t503);
t481 = t527 * t531 + t532 * t736;
t586 = -t481 * t534 + t530 * t735;
t713 = qJD(6) * t586 + t530 * t899 + t534 * t898;
t433 = -t481 * t530 - t534 * t735;
t712 = qJD(6) * t433 + t530 * t898 - t534 * t899;
t309 = -t354 * t530 - t534 * t894;
t182 = qJD(6) * t309 + t302 * t534 + t403 * t530;
t187 = t263 * t534 + t320 * t530;
t711 = t182 - t187;
t310 = t354 * t534 - t530 * t894;
t183 = -qJD(6) * t310 - t302 * t530 + t403 * t534;
t186 = -t263 * t530 + t320 * t534;
t710 = t183 - t186;
t379 = -t482 * t526 - t483 * t740;
t470 = t482 * pkin(2);
t700 = t379 * pkin(3) - t470;
t381 = -t484 * t740 - t526 * t572;
t472 = t572 * pkin(2);
t699 = t381 * pkin(3) - t472;
t696 = t538 * pkin(1) + pkin(10) * t663;
t687 = qJD(6) * t530;
t684 = qJD(6) * t535;
t681 = 0.2e1 * m(7);
t677 = mrSges(3,3) * t729;
t675 = pkin(12) * t688;
t670 = t524 * t735;
t665 = Ifges(3,5) * t469 + Ifges(3,6) * t468 + Ifges(3,3) * t507;
t647 = t854 * t680;
t645 = t688 / 0.2e1;
t644 = -t686 / 0.2e1;
t640 = -t681 / 0.2e1;
t265 = -t345 * mrSges(4,1) + t346 * mrSges(4,2);
t639 = pkin(10) * t658;
t637 = t523 * t668;
t636 = t527 * t668;
t633 = t455 * pkin(3) + pkin(11) * t636 + t698;
t627 = -pkin(1) * t803 + pkin(10) * t728;
t626 = pkin(5) * t531 - pkin(13) * t535;
t623 = mrSges(5,1) * t305 + mrSges(5,2) * t306;
t223 = t306 * t531 - t350 * t535;
t224 = t306 * t535 + t350 * t531;
t622 = mrSges(6,1) * t223 + mrSges(6,2) * t224;
t165 = -t224 * t530 + t305 * t534;
t166 = t224 * t534 + t305 * t530;
t620 = mrSges(7,1) * t165 - mrSges(7,2) * t166;
t617 = Ifges(6,1) * t535 - t777;
t616 = Ifges(7,1) * t534 - t774;
t615 = Ifges(7,1) * t530 + t773;
t614 = -Ifges(6,2) * t531 + t776;
t613 = -Ifges(7,2) * t530 + t773;
t612 = Ifges(7,2) * t534 + t774;
t611 = Ifges(6,5) * t535 - Ifges(6,6) * t531;
t610 = Ifges(7,5) * t534 - Ifges(7,6) * t530;
t609 = Ifges(7,5) * t530 + Ifges(7,6) * t534;
t53 = pkin(13) * t305 + t885;
t106 = -pkin(4) * t350 - t118;
t70 = pkin(5) * t223 - pkin(13) * t224 + t106;
t21 = t53 * t534 + t530 * t70;
t20 = -t53 * t530 + t534 * t70;
t76 = t185 * t535 - t531 * t99;
t57 = -t107 * t531 + t115 * t535;
t168 = t261 * t535 - t269 * t531;
t590 = -mrSges(6,1) + t619;
t587 = pkin(1) * (mrSges(3,1) * t533 + mrSges(3,2) * t537);
t331 = t483 * t733 - t754;
t332 = t484 * t733 - t753;
t41 = -pkin(5) * t281 - t43;
t583 = t41 * t618;
t582 = t87 * (mrSges(6,1) * t531 + mrSges(6,2) * t535);
t581 = t510 * (Ifges(3,5) * t537 - Ifges(3,6) * t533);
t580 = t533 * (Ifges(3,1) * t537 - t781);
t579 = (t772 + t781) * t525;
t15 = -t107 * t689 + t115 * t688 + t531 * t84 + t535 * t55;
t329 = -t454 * t724 + t455 * t532 - t536 * t637;
t330 = t455 * t536 + (t454 * t527 + t637) * t532;
t571 = t330 * pkin(4) + pkin(12) * t329 + t633;
t273 = -t378 * t724 + t379 * t532 - t483 * t670;
t274 = t379 * t536 + (t378 * t527 + t483 * t738) * t532;
t563 = t274 * pkin(4) + pkin(12) * t273 + t483 * t762 + t700;
t275 = -t380 * t724 + t381 * t532 - t484 * t670;
t276 = t381 * t536 + (t380 * t527 + t484 * t738) * t532;
t562 = t276 * pkin(4) + pkin(12) * t275 + t484 * t762 + t699;
t560 = -t483 * pkin(2) + qJ(3) * t430 + t627;
t56 = -pkin(4) * t377 - t61;
t552 = t484 * pkin(2) + qJ(3) * t897 + t696;
t551 = t357 * pkin(3) + pkin(11) * t312 + t560;
t543 = (-g(1) * t332 - g(2) * t331 + g(3) * t747) * pkin(11);
t542 = -g(1) * t897 + g(2) * t430 - g(3) * t475;
t358 = t484 * t526 + t522 * t553;
t540 = t358 * pkin(3) + pkin(11) * t541 + t552;
t238 = t358 * t532 + t536 * t881;
t239 = t358 * t536 - t532 * t881;
t539 = t239 * pkin(4) + t238 * pkin(12) + t540;
t495 = -pkin(5) * t535 - pkin(13) * t531 - pkin(4);
t491 = t626 * qJD(5);
t486 = -pkin(10) * t730 + t516;
t485 = (-mrSges(3,1) * t537 + mrSges(3,2) * t533) * t525;
t478 = -qJ(3) * t741 + t509;
t466 = t505 - t639;
t465 = t697 * qJD(1);
t464 = -pkin(10) * t660 + t504;
t463 = -mrSges(3,2) * t510 + mrSges(3,3) * t659;
t462 = mrSges(3,1) * t510 - mrSges(3,3) * t660;
t458 = pkin(12) * t718 + t495 * t530;
t457 = -pkin(12) * t720 + t495 * t534;
t438 = t510 * Ifges(3,6) + qJD(1) * t579;
t392 = -pkin(10) * t469 - t533 * t638 + t502;
t391 = -qJD(1) * t639 + t664;
t387 = t636 - t747;
t366 = -t495 * t687 + t491 * t534 + (t530 * t689 - t534 * t684) * pkin(12);
t365 = t495 * t685 + t491 * t530 + (-t530 * t684 - t534 * t689) * pkin(12);
t336 = mrSges(4,1) * t452 - mrSges(4,3) * t386;
t335 = -mrSges(4,2) * t452 + mrSges(4,3) * t385;
t318 = -mrSges(4,1) * t385 + mrSges(4,2) * t386;
t301 = mrSges(4,1) * t415 - mrSges(4,3) * t346;
t300 = -mrSges(4,2) * t415 + mrSges(4,3) * t345;
t289 = t386 * Ifges(4,1) + t385 * Ifges(4,4) + t452 * Ifges(4,5);
t288 = Ifges(4,4) * t386 + Ifges(4,2) * t385 + Ifges(4,6) * t452;
t287 = t766 + t768 + t769;
t285 = t330 * t535 + t387 * t531;
t284 = t330 * t531 - t387 * t535;
t278 = Ifges(5,4) * t282;
t236 = -t430 * t735 - t724 * t867 + t904;
t234 = t536 * t879 - t904;
t202 = -mrSges(5,2) * t328 + t785;
t197 = t276 * t535 + t332 * t531;
t196 = t276 * t531 - t332 * t535;
t195 = t274 * t535 + t331 * t531;
t194 = t274 * t531 - t331 * t535;
t192 = Ifges(6,4) * t199;
t184 = -mrSges(5,1) * t282 + mrSges(5,2) * t283;
t177 = t239 * t535 + t531 * t541;
t176 = t239 * t531 - t535 * t541;
t155 = pkin(5) * t894 - t168;
t150 = t283 * Ifges(5,1) + t328 * Ifges(5,5) + t278;
t149 = Ifges(5,2) * t282 + Ifges(5,6) * t328 + t779;
t148 = t283 * Ifges(5,5) + t282 * Ifges(5,6) + t328 * Ifges(5,3);
t146 = -mrSges(6,2) * t281 + t783;
t134 = -qJD(5) * t223 + t212 * t535 + t377 * t531;
t133 = qJD(5) * t224 + t212 * t531 - t377 * t535;
t132 = pkin(5) * t200 - pkin(13) * t199;
t130 = -mrSges(5,2) * t295 + mrSges(5,3) * t161;
t129 = mrSges(5,1) * t295 - mrSges(5,3) * t160;
t126 = t177 * t534 + t238 * t530;
t125 = -t177 * t530 + t238 * t534;
t112 = t200 * Ifges(6,1) + t281 * Ifges(6,5) + t192;
t110 = Ifges(6,5) * t200 + Ifges(6,6) * t199 + Ifges(6,3) * t281;
t109 = mrSges(7,1) * t193 - mrSges(7,3) * t145;
t108 = -mrSges(7,2) * t193 + mrSges(7,3) * t144;
t97 = -mrSges(5,1) * t161 + mrSges(5,2) * t160;
t78 = t282 * t626 + t100;
t72 = qJD(6) * t165 + t134 * t534 + t211 * t530;
t71 = -qJD(6) * t166 - t134 * t530 + t211 * t534;
t67 = -mrSges(6,2) * t159 + mrSges(6,3) * t95;
t65 = pkin(13) * t283 + t77;
t64 = -pkin(5) * t283 - t76;
t52 = -pkin(5) * t305 - t57;
t51 = -mrSges(6,1) * t95 + mrSges(6,2) * t94;
t33 = t132 * t530 + t43 * t534;
t32 = t132 * t534 - t43 * t530;
t26 = t530 * t78 + t534 * t65;
t25 = -t530 * t65 + t534 * t78;
t24 = pkin(5) * t133 - pkin(13) * t134 + t56;
t23 = -mrSges(7,2) * t93 + mrSges(7,3) * t50;
t22 = mrSges(7,1) * t93 - mrSges(7,3) * t49;
t13 = -pkin(5) * t211 - t16;
t12 = pkin(13) * t211 + t15;
t4 = -qJD(6) * t21 - t12 * t530 + t24 * t534;
t3 = qJD(6) * t20 + t12 * t534 + t24 * t530;
t11 = [t306 * t837 + t350 * t839 + t72 * t840 + t71 * t842 + t855 * t305 + t856 * t223 + (-Ifges(5,6) * t808 - Ifges(5,4) * t811 - Ifges(5,2) * t813 + Ifges(6,3) * t815 + Ifges(6,5) * t820 + Ifges(6,6) * t822 + t110 / 0.2e1 - t149 / 0.2e1 - t100 * mrSges(5,3) + t860) * t211 + (t1 * t165 - t166 * t2 - t18 * t72 + t19 * t71) * mrSges(7,3) + (Ifges(5,1) * t306 + Ifges(5,5) * t350) * t827 + (Ifges(5,1) * t212 + Ifges(5,5) * t377) * t811 - t6 * t620 + t29 * t622 + t85 * t623 + (-m(3) * t627 - m(4) * t560 - m(5) * t551 + mrSges(2,1) * t803 + t483 * mrSges(3,1) - t357 * mrSges(4,1) - t237 * mrSges(5,1) + t538 * mrSges(2,2) - t482 * mrSges(3,2) - t867 * mrSges(4,2) - mrSges(3,3) * t728 - t430 * mrSges(4,3) - t312 * mrSges(5,3) + t892 * (t237 * pkin(4) + t236 * pkin(12) + t551) - t887 * t908 + t577 * t236 + t865 * t907) * g(1) + (-m(5) * t540 - t239 * mrSges(5,1) - mrSges(5,3) * t541 - m(3) * t696 - t484 * mrSges(3,1) + mrSges(3,2) * t572 - mrSges(3,3) * t663 - m(4) * t552 - t358 * mrSges(4,1) + mrSges(4,2) * t546 - mrSges(4,3) * t897 - m(7) * (t177 * pkin(5) + t539) - t126 * mrSges(7,1) - t125 * mrSges(7,2) - m(6) * t539 - t177 * mrSges(6,1) - t538 * mrSges(2,1) + mrSges(2,2) * t803 + t889 * t238 + t865 * t176) * g(2) + (t106 * t29 + t15 * t44 + t16 * t43 + t56 * t87 + t57 * t8 + t7 * t885) * t682 / 0.2e1 + t885 * t67 + t857 * t224 + t166 * t850 + t165 * t851 + t469 * (Ifges(3,5) * t529 + (t533 * Ifges(3,1) + t780) * t525) / 0.2e1 + m(4) * (t170 * t293 + t171 * t294 + t240 * t340 + t250 * t266 + t251 * t267 + t319 * t322) + m(5) * (t100 * t60 + t118 * t35 + t119 * t34 + t139 * t142 + t153 * t85 + t61 * t99) + (t1 * t21 + t13 * t41 + t18 * t4 + t19 * t3 + t2 * t20 + t52 * t6) * t681 / 0.2e1 + ((t717 + t581) * t525 / 0.2e1 - t464 * t677) * qJD(2) + (Ifges(4,4) * t417 + Ifges(4,2) * t558 + Ifges(4,6) * t475) * t807 + t417 * t817 + t558 * t818 + t475 * t819 + (Ifges(3,4) * t469 + Ifges(3,2) * t468 + Ifges(3,6) * t507) * t729 / 0.2e1 + t507 * (Ifges(3,3) * t529 + (Ifges(3,5) * t533 + Ifges(3,6) * t537) * t525) / 0.2e1 + (-t438 / 0.2e1 - t465 * mrSges(3,3)) * t658 - (-mrSges(3,1) * t468 + mrSges(3,2) * t469) * t802 + t468 * (Ifges(3,6) * t529 + t579) / 0.2e1 + (-t100 * t377 + t139 * t212 - t34 * t350) * mrSges(5,2) - t587 * t647 - t485 * t673 + (Ifges(7,4) * t166 + Ifges(7,2) * t165) * t845 + (Ifges(7,4) * t72 + Ifges(7,2) * t71) * t831 + (Ifges(5,5) * t212 + Ifges(5,3) * t377) * t808 + (Ifges(5,5) * t306 + Ifges(5,3) * t350) * t810 + (Ifges(4,5) * t417 + Ifges(4,6) * t558 + Ifges(4,3) * t475) * t805 + (Ifges(4,1) * t417 + Ifges(4,4) * t558 + Ifges(4,5) * t475) * t806 + t171 * (-t475 * mrSges(4,2) + mrSges(4,3) * t558) + t240 * (-mrSges(4,1) * t558 + t417 * mrSges(4,2)) + (Ifges(3,1) * t469 + Ifges(3,4) * t468 + Ifges(3,5) * t507) * t730 / 0.2e1 + t392 * (mrSges(3,1) * t529 - mrSges(3,3) * t730) + t391 * (-mrSges(3,2) * t529 + t677) + t212 * t150 / 0.2e1 + Ifges(2,3) * qJDD(1) + t60 * t202 + t61 * t203 + t142 * t184 + (Ifges(5,4) * t306 + Ifges(5,6) * t350) * t826 + (Ifges(5,4) * t212 + Ifges(5,6) * t377) * t813 + t153 * t97 + t15 * t146 + t16 * t147 + t118 * t129 + t119 * t130 + t56 * t131 + t3 * t108 + t4 * t109 + t106 * t51 + t13 * t86 + t41 * (-mrSges(7,1) * t71 + mrSges(7,2) * t72) + t57 * t66 + t52 * t17 + t287 * t630 / 0.2e1 + t267 * (-mrSges(4,2) * t630 + mrSges(4,3) * t444) + t452 * (Ifges(4,5) * t445 + Ifges(4,6) * t444 + Ifges(4,3) * t630) / 0.2e1 + t385 * (Ifges(4,4) * t445 + Ifges(4,2) * t444 + Ifges(4,6) * t630) / 0.2e1 + t386 * (Ifges(4,1) * t445 + Ifges(4,4) * t444 + Ifges(4,5) * t630) / 0.2e1 + t266 * (mrSges(4,1) * t630 - mrSges(4,3) * t445) + (Ifges(6,5) * t815 + Ifges(6,1) * t820 + Ifges(6,4) * t822 - t43 * mrSges(6,3) + t890 + t112 / 0.2e1) * t134 + t486 * (mrSges(3,1) * t507 - mrSges(3,3) * t469) + t20 * t22 + t21 * t23 + t170 * (mrSges(4,1) * t475 - mrSges(4,3) * t417) + t466 * t463 - t467 * t462 + (Ifges(7,1) * t166 + Ifges(7,4) * t165) * t846 + (Ifges(7,1) * t72 + Ifges(7,4) * t71) * t829 + t444 * t288 / 0.2e1 + t322 * (-mrSges(4,1) * t444 + mrSges(4,2) * t445) + t445 * t289 / 0.2e1 + (t537 * (-Ifges(3,2) * t533 + t780) + t580) * t647 / 0.2e1 + t697 * (-mrSges(3,2) * t507 + mrSges(3,3) * t468) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t854 + t391 * t697 + t392 * t486 - t464 * t467 + t465 * t466) + t294 * t300 + t293 * t301 + t319 * t318 + t251 * t335 + t250 * t336 + t340 * t265 + t35 * (mrSges(5,1) * t350 - mrSges(5,3) * t306) + (-Ifges(6,4) * t820 + Ifges(7,5) * t829 - Ifges(6,2) * t822 - Ifges(6,6) * t815 + Ifges(7,6) * t831 + Ifges(7,3) * t824 - t858 + t874) * t133 + t99 * (mrSges(5,1) * t377 - mrSges(5,3) * t212) + t377 * t148 / 0.2e1 + t529 * t665 / 0.2e1 + (Ifges(7,5) * t166 + Ifges(7,6) * t165) * t836 + (Ifges(7,5) * t72 + Ifges(7,6) * t71) * t824; t419 * t837 + t474 * t839 + (Ifges(6,4) * t263 - Ifges(6,2) * t262 + Ifges(6,6) * t320) * t823 + (Ifges(7,5) * t182 + Ifges(7,6) * t183 + Ifges(7,3) * t303) * t824 + (Ifges(7,5) * t187 + Ifges(7,6) * t186 + Ifges(7,3) * t262) * t825 + (Ifges(7,1) * t182 + Ifges(7,4) * t183 + Ifges(7,5) * t303) * t829 + (Ifges(7,1) * t187 + Ifges(7,4) * t186 + Ifges(7,5) * t262) * t830 + (Ifges(7,4) * t182 + Ifges(7,2) * t183 + Ifges(7,6) * t303) * t831 + (Ifges(7,4) * t187 + Ifges(7,2) * t186 + Ifges(7,6) * t262) * t832 + (mrSges(6,1) * t29 + t856) * t353 + (t100 * t702 - t35 * t419 + t701 * t99) * mrSges(5,3) + (-t139 * t702 + t35 * t474 - t376 * t99) * mrSges(5,1) + (t44 * t702 - t703 * t87) * mrSges(6,2) - (mrSges(5,1) * t85 + t855) * t894 + (-g(1) * t562 - g(2) * t563 - g(3) * t571 + t168 * t8 + t268 * t29 + t43 * t763 + t44 * t764 + t7 * t882 + t705 * t87 + t543) * m(6) + t882 * t67 + (mrSges(6,2) * t29 + t857) * t354 + t310 * t850 + t309 * t851 + (Ifges(5,5) * t419 + Ifges(5,3) * t474) * t810 + (Ifges(5,5) * t402 - Ifges(5,6) * t403) * t808 + (Ifges(5,5) * t321 - Ifges(5,6) * t320 + Ifges(5,3) * t376) * t809 + (Ifges(5,1) * t402 - Ifges(5,4) * t403) * t811 + (Ifges(5,1) * t321 - Ifges(5,4) * t320 + Ifges(5,5) * t376) * t812 + (Ifges(5,4) * t402 - Ifges(5,2) * t403) * t813 + (Ifges(5,4) * t321 - Ifges(5,2) * t320 + Ifges(5,6) * t376) * t814 + (Ifges(6,5) * t302 - Ifges(6,6) * t303 + Ifges(6,3) * t403) * t815 + (Ifges(6,5) * t263 - Ifges(6,6) * t262 + Ifges(6,3) * t320) * t816 + (Ifges(6,1) * t302 - Ifges(6,4) * t303 + Ifges(6,5) * t403) * t820 + (Ifges(6,1) * t263 - Ifges(6,4) * t262 + Ifges(6,5) * t320) * t821 + (Ifges(6,4) * t302 - Ifges(6,2) * t303 + Ifges(6,6) * t403) * t822 + (Ifges(7,5) * t310 + Ifges(7,6) * t309) * t836 - t452 * (Ifges(4,5) * t446 + Ifges(4,6) * t443) / 0.2e1 + (-t321 / 0.2e1 + t402 / 0.2e1) * t150 + (t266 * t446 - t267 * t443) * mrSges(4,3) + (-t263 / 0.2e1 + t302 / 0.2e1) * t112 + t268 * t51 + (t100 * t376 - t139 * t701 - t34 * t474 + t419 * t85) * mrSges(5,2) + (-t186 / 0.2e1 + t183 / 0.2e1) * t74 + t665 + (mrSges(3,1) * t572 - t381 * mrSges(4,1) - t276 * mrSges(5,1) + t484 * mrSges(3,2) - t380 * mrSges(4,2) - t332 * mrSges(5,3) + t196 * t788 + t197 * t590 + t275 * t577) * g(1) + (mrSges(3,1) * t482 - t379 * mrSges(4,1) - t274 * mrSges(5,1) + mrSges(3,2) * t483 - t378 * mrSges(4,2) - t331 * mrSges(5,3) + t194 * t788 + t195 * t590 + t273 * t577) * g(2) + (-t455 * mrSges(4,1) - t330 * mrSges(5,1) - t454 * mrSges(4,2) - t387 * mrSges(5,3) + t284 * t788 + t285 * t590 + t329 * t577 + t485) * g(3) + t763 * t147 + t764 * t146 + t765 * t86 + t786 * t109 + (t1 * t102 + t101 * t2 + t155 * t6 - g(2) * (pkin(5) * t195 + pkin(13) * t194 + t563) - g(1) * (pkin(5) * t197 + pkin(13) * t196 + t562) - g(3) * (pkin(5) * t285 + pkin(13) * t284 + t571) + t765 * t41 + t787 * t19 + t786 * t18 + t543) * m(7) + t787 * t108 + (t1 * t309 - t18 * t711 + t19 * t710 - t2 * t310) * mrSges(7,3) + (Ifges(7,4) * t310 + Ifges(7,2) * t309) * t845 + (Ifges(7,1) * t310 + Ifges(7,4) * t309) * t846 + ((Ifges(4,4) * t522 + Ifges(4,2) * t526) * t807 + (Ifges(4,1) * t522 + Ifges(4,4) * t526) * t806 + (Ifges(4,5) * t522 + Ifges(4,6) * t526) * t805 + t240 * (-mrSges(4,1) * t526 + mrSges(4,2) * t522) - pkin(2) * t265 + t522 * t817 + t526 * t818 + (t526 * t335 + (t184 * t523 - t336) * t522) * qJD(3) + (-g(3) * t730 - t170 * t522 + t171 * t526 + t625) * mrSges(4,3)) * t524 + (-t187 / 0.2e1 + t182 / 0.2e1) * t75 + (Ifges(5,4) * t419 + Ifges(5,6) * t474) * t826 + (-t537 * t503 / 0.2e1 - t717 / 0.2e1 - t581 / 0.2e1 + t438 * t804 + (-t580 / 0.2e1 + t772 * t804 + t587) * t695 + (t464 * t537 + t465 * t533) * mrSges(3,3) + (-t766 / 0.2e1 - t266 * mrSges(4,1) - t769 / 0.2e1 - t768 / 0.2e1 + t267 * mrSges(4,2) - t287 / 0.2e1) * t731) * t695 - t167 * t184 + t168 * t66 + t155 * t17 + t102 * t23 + t101 * t22 - t386 * (Ifges(4,1) * t446 + Ifges(4,4) * t443) / 0.2e1 + (Ifges(5,1) * t419 + Ifges(5,5) * t474) * t827 + t478 * t301 + t479 * t300 - t464 * t463 + t465 * t462 - t322 * (-mrSges(4,1) * t443 + mrSges(4,2) * t446) - t446 * t289 / 0.2e1 - t443 * t288 / 0.2e1 - t385 * (Ifges(4,4) * t446 + Ifges(4,2) * t443) / 0.2e1 + (t149 - t110) * (t320 / 0.2e1 - t403 / 0.2e1) + (t111 - t73) * (t262 / 0.2e1 - t303 / 0.2e1) + (-g(3) * t698 + t170 * t478 + t171 * t479 + g(2) * t470 + g(1) * t472 - t266 * t298 - t267 * t299 - t322 * t349 + (-pkin(2) * t240 + (-t266 * t522 + t267 * t526) * qJD(3) + t625 * qJ(3)) * t524) * m(4) + (t170 * mrSges(4,1) - t171 * mrSges(4,2) + 0.2e1 * t819) * t528 + t290 * t129 + t291 * t130 + t6 * (-mrSges(7,1) * t309 + mrSges(7,2) * t310) + t858 * t906 - t299 * t335 - t298 * t336 + t338 * t97 - t349 * t318 + (-t139 * t167 - g(2) * (-pkin(11) * t754 + t700) - g(1) * (-pkin(11) * t753 + t699) + t290 * t35 + t291 * t34 + t338 * t85 - g(3) * (-pkin(11) * t747 + t633) + t706 * t99 + t707 * t100 + (t139 * t657 + t625 * t642) * t524) * m(5) - t376 * t148 / 0.2e1 - t391 * mrSges(3,2) + t392 * mrSges(3,1) + (-mrSges(6,1) * t702 + mrSges(6,3) * t703) * t43 + t705 * t131 + t706 * t203 + t707 * t202 + (-mrSges(7,1) * t710 + mrSges(7,2) * t711) * t41; -t280 * t202 + t433 * t22 - t586 * t23 - t360 * t335 - t359 * t336 + t481 * t67 + t527 * t97 + t888 * t480 - t714 * t279 - t899 * t146 + t713 * t109 + t712 * t108 + (t532 * t130 + t359 * t184 + (t129 - t51) * t536 + (t202 * t536 + t532 * t714) * qJD(4)) * t523 + (-t266 * t359 - t267 * t360 + t240 + t542) * m(4) + (t527 * t85 - t100 * t280 + t279 * t99 + (t139 * t359 + t34 * t532 + t35 * t536 + (t100 * t536 - t532 * t99) * qJD(4)) * t523 + t542) * m(5) + (-t1 * t586 + t18 * t713 + t19 * t712 + t2 * t433 + t480 * t6 + t542) * m(7) + (-t480 * t8 + t481 * t7 - t279 * t87 + (-t29 * t536 + t691 * t87) * t523 - t899 * t44 + t542) * m(6) + t265 + (m(6) * t43 - m(7) * t41 + t886) * (-qJD(5) * t481 + t359 * t736 + t531 * t905); (-mrSges(5,2) * t139 + Ifges(5,1) * t812 + Ifges(5,5) * t809 + t611 * t816 + t614 * t823 + t617 * t821 - t582) * t282 + (t645 - t757 / 0.2e1) * t112 + (Ifges(6,1) * t531 + t776) * t835 + (-Ifges(7,3) * t535 + t531 * t610) * t836 + t179 * t841 + (-Ifges(7,6) * t535 + t531 * t613) * t845 + (-Ifges(7,5) * t535 + t531 * t616) * t846 + t531 * t847 + (-t609 * t686 + (Ifges(7,3) * t531 + t535 * t610) * qJD(5)) * t824 + (Ifges(6,5) * t531 + Ifges(6,6) * t535) * t828 + (-t615 * t686 + (Ifges(7,5) * t531 + t535 * t616) * qJD(5)) * t829 + (-t612 * t686 + (Ifges(7,6) * t531 + t535 * t613) * qJD(5)) * t831 + (Ifges(6,2) * t535 + t777) * t834 + (Ifges(6,5) * t821 - Ifges(5,2) * t814 - Ifges(5,6) * t809 + Ifges(6,6) * t823 + Ifges(6,3) * t816 - t860) * t283 + (t238 * t895 + t239 * t864) * g(1) + (t785 - t202) * t99 + (Ifges(7,1) * t179 + Ifges(7,4) * t178) * t830 + t29 * t621 + (t859 - t874) * t758 + (t234 * t895 - t237 * t864) * g(2) + (-g(1) * t239 + g(2) * t237 - g(3) * t306 - t877 * t44 + (-t688 + t757) * t43 + t876) * mrSges(6,3) + (-t779 + t110) * t812 + t719 * t850 + (t644 * t74 + t645 * t75) * t534 + (Ifges(7,5) * t179 + Ifges(7,6) * t178) * t825 + t530 * t75 * t644 + t1 * (mrSges(7,2) * t535 - mrSges(7,3) * t721) + t2 * (-mrSges(7,1) * t535 - mrSges(7,3) * t719) - t10 * t721 / 0.2e1 + t149 * t811 + (t1 * t458 + t2 * t457) * t903 + (mrSges(7,1) * t877 + mrSges(7,3) * t883 + t25 * t640 + t366 * t903) * t18 + (-mrSges(7,2) * t877 - mrSges(7,3) * t884 + t26 * t640 + t365 * t903) * t19 + ((t41 * t688 + t793) * t903 + m(6) * ((-t43 * t535 - t44 * t531) * qJD(5) + t876) - t146 * t689 + t535 * t67 + t888 * t531) * pkin(12) + (t278 + t150) * t814 + t80 + (Ifges(7,4) * t179 + Ifges(7,2) * t178) * t832 + qJD(5) * t582 + (-t25 + t366) * t109 + (-t26 + t365) * t108 + t618 * t793 + (-t675 - t76) * t147 + (t43 * t76 + t44 * t77) * t891 + t896 * t843 - t77 * t146 + t535 * t37 / 0.2e1 - t535 * t9 / 0.2e1 + (t87 * t891 - t714 + t784) * t100 + t35 * mrSges(5,1) - t34 * mrSges(5,2) + (-m(6) * t29 - t51) * pkin(4) + t457 * t22 + t458 * t23 + (t675 - t64) * t86 + (t199 * t614 + t200 * t617 + t281 * t611) * qJD(5) / 0.2e1 + t874 * t689 + (-m(6) * t641 + m(7) * t297 + t305 * t880 - t306 * t573 + t623) * g(3) + (mrSges(7,1) * t884 - mrSges(7,2) * t883 + t64 * t640) * t41; (t112 + t192) * t823 + t609 * t836 + t685 * t840 + t759 * t841 + t760 * t842 + t687 * t843 + t612 * t845 + t615 * t846 + (-Ifges(6,2) * t823 - Ifges(6,6) * t816 + t859 + t861) * t200 + t6 * t619 + (-t778 + t73) * t821 + t530 * t850 + t534 * t851 + t111 * t820 + (((-t18 * t534 - t19 * t530) * qJD(6) + t875) * t903 - t109 * t685 - t108 * t687 - t530 * t22 + t534 * t23) * pkin(13) + (t783 - t146) * t43 + qJD(6) * t583 + (-t865 * t908 - t887 * t907) * g(2) + t36 + (-0.2e1 * t6 * t853 - t17) * pkin(5) + t863 + (t18 * t32 + t19 * t33) * t640 - t32 * t109 - t33 * t108 + (Ifges(6,1) * t821 + Ifges(6,5) * t816 + t610 * t825 + t613 * t832 + t616 * t830 - t583 - t890) * t199 + (t41 * t640 + t782 + t886) * t44 + (t176 * t887 + t177 * t865) * g(1) + (t144 * t613 + t145 * t616 + t193 * t610) * qJD(6) / 0.2e1 + ((-t687 + t760) * t19 + (-t685 + t759) * t18 + t875) * mrSges(7,3) + (t223 * t574 - t224 * t674 + t622) * g(3); t19 * t109 - t41 * (mrSges(7,1) * t145 + mrSges(7,2) * t144) + (Ifges(7,5) * t144 - Ifges(7,6) * t145) * t825 + (Ifges(7,1) * t144 - t775) * t830 + t74 * t829 - t18 * t108 - g(1) * (mrSges(7,1) * t125 - mrSges(7,2) * t126) - g(2) * ((t234 * t534 + t530 * t908) * mrSges(7,1) + (-t234 * t530 + t534 * t908) * mrSges(7,2)) - g(3) * t620 + (t144 * t18 + t145 * t19) * mrSges(7,3) + t9 + (-Ifges(7,2) * t145 + t143 + t75) * t832 + t862;];
tau  = t11;
