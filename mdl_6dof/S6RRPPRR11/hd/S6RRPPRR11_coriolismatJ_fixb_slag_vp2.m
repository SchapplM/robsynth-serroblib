% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRR11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:39:10
% EndTime: 2019-03-09 09:39:37
% DurationCPUTime: 15.27s
% Computational Cost: add. (34822->921), mult. (79839->1244), div. (0->0), fcn. (87195->10), ass. (0->459)
t545 = sin(qJ(6));
t772 = -t545 / 0.2e1;
t771 = t545 / 0.2e1;
t548 = cos(qJ(6));
t768 = t548 / 0.2e1;
t859 = Ifges(7,3) / 0.2e1;
t541 = sin(pkin(11));
t543 = cos(pkin(11));
t544 = cos(pkin(6));
t542 = sin(pkin(6));
t549 = cos(qJ(2));
t682 = t542 * t549;
t486 = t541 * t682 - t544 * t543;
t546 = sin(qJ(5));
t684 = t541 * t544;
t585 = t543 * t682 + t684;
t766 = cos(qJ(5));
t376 = t486 * t546 - t766 * t585;
t858 = -t376 / 0.2e1;
t547 = sin(qJ(2));
t683 = t542 * t547;
t617 = t766 * t683;
t651 = t546 * t683;
t457 = -t541 * t617 - t543 * t651;
t393 = t457 * t545 + t548 * t682;
t802 = t393 / 0.2e1;
t394 = -t457 * t548 + t545 * t682;
t801 = t394 / 0.2e1;
t570 = -t486 * t766 - t546 * t585;
t854 = -t570 / 0.2e1;
t857 = mrSges(6,1) * t854;
t765 = pkin(1) * t544;
t525 = t549 * t765;
t825 = pkin(3) + pkin(8);
t463 = -t683 * t825 + t525;
t759 = pkin(2) + qJ(4);
t418 = -t544 * t759 - t463;
t628 = -qJ(3) * t547 - pkin(1);
t439 = (-t549 * t759 + t628) * t542;
t281 = t543 * t418 - t439 * t541;
t225 = pkin(4) * t683 + pkin(9) * t486 + t281;
t282 = t541 * t418 + t543 * t439;
t232 = -pkin(9) * t585 + t282;
t107 = t225 * t766 - t546 * t232;
t102 = -pkin(5) * t683 - t107;
t108 = t546 * t225 + t232 * t766;
t320 = -t545 * t570 + t548 * t683;
t321 = t545 * t683 + t548 * t570;
t131 = Ifges(7,5) * t321 + Ifges(7,6) * t320 - Ifges(7,3) * t376;
t755 = mrSges(7,3) * t321;
t212 = -mrSges(7,1) * t376 - t755;
t751 = Ifges(6,4) * t570;
t220 = Ifges(6,2) * t376 + Ifges(6,6) * t683 + t751;
t369 = Ifges(6,4) * t376;
t221 = Ifges(6,1) * t570 + Ifges(6,5) * t683 + t369;
t710 = t545 * mrSges(7,3);
t230 = -mrSges(7,2) * t570 - t376 * t710;
t705 = t548 * mrSges(7,3);
t231 = mrSges(7,1) * t570 - t376 * t705;
t527 = t541 * pkin(4) + qJ(3);
t499 = t541 * t546 - t543 * t766;
t761 = pkin(10) * t499;
t500 = t541 * t766 + t546 * t543;
t762 = pkin(5) * t500;
t397 = t527 + t761 + t762;
t656 = -pkin(9) - t759;
t508 = t656 * t541;
t620 = t543 * t656;
t423 = t508 * t766 + t546 * t620;
t252 = t397 * t548 - t423 * t545;
t253 = t397 * t545 + t423 * t548;
t422 = t508 * t546 - t620 * t766;
t760 = pkin(10) * t500;
t763 = pkin(5) * t499;
t427 = t760 - t763;
t276 = t422 * t545 + t427 * t548;
t277 = -t422 * t548 + t427 * t545;
t741 = Ifges(7,6) * t545;
t746 = Ifges(7,5) * t548;
t607 = -t741 + t746;
t337 = t500 * Ifges(7,3) - t499 * t607;
t536 = Ifges(7,4) * t548;
t608 = Ifges(7,2) * t545 - t536;
t338 = -Ifges(7,6) * t499 + t500 * t608;
t749 = Ifges(7,4) * t545;
t610 = Ifges(7,1) * t548 - t749;
t340 = -Ifges(7,5) * t499 - t500 * t610;
t492 = pkin(8) * t682 + t547 * t765;
t522 = pkin(3) * t682;
t464 = t522 + t492;
t534 = t544 * qJ(3);
t445 = t534 + t464;
t360 = pkin(4) * t585 + t445;
t410 = mrSges(7,2) * t499 + t500 * t710;
t412 = -mrSges(7,1) * t499 + t500 * t705;
t750 = Ifges(6,4) * t499;
t425 = -t500 * Ifges(6,2) - t750;
t498 = Ifges(6,4) * t500;
t426 = -t499 * Ifges(6,1) - t498;
t616 = Ifges(7,3) / 0.4e1 - Ifges(6,1) / 0.4e1 + Ifges(6,2) / 0.4e1;
t494 = t500 * mrSges(6,2);
t627 = -t499 * mrSges(6,1) - t494;
t178 = -mrSges(7,1) * t320 + mrSges(7,2) * t321;
t727 = t570 * mrSges(6,3);
t343 = mrSges(6,1) * t683 - t727;
t633 = t178 / 0.2e1 - t343 / 0.2e1;
t650 = t750 / 0.2e1;
t103 = pkin(10) * t683 + t108;
t156 = -pkin(5) * t376 - pkin(10) * t570 + t360;
t68 = -t103 * t545 + t156 * t548;
t69 = t103 * t548 + t156 * t545;
t725 = t423 * mrSges(6,3);
t787 = t499 / 0.4e1;
t789 = -t499 / 0.4e1;
t853 = pkin(5) * t570;
t234 = -pkin(10) * t376 + t853;
t79 = -t107 * t545 + t234 * t548;
t795 = -t422 / 0.2e1;
t686 = t499 * t548;
t413 = mrSges(7,1) * t500 + mrSges(7,3) * t686;
t796 = t413 / 0.2e1;
t687 = t499 * t545;
t653 = mrSges(7,3) * t687;
t411 = -mrSges(7,2) * t500 + t653;
t797 = t411 / 0.2e1;
t706 = t548 * mrSges(7,2);
t712 = t545 * mrSges(7,1);
t513 = t706 + t712;
t400 = t499 * t513;
t798 = -t400 / 0.2e1;
t399 = t513 * t500;
t799 = -t399 / 0.2e1;
t80 = t107 * t548 + t234 * t545;
t815 = t321 / 0.4e1;
t817 = t320 / 0.4e1;
t819 = t276 / 0.2e1;
t820 = t252 / 0.2e1;
t756 = mrSges(7,3) * t320;
t211 = mrSges(7,2) * t376 + t756;
t821 = t211 / 0.2e1;
t830 = -mrSges(6,2) / 0.2e1;
t833 = m(7) / 0.2e1;
t856 = (-t221 / 0.4e1 - t369 / 0.4e1) * t500 - (t527 * t830 + t498 / 0.4e1 - t426 / 0.4e1 + mrSges(6,3) * t795 - t616 * t499) * t376 + (t527 * mrSges(6,1) / 0.2e1 + t650 - t425 / 0.4e1 + t337 / 0.4e1 - t725 / 0.2e1 + t616 * t500) * t570 + (t102 * t423 + t108 * t422 + t252 * t79 + t253 * t80 + t276 * t68 + t277 * t69) * t833 + t102 * t799 + t108 * t798 + t231 * t820 + t253 * t230 / 0.2e1 + t212 * t819 + t277 * t821 + t338 * t817 + t340 * t815 + t360 * t627 / 0.2e1 + t131 * t789 + t220 * t787 + t68 * t412 / 0.2e1 + t69 * t410 / 0.2e1 + t79 * t796 + t80 * t797 + t633 * t423;
t456 = t541 * t651 - t543 * t617;
t286 = -mrSges(7,2) * t456 + mrSges(7,3) * t393;
t287 = mrSges(7,1) * t456 - mrSges(7,3) * t394;
t846 = t286 * t768 + t287 * t772;
t584 = t211 * t771 + t212 * t768;
t855 = Ifges(6,2) / 0.2e1;
t539 = t545 ^ 2;
t540 = t548 ^ 2;
t852 = (t540 / 0.2e1 + t539 / 0.2e1) * mrSges(7,3);
t591 = t746 / 0.2e1 - t741 / 0.2e1;
t851 = t591 * t500;
t849 = t343 + t727;
t848 = Ifges(7,1) * t545 + t536;
t658 = t541 ^ 2 + t543 ^ 2;
t657 = t539 + t540;
t847 = -t608 + t848;
t661 = t548 * t413;
t673 = t545 * t411;
t845 = -t673 / 0.2e1 - t661 / 0.2e1;
t668 = t548 * t211;
t678 = t545 * t212;
t583 = t678 / 0.2e1 - t668 / 0.2e1;
t605 = -t545 * t79 + t548 * t80;
t451 = t543 * t464;
t521 = pkin(2) * t683;
t460 = t521 + (-qJ(3) * t549 + qJ(4) * t547) * t542;
t654 = pkin(9) * t683;
t267 = pkin(4) * t682 + t451 + (-t460 - t654) * t541;
t323 = t543 * t460 + t541 * t464;
t297 = t543 * t654 + t323;
t145 = t267 * t766 - t546 * t297;
t138 = -pkin(5) * t682 - t145;
t146 = t546 * t267 + t766 * t297;
t244 = -mrSges(7,1) * t393 + mrSges(7,2) * t394;
t748 = Ifges(6,5) * t457;
t778 = -t848 / 0.4e1;
t515 = Ifges(7,2) * t548 + t749;
t780 = -t515 / 0.4e1;
t707 = t548 * mrSges(7,1);
t711 = t545 * mrSges(7,2);
t512 = -t707 + t711;
t783 = -t512 / 0.2e1;
t832 = pkin(5) / 0.2e1;
t844 = t748 / 0.2e1 + t244 * t832 - t145 * mrSges(6,1) / 0.2e1 + t146 * mrSges(6,2) / 0.2e1 + t393 * t780 + t394 * t778 + (pkin(5) * t833 + t783) * t138;
t843 = -m(7) * pkin(5) - mrSges(6,1) + t512;
t647 = -t705 / 0.2e1;
t675 = t545 * t320;
t828 = mrSges(7,3) / 0.2e1;
t842 = t321 * t647 + t675 * t828 - t584;
t841 = Ifges(7,5) * t801 + Ifges(7,6) * t802 + t456 * t859;
t840 = t500 ^ 2;
t839 = m(4) / 0.2e1;
t838 = -m(5) / 0.2e1;
t837 = m(5) / 0.2e1;
t836 = -m(6) / 0.2e1;
t835 = m(6) / 0.2e1;
t834 = -m(7) / 0.2e1;
t831 = mrSges(7,1) / 0.2e1;
t829 = -mrSges(7,2) / 0.2e1;
t139 = pkin(10) * t682 + t146;
t414 = t525 + (-pkin(4) * t543 - t825) * t683;
t235 = pkin(5) * t456 + pkin(10) * t457 + t414;
t87 = -t139 * t545 + t235 * t548;
t827 = t87 / 0.2e1;
t88 = t139 * t548 + t235 * t545;
t826 = -t88 / 0.2e1;
t317 = Ifges(7,4) * t320;
t731 = t376 * Ifges(7,5);
t736 = t321 * Ifges(7,1);
t133 = t317 - t731 + t736;
t824 = -t133 / 0.4e1;
t735 = t321 * Ifges(7,4);
t181 = Ifges(7,1) * t320 - t735;
t822 = -t181 / 0.4e1;
t818 = -t320 / 0.2e1;
t816 = -t321 / 0.2e1;
t814 = t337 / 0.2e1;
t716 = t500 * Ifges(7,6);
t339 = t499 * t608 + t716;
t813 = t339 / 0.4e1;
t717 = t500 * Ifges(7,5);
t341 = -t499 * t610 + t717;
t812 = -t341 / 0.4e1;
t732 = t376 * mrSges(6,3);
t342 = -mrSges(6,2) * t683 + t732;
t811 = -t342 / 0.2e1;
t808 = t376 / 0.2e1;
t806 = -t376 / 0.4e1;
t803 = t570 / 0.2e1;
t398 = t499 * t512;
t800 = t398 / 0.2e1;
t794 = t422 / 0.2e1;
t793 = -t456 / 0.2e1;
t792 = -t457 / 0.2e1;
t791 = -t486 / 0.2e1;
t790 = -t499 / 0.2e1;
t788 = t499 / 0.2e1;
t786 = -t500 / 0.2e1;
t785 = t500 / 0.2e1;
t784 = t500 / 0.4e1;
t782 = t512 / 0.2e1;
t781 = -t513 / 0.2e1;
t779 = t515 / 0.4e1;
t777 = -t541 / 0.2e1;
t776 = t541 / 0.2e1;
t775 = -t543 / 0.2e1;
t774 = t543 / 0.2e1;
t773 = t544 / 0.2e1;
t770 = t545 / 0.4e1;
t769 = -t548 / 0.2e1;
t767 = t548 / 0.4e1;
t758 = mrSges(5,2) * t541;
t757 = mrSges(5,3) * t547;
t753 = Ifges(5,4) * t541;
t752 = Ifges(5,4) * t543;
t745 = Ifges(5,2) * t543;
t743 = Ifges(6,6) * t570;
t740 = Ifges(7,3) * t570;
t730 = t376 * Ifges(7,6);
t737 = t320 * Ifges(7,2);
t132 = -t730 + t735 + t737;
t199 = Ifges(7,4) * t394 + Ifges(7,2) * t393 + Ifges(7,6) * t456;
t200 = Ifges(7,1) * t394 + Ifges(7,4) * t393 + Ifges(7,5) * t456;
t728 = t570 * mrSges(6,2);
t734 = t376 * mrSges(6,1);
t233 = t728 - t734;
t302 = -Ifges(6,1) * t457 - Ifges(6,4) * t456 + Ifges(6,5) * t682;
t322 = -t460 * t541 + t451;
t722 = t457 * mrSges(6,2);
t724 = t456 * mrSges(6,1);
t324 = -t722 + t724;
t720 = t486 * mrSges(5,2);
t395 = mrSges(5,1) * t585 - t720;
t723 = t456 * mrSges(6,3);
t408 = -mrSges(6,2) * t682 - t723;
t721 = t457 * mrSges(6,3);
t409 = mrSges(6,1) * t682 + t721;
t430 = (Ifges(5,6) * t549 + (t745 + t753) * t547) * t542;
t431 = (t549 * Ifges(5,5) + (t541 * Ifges(5,1) + t752) * t547) * t542;
t448 = -mrSges(5,2) * t683 - mrSges(5,3) * t585;
t449 = mrSges(5,1) * t683 + mrSges(5,3) * t486;
t714 = t543 * mrSges(5,1);
t465 = (-t714 + t758) * t683;
t466 = -t534 - t492;
t467 = (-pkin(2) * t549 + t628) * t542;
t491 = pkin(8) * t683 - t525;
t475 = -pkin(2) * t544 + t491;
t484 = (mrSges(5,1) * t549 - t541 * t757) * t542;
t485 = (-mrSges(5,2) * t549 + t543 * t757) * t542;
t488 = (mrSges(4,2) * t549 - mrSges(4,3) * t547) * t542;
t489 = -qJ(3) * t682 + t521;
t503 = -mrSges(4,1) * t682 - t544 * mrSges(4,3);
t518 = Ifges(4,5) * t683;
t519 = Ifges(3,5) * t682;
t621 = Ifges(5,6) * t777 - Ifges(4,4);
t632 = Ifges(6,4) * t457 / 0.2e1 + t456 * t855 - Ifges(6,6) * t682 / 0.2e1 + t841;
t649 = Ifges(6,6) * t793;
t713 = t543 * Ifges(5,6);
t719 = t491 * mrSges(3,2);
t3 = ((-t467 * mrSges(4,3) + t430 * t775 + t475 * mrSges(4,1) + Ifges(3,4) * t682 + Ifges(6,5) * t803 + Ifges(6,6) * t808 + Ifges(5,5) * t791 + (-pkin(1) * mrSges(3,2) + (Ifges(4,6) - t713 / 0.2e1) * t549) * t542 + (Ifges(3,5) / 0.2e1 + t621) * t544) * t549 + (-t467 * mrSges(4,2) - Ifges(3,6) * t544 + (-Ifges(5,1) * t486 - Ifges(5,4) * t684) * t776 + (-Ifges(5,4) * t486 - Ifges(5,2) * t684) * t774 + Ifges(4,5) * t773 - t748 / 0.2e1 + t649 + (-pkin(1) * mrSges(3,1) + (Ifges(5,5) * t541 - Ifges(3,4) - Ifges(4,6) + t713) * t547) * t542 + (t492 + t466) * mrSges(4,1) + (Ifges(4,2) - Ifges(4,3) + Ifges(3,1) - Ifges(3,2) + Ifges(6,3) + Ifges(5,3) + (-t753 / 0.2e1 - t745 / 0.2e1) * t543) * t682) * t547) * t542 + (t518 / 0.2e1 + t519 / 0.2e1 + t430 * t777 + t719 + (mrSges(4,2) - mrSges(3,1)) * t492) * t544 + t491 * t503 + t281 * t484 + t282 * t485 + t489 * t488 + t463 * t395 + t445 * t465 + t323 * t448 + t322 * t449 + t414 * t233 + t108 * t408 + t107 * t409 + t360 * t324 + t146 * t342 + t145 * t343 + t320 * t199 / 0.2e1 + t321 * t200 / 0.2e1 + t69 * t286 + t68 * t287 + t102 * t244 + t88 * t211 + t87 * t212 + t138 * t178 + (t131 / 0.2e1 - t220 / 0.2e1) * t456 + t431 * t791 + t221 * t792 + t133 * t801 + t132 * t802 + t302 * t803 + m(7) * (t102 * t138 + t68 * t87 + t69 * t88) + m(6) * (t107 * t145 + t108 * t146 + t360 * t414) + m(5) * (t281 * t322 + t282 * t323 + t445 * t463) + m(4) * (t466 * t491 + t467 * t489 + t475 * t492) - t632 * t376;
t738 = t3 * qJD(1);
t733 = t376 * mrSges(6,2);
t729 = t570 * mrSges(6,1);
t150 = Ifges(7,6) * t570 - t376 * t608;
t151 = Ifges(7,5) * t570 + t376 * t610;
t229 = t513 * t376;
t644 = -t683 / 0.2e1;
t669 = t548 * t133;
t680 = t545 * t132;
t4 = (Ifges(6,5) * t376 - t743) * t644 + t150 * t818 + t151 * t816 - t102 * t229 - t79 * t212 - t68 * t231 - t80 * t211 - t69 * t230 - m(7) * (t68 * t79 + t69 * t80) - t360 * (t729 + t733) + t220 * t803 + (t376 * t607 + t680 + t740) * t808 + (Ifges(6,1) * t376 + t131 - t751) * t854 + (t732 - t342) * t107 + (-Ifges(6,2) * t570 + t221 + t369 + t669) * t858 + (-m(7) * t102 - t178 + t849) * t108;
t726 = t4 * qJD(1);
t718 = t499 * mrSges(6,3);
t715 = t541 * mrSges(5,1);
t709 = t545 * t68;
t704 = t548 * t69;
t177 = mrSges(7,1) * t321 + mrSges(7,2) * t320;
t179 = Ifges(7,5) * t320 - Ifges(7,6) * t321;
t180 = -Ifges(7,2) * t321 + t317;
t7 = t102 * t177 + t179 * t858 - t69 * t212 + t68 * t211 + (-t132 / 0.2e1 + t181 / 0.2e1 - t69 * mrSges(7,3)) * t321 + (t180 / 0.2e1 + t133 / 0.2e1 - t68 * mrSges(7,3)) * t320;
t702 = t7 * qJD(1);
t701 = t88 * t548;
t606 = -t704 + t709;
t660 = t178 - t343;
t17 = -t585 * t448 + t486 * t449 + t660 * t570 + (t342 + t668 - t678) * t376 + m(7) * (t102 * t570 - t376 * t606) + m(6) * (-t107 * t570 + t108 * t376) + m(5) * (t281 * t486 - t282 * t585);
t700 = qJD(1) * t17;
t406 = -t456 * t545 + t544 * t548;
t407 = t456 * t548 + t544 * t545;
t16 = t407 * t211 + t406 * t212 + t456 * t342 + t660 * t457 + (t233 + t395 - t503) * t544 + (-t448 * t543 + t449 * t541 - t488) * t683 + m(7) * (t102 * t457 + t406 * t68 + t407 * t69) + m(6) * (-t107 * t457 + t108 * t456 + t360 * t544) + m(5) * (t445 * t544 + (t281 * t541 - t282 * t543) * t683) + m(4) * (-t466 * t544 - t467 * t683);
t699 = t16 * qJD(1);
t698 = t253 * t548;
t697 = t323 * t541;
t696 = t570 * t422;
t695 = t406 * t545;
t694 = t407 * t548;
t693 = t422 * t457;
t692 = t422 * t499;
t691 = t456 * t500;
t690 = t457 * t499;
t689 = t499 * t177;
t688 = t499 * t570;
t685 = t500 * t512;
t681 = t543 * t486;
t677 = t545 * t231;
t674 = t545 * t339;
t672 = t545 * t412;
t671 = t545 * t413;
t670 = t545 * t515;
t666 = t548 * t230;
t664 = t548 * t341;
t663 = t548 * t410;
t662 = t548 * t411;
t659 = -Ifges(6,5) * t500 + Ifges(6,6) * t499;
t655 = t838 + t836;
t652 = -Ifges(7,2) / 0.4e1 + Ifges(7,1) / 0.4e1;
t648 = -t710 / 0.2e1;
t646 = 0.2e1 * t534 + t492;
t643 = t683 / 0.2e1;
t642 = t682 / 0.2e1;
t631 = -t536 / 0.4e1 + t778;
t629 = t543 * t759;
t626 = -t500 * mrSges(6,1) + t499 * mrSges(6,2);
t624 = t657 * t499;
t623 = t657 * t376;
t613 = t212 / 0.2e1 + t755 / 0.2e1;
t612 = -t756 / 0.2e1 + t821;
t611 = t547 * t659 / 0.4e1;
t511 = t543 * mrSges(5,2) + t715;
t514 = Ifges(7,5) * t545 + Ifges(7,6) * t548;
t604 = -t87 * t545 + t701;
t30 = t253 * t413 - t422 * t398 + (t339 * t769 + t341 * t772 - mrSges(7,3) * t698 + t514 * t786 + (-t670 / 0.2e1 + t848 * t768) * t499) * t499 + (-t411 + t653) * t252;
t555 = (t252 * t818 + t253 * t816) * mrSges(7,3) + (t514 * t806 + t848 * t815 + t320 * t779 + t548 * t822 + t132 * t767 + (t704 / 0.2e1 - t709 / 0.2e1) * mrSges(7,3) + (t180 + t133) * t770) * t499 + t102 * t800 + t211 * t820 - t253 * t212 / 0.2e1 + t341 * t817 - t321 * t339 / 0.4e1 + t177 * t794 + t179 * t784 + t68 * t797 - t69 * t413 / 0.2e1;
t565 = mrSges(7,1) * t827 + mrSges(7,2) * t826 + t841;
t6 = t555 - t565;
t603 = t6 * qJD(1) - t30 * qJD(2);
t574 = t811 + t583;
t576 = t541 * t585;
t581 = t662 / 0.2e1 - t671 / 0.2e1;
t597 = t252 * t545 - t698;
t553 = (t576 / 0.2e1 - t681 / 0.2e1) * mrSges(5,3) + t581 * t376 + (mrSges(6,3) * t854 - t633) * t499 + (mrSges(6,3) * t858 + t574) * t500 + (-t543 * t281 - t541 * t282 - t486 * t629 + t576 * t759) * t837 + (t107 * t499 - t108 * t500 + t376 * t423 + t696) * t835 + (-t102 * t499 - t376 * t597 + t500 * t606 + t696) * t833 + t570 * t798;
t559 = t463 * t838 + t414 * t836 + (t545 * t88 + t548 * t87) * t834 - t724 / 0.2e1 + t722 / 0.2e1 + t286 * t772 + t287 * t769;
t11 = (-t449 / 0.2e1 + mrSges(5,1) * t643) * t543 + (-t448 / 0.2e1 + mrSges(5,2) * t644) * t541 + t553 + t559;
t46 = (t400 + t718) * t499 + (mrSges(6,3) * t500 - t662 + t671) * t500 + m(7) * (t500 * t597 - t692) + m(6) * (-t423 * t500 - t692) + (m(5) * t759 + mrSges(5,3)) * t658;
t602 = qJD(1) * t11 + qJD(2) * t46;
t554 = (t715 / 0.2e1 + mrSges(4,3) + t511 / 0.2e1 - t626 / 0.2e1) * t544 + t646 * t839 + (t522 + t646) * t837 + (t423 * t456 + t527 * t544 + t360 + t693) * t835 + (t252 * t406 + t253 * t407 + t545 * t69 + t548 * t68 + t693) * t833 - t734 / 0.2e1 + t728 / 0.2e1 + t406 * t796 + t407 * t797 + t457 * t798 - t720 / 0.2e1 + t584;
t556 = -m(4) * t492 / 0.2e1 + (t322 * t543 + t697) * t838 + (-t145 * t499 + t146 * t500) * t836 + (t138 * t499 + t500 * t604) * t834 + t485 * t777;
t15 = (t409 / 0.2e1 - t244 / 0.2e1 - t721 / 0.2e1) * t499 + (mrSges(5,1) * t642 - t484 / 0.2e1) * t543 + t554 + (-t408 / 0.2e1 - t723 / 0.2e1 - t846) * t500 + t556;
t75 = t673 + t661 + mrSges(4,3) + (m(5) + m(4)) * qJ(3) + m(7) * (t252 * t548 + t253 * t545) + m(6) * t527 + t511 - t626;
t601 = qJD(1) * t15 + qJD(2) * t75;
t588 = t406 * t831 + t407 * t829;
t23 = -t689 / 0.2e1 + ((-t675 / 0.2e1 + t321 * t768) * mrSges(7,3) + t584) * t500 + t588;
t561 = (t499 * t852 + t845) * t500 + t398 * t788;
t587 = -t711 / 0.2e1 + t707 / 0.2e1;
t53 = t561 - t587;
t600 = -t23 * qJD(1) + t53 * qJD(2);
t563 = (t545 * t80 + t548 * t79) * t834 - t733 / 0.2e1 + t230 * t772 + t231 * t769;
t569 = (pkin(10) * t623 - t853) * t833 + t570 * t782;
t26 = 0.2e1 * t857 + (t830 + t852) * t376 + t563 + t569;
t564 = -t500 * t852 + (-t657 * t760 + t763) * t833;
t566 = (t276 * t548 + t277 * t545) * t833 + t410 * t771 + t412 * t768;
t57 = t494 + (t783 + mrSges(6,1)) * t499 + t564 - t566;
t599 = qJD(1) * t26 + qJD(2) * t57;
t28 = (mrSges(7,2) * t858 - t612) * t548 + (mrSges(7,1) * t858 + t613) * t545;
t586 = t706 / 0.2e1 + t712 / 0.2e1;
t578 = t586 * t500;
t81 = t578 - t581;
t598 = qJD(1) * t28 + qJD(2) * t81;
t596 = -t276 * t545 + t277 * t548;
t595 = t694 - t695;
t444 = t499 ^ 2;
t562 = (-t657 * t840 - t444) * t833 + (-t444 - t840) * t835 + t658 * t838;
t579 = t657 * t834 + t655;
t111 = t562 + t579;
t557 = (t500 * t623 + t688) * t833 + (t376 * t500 + t688) * t835 + (-t576 + t681) * t837;
t580 = m(7) * (t406 * t548 + t407 * t545);
t50 = t655 * t544 - t580 / 0.2e1 + t557;
t594 = qJD(1) * t50 + qJD(2) * t111;
t593 = t102 + t605;
t592 = t108 + t606;
t590 = pkin(5) * t800 + t422 * t781;
t589 = mrSges(7,1) * t819 + t277 * t829;
t582 = t670 / 0.2e1 + t848 * t769;
t19 = -m(7) * (t252 * t276 + t253 * t277) - t277 * t411 - t253 * t410 - t276 * t413 - t252 * t412 + t422 * t399 - t527 * t627 + (-t674 / 0.2e1 + t664 / 0.2e1 + t426 / 0.2e1 - t498 / 0.2e1 + t851) * t500 + (t338 * t772 + t340 * t768 + t725 + t814 - t425 / 0.2e1 + t650 + (t859 - Ifges(6,1) / 0.2e1 + t855) * t500) * t499 + (-m(7) * t422 + t400 - t718) * t423;
t567 = t151 * t789 + t500 * t824 - (-t717 / 0.2e1 + t812) * t376;
t568 = t150 * t787 + t132 * t784 - (t716 / 0.2e1 + t813) * t376;
t2 = (-t200 / 0.4e1 + mrSges(7,3) * t827 + (m(7) * t827 + t287 / 0.2e1) * pkin(10) + t568) * t545 + (-t199 / 0.4e1 + mrSges(7,3) * t826 + (m(7) * t826 - t286 / 0.2e1) * pkin(10) + t567) * t548 + (t229 / 0.2e1 + t811) * t422 + (-Ifges(6,3) * t549 / 0.2e1 + t611) * t542 + (-t514 / 0.4e1 + Ifges(6,6) / 0.2e1) * t456 + t844 + t856;
t34 = ((t422 + t596) * t833 + t663 / 0.2e1 - t672 / 0.2e1 + t798) * t500 + ((t423 + t597) * t833 + t799 - t581) * t499;
t577 = t2 * qJD(1) - t19 * qJD(2) + t34 * qJD(3);
t134 = m(7) * (0.1e1 - t657) * t500 * t499;
t558 = (-t695 / 0.2e1 + t694 / 0.2e1) * mrSges(7,3) + (-pkin(5) * t457 + pkin(10) * t595) * t833 + mrSges(6,2) * t793 + mrSges(6,1) * t792 + t457 * t782;
t8 = (t677 / 0.2e1 - t666 / 0.2e1 + t593 * t834 + t727 / 0.2e1 - t633) * t500 + (-t229 / 0.2e1 + t592 * t834 - t732 / 0.2e1 - t574) * t499 + t558;
t575 = -t8 * qJD(1) + t34 * qJD(2) + t134 * qJD(3);
t573 = t740 / 0.2e1 + t79 * t831 + t80 * t829;
t572 = t102 * t781 + t177 * t832 + t321 * t779;
t12 = t631 * t320 + (-0.3e1 / 0.4e1 * t730 + t132 / 0.4e1 + t822 + t737 / 0.4e1 + t735 / 0.4e1 + t612 * pkin(10)) * t545 + (0.3e1 / 0.4e1 * t731 - t180 / 0.4e1 + t824 - t736 / 0.4e1 + t613 * pkin(10)) * t548 + t572 + t573;
t292 = -pkin(5) * t513 - t608 * t768 + t610 * t771 - t582;
t299 = (t781 + t586) * t499;
t31 = (0.3e1 / 0.4e1 * t716 + pkin(10) * t797 + t813) * t545 + (-Ifges(7,3) / 0.2e1 - pkin(10) * t852 + (-t545 * t652 + t631) * t545) * t499 + (-0.3e1 / 0.4e1 * t717 + pkin(10) * t796 + t812 + (-0.3e1 / 0.4e1 * t749 + t780 + t652 * t548) * t499) * t548 + t589 + t590;
t571 = t12 * qJD(1) + t31 * qJD(2) + t299 * qJD(3) - t292 * qJD(5);
t300 = t499 * t586 + t513 * t788;
t110 = t562 - t579;
t82 = t578 + t581;
t65 = t499 * t783 + t564 + t566;
t54 = t561 + t587;
t49 = t580 / 0.2e1 + t557 + (m(6) + m(5)) * t773;
t33 = t34 * qJD(5);
t32 = -t674 / 0.4e1 + t664 / 0.4e1 + t607 * t784 + Ifges(7,3) * t790 - t851 + t589 - t590 + (t515 / 0.2e1 - t610 / 0.4e1) * t686 + t657 * t761 * t828 + t845 * pkin(10) + (t848 + t847) * t687 / 0.4e1;
t29 = t320 * t647 + t321 * t648 - t376 * t586 - t583;
t27 = t729 / 0.2e1 + mrSges(6,2) * t858 + t857 + t376 * t852 - t563 + t569;
t24 = t689 / 0.2e1 + t588 + t842 * t500;
t14 = (-t691 / 0.2e1 - t690 / 0.2e1) * mrSges(6,3) + t484 * t774 + t408 * t785 + t409 * t790 + t244 * t788 + t554 + (t714 / 0.2e1 + mrSges(4,1)) * t682 - t556 + t846 * t500;
t13 = -t680 / 0.4e1 + t180 * t767 + t669 / 0.4e1 + t181 * t770 + t607 * t806 + t610 * t815 + t591 * t376 - t572 + t573 + t847 * t817 + t842 * pkin(10);
t10 = t448 * t777 + t449 * t775 + t643 * t758 + t644 * t714 + t553 - t559;
t9 = t593 * t500 * t833 + t229 * t788 + t342 * t790 + t718 * t808 + t558 + (t178 + t666) * t785 + (t592 * t833 + t583) * t499 + (t677 + t849) * t786;
t5 = t555 + t565;
t1 = t542 * t611 + t199 * t767 + t200 * t770 + t456 * t514 / 0.4e1 + t649 + t229 * t794 + t342 * t795 + t568 * t545 + Ifges(6,3) * t642 + t701 * t828 + t87 * t648 + t567 * t548 + (t604 * t833 + t846) * pkin(10) - t844 + t856;
t18 = [qJD(2) * t3 + qJD(3) * t16 + qJD(4) * t17 - qJD(5) * t4 + qJD(6) * t7, t14 * qJD(3) + t10 * qJD(4) + t1 * qJD(5) + t5 * qJD(6) + t738 + (t527 * t324 + t463 * t511 - t491 * mrSges(4,3) + t492 * mrSges(4,2) - t492 * mrSges(3,1) + qJ(3) * t465 - t422 * t409 + t423 * t408 + t422 * t244 + t87 * t413 + t88 * t411 - t138 * t400 + t253 * t286 + t252 * t287 + 0.2e1 * (-t145 * t422 + t146 * t423 + t414 * t527) * t835 + 0.2e1 * (-pkin(2) * t492 - qJ(3) * t491) * t839 + 0.2e1 * (t138 * t422 + t252 * t87 + t253 * t88) * t833 + t719 + t456 * t814 + t426 * t792 + t425 * t793 + t341 * t801 + t339 * t802 + 0.2e1 * (qJ(3) * t463 - t322 * t629 - t697 * t759) * t837 + t518 + t519 - t414 * t626 + (t431 / 0.2e1 - t322 * mrSges(5,3) - t759 * t484) * t543 + (-t430 / 0.2e1 - t323 * mrSges(5,3) - t759 * t485) * t541 + (-t146 * mrSges(6,3) + t632) * t500 + (-t302 / 0.2e1 + t200 * t769 + t199 * t771 + t145 * mrSges(6,3)) * t499 + ((-pkin(2) * mrSges(4,1) + Ifges(5,5) * t774 + Ifges(6,5) * t790 + Ifges(6,6) * t786 + t621) * t549 + ((-Ifges(5,2) * t541 + t752) * t774 - qJ(3) * mrSges(4,1) + (Ifges(5,1) * t543 - t753) * t776 - Ifges(3,6)) * t547) * t542) * qJD(2), t699 + t14 * qJD(2) + 0.2e1 * ((t500 * t595 + t690) * t833 + (t690 + t691) * t835) * qJD(3) + t49 * qJD(4) + t9 * qJD(5) + t24 * qJD(6), qJD(2) * t10 + qJD(3) * t49 + qJD(5) * t27 + qJD(6) * t29 + t700, t1 * qJD(2) + t9 * qJD(3) + t27 * qJD(4) + t13 * qJD(6) - t726 + (t150 * t768 + t151 * t771 - pkin(5) * t229 + t514 * t803 - t107 * mrSges(6,2) - t743 - (-Ifges(6,5) + t582) * t376 + t843 * t108 + (m(7) * t605 + t666 - t677) * pkin(10) + t605 * mrSges(7,3)) * qJD(5), t702 + t5 * qJD(2) + t24 * qJD(3) + t29 * qJD(4) + t13 * qJD(5) + (-mrSges(7,1) * t69 - mrSges(7,2) * t68 + t179) * qJD(6); qJD(3) * t15 + qJD(4) * t11 + qJD(5) * t2 + qJD(6) * t6 - t738, qJD(3) * t75 + qJD(4) * t46 - qJD(5) * t19 - qJD(6) * t30, qJD(4) * t110 + qJD(6) * t54 + t33 + t601, qJD(3) * t110 + qJD(5) * t65 + qJD(6) * t82 + t602, t65 * qJD(4) + (t422 * mrSges(6,2) + t596 * mrSges(7,3) + pkin(5) * t399 + t338 * t768 + t340 * t771 + t582 * t500 + t514 * t790 + t659 + t843 * t423 + (m(7) * t596 + t663 - t672) * pkin(10)) * qJD(5) + t32 * qJD(6) + t577, t54 * qJD(3) + t82 * qJD(4) + t32 * qJD(5) + (-mrSges(7,1) * t253 - mrSges(7,2) * t252 + t499 * t514) * qJD(6) + t603; -qJD(2) * t15 + qJD(4) * t50 - qJD(5) * t8 - qJD(6) * t23 - t699, qJD(4) * t111 + qJD(6) * t53 + t33 - t601, t134 * qJD(5), t594 (t685 + m(7) * (-pkin(10) * t624 - t762) - mrSges(7,3) * t624 + t626) * qJD(5) + t300 * qJD(6) + t575, t300 * qJD(5) + qJD(6) * t685 + t600; -qJD(2) * t11 - qJD(3) * t50 - qJD(5) * t26 - qJD(6) * t28 - t700, -qJD(3) * t111 - qJD(5) * t57 - qJD(6) * t81 - t602, -t594, 0, -t599, -qJD(6) * t513 - t598; -qJD(2) * t2 + qJD(3) * t8 + qJD(4) * t26 - qJD(6) * t12 + t726, qJD(4) * t57 - qJD(6) * t31 - t577, -qJD(6) * t299 - t575, t599, t292 * qJD(6) (pkin(10) * t512 + t607) * qJD(6) - t571; -qJD(2) * t6 + qJD(3) * t23 + qJD(4) * t28 + qJD(5) * t12 - t702, -qJD(3) * t53 + qJD(4) * t81 + qJD(5) * t31 - t603, qJD(5) * t299 - t600, t598, t571, 0;];
Cq  = t18;
