% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:06:31
% EndTime: 2019-03-09 18:07:05
% DurationCPUTime: 20.85s
% Computational Cost: add. (65449->768), mult. (125533->1000), div. (0->0), fcn. (151137->10), ass. (0->457)
t529 = sin(qJ(3));
t530 = sin(qJ(2));
t533 = cos(qJ(3));
t534 = cos(qJ(2));
t495 = -t529 * t530 + t533 * t534;
t496 = -t529 * t534 - t533 * t530;
t526 = sin(pkin(11));
t752 = cos(pkin(11));
t425 = t752 * t495 + t496 * t526;
t527 = sin(qJ(6));
t528 = sin(qJ(5));
t531 = cos(qJ(6));
t532 = cos(qJ(5));
t494 = -t527 * t532 - t531 * t528;
t300 = t494 * t425;
t493 = -t527 * t528 + t531 * t532;
t302 = t493 * t425;
t580 = t526 * t495 - t496 * t752;
t134 = Ifges(7,4) * t302 + Ifges(7,2) * t300 + Ifges(7,6) * t580;
t136 = Ifges(7,1) * t302 + Ifges(7,4) * t300 + Ifges(7,5) * t580;
t521 = Ifges(6,4) * t532;
t894 = -Ifges(6,2) * t528 + t521;
t226 = Ifges(6,6) * t580 + t425 * t894;
t799 = Ifges(6,4) * t528;
t506 = Ifges(6,1) * t532 - t799;
t228 = Ifges(6,5) * t580 + t425 * t506;
t797 = Ifges(7,4) * t494;
t432 = Ifges(7,2) * t493 - t797;
t482 = Ifges(7,4) * t493;
t434 = -Ifges(7,1) * t494 + t482;
t833 = t493 / 0.2e1;
t846 = t580 / 0.2e1;
t503 = Ifges(6,2) * t532 + t799;
t505 = Ifges(6,1) * t528 + t521;
t930 = t532 / 0.2e1;
t972 = -t528 / 0.2e1;
t920 = t503 * t972 + t505 * t930;
t931 = t528 / 0.2e1;
t932 = -t494 / 0.2e1;
t558 = t134 * t833 + t136 * t932 + t228 * t931 + t226 * t930 + t300 * t432 / 0.2e1 + t302 * t434 / 0.2e1 - Ifges(5,6) * t580 + Ifges(4,6) * t496 + Ifges(4,5) * t495 + (Ifges(6,5) * t528 - Ifges(7,5) * t494 + Ifges(6,6) * t532 + Ifges(7,6) * t493) * t846 + (Ifges(5,5) + t920) * t425;
t863 = -pkin(8) - pkin(7);
t508 = t863 * t530;
t509 = t863 * t534;
t896 = t533 * t508 + t529 * t509;
t921 = t896 * mrSges(4,2);
t445 = t508 * t529 - t509 * t533;
t922 = t445 * mrSges(4,1);
t807 = mrSges(6,1) * t532;
t885 = mrSges(6,2) * t528 - t807;
t395 = qJ(4) * t495 + t445;
t917 = t496 * qJ(4) + t896;
t937 = t752 * t395 + t526 * t917;
t953 = t937 * t885;
t959 = t937 * mrSges(5,1);
t936 = -t526 * t395 + t752 * t917;
t961 = t936 * mrSges(5,2);
t430 = -mrSges(7,1) * t493 - mrSges(7,2) * t494;
t728 = t425 * t528;
t948 = pkin(5) * t728 + t937;
t967 = t948 * t430;
t974 = t558 - t921 - t922 + t953 - t959 - t961 + t967;
t973 = -t921 / 0.2e1 - t922 / 0.2e1 + t953 / 0.2e1 - t959 / 0.2e1 - t961 / 0.2e1 + t967 / 0.2e1;
t519 = -pkin(2) * t534 - pkin(1);
t451 = -t495 * pkin(3) + t519;
t971 = m(5) * t451;
t726 = t580 * t528;
t206 = pkin(5) * t726 - t936;
t970 = t206 * t948;
t824 = pkin(2) * t529;
t510 = t526 * t824;
t823 = pkin(2) * t533;
t518 = pkin(3) + t823;
t464 = t518 * t752 - t510;
t459 = -pkin(4) - t464;
t818 = t532 * pkin(5);
t448 = t459 - t818;
t969 = t448 * t948;
t644 = t752 * pkin(3);
t513 = -t644 - pkin(4);
t500 = t513 - t818;
t968 = t500 * t948;
t951 = -Ifges(5,1) * t580 - Ifges(5,4) * t425;
t301 = t494 * t580;
t299 = t493 * t580;
t798 = Ifges(7,4) * t299;
t135 = Ifges(7,2) * t301 - Ifges(7,6) * t425 + t798;
t296 = Ifges(7,4) * t301;
t137 = Ifges(7,1) * t299 - Ifges(7,5) * t425 + t296;
t269 = -pkin(4) * t425 - pkin(9) * t580 + t451;
t149 = t532 * t269 - t528 * t937;
t150 = t269 * t528 + t532 * t937;
t780 = t302 * mrSges(7,2);
t782 = t300 * mrSges(7,1);
t174 = t780 - t782;
t175 = -mrSges(7,1) * t301 + mrSges(7,2) * t299;
t216 = -mrSges(7,2) * t580 + mrSges(7,3) * t300;
t218 = mrSges(7,1) * t580 - mrSges(7,3) * t302;
t754 = t532 * mrSges(6,2);
t755 = t528 * mrSges(6,1);
t502 = t754 + t755;
t308 = t502 * t425;
t309 = t502 * t580;
t321 = -mrSges(6,2) * t580 - mrSges(6,3) * t728;
t727 = t425 * t532;
t323 = mrSges(6,1) * t580 - mrSges(6,3) * t727;
t417 = t425 * mrSges(5,2);
t725 = t580 * t532;
t121 = -pkin(10) * t725 + t149;
t110 = -pkin(5) * t425 + t121;
t122 = -pkin(10) * t726 + t150;
t750 = t122 * t527;
t67 = t110 * t531 - t750;
t749 = t122 * t531;
t68 = t110 * t527 + t749;
t756 = t496 * mrSges(4,3);
t769 = t580 * mrSges(5,1);
t905 = t580 * mrSges(5,3);
t966 = -t948 * t175 - t149 * t323 - t150 * t321 + t936 * t308 - t451 * (t417 + t769) - t445 * t756 - t519 * (-mrSges(4,1) * t496 + mrSges(4,2) * t495) - t67 * t218 - t68 * t216 - t206 * t174 - t299 * t136 / 0.2e1 - t300 * t135 / 0.2e1 - t301 * t134 / 0.2e1 - t302 * t137 / 0.2e1 + (-t309 + t905) * t937;
t963 = mrSges(7,1) / 0.2e1;
t962 = -mrSges(7,2) / 0.2e1;
t960 = t936 * mrSges(5,3);
t958 = t459 * t937;
t957 = t513 * t937;
t956 = t528 * t936;
t955 = t532 * t936;
t923 = t425 * mrSges(5,3);
t954 = t936 * t923;
t952 = t937 * t936;
t645 = pkin(2) * t752;
t511 = t529 * t645;
t465 = t526 * t518 + t511;
t460 = pkin(9) + t465;
t810 = pkin(10) + t460;
t440 = t810 * t528;
t441 = t810 * t532;
t340 = -t440 * t527 + t441 * t531;
t609 = -t531 * t440 - t441 * t527;
t950 = -t340 * mrSges(7,1) - t609 * mrSges(7,2);
t821 = pkin(3) * t526;
t512 = pkin(9) + t821;
t809 = pkin(10) + t512;
t475 = t809 * t528;
t476 = t809 * t532;
t415 = -t475 * t527 + t476 * t531;
t608 = -t531 * t475 - t476 * t527;
t949 = -t415 * mrSges(7,1) - t608 * mrSges(7,2);
t946 = t937 * mrSges(5,3) + t226 * t972;
t945 = -t464 * t937 + t465 * t936;
t944 = t526 * t936 - t752 * t937;
t79 = t121 * t531 - t750;
t943 = t67 - t79;
t685 = Ifges(7,5) * t493 + Ifges(7,6) * t494;
t80 = t685 + t950;
t942 = t80 * qJD(6);
t92 = t685 + t949;
t941 = t92 * qJD(6);
t819 = pkin(5) * t528;
t876 = t430 * t819 + t506 * t931 + t894 * t930 + t920;
t822 = pkin(3) * t496;
t282 = pkin(4) * t580 - pkin(9) * t425 - t822;
t155 = t532 * t282 - t956;
t604 = pkin(5) * t580 - pkin(10) * t727;
t114 = t155 + t604;
t156 = t528 * t282 + t955;
t680 = pkin(10) * t728;
t124 = -t680 + t156;
t597 = t114 * t527 + t124 * t531;
t598 = t114 * t531 - t124 * t527;
t927 = t597 * t962 + t598 * t963;
t523 = t530 * pkin(2);
t270 = t282 + t523;
t151 = t532 * t270 - t956;
t111 = t151 + t604;
t152 = t528 * t270 + t955;
t123 = -t680 + t152;
t69 = t111 * t531 - t123 * t527;
t70 = t111 * t527 + t123 * t531;
t928 = t69 * t963 + t70 * t962;
t173 = mrSges(7,1) * t299 + mrSges(7,2) * t301;
t219 = -mrSges(7,1) * t425 - t299 * mrSges(7,3);
t836 = t448 / 0.2e1;
t934 = -t340 / 0.2e1;
t938 = t173 * t836 + t219 * t934;
t851 = t340 / 0.2e1;
t847 = t415 / 0.2e1;
t933 = t425 / 0.2e1;
t831 = t494 / 0.2e1;
t929 = Ifges(4,2) - Ifges(4,1);
t78 = -t121 * t527 - t749;
t908 = t68 + t78;
t766 = Ifges(6,5) * t425;
t524 = t528 ^ 2;
t525 = t532 ^ 2;
t684 = t524 + t525;
t429 = -t494 * mrSges(7,1) + t493 * mrSges(7,2);
t916 = t429 * qJD(6);
t433 = Ifges(7,1) * t493 + t797;
t616 = t432 / 0.4e1 - t433 / 0.4e1;
t431 = Ifges(7,2) * t494 + t482;
t617 = t431 / 0.4e1 + t434 / 0.4e1;
t915 = -t936 * t502 / 0.2e1 + t617 * t301 - t616 * t299;
t781 = t301 * mrSges(7,3);
t217 = mrSges(7,2) * t425 + t781;
t914 = t219 * t932 - t493 * t217 / 0.2e1;
t792 = Ifges(7,6) * t300;
t795 = Ifges(7,5) * t302;
t888 = -t792 / 0.2e1 - t795 / 0.2e1;
t913 = -mrSges(4,3) * t445 - Ifges(4,4) * t496;
t910 = -m(6) / 0.2e1;
t859 = -t219 / 0.2e1;
t909 = Ifges(6,3) + Ifges(7,3);
t907 = Ifges(5,2) * t580;
t757 = t494 * mrSges(7,3);
t866 = m(7) / 0.2e1;
t902 = t500 * t866;
t472 = t533 * t645 - t510;
t579 = t493 * t472;
t573 = -t579 / 0.2e1;
t389 = t494 * t472;
t775 = t389 * mrSges(7,1);
t899 = t775 / 0.2e1 + mrSges(7,2) * t573;
t520 = Ifges(6,5) * t532;
t793 = Ifges(6,6) * t528;
t895 = t520 - t793;
t893 = mrSges(7,2) * t579 / 0.2e1 - t775 / 0.2e1;
t322 = mrSges(6,2) * t425 - mrSges(6,3) * t726;
t693 = t532 * t322;
t324 = -mrSges(6,1) * t425 - mrSges(6,3) * t725;
t701 = t528 * t324;
t892 = -t693 / 0.2e1 + t701 / 0.2e1;
t698 = t531 * t218;
t708 = t527 * t216;
t891 = -t698 / 0.2e1 - t708 / 0.2e1;
t850 = -t608 / 0.2e1;
t890 = -t299 * t847 + t301 * t850;
t889 = -t754 / 0.2e1 - t755 / 0.2e1;
t887 = -t520 / 0.2e1 + t793 / 0.2e1;
t311 = t505 * t580;
t865 = pkin(5) / 0.2e1;
t886 = t175 * t865 - t311 / 0.4e1;
t615 = -t503 / 0.4e1 + t506 / 0.4e1;
t841 = t430 / 0.2e1;
t618 = t841 + t885 / 0.2e1;
t883 = t175 + t309 + t905;
t789 = Ifges(7,3) * t580;
t882 = t789 / 0.2e1 - t888;
t881 = (-t299 * t831 + t301 * t833) * mrSges(7,3) + t914;
t880 = (t300 * t831 + t302 * t833) * mrSges(7,3);
t591 = t782 / 0.2e1 - t780 / 0.2e1;
t24 = t591 + t881;
t879 = t24 * qJD(1);
t649 = -t757 / 0.2e1;
t758 = t493 * mrSges(7,3);
t651 = -t758 / 0.2e1;
t878 = -t299 * t649 + t301 * t651 - t914;
t877 = t888 - t928;
t654 = -t789 / 0.2e1;
t875 = t654 + t888 - t927;
t177 = Ifges(7,1) * t301 - t798;
t843 = -t425 / 0.4e1;
t874 = t206 * t429 / 0.2e1 + t685 * t843 + (-t177 / 0.4e1 + t135 / 0.4e1) * t494;
t605 = t432 * t831 + t433 * t932 + (t431 + t434) * t833;
t606 = t216 * t932 + t218 * t833 + t321 * t931 + t323 * t930;
t203 = t206 * t819;
t829 = t500 / 0.2e1;
t629 = t173 * t829;
t849 = t608 / 0.2e1;
t307 = t885 * t580;
t858 = -t307 / 0.2e1;
t873 = (t908 * t608 + t203) * t866 + t513 * t858 + t629 + t217 * t849 + (-t866 * t943 + t859) * t415;
t635 = t727 / 0.2e1;
t636 = -t728 / 0.2e1;
t790 = Ifges(6,3) * t580;
t872 = Ifges(6,5) * t635 + Ifges(6,6) * t636 + t790 / 0.2e1 + t882 + (t698 + t708) * t865;
t870 = -m(5) / 0.2e1;
t869 = m(5) / 0.2e1;
t868 = m(6) / 0.2e1;
t867 = -m(7) / 0.2e1;
t864 = m(5) * pkin(3);
t860 = t217 / 0.2e1;
t310 = t580 * t503;
t857 = -t310 / 0.4e1;
t854 = -t609 / 0.2e1;
t853 = t609 / 0.2e1;
t848 = -t415 / 0.2e1;
t844 = -t425 / 0.2e1;
t835 = -t460 / 0.2e1;
t825 = -t512 / 0.2e1;
t820 = pkin(5) * t527;
t817 = t67 * mrSges(7,2);
t816 = t68 * mrSges(7,1);
t813 = t78 * mrSges(7,1);
t812 = t79 * mrSges(7,2);
t808 = m(7) * qJD(4);
t805 = mrSges(6,3) * t425;
t804 = mrSges(6,3) * t580;
t796 = Ifges(7,5) * t299;
t791 = Ifges(7,6) * t301;
t788 = pkin(5) * qJD(5);
t765 = Ifges(6,6) * t425;
t227 = t580 * t894 - t765;
t229 = t506 * t580 - t766;
t326 = -mrSges(5,1) * t425 + mrSges(5,2) * t580;
t458 = t523 - t822;
t483 = Ifges(4,4) * t495;
t600 = t792 + t795;
t633 = t725 / 0.2e1;
t3 = -t954 + (-t907 - 0.2e1 * t951) * t933 - (-t960 + t907 / 0.2e1 + t909 * t846) * t425 + (-Ifges(5,4) * t580 + t895 * t846 + t946) * t580 + (-mrSges(4,2) * t523 + t495 * t929 + t913) * t496 + t228 * t633 + t229 * t635 + t227 * t636 + (-pkin(2) * mrSges(4,1) * t495 - pkin(1) * mrSges(3,1) - Ifges(3,4) * t530) * t530 + (t791 + t796) * t846 + t152 * t322 + t151 * t324 + (t425 * t895 + t600 + t789 + t790) * t844 + m(4) * t519 * t523 + t70 * t217 + t69 * t219 - t966 + t483 * t495 + m(7) * (t67 * t69 + t68 * t70 + t970) + m(6) * (t149 * t151 + t150 * t152 - t952) + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t534 + (Ifges(3,1) - Ifges(3,2)) * t530) * t534 + (t971 + t326) * t458;
t783 = t3 * qJD(1);
t695 = t532 * t229;
t704 = t528 * t227;
t4 = (-t496 * t929 - t483) * t495 + (pkin(3) * t326 - t913) * t496 - m(7) * (t597 * t68 + t598 * t67 + t970) + t822 * t971 - m(6) * (t149 * t155 + t150 * t156 - t952) - t954 - t597 * t217 - t598 * t219 + t600 * t933 - t156 * t322 - t155 * t324 + (-t695 / 0.2e1 + t704 / 0.2e1 + t960 + t895 * t933 + t951) * t425 + (-t796 / 0.2e1 - t791 / 0.2e1 - (-Ifges(5,2) - t909) * t425 - t532 * t228 / 0.2e1 + (t887 + Ifges(5,4)) * t580 - t946) * t580 + t966;
t774 = t4 * qJD(1);
t176 = -Ifges(7,2) * t299 + t296;
t688 = Ifges(7,5) * t301 - Ifges(7,6) * t299;
t548 = t206 * t173 + (t137 / 0.2e1 + t176 / 0.2e1) * t301 - (-t177 / 0.2e1 + t68 * mrSges(7,3) + t135 / 0.2e1) * t299 - t67 * t781 + t688 * t844;
t9 = t78 * t219 + t79 * t217 + m(7) * (t67 * t78 + t68 * t79) + t149 * t322 - t150 * t324 + t936 * t307 + ((t149 * mrSges(6,3) + t766 / 0.2e1 - t229 / 0.2e1 + t310 / 0.2e1) * t528 + (-t150 * mrSges(6,3) + t765 / 0.2e1 - t311 / 0.2e1 - t227 / 0.2e1 + (m(7) * t206 + t175) * pkin(5)) * t532) * t580 + t548;
t753 = t9 * qJD(1);
t612 = t684 * t460;
t614 = t525 / 0.2e1 + t524 / 0.2e1;
t541 = t880 + t614 * t805 + (t425 * t465 - t464 * t580) * t869 + (t425 * t612 + t459 * t580) * t868 + (t300 * t609 + t302 * t340 + t448 * t580) * t866;
t552 = t458 * t869 + (t151 * t532 + t152 * t528) * t868 + (t493 * t69 - t494 * t70) * t866;
t559 = -t417 - t606;
t20 = (-mrSges(5,1) + t618) * t580 + t541 - t552 + t559;
t751 = qJD(1) * t20;
t14 = t67 * t217 - t68 * t219 + t548;
t748 = t14 * qJD(1);
t747 = t151 * t528;
t746 = t152 * t532;
t745 = t155 * t528;
t744 = t156 * t532;
t682 = m(7) * t865;
t545 = t889 * t425 + (t300 * t531 + t302 * t527) * t682 + t591;
t553 = (t908 * t493 + t494 * t943) * t867 + t892;
t16 = t545 + t553 + t881;
t743 = t16 * qJD(1);
t596 = -t149 * t528 + t150 * t532;
t734 = t936 * t580;
t19 = t302 * t217 + t300 * t219 + (t693 - t701 + t923) * t425 + t883 * t580 + m(7) * (t206 * t580 + t300 * t67 + t302 * t68) + m(6) * (t425 * t596 - t734) + m(5) * (t425 * t937 - t734);
t742 = t19 * qJD(1);
t471 = t526 * t823 + t511;
t733 = t936 * t471;
t732 = t609 * t218;
t731 = t340 * t216;
t730 = t608 * t218;
t729 = t415 * t216;
t724 = t448 * t174;
t723 = t459 * t308;
t714 = t500 * t174;
t713 = t500 * t429;
t712 = t512 * t528;
t711 = t513 * t308;
t710 = t513 * t502;
t707 = t527 * t299;
t702 = t528 * t323;
t697 = t531 * t301;
t694 = t532 * t321;
t683 = t864 / 0.2e1;
t681 = t448 * t819;
t224 = t493 * t389 - t494 * t579;
t671 = t224 * t866;
t679 = qJD(3) * t671;
t678 = mrSges(6,3) * t747;
t677 = mrSges(6,3) * t746;
t676 = mrSges(6,3) * t745;
t675 = mrSges(6,3) * t744;
t318 = t340 * t757;
t674 = t464 * t923;
t673 = t465 * t905;
t672 = t531 * t758;
t667 = -t67 / 0.2e1 + t79 / 0.2e1;
t666 = t78 / 0.2e1 + t68 / 0.2e1;
t665 = t460 * t702;
t664 = t460 * t694;
t663 = t512 * t694;
t652 = t765 / 0.4e1;
t650 = t758 / 0.2e1;
t648 = t757 / 0.2e1;
t628 = t324 * t825;
t627 = -t712 / 0.2e1;
t621 = t229 / 0.4e1 + t857;
t620 = t850 + t849;
t619 = t848 + t847;
t611 = t684 * t512;
t610 = t905 * t821;
t599 = t644 * t923;
t595 = t746 - t747;
t594 = t744 - t745;
t593 = t614 * t804;
t592 = -t151 * mrSges(6,1) / 0.2e1 + t152 * mrSges(6,2) / 0.2e1;
t350 = t448 * t429;
t590 = -t318 - t350 - t605;
t589 = m(7) * (t527 * t70 + t531 * t69);
t584 = t757 * t820 + t685 + t895;
t538 = (t300 * t608 + t302 * t415) * t866 + t684 * t805 / 0.2e1 + (t526 * t683 + t611 * t868) * t425 + (t513 * t868 - t683 * t752 + t618 + t902) * t580 + t880;
t544 = (t155 * t532 + t156 * t528) * t868 + (t493 * t598 - t494 * t597) * t866 + t822 * t870;
t22 = t538 - t544 + t559 - t769;
t581 = -t22 * qJD(1) + qJD(2) * t671;
t75 = t605 + t713;
t578 = t598 * t757;
t577 = t597 * t758;
t317 = t609 * t758;
t439 = t459 * t502;
t45 = -t317 - m(7) * t681 - t439 + (t340 * t494 + t493 * t609) * mrSges(7,3) + t590 - t876;
t549 = (t137 / 0.4e1 + t176 / 0.4e1) * t493 + t874;
t543 = t520 * t843 + t549 + t915;
t566 = t493 * t667 + t494 * t666;
t537 = (-t299 * t851 + t301 * t854 + t566) * mrSges(7,3) - t460 * t593 + (-t340 * t943 + t908 * t609 + t203) * t866 + t217 * t853 + t459 * t858 + t543 + t938;
t551 = t324 * t835 + ((m(7) * t836 + t841) * pkin(5) + t615) * t580 + t621;
t563 = (-t505 / 0.4e1 - t894 / 0.4e1) * t580 - t227 / 0.4e1 + t886;
t556 = t322 * t835 + t563;
t5 = t537 + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1) * t580 + (-t589 / 0.2e1 + t891) * pkin(5) + (-t766 / 0.2e1 + t551) * t532 + ((t933 + t425 / 0.4e1) * Ifges(6,6) + t556) * t528 + t592 + t877;
t569 = -t5 * qJD(1) + t45 * qJD(2);
t535 = (t460 * t594 - t733 + t958) * t910 + (-t733 + t945) * t870 - t731 / 0.2e1 + t673 / 0.2e1 + t676 / 0.2e1 - t724 / 0.2e1 - t577 / 0.2e1 + t217 * t573 + t389 * t859 - t732 / 0.2e1 + t674 / 0.2e1 - t723 / 0.2e1 - t578 / 0.2e1 - t883 * t471 / 0.2e1 + (t937 * t870 + t596 * t910 - t923 / 0.2e1 + t892) * t472 - t664 / 0.2e1 + t665 / 0.2e1 + (t471 * t206 + t340 * t597 + t389 * t67 + t579 * t68 + t598 * t609 + t969) * t867 - t675 / 0.2e1 - t973;
t536 = (t512 * t595 + t957) * t868 - t599 / 0.2e1 + t944 * t683 - t678 / 0.2e1 + t69 * t648 + t70 * t650 + t323 * t627 + t677 / 0.2e1 + t730 / 0.2e1 + t729 / 0.2e1 + t714 / 0.2e1 + t711 / 0.2e1 + t663 / 0.2e1 + (t415 * t70 + t608 * t69 + t968) * t866 - t610 / 0.2e1 + t973;
t2 = t535 + t536;
t550 = mrSges(4,1) * t824 + mrSges(4,2) * t823 - t389 * t757 - t579 * t758 + (mrSges(5,1) - t430 - t885) * t471 + (-mrSges(6,3) * t684 + mrSges(5,2)) * t472;
t47 = -m(7) * (t340 * t579 + t389 * t609 + t448 * t471) - m(6) * (t459 * t471 + t472 * t612) - m(5) * (-t464 * t471 + t465 * t472) + t550;
t568 = t2 * qJD(1) + t47 * qJD(2) - t224 * t808 / 0.2e1;
t539 = -(mrSges(7,3) * t851 + t616) * t299 + (mrSges(7,3) * t854 + t617) * t301 + t609 * t860 + t549 + t938;
t11 = t654 + t539 + t877;
t49 = t590 + t318;
t567 = -t11 * qJD(1) + t49 * qJD(2);
t565 = t318 / 0.2e1 + t350 / 0.2e1 + t713 / 0.2e1 + t605;
t473 = t500 * t819;
t542 = t317 / 0.2e1 + t439 / 0.2e1 + (t473 + t681) * t866 + t710 / 0.2e1 + t565 + t876;
t547 = (t531 * t389 + t527 * t579) * t682 + t889 * t472;
t35 = t608 * t650 + t542 - t547 + (t609 + t608) * t651 - t340 * t648 + t893;
t48 = m(7) * t473 + t710 + t75 + t876;
t560 = t155 * mrSges(6,1) / 0.2e1 - t156 * mrSges(6,2) / 0.2e1 + (t527 ^ 2 + t531 ^ 2) * t114 * t682;
t7 = (t628 + t857) * t532 + t873 + t908 * t648 + t915 + t874 + t875 + (t430 * t633 + t725 * t902 + t891) * pkin(5) + t79 * t650 + t67 * t651 + t322 * t627 - (t505 + t894) * t726 / 0.4e1 + t895 * t843 - t560 - t790 / 0.2e1 + t615 * t725 + t886 * t528 + t887 * t425 + t695 / 0.4e1 + t890 * mrSges(7,3) + t684 * t804 * t825 + (t176 + t137) * t493 / 0.4e1 - t704 / 0.4e1;
t562 = -t7 * qJD(1) - t35 * qJD(2) - t48 * qJD(3);
t540 = (mrSges(7,3) * t850 + t617) * t301 - (mrSges(7,3) * t847 + t616) * t299 + t608 * t860 + t219 * t848 + t629 + t549;
t13 = t540 + t875;
t557 = t340 * t649 + t565;
t39 = t557 + t893;
t561 = -t13 * qJD(1) - t39 * qJD(2) - t75 * qJD(3);
t106 = mrSges(7,1) * t619 + mrSges(7,2) * t620;
t554 = (t527 * t859 + t531 * t860 + (-t697 / 0.2e1 - t707 / 0.2e1) * mrSges(7,3)) * pkin(5);
t18 = -mrSges(7,1) * t666 + mrSges(7,2) * t667 + t554;
t499 = (t527 * mrSges(7,1) + t531 * mrSges(7,2)) * pkin(5);
t83 = (t854 + t853) * mrSges(7,2) + (t934 + t851) * mrSges(7,1);
t555 = -t18 * qJD(1) - t83 * qJD(2) - t106 * qJD(3) + t499 * qJD(5);
t492 = t499 * qJD(6);
t38 = t557 + t899;
t34 = ((t934 - t619) * t494 + (t854 + t620) * t493) * mrSges(7,3) + t542 + t547 + t899;
t25 = t591 + t878;
t23 = t538 + t544 + t606;
t21 = t580 * t618 + t541 + t552 + t606;
t17 = t545 - t553 + t878;
t15 = -t817 / 0.2e1 - t816 / 0.2e1 - t812 / 0.2e1 + t813 / 0.2e1 + t554 + t688;
t12 = t540 + t882 + t927;
t10 = t539 + t882 + t928;
t8 = t872 + (t566 + t890) * mrSges(7,3) + t560 + t543 - t512 * t593 + (t628 + ((m(7) * t829 + t841) * pkin(5) + t615) * t580 + t621) * t532 + (t322 * t825 + t563 + t652) * t528 + t873 + t927;
t6 = (t652 + t556) * t528 + t589 * t865 + t537 + t551 * t532 - t592 + t872 + t928;
t1 = t558 - t535 + t536;
t26 = [qJD(2) * t3 - qJD(3) * t4 + qJD(4) * t19 + qJD(5) * t9 + qJD(6) * t14, t783 + (t731 - t673 + t724 - t678 + t677 + Ifges(3,5) * t534 - Ifges(3,6) * t530 - t495 * mrSges(4,3) * t823 + t732 - t674 + t723 + (-mrSges(3,1) * t534 + mrSges(3,2) * t530) * pkin(7) + t664 - t665 + t69 * t757 + t70 * t758 + t756 * t824 + t974) * qJD(2) + t1 * qJD(3) + t21 * qJD(4) + t6 * qJD(5) + t10 * qJD(6) + 0.2e1 * ((t460 * t595 + t958) * t868 + t945 * t869 + (t340 * t70 + t609 * t69 + t969) * t866 + m(4) * (-t445 * t533 + t529 * t896) * pkin(2) / 0.2e1) * qJD(2), -t774 + t1 * qJD(2) + (m(6) * (t512 * t594 + t957) - t599 + t944 * t864 - t676 + t577 + t730 + t729 + t578 + t714 + t711 + t663 + m(7) * (t415 * t597 + t598 * t608 + t968) + t675 - t610 - t512 * t702 + t974) * qJD(3) + t23 * qJD(4) + t8 * qJD(5) + t12 * qJD(6), t742 + t21 * qJD(2) + t23 * qJD(3) + (t300 * t493 - t302 * t494) * t808 + t17 * qJD(5) + t25 * qJD(6), t753 + t6 * qJD(2) + t8 * qJD(3) + t17 * qJD(4) + (-t150 * mrSges(6,1) - t149 * mrSges(6,2) - Ifges(6,5) * t726 - Ifges(6,6) * t725 + t688 - t812 + t813) * qJD(5) + t15 * qJD(6) + (m(7) * (t527 * t79 + t531 * t78) + (-t697 - t707) * mrSges(7,3)) * t788, t748 + t10 * qJD(2) + t12 * qJD(3) + t25 * qJD(4) + t15 * qJD(5) + (t688 - t816 - t817) * qJD(6); -qJD(3) * t2 + qJD(4) * t20 + qJD(5) * t5 + qJD(6) * t11 - t783, -qJD(3) * t47 - qJD(5) * t45 - qJD(6) * t49 (m(7) * (t389 * t608 + t415 * t579 + t500 * t471) + m(6) * (t471 * t513 + t472 * t611) + (-t471 * t752 + t472 * t526) * t864 - t550) * qJD(3) + t34 * qJD(5) + t38 * qJD(6) - t568, t679 + t751, t34 * qJD(3) + (t460 * t885 + t584 + t950) * qJD(5) + t942 + (m(7) * (-t340 * t531 + t527 * t609) - t672) * t788 - t569, t38 * qJD(3) + t80 * qJD(5) - t567 + t942; qJD(2) * t2 + qJD(4) * t22 + qJD(5) * t7 + qJD(6) * t13 + t774, t35 * qJD(5) + t39 * qJD(6) + t568, qJD(5) * t48 + qJD(6) * t75, -t581 (mrSges(6,2) * t712 - t512 * t807 + t584 + t949) * qJD(5) + t941 + (m(7) * (-t415 * t531 + t527 * t608) - t672) * t788 - t562, t92 * qJD(5) - t561 + t941; -qJD(2) * t20 - qJD(3) * t22 - qJD(5) * t16 - qJD(6) * t24 - t742, t679 - t751, t581, 0, -t743 - t916 + (-t429 - t502 + (t531 * pkin(5) * t494 + t493 * t820) * m(7)) * qJD(5), -qJD(5) * t429 - t879 - t916; -qJD(2) * t5 - qJD(3) * t7 + qJD(4) * t16 + qJD(6) * t18 - t753, -t35 * qJD(3) + t83 * qJD(6) + t569, t106 * qJD(6) + t562, t743, -t492, -t492 - t555; -qJD(2) * t11 - qJD(3) * t13 + qJD(4) * t24 - qJD(5) * t18 - t748, -t39 * qJD(3) - t83 * qJD(5) + t567, -t106 * qJD(5) + t561, t879, t555, 0;];
Cq  = t26;
