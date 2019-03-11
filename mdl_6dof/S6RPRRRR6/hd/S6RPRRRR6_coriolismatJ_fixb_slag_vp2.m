% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:24
% EndTime: 2019-03-09 07:13:09
% DurationCPUTime: 28.30s
% Computational Cost: add. (66123->844), mult. (131178->1124), div. (0->0), fcn. (154263->10), ass. (0->452)
t971 = qJD(4) + qJD(5);
t587 = sin(qJ(4));
t907 = -pkin(9) - pkin(8);
t570 = t907 * t587;
t591 = cos(qJ(4));
t571 = t907 * t591;
t586 = sin(qJ(5));
t590 = cos(qJ(5));
t515 = t570 * t586 - t571 * t590;
t637 = t586 * t587 - t590 * t591;
t448 = -pkin(10) * t637 + t515;
t585 = sin(qJ(6));
t589 = cos(qJ(6));
t638 = t586 * t591 + t587 * t590;
t949 = t570 * t590 + t571 * t586;
t972 = -pkin(10) * t638 + t949;
t1005 = -t448 * t585 + t589 * t972;
t1009 = t1005 * mrSges(7,2);
t325 = t448 * t589 + t585 * t972;
t1015 = t325 * mrSges(7,1);
t1019 = -t1015 / 0.2e1 - t1009 / 0.2e1;
t666 = -t585 * t638 - t589 * t637;
t956 = Ifges(7,5) * t666;
t493 = -t585 * t637 + t589 * t638;
t976 = Ifges(7,6) * t493;
t981 = t956 - t976;
t1023 = 0.2e1 * t1019 + t981;
t1024 = t1023 * qJD(6);
t365 = -mrSges(7,1) * t666 + mrSges(7,2) * t493;
t578 = -pkin(4) * t591 - pkin(3);
t518 = pkin(5) * t637 + t578;
t1022 = m(7) * t518 + t365;
t1010 = t1005 / 0.2e1;
t893 = -t325 / 0.2e1;
t896 = -t1005 / 0.2e1;
t1011 = (t896 + t1010) * mrSges(7,2) + (t893 + t325 / 0.2e1) * mrSges(7,1);
t1021 = qJD(6) * t1011;
t548 = Ifges(6,6) * t638;
t549 = Ifges(6,5) * t637;
t660 = -t549 - t548 + t981;
t1000 = -mrSges(6,1) * t515 - mrSges(6,2) * t949 + t660;
t1006 = -t1015 - t1009;
t1020 = t1000 + t1006;
t975 = t493 * mrSges(7,1);
t695 = -t975 / 0.2e1;
t879 = t666 / 0.2e1;
t960 = -t666 / 0.2e1;
t622 = t695 + t975 / 0.2e1 + (t960 + t879) * mrSges(7,2);
t1018 = qJD(2) * t622 + qJD(3) * t1011;
t1014 = t325 * t493;
t1013 = t1005 * t585 - t325 * t589;
t583 = sin(pkin(11));
t584 = cos(pkin(11));
t588 = sin(qJ(3));
t850 = cos(qJ(3));
t553 = t583 * t588 - t584 * t850;
t555 = t583 * t850 + t584 * t588;
t449 = t637 * t555;
t451 = t638 * t555;
t667 = t449 * t585 - t451 * t589;
t273 = -mrSges(7,2) * t553 + mrSges(7,3) * t667;
t329 = t449 * t589 + t451 * t585;
t343 = -mrSges(6,1) * t449 - mrSges(6,2) * t451;
t829 = pkin(7) + qJ(2);
t563 = t829 * t583;
t564 = t829 * t584;
t506 = t563 * t850 + t564 * t588;
t745 = t555 * t587;
t423 = pkin(4) * t745 + t506;
t496 = mrSges(6,1) * t638 - mrSges(6,2) * t637;
t856 = t578 / 0.2e1;
t275 = mrSges(7,1) * t553 + mrSges(7,3) * t329;
t902 = t275 / 0.2e1;
t483 = Ifges(7,4) * t666;
t371 = Ifges(7,1) * t493 + t483;
t964 = -Ifges(7,2) * t493 + t483;
t979 = t371 / 0.4e1 + t964 / 0.4e1;
t824 = Ifges(7,4) * t329;
t172 = Ifges(7,2) * t667 + Ifges(7,6) * t553 - t824;
t319 = Ifges(7,4) * t667;
t173 = -Ifges(7,1) * t329 + Ifges(7,5) * t553 + t319;
t327 = pkin(5) * t451 + t423;
t965 = Ifges(7,2) * t329 + t319;
t967 = -mrSges(7,1) * t329 + mrSges(7,2) * t667;
t954 = t666 * mrSges(7,2);
t980 = t975 + t954;
t982 = Ifges(7,1) * t667 + t824;
t995 = t967 * t518 / 0.2e1 + (t173 / 0.4e1 + t965 / 0.4e1) * t666 + t980 * t327 / 0.2e1 + (t982 / 0.4e1 - t172 / 0.4e1) * t493;
t787 = t493 * Ifges(7,4);
t368 = Ifges(7,2) * t666 + t787;
t983 = Ifges(7,1) * t666 - t787;
t996 = -t983 / 0.4e1 + t368 / 0.4e1;
t1012 = t273 * t1010 - t325 * t902 + t343 * t856 + t423 * t496 / 0.2e1 + t979 * t667 + t996 * t329 + t995;
t841 = pkin(8) * t553;
t847 = pkin(3) * t555;
t495 = t841 + t847;
t384 = t495 * t591 + t506 * t587;
t746 = t553 * t591;
t282 = pkin(4) * t555 + pkin(9) * t746 + t384;
t385 = t495 * t587 - t506 * t591;
t747 = t553 * t587;
t342 = pkin(9) * t747 + t385;
t174 = t282 * t590 - t342 * t586;
t452 = t637 * t553;
t122 = pkin(5) * t555 - pkin(10) * t452 + t174;
t450 = t638 * t553;
t330 = t450 * t589 - t452 * t585;
t272 = -mrSges(7,2) * t555 + mrSges(7,3) * t330;
t333 = t450 * t585 + t452 * t589;
t274 = mrSges(7,1) * t555 - mrSges(7,3) * t333;
t843 = pkin(5) * t585;
t705 = t843 / 0.2e1;
t842 = pkin(5) * t589;
t917 = m(7) * pkin(5);
t636 = -Ifges(7,5) * t333 / 0.2e1 - Ifges(7,6) * t330 / 0.2e1;
t175 = t282 * t586 + t342 * t590;
t129 = pkin(10) * t450 + t175;
t641 = t122 * t585 + t129 * t589;
t642 = t122 * t589 - t129 * t585;
t603 = t636 + t641 * mrSges(7,2) / 0.2e1 - t642 * mrSges(7,1) / 0.2e1;
t860 = t555 / 0.2e1;
t654 = Ifges(7,3) * t860 - t603;
t945 = -Ifges(6,6) * t450 / 0.2e1 - Ifges(6,5) * t452 / 0.2e1;
t931 = t945 + t175 * mrSges(6,2) / 0.2e1 - t174 * mrSges(6,1) / 0.2e1;
t930 = Ifges(6,3) * t860 + t654 - t931;
t1008 = t272 * t705 + t274 * t842 / 0.2e1 + (t585 ^ 2 + t589 ^ 2) * t122 * t917 / 0.2e1 + t930;
t994 = -t917 / 0.2e1;
t1002 = t518 * t980;
t690 = -pkin(2) * t584 - pkin(1);
t468 = pkin(3) * t553 - pkin(8) * t555 + t690;
t507 = -t563 * t588 + t564 * t850;
t357 = t468 * t591 - t507 * t587;
t744 = t555 * t591;
t312 = -pkin(9) * t744 + t357;
t271 = pkin(4) * t553 + t312;
t358 = t468 * t587 + t507 * t591;
t313 = -pkin(9) * t745 + t358;
t292 = t586 * t313;
t154 = t271 * t590 - t292;
t182 = t312 * t590 - t292;
t1001 = t154 - t182;
t797 = t333 * mrSges(7,2);
t798 = t330 * mrSges(7,1);
t722 = t798 / 0.2e1 - t797 / 0.2e1;
t792 = t452 * mrSges(6,2);
t794 = t450 * mrSges(6,1);
t662 = t722 - t792 / 0.2e1 + t794 / 0.2e1;
t999 = (t330 * t589 + t333 * t585) * t994 - t662;
t442 = t449 * pkin(10);
t126 = t154 + t442;
t119 = pkin(5) * t553 + t126;
t294 = t590 * t313;
t155 = t271 * t586 + t294;
t840 = pkin(10) * t451;
t127 = t155 - t840;
t769 = t127 * t589;
t76 = t119 * t585 + t769;
t997 = t76 * mrSges(7,3) - t982 / 0.2e1;
t993 = t965 / 0.2e1;
t990 = t327 * t967;
t779 = t638 * mrSges(6,3);
t988 = t515 * t779;
t181 = -t312 * t586 - t294;
t985 = t155 + t181;
t579 = Ifges(5,5) * t591;
t948 = -Ifges(5,6) * t587 + t579;
t955 = Ifges(7,5) * t667;
t977 = Ifges(7,6) * t329;
t721 = t955 + t977;
t671 = t977 / 0.2e1 + t955 / 0.2e1;
t891 = t329 / 0.2e1;
t793 = t451 * mrSges(6,3);
t396 = -mrSges(6,2) * t553 - t793;
t884 = t396 / 0.2e1;
t796 = t449 * mrSges(6,3);
t398 = mrSges(6,1) * t553 + t796;
t883 = -t398 / 0.2e1;
t876 = t493 / 0.2e1;
t857 = t638 / 0.2e1;
t973 = Ifges(4,4) - t948;
t970 = t1013 * t994 - t1019;
t134 = t181 + t840;
t135 = t442 + t182;
t92 = t134 * t585 + t135 * t589;
t830 = t92 * mrSges(7,2);
t91 = t134 * t589 - t135 * t585;
t831 = t91 * mrSges(7,1);
t827 = t831 / 0.2e1 - t830 / 0.2e1;
t969 = (t585 * t92 + t589 * t91) * t994 - t827;
t770 = t127 * t585;
t75 = t119 * t589 - t770;
t908 = -t76 / 0.2e1;
t615 = t325 * t891 + t493 * t908 + t667 * t896 + t75 * t960;
t82 = -t126 * t585 - t769;
t83 = t126 * t589 - t770;
t838 = t449 * pkin(5);
t844 = pkin(5) * t638;
t778 = t638 * Ifges(6,4);
t500 = -Ifges(6,1) * t637 - t778;
t872 = -t500 / 0.4e1;
t499 = -Ifges(6,2) * t637 + t778;
t873 = t499 / 0.4e1;
t877 = -t493 / 0.2e1;
t888 = -t365 / 0.2e1;
t795 = t449 * Ifges(6,4);
t346 = -Ifges(6,1) * t451 + t795;
t889 = -t346 / 0.4e1;
t283 = -Ifges(6,2) * t451 + Ifges(6,6) * t553 - t795;
t900 = t283 / 0.4e1;
t921 = m(7) / 0.2e1;
t861 = t553 / 0.4e1;
t947 = t660 * t861;
t185 = -mrSges(7,1) * t667 - mrSges(7,2) * t329;
t961 = t185 / 0.2e1;
t966 = (t82 * t877 + t83 * t879 + t615) * mrSges(7,3) + t844 * t961 + ((t327 * t638 - t449 * t518) * pkin(5) + (t82 + t76) * t1005 + (-t75 + t83) * t325) * t921 + t449 * t872 + t449 * t873 + t949 * t884 + t838 * t888 - t638 * t889 - t638 * t900 + t947 + t1012;
t775 = t591 * mrSges(5,2);
t777 = t587 * mrSges(5,1);
t565 = t775 + t777;
t953 = t506 * t565;
t550 = Ifges(6,4) * t637;
t498 = -Ifges(6,2) * t638 - t550;
t501 = Ifges(6,1) * t638 - t550;
t669 = t498 / 0.4e1 + t501 / 0.4e1;
t952 = t669 * t451;
t439 = Ifges(6,4) * t451;
t284 = -Ifges(6,1) * t449 + Ifges(6,5) * t553 - t439;
t345 = Ifges(6,2) * t449 - t439;
t951 = t345 + t284;
t870 = -t949 / 0.2e1;
t882 = t449 / 0.2e1;
t946 = -t451 * t870 + t515 * t882;
t661 = -Ifges(6,5) * t451 + Ifges(6,6) * t449 + t721;
t944 = -t384 * t587 + t385 * t591;
t580 = Ifges(5,4) * t591;
t820 = Ifges(5,2) * t587;
t943 = t820 - t580;
t776 = t587 * mrSges(5,2);
t826 = mrSges(5,1) * t591;
t942 = -t826 + t776;
t926 = t591 ^ 2;
t927 = t587 ^ 2;
t940 = t926 + t927;
t934 = qJD(6) * t622;
t624 = -t954 + 0.2e1 * t695;
t933 = t624 * qJD(6);
t835 = t76 * mrSges(7,1);
t837 = t75 * mrSges(7,2);
t932 = -t835 / 0.2e1 - t837 / 0.2e1 + t671;
t929 = t273 * t960 + t275 * t876 + t398 * t857 + t637 * t884;
t925 = m(5) / 0.2e1;
t924 = -m(6) / 0.2e1;
t923 = m(6) / 0.2e1;
t922 = -m(7) / 0.2e1;
t919 = -pkin(8) / 0.2e1;
t918 = m(6) * pkin(4);
t916 = -mrSges(5,1) / 0.2e1;
t915 = mrSges(5,2) / 0.2e1;
t914 = -mrSges(5,3) / 0.2e1;
t911 = -mrSges(7,3) / 0.2e1;
t910 = mrSges(7,3) / 0.2e1;
t909 = -t75 / 0.2e1;
t906 = -t172 / 0.2e1;
t905 = -t173 / 0.2e1;
t903 = t273 / 0.2e1;
t901 = -t283 / 0.2e1;
t899 = -t284 / 0.2e1;
t890 = -t667 / 0.2e1;
t887 = -t368 / 0.2e1;
t886 = t371 / 0.2e1;
t874 = t499 / 0.2e1;
t869 = t949 / 0.2e1;
t868 = -t515 / 0.2e1;
t577 = pkin(4) * t590 + pkin(5);
t736 = t585 * t586;
t538 = -pkin(4) * t736 + t577 * t589;
t866 = -t538 / 0.2e1;
t732 = t586 * t589;
t539 = pkin(4) * t732 + t577 * t585;
t865 = -t539 / 0.2e1;
t545 = (-t585 * t590 - t732) * pkin(4);
t864 = -t545 / 0.2e1;
t546 = (t589 * t590 - t736) * pkin(4);
t863 = t546 / 0.2e1;
t862 = t553 / 0.2e1;
t859 = -t637 / 0.2e1;
t855 = -t585 / 0.2e1;
t854 = -t587 / 0.2e1;
t853 = t587 / 0.2e1;
t852 = -t591 / 0.2e1;
t851 = t591 / 0.2e1;
t849 = m(6) * t423;
t848 = (-t493 * t589 + t585 * t666) * t917;
t846 = pkin(3) * t565;
t845 = pkin(4) * t587;
t836 = t75 * mrSges(7,3);
t833 = t82 * mrSges(7,1);
t832 = t83 * mrSges(7,2);
t828 = t833 / 0.2e1 - t832 / 0.2e1;
t825 = Ifges(5,4) * t587;
t823 = Ifges(5,5) * t553;
t819 = Ifges(5,6) * t553;
t813 = pkin(4) * qJD(4);
t812 = pkin(5) * qJD(5);
t811 = t154 * mrSges(6,2);
t810 = t155 * mrSges(6,1);
t807 = t181 * mrSges(6,1);
t806 = t182 * mrSges(6,2);
t799 = t325 * mrSges(7,3);
t790 = t666 * mrSges(7,3);
t424 = -pkin(4) * t747 + t507;
t328 = -pkin(5) * t450 + t424;
t344 = mrSges(6,1) * t451 - mrSges(6,2) * t449;
t395 = -mrSges(6,2) * t555 + mrSges(6,3) * t450;
t397 = mrSges(6,1) * t555 - mrSges(6,3) * t452;
t486 = -mrSges(5,2) * t555 + mrSges(5,3) * t747;
t707 = mrSges(5,3) * t745;
t487 = -mrSges(5,2) * t553 - t707;
t488 = mrSges(5,1) * t555 + mrSges(5,3) * t746;
t489 = mrSges(5,1) * t553 - mrSges(5,3) * t744;
t544 = t555 * mrSges(4,1);
t647 = Ifges(7,4) * t333 + Ifges(7,2) * t330;
t648 = Ifges(6,4) * t452 + Ifges(6,2) * t450;
t649 = Ifges(7,1) * t333 + Ifges(7,4) * t330;
t650 = Ifges(6,1) * t452 + Ifges(6,4) * t450;
t651 = t797 - t798;
t652 = t792 - t794;
t5 = t649 * t891 - m(5) * (t357 * t384 + t358 * t385 + t506 * t507) - m(6) * (t154 * t174 + t155 * t175 + t423 * t424) + ((Ifges(5,1) * t926 + Ifges(4,1) - Ifges(4,2) - Ifges(5,3) - Ifges(6,3) - Ifges(7,3) + (t820 - 0.2e1 * t580) * t587) * t553 + Ifges(6,5) * t449 + Ifges(7,5) * t329 + Ifges(6,6) * t451 - Ifges(7,6) * t667 - t565 * t507 + t973 * t555) * t555 + (t690 * mrSges(4,2) - t973 * t553 + t636 + t945 + t953) * t553 - t358 * t486 - t385 * t487 - t357 * t488 - t384 * t489 - t424 * t344 - t175 * t396 - t154 * t397 - t174 * t398 - t155 * t395 - t328 * t185 - t75 * t274 - t76 * t272 - t327 * t651 - t423 * t652 + t451 * t648 / 0.2e1 - t690 * t544 - t641 * t273 - t642 * t275 - m(7) * (t327 * t328 + t641 * t76 + t642 * t75) + t650 * t882 + t647 * t890 + t452 * t899 + t450 * t901 + t333 * t905 + t330 * t906;
t786 = t5 * qJD(1);
t781 = t545 * mrSges(7,1);
t780 = t546 * mrSges(7,2);
t709 = pkin(4) * t744;
t394 = t709 - t838;
t602 = t155 * t796 + t172 * t891 + t329 * t997 + t423 * t343 + t661 * t862 + t667 * t993 + t990;
t663 = t905 + t836;
t753 = t506 * t555;
t6 = -m(6) * (t154 * t181 + t155 * t182) - m(7) * (t327 * t394 + t75 * t91 + t76 * t92) + t358 * t489 - t182 * t396 - t181 * t398 - (-t345 / 0.2e1 + t899 + t154 * mrSges(6,3)) * t451 - t394 * t185 + (t901 + t346 / 0.2e1) * t449 + t663 * t667 + t942 * t753 - t91 * t275 - t92 * t273 - t602 - t344 * t709 + ((-Ifges(5,4) * t745 + t823) * t587 + (t819 - pkin(4) * t849 + t358 * mrSges(5,3) + (t580 + (Ifges(5,1) - Ifges(5,2)) * t587) * t555) * t591) * t555 + (-t487 - t707) * t357;
t774 = t6 * qJD(1);
t7 = t283 * t882 - t449 * t346 / 0.2e1 - t155 * t398 + t83 * t273 + t82 * t275 + t602 - t185 * t838 + (t396 + t793) * t154 + m(7) * (-t327 * t838 + t75 * t82 + t76 * t83) - t951 * t451 / 0.2e1 + (t173 / 0.2e1 - t836) * t667;
t773 = t7 * qJD(1);
t497 = mrSges(6,1) * t637 + mrSges(6,2) * t638;
t594 = (t384 * t591 + t385 * t587) * t925 + (-t174 * t637 + t175 * t638) * t923 + (t493 * t641 + t642 * t666) * t921 + t274 * t879 + t272 * t876 + t397 * t859 + t395 * t857 + t486 * t853 + t488 * t851;
t598 = -m(5) * (-t841 * t940 - t847) / 0.2e1 + (t450 * t949 + t452 * t515 + t555 * t578) * t924 + (t1005 * t330 + t325 * t333 + t518 * t555) * t922;
t740 = t638 * t450;
t742 = t637 * t452;
t754 = t493 * t330;
t756 = t666 * t333;
t21 = (-mrSges(4,2) + (t926 / 0.2e1 + t927 / 0.2e1) * mrSges(5,3)) * t553 + (t888 - t497 / 0.2e1 + t826 / 0.2e1 - t776 / 0.2e1) * t555 + (t754 / 0.2e1 - t756 / 0.2e1) * mrSges(7,3) + (t740 / 0.2e1 + t742 / 0.2e1) * mrSges(6,3) + t594 + t544 + t598;
t772 = qJD(1) * t21;
t12 = t75 * t273 + t990 - t76 * t275 + t721 * t862 - (t906 - t997) * t329 + (t993 - t663) * t667;
t771 = t12 * qJD(1);
t725 = t591 * t487;
t731 = t587 * t489;
t22 = t333 * t273 + t330 * t275 + t452 * t396 + t450 * t398 + (mrSges(4,3) * t553 - t725 + t731) * t553 + (t185 + t344 + (mrSges(4,3) + t565) * t555) * t555 + m(7) * (t327 * t555 + t330 * t75 + t333 * t76) + m(6) * (t154 * t450 + t155 * t452 + t423 * t555) + m(5) * (t753 + (t357 * t587 - t358 * t591) * t553) + m(4) * (-t507 * t553 + t753) + (m(3) * qJ(2) + mrSges(3,3)) * (t583 ^ 2 + t584 ^ 2);
t768 = t22 * qJD(1);
t611 = (-t329 * t877 + t667 * t960) * mrSges(7,3) + t273 * t879 + t275 * t877;
t25 = t611 - t722;
t767 = t25 * qJD(1);
t765 = t329 * t493;
t764 = t667 * t666;
t759 = t449 * t638;
t758 = t451 * t637;
t751 = t538 * t667;
t750 = t538 * t666;
t749 = t539 * t329;
t748 = t539 * t493;
t738 = t585 * t329;
t737 = t585 * t493;
t735 = t586 * t395;
t734 = t586 * t449;
t730 = t589 * t667;
t729 = t589 * t666;
t728 = t590 * t397;
t727 = t590 * t451;
t252 = -mrSges(7,1) * t539 - t538 * mrSges(7,2);
t713 = qJD(6) * t252;
t711 = t918 / 0.2e1;
t708 = pkin(8) * t914;
t706 = t848 / 0.2e1;
t704 = -t842 / 0.2e1;
t701 = -Ifges(5,1) / 0.4e1 + Ifges(5,2) / 0.4e1;
t673 = -t154 / 0.2e1 + t182 / 0.2e1;
t672 = t181 / 0.2e1 + t155 / 0.2e1;
t665 = mrSges(7,3) * t705;
t664 = mrSges(7,3) * t704;
t659 = -t843 / 0.2e1 + t865;
t658 = t704 + t866;
t657 = (t345 / 0.4e1 + t284 / 0.4e1) * t637;
t569 = Ifges(5,1) * t591 - t825;
t646 = -t493 * t75 + t666 * t76;
t40 = t1002 + (t983 / 0.2e1 + t887) * t493 + (t886 + t964 / 0.2e1) * t666;
t595 = (mrSges(7,3) * t896 + t979) * t667 - (-t799 / 0.2e1 - t996) * t329 + t1005 * t903 + t275 * t893 + t981 * t861 + t995;
t8 = t595 - t654;
t645 = qJD(1) * t8 + qJD(3) * t40;
t599 = (-t1001 * t638 - t637 * t985) * t924 + (t493 * t92 + t666 * t91 + t646) * t922 + t731 / 0.2e1 - t725 / 0.2e1;
t600 = (t775 / 0.2e1 + t777 / 0.2e1) * t553 + (t330 * t538 + t333 * t539) * t921 + (t450 * t590 + t452 * t586) * t711 + t662;
t15 = (-t765 / 0.2e1 + t764 / 0.2e1) * mrSges(7,3) + (t758 / 0.2e1 - t759 / 0.2e1) * mrSges(6,3) + t599 + t600 + t929;
t644 = t15 * qJD(1);
t623 = t764 * t911 + t765 * t910 - t929 - (t758 - t759) * mrSges(6,3) / 0.2e1;
t604 = (t493 * t83 + t666 * t82 + t646) * t921 + t623;
t18 = t604 + t999;
t643 = t18 * qJD(1);
t639 = -t493 * t538 + t539 * t666;
t634 = t749 / 0.2e1 - t751 / 0.2e1;
t629 = -t980 - t496;
t526 = t844 + t845;
t592 = -(t900 + t889) * t638 + (t873 + t872) * t449 + (t877 * t91 + t879 * t92 + t615) * mrSges(7,3) + (-t637 * t673 - t638 * t672 + t946) * mrSges(6,3) - t952 + t953 / 0.2e1 + t526 * t961 + t515 * t883 + t394 * t365 / 0.2e1 - t657 + t396 * t869 + (t327 * t526 + t394 * t518 + (t76 + t91) * t1005 + (-t75 + t92) * t325) * t921 + t1012;
t596 = (t538 * t642 + t539 * t641) * t922 + t384 * t916 + t385 * t915 + t274 * t866 + t272 * t865 - (t174 * t590 + t175 * t586) * t918 / 0.2e1;
t566 = Ifges(5,2) * t591 + t825;
t605 = pkin(3) * t916 + t569 / 0.4e1 - t566 / 0.4e1 + (t708 - t701) * t591 + (m(6) * t856 + t497 / 0.2e1) * pkin(4);
t568 = Ifges(5,1) * t587 + t580;
t612 = pkin(3) * t915 - t568 / 0.4e1 + t943 / 0.4e1 + (t708 + t701) * t587;
t621 = -t1001 * t515 + t985 * t949;
t2 = (t423 * t845 + t621) * t923 + (-t728 / 0.2e1 + t344 * t853 - t735 / 0.2e1) * pkin(4) + t596 + (-t549 / 0.4e1 - t548 / 0.4e1 + t956 / 0.4e1 - t976 / 0.4e1 + t948) * t553 + t592 + (-Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t612 * t587 + (t605 - t825) * t591) * t555 + (t487 * t854 + t489 * t852) * pkin(8) + t603 + t931;
t608 = -t1002 + t983 * t877 - t578 * t496 + t988 + (t799 - t887) * t493 + (t964 + t371) * t960;
t34 = -t608 - t846 - t638 * t874 + t500 * t857 + t988 + (-t566 / 0.2e1 + t569 / 0.2e1 + pkin(4) * t497) * t587 + m(6) * t578 * t845 + (-t943 / 0.2e1 + t568 / 0.2e1) * t591 + mrSges(7,3) * t1014 + (t501 + t498) * t859 + t1022 * t526;
t626 = t2 * qJD(1) + t34 * qJD(3);
t3 = t398 * t868 - t951 * t637 / 0.4e1 - t952 + t946 * mrSges(6,3) + t966 - t1008;
t37 = -(t500 / 0.2e1 - t499 / 0.2e1 + t515 * mrSges(6,3) + t1022 * pkin(5)) * t638 - (-t501 / 0.2e1 - t498 / 0.2e1) * t637 + (t1005 * t666 - t1014) * mrSges(7,3) + t608 - t1005 * t790;
t625 = qJD(1) * t3 - qJD(3) * t37;
t607 = (-t329 * t865 + t667 * t866) * mrSges(7,3) + t538 * t903 + t275 * t865 + t671;
t13 = (t909 + t92 / 0.2e1) * mrSges(7,2) + (t908 - t91 / 0.2e1) * mrSges(7,1) + t607 - t671;
t620 = qJD(1) * t13 + qJD(4) * t252 + t1018;
t617 = m(7) * (t493 * t546 + t545 * t666 + t639);
t104 = t706 - t617 / 0.2e1;
t597 = (t586 * t883 + t590 * t884 + (t727 / 0.2e1 + t734 / 0.2e1) * mrSges(6,3)) * pkin(4) + (t538 * t82 + t539 * t83 + t545 * t75 + t546 * t76) * t921 + t545 * t902 + t273 * t863 + t828;
t11 = t673 * mrSges(6,2) - t672 * mrSges(6,1) + ((-t738 / 0.2e1 + t730 / 0.2e1) * pkin(5) + t634) * mrSges(7,3) + t597 + t969;
t616 = t781 + (-mrSges(6,1) * t586 - mrSges(6,2) * t590) * pkin(4) - t780;
t250 = m(7) * (t538 * t545 + t539 * t546) + t616;
t610 = ((t539 + t545) * t1005 + (-t538 + t546) * t325) * t921 + t1019;
t614 = -t748 / 0.2e1 - t750 / 0.2e1 + t666 * t863 + t493 * t864;
t39 = (t869 + t870) * mrSges(6,2) + (t515 / 0.2e1 + t868) * mrSges(6,1) + ((t737 / 0.2e1 + t729 / 0.2e1) * pkin(5) + t614) * mrSges(7,3) + t610 + t970;
t619 = t11 * qJD(1) - t104 * qJD(2) + t39 * qJD(3) + t250 * qJD(4);
t606 = (t589 * t903 + t275 * t855 + (-t329 * t855 + t589 * t890) * mrSges(7,3)) * pkin(5) + t671;
t19 = (t909 + t83 / 0.2e1) * mrSges(7,2) + (t908 - t82 / 0.2e1) * mrSges(7,1) + t606 - t671;
t223 = (t863 + t658) * mrSges(7,2) + (t864 + t659) * mrSges(7,1);
t562 = (t585 * mrSges(7,1) + t589 * mrSges(7,2)) * pkin(5);
t613 = -qJD(1) * t19 - qJD(4) * t223 + qJD(5) * t562 - t1018;
t556 = t562 * qJD(6);
t224 = -t780 / 0.2e1 + t781 / 0.2e1 + t658 * mrSges(7,2) + t659 * mrSges(7,1);
t95 = t617 / 0.2e1 + t706 + t629;
t38 = mrSges(7,3) * t614 - t493 * t665 + t664 * t666 + t1000 + t610 - t970;
t24 = t611 + t722;
t23 = t754 * t911 + t756 * t910 + t594 - t598 + (-t740 - t742) * mrSges(6,3) / 0.2e1 + t940 * t553 * t914 + (t942 + t497 + t365) * t860;
t20 = t606 + t828 + t932;
t17 = t604 - t999;
t16 = -t599 + t600 + t623;
t14 = t607 + t827 + t932;
t10 = t597 + t807 / 0.2e1 - t806 / 0.2e1 - t811 / 0.2e1 - t810 / 0.2e1 + t634 * mrSges(7,3) + t329 * t665 + t667 * t664 + t661 - t969;
t9 = t595 + t654;
t4 = -t657 + (t883 + t796 / 0.2e1) * t515 - (mrSges(6,3) * t870 + t669) * t451 + t966 + t1008;
t1 = (t489 * t919 + t823 / 0.4e1 + t605 * t555) * t591 + (-t819 / 0.2e1 + t487 * t919 + (t849 / 0.2e1 + t344 / 0.2e1) * pkin(4) + (t612 - t580) * t555) * t587 + Ifges(5,6) * t747 / 0.2e1 - Ifges(5,5) * t746 / 0.2e1 + t930 - t596 + (t728 + t735) * pkin(4) / 0.2e1 + t592 + Ifges(5,3) * t860 + t579 * t861 + t621 * t923 + t947;
t26 = [qJD(2) * t22 - qJD(3) * t5 - qJD(4) * t6 + qJD(5) * t7 + qJD(6) * t12, t768 + 0.2e1 * ((t330 * t666 + t333 * t493) * t921 + (-t450 * t637 + t452 * t638) * t923) * qJD(2) + t23 * qJD(3) + t16 * qJD(4) + t17 * qJD(5) + t24 * qJD(6), t23 * qJD(2) + t1 * qJD(4) + t4 * qJD(5) + t9 * qJD(6) - t786 + (-t175 * t637 * mrSges(6,3) + t641 * t790 - t174 * t779 + t515 * t395 + t506 * mrSges(4,2) - t507 * mrSges(4,1) + t424 * t497 + t452 * t501 / 0.2e1 + t328 * t365 + t330 * t368 / 0.2e1 + t325 * t272 + (Ifges(5,5) * t587 + Ifges(6,5) * t638 + Ifges(7,5) * t493 + Ifges(5,6) * t591 - Ifges(6,6) * t637 + Ifges(7,6) * t666 - Ifges(4,6)) * t555 + t949 * t397 + 0.2e1 * (t174 * t949 + t175 * t515 + t424 * t578) * t923 + t518 * t651 + t578 * t652 + t591 * pkin(8) * t486 - t587 * pkin(8) * t488 + (t566 * t853 + t568 * t852 + t569 * t854 + t851 * t943 - Ifges(4,5) + t846) * t553 + t944 * mrSges(5,3) + 0.2e1 * (-pkin(3) * t507 + pkin(8) * t944) * t925 - t642 * t493 * mrSges(7,3) + t507 * t942 + t1005 * t274 + 0.2e1 * (t1005 * t642 + t325 * t641 + t518 * t328) * t921 + t650 * t857 + t648 * t859 + t450 * t874 + t649 * t876 + t647 * t879 + t333 * t886) * qJD(3), -t774 + t16 * qJD(2) + t1 * qJD(3) + (-Ifges(5,5) * t745 - Ifges(5,6) * t744 + m(7) * (t538 * t91 + t539 * t92) - t830 + t831 + t807 - t806 - t357 * mrSges(5,2) - t358 * mrSges(5,1) + t661 + (t749 - t751) * mrSges(7,3)) * qJD(4) + t10 * qJD(5) + t14 * qJD(6) + (m(6) * (t181 * t590 + t182 * t586) + (t727 + t734) * mrSges(6,3)) * t813, t773 + t17 * qJD(2) + t4 * qJD(3) + t10 * qJD(4) + (t661 - t810 - t811 - t832 + t833) * qJD(5) + t20 * qJD(6) + (m(7) * (t585 * t83 + t589 * t82) + (-t730 + t738) * mrSges(7,3)) * t812, t771 + t24 * qJD(2) + t9 * qJD(3) + t14 * qJD(4) + t20 * qJD(5) + (t721 - t835 - t837) * qJD(6); qJD(3) * t21 - qJD(4) * t15 + qJD(5) * t18 + qJD(6) * t25 - t768, 0, t772 (-t565 + t629) * qJD(4) + t95 * qJD(5) + 0.2e1 * (t639 * t921 + (-t586 * t637 - t590 * t638) * t711) * qJD(4) - t644 + t933, t95 * qJD(4) + (t629 + t848) * qJD(5) + t643 + t933, -qJD(6) * t980 + t624 * t971 + t767; -qJD(2) * t21 + qJD(4) * t2 + qJD(5) * t3 + qJD(6) * t8 + t786, -t772, qJD(4) * t34 - qJD(5) * t37 + qJD(6) * t40 (m(7) * (t1005 * t539 - t325 * t538) + t942 * pkin(8) + (-t748 - t750) * mrSges(7,3) + t948 + t1020) * qJD(4) + t38 * qJD(5) + t1024 + (m(6) * (-t515 * t590 + t586 * t949) + (-t586 * t638 + t590 * t637) * mrSges(6,3)) * t813 + t626, t38 * qJD(4) + t1020 * qJD(5) + t1024 + (m(7) * t1013 + (-t729 - t737) * mrSges(7,3)) * t812 + t625 (t981 + t1006) * qJD(6) + t645 + t971 * t1023; qJD(2) * t15 - qJD(3) * t2 + qJD(5) * t11 + qJD(6) * t13 + t774, -qJD(5) * t104 + t644 + t934, qJD(5) * t39 + t1021 - t626, qJD(5) * t250 + t713 ((t545 * t589 + t546 * t585) * t917 + t616) * qJD(5) + t224 * qJD(6) + t619, qJD(5) * t224 + t620 + t713; -qJD(2) * t18 - qJD(3) * t3 - qJD(4) * t11 + qJD(6) * t19 - t773, qJD(4) * t104 - t643 + t934, -qJD(4) * t39 + t1021 - t625, qJD(6) * t223 - t619, -t556, -t556 - t613; -qJD(2) * t25 - qJD(3) * t8 - qJD(4) * t13 - qJD(5) * t19 - t771, -t622 * t971 - t767, -t1011 * t971 - t645, -qJD(5) * t223 - t620, t613, 0;];
Cq  = t26;
