% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRR11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR11_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR11_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:29:22
% EndTime: 2019-03-09 14:30:04
% DurationCPUTime: 26.04s
% Computational Cost: add. (44279->835), mult. (82361->1096), div. (0->0), fcn. (86637->8), ass. (0->448)
t947 = qJD(4) + qJD(5);
t568 = sin(qJ(4));
t573 = -pkin(2) - pkin(8);
t721 = t568 * t573;
t529 = -t568 * pkin(9) + t721;
t571 = cos(qJ(4));
t530 = (-pkin(9) + t573) * t571;
t567 = sin(qJ(5));
t845 = cos(qJ(5));
t437 = t845 * t529 + t567 * t530;
t622 = t567 * t571 + t568 * t845;
t356 = -t622 * pkin(10) + t437;
t566 = sin(qJ(6));
t570 = cos(qJ(6));
t677 = t845 * t571;
t517 = t567 * t568 - t677;
t922 = -t567 * t529 + t845 * t530;
t948 = t517 * pkin(10) + t922;
t1008 = -t356 * t566 + t570 * t948;
t1012 = t1008 * mrSges(7,2);
t208 = t356 * t570 + t566 * t948;
t1016 = t208 * mrSges(7,1);
t1019 = -t1016 / 0.2e1 - t1012 / 0.2e1;
t410 = t517 * t570 + t566 * t622;
t657 = t517 * t566 - t570 * t622;
t971 = Ifges(7,5) * t657 + Ifges(7,6) * t410;
t1024 = 0.2e1 * t1019 + t971;
t1025 = t1024 * qJD(6);
t889 = -t208 / 0.2e1;
t955 = t410 * mrSges(7,2);
t990 = t657 * mrSges(7,1);
t1001 = t990 + t955;
t552 = t568 * pkin(4) + qJ(3);
t463 = pkin(5) * t622 + t552;
t843 = m(7) * t463;
t1023 = -t1001 + t843;
t569 = sin(qJ(2));
t572 = cos(qJ(2));
t482 = t517 * t572;
t483 = t622 * t572;
t351 = t482 * t566 - t483 * t570;
t800 = mrSges(7,3) * t351;
t301 = mrSges(7,1) * t569 - t800;
t1004 = t301 * t889;
t780 = t483 * mrSges(6,3);
t446 = mrSges(6,1) * t569 + t780;
t1022 = t1004 - t437 * t446 / 0.2e1;
t993 = t1008 / 0.2e1;
t1013 = (-t1008 / 0.2e1 + t993) * mrSges(7,2) + (t889 + t208 / 0.2e1) * mrSges(7,1);
t1021 = qJD(6) * t1013;
t650 = -Ifges(6,5) * t622 + Ifges(6,6) * t517 + t971;
t1002 = -t437 * mrSges(6,1) - t922 * mrSges(6,2) + t650;
t1009 = -t1016 - t1012;
t1020 = t1002 + t1009;
t681 = t990 / 0.2e1;
t964 = t410 / 0.2e1;
t965 = -t410 / 0.2e1;
t603 = t681 - t990 / 0.2e1 + (t964 + t965) * mrSges(7,2);
t1018 = qJD(2) * t1013 + qJD(3) * t603;
t1015 = t208 * t410;
t1014 = t1008 * t566 - t208 * t570;
t658 = t570 * t482 + t483 * t566;
t796 = t351 * Ifges(7,4);
t201 = Ifges(7,2) * t658 + Ifges(7,6) * t569 + t796;
t332 = Ifges(7,4) * t658;
t203 = Ifges(7,1) * t351 + t569 * Ifges(7,5) + t332;
t891 = pkin(3) + pkin(7);
t537 = t891 * t572;
t713 = t571 * t572;
t499 = pkin(4) * t713 + t537;
t396 = -t482 * pkin(5) + t499;
t402 = Ifges(7,4) * t657;
t260 = -Ifges(7,1) * t410 + t402;
t925 = Ifges(7,2) * t410 + t260 + t402;
t940 = -Ifges(7,2) * t351 + t332;
t942 = -mrSges(7,1) * t410 + mrSges(7,2) * t657;
t943 = mrSges(7,1) * t351 + mrSges(7,2) * t658;
t974 = Ifges(7,1) * t658 - t796;
t1007 = (-t974 / 0.4e1 + t201 / 0.4e1) * t410 + t942 * t396 / 0.2e1 + t943 * t463 / 0.2e1 + (t203 + t940) * t657 / 0.4e1 + t925 * t658 / 0.4e1;
t821 = Ifges(7,4) * t410;
t258 = Ifges(7,2) * t657 - t821;
t973 = Ifges(7,1) * t657 + t821;
t1000 = t973 / 0.4e1 - t258 / 0.4e1;
t798 = t658 * mrSges(7,3);
t299 = -mrSges(7,2) * t569 + t798;
t994 = t657 / 0.2e1;
t998 = t299 * t965 + t301 * t994;
t995 = -t657 / 0.2e1;
t33 = t258 * t965 - t463 * t942 + t925 * t995 + t964 * t973;
t997 = -mrSges(6,2) / 0.2e1;
t996 = -Ifges(7,3) / 0.2e1;
t992 = m(6) * t552;
t954 = t410 * mrSges(7,3);
t987 = t396 * t943;
t775 = t517 * mrSges(6,3);
t986 = t437 * t775;
t769 = qJ(3) * t569;
t502 = t572 * t573 - pkin(1) - t769;
t536 = t891 * t569;
t420 = -t502 * t568 + t571 * t536;
t722 = t568 * t572;
t380 = pkin(9) * t722 + t420;
t358 = pkin(4) * t569 + t380;
t421 = t502 * t571 + t536 * t568;
t381 = -pkin(9) * t713 + t421;
t362 = t845 * t381;
t226 = t567 * t358 + t362;
t237 = -t380 * t567 - t362;
t980 = t226 + t237;
t468 = Ifges(6,4) * t482;
t339 = -Ifges(6,1) * t483 + Ifges(6,5) * t569 + t468;
t367 = Ifges(6,2) * t483 + t468;
t979 = t367 + t339;
t508 = Ifges(6,4) * t622;
t431 = Ifges(6,2) * t517 - t508;
t434 = -Ifges(6,1) * t517 - t508;
t978 = t434 + t431;
t715 = t570 * t657;
t729 = t566 * t410;
t976 = t715 - t729;
t696 = t845 * pkin(4);
t556 = t696 + pkin(5);
t725 = t567 * t570;
t497 = pkin(4) * t725 + t556 * t566;
t739 = t497 * t410;
t728 = t566 * t567;
t496 = -pkin(4) * t728 + t556 * t570;
t741 = t496 * t657;
t975 = t739 - t741;
t723 = t568 * t569;
t481 = -t567 * t723 + t569 * t677;
t484 = t622 * t569;
t455 = t570 * t484;
t350 = t481 * t566 + t455;
t881 = t350 / 0.2e1;
t347 = t481 * t570 - t484 * t566;
t885 = t347 / 0.2e1;
t708 = mrSges(7,1) * t881 + mrSges(7,2) * t885;
t972 = -t484 * mrSges(6,1) / 0.2e1 + t481 * t997 - t708;
t626 = Ifges(7,5) * t881 + Ifges(7,6) * t885;
t659 = t569 * pkin(2) - qJ(3) * t572;
t515 = pkin(8) * t569 + t659;
t521 = t571 * t537;
t359 = pkin(4) * t572 + t521 + (-pkin(9) * t569 - t515) * t568;
t427 = t571 * t515 + t568 * t537;
t718 = t569 * t571;
t383 = pkin(9) * t718 + t427;
t229 = t845 * t359 - t383 * t567;
t164 = pkin(5) * t572 - pkin(10) * t484 + t229;
t230 = t567 * t359 + t845 * t383;
t173 = pkin(10) * t481 + t230;
t93 = t164 * t570 - t173 * t566;
t94 = t164 * t566 + t173 * t570;
t920 = t94 * mrSges(7,2) / 0.2e1 - t93 * mrSges(7,1) / 0.2e1;
t646 = -t572 * t996 + t626 - t920;
t928 = Ifges(7,5) * t658;
t957 = Ifges(7,6) * t351;
t707 = t928 - t957;
t665 = t957 / 0.2e1 - t928 / 0.2e1;
t961 = t569 / 0.4e1;
t927 = t657 * mrSges(7,3);
t682 = -t927 / 0.2e1;
t950 = t350 * t410;
t503 = (-t566 * t845 - t725) * pkin(4);
t504 = (t570 * t845 - t728) * pkin(4);
t360 = t567 * t381;
t225 = t845 * t358 - t360;
t472 = t483 * pkin(10);
t170 = t225 + t472;
t837 = t482 * pkin(10);
t171 = t226 + t837;
t767 = t171 * t566;
t96 = t170 * t570 - t767;
t830 = t96 * mrSges(7,2);
t766 = t171 * t570;
t95 = -t170 * t566 - t766;
t831 = t95 * mrSges(7,1);
t828 = t831 / 0.2e1 - t830 / 0.2e1;
t860 = -t503 / 0.2e1;
t162 = pkin(5) * t569 + t170;
t87 = t162 * t570 - t767;
t88 = t162 * t566 + t766;
t888 = -t299 / 0.2e1;
t899 = -m(7) / 0.2e1;
t949 = (t496 * t95 + t497 * t96 + t503 * t87 + t504 * t88) * t899 + t301 * t860 + t504 * t888 - t828;
t895 = m(7) * pkin(5);
t946 = -t1014 * t895 / 0.2e1 - t1019;
t944 = t230 * t997 + t229 * mrSges(6,1) / 0.2e1;
t781 = t483 * mrSges(6,2);
t783 = t482 * mrSges(6,1);
t797 = t351 * mrSges(7,2);
t799 = t658 * mrSges(7,1);
t854 = -t537 / 0.2e1;
t901 = -m(6) / 0.2e1;
t941 = -(t950 / 0.2e1 + t347 * t995) * mrSges(7,3) + m(5) * t854 + (-t437 * t481 + t484 * t922 + t499) * t901 + (t1008 * t350 - t208 * t347 + t396) * t899 + t799 / 0.2e1 - t797 / 0.2e1 + t783 / 0.2e1 + t781 / 0.2e1;
t185 = t237 - t837;
t238 = t845 * t380 - t360;
t186 = t472 + t238;
t101 = t185 * t570 - t186 * t566;
t937 = -t101 / 0.2e1;
t936 = -t201 / 0.2e1;
t883 = t658 / 0.2e1;
t934 = -Ifges(3,4) - Ifges(4,6);
t676 = t658 * t964;
t864 = t484 / 0.2e1;
t867 = t481 / 0.2e1;
t923 = Ifges(6,5) * t864 + Ifges(6,6) * t867;
t829 = -t88 * t410 + t657 * t87;
t898 = m(7) / 0.2e1;
t921 = (-t410 * t95 - t657 * t96 + t829) * t898;
t466 = Ifges(6,6) * t483;
t467 = Ifges(6,5) * t482;
t651 = t467 + t466 + t707;
t919 = -t229 * t517 + t230 * t622;
t827 = mrSges(5,3) * t572;
t903 = t571 ^ 2;
t904 = t568 ^ 2;
t914 = (t903 + t904) * t827;
t782 = t482 * mrSges(6,3);
t444 = -mrSges(6,2) * t569 + t782;
t679 = t775 / 0.2e1;
t683 = t954 / 0.2e1;
t684 = -t798 / 0.2e1;
t685 = -t800 / 0.2e1;
t823 = Ifges(6,4) * t483;
t337 = Ifges(6,2) * t482 + t569 * Ifges(6,6) - t823;
t364 = -mrSges(6,1) * t483 + mrSges(6,2) * t482;
t368 = Ifges(6,1) * t482 + t823;
t428 = -mrSges(6,1) * t517 - mrSges(6,2) * t622;
t822 = Ifges(6,4) * t517;
t432 = -Ifges(6,2) * t622 - t822;
t433 = -Ifges(6,1) * t622 + t822;
t905 = t437 * t780 / 0.2e1 - t922 * t782 / 0.2e1 + t552 * t364 / 0.2e1 + t499 * t428 / 0.2e1 + t978 * t482 / 0.4e1 + t650 * t961 - t979 * t622 / 0.4e1 + (-t368 / 0.4e1 + t337 / 0.4e1) * t517 + (-t433 / 0.4e1 + t432 / 0.4e1) * t483 + t1000 * t351 + t1007;
t583 = t1008 * t684 + t208 * t685 + t226 * t679 + t87 * t682 + t88 * t683 + t905;
t680 = -t775 / 0.2e1;
t653 = t226 * t680;
t887 = t299 / 0.2e1;
t913 = (t95 * t964 + t96 * t994) * mrSges(7,3) + t922 * t444 / 0.2e1 + t653 + t583 + ((t88 + t95) * t898 + t887) * t1008 + (-t87 + t96) * t898 * t208 + t1022;
t912 = qJD(6) * t603;
t608 = t955 + 0.2e1 * t681;
t911 = t608 * qJD(6);
t834 = t88 * mrSges(7,1);
t835 = t87 * mrSges(7,2);
t910 = -t834 / 0.2e1 - t835 / 0.2e1 - t665;
t102 = t185 * t566 + t186 * t570;
t187 = t622 * t225;
t692 = mrSges(5,3) * t713;
t525 = -t569 * mrSges(5,2) - t692;
t714 = t571 * t525;
t693 = mrSges(5,3) * t722;
t523 = mrSges(5,1) * t569 + t693;
t724 = t568 * t523;
t900 = m(6) / 0.2e1;
t909 = (t238 * t622 - t517 * t980 - t187) * t900 + (-t101 * t410 - t102 * t657 + t829) * t898 + t714 / 0.2e1 - t724 / 0.2e1;
t213 = t797 - t799;
t858 = -t517 / 0.2e1;
t865 = -t483 / 0.2e1;
t908 = (-t396 * t517 - t463 * t483) * t898 + t213 * t858 - t1001 * t865;
t846 = t572 / 0.2e1;
t907 = Ifges(6,3) * t846 + t646 + t923 + t944;
t855 = -t622 / 0.2e1;
t856 = t622 / 0.2e1;
t857 = t517 / 0.2e1;
t906 = t444 * t858 + t446 * t855 + (t482 * t857 + t483 * t856) * mrSges(6,3) + (-t351 * t995 + t676) * mrSges(7,3) + t998;
t902 = m(5) / 0.2e1;
t896 = m(6) * pkin(4);
t894 = mrSges(5,1) / 0.2e1;
t893 = -t87 / 0.2e1;
t892 = -t88 / 0.2e1;
t884 = -t658 / 0.2e1;
t880 = t351 / 0.2e1;
t869 = -t444 / 0.2e1;
t866 = t482 / 0.2e1;
t863 = -t496 / 0.2e1;
t862 = t496 / 0.2e1;
t861 = -t497 / 0.2e1;
t859 = t504 / 0.2e1;
t853 = -t566 / 0.2e1;
t852 = t568 / 0.2e1;
t851 = t568 / 0.4e1;
t850 = t569 / 0.2e1;
t849 = -t571 / 0.2e1;
t848 = t571 / 0.2e1;
t847 = t571 / 0.4e1;
t844 = t976 * t895;
t842 = pkin(4) * t567;
t841 = pkin(4) * t571;
t840 = pkin(5) * t483;
t839 = pkin(5) * t566;
t838 = pkin(5) * t570;
t826 = Ifges(5,1) * t568;
t825 = Ifges(5,4) * t568;
t824 = Ifges(5,4) * t571;
t820 = Ifges(5,5) * t568;
t819 = Ifges(5,5) * t571;
t818 = Ifges(5,6) * t568;
t817 = Ifges(5,6) * t571;
t815 = pkin(5) * qJD(5);
t814 = t101 * mrSges(7,1);
t813 = t102 * mrSges(7,2);
t806 = t225 * mrSges(6,2);
t805 = t226 * mrSges(6,1);
t802 = t237 * mrSges(6,1);
t801 = t238 * mrSges(6,2);
t200 = Ifges(7,4) * t350 + Ifges(7,2) * t347 + Ifges(7,6) * t572;
t202 = Ifges(7,1) * t350 + Ifges(7,4) * t347 + Ifges(7,5) * t572;
t212 = -mrSges(7,1) * t347 + mrSges(7,2) * t350;
t298 = -mrSges(7,2) * t572 + mrSges(7,3) * t347;
t300 = mrSges(7,1) * t572 - mrSges(7,3) * t350;
t336 = Ifges(6,4) * t484 + Ifges(6,2) * t481 + Ifges(6,6) * t572;
t338 = Ifges(6,1) * t484 + Ifges(6,4) * t481 + Ifges(6,5) * t572;
t365 = -mrSges(6,1) * t481 + mrSges(6,2) * t484;
t366 = -t781 - t783;
t498 = (-t841 - t891) * t569;
t395 = -pkin(5) * t481 + t498;
t426 = -t515 * t568 + t521;
t443 = -mrSges(6,2) * t572 + t481 * mrSges(6,3);
t445 = mrSges(6,1) * t572 - mrSges(6,3) * t484;
t640 = Ifges(5,2) * t571 + t825;
t477 = Ifges(5,6) * t572 + t569 * t640;
t642 = t824 + t826;
t478 = Ifges(5,5) * t572 + t569 * t642;
t645 = t571 * mrSges(5,1) - t568 * mrSges(5,2);
t501 = t645 * t569;
t774 = t572 * mrSges(5,1);
t522 = -mrSges(5,3) * t723 + t774;
t773 = t572 * mrSges(5,2);
t524 = mrSges(5,3) * t718 - t773;
t638 = -pkin(2) * t572 - t769;
t531 = -pkin(1) + t638;
t596 = t626 + t923;
t639 = -t817 - t820;
t532 = t572 * mrSges(4,2) - t569 * mrSges(4,3);
t662 = m(4) * t531 + t532;
t5 = m(5) * (t420 * t426 + t421 * t427 - t536 * t537) + m(6) * (t225 * t229 + t226 * t230 + t498 * t499) + m(7) * (t395 * t396 + t87 * t93 + t88 * t94) + t662 * t659 - t537 * t501 + t420 * t522 + t426 * t523 + t421 * t524 + t427 * t525 + t498 * t366 + t499 * t365 + t226 * t443 + t230 * t444 + t225 * t445 + t229 * t446 + t395 * t213 + t396 * t212 + t94 * t299 + t87 * t300 + t93 * t301 + t88 * t298 + t203 * t881 + t200 * t883 + t201 * t885 + t337 * t867 + t202 * t880 + t339 * t864 + t338 * t865 + t336 * t866 + (Ifges(7,5) * t880 + Ifges(7,6) * t883 + Ifges(6,5) * t865 + Ifges(6,6) * t866 - t536 * t645 + t477 * t849 - t568 * t478 / 0.2e1 - pkin(1) * mrSges(3,2) - t531 * mrSges(4,3) + (-t820 / 0.2e1 - t817 / 0.2e1 - t934) * t572) * t572 + (-pkin(1) * mrSges(3,1) - t531 * mrSges(4,2) + (-t639 + t934) * t569 + (Ifges(5,3) + Ifges(4,2) - Ifges(4,3) + Ifges(3,1) - Ifges(3,2) + Ifges(6,3) + Ifges(7,3) - t903 * Ifges(5,2) / 0.2e1 + (-t824 - t826 / 0.2e1) * t568) * t572 + t596) * t569;
t778 = t5 * qJD(1);
t777 = t503 * mrSges(7,1);
t776 = t504 * mrSges(7,2);
t591 = t226 * t780 + t351 * t936 + t499 * t364 - t88 * t800 + t880 * t974 + t883 * t940 + t987;
t644 = -t568 * mrSges(5,1) - t571 * mrSges(5,2);
t620 = t644 * t572;
t697 = pkin(4) * t722;
t627 = -t697 - t840;
t655 = -t203 / 0.2e1 + t87 * mrSges(7,3);
t6 = -t627 * t213 + (-t467 / 0.2e1 - t466 / 0.2e1 + t572 * (-t818 + t819) + t665) * t569 - t238 * t444 - t237 * t446 + (-t337 / 0.2e1 + t368 / 0.2e1) * t483 - t102 * t299 - t101 * t301 - t591 + (Ifges(5,4) * t722 + pkin(4) * t366) * t722 - m(6) * (t225 * t237 + t226 * t238 - t499 * t697) - m(7) * (t87 * t101 + t88 * t102 + t396 * t627) + (-t824 + (-Ifges(5,1) + Ifges(5,2)) * t568) * t572 ^ 2 * t571 + (-t367 / 0.2e1 - t339 / 0.2e1 + t225 * mrSges(6,3)) * t482 + t655 * t658 - t537 * t620 + (t523 - t693) * t421 + (-t525 - t692) * t420;
t772 = t6 * qJD(1);
t9 = m(7) * (-t396 * t840 + t87 * t95 + t88 * t96) + (t444 - t782) * t225 + t483 * t337 / 0.2e1 + t368 * t865 - t226 * t446 + t203 * t883 + t96 * t299 + t95 * t301 + t591 - t213 * t840 - t87 * t798 + t979 * t866 + t651 * t850;
t771 = t9 * qJD(1);
t770 = -t813 / 0.2e1 + t814 / 0.2e1;
t12 = t707 * t850 + t987 - t88 * t301 + t87 * t299 + (t936 + t974 / 0.2e1 - t88 * mrSges(7,3)) * t351 + (t940 / 0.2e1 - t655) * t658;
t768 = t12 * qJD(1);
t589 = (t351 * t994 + t676) * mrSges(7,3) + t998;
t28 = t589 - t708;
t762 = t28 * qJD(1);
t32 = t350 * t301 - t347 * t299 + t484 * t446 - t481 * t444 + m(7) * (-t347 * t88 + t350 * t87) + m(6) * (t225 * t484 - t226 * t481) + (-t714 + t724 + m(5) * (t420 * t568 - t421 * t571) - t662) * t569;
t761 = t32 * qJD(1);
t748 = t481 * t622;
t743 = t484 * t517;
t742 = t496 * t658;
t740 = t497 * t351;
t731 = t566 * t298;
t730 = t566 * t351;
t727 = t567 * t446;
t717 = t570 * t300;
t716 = t570 * t658;
t712 = t571 * t573;
t163 = -t497 * mrSges(7,1) - t496 * mrSges(7,2);
t702 = qJD(6) * t163;
t700 = mrSges(6,3) * t842;
t699 = t896 / 0.2e1;
t691 = t844 / 0.2e1;
t678 = mrSges(6,3) * t856;
t669 = t723 / 0.2e1;
t668 = -t722 / 0.2e1;
t667 = t718 / 0.2e1;
t429 = mrSges(6,1) * t622 - t517 * mrSges(6,2);
t660 = t483 * t700;
t656 = mrSges(6,3) * t696;
t654 = t696 / 0.2e1;
t649 = -t839 / 0.2e1 + t861;
t648 = -t838 / 0.2e1 + t863;
t647 = t482 * t656;
t535 = Ifges(5,1) * t571 - t825;
t534 = -t568 * Ifges(5,2) + t824;
t575 = t1004 + t971 * t961 + (mrSges(7,3) * t884 + t887) * t1008 + (t889 * mrSges(7,3) + t1000) * t351 + t1007;
t8 = t575 - t646;
t636 = t8 * qJD(1) - t33 * qJD(2);
t578 = t906 + t972;
t582 = (-t347 * t497 + t350 * t496) * t898 + (-t481 * t567 + t484 * t845) * t699 + mrSges(5,1) * t669 + mrSges(5,2) * t667;
t14 = t578 - t582 + t914 / 0.2e1 + t909;
t634 = t14 * qJD(1);
t613 = (t484 * t566 ^ 2 + t455 * t570) * t895;
t19 = -t613 / 0.2e1 + t578 + t921;
t633 = t19 * qJD(1);
t629 = t426 * t571 + t427 * t568;
t615 = t1001 - t429;
t597 = t615 + t644;
t115 = mrSges(4,3) + (m(5) + m(4)) * qJ(3) + t992 + t843 - t597;
t580 = t629 * t902 + t919 * t900 + (-t410 * t93 - t657 * t94) * t898 + t300 * t965 + t298 * t995 + t445 * t858 + t443 * t856;
t22 = (t522 / 0.2e1 - t774 / 0.2e1) * t571 + (t524 / 0.2e1 + t773 / 0.2e1) * t568 + (-t743 / 0.2e1 - t748 / 0.2e1) * mrSges(6,3) + t580 + t941;
t628 = -qJD(1) * t22 + qJD(2) * t115;
t625 = m(7) * (t566 * t94 + t570 * t93);
t624 = -t730 / 0.2e1 - t716 / 0.2e1;
t623 = t534 * t848 + t535 * t852;
t612 = (t101 * t570 + t102 * t566) * t898;
t479 = -pkin(5) * t517 + t841;
t574 = -mrSges(6,3) * t187 / 0.2e1 + t523 * t721 / 0.2e1 - t525 * t712 / 0.2e1 + t535 * t713 / 0.2e1 + t429 * t697 / 0.2e1 + t102 * t682 - t1022 - t366 * t841 / 0.2e1 + t922 * t869 - 0.2e1 * (t640 * t847 + t642 * t851) * t572 + ((t499 * t571 - t552 * t722) * pkin(4) + t980 * t922 + (-t225 + t238) * t437) * t901 + (Ifges(5,6) * t847 + Ifges(5,5) * t851 - t639 / 0.4e1) * t569 - t479 * t213 / 0.2e1 + t645 * t854 + t627 * t1001 / 0.2e1 + t954 * t937 - qJ(3) * t620 / 0.2e1 + ((t101 + t88) * t899 + t888) * t1008 + t534 * t668 + t238 * t678 + t237 * t680 - t573 * t914 / 0.2e1 + (t479 * t396 + t463 * t627 + (t102 - t87) * t208) * t899;
t576 = (t496 * t93 + t497 * t94) * t898 + Ifges(5,3) * t846 + t300 * t862 + t497 * t298 / 0.2e1 + t426 * t894 - t427 * mrSges(5,2) / 0.2e1 + (t229 * t845 + t230 * t567) * t699 + Ifges(5,6) * t667 + Ifges(5,5) * t669 + t443 * t842 / 0.2e1 + t445 * t654 + t907;
t1 = t574 + t798 * t993 + t653 + t576 - t800 * t889 + t954 * t892 - t905 + t87 * t927 / 0.2e1;
t587 = -t208 * t954 - t552 * t428 + t33 - t986;
t17 = -t587 + qJ(3) * t645 + t640 * t852 + t642 * t849 + t432 * t857 + t433 * t858 - t986 - t1015 * mrSges(7,3) - t623 + t978 * t855 + (t992 + t429) * t841 + t1023 * t479;
t611 = -t1 * qJD(1) + t17 * qJD(2);
t20 = (t437 * mrSges(6,3) - t432 / 0.2e1 + t433 / 0.2e1 + t1023 * pkin(5)) * t517 - (-t431 / 0.2e1 - t434 / 0.2e1) * t622 + (t1008 * t657 + t1015) * mrSges(7,3) + t587 - t1008 * t927;
t4 = (-t717 / 0.2e1 - t731 / 0.2e1 - t625 / 0.2e1 + t908) * pkin(5) + (-Ifges(6,3) / 0.2e1 + t996) * t572 - t596 + t913 + t920 - t944;
t610 = t4 * qJD(1) - t20 * qJD(2);
t11 = (-t238 / 0.2e1 + t225 / 0.2e1) * mrSges(6,2) + (t237 / 0.2e1 + t226 / 0.2e1) * mrSges(6,1) + pkin(5) * t612 + (t727 / 0.2e1 + t845 * t869 + (t567 * t865 + t845 * t866) * mrSges(6,3)) * pkin(4) + (t740 / 0.2e1 + t742 / 0.2e1 + t624 * pkin(5)) * mrSges(7,3) + t770 + t949;
t593 = t777 + (-t567 * mrSges(6,1) - mrSges(6,2) * t845) * pkin(4) - t776;
t160 = m(7) * (t496 * t503 + t497 * t504) + t593;
t588 = ((t497 + t503) * t1008 + (-t496 + t504) * t208) * t898 + t1019;
t594 = t739 / 0.2e1 - t741 / 0.2e1 + t657 * t859 - t410 * t860;
t26 = ((-t729 / 0.2e1 + t715 / 0.2e1) * pkin(5) + t594) * mrSges(7,3) + t588 + t946;
t598 = m(7) * (-t410 * t503 - t504 * t657 - t975);
t70 = t691 - t598 / 0.2e1;
t601 = -t11 * qJD(1) + t26 * qJD(2) - t70 * qJD(3) + t160 * qJD(4);
t586 = (t351 * t861 + t658 * t863) * mrSges(7,3) + t299 * t862 + t301 * t861 - t665;
t16 = (t893 + t102 / 0.2e1) * mrSges(7,2) + (t892 + t937) * mrSges(7,1) + t586 + t665;
t600 = t16 * qJD(1) + t163 * qJD(4) + t1018;
t131 = (t859 + t648) * mrSges(7,2) + (t860 + t649) * mrSges(7,1);
t584 = (t570 * t887 + t301 * t853 + (t351 * t853 + t570 * t884) * mrSges(7,3)) * pkin(5) - t665;
t24 = (t893 + t96 / 0.2e1) * mrSges(7,2) + (t892 - t95 / 0.2e1) * mrSges(7,1) + t584 + t665;
t526 = (t566 * mrSges(7,1) + t570 * mrSges(7,2)) * pkin(5);
t592 = -qJD(1) * t24 - qJD(4) * t131 + qJD(5) * t526 - t1018;
t581 = t906 - t972;
t516 = t526 * qJD(6);
t132 = -t776 / 0.2e1 + t777 / 0.2e1 + t648 * mrSges(7,2) + t649 * mrSges(7,1);
t62 = t598 / 0.2e1 + t691 + t615;
t27 = t589 + t708;
t25 = t594 * mrSges(7,3) + t682 * t838 + t683 * t839 + t1002 + t588 - t946;
t23 = t584 + t828 + t910;
t21 = t522 * t848 + t524 * t852 + t580 + (m(4) * pkin(7) + mrSges(4,1)) * t572 + t713 * t894 + mrSges(5,2) * t668 + t481 * t678 + t484 * t679 - t941;
t18 = t613 / 0.2e1 + t581 + t921;
t15 = t586 + t770 + t910;
t13 = t581 + (t904 / 0.2e1 + t903 / 0.2e1) * t827 + t582 + t909;
t10 = -t647 / 0.2e1 + t802 / 0.2e1 - t801 / 0.2e1 - t805 / 0.2e1 - t806 / 0.2e1 + (mrSges(7,3) * t624 + t612) * pkin(5) - pkin(4) * t727 / 0.2e1 + t496 * t684 + t497 * t685 + t444 * t654 + t660 / 0.2e1 + t651 + t770 - t949;
t7 = t575 + t646;
t3 = t908 * pkin(5) + (t625 + t717 + t731) * pkin(5) / 0.2e1 + t907 + t913;
t2 = -t574 + t576 + t583;
t29 = [qJD(2) * t5 + qJD(3) * t32 - qJD(4) * t6 + qJD(5) * t9 + qJD(6) * t12, t21 * qJD(3) + t2 * qJD(4) + t3 * qJD(5) + t7 * qJD(6) + t778 + ((-t536 * mrSges(5,2) + t478 / 0.2e1 + t573 * t522 - t426 * mrSges(5,3)) * t571 + t922 * t445 + 0.2e1 * (t229 * t922 + t230 * t437 + t498 * t552) * t900 + t94 * t927 + t552 * t365 + t498 * t429 - qJ(3) * t501 + t463 * t212 + t437 * t443 + t208 * t298 + 0.2e1 * (-qJ(3) * t536 + t573 * t629) * t902 + t260 * t881 + t258 * t885 + t200 * t994 + (-Ifges(4,4) + Ifges(3,5) + t819 / 0.2e1 - t818 / 0.2e1 + Ifges(6,5) * t858 + Ifges(6,6) * t855 + Ifges(7,5) * t965 + Ifges(7,6) * t994 - pkin(2) * mrSges(4,1)) * t572 + t336 * t855 + t338 * t858 + t434 * t864 + t432 * t867 - t395 * t1001 + t1008 * t300 + 0.2e1 * (t1008 * t93 + t208 * t94 + t395 * t463) * t898 + t93 * t954 + t202 * t965 + (-t536 * mrSges(5,1) - t477 / 0.2e1 + t573 * t524 - t427 * mrSges(5,3)) * t568 + (-qJ(3) * mrSges(4,1) + Ifges(4,5) - Ifges(3,6) + t623) * t569 + (m(4) * t638 - t572 * mrSges(3,1) + t569 * mrSges(3,2) + t532) * pkin(7) - t919 * mrSges(6,3)) * qJD(2), t761 + t21 * qJD(2) + 0.2e1 * ((t347 * t657 - t950) * t898 + (-t743 - t748) * t900) * qJD(3) + t13 * qJD(4) + t18 * qJD(5) + t27 * qJD(6), -t772 + t2 * qJD(2) + t13 * qJD(3) + (-t647 + t660 + t814 - t813 + m(7) * (t101 * t496 + t102 * t497) + t802 - t801 + (t237 * t845 + t238 * t567) * t896 - Ifges(5,5) * t713 + Ifges(5,6) * t722 - t420 * mrSges(5,2) - t421 * mrSges(5,1) + t651 + (-t740 - t742) * mrSges(7,3)) * qJD(4) + t10 * qJD(5) + t15 * qJD(6), t771 + t3 * qJD(2) + t18 * qJD(3) + t10 * qJD(4) + (t651 - t805 - t806 - t830 + t831) * qJD(5) + t23 * qJD(6) + (m(7) * (t566 * t96 + t570 * t95) + (-t716 - t730) * mrSges(7,3)) * t815, t768 + t7 * qJD(2) + t27 * qJD(3) + t15 * qJD(4) + t23 * qJD(5) + (t707 - t834 - t835) * qJD(6); -qJD(3) * t22 - qJD(4) * t1 + qJD(5) * t4 + qJD(6) * t8 - t778, qJD(3) * t115 + qJD(4) * t17 - qJD(5) * t20 - qJD(6) * t33, t628 (m(7) * (t1008 * t497 - t208 * t496) + (-t437 * t845 + t567 * t922) * t896 + t517 * t700 + t622 * t656 - mrSges(5,2) * t712 - mrSges(5,1) * t721 + t639 + t975 * mrSges(7,3) + t1020) * qJD(4) + t25 * qJD(5) + t1025 + t611, t25 * qJD(4) + t1020 * qJD(5) + t1025 + (m(7) * t1014 - t976 * mrSges(7,3)) * t815 + t610 (t971 + t1009) * qJD(6) + t636 + t947 * t1024; qJD(2) * t22 + qJD(4) * t14 + qJD(5) * t19 + qJD(6) * t28 - t761, -t628, 0 ((-t517 * t567 - t622 * t845) * t896 - m(7) * t975 + t597) * qJD(4) + t62 * qJD(5) + t634 + t911, t62 * qJD(4) + (t615 + t844) * qJD(5) + t633 + t911, qJD(6) * t1001 + t947 * t608 + t762; qJD(2) * t1 - qJD(3) * t14 - qJD(5) * t11 + qJD(6) * t16 + t772, qJD(5) * t26 + t1021 - t611, -qJD(5) * t70 - t634 + t912, qJD(5) * t160 + t702 ((t503 * t570 + t504 * t566) * t895 + t593) * qJD(5) + t132 * qJD(6) + t601, t132 * qJD(5) + t600 + t702; -qJD(2) * t4 - qJD(3) * t19 + qJD(4) * t11 + qJD(6) * t24 - t771, -qJD(4) * t26 + t1021 - t610, qJD(4) * t70 - t633 + t912, qJD(6) * t131 - t601, -t516, -t516 - t592; -qJD(2) * t8 - qJD(3) * t28 - qJD(4) * t16 - qJD(5) * t24 - t768, -t1013 * t947 - t636, -t603 * t947 - t762, -qJD(5) * t131 - t600, t592, 0;];
Cq  = t29;
