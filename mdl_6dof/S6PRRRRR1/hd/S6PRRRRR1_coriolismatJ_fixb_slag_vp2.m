% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:37:51
% EndTime: 2019-03-09 00:38:27
% DurationCPUTime: 23.40s
% Computational Cost: add. (39208->746), mult. (85660->973), div. (0->0), fcn. (100665->12), ass. (0->445)
t470 = sin(qJ(3));
t472 = cos(qJ(3));
t467 = sin(pkin(6));
t785 = sin(qJ(2));
t633 = t467 * t785;
t727 = cos(pkin(6));
t522 = t470 * t727 + t472 * t633;
t523 = t470 * t633 - t472 * t727;
t784 = sin(qJ(4));
t787 = cos(qJ(4));
t828 = -t522 * t784 - t787 * t523;
t900 = t828 * mrSges(5,2);
t354 = t522 * t787 - t784 * t523;
t924 = t354 * mrSges(5,1);
t469 = sin(qJ(5));
t786 = cos(qJ(5));
t234 = t354 * t469 - t786 * t828;
t468 = sin(qJ(6));
t471 = cos(qJ(6));
t578 = mrSges(7,1) * t471 - mrSges(7,2) * t468;
t909 = t354 * t786 + t469 * t828;
t936 = t909 * t578;
t946 = t909 * mrSges(6,1);
t550 = -t946 / 0.2e1 - t936 / 0.2e1;
t465 = t471 ^ 2;
t734 = t465 * mrSges(7,3);
t463 = t468 ^ 2;
t735 = t463 * mrSges(7,3);
t843 = t734 / 0.2e1 + t735 / 0.2e1;
t962 = t234 * mrSges(6,2);
t982 = -t234 * t843 + t962 / 0.2e1 + t550;
t987 = -t900 / 0.2e1 - t924 / 0.2e1 + t982;
t986 = qJD(3) + qJD(4);
t782 = pkin(4) * t469;
t454 = pkin(11) + t782;
t662 = t463 + t465;
t597 = t662 * t454;
t812 = m(6) * pkin(4);
t660 = t812 / 0.2e1;
t814 = m(7) / 0.2e1;
t657 = t786 * pkin(4);
t455 = -t657 - pkin(5);
t939 = t455 * t909;
t961 = t234 * t469;
t985 = -(-t234 * t597 + t939) * t814 - (-t786 * t909 - t961) * t660 - t987;
t777 = mrSges(7,2) * t471;
t779 = mrSges(7,1) * t468;
t577 = t777 + t779;
t546 = t234 * t577;
t555 = t777 / 0.2e1 + t779 / 0.2e1;
t960 = t234 * t555;
t969 = -t546 / 0.2e1 + t960;
t473 = cos(qJ(2));
t682 = t467 * t473;
t174 = -t468 * t682 + t471 * t909;
t719 = t174 * t471;
t173 = -t468 * t909 - t471 * t682;
t721 = t173 * t468;
t560 = t719 - t721;
t973 = m(7) * (t909 - t560) * t234;
t980 = t973 * qJD(1);
t984 = -qJD(6) * t969 - t980;
t733 = t468 * mrSges(7,3);
t637 = t733 / 0.2e1;
t638 = -t733 / 0.2e1;
t885 = t637 + t638;
t967 = t960 + t546 / 0.2e1 + t885 * t174;
t968 = t967 * qJD(6) + t980;
t983 = qJD(5) + t986;
t601 = t662 * t234;
t949 = m(7) * (-pkin(5) * t909 - pkin(11) * t601);
t981 = t949 / 0.2e1 + t982;
t912 = t962 - t936 - t946;
t979 = t912 - t900 - t924;
t976 = qJD(1) * t969;
t972 = -t949 / 0.2e1;
t627 = t784 * t472;
t434 = -t470 * t787 - t627;
t540 = t470 * t784 - t472 * t787;
t377 = t434 * t469 - t786 * t540;
t373 = Ifges(6,5) * t377;
t893 = t577 * t377;
t801 = t893 / 0.2e1;
t441 = t468 * Ifges(7,5) + t471 * Ifges(7,6);
t500 = -t434 * t786 - t469 * t540;
t854 = t500 * t441;
t860 = Ifges(6,6) * t500;
t878 = t854 / 0.4e1 - t860 / 0.2e1;
t517 = -t373 / 0.2e1 + pkin(5) * t801 - t878;
t769 = Ifges(7,6) * t500;
t459 = Ifges(7,4) * t471;
t847 = -Ifges(7,2) * t468 + t459;
t161 = t377 * t847 + t769;
t894 = t471 * t161;
t773 = Ifges(7,4) * t468;
t446 = Ifges(7,1) * t471 - t773;
t772 = Ifges(7,5) * t500;
t163 = t377 * t446 + t772;
t896 = t468 * t163;
t880 = -t894 / 0.4e1 - t896 / 0.4e1;
t899 = t377 * t468;
t911 = -mrSges(7,2) * t500 - mrSges(7,3) * t899;
t937 = t471 * t911;
t898 = t377 * t471;
t910 = mrSges(7,1) * t500 - mrSges(7,3) * t898;
t938 = t468 * t910;
t913 = t937 / 0.2e1 - t938 / 0.2e1;
t964 = t913 * pkin(11) - t517 - t880;
t963 = -t234 / 0.2e1;
t804 = t234 / 0.2e1;
t807 = -pkin(9) - pkin(8);
t850 = t787 * t807;
t566 = t470 * t850;
t392 = t627 * t807 + t566;
t433 = t434 * pkin(10);
t849 = t433 + t392;
t588 = t784 * t807;
t494 = t472 * (-t787 * pkin(10) + t850) - (-pkin(10) * t784 + t588) * t470;
t956 = t494 * t786;
t876 = t469 * t849 - t956;
t958 = t234 * t876;
t957 = t469 * t494;
t393 = t472 * t588 + t566;
t954 = t392 - t393;
t840 = -t468 * Ifges(7,1) - t459;
t663 = t471 * t840;
t443 = t471 * Ifges(7,2) + t773;
t674 = t468 * t443;
t548 = t674 / 0.2e1 + t663 / 0.2e1;
t908 = t854 / 0.2e1 - t860 + t894 / 0.2e1 + t896 / 0.2e1;
t518 = -Ifges(5,5) * t540 + Ifges(5,6) * t434 - t377 * t548 + t373 + t908;
t391 = -t470 * t588 + t472 * t850;
t947 = t391 * mrSges(5,1);
t953 = t518 + t947;
t789 = -t471 / 0.2e1;
t790 = -t468 / 0.2e1;
t591 = -t446 * t790 - t789 * t847 - t548;
t917 = t876 * t578;
t923 = t876 * mrSges(6,1);
t229 = t786 * t849 + t957;
t948 = t229 * mrSges(6,2);
t952 = -t917 - t923 - t948;
t841 = t937 - t938;
t824 = t948 / 0.2e1 + t917 / 0.2e1 + t923 / 0.2e1;
t945 = t173 * t910;
t659 = t787 * pkin(3);
t456 = t659 + pkin(4);
t586 = t786 * t784;
t420 = pkin(3) * t586 + t469 * t456;
t944 = t229 * t420;
t943 = t229 * t468;
t942 = t229 * t469;
t941 = t229 * t471;
t718 = t229 * t876;
t856 = t229 * t909;
t736 = t434 * mrSges(5,3);
t940 = t354 * t736;
t457 = -t472 * pkin(3) - pkin(2);
t407 = pkin(4) * t540 + t457;
t877 = mrSges(6,1) * t500 + mrSges(6,2) * t377;
t224 = -pkin(5) * t377 - pkin(11) * t500 + t407;
t94 = t224 * t471 - t468 * t876;
t95 = t224 * t468 + t471 * t876;
t928 = t407 * t877 + t94 * t910 + t95 * t911;
t927 = -t229 / 0.2e1;
t926 = -t377 / 0.4e1;
t902 = mrSges(6,3) * t377;
t922 = mrSges(6,1) + t578;
t921 = t940 / 0.2e1;
t628 = t784 * t469;
t419 = -pkin(3) * t628 + t456 * t786;
t413 = -pkin(5) - t419;
t919 = t413 * t876;
t268 = t577 * t500;
t918 = t876 * t268;
t906 = -t911 / 0.2e1;
t795 = t377 / 0.2e1;
t858 = t500 * mrSges(6,3);
t903 = -t858 / 0.2e1;
t646 = t858 / 0.2e1;
t897 = t467 ^ 2 * t785;
t458 = Ifges(7,5) * t471;
t768 = Ifges(7,6) * t468;
t852 = (t458 / 0.2e1 - t768 / 0.2e1) * t377;
t846 = t663 / 0.4e1 + t674 / 0.4e1;
t891 = t846 * t377;
t379 = -t434 * mrSges(5,1) - mrSges(5,2) * t540;
t890 = t379 + t877;
t414 = pkin(11) + t420;
t600 = t662 * t414;
t889 = t420 - t600;
t887 = t470 ^ 2 + t472 ^ 2;
t848 = t458 - t768;
t871 = t500 / 0.2e1;
t872 = -t377 / 0.2e1;
t886 = t848 * t871 - Ifges(6,4) * t500 + (Ifges(6,2) + Ifges(7,3)) * t872;
t703 = t909 * t413;
t884 = t560 * t419 + t703;
t287 = pkin(5) * t500 - pkin(11) * t377;
t691 = t500 * t468;
t275 = mrSges(7,2) * t377 - mrSges(7,3) * t691;
t690 = t500 * t471;
t278 = -mrSges(7,1) * t377 - mrSges(7,3) * t690;
t549 = t275 * t790 + t278 * t789;
t861 = mrSges(7,3) * (t465 / 0.2e1 + t463 / 0.2e1);
t875 = t500 * t861 - t549;
t269 = t443 * t500;
t270 = t840 * t500;
t164 = -Ifges(7,5) * t377 + t446 * t500;
t670 = t471 * t164;
t162 = -Ifges(7,6) * t377 + t500 * t847;
t681 = t468 * t162;
t874 = t577 * t927 - t471 * t269 / 0.4e1 + t670 / 0.4e1 + t468 * t270 / 0.4e1 - t681 / 0.4e1 + t848 * t926 - (-t840 + t847) * t691 / 0.4e1;
t767 = Ifges(7,3) * t500;
t870 = -t767 / 0.2e1;
t868 = m(7) * t413;
t867 = m(7) * t455;
t830 = mrSges(7,3) * t662;
t666 = t471 * t275;
t676 = t468 * t278;
t845 = -t666 / 0.2e1 + t676 / 0.2e1;
t112 = t287 * t471 - t943;
t113 = t287 * t468 + t941;
t561 = -t112 * t468 + t113 * t471;
t783 = pkin(4) * t434;
t238 = t287 - t783;
t461 = t470 * pkin(3);
t225 = t238 + t461;
t350 = t433 + t393;
t230 = t350 * t786 + t957;
t96 = t225 * t471 - t230 * t468;
t97 = t225 * t468 + t230 * t471;
t571 = -t96 * t468 + t97 * t471;
t440 = t470 * mrSges(4,1) + t472 * mrSges(4,2);
t839 = t734 + t735;
t838 = -mrSges(4,1) * t472 + mrSges(4,2) * t470;
t837 = -t268 / 0.2e1 + t903;
t836 = Ifges(6,4) * t377 + Ifges(6,1) * t871 - t681 / 0.2e1 + t670 / 0.2e1;
t395 = t434 * t682;
t396 = t540 * t682;
t299 = -t395 * t786 - t396 * t469;
t834 = -t922 * t299 / 0.2e1;
t833 = t396 * mrSges(5,2) / 0.2e1 + t395 * mrSges(5,1) / 0.2e1;
t832 = (-t784 * mrSges(5,1) - t787 * mrSges(5,2)) * pkin(3);
t284 = -mrSges(6,1) * t377 + mrSges(6,2) * t500;
t831 = -m(6) * t407 - t284;
t829 = t473 * t887;
t827 = t893 * t963 + t174 * t906 - t945 / 0.2e1;
t823 = -t786 * t812 + t867;
t300 = t469 * t395 - t396 * t786;
t289 = t471 * t300 + t468 * t633;
t700 = t289 * t471;
t288 = -t468 * t300 + t471 * t633;
t701 = t288 * t468;
t753 = t300 * mrSges(6,2);
t822 = t753 / 0.2e1 + (-t700 / 0.2e1 + t701 / 0.2e1) * mrSges(7,3);
t821 = m(7) * t597 + t469 * t812 + t839;
t572 = t468 * t94 - t471 * t95;
t815 = -m(7) / 0.2e1;
t817 = -m(6) / 0.2e1;
t820 = t876 * t817 - t572 * t815 - t902 / 0.2e1 + t845;
t818 = (t663 + t674) * t926 + t964;
t816 = m(6) / 0.2e1;
t813 = m(5) * pkin(3);
t811 = m(7) * pkin(4);
t810 = -mrSges(7,2) / 0.2e1;
t809 = t96 / 0.2e1;
t808 = -t97 / 0.2e1;
t806 = -t112 / 0.2e1;
t805 = t113 / 0.2e1;
t799 = t268 / 0.2e1;
t798 = t911 / 0.2e1;
t797 = t275 / 0.2e1;
t796 = -t278 / 0.2e1;
t794 = -t414 / 0.2e1;
t793 = t443 / 0.4e1;
t792 = t840 / 0.4e1;
t791 = t455 / 0.2e1;
t788 = t471 / 0.2e1;
t775 = mrSges(7,3) * t500;
t227 = t350 * t469 - t956;
t765 = t227 * mrSges(6,1);
t762 = t230 * mrSges(6,2);
t745 = t392 * mrSges(5,2);
t744 = t393 * mrSges(5,2);
t740 = t419 * mrSges(6,2);
t739 = t420 * mrSges(6,1);
t424 = (t469 * t787 + t586) * pkin(3);
t738 = t424 * mrSges(6,1);
t425 = (t786 * t787 - t628) * pkin(3);
t737 = t425 * mrSges(6,2);
t100 = t238 * t471 - t943;
t726 = t100 * t468;
t101 = t238 * t468 + t941;
t725 = t101 * t471;
t716 = t227 * t229;
t715 = t227 * t578;
t714 = t229 * t299;
t713 = t229 * t424;
t705 = t234 * t299;
t704 = t234 * t424;
t587 = t473 * t897;
t39 = m(7) * (t173 * t288 + t174 * t289 + t705) + m(6) * (t300 * t909 - t587 + t705) + m(5) * (-t354 * t396 + t395 * t828 - t587) + m(4) * (t829 * t897 - t587);
t689 = t39 * qJD(1);
t265 = t578 * t500;
t687 = t413 * t265;
t686 = t413 * t893;
t685 = t420 * t578;
t684 = t424 * t578;
t683 = t455 * t893;
t661 = mrSges(6,1) * t782;
t655 = mrSges(7,3) * t725;
t654 = t782 / 0.2e1;
t653 = Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t650 = t902 / 0.2e1;
t636 = mrSges(7,3) * t788;
t634 = t767 / 0.2e1 + t852;
t632 = t468 * t786;
t631 = t471 * t786;
t629 = t784 * t393;
t622 = -t690 / 0.4e1;
t621 = t690 / 0.4e1;
t620 = -t682 / 0.2e1;
t619 = t682 / 0.2e1;
t599 = t662 * t419;
t598 = t662 * t425;
t596 = t858 * t782;
t594 = mrSges(6,2) * t657;
t593 = t657 / 0.2e1;
t585 = -t632 / 0.2e1;
t582 = t657 * t902;
t581 = t100 * t638 + t101 * t636 - t824;
t580 = t97 * t636 + t96 * t638 - t715 / 0.2e1 - t765 / 0.2e1 - t762 / 0.2e1;
t412 = t461 - t783;
t527 = t234 * t650 + t620 * t890 + t903 * t909 - t827 + t921;
t558 = t227 * t234 - t856;
t474 = m(5) * (-t354 * t954 - t461 * t682) / 0.2e1 + t527 + (t173 * t96 + t174 * t97 + t234 * t572 + t558) * t814 + t440 * t620 + t666 * t963 - t921 + (t230 * t909 - t412 * t682 + t558 - t958) * t816 + (t676 - t902) * t804 + (t799 + t646) * t909;
t556 = t700 - t701;
t484 = (-t299 * t419 + t300 * t420) * t817 + (t299 * t413 + t414 * t556) * t815 - (t395 * t787 - t396 * t784) * t813 / 0.2e1 + t440 * t619;
t2 = t474 + t484 + t822 - t833 - t834;
t380 = mrSges(5,1) * t540 - t434 * mrSges(5,2);
t428 = Ifges(5,4) * t540;
t477 = t457 * t379 + (t848 * t872 + t836) * t377 + (-t377 * t653 + Ifges(6,1) * t795 + (Ifges(7,1) * t898 + t772) * t788 + (-Ifges(7,2) * t899 + t769) * t790 - t898 * t773 + t886) * t500 - (t902 + t893) * t229 + t928;
t493 = -Ifges(5,4) * t434 + (Ifges(5,1) - Ifges(5,2)) * t540;
t3 = m(5) * (t391 * t954 + t457 * t461) + m(6) * (t230 * t876 + t407 * t412 - t716) + m(7) * (t94 * t96 + t95 * t97 - t716) + (t230 * t377 + (t227 - t876) * t500) * mrSges(6,3) - pkin(2) * t440 + t412 * t284 + t477 + t97 * t275 + t96 * t278 + t227 * t268 + (Ifges(4,4) * t472 + (Ifges(4,1) - Ifges(4,2)) * t470 + (-mrSges(5,3) * t954 - t428) * t787) * t472 + t493 * t434 + (t784 * t428 - Ifges(4,4) * t470 + pkin(3) * t380 + (t392 * t784 - t629) * mrSges(5,3)) * t470;
t574 = t2 * qJD(1) + t3 * qJD(2);
t5 = m(7) * (t100 * t94 + t101 * t95 - t718) + t540 * t428 + t229 * t902 + t477 + t101 * t275 + t100 * t278 + t918 + (pkin(4) * t831 + t493) * t434;
t559 = -t856 + t958;
t480 = (t682 * t783 + t559 + t856) * t817 + (t100 * t173 + t101 * t174 + t559) * t815 + t837 * t909 - t820 * t234;
t568 = t288 * t638 + t289 * t636 - t753 / 0.2e1 + t834;
t543 = t568 + t833;
t488 = (t299 * t455 + t454 * t556) * t814 + (-t299 * t786 + t300 * t469) * t660 + t543;
t8 = t890 * t619 + t646 * t909 + t902 * t963 + t480 + t488 + t827;
t573 = -t8 * qJD(1) + t5 * qJD(2);
t551 = pkin(5) * t577;
t567 = -t551 / 0.2e1 + t591;
t13 = m(7) * (t112 * t94 + t113 * t95 - t718) + t918 - t229 * t893 + t113 * t275 + t112 * t278 + (t161 * t790 + t163 * t788 + t886) * t500 - (t852 + (-Ifges(6,1) / 0.2e1 + t653) * t500 - t836) * t377 + t928;
t482 = (t801 + (t876 + t572) * t814 + t845) * t234 + (t112 * t173 + t113 * t174 - t856) * t814 + t945 / 0.2e1 + t174 * t798 + t909 * t799 + t877 * t620;
t513 = m(7) * (-pkin(5) * t299 + pkin(11) * t556);
t15 = (t578 / 0.2e1 + mrSges(6,1) / 0.2e1) * t299 - t513 / 0.2e1 + t482 + t822;
t564 = t15 * qJD(1) + t13 * qJD(2);
t18 = -t229 * t265 + t94 * t275 - t95 * t278 + (t572 * mrSges(7,3) + t162 * t789 + t270 * t788 + t441 * t795 + (t164 - t269) * t790) * t500;
t496 = (-t719 / 0.2e1 + t721 / 0.2e1) * t775 + t173 * t797 + t174 * t796 + t265 * t804;
t552 = t288 * mrSges(7,1) / 0.2e1 + t289 * t810;
t27 = t496 - t552;
t563 = t27 * qJD(1) + t18 * qJD(2);
t562 = t725 - t726;
t557 = t919 - t944;
t553 = -t100 * mrSges(7,1) / 0.2e1 + t101 * mrSges(7,2) / 0.2e1;
t545 = t413 * t577;
t544 = t455 * t577;
t542 = t551 / 0.2e1;
t541 = t662 * t786;
t537 = t555 * t419;
t536 = t555 * t425;
t534 = -t545 / 0.2e1;
t24 = (t804 + t963) * mrSges(6,2) + (t234 * t889 + t884) * t814 + t972;
t506 = (t441 / 0.4e1 - Ifges(6,6) / 0.2e1) * t500 + t413 * t801 + t420 * t799;
t516 = m(7) * (-pkin(5) * t227 + pkin(11) * t571);
t530 = t163 / 0.4e1 + t910 * t794 + t419 * t796;
t531 = t161 / 0.4e1 + t414 * t798 + t419 * t797;
t6 = (t230 / 0.2e1 + t927) * mrSges(6,2) + (-t769 / 0.4e1 + pkin(11) * t906 + (t792 - t459 / 0.4e1) * t377 + t531) * t471 + (-t772 / 0.4e1 + pkin(11) * t910 / 0.2e1 + (t793 + t773 / 0.4e1 + (Ifges(7,2) / 0.4e1 - Ifges(7,1) / 0.4e1) * t471) * t377 + t530) * t468 - t516 / 0.2e1 + (t414 * t561 - t419 * t572 + t557) * t814 - (-Ifges(6,5) / 0.2e1 + t846) * t377 + ((t805 + t808) * t471 + (t806 + t809) * t468) * mrSges(7,3) + t506 + t517 - t922 * (t876 / 0.2e1 - t227 / 0.2e1);
t501 = mrSges(7,3) * t599 - t685 - t739 - t740;
t67 = m(7) * (t413 * t420 + t414 * t599) + t501;
t533 = t24 * qJD(1) + t6 * qJD(2) + t67 * qJD(3);
t476 = (-t419 * t876 - t713 + t944) * t817 + (-t713 + t919) * t815 + t745 / 0.2e1 - t947 / 0.2e1 - t686 / 0.2e1 + t419 * t650 + t420 * t646 + t837 * t424 + (t562 * t815 - t913) * t414 + t820 * t425;
t481 = t227 * t455 * t814 + t947 / 0.2e1 - t744 / 0.2e1 + t683 / 0.2e1 + (-t227 * t786 + t230 * t469) * t660 - t596 / 0.2e1 - t582 / 0.2e1 + t580 + (t571 * t814 + t913) * t454;
t498 = t100 * t637 - t655 / 0.2e1 + t824;
t10 = t476 + t481 + t498;
t478 = (-t234 * t420 + t704 + (-t419 + t425) * t909) * t816 + (-t234 * t600 + t425 * t560 + t703 + t704) * t814 + t987;
t17 = t478 + t985;
t69 = -t922 * t424 + t832 + (-mrSges(6,2) + t830) * t425 + m(7) * (t413 * t424 + t414 * t598) + m(6) * (-t419 * t424 + t420 * t425);
t532 = t17 * qJD(1) - t10 * qJD(2) + t69 * qJD(3);
t521 = mrSges(7,1) * t809 + mrSges(7,2) * t808 + t634;
t19 = -t687 / 0.2e1 + t446 * t622 + t443 * t621 + t875 * t414 + t521 - t874;
t248 = t545 + t591;
t528 = -t19 * qJD(2) + t248 * qJD(3) - t976;
t520 = mrSges(5,3) * t540;
t519 = mrSges(7,1) * t806 + mrSges(7,2) * t805 + t870;
t503 = Ifges(6,5) * t795 - t824;
t475 = (t455 * t876 + (t631 * t95 - t632 * t94 - t942) * pkin(4)) * t814 + t893 * t791 + t112 * t638 + t113 * t636 + t268 * t654 + pkin(4) * t278 * t585 + t593 * t666 + t503 + (t561 * t814 + t913) * t454 - t891 + t878 - t880;
t514 = m(7) * (-pkin(5) * t876 + pkin(11) * t562);
t11 = t475 - t514 / 0.2e1 + t498 + t891 - t964;
t492 = -t578 * t782 + t657 * t830 - t594 - t661;
t249 = (t454 * t541 + t455 * t469) * t811 + t492;
t483 = (t939 - t454 * t601 + (-t173 * t632 + t174 * t631 + t961) * pkin(4)) * t814 + t982;
t30 = t972 + t483 - t982;
t479 = (t455 * t420 + (t413 * t469 + t414 * t541) * pkin(4)) * t814 - t740 / 0.2e1 - t739 / 0.2e1 - t685 / 0.2e1 - t661 / 0.2e1 - t578 * t654 - t594 / 0.2e1 + t593 * t830 + (t597 * t814 + t843) * t419;
t497 = (-pkin(5) * t424 + pkin(11) * t598) * t814 - t738 / 0.2e1 - t684 / 0.2e1 - t737 / 0.2e1;
t54 = t425 * t843 - t479 + t497;
t515 = t30 * qJD(1) + t11 * qJD(2) - t54 * qJD(3) + t249 * qJD(4);
t495 = -t544 / 0.2e1 - t591;
t134 = t534 - t536 + t495;
t504 = t443 * t622 + t446 * t621 + t885 * t95 + t874;
t486 = t265 * t791 - t454 * t875 + t504;
t22 = t486 + t553 + t870 - t852;
t281 = t544 + t591;
t512 = t22 * qJD(2) - t134 * qJD(3) + t281 * qJD(4) - t976;
t136 = t534 + t542 - t537 - t591;
t505 = (mrSges(7,1) * t585 + t631 * t810) * pkin(4);
t171 = t542 + t505 + t495;
t487 = -t875 * pkin(11) + t504 - pkin(5) * t265 / 0.2e1;
t26 = t487 + t519 - t852;
t291 = -t551 + t591;
t511 = t26 * qJD(2) - t136 * qJD(3) - t171 * qJD(4) + t291 * qJD(5);
t417 = t544 / 0.2e1;
t381 = t545 / 0.2e1;
t172 = t417 + t505 + t567;
t137 = t381 - t537 + t567;
t135 = t381 + t417 - t536 + t591;
t55 = t425 * t861 + t479 + t497;
t29 = t483 + t981;
t28 = t496 + t552;
t25 = Ifges(7,5) * t898 / 0.2e1 - Ifges(7,6) * t899 / 0.2e1 + t487 - t519;
t23 = t884 * t814 + (t889 * t814 + mrSges(6,2) / 0.2e1 - t861) * t234 + t550 + t981;
t21 = t486 - t553 + t634;
t20 = t504 + t687 / 0.2e1 + t521 + t662 * t775 * t794 + t549 * t414;
t16 = t478 - t985;
t14 = t513 / 0.2e1 + t482 + t568;
t12 = t475 + t514 / 0.2e1 + t581 + t818;
t9 = -t940 / 0.2e1 + t488 - t480 + t527;
t7 = t580 + t516 / 0.2e1 + t557 * t814 + t506 + ((-t112 * t414 - t419 * t94) * t814 + mrSges(7,3) * t806 - t377 * t793 + t530) * t468 + (mrSges(7,3) * t805 + (t113 * t414 + t419 * t95) * t814 - t377 * t792 + t531) * t471 + t503 + t818;
t4 = t518 - t476 + t481 + t581;
t1 = t474 - t484 + t543;
t31 = [t39 * qJD(2) + t973 * t983, t1 * qJD(3) + t9 * qJD(4) + t14 * qJD(5) + t28 * qJD(6) + t689 + (t289 * t275 + t288 * t278 + m(7) * (t288 * t94 + t289 * t95 - t714) + t300 * t902 + m(6) * (t300 * t876 - t714) + t395 * t736 + t396 * t520 + m(5) * (t391 * t396 + t392 * t395) + m(4) * (-t785 * pkin(2) + pkin(8) * t829) * t467 + (m(5) * t457 - mrSges(3,1) + t380 - t831 + t838) * t633 + (t268 + t858) * t299 + (mrSges(4,3) * t887 - mrSges(3,2)) * t682) * qJD(2), t1 * qJD(2) + ((-t354 * t787 + t784 * t828) * t813 + t523 * mrSges(4,2) - t522 * mrSges(4,1) + (-m(6) * t419 + t868) * t909 - (m(6) * t420 + m(7) * t600 + t839) * t234 + t979) * qJD(3) + t16 * qJD(4) + t23 * qJD(5) + t968, t9 * qJD(2) + t16 * qJD(3) + (-t234 * t821 + t823 * t909 + t979) * qJD(4) + t29 * qJD(5) + t968, t14 * qJD(2) + t23 * qJD(3) + t29 * qJD(4) + (-mrSges(7,3) * t601 + t912 + t949) * qJD(5) + t968, t28 * qJD(2) + (-mrSges(7,1) * t174 - mrSges(7,2) * t173) * qJD(6) + t983 * t967; qJD(3) * t2 - qJD(4) * t8 + qJD(5) * t15 + qJD(6) * t27 - t689, qJD(3) * t3 + qJD(4) * t5 + qJD(5) * t13 + qJD(6) * t18 ((t391 * t787 + t629) * t813 + m(6) * (-t227 * t419 + t230 * t420) + t520 * t659 + t784 * pkin(3) * t736 + Ifges(4,5) * t472 - Ifges(4,6) * t470 - t715 + t686 - t744 - t762 - t765 + t227 * t868 - t420 * t858 - t419 * t902 + (m(7) * t571 + t841) * t414 + t838 * pkin(8) + t571 * mrSges(7,3) + t953) * qJD(3) + t4 * qJD(4) + t7 * qJD(5) + t20 * qJD(6) + t574, t4 * qJD(3) + ((-t786 * t876 + t942) * t812 - t596 - t582 + t683 - t745 + t876 * t867 + t655 - mrSges(7,3) * t726 + (m(7) * t562 + t841) * t454 + t952 + t953) * qJD(4) + t12 * qJD(5) + t21 * qJD(6) + t573, t7 * qJD(3) + t12 * qJD(4) + t25 * qJD(6) + t564 + (-(-Ifges(6,5) + t548) * t377 + (-m(7) * t876 - t893) * pkin(5) + (m(7) * t561 + t841) * pkin(11) + t561 * mrSges(7,3) + t908 + t952) * qJD(5), t20 * qJD(3) + t21 * qJD(4) + t25 * qJD(5) + (-mrSges(7,1) * t95 - mrSges(7,2) * t94 - t854) * qJD(6) + t563; -qJD(2) * t2 + qJD(4) * t17 + qJD(5) * t24 + t984, -qJD(4) * t10 + qJD(5) * t6 - qJD(6) * t19 - t574, qJD(4) * t69 + qJD(5) * t67 + qJD(6) * t248 (t424 * t823 + t425 * t821 - t684 - t737 - t738 + t832) * qJD(4) + t55 * qJD(5) + t135 * qJD(6) + t532, t55 * qJD(4) + (m(7) * (-pkin(5) * t420 + pkin(11) * t599) + t501) * qJD(5) + t137 * qJD(6) + t533, t135 * qJD(4) + t137 * qJD(5) + (-t414 * t578 + t848) * qJD(6) + t528; qJD(2) * t8 - qJD(3) * t17 + qJD(5) * t30 + t984, qJD(3) * t10 + qJD(5) * t11 + qJD(6) * t22 - t573, -qJD(5) * t54 - qJD(6) * t134 - t532, qJD(5) * t249 + qJD(6) * t281 ((-pkin(5) * t469 + pkin(11) * t541) * t811 + t492) * qJD(5) + t172 * qJD(6) + t515, t172 * qJD(5) + (-t454 * t578 + t848) * qJD(6) + t512; -qJD(2) * t15 - qJD(3) * t24 - qJD(4) * t30 - t980, -qJD(3) * t6 - qJD(4) * t11 + qJD(6) * t26 - t564, qJD(4) * t54 - qJD(6) * t136 - t533, -qJD(6) * t171 - t515, t291 * qJD(6) (-pkin(11) * t578 + t848) * qJD(6) + t511; -t27 * qJD(2) + t969 * t986, qJD(3) * t19 - qJD(4) * t22 - qJD(5) * t26 - t563, qJD(4) * t134 + qJD(5) * t136 - t528, qJD(5) * t171 - t512, -t511, 0;];
Cq  = t31;
