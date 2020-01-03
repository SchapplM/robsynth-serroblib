% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:48
% EndTime: 2019-12-05 16:48:46
% DurationCPUTime: 44.45s
% Computational Cost: add. (80389->867), mult. (112721->1255), div. (0->0), fcn. (121234->8), ass. (0->503)
t968 = Icges(5,4) + Icges(6,4);
t955 = Icges(5,5) + Icges(6,5);
t950 = Icges(5,6) + Icges(6,6);
t956 = Icges(5,1) + Icges(6,1);
t967 = Icges(5,2) + Icges(6,2);
t961 = Icges(5,3) + Icges(6,3);
t586 = qJ(3) + qJ(4);
t579 = sin(t586);
t580 = cos(t586);
t588 = cos(pkin(8));
t587 = sin(pkin(8));
t592 = cos(qJ(2));
t769 = t587 * t592;
t518 = -t579 * t769 - t580 * t588;
t966 = t968 * t518;
t764 = t588 * t592;
t520 = -t579 * t764 + t580 * t587;
t965 = t968 * t520;
t964 = -t950 * t579 + t955 * t580;
t521 = t579 * t587 + t580 * t764;
t963 = t968 * t521;
t519 = -t588 * t579 + t580 * t769;
t962 = t968 * t519;
t590 = sin(qJ(2));
t770 = t587 * t590;
t937 = t967 * t518 + t950 * t770 + t962;
t765 = t588 * t590;
t936 = t967 * t520 + t950 * t765 + t963;
t935 = t956 * t519 + t955 * t770 + t966;
t934 = t956 * t521 + t955 * t765 + t965;
t960 = t968 * t580;
t959 = t968 * t579;
t929 = t964 * t590 - t961 * t592;
t958 = -t579 * t967 + t960;
t957 = t956 * t580 - t959;
t939 = t950 * t518 + t955 * t519 + t961 * t770;
t938 = t950 * t520 + t955 * t521 + t961 * t765;
t953 = -t936 * t579 + t934 * t580;
t952 = -t937 * t579 + t935 * t580;
t928 = t958 * t590 - t950 * t592;
t927 = t957 * t590 - t955 * t592;
t949 = -t929 * t587 - t952;
t948 = t929 * t588 + t953;
t947 = t955 * t518 - t950 * t519;
t946 = t955 * t520 - t950 * t521;
t945 = -t521 * t967 + t934 + t965;
t944 = -t519 * t967 + t935 + t966;
t943 = t956 * t520 - t936 - t963;
t942 = t956 * t518 - t937 - t962;
t941 = t938 * t592;
t940 = t939 * t592;
t933 = t928 * t587;
t932 = t928 * t588;
t931 = t927 * t587;
t930 = t927 * t588;
t926 = (-t579 * t955 - t950 * t580) * t590;
t925 = -t948 * t590 + t941;
t924 = t949 * t590 + t940;
t923 = t928 * t579 - t927 * t580;
t585 = t590 ^ 2;
t869 = t592 ^ 2;
t760 = t590 * t592;
t922 = -t944 * t518 - t942 * t519 - t947 * t770;
t921 = t945 * t518 + t943 * t519 + t946 * t770;
t920 = t944 * t520 + t942 * t521 + t947 * t765;
t919 = t945 * t520 + t943 * t521 + t946 * t765;
t918 = t952 * t590 - t940;
t917 = t953 * t590 - t941;
t915 = t933 * t518 + t931 * t519 - t924 * t587;
t914 = -t932 * t518 - t930 * t519 + t925 * t587;
t913 = -t933 * t520 - t931 * t521 + t924 * t588;
t912 = -t932 * t520 - t930 * t521 + t925 * t588;
t911 = t937 * t518 + t935 * t519 + t939 * t770;
t910 = t936 * t518 + t934 * t519 + t938 * t770;
t909 = t937 * t520 + t935 * t521 + t939 * t765;
t908 = t936 * t520 + t934 * t521 + t938 * t765;
t907 = t928 * t518 + t927 * t519 + t929 * t770;
t906 = t928 * t520 + t927 * t521 + t929 * t765;
t591 = cos(qJ(3));
t581 = t591 * pkin(3);
t564 = pkin(4) * t580 + t581;
t710 = t564 - t581;
t425 = qJ(5) * t590 + t592 * t710;
t825 = pkin(4) * t579;
t589 = sin(qJ(3));
t826 = pkin(3) * t589;
t563 = t825 + t826;
t680 = -t563 + t826;
t752 = rSges(6,1) * t519 + rSges(6,2) * t518 + rSges(6,3) * t770 + t425 * t587 + t588 * t680;
t665 = rSges(6,1) * t580 - rSges(6,2) * t579;
t879 = (-rSges(6,3) - qJ(5)) * t592 + (t665 + t710) * t590;
t905 = t950 * t590 + t958 * t592;
t904 = -t590 * t955 - t957 * t592;
t903 = (t961 * t590 + t964 * t592 + t923) * t592;
t902 = (-t580 * t967 - t959) * t590 + t927;
t901 = (t956 * t579 + t960) * t590 + t928;
t900 = t926 * t592;
t899 = t929 * t592;
t898 = t921 * t587 + t588 * t922;
t897 = t587 * t919 - t588 * t920;
t896 = t948 * t592 + (t932 * t579 - t930 * t580 + t938) * t590;
t895 = -t949 * t592 + (t933 * t579 - t931 * t580 + t939) * t590;
t894 = t923 * t590 + t899;
t893 = t587 * t918 + t588 * t917;
t892 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t760 + (0.2e1 * t869 - 0.2e1 * t585) * Icges(3,4);
t861 = m(5) / 0.2e1;
t860 = m(6) / 0.2e1;
t838 = t587 / 0.2e1;
t837 = -t588 / 0.2e1;
t835 = t590 / 0.2e1;
t834 = -t592 / 0.2e1;
t891 = (-t905 * t518 + t904 * t519) * t592 + t907 * t590 + (t914 * t590 + t910 * t592) * t588 + ((-t899 + t911) * t592 + (-t903 - t915) * t590) * t587;
t890 = (-t905 * t520 + t904 * t521) * t592 + t906 * t590 + ((-t899 + t908) * t592 + (-t903 + t912) * t590) * t588 + (t913 * t590 + t909 * t592) * t587;
t889 = (-t902 * t518 + t901 * t519) * t592 + (t921 * t588 + (-t900 - t922) * t587) * t590;
t888 = (-t902 * t520 + t901 * t521) * t592 + ((-t900 + t919) * t588 + t920 * t587) * t590;
t653 = -Icges(3,5) * t590 - Icges(3,6) * t592;
t551 = t653 * t587;
t552 = t653 * t588;
t887 = t914 * t587 + t915 * t588;
t886 = t912 * t587 - t913 * t588;
t885 = -t907 * t592 + (t911 * t587 + t910 * t588) * t590;
t884 = -t906 * t592 + (t909 * t587 + t908 * t588) * t590;
t883 = -t947 * t592 + (-t944 * t579 + t942 * t580) * t590;
t882 = -t946 * t592 + (-t945 * t579 + t943 * t580) * t590;
t881 = t752 * t592;
t751 = rSges(6,1) * t521 + rSges(6,2) * t520 + rSges(6,3) * t765 + t425 * t588 - t587 * t680;
t880 = t879 * t587;
t878 = rSges(6,3) * t590 + t592 * t665 + t425;
t583 = t587 ^ 2;
t584 = t588 ^ 2;
t874 = t583 + t584;
t484 = pkin(7) * t590 + t581 * t592;
t766 = t588 * t589;
t703 = pkin(3) * t766;
t415 = t484 * t587 - t703;
t376 = t415 * t765;
t771 = t587 * t589;
t704 = pkin(3) * t771;
t416 = t484 * t588 + t704;
t697 = -t416 - t751;
t754 = t752 * t765;
t176 = t697 * t770 + t376 + t754;
t246 = t879 * t770 + t881;
t387 = t592 * t415;
t620 = pkin(7) * t592 - t581 * t590;
t734 = -t620 * t770 + t387;
t215 = t246 + t734;
t669 = t697 * t592;
t693 = t620 - t879;
t216 = t693 * t765 + t669;
t367 = rSges(5,1) * t519 + rSges(5,2) * t518 + rSges(5,3) * t770;
t343 = t367 * t765;
t369 = rSges(5,1) * t521 + rSges(5,2) * t520 + rSges(5,3) * t765;
t740 = -t369 - t416;
t242 = t740 * t770 + t343 + t376;
t351 = t592 * t367;
t666 = rSges(5,1) * t580 - rSges(5,2) * t579;
t501 = -rSges(5,3) * t592 + t590 * t666;
t320 = t501 * t770 + t351;
t272 = t320 + t734;
t678 = t740 * t592;
t724 = t620 - t501;
t273 = t724 * t765 + t678;
t509 = t518 * pkin(4);
t510 = t520 * pkin(4);
t409 = rSges(6,1) * t518 - rSges(6,2) * t519;
t411 = rSges(6,1) * t520 - rSges(6,2) * t521;
t736 = t587 * t409 + t588 * t411;
t285 = t509 * t587 + t510 * t588 + t736;
t410 = rSges(5,1) * t518 - rSges(5,2) * t519;
t412 = rSges(5,1) * t520 - rSges(5,2) * t521;
t315 = t587 * t410 + t588 * t412;
t540 = (-rSges(6,1) * t579 - rSges(6,2) * t580) * t590;
t670 = t590 * t825 - t540;
t442 = t670 * t587;
t443 = t670 * t588;
t541 = (-rSges(5,1) * t579 - rSges(5,2) * t580) * t590;
t571 = t592 * pkin(2) + t590 * pkin(6);
t711 = t874 * t571;
t674 = t587 * t415 + t588 * t416 + t711;
t231 = t367 * t587 + t369 * t588 + t674;
t371 = t410 * t765;
t310 = -t412 * t770 + t371;
t330 = t592 * t410 + t541 * t770;
t701 = t541 * t765;
t331 = -t412 * t592 - t701;
t570 = t590 * pkin(2) - t592 * pkin(6);
t690 = -t570 + t724;
t334 = t690 * t587;
t336 = t690 * t588;
t699 = t310 * t231 + t330 * t336 + t331 * t334;
t166 = t587 * t752 + t588 * t751 + t674;
t370 = t409 * t765;
t729 = -t411 - t510;
t276 = t370 + (t509 * t588 + t587 * t729) * t590;
t705 = t585 * t825;
t735 = t592 * t409 + t540 * t770;
t300 = t509 * t592 - t587 * t705 + t735;
t301 = t729 * t592 + (-t540 * t590 + t705) * t588;
t672 = -t570 + t693;
t302 = t672 * t587;
t304 = t672 * t588;
t700 = t276 * t166 + t300 * t304 + t301 * t302;
t821 = (t176 * t285 + t215 * t443 + t216 * t442 + t700) * t860 + (t242 * t315 + (-t272 * t588 - t273 * t587) * t541 + t699) * t861;
t217 = -t751 * t770 + t754;
t763 = t589 * t592;
t631 = t587 * t763 + t588 * t591;
t545 = t631 * pkin(3);
t372 = -t563 * t769 - t564 * t588 + t545;
t632 = -t587 * t591 + t588 * t763;
t546 = t632 * pkin(3);
t373 = -t563 * t764 + t564 * t587 + t546;
t717 = -t587 * t545 - t588 * t546;
t230 = t372 * t587 + t373 * t588 + t717 + t736;
t739 = -t373 - t411;
t749 = t372 * t765 + t370;
t243 = t739 * t770 + t749;
t679 = t751 * t592;
t247 = -t765 * t879 - t679;
t504 = t680 * t590;
t277 = t592 * t372 + t504 * t770 + t735;
t718 = -t504 - t540;
t677 = t718 * t588;
t668 = t590 * t677;
t278 = t592 * t739 + t668;
t289 = t717 + t315;
t349 = t590 * t369;
t297 = -t349 * t587 + t343;
t785 = t369 * t592;
t321 = -t501 * t765 - t785;
t565 = t590 * t704;
t340 = t587 * t718 + t565;
t567 = t590 * t703;
t341 = t567 + t677;
t444 = -t541 * t587 + t565;
t445 = -t541 * t588 + t567;
t822 = (t166 * t243 + t217 * t230 + t246 * t341 + t247 * t340 + t277 * t304 + t278 * t302) * t860 + (t289 * t297 + t320 * t445 + t321 * t444 + t699) * t861;
t875 = t821 - t822;
t775 = (-Icges(4,5) * t589 - Icges(4,6) * t591) * t760;
t808 = Icges(4,4) * t591;
t656 = -Icges(4,2) * t589 + t808;
t524 = -Icges(4,6) * t592 + t590 * t656;
t716 = -t524 + (-Icges(4,1) * t589 - t808) * t590;
t809 = Icges(4,4) * t589;
t661 = Icges(4,1) * t591 - t809;
t526 = -Icges(4,5) * t592 + t590 * t661;
t715 = t526 + (-Icges(4,2) * t591 - t809) * t590;
t675 = t881 * t588 - t880 * t765;
t439 = t620 * t588;
t738 = t879 * t588;
t695 = -t439 + t738;
t438 = t620 * t587;
t737 = t588 * t387 + t438 * t765;
t118 = (t590 * t695 + t669) * t587 + t675 + t737;
t634 = -t880 * t592 + t879 * t769 + t878 * t770;
t694 = t592 * t438 + t484 * t770 - t620 * t769;
t147 = (-t415 - t752) * t590 + t634 + t694;
t386 = t590 * t416;
t692 = -t484 - t878;
t753 = t751 * t590;
t148 = t386 + t695 * t592 + (t590 * t692 + t592 * t693) * t588 + t753;
t157 = (t590 * t738 - t679) * t587 + t675;
t461 = t501 * t588;
t725 = -t439 + t461;
t459 = t501 * t587;
t750 = t588 * t351 - t459 * t765;
t175 = (t590 * t725 + t678) * t587 + t737 + t750;
t177 = -t590 * t752 + t634;
t178 = t738 * t592 + (-t590 * t878 - t592 * t879) * t588 + t753;
t503 = rSges(5,3) * t590 + t592 * t666;
t691 = -t592 * t459 + t501 * t769 + t503 * t770;
t203 = (-t367 - t415) * t590 + t691 + t694;
t723 = -t484 - t503;
t204 = t349 + t386 + t725 * t592 + (t590 * t723 + t592 * t724) * t588;
t257 = (t461 * t590 - t785) * t587 + t750;
t280 = -t367 * t590 + t691;
t282 = t461 * t592 + t349 + (-t501 * t592 - t503 * t590) * t588;
t873 = (t118 * t217 + t147 * t246 + t148 * t247 + t157 * t176 + t177 * t215 + t178 * t216) * t860 + (t175 * t297 + t203 * t320 + t204 * t321 + t242 * t257 + t272 * t280 + t273 * t282) * t861;
t652 = Icges(4,5) * t591 - Icges(4,6) * t589;
t522 = -Icges(4,3) * t592 + t590 * t652;
t667 = rSges(4,1) * t591 - rSges(4,2) * t589;
t528 = -rSges(4,3) * t592 + t590 * t667;
t871 = (t893 * t590 + t894 * t592) * t835 + ((t893 + t903) * t592 + ((t905 * t579 + t904 * t580 - t929) * t592 + t896 * t588 + t895 * t587 - t894) * t590) * t834;
t676 = t898 * t837 + t897 * t838;
t868 = 2 * qJD(2);
t867 = 4 * qJD(2);
t866 = 2 * qJD(3);
t865 = 4 * qJD(3);
t864 = 2 * qJD(4);
t863 = 4 * qJD(4);
t862 = m(4) / 0.2e1;
t478 = t528 * t587;
t479 = t528 * t588;
t759 = t591 * t592;
t550 = t588 * t759 + t771;
t414 = rSges(4,1) * t550 - rSges(4,2) * t632 + rSges(4,3) * t765;
t781 = t414 * t592;
t548 = t587 * t759 - t766;
t413 = rSges(4,1) * t548 - rSges(4,2) * t631 + rSges(4,3) * t770;
t782 = t413 * t592;
t270 = (-t478 * t590 + t782) * t588 + (t479 * t590 - t781) * t587;
t529 = rSges(4,3) * t590 + t592 * t667;
t637 = t528 * t592 + t529 * t590;
t291 = -t413 * t590 - t478 * t592 + t587 * t637;
t292 = t414 * t590 + t479 * t592 - t588 * t637;
t312 = (t413 * t588 - t414 * t587) * t590;
t328 = t528 * t770 + t782;
t329 = -t528 * t765 - t781;
t859 = m(4) * (t270 * t312 + t291 * t328 + t292 * t329);
t857 = m(5) * (t175 * t242 + t203 * t272 + t204 * t273);
t84 = t242 * t310 + t272 * t330 + t273 * t331;
t854 = m(5) * t84;
t853 = m(5) * (t257 * t297 + t280 * t320 + t282 * t321);
t757 = t215 * t764 + t216 * t769;
t851 = m(6) * (-t118 * t592 + (t147 * t588 + t148 * t587 + t176) * t590 + t757);
t850 = m(6) * (t118 * t176 + t147 * t215 + t148 * t216);
t756 = t246 * t764 + t247 * t769;
t847 = m(6) * (-t157 * t592 + (t177 * t588 + t178 * t587 + t217) * t590 + t756);
t846 = m(6) * (t157 * t217 + t177 * t246 + t178 * t247);
t49 = t176 * t276 + t215 * t300 + t216 * t301;
t845 = m(6) * t49;
t844 = m(6) * (t217 * t243 + t246 * t277 + t247 * t278);
t843 = m(6) * (t166 * t230 + t302 * t340 + t304 * t341);
t561 = t874 * t590;
t842 = m(6) * (t176 * t561 + t757);
t841 = m(6) * (t166 * t285 + t302 * t442 + t304 * t443);
t840 = m(6) * (t217 * t561 + t756);
t839 = t561 / 0.2e1;
t836 = t588 / 0.2e1;
t290 = t413 * t587 + t414 * t588 + t711;
t436 = -rSges(4,1) * t631 - rSges(4,2) * t548;
t437 = -rSges(4,1) * t632 - rSges(4,2) * t550;
t322 = t436 * t587 + t437 * t588;
t714 = -t528 - t570;
t432 = t714 * t587;
t434 = t714 * t588;
t560 = (-rSges(4,1) * t589 - rSges(4,2) * t591) * t590;
t833 = m(4) * (t290 * t322 + (-t432 * t587 - t434 * t588) * t560);
t832 = m(5) * (t231 * t289 + t334 * t444 + t336 * t445);
t144 = t297 * t310 + t320 * t330 + t321 * t331;
t143 = m(5) * t144;
t831 = m(5) * (t231 * t315 + (-t334 * t587 - t336 * t588) * t541);
t830 = m(6) * (-t276 * t592 + (t300 * t588 + t301 * t587) * t590);
t829 = m(6) * (-t230 * t592 + (t340 * t587 + t341 * t588) * t590);
t828 = m(6) * (-t285 * t592 + (t442 * t587 + t443 * t588) * t590);
t827 = m(6) * t276;
t820 = m(6) * qJD(2);
t819 = m(6) * qJD(5);
t811 = Icges(4,4) * t548;
t810 = Icges(4,4) * t550;
t30 = 0.2e1 * (t118 / 0.4e1 - t230 / 0.4e1) * m(6) + 0.2e1 * (t175 / 0.4e1 - t289 / 0.4e1) * m(5) + 0.2e1 * (t270 / 0.4e1 - t322 / 0.4e1) * m(4);
t790 = t30 * qJD(1);
t403 = Icges(4,5) * t548 - Icges(4,6) * t631 + Icges(4,3) * t770;
t784 = t403 * t592;
t404 = Icges(4,5) * t550 - Icges(4,6) * t632 + Icges(4,3) * t765;
t783 = t404 * t592;
t778 = t522 * t592;
t69 = 0.2e1 * (t157 / 0.4e1 - t285 / 0.4e1) * m(6) + 0.2e1 * (t257 / 0.4e1 - t315 / 0.4e1) * m(5);
t758 = t69 * qJD(1);
t755 = t302 * t769 + t304 * t764;
t405 = -Icges(4,2) * t631 + Icges(4,6) * t770 + t811;
t733 = -Icges(4,1) * t631 - t405 - t811;
t406 = -Icges(4,2) * t632 + Icges(4,6) * t765 + t810;
t732 = -Icges(4,1) * t632 - t406 - t810;
t543 = Icges(4,4) * t631;
t407 = Icges(4,1) * t548 + Icges(4,5) * t770 - t543;
t731 = -Icges(4,2) * t548 + t407 - t543;
t544 = Icges(4,4) * t632;
t408 = Icges(4,1) * t550 + Icges(4,5) * t765 - t544;
t730 = -Icges(4,2) * t550 + t408 - t544;
t728 = -t412 + t546;
t713 = -t529 - t571;
t712 = t874 * t570;
t709 = t874 * t760;
t707 = qJD(3) * t590;
t499 = (t839 - t590 / 0.2e1) * m(6);
t706 = t499 * qJD(1);
t702 = t830 / 0.2e1;
t698 = t217 * t276 + t246 * t300 + t247 * t301;
t696 = t546 + t739;
t689 = -t571 + t723;
t687 = t770 / 0.2e1;
t685 = t769 / 0.2e1;
t683 = t765 / 0.2e1;
t681 = t764 / 0.2e1;
t673 = t587 * t438 + t588 * t439 - t712;
t671 = -t571 + t692;
t568 = rSges(3,1) * t590 + rSges(3,2) * t592;
t664 = (t926 * t869 + ((t902 * t579 + t901 * t580) * t592 + t882 * t588 + t883 * t587) * t590) * t834 + t889 * t687 + t888 * t683;
t642 = -t405 * t589 + t407 * t591;
t283 = t590 * t642 - t784;
t641 = -t406 * t589 + t408 * t591;
t284 = t590 * t641 - t783;
t647 = t283 * t587 + t284 * t588;
t638 = -t524 * t589 + t526 * t591;
t173 = 0.2e1 * (t243 / 0.4e1 - t276 / 0.4e1) * m(6);
t598 = m(6) * (-t592 * t243 + (t277 * t588 + t278 * t587) * t590);
t68 = t702 - t598 / 0.2e1;
t636 = -t173 * qJD(1) + t68 * qJD(5);
t635 = -t592 * t545 - t585 * t704;
t633 = t143 + t664;
t626 = -t522 * t587 - t642;
t625 = -t522 * t588 - t641;
t426 = -Icges(4,5) * t631 - Icges(4,6) * t548;
t222 = t426 * t770 + t548 * t733 - t631 * t731;
t427 = -Icges(4,5) * t632 - Icges(4,6) * t550;
t223 = t427 * t770 + t548 * t732 - t631 * t730;
t127 = -t222 * t588 + t223 * t587;
t224 = t426 * t765 + t550 * t733 - t632 * t731;
t225 = t427 * t765 + t550 * t732 - t632 * t730;
t128 = -t224 * t588 + t225 * t587;
t624 = t127 * t837 + t128 * t838;
t621 = (Icges(4,3) * t590 + t592 * t652 - t638) * t592;
t619 = t892 * t587 + t552;
t618 = -t892 * t588 + t551;
t600 = t590 * t626 + t784;
t599 = t590 * t625 + t783;
t597 = t884 * t681 + t890 * t683 + t885 * t685 + t891 * t687 + t871;
t596 = t597 + t873;
t595 = -t676 + t890 * t838 + t891 * t837 + (t587 * t917 - t588 * t918) * t835 + (t896 * t587 - t895 * t588) * t834 + t887 * t687 + (t910 * t587 - t911 * t588) * t685 + t886 * t683 + (t908 * t587 - t909 * t588) * t681;
t594 = -t871 + t888 * t838 + t889 * t837 + (t882 * t587 - t883 * t588) * t834 - t891 * t770 / 0.2e1 + t898 * t687 - t885 * t769 / 0.2e1 - t890 * t765 / 0.2e1 + t897 * t683 - t884 * t764 / 0.2e1;
t566 = t585 * t703;
t527 = Icges(4,5) * t590 + t592 * t661;
t525 = Icges(4,6) * t590 + t592 * t656;
t498 = (t835 + t839) * m(6);
t485 = t545 * t765;
t482 = t709 - t760;
t476 = t526 * t588;
t475 = t526 * t587;
t474 = t524 * t588;
t473 = t524 * t587;
t435 = t713 * t588;
t433 = t713 * t587;
t421 = t874 * t568;
t339 = -t437 * t592 - t560 * t765;
t338 = t436 * t592 + t560 * t770;
t337 = t689 * t588;
t335 = t689 * t587;
t327 = t590 * t638 - t778;
t318 = (t436 * t588 - t437 * t587) * t590;
t314 = t522 * t765 - t524 * t632 + t526 * t550;
t313 = t522 * t770 - t524 * t631 + t526 * t548;
t309 = -t478 * t587 - t479 * t588 - t712;
t308 = t592 * t728 + t566 - t701;
t307 = t635 + t330;
t306 = m(5) * t310;
t305 = t671 * t588;
t303 = t671 * t587;
t281 = t728 * t770 + t371 - t485;
t269 = t404 * t765 - t406 * t632 + t408 * t550;
t268 = t403 * t765 - t405 * t632 + t407 * t550;
t267 = t404 * t770 - t406 * t631 + t408 * t548;
t266 = t403 * t770 - t405 * t631 + t407 * t548;
t260 = t592 * t696 + t566 + t668;
t259 = t635 + t277;
t258 = -t459 * t587 - t461 * t588 + t673;
t237 = -t427 * t592 + (-t589 * t730 + t591 * t732) * t590;
t236 = -t426 * t592 + (-t589 * t731 + t591 * t733) * t590;
t232 = t828 / 0.2e1;
t229 = -t625 * t592 + (t474 * t589 - t476 * t591 + t404) * t590;
t228 = -t626 * t592 + (t473 * t589 - t475 * t591 + t403) * t590;
t226 = t696 * t770 - t485 + t749;
t214 = t474 * t632 - t476 * t550 + t588 * t599;
t213 = t473 * t632 - t475 * t550 + t588 * t600;
t212 = t474 * t631 - t476 * t548 + t587 * t599;
t211 = t473 * t631 - t475 * t548 + t587 * t600;
t202 = -t880 * t587 - t738 * t588 + t673;
t171 = t829 / 0.2e1;
t170 = -t327 * t592 + t590 * t647;
t156 = -t314 * t592 + (t268 * t587 + t588 * t269) * t590;
t155 = -t313 * t592 + (t587 * t266 + t267 * t588) * t590;
t149 = t306 + t827 / 0.2e1 + t243 * t860;
t124 = -t213 * t588 + t214 * t587;
t123 = -t211 * t588 + t212 * t587;
t119 = t166 * t561 + t755;
t88 = t840 / 0.2e1;
t87 = -(t550 * t716 - t632 * t715) * t592 + (t224 * t587 + (t225 - t775) * t588) * t590;
t86 = -(t548 * t716 - t631 * t715) * t592 + (t223 * t588 + (t222 - t775) * t587) * t590;
t70 = (t257 + t315) * t861 + (t157 + t285) * t860;
t67 = t702 + t598 / 0.2e1;
t65 = t842 / 0.2e1;
t63 = (t621 + t647) * t592 + (t229 * t588 + t228 * t587 - (-t525 * t589 + t527 * t591 + t522) * t592 + t327) * t590;
t59 = -(-t525 * t632 + t527 * t550) * t592 + t314 * t590 + (t213 * t590 + t268 * t592) * t587 + ((t269 - t778) * t592 + (t214 - t621) * t590) * t588;
t58 = -(-t525 * t631 + t527 * t548) * t592 + t313 * t590 + (t212 * t590 + t267 * t592) * t588 + ((t266 - t778) * t592 + (t211 - t621) * t590) * t587;
t31 = (t270 + t322) * t862 + (t175 + t289) * t861 + (t118 + t230) * t860;
t27 = t847 / 0.2e1;
t18 = t851 / 0.2e1;
t17 = t88 + t232 - t847 / 0.2e1;
t16 = t88 + t27 - t828 / 0.2e1;
t15 = t232 + t27 - t840 / 0.2e1;
t14 = t65 + t171 - t851 / 0.2e1;
t13 = t65 + t18 - t829 / 0.2e1;
t12 = t171 + t18 - t842 / 0.2e1;
t11 = t676 + t831 + t841;
t8 = t624 + t676 + t832 + t833 + t843;
t7 = t633 + t844;
t6 = t664 + t845 + t854;
t5 = t597 + t846 + t853;
t4 = t597 + t859 + t857 + t850 + (-t63 / 0.2e1 + t156 * t836 + t155 * t838) * t592 + (t170 / 0.2e1 + t59 * t836 + t58 * t838) * t590;
t3 = t596 - t875;
t2 = t596 + t875;
t1 = t594 + t821 + t822 - t873;
t9 = [0, t31 * qJD(3) + t70 * qJD(4) + t498 * qJD(5) + (-m(3) * t421 / 0.2e1 + t309 * t862 + t258 * t861 + t202 * t860) * t868, t31 * qJD(2) + t149 * qJD(4) + (t226 * t860 + t281 * t861 + t318 * t862) * t866, t70 * qJD(2) + t149 * qJD(3) + (t306 + t827) * qJD(4), t498 * qJD(2); -qJD(3) * t30 - qJD(4) * t69 + qJD(5) * t499, t8 * qJD(3) + t11 * qJD(4) + t119 * t819 + (m(5) * (t231 * t258 + t334 * t335 + t336 * t337) + m(6) * (t166 * t202 + t302 * t303 + t304 * t305) + m(4) * (t290 * t309 + t432 * t433 + t434 * t435) + m(3) * (-t421 + t568) * t874 * (rSges(3,1) * t592 - rSges(3,2) * t590) + (t124 + t583 * t552 + (t619 * t588 + (-t551 + t618) * t587) * t588 + t886) * t838 + (t123 + t584 * t551 + (t618 * t587 + (-t552 + t619) * t588) * t587 + t887) * t837) * qJD(2), -t790 + t8 * qJD(2) + t1 * qJD(4) + t14 * qJD(5) + (-t859 / 0.4e1 - t857 / 0.4e1 - t850 / 0.4e1) * t865 + ((t290 * t318 + t312 * t322 + t338 * t434 + t339 * t432 + (-t328 * t588 - t329 * t587) * t560) * t862 + (t231 * t281 + t242 * t289 + t272 * t445 + t273 * t444 + t307 * t336 + t308 * t334) * t861 + (t166 * t226 + t176 * t230 + t215 * t341 + t216 * t340 + t259 * t304 + t260 * t302) * t860) * t866 + (-t170 / 0.2e1 + (-t59 / 0.2e1 + t128 / 0.2e1) * t588 + (-t58 / 0.2e1 + t127 / 0.2e1) * t587) * t707 + (t86 * t837 + t87 * t838 + t594 + (t63 / 0.2e1 + (-t156 / 0.2e1 + t236 / 0.2e1) * t588 + (-t155 / 0.2e1 - t237 / 0.2e1) * t587) * t592) * qJD(3), -t758 + t11 * qJD(2) + t1 * qJD(3) + t594 * qJD(4) + t17 * qJD(5) + (-t853 / 0.4e1 - t846 / 0.4e1) * t863 + ((t297 * t315 + (-t320 * t588 - t321 * t587) * t541 + t699) * t861 + (t217 * t285 + t246 * t443 + t247 * t442 + t700) * t860) * t864, t706 + t119 * t820 + t14 * qJD(3) + t17 * qJD(4) + (-t561 * t592 - t482 + t709) * t819; qJD(2) * t30 - qJD(4) * t173, t790 + ((-t266 * t588 + t267 * t587) * t685 + (-t268 * t588 + t269 * t587) * t681 + (-t228 * t588 + t229 * t587) * t834 + (-t283 * t588 + t284 * t587) * t835 + t58 * t837 + t124 * t683 + t123 * t687 + t59 * t838 + t595 - t624) * qJD(2) + t4 * qJD(3) + t2 * qJD(4) + t13 * qJD(5) + (-t833 / 0.4e1 - t832 / 0.4e1 - t843 / 0.4e1) * t867 + ((t270 * t290 + t291 * t434 + t292 * t432 + t309 * t312 + t328 * t435 + t329 * t433) * t862 + (t175 * t231 + t203 * t336 + t204 * t334 + t242 * t258 + t272 * t337 + t273 * t335) * t861 + (t118 * t166 + t147 * t304 + t148 * t302 + t176 * t202 + t215 * t305 + t216 * t303) * t860) * t868, t4 * qJD(2) + (-t869 * t775 / 0.2e1 + t664) * qJD(3) + t6 * qJD(4) + (m(6) * (t176 * t226 + t215 * t259 + t216 * t260) / 0.4e1 + m(5) * (t242 * t281 + t272 * t307 + t273 * t308) / 0.4e1 + m(4) * (t312 * t318 + t328 * t338 + t329 * t339) / 0.4e1) * t865 + ((t237 * t588 + t236 * t587 - (-t715 * t589 + t716 * t591) * t592) * t834 + t87 * t836 + t86 * t838) * t707, t2 * qJD(2) + t6 * qJD(3) + t664 * qJD(4) + (-t844 / 0.4e1 - t143 / 0.4e1) * t863 + ((t698 + t49) * t860 + (t144 + t84) * t861) * t864 + t636, qJD(2) * t13 + qJD(4) * t68; qJD(2) * t69 + qJD(3) * t173, t758 + t595 * qJD(2) + t3 * qJD(3) + t5 * qJD(4) + t16 * qJD(5) + (-t841 / 0.4e1 - t831 / 0.4e1) * t867 + ((t231 * t257 + t258 * t297 + t280 * t336 + t282 * t334 + t320 * t337 + t321 * t335) * t861 + (t157 * t166 + t177 * t304 + t178 * t302 + t202 * t217 + t246 * t305 + t247 * t303) * t860) * t868, t3 * qJD(2) + t664 * qJD(3) + t7 * qJD(4) + (-t845 / 0.4e1 - t854 / 0.4e1) * t865 + ((t176 * t243 + t215 * t277 + t216 * t278 + t217 * t226 + t246 * t259 + t247 * t260) * t860 + (t281 * t297 + t307 * t320 + t308 * t321 + t84) * t861) * t866 - t636, t5 * qJD(2) + t7 * qJD(3) + (m(6) * t698 + t633) * qJD(4), qJD(2) * t16 - qJD(3) * t68; -t499 * qJD(2), -t706 + (-t202 * t592 + (t303 * t587 + t305 * t588 + t166) * t590 - t119 + t755) * t820 + t12 * qJD(3) + t15 * qJD(4) + t482 * t819, t12 * qJD(2) + m(6) * (-t226 * t592 + (t259 * t588 + t260 * t587) * t590) * qJD(3) + t67 * qJD(4), t15 * qJD(2) + t67 * qJD(3) + qJD(4) * t830, t482 * t820;];
Cq = t9;