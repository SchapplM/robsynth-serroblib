% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:01
% EndTime: 2019-12-31 21:11:25
% DurationCPUTime: 18.39s
% Computational Cost: add. (62461->619), mult. (73899->770), div. (0->0), fcn. (80259->8), ass. (0->388)
t760 = Icges(6,1) - Icges(6,2);
t1002 = 2 * Icges(6,4);
t901 = m(6) / 0.2e1;
t588 = qJ(1) + qJ(2);
t584 = sin(t588);
t585 = cos(t588);
t591 = cos(qJ(3));
t589 = sin(qJ(3));
t812 = qJ(4) * t589;
t900 = pkin(3) + pkin(4);
t915 = t900 * t591 + pkin(2) + t812;
t882 = sin(qJ(5));
t883 = cos(qJ(5));
t536 = t589 * t883 - t591 * t882;
t494 = t536 * t585;
t619 = t589 * t882 + t591 * t883;
t495 = t619 * t585;
t387 = t495 * rSges(6,1) + t494 * rSges(6,2) - t584 * rSges(6,3);
t926 = -t584 * pkin(8) + t387;
t311 = t584 * pkin(7) + t915 * t585 + t926;
t840 = cos(qJ(1)) * pkin(1);
t308 = t311 + t840;
t771 = t585 * t589;
t284 = t308 * t771;
t829 = rSges(5,3) + qJ(4);
t884 = rSges(5,1) + pkin(3);
t912 = t829 * t589 + t884 * t591 + pkin(2);
t379 = (rSges(5,2) + pkin(7)) * t584 + t912 * t585;
t366 = t379 + t840;
t347 = t366 * t771;
t576 = t585 * rSges(5,2);
t578 = t585 * pkin(7);
t378 = -t912 * t584 + t576 + t578;
t841 = sin(qJ(1)) * pkin(1);
t365 = t378 - t841;
t695 = t379 * t771;
t696 = t311 * t771;
t779 = t584 * t589;
t902 = m(5) / 0.2e1;
t492 = t536 * t584;
t493 = t619 * t584;
t925 = -t493 * rSges(6,1) - t492 * rSges(6,2);
t954 = t578 + (-rSges(6,3) - pkin(8)) * t585 - t915 * t584 + t925;
t970 = t954 - t841;
t758 = (t284 + t696 + (-t970 - t954) * t779) * t901 + (t347 + t695 + (-t365 - t378) * t779) * t902;
t199 = -t779 * t970 + t284;
t744 = t779 * t954 - t696;
t978 = t365 - t378;
t759 = (t199 + t744) * t901 + (-t978 * t779 + t347 - t695) * t902;
t15 = t759 - t758;
t1001 = t15 * qJD(1);
t722 = t536 * t1002 + t760 * t619;
t724 = t619 * t1002 - t760 * t536;
t613 = Icges(6,5) * t495 + Icges(6,6) * t494 - Icges(6,3) * t584;
t604 = t585 * t613;
t487 = Icges(6,4) * t495;
t382 = Icges(6,2) * t494 - Icges(6,6) * t584 + t487;
t486 = Icges(6,4) * t494;
t384 = Icges(6,1) * t495 - Icges(6,5) * t584 + t486;
t634 = -t492 * t382 - t493 * t384;
t210 = t604 - t634;
t583 = t585 ^ 2;
t485 = Icges(6,4) * t493;
t701 = 0.2e1 * t485;
t815 = Icges(6,6) * t585;
t644 = t701 + 0.2e1 * t815;
t820 = Icges(6,5) * t585;
t700 = 0.2e1 * t820;
t818 = Icges(6,2) * t492;
t825 = Icges(6,1) * t493;
t595 = t492 * (t644 + t818) + t493 * (t700 + t825) + Icges(6,3) * t583;
t931 = t585 * t595;
t109 = t210 * t584 - t931;
t211 = t494 * t382 + t495 * t384 - t584 * t613;
t641 = Icges(6,5) * t493 + Icges(6,6) * t492;
t614 = Icges(6,3) * t585 + t641;
t645 = t818 + t485;
t616 = t815 + t645;
t484 = Icges(6,4) * t492;
t649 = t484 + t825;
t617 = t820 + t649;
t110 = t211 * t584 - (t494 * t616 + t495 * t617 - t584 * t614) * t585;
t927 = t614 - t641;
t40 = t931 + (t927 * t584 + (-t617 + t649) * t495 - (t616 - t645) * t494 + t634) * t584;
t41 = (-t494 * t645 - t495 * t649 + t210 - 0.2e1 * t604 + t634) * t585 + (-t492 * t616 - t493 * t617 - t927 * t585 + t211 + t595) * t584;
t886 = -t585 / 0.2e1;
t888 = -t584 / 0.2e1;
t965 = t585 / 0.2e1;
t914 = t110 * t886 + t41 * t965 + (t109 + t40) * t888;
t733 = -Icges(6,1) * t494 + t382 + t487;
t917 = -Icges(6,2) * t495 + t384 + t486;
t987 = (t536 * t733 + t619 * t917) * t584;
t449 = -Icges(6,5) * t619 - Icges(6,6) * t536;
t996 = t585 * (-t449 * t585 + t492 * t724 + t493 * t722);
t997 = t584 * (t449 * t584 + t494 * t724 + t495 * t722);
t1000 = t914 + t987 / 0.2e1 - t996 / 0.2e1 + t997 / 0.2e1;
t999 = -t987 / 0.4e1 + t996 / 0.4e1 - t997 / 0.4e1;
t396 = -rSges(6,1) * t492 + rSges(6,2) * t493;
t399 = t494 * rSges(6,1) - t495 * rSges(6,2);
t150 = t308 * t399 + t396 * t970;
t151 = t311 * t399 + t396 * t954;
t660 = -t722 * t536 / 0.2e1 + t724 * t619 / 0.2e1;
t653 = (t151 + t150) * t901 + t660;
t811 = qJ(4) * t591;
t550 = pkin(3) * t589 - t811;
t837 = rSges(5,1) * t589;
t551 = -rSges(5,3) * t591 + t837;
t708 = t550 + t551;
t445 = t708 * t584;
t447 = t708 * t585;
t656 = -t445 * t771 + t447 * t779;
t459 = rSges(6,1) * t536 - rSges(6,2) * t619;
t839 = pkin(4) * t589;
t659 = t459 + t550 + t839;
t367 = t659 * t584;
t369 = t659 * t585;
t657 = -t367 * t771 + t369 * t779;
t770 = t585 * t591;
t778 = t584 * t591;
t737 = t365 * t770 + t366 * t778;
t747 = t308 * t778 + t770 * t970;
t755 = (t657 + t747) * t901 + (t656 + t737) * t902;
t566 = pkin(3) * t779;
t428 = t566 + (-t829 * t591 + t837) * t584;
t558 = qJ(4) * t770;
t565 = rSges(5,3) * t770;
t429 = -t884 * t771 + t558 + t565;
t729 = t428 * t771 + t429 * t779;
t341 = t566 + (-t811 + t839) * t584 - t396;
t342 = -t900 * t771 - t399 + t558;
t739 = t341 * t771 + t342 * t779;
t757 = (t739 + t747) * t901 + (t729 + t737) * t902;
t947 = t755 - t757;
t735 = t378 * t770 + t379 * t778;
t743 = t311 * t778 + t770 * t954;
t754 = (t657 + t743) * t901 + (t656 + t735) * t902;
t756 = (t739 + t743) * t901 + (t729 + t735) * t902;
t948 = t754 - t756;
t995 = -t947 * qJD(1) - t948 * qJD(2);
t586 = Icges(5,5) * t589;
t642 = Icges(5,3) * t591 - t586;
t822 = Icges(4,4) * t589;
t983 = Icges(4,2) * t591 + t642 + t822;
t587 = Icges(4,4) * t591;
t548 = Icges(4,1) * t589 + t587;
t821 = Icges(5,5) * t591;
t993 = Icges(5,1) * t589 + t548 - t821;
t670 = t760 * t493;
t622 = t670 + t820;
t941 = -0.2e1 * t484;
t976 = t622 - t941;
t991 = t976 * t494 + t584 * (-Icges(6,5) * t492 + Icges(6,6) * t493);
t385 = rSges(6,3) * t585 - t925;
t276 = -t584 * t385 - t585 * t387;
t969 = -t396 * t584 + t399 * t585;
t988 = t276 * t969;
t390 = -Icges(6,5) * t494 + Icges(6,6) * t495;
t986 = (-t390 * t584 - t494 * t917 + t495 * t733) * t584;
t985 = (t390 * t585 - t492 * t917 + t493 * t733) * t584;
t984 = t969 * t591;
t545 = -Icges(4,2) * t589 + t587;
t982 = t545 + t993;
t549 = Icges(4,1) * t591 - t822;
t922 = Icges(5,1) * t591 + t586;
t981 = t549 + t922;
t541 = Icges(5,3) * t589 + t821;
t601 = (-t541 / 0.2e1 + t982 / 0.2e1) * t591 + (-t983 / 0.2e1 + t981 / 0.2e1) * t589 - t660;
t117 = t308 * t954 - t311 * t970;
t458 = -rSges(6,1) * t619 - rSges(6,2) * t536;
t582 = t584 ^ 2;
t705 = t582 + t583;
t929 = t705 * t589;
t842 = m(6) * (t458 * t929 + t984);
t975 = (-Icges(4,6) + Icges(5,6)) * t591 + (-Icges(5,4) - Icges(4,5)) * t589;
t473 = Icges(5,4) * t584 + t585 * t922;
t475 = Icges(4,5) * t584 + t549 * t585;
t974 = -t983 * t585 + t473 + t475;
t472 = -Icges(5,4) * t585 + t584 * t922;
t562 = Icges(4,4) * t779;
t474 = Icges(4,1) * t778 - Icges(4,5) * t585 - t562;
t973 = -Icges(4,2) * t778 - t642 * t584 + t472 + t474 - t562;
t561 = Icges(5,5) * t770;
t465 = Icges(5,6) * t584 + Icges(5,3) * t771 + t561;
t471 = Icges(4,6) * t584 + t545 * t585;
t972 = -Icges(5,1) * t771 - t548 * t585 + t465 - t471 + t561;
t464 = -Icges(5,6) * t585 + t541 * t584;
t470 = Icges(4,4) * t778 - Icges(4,2) * t779 - Icges(4,6) * t585;
t971 = t993 * t584 - t464 + t470;
t903 = m(4) / 0.2e1;
t887 = t584 / 0.2e1;
t881 = m(3) * (t840 * (-rSges(3,1) * t584 - rSges(3,2) * t585) + (t585 * rSges(3,1) - t584 * rSges(3,2)) * t841);
t543 = Icges(5,4) * t591 + Icges(5,6) * t589;
t788 = t543 * t584;
t468 = -Icges(5,2) * t585 + t788;
t443 = t584 * t468;
t316 = t464 * t771 + t472 * t770 + t443;
t961 = t316 * t585;
t937 = (t311 * t584 + t585 * t954) * t458;
t936 = (t308 * t584 + t585 * t970) * t458;
t935 = t459 * t705;
t953 = (t981 - t983) * t591 + (t541 - t982) * t589;
t669 = t760 * t492;
t950 = t669 - t815;
t603 = t701 - t950;
t808 = (t603 * t536 + t619 * t976) * t585;
t661 = t808 / 0.4e1 + t999;
t665 = -0.2e1 * t485;
t602 = t665 + t950;
t810 = (t602 * t536 - (t622 + 0.2e1 * t484) * t619) * t585;
t662 = -t810 / 0.4e1 + t999;
t951 = t661 - t662;
t933 = (t464 * t589 + t472 * t591) * t584;
t765 = t589 * t591;
t709 = t705 * t765;
t930 = (m(5) / 0.4e1 + m(6) / 0.4e1) * (t709 - t765);
t196 = -t379 * t365 + t366 * t378;
t838 = rSges(4,1) * t591;
t673 = pkin(2) + t838;
t706 = rSges(4,2) * t779 + t585 * rSges(4,3);
t426 = -t584 * t673 + t578 + t706;
t419 = t426 - t841;
t564 = rSges(4,2) * t771;
t427 = -t564 + t673 * t585 + (rSges(4,3) + pkin(7)) * t584;
t420 = t427 + t840;
t245 = -t427 * t419 + t420 * t426;
t924 = t975 * t585;
t923 = t975 * t584;
t921 = t973 * t589 + t971 * t591;
t920 = -t974 * t589 + t972 * t591;
t703 = qJD(1) + qJD(2);
t625 = (-t986 - (t602 * t495 + t991) * t585) * t887 + (-t985 - ((t665 - 0.2e1 * t815) * t493 - (-0.2e1 * t670 + t941 - 0.2e1 * t820) * t492) * t585) * t886;
t626 = (t986 - (t603 * t495 - t991) * t585) * t887 + (t985 - (-(-t941 + t700) * t492 + (t644 - 0.2e1 * t669) * t493) * t585) * t886;
t552 = rSges(4,1) * t589 + rSges(4,2) * t591;
t698 = ((-t366 + t379) * t447 + t978 * t445) * t902 + ((-t420 + t427) * t585 + (t419 - t426) * t584) * t552 * t903 + ((-t308 + t311) * t369 + (-t954 + t970) * t367) * t901;
t131 = t308 * t342 + t341 * t970;
t144 = t311 * t342 + t341 * t954;
t233 = t365 * t428 + t366 * t429;
t237 = t378 * t428 + t379 * t429;
t508 = t552 * t584;
t510 = t552 * t585;
t275 = t419 * t508 - t420 * t510;
t281 = t426 * t508 - t427 * t510;
t913 = (t237 + t233) * t902 + (t281 + t275) * t903 + (t144 + t131) * t901;
t908 = 0.4e1 * qJD(1);
t906 = 0.4e1 * qJD(2);
t905 = 2 * qJD(3);
t799 = t459 * t585;
t263 = t308 * t799;
t800 = t459 * t584;
t748 = -t311 * t799 + t800 * t954;
t894 = m(6) * (-t800 * t970 + t263 + t748);
t83 = (t311 * t585 - t584 * t954) * t459 + t748;
t893 = m(6) * t83;
t746 = t341 * t799 + t342 * t800;
t892 = m(6) * (t746 - t936);
t745 = -t367 * t399 - t369 * t396;
t891 = m(6) * (t745 - t936);
t890 = m(6) * (t746 - t937);
t889 = m(6) * (t745 - t937);
t877 = m(4) * t245;
t875 = m(4) * t275;
t874 = m(4) * t281;
t868 = m(5) * t196;
t862 = m(5) * t233;
t860 = m(5) * t237;
t859 = m(5) * (-t365 * t779 + t347);
t858 = m(5) * (-t378 * t779 + t695);
t855 = m(6) * t117;
t667 = t705 * t458;
t727 = t591 * t935;
t852 = m(6) * (-t984 + (t276 - t667) * t589 + t727);
t850 = m(6) * t131;
t848 = m(6) * t144;
t845 = m(6) * t199;
t844 = m(6) * t744;
t843 = m(6) * (t276 * t929 + t727);
t687 = t40 / 0.2e1 + t109 / 0.2e1;
t688 = t110 / 0.2e1 - t41 / 0.2e1;
t13 = t584 * t687 + t585 * t688;
t624 = (t396 * t585 + t399 * t584) * t589;
t192 = -t624 * m(6) / 0.2e1;
t704 = t192 * qJD(4);
t828 = t13 * qJD(5) + t704;
t797 = t470 * t589;
t542 = Icges(4,5) * t591 - Icges(4,6) * t589;
t789 = t542 * t585;
t787 = t543 * t585;
t773 = t585 * t468;
t749 = -t369 * t341 - t367 * t342;
t738 = -t447 * t428 - t445 * t429;
t736 = -t367 * t778 - t369 * t770;
t728 = -t445 * t778 - t447 * t770;
t466 = Icges(4,5) * t778 - Icges(4,6) * t779 - Icges(4,3) * t585;
t726 = -t584 * t466 - t474 * t770;
t467 = Icges(4,3) * t584 + t789;
t725 = t584 * t467 + t475 * t770;
t712 = t584 * (qJ(4) * t778 - t566) + t585 * (-pkin(3) * t771 + t558);
t554 = pkin(3) * t591 + t812;
t711 = t705 * t554;
t555 = rSges(5,1) * t591 + rSges(5,3) * t589;
t707 = -t554 - t555;
t469 = Icges(5,2) * t584 + t787;
t317 = t465 * t771 + t584 * t469 + t473 * t770;
t671 = -t626 - t625;
t432 = t475 * t778;
t664 = t585 * t467 - t432;
t663 = t471 * t589 - t466;
t658 = -pkin(4) * t591 + t458 - t554;
t218 = (pkin(4) * t770 + t926) * t585 + (pkin(4) * t778 + pkin(8) * t585 + t385) * t584 + t711;
t294 = t584 * (t555 * t584 - t576) + (t584 * rSges(5,2) + t555 * t585) * t585 + t711;
t79 = m(5) * (t294 * t929 + t728) + m(6) * (t218 * t929 + t736);
t655 = t396 * t799 + t399 * t800;
t654 = -t465 * t779 + t585 * t469 - t473 * t778;
t636 = t367 * t584 + t369 * t585;
t368 = t658 * t584;
t370 = t658 * t585;
t635 = t368 * t584 + t370 * t585;
t628 = t317 + t773;
t623 = (-t508 * t585 + t510 * t584) * t552;
t621 = t308 * t585;
t618 = t810 / 0.2e1 + t1000;
t599 = t601 + t913;
t598 = t914 + t951;
t597 = t914 - t951;
t596 = t661 + t662 - t914;
t594 = ((t464 + t470) * t591 + (-t472 + t474) * t589) * (t887 + t888) - t601;
t312 = -t773 + t933;
t221 = -t312 * t585 - t584 * t654;
t315 = -t471 * t779 - t664;
t222 = -(-t584 * (-t474 * t591 + t797) - t585 * t466) * t585 + t315 * t584;
t223 = t317 * t584 - t961;
t318 = -t470 * t771 - t726;
t319 = -t471 * t771 + t725;
t224 = -t318 * t585 + t319 * t584;
t73 = (t317 - t628) * t585 + (t316 + t654 - t443) * t584;
t74 = (t585 * t663 + t319 - t725) * t585 + (t584 * t663 + t318 + t664) * t584;
t75 = -t961 + (t312 + t628 - t933) * t584;
t76 = (t315 - t432 + (t467 + t797) * t585 + t726) * t585 + t725 * t584;
t593 = (-t808 / 0.2e1 + (t76 + t75) * t965 + (t222 + t221 + t74 + t73) * t888 + (t542 * t584 + t953 * t585 + t972 * t589 + t974 * t591 + t788) * t887 + (t953 * t584 - t971 * t589 + t973 * t591 + t223 + t224 - t787 - t789) * t886 + t1000) * qJD(3);
t556 = -rSges(4,2) * t589 + t838;
t448 = t707 * t585;
t446 = t707 * t584;
t355 = 0.4e1 * t930;
t322 = t585 * (-rSges(5,1) * t771 + t565) - t551 * t582 + t712;
t242 = -pkin(4) * t929 + t712 - t969;
t234 = -t842 / 0.2e1;
t206 = t843 / 0.2e1;
t193 = t624 * t901;
t191 = t193 * qJD(5);
t189 = t192 * qJD(5);
t142 = -t844 + t858;
t133 = -t459 * t667 + t988;
t127 = t845 + t859;
t121 = t852 / 0.2e1;
t94 = m(6) * t151 + t660;
t92 = m(6) * t150 + t660;
t90 = t889 / 0.2e1;
t89 = t890 / 0.2e1;
t86 = t891 / 0.2e1;
t85 = t892 / 0.2e1;
t82 = t893 / 0.2e1;
t80 = t894 / 0.2e1;
t52 = t206 + t121 + t842 / 0.2e1;
t51 = t234 + t206 - t852 / 0.2e1;
t50 = t234 + t121 - t843 / 0.2e1;
t44 = t855 + t868 + t877 + t881;
t43 = t601 + t848 + t860 + t874;
t42 = t601 + t850 + t862 + t875;
t26 = t754 + t756;
t23 = t755 + t757;
t21 = m(6) * t133 + t626;
t20 = t82 - t894 / 0.2e1 + t653;
t19 = t80 - t893 / 0.2e1 + t653;
t18 = m(6) * (t218 * t969 + t458 * t636) + t625;
t16 = t758 + t759;
t14 = t80 + t82 - t653;
t11 = t599 + t698;
t10 = t599 - t698;
t9 = t594 + t698 - t913;
t8 = (t224 / 0.2e1 - t75 / 0.2e1 - t76 / 0.2e1 + t223 / 0.2e1 + t688) * t585 + (t74 / 0.2e1 + t221 / 0.2e1 + t222 / 0.2e1 + t73 / 0.2e1 + t687) * t584;
t7 = t8 * qJD(3);
t6 = t89 + t598 - t889 / 0.2e1;
t5 = t89 + t596 + t90;
t4 = t597 + t90 - t890 / 0.2e1;
t3 = t598 + t85 - t891 / 0.2e1;
t2 = t86 + t596 + t85;
t1 = t86 + t597 - t892 / 0.2e1;
t12 = [t44 * qJD(2) + t42 * qJD(3) + t127 * qJD(4) + t92 * qJD(5), t44 * qJD(1) + t11 * qJD(3) + t16 * qJD(4) + t19 * qJD(5) + 0.2e1 * (t117 * t901 + t196 * t902 + t245 * t903 + t881 / 0.2e1) * qJD(2), t42 * qJD(1) + t11 * qJD(2) + t593 + t23 * qJD(4) + t2 * qJD(5) + (((-t419 * t585 - t420 * t584) * t556 + t623) * t903 + (t365 * t448 + t366 * t446 + t738) * t902 + (t308 * t368 + t370 * t970 + t749) * t901) * t905, qJD(1) * t127 + qJD(2) * t16 + qJD(3) * t23 - t189, t92 * qJD(1) + t19 * qJD(2) + t2 * qJD(3) - t704 + ((t655 + t936) * m(6) + t618) * qJD(5); t10 * qJD(3) - t15 * qJD(4) + t20 * qJD(5) + (-t855 / 0.4e1 - t868 / 0.4e1 - t877 / 0.4e1 - t881 / 0.4e1) * t908, qJD(3) * t43 + qJD(4) * t142 + qJD(5) * t94, t10 * qJD(1) + t43 * qJD(2) + t593 + t26 * qJD(4) + t5 * qJD(5) + (((-t426 * t585 - t427 * t584) * t556 + t623) * t903 + (t378 * t448 + t379 * t446 + t738) * t902 + (t311 * t368 + t370 * t954 + t749) * t901) * t905, qJD(2) * t142 + qJD(3) * t26 - t1001 - t189, t20 * qJD(1) + t94 * qJD(2) + t5 * qJD(3) - t704 + ((t655 + t937) * m(6) + t618) * qJD(5); t594 * qJD(1) + t9 * qJD(2) + t7 + t947 * qJD(4) + t1 * qJD(5) + (-t875 / 0.4e1 - t862 / 0.4e1 - t850 / 0.4e1) * t908, t9 * qJD(1) + t594 * qJD(2) + t7 + t948 * qJD(4) + t4 * qJD(5) + (-t874 / 0.4e1 - t860 / 0.4e1 - t848 / 0.4e1) * t906, (m(6) * (t218 * t242 - t367 * t368 - t369 * t370) + m(5) * (t294 * t322 - t445 * t446 - t447 * t448) + m(4) * ((t584 * (rSges(4,1) * t778 - t706) + t585 * (rSges(4,1) * t770 + t584 * rSges(4,3) - t564)) * (-t508 * t584 - t510 * t585) + t705 * t556 * t552) + t626 + (t924 * t582 + (t921 * t585 + (t920 - t923) * t584) * t585) * t887 + (t923 * t583 + (t920 * t584 + (t921 - t924) * t585) * t584) * t886) * qJD(3) + t79 * qJD(4) + t18 * qJD(5) + t703 * t8, t79 * qJD(3) + (-0.4e1 * t930 + 0.2e1 * (t901 + t902) * (-t591 * t929 + t709)) * qJD(4) + t51 * qJD(5) - t995, t1 * qJD(1) + t4 * qJD(2) + t18 * qJD(3) + t51 * qJD(4) + ((-t133 + (-t218 + t276) * t969 + (-t636 - t935) * t458) * m(6) - t625 + t671) * qJD(5); t15 * qJD(2) - t947 * qJD(3) + t191 + (-t845 / 0.4e1 - t859 / 0.4e1) * t908 + (-t589 * t621 + t284) * m(6) * qJD(1), t1001 - t948 * qJD(3) + t191 + (t844 / 0.4e1 - t858 / 0.4e1) * t906, (m(6) * (-t242 * t591 + t736) + m(5) * (-t322 * t591 + t728) + 0.2e1 * ((t218 + t635) * t901 + (t446 * t584 + t448 * t585 + t294) * t902) * t589 - t79) * qJD(3) + t355 * qJD(4) + t50 * qJD(5) + t995, t355 * qJD(3), t50 * qJD(3) + qJD(5) * t842 + t703 * t193; ((-t459 * t621 - t150 + t263) * m(6) - t660) * qJD(1) + t14 * qJD(2) + t3 * qJD(3) + t828, t14 * qJD(1) + ((t83 - t151) * m(6) - t660) * qJD(2) + t6 * qJD(3) + t828, t3 * qJD(1) + t6 * qJD(2) + ((t242 * t276 + t459 * t635) * m(6) - t626 + t671) * qJD(3) + t52 * qJD(4) + t21 * qJD(5), qJD(3) * t52 + t192 * t703, t21 * qJD(3) + (m(6) * (t458 * t935 - t988) + t625) * qJD(5) + t703 * t13;];
Cq = t12;
