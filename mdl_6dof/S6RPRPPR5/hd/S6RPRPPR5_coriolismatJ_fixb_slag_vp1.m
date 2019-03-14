% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPPR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:07
% EndTime: 2019-03-09 02:50:48
% DurationCPUTime: 33.21s
% Computational Cost: add. (62787->914), mult. (62661->1271), div. (0->0), fcn. (65414->9), ass. (0->525)
t597 = pkin(9) + qJ(3);
t589 = sin(t597);
t604 = sin(qJ(1));
t801 = t589 * t604;
t578 = pkin(3) * t801;
t591 = cos(t597);
t842 = rSges(5,3) + qJ(4);
t850 = rSges(5,2) * t589;
t417 = t578 + (-t591 * t842 - t850) * t604;
t605 = cos(qJ(1));
t795 = t591 * t605;
t563 = qJ(4) * t795;
t800 = t589 * t605;
t679 = -pkin(3) * t800 + t563;
t738 = rSges(5,2) * t800 + rSges(5,3) * t795;
t418 = t679 + t738;
t600 = sin(pkin(10));
t601 = cos(pkin(10));
t674 = rSges(6,1) * t600 + rSges(6,2) * t601;
t841 = rSges(6,3) + qJ(5);
t342 = t578 + (t841 * t589 + (-qJ(4) - t674) * t591) * t604;
t726 = pkin(3) + t841;
t794 = t600 * t605;
t716 = t591 * t794;
t793 = t601 * t605;
t742 = t591 * rSges(6,2) * t793 + rSges(6,1) * t716;
t343 = -t726 * t800 + t563 + t742;
t651 = t342 * t605 + t343 * t604;
t603 = -pkin(8) - qJ(5);
t596 = pkin(10) + qJ(6);
t588 = sin(t596);
t590 = cos(t596);
t671 = rSges(7,1) * t588 + rSges(7,2) * t590;
t856 = pkin(5) * t600;
t702 = qJ(4) + t856;
t318 = t578 + ((rSges(7,3) - t603) * t589 + (-t671 - t702) * t591) * t604;
t792 = t603 * t605;
t741 = pkin(5) * t716 + t589 * t792;
t745 = t671 * t795;
t908 = rSges(7,3) + pkin(3);
t319 = -t800 * t908 + t563 + t741 + t745;
t652 = t318 * t605 + t319 * t604;
t572 = rSges(5,3) * t801;
t584 = cos(pkin(9)) * pkin(2) + pkin(1);
t804 = t589 * qJ(4);
t692 = -t584 - t804;
t854 = pkin(7) + qJ(2);
t358 = -t572 + (rSges(5,1) + t854) * t605 + ((rSges(5,2) - pkin(3)) * t591 + t692) * t604;
t585 = t604 * t854;
t700 = t604 * rSges(5,1) - rSges(5,2) * t795;
t857 = pkin(3) * t591;
t359 = t585 + (t589 * t842 + t584 + t857) * t605 + t700;
t796 = t591 * t604;
t770 = t358 * t795 + t359 * t796;
t789 = t604 * t600;
t522 = t589 * t793 - t789;
t788 = t604 * t601;
t523 = t589 * t794 + t788;
t676 = t523 * rSges(6,1) + t522 * rSges(6,2);
t855 = t604 * pkin(4);
t286 = t855 + t585 + (t591 * t726 - t692) * t605 + t676;
t562 = qJ(5) * t796;
t524 = t589 * t788 + t794;
t715 = t589 * t789;
t525 = -t715 + t793;
t960 = t525 * rSges(6,1) - t524 * rSges(6,2);
t984 = (pkin(4) + t854) * t605 + ((-rSges(6,3) - pkin(3)) * t591 + t692) * t604 - t562 + t960;
t774 = t286 * t796 + t795 * t984;
t583 = pkin(5) * t601 + pkin(4);
t791 = t604 * t588;
t797 = t590 * t605;
t491 = t589 * t797 - t791;
t790 = t604 * t590;
t492 = t588 * t800 + t790;
t673 = t492 * rSges(7,1) + t491 * rSges(7,2);
t639 = -t591 * t792 + t673;
t950 = t589 * t702 + t591 * t908 + t584;
t250 = t604 * t583 + t605 * t950 + t585 + t639;
t704 = t605 * t854;
t740 = t605 * t583 + t603 * t796;
t493 = t588 * t605 + t589 * t790;
t494 = -t589 * t791 + t797;
t961 = t494 * rSges(7,1) - t493 * rSges(7,2);
t985 = -t604 * t950 + t704 + t740 + t961;
t777 = t250 * t796 + t795 * t985;
t940 = m(7) / 0.2e1;
t941 = m(6) / 0.2e1;
t942 = m(5) / 0.2e1;
t714 = (t589 * t652 + t777) * t940 + (t589 * t651 + t774) * t941 + ((t417 * t605 + t418 * t604) * t589 + t770) * t942;
t545 = pkin(3) * t589 - qJ(4) * t591;
t824 = qJ(5) * t589;
t691 = t545 + t824;
t724 = t591 * t856;
t784 = -qJ(5) - t603;
t954 = -rSges(7,3) * t589 + t591 * t671;
t640 = t589 * t784 + t691 - t724 - t954;
t298 = t640 * t604;
t300 = t640 * t605;
t953 = -rSges(6,3) * t589 + t591 * t674;
t678 = t691 - t953;
t362 = t678 * t604;
t364 = t678 * t605;
t670 = rSges(5,3) * t591 + t850;
t744 = t545 - t670;
t420 = t744 * t604;
t423 = t744 * t605;
t718 = (-t362 * t800 + t364 * t801 + t774) * t941 + (-t420 * t800 + t423 * t801 + t770) * t942 + (-t298 * t800 + t300 * t801 + t777) * t940;
t9 = t718 - t714;
t1016 = t9 * qJD(1);
t963 = -t286 * t604 - t605 * t984;
t628 = t963 * t589;
t656 = -t250 * t604 - t605 * t985;
t629 = t656 * t589;
t779 = (t591 * t652 + t629) * t940 + (t591 * t651 + t628) * t941;
t840 = (-t362 * t795 + t364 * t796 + t628) * t941 + (-t298 * t795 + t300 * t796 + t629) * t940;
t12 = t840 - t779;
t1015 = t12 * qJD(1);
t1014 = Icges(5,1) + Icges(4,3);
t1013 = Icges(5,4) - Icges(4,5);
t1012 = Icges(5,5) - Icges(4,6);
t661 = Icges(7,5) * t588 + Icges(7,6) * t590;
t617 = -Icges(7,3) * t589 + t591 * t661;
t835 = Icges(7,4) * t588;
t664 = Icges(7,2) * t590 + t835;
t619 = -Icges(7,6) * t589 + t591 * t664;
t834 = Icges(7,4) * t590;
t667 = Icges(7,1) * t588 + t834;
t621 = -Icges(7,5) * t589 + t591 * t667;
t229 = -t493 * t619 + t494 * t621 - t617 * t796;
t349 = Icges(7,5) * t492 + Icges(7,6) * t491 + Icges(7,3) * t795;
t836 = Icges(7,4) * t492;
t352 = Icges(7,2) * t491 + Icges(7,6) * t795 + t836;
t477 = Icges(7,4) * t491;
t355 = Icges(7,1) * t492 + Icges(7,5) * t795 + t477;
t176 = t349 * t796 + t493 * t352 - t494 * t355;
t351 = -Icges(7,5) * t494 + Icges(7,6) * t493 + Icges(7,3) * t796;
t479 = Icges(7,4) * t494;
t354 = Icges(7,2) * t493 + Icges(7,6) * t796 - t479;
t478 = Icges(7,4) * t493;
t356 = Icges(7,1) * t494 - Icges(7,5) * t796 - t478;
t177 = t351 * t796 + t354 * t493 + t356 * t494;
t658 = t176 * t605 + t177 * t604;
t70 = t229 * t589 + t591 * t658;
t537 = Icges(4,5) * t591 - Icges(4,6) * t589;
t538 = -Icges(5,4) * t591 + Icges(5,5) * t589;
t1010 = (t537 + t538) * t605 + t1014 * t604;
t1002 = -t1012 * t801 + t1013 * t796 + t1014 * t605;
t227 = -t491 * t619 - t492 * t621 - t617 * t795;
t174 = t349 * t795 + t491 * t352 + t492 * t355;
t175 = t351 * t795 + t491 * t354 - t492 * t356;
t659 = t174 * t605 + t604 * t175;
t1009 = t227 * t589 + t591 * t659;
t817 = t351 * t589;
t997 = t354 * t590 - t356 * t588;
t196 = t591 * t997 - t817;
t361 = rSges(7,3) * t796 - t961;
t1004 = -t361 * t589 - t954 * t796;
t624 = rSges(7,3) * t795 + t673;
t276 = -t589 * t624 - t795 * t954;
t655 = -t1004 * t605 + t276 * t604;
t547 = rSges(4,1) * t589 + rSges(4,2) * t591;
t519 = t547 * t604;
t521 = t547 * t605;
t392 = t604 * t519 + t521 * t605;
t682 = (t604 * t318 - t319 * t605) * t940 + (t604 * t342 - t343 * t605) * t941 + (t604 * t417 - t418 * t605) * t942 + m(4) * t392 / 0.2e1;
t648 = t362 * t604 + t364 * t605;
t654 = t298 * t604 + t300 * t605;
t598 = t604 ^ 2;
t599 = t605 ^ 2;
t736 = t598 + t599;
t966 = t736 * t547;
t996 = -m(7) / 0.2e1;
t683 = t654 * t996 - m(6) * t648 / 0.2e1 + (-t420 * t604 - t423 * t605) * t942 - m(4) * t966 / 0.2e1;
t36 = t683 - t682;
t1008 = t36 * qJD(1);
t837 = Icges(4,4) * t589;
t542 = Icges(4,1) * t591 - t837;
t467 = Icges(4,5) * t604 + t542 * t605;
t829 = Icges(5,6) * t591;
t534 = Icges(5,3) * t589 - t829;
t468 = Icges(5,5) * t604 + t534 * t605;
t1006 = -t467 * t796 - t468 * t801;
t1005 = t174 * t604 - t175 * t605;
t975 = m(6) + m(7);
t1003 = t975 / 0.2e1;
t415 = t954 * t604;
t1001 = t1012 * t591 + t1013 * t589;
t1000 = t1010 * t605 + t1006;
t981 = -t1010 * t604 - t467 * t795 - t468 * t800;
t570 = Icges(4,4) * t801;
t466 = Icges(4,1) * t796 - Icges(4,5) * t605 - t570;
t469 = Icges(5,5) * t605 + Icges(5,6) * t796 - Icges(5,3) * t801;
t999 = t1002 * t604 - t466 * t795 + t469 * t800;
t464 = Icges(4,4) * t796 - Icges(4,2) * t801 - Icges(4,6) * t605;
t565 = Icges(5,6) * t801;
t471 = Icges(5,4) * t605 + Icges(5,2) * t796 - t565;
t998 = t464 * t589 - t471 * t591;
t917 = t589 / 0.2e1;
t582 = Icges(4,4) * t591;
t831 = Icges(4,2) * t589;
t465 = Icges(4,6) * t604 + (t582 - t831) * t605;
t566 = Icges(5,6) * t800;
t470 = Icges(5,4) * t604 - Icges(5,2) * t795 + t566;
t986 = -t465 * t800 - t470 * t795 - t981;
t983 = t465 * t589 + t470 * t591 + t1002;
t979 = -t464 * t800 - t465 * t801 - t470 * t796 + t471 * t795 - t1000 - t999;
t374 = -Icges(6,4) * t525 + Icges(6,2) * t524 + Icges(6,6) * t796;
t377 = -Icges(6,1) * t525 + Icges(6,4) * t524 + Icges(6,5) * t796;
t647 = t374 * t524 - t377 * t525;
t964 = t736 * t591;
t419 = (t964 - t591) * t589;
t978 = -0.4e1 * t419;
t965 = t736 * t589;
t977 = -0.2e1 * t965;
t728 = m(6) / 0.4e1 + m(7) / 0.4e1;
t976 = -0.2e1 * t728;
t916 = t591 / 0.2e1;
t913 = t604 / 0.2e1;
t911 = -t605 / 0.2e1;
t909 = t605 / 0.2e1;
t438 = Icges(7,3) * t591 + t589 * t661;
t539 = Icges(4,2) * t591 + t837;
t665 = Icges(6,4) * t600 + Icges(6,2) * t601;
t620 = -Icges(6,6) * t589 + t591 * t665;
t668 = Icges(6,1) * t600 + Icges(6,4) * t601;
t622 = -Icges(6,5) * t589 + t591 * t668;
t662 = Icges(6,5) * t600 + Icges(6,6) * t601;
t798 = t590 * t619;
t805 = t588 * t621;
t914 = t601 / 0.2e1;
t915 = t600 / 0.2e1;
t969 = t589 * (-t622 * t915 - t620 * t914 - t805 / 0.2e1 - t798 / 0.2e1 + t542 / 0.2e1 - t539 / 0.2e1 - Icges(5,6) * t589 - Icges(5,3) * t591 / 0.2e1 + t662 * t917 + t438 / 0.2e1 + (Icges(5,2) + Icges(6,3)) * t916);
t968 = t591 * t662;
t389 = rSges(7,1) * t491 - rSges(7,2) * t492;
t390 = rSges(7,1) * t493 + rSges(7,2) * t494;
t646 = t389 * t604 - t390 * t605;
t259 = t646 * t591;
t846 = rSges(7,3) * t591;
t444 = t589 * t671 + t846;
t725 = t589 * t856;
t962 = -t591 * t784 - t444 - t725;
t959 = t1001 * t604;
t958 = t1001 * t605;
t952 = -t605 * t361 + t604 * t624;
t475 = (Icges(7,2) * t588 - t834) * t591;
t476 = (-Icges(7,1) * t590 + t835) * t591;
t951 = -t588 * (t476 / 0.2e1 + t619 / 0.2e1) - t590 * (-t621 / 0.2e1 + t475 / 0.2e1);
t369 = Icges(6,5) * t523 + Icges(6,6) * t522 + Icges(6,3) * t795;
t371 = -Icges(6,5) * t525 + Icges(6,6) * t524 + Icges(6,3) * t796;
t826 = Icges(6,3) * t589;
t618 = -t826 + t968;
t425 = t618 * t604;
t426 = t618 * t605;
t507 = Icges(5,3) * t796 + t565;
t508 = Icges(5,3) * t795 + t566;
t660 = Icges(5,2) * t589 + t829;
t509 = t660 * t604;
t510 = t660 * t605;
t515 = -Icges(4,2) * t796 - t570;
t516 = t539 * t605;
t838 = Icges(4,1) * t589;
t669 = -t582 - t838;
t517 = t669 * t604;
t518 = t669 * t605;
t812 = t377 * t600;
t375 = Icges(6,1) * t523 + Icges(6,4) * t522 + Icges(6,5) * t795;
t813 = t375 * t600;
t814 = t374 * t601;
t372 = Icges(6,4) * t523 + Icges(6,2) * t522 + Icges(6,6) * t795;
t815 = t372 * t601;
t949 = ((-t425 - t812 - t814 + t469 + t509 + t464 - t517) * t605 + (t426 + t813 + t815 + t468 - t510 - t465 + t518) * t604) * t591 + ((t371 + t471 - t507 + t466 + t515) * t605 + (-t369 + t470 + t508 - t467 + t516) * t604) * t589;
t948 = 0.4e1 * t419;
t946 = 0.2e1 * qJD(1);
t945 = 0.4e1 * qJD(1);
t944 = 2 * qJD(3);
t939 = t70 / 0.2e1;
t216 = (t444 * t604 - t361) * t591;
t416 = -rSges(7,3) * t800 + t745;
t217 = ((-t444 + t846) * t605 + t673) * t591 + (-t605 * t954 + t416) * t589;
t936 = m(7) * (t1004 * t318 + t216 * t985 + t217 * t250 - t276 * t319);
t180 = (t415 * t605 - t416 * t604) * t591 + t952 * t589;
t238 = t952 * t591;
t635 = t216 * t605 + t217 * t604 - t238;
t775 = t1004 * t795 - t276 * t796;
t935 = m(7) * (-t180 * t591 + t589 * t635 + t775);
t934 = m(7) * (t635 * t591 + (t180 + t655) * t589);
t263 = t276 * t800;
t626 = t276 * t605;
t932 = m(7) * (t589 * t626 - t263);
t264 = t276 * t795;
t931 = m(7) * (t591 * t626 - t264);
t480 = (-rSges(7,1) * t590 + rSges(7,2) * t588) * t591;
t930 = m(7) * (-t298 * t389 + t300 * t390 + t480 * t656);
t532 = pkin(4) * t605 - t562;
t548 = t804 + t857;
t747 = t736 * t548;
t681 = -t604 * t532 + t605 * (qJ(5) * t795 + t855) + t747;
t823 = qJ(5) * t591;
t155 = (pkin(5) * t715 + t361 + t532 - t740) * t604 + ((t725 - t823 + t846) * t605 + (-pkin(4) + t583) * t604 + t639) * t605 + t681;
t772 = -t298 * t796 - t300 * t795;
t929 = m(7) * (t155 * t965 + t772);
t383 = Icges(7,5) * t491 - Icges(7,6) * t492;
t766 = -Icges(7,2) * t492 + t355 + t477;
t768 = -Icges(7,1) * t491 + t352 + t836;
t153 = t383 * t589 + (t588 * t768 - t590 * t766) * t591;
t924 = t153 / 0.2e1;
t456 = Icges(6,5) * t591 + t589 * t668;
t920 = t456 / 0.2e1;
t912 = t604 / 0.4e1;
t910 = -t605 / 0.4e1;
t907 = m(3) * t736 * (rSges(3,3) + qJ(2));
t852 = rSges(4,1) * t591;
t701 = t584 + t852;
t739 = rSges(4,2) * t801 + t605 * rSges(4,3);
t405 = -t604 * t701 + t704 + t739;
t699 = -rSges(4,2) * t800 + t604 * rSges(4,3);
t406 = t605 * t701 + t585 + t699;
t906 = m(4) * (t405 * t519 - t406 * t521);
t905 = m(4) * (t405 * t605 + t406 * t604);
t279 = -t604 * (rSges(5,1) * t605 + rSges(5,2) * t796 - t572) + t605 * (rSges(5,3) * t800 + t700) + t747;
t760 = -t420 * t796 - t423 * t795;
t899 = m(5) * (t279 * t965 + t760);
t897 = m(5) * (t358 * t417 + t359 * t418);
t896 = m(5) * (-t358 * t801 + t359 * t800);
t895 = m(5) * (t358 * t605 + t359 * t604);
t212 = t604 * (rSges(6,3) * t796 - t960) + t605 * (rSges(6,3) * t795 + t676) + t681;
t769 = -t362 * t796 - t364 * t795;
t887 = m(6) * (t212 * t965 + t769);
t884 = m(6) * (t286 * t343 + t342 * t984);
t268 = t286 * t800;
t883 = m(6) * (-t801 * t984 + t268);
t269 = t286 * t795;
t882 = m(6) * (-t796 * t984 + t269);
t881 = m(6) * t963;
t874 = m(7) * (t250 * t319 + t318 * t985);
t873 = m(7) * (-t238 * t964 + t589 * t655);
t872 = m(7) * (-t238 * t965 + t775);
t235 = t250 * t800;
t871 = m(7) * (-t801 * t985 + t235);
t236 = t250 * t795;
t870 = m(7) * (-t796 * t985 + t236);
t869 = m(7) * t656;
t868 = m(7) * (-t1004 * t801 - t263);
t867 = m(7) * (-t1004 * t796 - t264);
t866 = m(7) * t655;
t272 = t389 * t605 + t604 * t390;
t863 = m(7) * (-t272 * t591 - t480 * t965);
t862 = m(7) * (t272 * t589 - t480 * t964);
t860 = t646 * m(7) * t589;
t859 = m(7) * t259;
t858 = m(7) * t272;
t853 = m(7) * qJD(6);
t845 = t604 * t70;
t844 = t605 * t1009;
t821 = t250 * t605;
t820 = t286 * t605;
t818 = t349 * t589;
t811 = t617 * t589;
t442 = Icges(7,5) * t591 + t589 * t667;
t806 = t588 * t442;
t474 = (-Icges(7,5) * t590 + Icges(7,6) * t588) * t591;
t802 = t589 * t474;
t440 = Icges(7,6) * t591 + t589 * t664;
t799 = t590 * t440;
t767 = -Icges(7,1) * t493 + t354 - t479;
t765 = Icges(7,2) * t494 - t356 + t478;
t586 = t589 ^ 2;
t587 = t591 ^ 2;
t737 = t736 * t587;
t393 = -t587 + (0.1e1 - t736) * t586 + t737;
t763 = t393 * t1003;
t403 = -t586 * t736 - t591 * t964;
t762 = t403 * t1003;
t404 = t589 * t965 + t737;
t761 = t404 * t1003;
t755 = -t619 - t476;
t754 = -t621 + t475;
t748 = t604 * (qJ(4) * t796 - t578) + t605 * t679;
t743 = rSges(5,2) * t591 - rSges(5,3) * t589 - t548;
t735 = qJD(1) * t591;
t734 = qJD(3) * t589;
t733 = qJD(3) * t591;
t732 = qJD(6) * t591;
t121 = t216 * t604 - t217 * t605;
t731 = t121 * qJD(2);
t727 = t940 + t941;
t241 = (t942 + t727) * t977;
t730 = t241 * qJD(1);
t306 = -0.2e1 * t727 * t964;
t729 = t306 * qJD(1);
t721 = t121 * t940;
t720 = -t70 / 0.2e1 + t939;
t199 = t369 * t795 + t522 * t372 + t523 * t375;
t710 = t796 / 0.4e1;
t706 = t538 / 0.2e1 + t537 / 0.2e1;
t690 = -t548 - t823;
t685 = m(5) / 0.4e1 + t728;
t56 = m(6) * (t212 * t964 + t589 * t648) + m(7) * (t155 * t964 + t589 * t654);
t384 = Icges(7,5) * t493 + Icges(7,6) * t494;
t154 = t384 * t589 + (t588 * t767 - t590 * t765) * t591;
t171 = t474 * t795 + t491 * t754 - t492 * t755;
t172 = t474 * t796 + t493 * t754 + t494 * t755;
t677 = t930 / 0.2e1 + (t153 + t171) * t912 + (t154 + t172) * t910;
t650 = -t352 * t590 - t355 * t588;
t195 = t591 * t650 + t818;
t657 = t195 * t605 - t196 * t604;
t303 = -t390 * t589 + t480 * t796;
t304 = t589 * t389 - t480 * t795;
t653 = t303 * t605 + t304 * t604;
t645 = t798 + t805;
t641 = -m(7) * (t250 * t389 - t390 * t985) - t802 / 0.2e1;
t126 = t383 * t795 + t491 * t766 - t492 * t768;
t127 = t384 * t795 + t491 * t765 - t492 * t767;
t63 = t126 * t604 - t127 * t605;
t128 = t383 * t796 + t493 * t766 + t494 * t768;
t129 = t384 * t796 + t493 * t765 + t494 * t767;
t64 = t128 * t604 - t129 * t605;
t638 = t63 * t913 + t64 * t911;
t634 = t617 * t605 - t650;
t633 = t617 * t604 + t997;
t631 = t438 - t645;
t623 = -t70 * t912 + t1009 * t910 + t845 / 0.4e1 + t844 / 0.4e1 + (t710 - t796 / 0.4e1) * t1005;
t616 = t591 * t634 - t818;
t615 = t591 * t633 - t817;
t614 = t591 * t631 + t811;
t410 = t619 * t605;
t412 = t621 * t605;
t139 = (-t410 * t590 - t412 * t588 + t349) * t591 + t634 * t589;
t409 = t619 * t604;
t411 = t621 * t604;
t140 = (-t409 * t590 - t411 * t588 + t351) * t591 + t633 * t589;
t162 = t440 * t493 - t442 * t494 + t604 * t614;
t163 = t491 * t440 + t492 * t442 + t605 * t614;
t168 = (-t617 - t799 - t806) * t591 + t631 * t589;
t254 = t591 * t645 - t811;
t613 = t168 * t917 + t254 * t916 + t936 / 0.2e1 - (-t196 + t229) * t801 / 0.4e1 - (t195 + t227) * t800 / 0.4e1 + (t140 + t162) * t710 + (t139 + t163) * t795 / 0.4e1;
t612 = -(t199 - t981) * t604 / 0.2e1 + (t199 + t986) * t913 + ((t998 + t1010) * t605 + t979 + t999 + t1006) * t911;
t611 = (t371 * t796 + t647 + t1002 * t605 + (t466 * t591 - t469 * t589 - t998) * t604) * t911 + (t605 * t983 + t647 + t981 + t986) * t909 + (t371 * t795 + t604 * t983 + t1000 + t979) * t913;
t299 = -t562 + (-t548 + t962) * t604;
t301 = (t690 + t962) * t605;
t458 = rSges(6,3) * t591 + t589 * t674;
t363 = -t562 + (-t458 - t548) * t604;
t365 = (-t458 + t690) * t605;
t609 = (t299 * t604 + t301 * t605 + t155) * t940 + (t363 * t604 + t365 * t605 + t212) * t941;
t454 = Icges(6,6) * t591 + t589 * t665;
t607 = t456 * t915 + t454 * t914 + t617 / 0.2e1 - t582 - t838 / 0.2e1 + t831 / 0.2e1 - t660 / 0.2e1 + t534 / 0.2e1 + t968 / 0.2e1 - t826 / 0.2e1 + t806 / 0.2e1 + t799 / 0.2e1;
t550 = -rSges(4,2) * t589 + t852;
t430 = t622 * t605;
t429 = t622 * t604;
t428 = t620 * t605;
t427 = t620 * t604;
t424 = t743 * t605;
t421 = t743 * t604;
t308 = t728 * t978;
t305 = (t976 + t1003) * t964;
t297 = t598 * t670 + t605 * t738 + t748;
t265 = -t858 / 0.2e1;
t257 = t859 / 0.2e1;
t256 = t860 / 0.2e1;
t247 = t685 * t948;
t240 = t685 * t977 + (m(5) + t975) * t965 / 0.2e1;
t230 = t605 * (-rSges(6,3) * t800 + t742) - qJ(5) * t965 + t953 * t598 + t748;
t219 = t862 / 0.2e1;
t218 = t863 / 0.2e1;
t211 = (t802 + (t588 * t755 - t590 * t754) * t591) * t589;
t187 = -t866 / 0.2e1;
t185 = (t416 + t741) * t605 + (-t736 + t599) * t824 + (t415 + (t589 * t603 + t724 + t824) * t604) * t604 + t748;
t179 = t867 / 0.2e1;
t178 = t868 / 0.2e1;
t161 = t403 * t976 + t761 + t763;
t160 = t404 * t976 + t762 + t763;
t159 = t393 * t976 + t761 + t762;
t123 = t872 / 0.2e1;
t122 = t873 / 0.2e1;
t120 = t491 * t409 + t492 * t411 + t605 * t615;
t119 = t491 * t410 + t492 * t412 + t605 * t616;
t118 = t409 * t493 - t411 * t494 + t604 * t615;
t117 = t410 * t493 - t412 * t494 + t604 * t616;
t116 = qJD(3) * t721;
t91 = t870 + t882;
t84 = t591 * t951 - t641;
t83 = t254 * t589 + t591 * t657;
t74 = t931 / 0.2e1;
t73 = t932 / 0.2e1;
t71 = t871 + t883 + t896;
t66 = t155 * t272 + t480 * t654;
t62 = -t869 - t881 + t895 + t905 + t907;
t61 = t187 + t858 / 0.2e1;
t60 = t265 + t187;
t59 = t265 + t866 / 0.2e1;
t58 = t119 * t604 - t120 * t605;
t57 = t117 * t604 - t118 * t605;
t55 = t179 + t74 - t859 / 0.2e1;
t54 = t178 + t73 - t860 / 0.2e1;
t53 = t257 + t179 - t931 / 0.2e1;
t52 = t257 + t74 - t867 / 0.2e1;
t51 = t256 + t178 - t932 / 0.2e1;
t50 = t256 + t73 - t868 / 0.2e1;
t49 = t1004 * t216 - t180 * t238 - t217 * t276;
t45 = t172 * t589 + (t128 * t605 + t129 * t604) * t591;
t44 = t171 * t589 + (t126 * t605 + t127 * t604) * t591;
t43 = t934 / 0.2e1;
t42 = t935 / 0.2e1;
t39 = t887 + t899 + t929;
t34 = t682 + t683;
t29 = (t139 * t605 + t140 * t604 + t254) * t591 + (t168 - t657) * t589;
t28 = t123 + t42 - t863 / 0.2e1;
t27 = t122 + t43 - t862 / 0.2e1;
t26 = t219 + t122 - t934 / 0.2e1;
t25 = t219 + t43 - t873 / 0.2e1;
t24 = t218 + t123 - t935 / 0.2e1;
t23 = t218 + t42 - t872 / 0.2e1;
t16 = -t591 * t607 + t874 + t884 + t897 + t906 + t969;
t15 = (t119 * t605 + t120 * t604 + t227) * t591 + (t163 - t659) * t589;
t14 = (t117 * t605 + t118 * t604 + t229) * t591 + (t162 - t658) * t589;
t13 = m(7) * t66 + t638;
t11 = t779 + t840;
t8 = t714 + t718;
t6 = t720 * t795;
t5 = m(7) * t49 + (t15 * t909 + t14 * t913 + t83 / 0.2e1) * t591 + (-t844 / 0.2e1 - t845 / 0.2e1 + t29 / 0.2e1) * t589;
t4 = t604 * t611 + t605 * t612;
t3 = (-t168 / 0.2e1 + (t227 / 0.4e1 + t195 / 0.4e1) * t605 + (-t196 / 0.4e1 + t229 / 0.4e1) * t604) * t589 + t623 - t936 / 0.2e1 + (-t254 / 0.2e1 + (-t163 / 0.4e1 - t139 / 0.4e1) * t605 + (-t162 / 0.4e1 - t140 / 0.4e1) * t604) * t591 + t677;
t2 = t613 + (-t153 / 0.4e1 - t171 / 0.4e1) * t604 + t623 - t930 / 0.2e1 + (t172 / 0.4e1 + t154 / 0.4e1) * t605;
t1 = t613 + t677;
t7 = [t62 * qJD(2) + t16 * qJD(3) + t71 * qJD(4) + t91 * qJD(5) + t84 * qJD(6), qJD(1) * t62 + qJD(3) * t34 + qJD(4) * t240 + qJD(5) * t305 + qJD(6) * t60, t16 * qJD(1) + t34 * qJD(2) + t8 * qJD(4) + t11 * qJD(5) + t1 * qJD(6) + ((t286 * t363 - t342 * t364 - t343 * t362 + t365 * t984) * t941 + (t250 * t299 - t298 * t319 - t300 * t318 + t301 * t985) * t940 + (t358 * t424 + t359 * t421 - t417 * t423 - t418 * t420) * t942) * t944 + ((-t454 * t524 / 0.2e1 + t525 * t920 - t162 / 0.2e1 + m(4) * (-t405 * t550 - t519 * t547) - t140 / 0.2e1 + t706 * t605 - t612) * qJD(3) + (t427 * t914 + t429 * t915 - t371 / 0.2e1 - t471 / 0.2e1 + t507 / 0.2e1 - t466 / 0.2e1 - t515 / 0.2e1) * t733 + (-t814 / 0.2e1 - t812 / 0.2e1 - t425 / 0.2e1 + t469 / 0.2e1 + t509 / 0.2e1 + t464 / 0.2e1 - t517 / 0.2e1) * t734) * t605 + ((t163 / 0.2e1 + t139 / 0.2e1 + m(4) * (-t406 * t550 + t521 * t547) + t522 * t454 / 0.2e1 + t523 * t920 + t706 * t604 - t611) * qJD(3) + (-t428 * t601 / 0.2e1 - t430 * t600 / 0.2e1 + t369 / 0.2e1 - t470 / 0.2e1 - t508 / 0.2e1 + t467 / 0.2e1 - t516 / 0.2e1) * t733 + (t815 / 0.2e1 + t813 / 0.2e1 + t426 / 0.2e1 + t468 / 0.2e1 - t510 / 0.2e1 - t465 / 0.2e1 + t518 / 0.2e1) * t734) * t604, qJD(1) * t71 + qJD(2) * t240 + qJD(3) * t8 + qJD(6) * t51, qJD(1) * t91 + qJD(2) * t305 + qJD(3) * t11 + qJD(6) * t53, t84 * qJD(1) + t60 * qJD(2) + t1 * qJD(3) + t51 * qJD(4) + t53 * qJD(5) + (t211 + m(7) * (-t1004 * t390 + t250 * t304 - t276 * t389 + t303 * t985)) * qJD(6) + ((t924 + t171 / 0.2e1 - t720) * t605 + (t154 / 0.2e1 + t172 / 0.2e1) * t604) * t732; -t36 * qJD(3) + t241 * qJD(4) + t306 * qJD(5) + t59 * qJD(6) + (t869 / 0.4e1 - t895 / 0.4e1 + t881 / 0.4e1 - t905 / 0.4e1 - t907 / 0.4e1) * t945, 0, -t1008 + ((-t421 * t605 + t424 * t604) * t942 + (-t363 * t605 + t365 * t604) * t941 + (-t299 * t605 + t301 * t604) * t940) * t944 + qJD(6) * t721, t730, t729, t59 * qJD(1) + t116 + (t303 * t604 - t304 * t605) * t853; t36 * qJD(2) + t4 * qJD(3) + t9 * qJD(4) + t12 * qJD(5) + t3 * qJD(6) + (-t906 / 0.4e1 - t897 / 0.4e1 - t884 / 0.4e1 - t874 / 0.4e1) * t945 + t607 * t735 - t969 * qJD(1), t1008 - t121 * t853 / 0.2e1, t4 * qJD(1) + t39 * qJD(4) + t56 * qJD(5) + t13 * qJD(6) + (m(7) * (t155 * t185 - t298 * t299 - t300 * t301) + m(5) * (t279 * t297 - t420 * t421 - t423 * t424) + m(6) * (t212 * t230 - t362 * t363 - t364 * t365) + m(4) * (t550 * t966 - (t604 * (rSges(4,1) * t796 - t739) + t605 * (rSges(4,1) * t795 + t699)) * t392) + (t58 + (t522 * t428 + t523 * t430) * t604 + (-t522 * t427 - t523 * t429 - t604 * t959 + t949) * t605 + t958 * t598) * t913 + (t57 - (t427 * t524 - t429 * t525) * t605 + t959 * t599 + (t428 * t524 - t430 * t525 - t605 * t958 + t949) * t604) * t911) * qJD(3), qJD(4) * t685 * t978 + t39 * qJD(3) + t159 * qJD(5) + t24 * qJD(6) + t1016, t728 * t948 * qJD(5) + t56 * qJD(3) + t159 * qJD(4) + t26 * qJD(6) + t1015, t3 * qJD(1) + t13 * qJD(3) + t24 * qJD(4) + t26 * qJD(5) + t731 * t996 + (-t83 / 0.2e1 + (t63 / 0.2e1 - t15 / 0.2e1) * t605 + (t64 / 0.2e1 - t14 / 0.2e1) * t604) * t732 + (t44 * t913 + t45 * t911 + (-t29 / 0.2e1 + (-t154 / 0.2e1 + t1009 / 0.2e1) * t605 + (t924 + t939) * t604) * t589 + (-t259 * t155 - t238 * t272 - t304 * t298 - t303 * t300 + t480 * t655 - t49) * m(7)) * qJD(6); -t241 * qJD(2) - t9 * qJD(3) + t50 * qJD(6) + (-t871 / 0.4e1 - t896 / 0.4e1 - t883 / 0.4e1) * t945 + ((-t589 * t821 + t235) * t940 + (-t589 * t820 + t268) * t941) * t946, -t730, -t1016 + t247 * qJD(4) + t160 * qJD(5) + t23 * qJD(6) + 0.4e1 * (-t929 / 0.4e1 - t899 / 0.4e1 - t887 / 0.4e1) * qJD(3) + ((-t591 * t185 + t772) * t940 + (-t591 * t297 + t760) * t942 + (-t591 * t230 + t769) * t941 + ((t421 * t604 + t424 * t605 + t279) * t942 + t609) * t589) * t944, t247 * qJD(3), t160 * qJD(3), t50 * qJD(1) + t23 * qJD(3) + (t259 * t591 + t589 * t653) * t853; -t306 * qJD(2) - t12 * qJD(3) + t52 * qJD(6) + (-t870 / 0.4e1 - t882 / 0.4e1) * t945 + (t236 * t940 + t269 * t941 + (-t820 * t941 - t821 * t940) * t591) * t946, -t729, -t1015 + (0.2e1 * t609 * t591 + 0.2e1 * ((t185 + t654) * t940 + (t230 + t648) * t941) * t589 - t56) * qJD(3) + t161 * qJD(4) + t308 * qJD(5) + t25 * qJD(6), t161 * qJD(3), t308 * qJD(3), t52 * qJD(1) + t25 * qJD(3) + (-t259 * t589 + t591 * t653) * t853; t641 * qJD(1) + t61 * qJD(2) + t2 * qJD(3) + t54 * qJD(4) + t55 * qJD(5) + t6 * qJD(6) - t735 * t951, qJD(1) * t61 + t116, t2 * qJD(1) + (-(t176 * t604 - t177 * t605) * t801 / 0.2e1 + t57 * t796 / 0.2e1 - t1005 * t800 / 0.2e1 + t58 * t795 / 0.2e1 + t14 * t911 + t15 * t913 + (t195 * t604 + t196 * t605) * t916 + (t139 * t604 - t140 * t605) * t917 - t638) * qJD(3) + t28 * qJD(4) + t27 * qJD(5) + t5 * qJD(6) + (t731 / 0.2e1 + (t1004 * t301 + t155 * t180 - t185 * t238 - t216 * t300 - t217 * t298 - t276 * t299 - t66) * qJD(3)) * m(7), qJD(1) * t54 + qJD(3) * t28, qJD(1) * t55 + qJD(3) * t27, t6 * qJD(1) + t5 * qJD(3) + (t211 * t917 + m(7) * (t1004 * t303 + t238 * t259 - t276 * t304) + (t44 * t909 + t45 * t913 + (t153 * t605 + t154 * t604) * t917) * t591) * qJD(6);];
Cq  = t7;