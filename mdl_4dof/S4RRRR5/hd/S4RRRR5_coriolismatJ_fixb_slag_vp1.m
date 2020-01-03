% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:24
% EndTime: 2019-12-31 17:28:03
% DurationCPUTime: 32.42s
% Computational Cost: add. (83680->905), mult. (127501->1255), div. (0->0), fcn. (141365->8), ass. (0->533)
t606 = sin(qJ(2));
t607 = sin(qJ(1));
t809 = t606 * t607;
t605 = sin(qJ(3));
t608 = cos(qJ(3));
t610 = cos(qJ(1));
t799 = t608 * t610;
t609 = cos(qJ(2));
t802 = t607 * t609;
t550 = t605 * t802 + t799;
t783 = t610 * t605;
t551 = t608 * t802 - t783;
t445 = -Icges(4,5) * t550 - Icges(4,6) * t551;
t543 = Icges(4,4) * t550;
t431 = -Icges(4,1) * t551 - Icges(4,5) * t809 + t543;
t754 = -Icges(4,2) * t551 - t431 - t543;
t544 = Icges(4,4) * t551;
t427 = -Icges(4,2) * t550 + Icges(4,6) * t809 + t544;
t756 = -Icges(4,1) * t550 - t427 - t544;
t221 = t445 * t809 - t550 * t754 + t551 * t756;
t552 = t607 * t608 - t609 * t783;
t791 = t609 * t610;
t812 = t605 * t607;
t553 = t608 * t791 + t812;
t446 = Icges(4,5) * t552 - Icges(4,6) * t553;
t545 = Icges(4,4) * t552;
t808 = t606 * t610;
t432 = Icges(4,1) * t553 + Icges(4,5) * t808 + t545;
t753 = -Icges(4,2) * t553 + t432 + t545;
t848 = Icges(4,4) * t553;
t429 = Icges(4,2) * t552 + Icges(4,6) * t808 + t848;
t755 = Icges(4,1) * t552 - t429 - t848;
t222 = t446 * t809 - t550 * t753 + t551 * t755;
t132 = -t221 * t610 + t222 * t607;
t223 = t445 * t808 + t552 * t754 + t553 * t756;
t224 = t446 * t808 + t552 * t753 + t553 * t755;
t133 = -t223 * t610 + t224 * t607;
t863 = t609 * pkin(2);
t864 = t606 * pkin(6);
t580 = t863 + t864;
t915 = t610 ^ 2;
t916 = t607 ^ 2;
t967 = t915 + t916;
t736 = t967 * t580;
t604 = qJ(3) + qJ(4);
t596 = cos(t604);
t595 = sin(t604);
t784 = t610 * t595;
t531 = t596 * t607 - t609 * t784;
t532 = t595 * t607 + t596 * t791;
t680 = rSges(5,1) * t532 + rSges(5,2) * t531;
t405 = rSges(5,3) * t808 + t680;
t611 = -pkin(7) - pkin(6);
t581 = t611 * t808;
t594 = pkin(3) * t608 + pkin(2);
t862 = pkin(2) - t594;
t694 = t862 * t609;
t726 = pkin(3) * t812;
t437 = t726 - t581 + (-t694 - t864) * t610;
t758 = t405 + t437;
t529 = t595 * t802 + t596 * t610;
t530 = t596 * t802 - t784;
t927 = -t530 * rSges(5,1) + t529 * rSges(5,2);
t403 = rSges(5,3) * t809 - t927;
t792 = t609 * t594;
t687 = pkin(3) * t783 - t607 * t792;
t861 = pkin(6) + t611;
t930 = t861 * t606;
t436 = (t863 + t930) * t607 + t687;
t759 = t403 - t436;
t230 = t607 * t759 + t610 * t758 + t736;
t422 = -rSges(5,1) * t529 - rSges(5,2) * t530;
t423 = t531 * rSges(5,1) - t532 * rSges(5,2);
t340 = t607 * t422 + t610 * t423;
t546 = t550 * pkin(3);
t547 = t552 * pkin(3);
t310 = -t546 * t607 + t547 * t610 + t340;
t926 = -t551 * rSges(4,1) + t550 * rSges(4,2);
t433 = rSges(4,3) * t809 - t926;
t683 = rSges(4,1) * t553 + rSges(4,2) * t552;
t435 = rSges(4,3) * t808 + t683;
t313 = t433 * t607 + t435 * t610 + t736;
t451 = -rSges(4,1) * t550 - rSges(4,2) * t551;
t452 = rSges(4,1) * t552 - rSges(4,2) * t553;
t351 = t451 * t607 + t452 * t610;
t579 = t606 * pkin(2) - t609 * pkin(6);
t928 = t862 * t606;
t929 = t861 * t609;
t493 = -t928 + t929;
t858 = rSges(5,1) * t596;
t679 = -rSges(5,2) * t595 + t858;
t503 = -rSges(5,3) * t609 + t606 * t679;
t748 = t493 + t503;
t718 = t579 + t748;
t368 = t718 * t607;
t370 = t718 * t610;
t859 = rSges(4,1) * t608;
t682 = -rSges(4,2) * t605 + t859;
t641 = t682 * t606;
t946 = rSges(4,3) * t609 - t641;
t739 = -t946 + t579;
t441 = t739 * t607;
t443 = t739 * t610;
t536 = (-rSges(5,1) * t595 - rSges(5,2) * t596) * t606;
t865 = pkin(3) * t605;
t688 = t606 * t865 - t536;
t460 = t688 * t607;
t461 = t688 * t610;
t563 = (-rSges(4,1) * t605 - rSges(4,2) * t608) * t606;
t883 = -t610 / 0.2e1;
t886 = t607 / 0.2e1;
t984 = t132 * t883 + t133 * t886 + m(5) * (t230 * t310 - t368 * t460 - t370 * t461) + m(4) * (t313 * t351 + (t441 * t607 + t443 * t610) * t563);
t711 = -t809 / 0.2e1;
t670 = Icges(5,5) * t596 - Icges(5,6) * t595;
t497 = -Icges(5,3) * t609 + t606 * t670;
t843 = Icges(5,4) * t596;
t673 = -Icges(5,2) * t595 + t843;
t499 = -Icges(5,6) * t609 + t606 * t673;
t844 = Icges(5,4) * t595;
t675 = Icges(5,1) * t596 - t844;
t501 = -Icges(5,5) * t609 + t606 * t675;
t323 = t497 * t808 + t499 * t531 + t501 * t532;
t394 = Icges(5,5) * t530 - Icges(5,6) * t529 + Icges(5,3) * t809;
t508 = Icges(5,4) * t530;
t397 = -Icges(5,2) * t529 + Icges(5,6) * t809 + t508;
t507 = Icges(5,4) * t529;
t401 = -Icges(5,1) * t530 - Icges(5,5) * t809 + t507;
t253 = t394 * t808 + t531 * t397 - t401 * t532;
t396 = Icges(5,5) * t532 + Icges(5,6) * t531 + Icges(5,3) * t808;
t845 = Icges(5,4) * t532;
t399 = Icges(5,2) * t531 + Icges(5,6) * t808 + t845;
t509 = Icges(5,4) * t531;
t402 = Icges(5,1) * t532 + Icges(5,5) * t808 + t509;
t254 = t396 * t808 + t531 * t399 + t532 * t402;
t668 = t253 * t607 + t254 * t610;
t976 = -t323 * t609 + t606 * t668;
t983 = t976 * t711;
t708 = t809 / 0.4e1;
t710 = -t809 / 0.4e1;
t982 = t710 + t708;
t671 = Icges(4,5) * t608 - Icges(4,6) * t605;
t514 = -Icges(4,3) * t609 + t606 * t671;
t846 = Icges(4,4) * t608;
t674 = -Icges(4,2) * t605 + t846;
t518 = -Icges(4,6) * t609 + t606 * t674;
t847 = Icges(4,4) * t605;
t676 = Icges(4,1) * t608 - t847;
t522 = -Icges(4,5) * t609 + t606 * t676;
t331 = t514 * t809 - t518 * t550 + t522 * t551;
t424 = Icges(4,5) * t551 - Icges(4,6) * t550 + Icges(4,3) * t809;
t279 = t424 * t809 - t427 * t550 - t431 * t551;
t426 = Icges(4,5) * t553 + Icges(4,6) * t552 + Icges(4,3) * t808;
t280 = t426 * t809 - t550 * t429 + t551 * t432;
t666 = t279 * t607 + t280 * t610;
t67 = t609 * t331 - t606 * t666;
t321 = t497 * t809 - t499 * t529 + t501 * t530;
t251 = t394 * t809 - t397 * t529 - t401 * t530;
t252 = t396 * t809 - t529 * t399 + t530 * t402;
t669 = t251 * t607 + t252 * t610;
t134 = -t321 * t609 + t606 * t669;
t602 = t610 * pkin(5);
t851 = -rSges(5,3) + t611;
t344 = t602 + (t851 * t606 - pkin(1)) * t607 + t687 + t927;
t856 = rSges(5,3) * t606;
t345 = -t581 + (pkin(5) + t865) * t607 + (pkin(1) + t792 + t856) * t610 + t680;
t981 = m(5) * (-t344 * t422 + t345 * t423);
t471 = t503 * t809;
t346 = t609 * t403 + t471;
t470 = t493 * t809;
t824 = t436 * t609;
t287 = t346 + t470 - t824;
t772 = t759 * t609 - t287 + t470 + t471;
t980 = m(5) * t772;
t333 = t514 * t808 + t518 * t552 + t522 * t553;
t281 = t424 * t808 + t552 * t427 - t431 * t553;
t282 = t426 * t808 + t552 * t429 + t553 * t432;
t665 = t281 * t607 + t282 * t610;
t977 = -t609 * t333 + t606 * t665;
t826 = t424 * t609;
t965 = t427 * t605 + t431 * t608;
t305 = t965 * t606 + t826;
t829 = t394 * t609;
t966 = t397 * t595 + t401 * t596;
t277 = t966 * t606 + t829;
t972 = -t253 * t610 + t254 * t607;
t975 = t982 * t972;
t486 = t946 * t607;
t971 = -t281 * t610 + t282 * t607;
t359 = t433 * t609 - t946 * t809;
t908 = m(4) / 0.2e1;
t969 = -m(5) / 0.2e1;
t906 = m(5) / 0.2e1;
t701 = t802 / 0.4e1;
t745 = t501 + (-Icges(5,2) * t596 - t844) * t606;
t968 = t595 * t745;
t887 = t606 / 0.2e1;
t885 = t607 / 0.4e1;
t884 = -t609 / 0.2e1;
t882 = -t610 / 0.4e1;
t464 = t499 * t607;
t466 = t501 * t607;
t649 = -t497 * t607 + t966;
t206 = -t649 * t609 + (t464 * t595 - t466 * t596 + t394) * t606;
t500 = Icges(5,6) * t606 + t609 * t673;
t502 = Icges(5,5) * t606 + t609 * t675;
t498 = Icges(5,3) * t606 + t609 * t670;
t817 = t596 * t501;
t820 = t595 * t499;
t658 = t817 - t820;
t645 = t498 - t658;
t823 = t497 * t609;
t620 = t606 * t645 + t823;
t242 = -t500 * t529 + t502 * t530 + t607 * t620;
t952 = t206 + t242;
t465 = t499 * t610;
t467 = t501 * t610;
t662 = -t399 * t595 + t402 * t596;
t648 = -t497 * t610 - t662;
t207 = -t648 * t609 + (t465 * t595 - t467 * t596 + t396) * t606;
t243 = t500 * t531 + t502 * t532 + t610 * t620;
t951 = t207 + t243;
t949 = -t277 + t321;
t828 = t396 * t609;
t278 = t606 * t662 - t828;
t948 = t278 + t323;
t520 = Icges(3,4) * t802 - Icges(3,2) * t809 - Icges(3,6) * t610;
t600 = Icges(3,4) * t609;
t840 = Icges(3,2) * t606;
t521 = Icges(3,6) * t607 + (t600 - t840) * t610;
t849 = Icges(3,4) * t606;
t573 = Icges(3,1) * t609 - t849;
t525 = Icges(3,5) * t607 + t573 * t610;
t488 = t525 * t802;
t569 = Icges(3,5) * t609 - Icges(3,6) * t606;
t517 = Icges(3,3) * t607 + t569 * t610;
t691 = t517 * t610 - t488;
t516 = Icges(3,5) * t802 - Icges(3,6) * t809 - Icges(3,3) * t610;
t585 = Icges(3,4) * t809;
t524 = Icges(3,1) * t802 - Icges(3,5) * t610 - t585;
t750 = -t607 * t516 - t524 * t791;
t947 = -t520 * t808 - t521 * t809 - t691 - t750;
t880 = rSges(4,3) + pkin(6);
t923 = t880 * t606 + pkin(1) + t863;
t373 = -t923 * t607 + t602 + t926;
t374 = pkin(5) * t607 + t923 * t610 + t683;
t377 = -t422 + t546;
t757 = -t423 - t547;
t945 = (t344 * t461 + t345 * t460 + t368 * t757 - t370 * t377) * t906 + (-t441 * t452 + t443 * t451 + (-t373 * t610 - t374 * t607) * t563) * t908;
t590 = pkin(2) * t809;
t453 = t590 + (-t606 * t594 - t929) * t607;
t494 = -t694 - t930;
t468 = t503 * t607;
t504 = t609 * t679 + t856;
t719 = -t609 * t468 + t503 * t802 + t504 * t809;
t204 = (t493 * t607 + t453) * t609 + (t494 * t607 - t759) * t606 + t719;
t390 = t606 * t405;
t747 = -t494 - t504;
t733 = t606 * rSges(5,2) * t784 + rSges(5,3) * t791;
t469 = -t808 * t858 + t733;
t592 = pkin(6) * t791;
t790 = t609 * t611;
t751 = -t592 + (-t790 + t928) * t610 + t469;
t205 = t390 + (t610 * t747 + t437) * t606 + (-t610 * t748 - t751) * t609;
t932 = t609 * t758;
t288 = t748 * t808 + t932;
t528 = rSges(4,3) * t606 + t609 * t682;
t311 = (t528 * t607 - t433) * t606;
t732 = t606 * rSges(4,2) * t783 + rSges(4,3) * t791;
t487 = -rSges(4,1) * t606 * t799 + t732;
t312 = (t610 * t946 - t487) * t609 + (-t528 * t610 + t435) * t606;
t361 = t435 * t609 - t808 * t946;
t392 = (t851 * t609 + (t594 + t679) * t606) * t607;
t393 = (-t790 + (-t594 - t858) * t606) * t610 + t733;
t438 = t590 + (-t880 * t609 + t641) * t607;
t439 = t592 + (-pkin(2) - t859) * t808 + t732;
t944 = -(t311 * t373 + t312 * t374 + t359 * t438 - t361 * t439) * m(4) / 0.2e1 + (t204 * t344 + t205 * t345 + t287 * t392 - t288 * t393) * t969;
t943 = t772 * t368 * t969;
t533 = (-Icges(5,5) * t595 - Icges(5,6) * t596) * t606;
t794 = t609 * t533;
t942 = -t794 / 0.2e1 - t606 * t968 / 0.2e1;
t881 = t610 / 0.2e1;
t388 = t403 * t808;
t247 = t388 + (-t436 * t610 - t607 * t758) * t606;
t406 = t422 * t808;
t334 = -t423 * t809 + t406;
t357 = t609 * t422 + t536 * t809;
t358 = -t423 * t609 - t536 * t808;
t416 = -Icges(5,5) * t529 - Icges(5,6) * t530;
t762 = -Icges(5,2) * t530 - t401 - t507;
t764 = -Icges(5,1) * t529 - t397 - t508;
t219 = -t416 * t609 + (-t762 * t595 + t596 * t764) * t606;
t417 = Icges(5,5) * t531 - Icges(5,6) * t532;
t761 = -Icges(5,2) * t532 + t402 + t509;
t763 = Icges(5,1) * t531 - t399 - t845;
t220 = -t417 * t609 + (-t595 * t761 + t596 * t763) * t606;
t705 = t808 / 0.2e1;
t709 = t809 / 0.2e1;
t535 = (-Icges(5,1) * t595 - t843) * t606;
t746 = -t499 + t535;
t307 = -t794 + (t596 * t746 - t968) * t606;
t832 = t307 * t609;
t188 = t416 * t809 - t762 * t529 + t530 * t764;
t189 = t417 * t809 - t529 * t761 + t530 * t763;
t259 = -t529 * t745 + t530 * t746 + t533 * t809;
t84 = -t259 * t609 + (t188 * t607 + t189 * t610) * t606;
t190 = t416 * t808 + t762 * t531 + t532 * t764;
t191 = t417 * t808 + t531 * t761 + t532 * t763;
t260 = t531 * t745 + t532 * t746 + t533 * t808;
t85 = -t260 * t609 + (t190 * t607 + t191 * t610) * t606;
t727 = t85 * t705 + t84 * t709 + (-t832 + (t219 * t607 + t220 * t610) * t606) * t884;
t19 = t727 + m(5) * (t247 * t334 + t287 * t357 - t288 * t358);
t941 = t19 * qJD(4);
t515 = Icges(4,3) * t606 + t609 * t671;
t801 = t608 * t522;
t814 = t605 * t518;
t850 = Icges(3,1) * t606;
t933 = t609 * (t817 / 0.2e1 - t820 / 0.2e1 - t498 / 0.2e1 + t801 / 0.2e1 - t814 / 0.2e1 - t515 / 0.2e1 + t600 + t850 / 0.2e1 - t840 / 0.2e1);
t108 = -t188 * t610 + t189 * t607;
t109 = -t190 * t610 + t191 * t607;
t782 = t108 * t883 + t109 * t886;
t816 = t596 * t502;
t819 = t595 * t500;
t267 = -t645 * t609 + (t497 + t816 - t819) * t606;
t341 = t606 * t658 - t823;
t667 = -t277 * t607 + t278 * t610;
t925 = (-t341 * t609 + t606 * t667) * t887 + ((-t267 + t667) * t609 + (t206 * t607 + t207 * t610 + t341) * t606) * t884;
t924 = t267 * t884 + t341 * t887;
t557 = (-Icges(4,2) * t608 - t847) * t606;
t560 = (-Icges(4,1) * t605 - t846) * t606;
t922 = -t605 * (t522 / 0.2e1 + t557 / 0.2e1) + t608 * (t560 / 0.2e1 - t518 / 0.2e1);
t240 = -t445 * t609 + (-t605 * t754 + t608 * t756) * t606;
t241 = -t446 * t609 + (-t605 * t753 + t608 * t755) * t606;
t554 = (-Icges(4,5) * t605 - Icges(4,6) * t608) * t606;
t741 = t522 + t557;
t743 = -t518 + t560;
t293 = -t550 * t741 + t551 * t743 + t554 * t809;
t294 = t552 * t741 + t553 * t743 + t554 * t808;
t921 = t945 + (t241 + t294) * t885 + (t240 + t293) * t882;
t834 = t220 * t607;
t835 = t219 * t610;
t689 = -t835 / 0.4e1 + t834 / 0.4e1 + t260 * t885 + t259 * t882;
t788 = t610 * t976;
t920 = t788 / 0.4e1 + t976 * t882;
t787 = t610 * t977;
t805 = t607 * t67;
t919 = -t805 / 0.4e1 + t787 / 0.4e1 + t67 * t885 + t977 * t882 + t943 + t982 * t971;
t558 = -Icges(3,2) * t802 - t585;
t570 = Icges(3,2) * t609 + t849;
t559 = t570 * t610;
t677 = -t600 - t850;
t561 = t677 * t607;
t562 = t677 * t610;
t918 = (t610 * (t520 - t561) + (-t521 + t562) * t607) * t609 + ((t524 + t558) * t610 + (-t525 + t559) * t607) * t606;
t482 = t518 * t607;
t484 = t522 * t607;
t647 = -t514 * t607 + t965;
t228 = -t647 * t609 + (t482 * t605 - t484 * t608 + t424) * t606;
t483 = t518 * t610;
t485 = t522 * t610;
t660 = -t429 * t605 + t432 * t608;
t646 = -t514 * t610 - t660;
t229 = -t646 * t609 + (t483 * t605 - t485 * t608 + t426) * t606;
t519 = Icges(4,6) * t606 + t609 * t674;
t523 = Icges(4,5) * t606 + t609 * t676;
t657 = t801 - t814;
t644 = t515 - t657;
t822 = t514 * t609;
t619 = t606 * t644 + t822;
t268 = -t519 * t550 + t523 * t551 + t607 * t619;
t269 = t519 * t552 + t523 * t553 + t610 * t619;
t800 = t608 * t523;
t813 = t605 * t519;
t300 = -t644 * t609 + (t514 + t800 - t813) * t606;
t825 = t426 * t609;
t306 = t606 * t660 - t825;
t352 = t606 * t657 - t822;
t698 = t791 / 0.4e1;
t704 = t808 / 0.4e1;
t917 = t300 * t884 + t352 * t887 - t944 + (t228 + t268) * t708 + (t229 + t269) * t704 + (-t305 + t331) * t701 + (t306 + t333) * t698;
t913 = 4 * qJD(1);
t912 = 2 * qJD(2);
t910 = 2 * qJD(3);
t909 = 4 * qJD(3);
t659 = t433 * t610 - t435 * t607;
t283 = t659 * t609 + (t486 * t610 - t487 * t607) * t606;
t337 = t659 * t606;
t903 = m(4) * (t283 * t337 + t311 * t359 - t312 * t361);
t765 = t403 * t791 - t468 * t808;
t167 = (t453 * t606 - t824) * t610 + (-t606 * t751 - t932) * t607 + t765;
t827 = t405 * t609;
t263 = (-t469 * t606 - t827) * t607 + t765;
t301 = -t403 * t606 + t719;
t302 = -t504 * t808 + t390 + (-t503 * t610 - t469) * t609;
t326 = -t405 * t809 + t388;
t348 = t503 * t808 + t827;
t902 = m(5) * (t167 * t326 + t204 * t346 - t205 * t348 + t247 * t263 + t287 * t301 - t288 * t302);
t900 = t348 * t980;
t899 = m(5) * (t167 * t247 + t204 * t287 - t205 * t288);
t720 = t334 * t230 - t357 * t370 - t358 * t368;
t897 = m(5) * (t247 * t340 + (-t287 * t610 + t288 * t607) * t536 + t720);
t895 = m(5) * (t310 * t326 + t346 * t461 - t348 * t460 + t720);
t894 = t288 * t980;
t891 = m(5) * (t301 * t344 + t302 * t345 + t346 * t392 - t348 * t393);
t771 = t357 * t344 + t358 * t345;
t890 = m(5) * (-t287 * t422 - t288 * t423 + t771);
t889 = -t977 / 0.2e1;
t888 = t240 / 0.2e1;
t860 = rSges(3,1) * t609;
t695 = pkin(1) + t860;
t731 = rSges(3,2) * t809 + t610 * rSges(3,3);
t478 = -t607 * t695 + t602 + t731;
t587 = rSges(3,2) * t808;
t479 = -t587 + t695 * t610 + (rSges(3,3) + pkin(5)) * t607;
t575 = rSges(3,1) * t606 + rSges(3,2) * t609;
t564 = t575 * t607;
t565 = t575 * t610;
t879 = m(3) * (t478 * t564 - t479 * t565);
t873 = m(4) * (t373 * t438 + t374 * t439);
t872 = m(4) * (-t373 * t451 + t374 * t452);
t871 = m(5) * (t346 * t377 + t348 * t757 + t771);
t868 = m(5) * (-t368 * t423 + t370 * t422 + (-t344 * t610 - t345 * t607) * t536);
t867 = m(5) * (t344 * t377 - t345 * t757);
t866 = m(5) * (t344 * t392 + t345 * t393);
t815 = t596 * t606;
t810 = t606 * t520;
t793 = t609 * t554;
t749 = t607 * t517 + t525 * t791;
t738 = -t528 - t580;
t737 = t607 * (pkin(6) * t802 - t590) + t610 * (-pkin(2) * t808 + t592);
t730 = qJD(1) * t606;
t729 = qJD(3) * t606;
t728 = qJD(3) * t609;
t722 = t889 + t977 / 0.2e1;
t717 = -t580 + t747;
t713 = -t815 / 0.2e1;
t712 = t815 / 0.2e1;
t702 = t802 / 0.2e1;
t699 = t791 / 0.2e1;
t690 = t606 * t521 - t516;
t685 = t976 * t709 + t983;
t678 = t499 * t713 + t535 * t712 + t942;
t672 = -Icges(3,5) * t606 - Icges(3,6) * t609;
t664 = -t305 * t607 + t306 * t610;
t655 = -t900 / 0.2e1 + t685;
t624 = t606 * t649 + t829;
t184 = t464 * t529 - t466 * t530 + t607 * t624;
t623 = t606 * t648 + t828;
t185 = t465 * t529 - t467 * t530 + t607 * t623;
t36 = (-t242 + t669) * t609 + (t184 * t607 + t185 * t610 + t321) * t606;
t186 = -t464 * t531 - t466 * t532 + t610 * t624;
t187 = -t465 * t531 - t467 * t532 + t610 * t623;
t37 = (-t243 + t668) * t609 + (t186 * t607 + t187 * t610 + t323) * t606;
t654 = t134 * t702 + t36 * t709 + t37 * t705 + t699 * t976 + t925;
t640 = t902 / 0.2e1 + t654;
t639 = t920 + t975;
t638 = t983 + (t220 + t260) * t705 + (t976 + t219 + t259) * t709;
t637 = t499 * t712 + t535 * t713 - t942;
t356 = -t521 * t808 + t749;
t636 = (t610 * t690 + t356 - t749) * t881 + (-t607 * (-t609 * t524 + t810) - t516 * t610) * t883 + (t607 * t690 + t691 + t947) * t886;
t635 = t356 * t886 - t749 * t607 / 0.2e1 + (-t488 + (t517 + t810) * t610 + t750 + t947) * t883;
t102 = -t184 * t610 + t185 * t607;
t103 = -t186 * t610 + t187 * t607;
t634 = t102 * t709 + t103 * t705 + t36 * t883 + t37 * t886 + (-t206 * t610 + t207 * t607) * t884 + (-t251 * t610 + t252 * t607) * t702 + t972 * t699 + (t277 * t610 + t278 * t607) * t887 - t782;
t626 = t948 * t698 + t949 * t701 + t951 * t704 + t952 * t708 + t924;
t625 = t36 * t711 - t37 * t808 / 0.2e1 + t108 * t709 + t109 * t705 - t134 * t802 / 0.2e1 + t85 * t886 + t84 * t883 - t925 + (t834 - t835 + t788) * t884;
t622 = t606 * t647 + t826;
t621 = t606 * t646 + t825;
t618 = t638 - t832;
t615 = -t800 / 0.2e1 + t813 / 0.2e1 - t816 / 0.2e1 + t819 / 0.2e1 - t514 / 0.2e1 - t497 / 0.2e1 - t573 / 0.2e1 + t570 / 0.2e1;
t614 = t626 + t639 - t689;
t613 = t639 + t689 - t924 + t952 * t710 - t951 * t808 / 0.4e1 - t949 * t802 / 0.4e1 - t948 * t791 / 0.4e1;
t612 = t626 + t689 - t920 + t975;
t603 = t606 ^ 2;
t577 = -rSges(3,2) * t606 + t860;
t556 = t672 * t610;
t555 = t672 * t607;
t444 = t738 * t610;
t442 = t738 * t607;
t376 = -t452 * t609 - t563 * t808;
t375 = t451 * t609 + t563 * t809;
t371 = t717 * t610;
t369 = t717 * t607;
t349 = (t451 * t610 - t452 * t607) * t606;
t335 = t486 * t607 + t487 * t610 + t737;
t330 = (-t536 * t606 + t603 * t865) * t610 + t757 * t609;
t329 = -t546 * t609 - t603 * t726 + t357;
t324 = -t793 + (-t605 * t741 + t608 * t743) * t606;
t303 = t406 + (-t546 * t610 + t607 * t757) * t606;
t266 = t751 * t610 + (t453 - t468) * t607 + t737;
t214 = -t483 * t552 - t485 * t553 + t610 * t621;
t213 = -t482 * t552 - t484 * t553 + t610 * t622;
t212 = t483 * t550 - t485 * t551 + t607 * t621;
t211 = t482 * t550 - t484 * t551 + t607 * t622;
t166 = -t352 * t609 + t606 * t664;
t158 = t678 + t981;
t150 = t868 / 0.2e1;
t137 = t230 * t340 + (t368 * t607 + t370 * t610) * t536;
t118 = -t213 * t610 + t214 * t607;
t117 = -t211 * t610 + t212 * t607;
t112 = t871 / 0.2e1;
t95 = -t294 * t609 + (t223 * t607 + t224 * t610) * t606;
t94 = -t293 * t609 + (t221 * t607 + t222 * t610) * t606;
t92 = t890 / 0.2e1;
t89 = t891 / 0.2e1;
t88 = -t793 / 0.2e1 + t872 + t867 + t922 * t606 + t678;
t87 = t263 * t326 + t301 * t346 - t302 * t348;
t77 = -t606 * t615 + t866 + t873 + t879 + t933;
t72 = (-t300 + t664) * t609 + (t228 * t607 + t229 * t610 + t352) * t606;
t69 = t895 / 0.2e1;
t62 = (-t269 + t665) * t609 + (t213 * t607 + t214 * t610 + t333) * t606;
t61 = (-t268 + t666) * t609 + (t211 * t607 + t212 * t610 + t331) * t606;
t50 = t897 / 0.2e1;
t24 = m(5) * t137 + t782;
t21 = m(5) * (t326 * t334 + t346 * t357 - t348 * t358) + t727;
t20 = t21 * qJD(4);
t18 = t782 + t984;
t16 = m(5) * t87 + t654;
t15 = t112 - t890 / 0.2e1 + t655;
t14 = t92 - t871 / 0.2e1 + t655;
t13 = t92 + t112 + t900 / 0.2e1 + t618;
t12 = t607 * t636 + t610 * t635;
t11 = t722 * t809 + t685 - t894;
t10 = t69 - t897 / 0.2e1 + t640;
t9 = t50 - t895 / 0.2e1 + t640;
t8 = t903 + t899 + (-t72 / 0.2e1 + t787 / 0.2e1 - t805 / 0.2e1) * t609 + (t166 / 0.2e1 + t62 * t881 + t61 * t886) * t606 + t654;
t7 = t614 + t89 - t868 / 0.2e1;
t6 = t150 + t613 - t891 / 0.2e1;
t5 = t150 + t612 + t89;
t4 = t69 + t625 + t50 - t902 / 0.2e1;
t3 = t919 + t917 + t614 + (t293 / 0.4e1 + t240 / 0.4e1) * t610 + (-t294 / 0.4e1 - t241 / 0.4e1) * t607 - t945;
t2 = t919 + t613 + (t300 / 0.2e1 + (-t306 / 0.4e1 - t333 / 0.4e1) * t610 + (t305 / 0.4e1 - t331 / 0.4e1) * t607) * t609 + (-t352 / 0.2e1 + (-t269 / 0.4e1 - t229 / 0.4e1) * t610 + (-t268 / 0.4e1 - t228 / 0.4e1) * t607) * t606 + t921 + t944;
t1 = t917 + t612 + t921 - t943;
t17 = [t77 * qJD(2) + t88 * qJD(3) + t158 * qJD(4), t77 * qJD(1) + t1 * qJD(3) + t5 * qJD(4) + ((t373 * t444 + t374 * t442 - t438 * t443 - t439 * t441) * t908 + (t344 * t371 + t345 * t369 - t368 * t393 - t370 * t392) * t906) * t912 + ((m(3) * (-t478 * t577 - t564 * t575) - t228 / 0.2e1 - t242 / 0.2e1 - t268 / 0.2e1 - t206 / 0.2e1 + t569 * t881 + (-t524 / 0.2e1 - t558 / 0.2e1) * t609 + (t520 / 0.2e1 - t561 / 0.2e1) * t606 - t635) * t610 + (m(3) * (-t479 * t577 + t565 * t575) + t243 / 0.2e1 + t269 / 0.2e1 + t207 / 0.2e1 + t229 / 0.2e1 + t569 * t886 + (t525 / 0.2e1 - t559 / 0.2e1) * t609 + (-t521 / 0.2e1 + t562 / 0.2e1) * t606 - t636) * t607) * qJD(2), t88 * qJD(1) + t1 * qJD(2) + t638 * qJD(3) + t13 * qJD(4) + t894 * t909 / 0.4e1 + ((t287 * t377 + t288 * t757 + t329 * t344 + t330 * t345) * t906 + (-t359 * t451 - t361 * t452 + t373 * t375 + t374 * t376) * t908) * t910 + (-t307 - t324) * t728 + ((t241 / 0.2e1 + t294 / 0.2e1) * t610 + (t293 / 0.2e1 + t888 - t722) * t607) * t729, t158 * qJD(1) + t5 * qJD(2) + t13 * qJD(3) + ((-t346 * t422 - t348 * t423 + t771) * m(5) + t618) * qJD(4); t12 * qJD(2) + t2 * qJD(3) + t6 * qJD(4) + (-t879 / 0.4e1 - t873 / 0.4e1 - t866 / 0.4e1) * t913 + t615 * t730 - t933 * qJD(1), t12 * qJD(1) + t18 * qJD(3) + t24 * qJD(4) + (m(5) * (t230 * t266 - t368 * t369 - t370 * t371) + m(4) * (t313 * t335 - t441 * t442 - t443 * t444) + m(3) * ((t607 * (rSges(3,1) * t802 - t731) + t610 * (rSges(3,1) * t791 + rSges(3,3) * t607 - t587)) * (-t564 * t607 - t565 * t610) + t967 * t577 * t575) + (t103 + t118 + t916 * t556 + (-t607 * t555 + t918) * t610) * t886 + (t102 + t117 + t915 * t555 + (-t610 * t556 + t918) * t607) * t883) * qJD(2), t2 * qJD(1) + t18 * qJD(2) + (t94 * t883 + t95 * t886 + t625) * qJD(3) + t4 * qJD(4) + (-t899 / 0.4e1 - t903 / 0.4e1) * t909 + ((t313 * t349 + t337 * t351 - t375 * t443 - t376 * t441 + (-t359 * t610 + t361 * t607) * t563) * t908 + (t230 * t303 + t247 * t310 + t287 * t461 - t288 * t460 - t329 * t370 - t330 * t368) * t906) * t910 + (t72 / 0.2e1 + (t889 + t888) * t610 + (t67 / 0.2e1 - t241 / 0.2e1) * t607) * t728 + (-t166 / 0.2e1 + (-t62 / 0.2e1 + t133 / 0.2e1) * t610 + (-t61 / 0.2e1 + t132 / 0.2e1) * t607) * t729, t6 * qJD(1) + t24 * qJD(2) + t4 * qJD(3) + ((t326 * t340 + (-t346 * t610 + t348 * t607) * t536 - t87 + t720) * m(5) + t625) * qJD(4); t3 * qJD(2) + t11 * qJD(3) + t14 * qJD(4) + (-t867 / 0.4e1 - t872 / 0.4e1) * t913 - t922 * t730 + (t793 / 0.2e1 + t637 + 0.2e1 * t772 * t345 * t906) * qJD(1), t3 * qJD(1) + t8 * qJD(3) + t9 * qJD(4) + ((t283 * t313 - t311 * t443 - t312 * t441 + t335 * t337 + t359 * t444 - t361 * t442) * t908 + (t167 * t230 - t204 * t370 - t205 * t368 + t247 * t266 + t287 * t371 - t288 * t369) * t906) * t912 + (t634 + t118 * t705 + t971 * t699 + (-t279 * t610 + t280 * t607) * t702 + t117 * t709 + t61 * t883 + (-t228 * t610 + t229 * t607) * t884 + (t305 * t610 + t306 * t607) * t887 + t62 * t886 - t984) * qJD(2), t11 * qJD(1) + t8 * qJD(2) + (t95 * t705 + t94 * t709 + (-t324 * t609 + (t240 * t607 + t241 * t610) * t606) * t884 + t727) * qJD(3) + t941 + ((t247 * t303 + t287 * t329 - t288 * t330) * m(5) / 0.4e1 + (t337 * t349 + t359 * t375 - t361 * t376) * m(4) / 0.4e1) * t909, t14 * qJD(1) + t9 * qJD(2) + t19 * qJD(3) + t941; (t637 - t981) * qJD(1) + t7 * qJD(2) + t15 * qJD(3) + t685 * qJD(4), t7 * qJD(1) + ((t230 * t263 + t266 * t326 - t301 * t370 - t302 * t368 + t346 * t371 - t348 * t369 - t137) * m(5) + t634) * qJD(2) + t10 * qJD(3) + t16 * qJD(4), t15 * qJD(1) + t10 * qJD(2) + ((t303 * t326 + t329 * t346 - t330 * t348) * m(5) + t727) * qJD(3) + t20, qJD(1) * t685 + qJD(2) * t16 + qJD(3) * t21 + t20;];
Cq = t17;