% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP11_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP11_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:29
% EndTime: 2019-12-31 18:54:14
% DurationCPUTime: 37.00s
% Computational Cost: add. (64565->829), mult. (86234->1167), div. (0->0), fcn. (95014->7), ass. (0->486)
t559 = pkin(8) + qJ(3);
t554 = cos(t559);
t563 = sin(qJ(4));
t566 = cos(qJ(1));
t748 = t566 * t563;
t564 = sin(qJ(1));
t565 = cos(qJ(4));
t751 = t564 * t565;
t521 = t554 * t751 - t748;
t499 = Icges(6,5) * t521;
t750 = t565 * t566;
t752 = t564 * t563;
t520 = t554 * t752 + t750;
t553 = sin(t559);
t758 = t553 * t564;
t360 = -Icges(6,6) * t758 - Icges(6,3) * t520 - t499;
t362 = Icges(5,5) * t521 - Icges(5,6) * t520 + Icges(5,3) * t758;
t365 = Icges(6,4) * t521 + Icges(6,2) * t758 + Icges(6,6) * t520;
t502 = Icges(5,4) * t521;
t368 = -Icges(5,2) * t520 + Icges(5,6) * t758 + t502;
t498 = Icges(6,5) * t520;
t371 = Icges(6,1) * t521 + Icges(6,4) * t758 + t498;
t501 = Icges(5,4) * t520;
t375 = -Icges(5,1) * t521 - Icges(5,5) * t758 + t501;
t979 = (t362 + t365) * t758 + (t371 - t375) * t521 + (-t360 - t368) * t520;
t523 = t554 * t750 + t752;
t500 = Icges(6,5) * t523;
t522 = t554 * t748 - t751;
t757 = t553 * t566;
t361 = Icges(6,6) * t757 + Icges(6,3) * t522 + t500;
t364 = Icges(5,5) * t523 - Icges(5,6) * t522 + Icges(5,3) * t757;
t367 = Icges(6,4) * t523 + Icges(6,2) * t757 + Icges(6,6) * t522;
t787 = Icges(5,4) * t523;
t370 = -Icges(5,2) * t522 + Icges(5,6) * t757 + t787;
t782 = Icges(6,5) * t522;
t373 = Icges(6,1) * t523 + Icges(6,4) * t757 + t782;
t503 = Icges(5,4) * t522;
t376 = Icges(5,1) * t523 + Icges(5,5) * t757 - t503;
t978 = (t364 + t367) * t758 + (t373 + t376) * t521 + (t361 - t370) * t520;
t780 = Icges(6,5) * t565;
t538 = t553 * t780;
t759 = t553 * t563;
t776 = Icges(6,6) * t554;
t452 = Icges(6,3) * t759 + t538 - t776;
t634 = Icges(6,4) * t565 + Icges(6,6) * t563;
t456 = -Icges(6,2) * t554 + t553 * t634;
t781 = Icges(6,5) * t563;
t636 = Icges(6,1) * t565 + t781;
t460 = -Icges(6,4) * t554 + t553 * t636;
t270 = t452 * t520 + t456 * t758 + t460 * t521;
t632 = Icges(5,5) * t565 - Icges(5,6) * t563;
t454 = -Icges(5,3) * t554 + t553 * t632;
t785 = Icges(5,4) * t565;
t635 = -Icges(5,2) * t563 + t785;
t458 = -Icges(5,6) * t554 + t553 * t635;
t786 = Icges(5,4) * t563;
t637 = Icges(5,1) * t565 - t786;
t462 = -Icges(5,5) * t554 + t553 * t637;
t271 = t454 * t758 - t458 * t520 + t462 * t521;
t914 = t270 + t271;
t977 = t979 * t564 + t978 * t566;
t846 = t553 / 0.2e1;
t845 = -t554 / 0.2e1;
t791 = -qJ(5) - rSges(6,3);
t835 = -pkin(4) - rSges(6,1);
t974 = -t791 * t563 - t835 * t565;
t940 = t974 * t553;
t843 = t564 / 0.2e1;
t841 = -t566 / 0.2e1;
t320 = t835 * t522 - t791 * t523;
t716 = t835 * t520 - t791 * t521;
t231 = t320 * t566 + t564 * t716;
t415 = -rSges(5,1) * t520 - rSges(5,2) * t521;
t420 = -rSges(5,1) * t522 - rSges(5,2) * t523;
t309 = t564 * t415 + t420 * t566;
t865 = m(6) / 0.2e1;
t868 = -m(5) / 0.2e1;
t735 = -t231 * t865 + t309 * t868;
t706 = -rSges(6,2) * t554 + t940;
t939 = -t791 * t522 - t835 * t523;
t717 = rSges(6,2) * t757 + t939;
t247 = t554 * t717 + t706 * t757;
t891 = t791 * t520 + t835 * t521;
t908 = rSges(6,2) * t758 - t891;
t938 = t554 * t908 + t706 * t758;
t624 = t247 * t564 - t566 * t938;
t642 = t523 * rSges(5,1) - t522 * rSges(5,2);
t382 = rSges(5,3) * t757 + t642;
t796 = rSges(5,1) * t565;
t641 = -rSges(5,2) * t563 + t796;
t593 = t641 * t553;
t465 = -rSges(5,3) * t554 + t593;
t306 = t554 * t382 + t465 * t757;
t889 = -t521 * rSges(5,1) + t520 * rSges(5,2);
t378 = rSges(5,3) * t758 - t889;
t949 = t378 * t554 + t465 * t758;
t893 = t306 * t564 - t566 * t949;
t935 = -m(6) / 0.2e1;
t739 = t624 * t935 + t868 * t893;
t17 = t739 - t735;
t976 = qJD(1) * t17;
t803 = t554 * pkin(3);
t532 = pkin(7) * t553 + t803;
t560 = t564 ^ 2;
t561 = t566 ^ 2;
t941 = t560 + t561;
t699 = t941 * t532;
t193 = t564 * t908 + t566 * t717 + t699;
t253 = t564 * t378 + t382 * t566 + t699;
t531 = pkin(3) * t553 - pkin(7) * t554;
t675 = t531 + t706;
t328 = t675 * t564;
t330 = t675 * t566;
t697 = (t835 * t563 - t791 * t565) * t553;
t391 = t697 * t564;
t393 = t697 * t566;
t506 = (-rSges(5,1) * t563 - rSges(5,2) * t565) * t553;
t705 = t465 + t531;
t387 = t705 * t564;
t389 = t705 * t566;
t887 = t387 * t564 + t389 * t566;
t975 = -m(5) * (t253 * t309 + t887 * t506) - m(6) * (t193 * t231 + t328 * t391 + t330 * t393);
t668 = t914 * t845 + t977 * t846;
t802 = -pkin(6) - qJ(2);
t549 = t564 * t802;
t548 = cos(pkin(8)) * pkin(2) + pkin(1);
t651 = t548 + t803;
t837 = rSges(6,2) + pkin(7);
t880 = t837 * t553 + t651;
t279 = t566 * t880 - t549 + t939;
t653 = t566 * t802;
t952 = -t564 * t880 - t653 + t891;
t959 = -t279 * t522 + t520 * t952;
t973 = m(6) * qJD(1) * t959;
t274 = t522 * t452 + t456 * t757 + t523 * t460;
t216 = -t522 * t360 + t365 * t757 + t523 * t371;
t217 = t522 * t361 + t367 * t757 + t523 * t373;
t628 = t564 * t216 + t217 * t566;
t961 = -t274 * t554 + t553 * t628;
t275 = t454 * t757 - t522 * t458 + t523 * t462;
t218 = t362 * t757 - t522 * t368 - t523 * t375;
t219 = t364 * t757 - t522 * t370 + t523 * t376;
t627 = t564 * t218 + t219 * t566;
t962 = -t275 * t554 + t553 * t627;
t667 = t961 / 0.2e1 + t962 / 0.2e1;
t401 = -Icges(5,5) * t520 - Icges(5,6) * t521;
t403 = -Icges(6,4) * t520 + Icges(6,6) * t521;
t969 = t401 + t403;
t402 = -Icges(5,5) * t522 - Icges(5,6) * t523;
t404 = -Icges(6,4) * t522 + Icges(6,6) * t523;
t968 = t402 + t404;
t721 = Icges(5,2) * t523 - t376 + t503;
t723 = Icges(6,3) * t523 - t373 - t782;
t966 = t721 + t723;
t722 = Icges(5,2) * t521 + t375 + t501;
t724 = Icges(6,3) * t521 - t371 - t498;
t965 = t722 + t724;
t725 = -Icges(5,1) * t522 - t370 - t787;
t727 = -Icges(6,1) * t522 + t361 + t500;
t964 = t725 + t727;
t726 = -Icges(5,1) * t520 - t368 - t502;
t728 = -Icges(6,1) * t520 - t360 + t499;
t963 = t726 + t728;
t894 = -t279 * t564 - t566 * t952;
t768 = t365 * t554;
t937 = t360 * t563 - t371 * t565;
t233 = t937 * t553 + t768;
t770 = t362 * t554;
t936 = t368 * t563 + t375 * t565;
t236 = t936 * t553 + t770;
t623 = t328 * t564 + t330 * t566;
t529 = rSges(4,1) * t553 + rSges(4,2) * t554;
t895 = t941 * t529;
t676 = t623 * t935 + t887 * t868 - m(4) * t895 / 0.2e1;
t755 = t554 * t566;
t543 = rSges(6,2) * t755;
t545 = pkin(7) * t755;
t321 = t543 + t545 + (-pkin(3) - t974) * t757;
t544 = pkin(3) * t758;
t322 = t544 + (-t837 * t554 + t940) * t564;
t836 = rSges(5,3) + pkin(7);
t355 = t544 + (-t836 * t554 + t593) * t564;
t694 = t553 * rSges(5,2) * t748 + rSges(5,3) * t755;
t356 = t545 + (-pkin(3) - t796) * t757 + t694;
t508 = t529 * t564;
t510 = t529 * t566;
t392 = t564 * t508 + t510 * t566;
t867 = m(5) / 0.2e1;
t677 = (-t321 * t566 + t564 * t322) * t865 + (t564 * t355 - t356 * t566) * t867 + m(4) * t392 / 0.2e1;
t61 = t677 - t676;
t956 = t61 * qJD(1);
t954 = -t233 - t236;
t621 = t361 * t563 + t373 * t565;
t767 = t367 * t554;
t234 = t553 * t621 - t767;
t619 = -t370 * t563 + t376 * t565;
t769 = t364 * t554;
t237 = t553 * t619 - t769;
t953 = t234 + t237;
t951 = -t216 * t566 + t217 * t564;
t950 = -t218 * t566 + t219 * t564;
t395 = t564 * t520 + t522 * t566;
t806 = m(6) * t395;
t385 = -t806 / 0.2e1;
t383 = t806 / 0.2e1;
t327 = (t554 * t563 - t395) * t759;
t798 = m(6) * qJD(5);
t948 = t327 * t798;
t947 = t965 * t520 + t963 * t521 + t969 * t758;
t946 = t966 * t520 + t964 * t521 + t968 * t758;
t945 = t965 * t522 + t963 * t523 + t969 * t757;
t944 = t966 * t522 + t964 * t523 + t968 * t757;
t489 = (-Icges(6,4) * t563 + Icges(6,6) * t565) * t553;
t485 = (Icges(6,3) * t565 - t781) * t553;
t708 = -t460 + t485;
t493 = -Icges(6,1) * t759 + t538;
t710 = t452 + t493;
t206 = t489 * t758 + t520 * t708 + t521 * t710;
t486 = (-Icges(5,5) * t563 - Icges(5,6) * t565) * t553;
t490 = (-Icges(5,2) * t565 - t786) * t553;
t707 = -t462 - t490;
t494 = (-Icges(5,1) * t563 - t785) * t553;
t709 = -t458 + t494;
t207 = t486 * t758 + t520 * t707 + t521 * t709;
t943 = -t207 - t206;
t208 = t489 * t757 + t522 * t708 + t523 * t710;
t209 = t486 * t757 + t522 * t707 + t523 * t709;
t942 = -t208 - t209;
t881 = t836 * t553 + t651;
t317 = t566 * t881 - t549 + t642;
t910 = -t564 * t881 - t653 + t889;
t892 = -t317 * t564 - t566 * t910;
t923 = t465 * t564;
t922 = t564 * t706;
t631 = Icges(6,3) * t563 + t780;
t581 = -t553 * t631 + t776;
t425 = t581 * t564;
t433 = t460 * t564;
t605 = t456 * t564 - t937;
t167 = t605 * t554 + (t425 * t563 - t433 * t565 + t365) * t553;
t431 = t458 * t564;
t435 = t462 * t564;
t603 = -t454 * t564 + t936;
t169 = -t603 * t554 + (t431 * t563 - t435 * t565 + t362) * t553;
t921 = t167 + t169;
t426 = t581 * t566;
t434 = t460 * t566;
t604 = t456 * t566 + t621;
t168 = t604 * t554 + (t426 * t563 - t434 * t565 + t367) * t553;
t432 = t458 * t566;
t436 = t462 * t566;
t602 = -t454 * t566 - t619;
t170 = -t602 * t554 + (t432 * t563 - t436 * t565 + t364) * t553;
t920 = t168 + t170;
t179 = -t403 * t554 + (t563 * t724 + t565 * t728) * t553;
t181 = -t401 * t554 + (t563 * t722 + t565 * t726) * t553;
t919 = t179 + t181;
t180 = -t404 * t554 + (t563 * t723 + t565 * t727) * t553;
t182 = -t402 * t554 + (t563 * t721 + t565 * t725) * t553;
t918 = t180 + t182;
t453 = Icges(6,6) * t553 + t554 * t631;
t461 = Icges(6,4) * t553 + t554 * t636;
t457 = Icges(6,2) * t553 + t554 * t634;
t617 = t563 * t452 + t565 * t460;
t601 = -t457 + t617;
t764 = t456 * t554;
t573 = -t553 * t601 + t764;
t189 = t453 * t520 + t461 * t521 + t564 * t573;
t459 = Icges(5,6) * t553 + t554 * t635;
t463 = Icges(5,5) * t553 + t554 * t637;
t455 = Icges(5,3) * t553 + t554 * t632;
t616 = -t563 * t458 + t565 * t462;
t600 = t455 - t616;
t765 = t454 * t554;
t574 = t553 * t600 + t765;
t190 = -t459 * t520 + t463 * t521 + t564 * t574;
t917 = t189 + t190;
t191 = t522 * t453 + t523 * t461 + t566 * t573;
t192 = -t522 * t459 + t523 * t463 + t566 * t574;
t916 = t191 + t192;
t199 = t601 * t554 + (t563 * t453 + t565 * t461 + t456) * t553;
t200 = -t600 * t554 + (-t563 * t459 + t565 * t463 + t454) * t553;
t915 = t199 + t200;
t913 = t274 + t275;
t292 = t553 * t617 - t764;
t293 = t553 * t616 - t765;
t912 = t292 + t293;
t756 = t554 * t564;
t472 = Icges(4,4) * t756 - Icges(4,2) * t758 - Icges(4,6) * t566;
t547 = Icges(4,4) * t554;
t779 = Icges(4,2) * t553;
t473 = Icges(4,6) * t564 + (t547 - t779) * t566;
t788 = Icges(4,4) * t553;
t528 = Icges(4,1) * t554 - t788;
t475 = Icges(4,5) * t564 + t528 * t566;
t443 = t475 * t756;
t524 = Icges(4,5) * t554 - Icges(4,6) * t553;
t471 = Icges(4,3) * t564 + t524 * t566;
t646 = t471 * t566 - t443;
t470 = Icges(4,5) * t756 - Icges(4,6) * t758 - Icges(4,3) * t566;
t539 = Icges(4,4) * t758;
t474 = Icges(4,1) * t756 - Icges(4,5) * t566 - t539;
t712 = -t564 * t470 - t474 * t755;
t911 = -t472 * t757 - t473 * t758 - t646 - t712;
t909 = t954 * t564 + t953 * t566;
t906 = -t553 / 0.2e1;
t905 = t554 / 0.2e1;
t904 = -t564 / 0.2e1;
t839 = t566 / 0.2e1;
t902 = t943 * t554 + (t947 * t564 + t946 * t566) * t553;
t901 = t942 * t554 + (t945 * t564 + t944 * t566) * t553;
t886 = t167 / 0.2e1 + t169 / 0.2e1;
t885 = -t168 / 0.2e1 - t170 / 0.2e1;
t884 = (-t387 * t420 + t389 * t415 + t892 * t506) * t868 + (-t279 * t391 - t320 * t328 + t330 * t716 - t393 * t952) * t935;
t704 = rSges(6,2) * t553 + t554 * t974;
t171 = (t564 * t704 - t908) * t553;
t713 = -t940 * t566 + t543;
t172 = (-t566 * t706 - t713) * t554 + (-t566 * t704 + t717) * t553;
t467 = rSges(5,3) * t553 + t554 * t641;
t243 = (t467 * t564 - t378) * t553;
t442 = -rSges(5,1) * t553 * t750 + t694;
t244 = (-t465 * t566 - t442) * t554 + (-t467 * t566 + t382) * t553;
t883 = (t171 * t952 + t172 * t279 - t247 * t321 + t322 * t938) * t935 + (t243 * t910 + t244 * t317 - t306 * t356 + t355 * t949) * t868;
t657 = -t462 / 0.2e1 - t460 / 0.2e1;
t659 = t458 / 0.2e1 - t452 / 0.2e1;
t879 = -t563 * (t490 / 0.2e1 - t485 / 0.2e1 - t657) - t565 * (-t494 / 0.2e1 - t493 / 0.2e1 + t659);
t789 = Icges(4,1) * t553;
t877 = t563 * t659 + t565 * t657 + t455 / 0.2e1 + t457 / 0.2e1 - t547 + t779 / 0.2e1 - t789 / 0.2e1;
t525 = Icges(4,2) * t554 + t788;
t876 = t563 * (t459 / 0.2e1 - t453 / 0.2e1) + t565 * (-t463 / 0.2e1 - t461 / 0.2e1) - t454 / 0.2e1 - t456 / 0.2e1 + t525 / 0.2e1 - t528 / 0.2e1;
t491 = -Icges(4,2) * t756 - t539;
t492 = t525 * t566;
t638 = -t547 - t789;
t495 = t638 * t564;
t496 = t638 * t566;
t875 = (t566 * (t472 - t495) + (-t473 + t496) * t564) * t554 + (t566 * (t474 + t491) + (-t475 + t492) * t564) * t553;
t873 = 0.4e1 * qJD(1);
t872 = 2 * qJD(3);
t870 = 2 * qJD(4);
t869 = 4 * qJD(4);
t618 = t378 * t566 - t382 * t564;
t201 = t618 * t554 + (-t442 * t564 - t566 * t923) * t553;
t285 = t618 * t553;
t863 = m(5) * (t201 * t285 + t243 * t949 - t244 * t306);
t588 = -t564 * t717 + t566 * t908;
t137 = t588 * t554 + (-t564 * t713 - t566 * t922) * t553;
t198 = t588 * t553;
t859 = m(6) * (t522 * t171 + t520 * t172 + (t198 * t554 + (t137 + t624) * t553) * t563);
t858 = m(6) * (t137 * t198 + t171 * t938 - t172 * t247);
t354 = (t520 * t566 - t522 * t564) * t553;
t552 = t553 ^ 2;
t423 = t520 * t554 + t552 * t752;
t424 = -t554 * t522 - t552 * t748;
t856 = m(6) * (t354 * t193 + t198 * t395 - t424 * t328 - t423 * t330 + t624 * t759);
t852 = m(6) * (-t328 * t521 - t330 * t523 - t391 * t520 - t393 * t522 + (t193 * t565 + t231 * t563) * t553);
t851 = m(6) * (-t522 * t247 + t279 * t424 + t423 * t952 - t520 * t938);
t840 = -t566 / 0.4e1;
t834 = m(3) * t941 * (rSges(3,3) + qJ(2));
t797 = rSges(4,1) * t554;
t650 = t548 + t797;
t693 = rSges(4,2) * t758 + t566 * rSges(4,3);
t411 = -t564 * t650 - t653 + t693;
t649 = -rSges(4,2) * t757 + t564 * rSges(4,3);
t412 = t566 * t650 - t549 + t649;
t833 = m(4) * (t411 * t508 - t412 * t510);
t832 = m(4) * (t411 * t566 + t412 * t564);
t825 = m(5) * (t317 * t356 + t355 * t910);
t824 = m(5) * (t317 * t420 - t415 * t910);
t822 = m(5) * t892;
t817 = m(6) * (t279 * t521 + t320 * t520 - t522 * t716 + t523 * t952);
t589 = t894 * t759;
t816 = m(6) * (t520 * t321 + t522 * t322 + t589);
t815 = m(6) * (-t522 * t328 + t520 * t330 + t589);
t814 = m(6) * (t279 * t320 - t716 * t952);
t813 = m(6) * (t279 * t321 + t322 * t952);
t811 = m(6) * t894;
t800 = m(6) * qJD(3);
t799 = m(6) * qJD(4);
t762 = t472 * t553;
t569 = (t243 * t564 - t244 * t566) * t868 + (t171 * t564 - t172 * t566) * t935;
t591 = (t391 * t566 - t393 * t564) * t865;
t56 = t591 + t569;
t754 = t56 * qJD(2);
t711 = t564 * t471 + t475 * t755;
t703 = -t467 - t532;
t700 = t564 * (pkin(7) * t756 - t544) + t566 * (-pkin(3) * t757 + t545);
t691 = qJD(1) * t553;
t690 = qJD(1) * t554;
t689 = qJD(3) * t564;
t688 = qJD(3) * t566;
t687 = qJD(4) * t553;
t686 = qJD(4) * t554;
t590 = (t423 * t564 - t424 * t566) * t865;
t595 = m(6) * (-t521 * t566 + t523 * t564);
t269 = t590 - t595 / 0.2e1;
t685 = t269 * qJD(2);
t576 = -t553 * t605 + t768;
t147 = t425 * t520 - t433 * t521 + t564 * t576;
t575 = -t553 * t604 + t767;
t148 = t426 * t520 - t434 * t521 + t564 * t575;
t578 = t553 * t603 + t770;
t149 = t431 * t520 - t435 * t521 + t564 * t578;
t577 = t553 * t602 + t769;
t150 = t432 * t520 - t436 * t521 + t564 * t577;
t682 = (-t917 + t977) * t905 + ((t148 + t150) * t566 + (t147 + t149) * t564 + t914) * t846;
t151 = t522 * t425 - t523 * t433 + t566 * t576;
t152 = t522 * t426 - t523 * t434 + t566 * t575;
t153 = t522 * t431 - t523 * t435 + t566 * t578;
t154 = t522 * t432 - t523 * t436 + t566 * t577;
t681 = (t628 + t627 - t916) * t905 + ((t152 + t154) * t566 + (t151 + t153) * t564 + t913) * t846;
t680 = (t921 * t564 + t920 * t566 + t912) * t906 + (t909 - t915) * t845;
t679 = t947 * t839 + t946 * t904;
t678 = t945 * t841 + t944 * t843;
t674 = -t532 - t704;
t672 = t758 / 0.4e1;
t666 = t912 * t905 + t909 * t906;
t662 = t179 / 0.2e1 + t181 / 0.2e1;
t661 = -t180 / 0.2e1 - t182 / 0.2e1;
t655 = -t489 / 0.2e1 - t486 / 0.2e1;
t645 = t473 * t553 - t470;
t633 = -Icges(4,5) * t553 - Icges(4,6) * t554;
t613 = -t679 - t682;
t612 = -t678 + t681;
t300 = -t473 * t757 + t711;
t587 = t300 * t843 + t711 * t904 + (-t443 + (t471 + t762) * t566 + t712 + t911) * t841;
t586 = (t566 * t645 + t300 - t711) * t839 + (-t564 * (-t474 * t554 + t762) - t470 * t566) * t841 + (t564 * t645 + t646 + t911) * t843;
t579 = -t884 + (t918 - t942) * t564 / 0.4e1 + (t919 - t943) * t840;
t568 = (-t758 / 0.4e1 + t672) * (t951 + t950) + (t840 + t566 / 0.4e1) * (t961 + t962);
t567 = -t883 + t912 * t846 + t915 * t845 + (t917 + t921) * t672 + (t916 + t920) * t757 / 0.4e1 + (t914 + t954) * t756 / 0.4e1 + (t913 + t953) * t755 / 0.4e1;
t530 = -rSges(4,2) * t553 + t797;
t488 = t633 * t566;
t487 = t633 * t564;
t390 = t703 * t566;
t388 = t703 * t564;
t331 = t674 * t566;
t329 = t674 * t564;
t324 = -t554 * t420 - t506 * t757;
t323 = t415 * t554 + t506 * t758;
t314 = t552 * t563 * t565 + t520 * t521 + t522 * t523;
t296 = (t415 * t566 - t420 * t564) * t553;
t281 = t442 * t566 - t564 * t923 + t700;
t268 = t590 + t595 / 0.2e1;
t261 = 0.2e1 * t385;
t260 = t385 + t383;
t259 = 0.2e1 * t383;
t258 = -t320 * t554 - t393 * t553;
t257 = t554 * t716 + t697 * t758;
t252 = -t554 * t486 + (t563 * t707 + t565 * t709) * t553;
t251 = -t554 * t489 + (t563 * t708 + t565 * t710) * t553;
t229 = -t564 * t922 + t566 * t713 + t700;
t220 = (-t320 * t564 + t566 * t716) * t553;
t113 = t193 * t395 + t623 * t759;
t111 = t815 / 0.2e1;
t109 = t816 / 0.2e1;
t105 = t817 / 0.2e1;
t92 = -t811 - t822 + t832 + t834;
t81 = t198 * t354 - t247 * t424 + t423 * t938;
t80 = -t153 * t566 + t154 * t564;
t79 = -t151 * t566 + t152 * t564;
t78 = -t149 * t566 + t150 * t564;
t77 = -t147 * t566 + t148 * t564;
t73 = t851 / 0.2e1;
t71 = t852 / 0.2e1;
t60 = t676 + t677;
t57 = t553 * t879 + t655 * t554 + t814 + t824;
t55 = t591 - t569;
t50 = t856 / 0.2e1;
t39 = -t553 * t876 - t554 * t877 + t813 + t825 + t833;
t23 = t109 - t815 / 0.2e1;
t22 = t111 + t109;
t21 = t111 - t816 / 0.2e1;
t19 = t859 / 0.2e1;
t15 = t735 + t739;
t13 = t105 - t851 / 0.2e1;
t12 = t73 + t105;
t11 = t73 - t817 / 0.2e1;
t10 = t71 + t19 - t856 / 0.2e1;
t9 = t50 + t71 - t859 / 0.2e1;
t8 = t50 + t19 - t852 / 0.2e1;
t7 = t564 * t678 + t566 * t679 - t975;
t6 = t564 * t586 + t566 * t587;
t4 = t863 + t858 + (t564 * t668 + t566 * t667 + t680) * t554 + (t564 * t682 + t566 * t681 - t666) * t553;
t3 = t579 + t567;
t2 = t579 + (-t293 / 0.2e1 - t292 / 0.2e1 + (-t192 / 0.4e1 - t191 / 0.4e1 - t170 / 0.4e1 - t168 / 0.4e1) * t566 + (-t190 / 0.4e1 - t189 / 0.4e1 - t169 / 0.4e1 - t167 / 0.4e1) * t564) * t553 + (t200 / 0.2e1 + t199 / 0.2e1 + (-t275 / 0.4e1 - t274 / 0.4e1 - t237 / 0.4e1 - t234 / 0.4e1) * t566 + (t236 / 0.4e1 + t233 / 0.4e1 - t271 / 0.4e1 - t270 / 0.4e1) * t564) * t554 + t568 + t883;
t1 = (t181 / 0.4e1 + t179 / 0.4e1 + t207 / 0.4e1 + t206 / 0.4e1) * t566 + (-t209 / 0.4e1 - t208 / 0.4e1 - t182 / 0.4e1 - t180 / 0.4e1) * t564 + t568 + t567 + t884;
t5 = [t92 * qJD(2) + t39 * qJD(3) + t57 * qJD(4) - t798 * t959, qJD(1) * t92 + qJD(3) * t60 + qJD(4) * t15 + qJD(5) * t260, t39 * qJD(1) + t60 * qJD(2) + t3 * qJD(4) + t22 * qJD(5) + ((t317 * t388 - t355 * t389 - t356 * t387 + t390 * t910) * t867 + (t279 * t329 - t321 * t328 - t322 * t330 + t331 * t952) * t865) * t872 + (m(4) * (-t411 * t530 - t508 * t529) - t189 / 0.2e1 - t190 / 0.2e1 + t524 * t839 + (-t474 / 0.2e1 - t491 / 0.2e1) * t554 + (t472 / 0.2e1 - t495 / 0.2e1) * t553 - t587 - t886) * t688 + (m(4) * (-t412 * t530 + t510 * t529) + t191 / 0.2e1 + t192 / 0.2e1 + t524 * t843 + (t475 / 0.2e1 - t492 / 0.2e1) * t554 + (-t473 / 0.2e1 + t496 / 0.2e1) * t553 - t586 - t885) * t689, t57 * qJD(1) + t15 * qJD(2) + t3 * qJD(3) + t12 * qJD(5) + ((-t247 * t320 + t257 * t952 + t258 * t279 - t716 * t938) * t865 + (-t306 * t420 + t317 * t324 + t323 * t910 - t415 * t949) * t867) * t870 + (-t251 - t252) * t686 + ((t209 / 0.2e1 + t208 / 0.2e1 - t661) * t566 + (t206 / 0.2e1 + t207 / 0.2e1 + t662) * t564) * t687, t260 * qJD(2) + t22 * qJD(3) + t12 * qJD(4) - t973; t61 * qJD(3) - t17 * qJD(4) + t261 * qJD(5) + (t811 / 0.4e1 + t822 / 0.4e1 - t832 / 0.4e1 - t834 / 0.4e1) * t873, 0, t956 + ((-t388 * t566 + t390 * t564) * t867 + (-t329 * t566 + t331 * t564) * t865) * t872 + t55 * qJD(4), -t976 + t55 * qJD(3) + ((t323 * t564 - t324 * t566) * t867 + (t257 * t564 - t258 * t566) * t865) * t870 + t268 * qJD(5), qJD(1) * t261 + qJD(4) * t268; -t61 * qJD(2) + t6 * qJD(3) + t2 * qJD(4) + t21 * qJD(5) + (-t833 / 0.4e1 - t825 / 0.4e1 - t813 / 0.4e1) * t873 + t877 * t690 + t876 * t691, qJD(4) * t56 - t956, t6 * qJD(1) + t7 * qJD(4) + t113 * t798 + (m(5) * (t253 * t281 - t387 * t388 - t389 * t390) + m(6) * (t193 * t229 - t328 * t329 - t330 * t331) + m(4) * (t530 * t895 - (t564 * (rSges(4,1) * t756 - t693) + t566 * (rSges(4,1) * t755 + t649)) * t392) + (t79 + t80 + t560 * t488 + (-t564 * t487 + t875) * t566) * t843 + (t77 + t78 + t561 * t487 + (-t566 * t488 + t875) * t564) * t841) * qJD(3), t2 * qJD(1) + t754 + t7 * qJD(3) + (t902 * t841 + t901 * t843) * qJD(4) + t9 * qJD(5) + (-t863 / 0.4e1 - t858 / 0.4e1) * t869 + ((t296 * t253 + t285 * t309 - t323 * t389 - t324 * t387 + t893 * t506) * t867 + (t193 * t220 + t198 * t231 + t247 * t391 - t257 * t330 - t258 * t328 - t393 * t938) * t865) * t870 + ((t662 - t667) * t566 + (t661 - t668) * t564 - t680) * t686 + (t564 * t613 - t566 * t612 + t666) * t687, t21 * qJD(1) + t9 * qJD(4) + t113 * t800 - t948; t17 * qJD(2) + t1 * qJD(3) + t11 * qJD(5) - t655 * t690 + (-t814 / 0.4e1 - t824 / 0.4e1) * t873 - t879 * t691, -qJD(3) * t56 + qJD(5) * t269 + t976, t1 * qJD(1) - t754 + t4 * qJD(4) + t8 * qJD(5) + ((t201 * t253 - t243 * t389 - t244 * t387 + t281 * t285 - t306 * t388 + t390 * t949) * t867 + (t137 * t193 - t171 * t330 - t172 * t328 + t198 * t229 - t247 * t329 + t331 * t938) * t865) * t872 + t613 * t688 + t612 * t689 + (((t951 / 0.2e1 + t950 / 0.2e1 + t886) * t566 + (t979 * t841 + t978 * t843 + t885) * t564) * t554 + ((t80 / 0.2e1 + t79 / 0.2e1 + t233 / 0.2e1 + t236 / 0.2e1) * t566 + (t237 / 0.2e1 + t234 / 0.2e1 + t77 / 0.2e1 + t78 / 0.2e1) * t564) * t553 + t975) * qJD(3), t4 * qJD(3) + (t252 / 0.2e1 + t251 / 0.2e1) * qJD(4) * t554 ^ 2 + (m(6) * (t198 * t220 - t247 * t258 + t257 * t938) / 0.4e1 + (t285 * t296 - t306 * t324 + t323 * t949) * m(5) / 0.4e1) * t869 + t81 * t798 + ((t919 * t564 + t918 * t566) * t845 + t902 * t843 + t901 * t839) * t687, t11 * qJD(1) + t685 + t8 * qJD(3) + t81 * t799 + (t354 * t759 + t423 * t522 + t424 * t520 - t314) * t798; t259 * qJD(2) + t23 * qJD(3) + t13 * qJD(4) + t973, qJD(1) * t259 - qJD(4) * t269, t23 * qJD(1) + (t520 * t329 + t522 * t331 - t113 + (t193 * t554 + (t229 + t623) * t553) * t563) * t800 + t10 * qJD(4) + t948, t13 * qJD(1) - t685 + t10 * qJD(3) + (t938 * t523 - t247 * t521 + t257 * t522 + t258 * t520 + (t198 * t565 + t220 * t563) * t553 - t81) * t799 + t314 * t798, 0.4e1 * (t327 * qJD(3) / 0.4e1 + t314 * qJD(4) / 0.4e1) * m(6);];
Cq = t5;
