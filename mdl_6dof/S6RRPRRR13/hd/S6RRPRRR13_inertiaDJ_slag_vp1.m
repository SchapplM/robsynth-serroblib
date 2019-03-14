% Calculate time derivative of joint inertia matrix for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR13_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:45:11
% EndTime: 2019-03-09 14:46:11
% DurationCPUTime: 34.85s
% Computational Cost: add. (115434->1505), mult. (272147->1977), div. (0->0), fcn. (305733->12), ass. (0->605)
t593 = sin(qJ(1));
t589 = cos(pkin(6));
t645 = qJD(2) * t589 + qJD(1);
t592 = sin(qJ(2));
t736 = t592 * t593;
t688 = t589 * t736;
t694 = qJD(2) * t592;
t595 = cos(qJ(2));
t596 = cos(qJ(1));
t733 = t595 * t596;
t489 = -qJD(1) * t688 - t593 * t694 + t645 * t733;
t590 = sin(qJ(5));
t687 = t589 * t733;
t556 = -t687 + t736;
t591 = sin(qJ(4));
t588 = sin(pkin(6));
t763 = cos(qJ(4));
t671 = t588 * t763;
t641 = t596 * t671;
t514 = t556 * t591 - t641;
t734 = t593 * t595;
t735 = t592 * t596;
t557 = t589 * t735 + t734;
t594 = cos(qJ(5));
t628 = t514 * t590 - t557 * t594;
t780 = qJD(5) * t628 - t489 * t590;
t800 = pkin(5) * t780;
t558 = t589 * t734 + t735;
t559 = -t688 + t733;
t740 = t588 * t593;
t459 = Icges(4,1) * t740 - Icges(4,4) * t559 + Icges(4,5) * t558;
t462 = Icges(3,5) * t559 - Icges(3,6) * t558 + Icges(3,3) * t740;
t787 = t459 + t462;
t738 = t588 * t596;
t460 = -Icges(4,1) * t738 - Icges(4,4) * t557 + Icges(4,5) * t556;
t461 = Icges(3,5) * t557 - Icges(3,6) * t556 - Icges(3,3) * t738;
t786 = t460 + t461;
t456 = -Icges(4,5) * t738 - Icges(4,6) * t557 + Icges(4,3) * t556;
t463 = Icges(3,4) * t557 - Icges(3,2) * t556 - Icges(3,6) * t738;
t785 = -t463 + t456;
t455 = Icges(4,5) * t740 - Icges(4,6) * t559 + Icges(4,3) * t558;
t464 = Icges(3,4) * t559 - Icges(3,2) * t558 + Icges(3,6) * t740;
t784 = t464 - t455;
t458 = -Icges(4,4) * t738 - Icges(4,2) * t557 + Icges(4,6) * t556;
t465 = Icges(3,1) * t557 - Icges(3,4) * t556 - Icges(3,5) * t738;
t783 = t465 - t458;
t457 = Icges(4,4) * t740 - Icges(4,2) * t559 + Icges(4,6) * t558;
t466 = Icges(3,1) * t559 - Icges(3,4) * t558 + Icges(3,5) * t740;
t782 = -t466 + t457;
t753 = Icges(3,4) * t592;
t527 = Icges(3,6) * t589 + (Icges(3,2) * t595 + t753) * t588;
t752 = Icges(3,4) * t595;
t528 = Icges(3,5) * t589 + (Icges(3,1) * t592 + t752) * t588;
t751 = Icges(4,6) * t592;
t529 = Icges(4,5) * t589 + (-Icges(4,3) * t595 - t751) * t588;
t695 = qJD(2) * t588;
t539 = (Icges(3,5) * t595 - Icges(3,6) * t592) * t695;
t540 = (-Icges(4,4) * t595 + Icges(4,5) * t592) * t695;
t541 = (-Icges(3,2) * t592 + t752) * t695;
t542 = (Icges(3,1) * t595 - t753) * t695;
t693 = qJD(2) * t595;
t665 = t588 * t693;
t666 = t588 * t694;
t739 = t588 * t595;
t741 = t588 * t592;
t799 = t528 * t665 + t541 * t739 + t542 * t741 + (-t527 + t529) * t666 + (t539 + t540) * t589;
t512 = t558 * t591 + t593 * t671;
t443 = -t512 * t590 + t559 * t594;
t487 = qJD(1) * t557 + qJD(2) * t558;
t798 = qJD(5) * t443 - t487 * t590;
t696 = qJD(1) * t596;
t668 = t588 * t696;
t797 = -t558 * t784 - t559 * t782 + t740 * t787;
t796 = -t556 * t785 - t557 * t783 + t738 * t786;
t795 = t556 * t784 + t557 * t782 + t738 * t787;
t794 = t558 * t785 + t559 * t783 + t740 * t786;
t750 = Icges(4,6) * t595;
t530 = Icges(4,4) * t589 + (-Icges(4,2) * t592 - t750) * t588;
t537 = (Icges(4,3) * t592 - t750) * t695;
t538 = (-Icges(4,2) * t595 + t751) * t695;
t737 = t592 * t538;
t793 = ((-t737 + (-qJD(2) * t530 - t537) * t595) * t588 + t799) * t589;
t486 = -qJD(1) * t687 - t596 * t693 + t645 * t736;
t372 = Icges(4,1) * t668 + Icges(4,4) * t487 - Icges(4,5) * t486;
t373 = -Icges(3,5) * t487 + Icges(3,6) * t486 + Icges(3,3) * t668;
t792 = t373 + t372;
t368 = Icges(4,5) * t668 + Icges(4,6) * t487 - Icges(4,3) * t486;
t375 = -Icges(3,4) * t487 + Icges(3,2) * t486 + Icges(3,6) * t668;
t791 = t375 - t368;
t488 = qJD(1) * t558 + qJD(2) * t557;
t697 = qJD(1) * t593;
t669 = t588 * t697;
t367 = Icges(4,5) * t669 - Icges(4,6) * t489 + Icges(4,3) * t488;
t376 = Icges(3,4) * t489 - Icges(3,2) * t488 + Icges(3,6) * t669;
t790 = t376 - t367;
t370 = Icges(4,4) * t668 + Icges(4,2) * t487 - Icges(4,6) * t486;
t377 = -Icges(3,1) * t487 + Icges(3,4) * t486 + Icges(3,5) * t668;
t789 = t377 - t370;
t369 = Icges(4,4) * t669 - Icges(4,2) * t489 + Icges(4,6) * t488;
t378 = Icges(3,1) * t489 - Icges(3,4) * t488 + Icges(3,5) * t669;
t788 = t378 - t369;
t555 = t589 * t763 - t591 * t739;
t371 = Icges(4,1) * t669 - Icges(4,4) * t489 + Icges(4,5) * t488;
t374 = Icges(3,5) * t489 - Icges(3,6) * t488 + Icges(3,3) * t669;
t781 = (-t374 - t371) * t596;
t585 = t588 ^ 2;
t779 = m(7) / 0.2e1;
t390 = -t488 * t763 - qJD(4) * t641 + (qJD(4) * t556 + t669) * t591;
t778 = t390 / 0.2e1;
t392 = qJD(4) * t512 + t486 * t763 + t591 * t668;
t777 = t392 / 0.2e1;
t776 = -t487 / 0.2e1;
t775 = t489 / 0.2e1;
t507 = qJD(4) * t555 - t666 * t763;
t774 = t507 / 0.2e1;
t511 = -t558 * t763 + t591 * t740;
t773 = t511 / 0.2e1;
t513 = t556 * t763 + t591 * t738;
t772 = -t513 / 0.2e1;
t619 = -t589 * t591 - t595 * t671;
t771 = -t619 / 0.2e1;
t770 = t557 / 0.2e1;
t769 = t559 / 0.2e1;
t768 = t588 / 0.2e1;
t767 = t589 / 0.2e1;
t766 = t593 / 0.2e1;
t765 = rSges(4,2) - pkin(2);
t764 = -rSges(5,3) - pkin(2);
t762 = t489 * pkin(2);
t761 = t557 * pkin(2);
t584 = t596 * pkin(1);
t581 = pkin(5) * t594 + pkin(4);
t760 = -pkin(4) + t581;
t597 = -pkin(11) - pkin(10);
t759 = -pkin(10) - t597;
t758 = rSges(4,3) * t488;
t587 = qJ(5) + qJ(6);
t582 = sin(t587);
t583 = cos(t587);
t586 = qJD(5) + qJD(6);
t648 = -t512 * t586 - t487;
t640 = qJD(1) * t671;
t393 = -qJD(4) * t511 - t486 * t591 + t596 * t640;
t649 = t559 * t586 + t393;
t268 = -t582 * t649 + t583 * t648;
t269 = t582 * t648 + t583 * t649;
t171 = Icges(7,5) * t269 + Icges(7,6) * t268 + Icges(7,3) * t392;
t173 = Icges(7,4) * t269 + Icges(7,2) * t268 + Icges(7,6) * t392;
t175 = Icges(7,1) * t269 + Icges(7,4) * t268 + Icges(7,5) * t392;
t431 = -t512 * t582 + t559 * t583;
t432 = t512 * t583 + t559 * t582;
t315 = Icges(7,5) * t432 + Icges(7,6) * t431 + Icges(7,3) * t511;
t317 = Icges(7,4) * t432 + Icges(7,2) * t431 + Icges(7,6) * t511;
t319 = Icges(7,1) * t432 + Icges(7,4) * t431 + Icges(7,5) * t511;
t622 = -t555 * t586 + t665;
t506 = qJD(4) * t619 + t591 * t666;
t636 = t586 * t741 + t506;
t365 = -t582 * t636 + t583 * t622;
t366 = t582 * t622 + t583 * t636;
t494 = -t555 * t582 + t583 * t741;
t495 = t555 * t583 + t582 * t741;
t56 = -t171 * t619 + t173 * t494 + t175 * t495 + t315 * t507 + t317 * t365 + t319 * t366;
t756 = t56 * t511;
t647 = -t514 * t586 + t489;
t391 = qJD(4) * t513 + t488 * t591 + t593 * t640;
t650 = t557 * t586 + t391;
t266 = -t582 * t650 + t583 * t647;
t267 = t582 * t647 + t583 * t650;
t170 = Icges(7,5) * t267 + Icges(7,6) * t266 + Icges(7,3) * t390;
t172 = Icges(7,4) * t267 + Icges(7,2) * t266 + Icges(7,6) * t390;
t174 = Icges(7,1) * t267 + Icges(7,4) * t266 + Icges(7,5) * t390;
t433 = -t514 * t582 + t557 * t583;
t434 = t514 * t583 + t557 * t582;
t316 = Icges(7,5) * t434 + Icges(7,6) * t433 - Icges(7,3) * t513;
t318 = Icges(7,4) * t434 + Icges(7,2) * t433 - Icges(7,6) * t513;
t320 = Icges(7,1) * t434 + Icges(7,4) * t433 - Icges(7,5) * t513;
t57 = -t170 * t619 + t172 * t494 + t174 * t495 + t316 * t507 + t318 * t365 + t320 * t366;
t755 = t57 * t513;
t754 = rSges(7,3) - t597;
t162 = -t315 * t619 + t317 * t494 + t319 * t495;
t749 = t162 * t392;
t163 = -t316 * t619 + t318 * t494 + t320 * t495;
t748 = t163 * t390;
t744 = t557 * t590;
t743 = t559 * t590;
t246 = Icges(7,5) * t366 + Icges(7,6) * t365 + Icges(7,3) * t507;
t247 = Icges(7,4) * t366 + Icges(7,2) * t365 + Icges(7,6) * t507;
t248 = Icges(7,1) * t366 + Icges(7,4) * t365 + Icges(7,5) * t507;
t359 = Icges(7,5) * t495 + Icges(7,6) * t494 - Icges(7,3) * t619;
t360 = Icges(7,4) * t495 + Icges(7,2) * t494 - Icges(7,6) * t619;
t361 = Icges(7,1) * t495 + Icges(7,4) * t494 - Icges(7,5) * t619;
t112 = -t246 * t619 + t494 * t247 + t495 * t248 + t507 * t359 + t365 * t360 + t366 * t361;
t211 = -t359 * t619 + t360 * t494 + t361 * t495;
t732 = -t112 * t619 + t211 * t507;
t731 = t112 * t741 + t211 * t665;
t690 = t590 * t741;
t509 = t555 * t594 + t690;
t409 = -qJD(5) * t509 - t506 * t590 + t594 * t665;
t508 = -t555 * t590 + t594 * t741;
t600 = qJD(5) * t508 + t590 * t665;
t410 = t506 * t594 + t600;
t270 = Icges(6,5) * t410 + Icges(6,6) * t409 + Icges(6,3) * t507;
t271 = Icges(6,4) * t410 + Icges(6,2) * t409 + Icges(6,6) * t507;
t272 = Icges(6,1) * t410 + Icges(6,4) * t409 + Icges(6,5) * t507;
t394 = Icges(6,5) * t509 + Icges(6,6) * t508 - Icges(6,3) * t619;
t395 = Icges(6,4) * t509 + Icges(6,2) * t508 - Icges(6,6) * t619;
t396 = Icges(6,1) * t509 + Icges(6,4) * t508 - Icges(6,5) * t619;
t123 = -t270 * t619 + t508 * t271 + t509 * t272 + t507 * t394 + t409 * t395 + t410 * t396;
t215 = -t394 * t619 + t395 * t508 + t396 * t509;
t730 = -t123 * t619 + t215 * t507;
t729 = t123 * t741 + t215 * t665;
t631 = -t267 * rSges(7,1) - t266 * rSges(7,2);
t176 = t390 * rSges(7,3) - t631;
t630 = -t434 * rSges(7,1) - t433 * rSges(7,2);
t324 = -t513 * rSges(7,3) - t630;
t728 = t511 * t176 + t392 * t324;
t177 = t269 * rSges(7,1) + t268 * rSges(7,2) + t392 * rSges(7,3);
t323 = t432 * rSges(7,1) + t431 * rSges(7,2) + t511 * rSges(7,3);
t727 = -t177 * t619 + t507 * t323;
t388 = t390 * pkin(10);
t194 = -t390 * t597 + t391 * t760 - t388 - t800;
t726 = t176 + t194;
t303 = t393 * pkin(4) + pkin(10) * t392;
t614 = pkin(5) * t798 - t392 * t597 + t393 * t581;
t195 = -t303 + t614;
t725 = t177 + t195;
t416 = Icges(5,5) * t506 - Icges(5,6) * t507 + Icges(5,3) * t665;
t417 = Icges(5,4) * t506 - Icges(5,2) * t507 + Icges(5,6) * t665;
t418 = Icges(5,1) * t506 - Icges(5,4) * t507 + Icges(5,5) * t665;
t450 = Icges(5,5) * t555 + Icges(5,6) * t619 + Icges(5,3) * t741;
t451 = Icges(5,4) * t555 + Icges(5,2) * t619 + Icges(5,6) * t741;
t452 = Icges(5,1) * t555 + Icges(5,4) * t619 + Icges(5,5) * t741;
t183 = t416 * t741 + t417 * t619 + t555 * t418 + t450 * t665 - t507 * t451 + t506 * t452;
t291 = t450 * t741 + t451 * t619 + t452 * t555;
t724 = t183 * t741 + t291 * t665;
t446 = t514 * t594 + t744;
t279 = -qJD(5) * t446 - t391 * t590 + t489 * t594;
t280 = t391 * t594 - t780;
t190 = rSges(6,1) * t280 + rSges(6,2) * t279 + rSges(6,3) * t390;
t302 = t391 * pkin(4) + t388;
t723 = -t190 - t302;
t444 = t512 * t594 + t743;
t281 = -qJD(5) * t444 - t393 * t590 - t487 * t594;
t282 = t393 * t594 + t798;
t191 = t282 * rSges(6,1) + t281 * rSges(6,2) + t392 * rSges(6,3);
t722 = -t191 - t303;
t249 = rSges(7,1) * t366 + rSges(7,2) * t365 + rSges(7,3) * t507;
t287 = pkin(5) * t600 + t506 * t760 + t507 * t759;
t721 = t249 + t287;
t273 = rSges(6,1) * t410 + rSges(6,2) * t409 + rSges(6,3) * t507;
t426 = t506 * pkin(4) + t507 * pkin(10);
t720 = -t273 - t426;
t505 = t513 * pkin(10);
t428 = t514 * pkin(4) - t505;
t719 = t559 * t302 - t487 * t428;
t427 = t512 * pkin(4) + pkin(10) * t511;
t718 = t303 * t741 + t427 * t665;
t363 = rSges(7,1) * t495 + rSges(7,2) * t494 - rSges(7,3) * t619;
t717 = -t513 * t249 + t390 * t363;
t672 = pkin(5) * t743 - t511 * t597 + t512 * t581;
t309 = -t427 + t672;
t716 = t309 + t323;
t310 = pkin(5) * t744 + t513 * t597 + t514 * t760 + t505;
t715 = t310 + t324;
t333 = t444 * rSges(6,1) + t443 * rSges(6,2) + t511 * rSges(6,3);
t714 = -t333 - t427;
t334 = rSges(6,1) * t446 - rSges(6,2) * t628 - rSges(6,3) * t513;
t713 = -t334 - t428;
t705 = -t488 * qJ(3) - t556 * qJD(3);
t354 = -t705 + t762;
t485 = t489 * pkin(9);
t448 = pkin(3) * t669 + t485;
t712 = -t354 - t448;
t364 = pkin(5) * t690 + t555 * t760 - t619 * t759;
t711 = t363 + t364;
t397 = rSges(6,1) * t509 + rSges(6,2) * t508 - rSges(6,3) * t619;
t490 = t555 * pkin(4) - pkin(10) * t619;
t710 = t397 + t490;
t709 = t557 * t426 + t489 * t490;
t419 = rSges(5,1) * t506 - rSges(5,2) * t507 + rSges(5,3) * t665;
t517 = (-qJD(3) * t595 + (pkin(2) * t595 + qJ(3) * t592) * qJD(2)) * t588;
t708 = -t419 - t517;
t546 = t556 * qJ(3);
t491 = t546 + t761;
t492 = t559 * pkin(2) + qJ(3) * t558;
t707 = t491 * t740 + t492 * t738;
t477 = t589 * t492;
t518 = pkin(3) * t740 + pkin(9) * t559;
t706 = t589 * t518 + t477;
t519 = -pkin(3) * t738 + t557 * pkin(9);
t704 = -t491 - t519;
t703 = -t492 - t518;
t701 = -t517 - (-rSges(4,2) * t595 + rSges(4,3) * t592) * t695;
t533 = rSges(4,1) * t589 + (-rSges(4,2) * t592 - rSges(4,3) * t595) * t588;
t560 = (pkin(2) * t592 - qJ(3) * t595) * t588;
t700 = -t533 - t560;
t561 = pkin(3) * t589 + pkin(9) * t741;
t699 = -t560 - t561;
t698 = pkin(8) * t740 + t584;
t101 = -t270 * t513 - t271 * t628 + t272 * t446 + t279 * t395 + t280 * t396 + t390 * t394;
t184 = Icges(6,5) * t280 + Icges(6,6) * t279 + Icges(6,3) * t390;
t186 = Icges(6,4) * t280 + Icges(6,2) * t279 + Icges(6,6) * t390;
t188 = Icges(6,1) * t280 + Icges(6,4) * t279 + Icges(6,5) * t390;
t328 = Icges(6,5) * t446 - Icges(6,6) * t628 - Icges(6,3) * t513;
t330 = Icges(6,4) * t446 - Icges(6,2) * t628 - Icges(6,6) * t513;
t332 = Icges(6,1) * t446 - Icges(6,4) * t628 - Icges(6,5) * t513;
t63 = -t184 * t619 + t186 * t508 + t188 * t509 + t328 * t507 + t330 * t409 + t332 * t410;
t686 = -t101 / 0.2e1 - t63 / 0.2e1;
t102 = t270 * t511 + t271 * t443 + t272 * t444 + t281 * t395 + t282 * t396 + t392 * t394;
t185 = Icges(6,5) * t282 + Icges(6,6) * t281 + Icges(6,3) * t392;
t187 = Icges(6,4) * t282 + Icges(6,2) * t281 + Icges(6,6) * t392;
t189 = Icges(6,1) * t282 + Icges(6,4) * t281 + Icges(6,5) * t392;
t327 = Icges(6,5) * t444 + Icges(6,6) * t443 + Icges(6,3) * t511;
t329 = Icges(6,4) * t444 + Icges(6,2) * t443 + Icges(6,6) * t511;
t331 = Icges(6,1) * t444 + Icges(6,4) * t443 + Icges(6,5) * t511;
t62 = -t185 * t619 + t187 * t508 + t189 * t509 + t327 * t507 + t329 * t409 + t331 * t410;
t685 = t62 / 0.2e1 + t102 / 0.2e1;
t684 = -t426 - t721;
t683 = -t517 + t720;
t682 = -t302 + t712;
t681 = -t427 - t716;
t680 = -t428 - t715;
t353 = -t487 * pkin(2) - qJ(3) * t486 + qJD(3) * t558;
t679 = t353 * t738 + t354 * t740 + t491 * t668;
t678 = t490 + t711;
t259 = t393 * rSges(5,1) - t392 * rSges(5,2) - t487 * rSges(5,3);
t677 = t589 * t427 + t706;
t676 = -t427 + t703;
t675 = -t428 + t704;
t453 = rSges(5,1) * t555 + rSges(5,2) * t619 + rSges(5,3) * t741;
t674 = -t453 + t699;
t381 = -t487 * rSges(3,1) + t486 * rSges(3,2) + rSges(3,3) * t668;
t673 = -t490 + t699;
t404 = t512 * rSges(5,1) - t511 * rSges(5,2) + t559 * rSges(5,3);
t471 = t559 * rSges(3,1) - t558 * rSges(3,2) + rSges(3,3) * t740;
t380 = rSges(4,1) * t668 + t487 * rSges(4,2) - t486 * rSges(4,3);
t468 = rSges(4,1) * t740 - t559 * rSges(4,2) + t558 * rSges(4,3);
t667 = t585 * t693;
t165 = -t327 * t619 + t329 * t508 + t331 * t509;
t203 = t394 * t511 + t395 * t443 + t396 * t444;
t664 = t165 / 0.2e1 + t203 / 0.2e1;
t166 = -t328 * t619 + t330 * t508 + t332 * t509;
t204 = -t394 * t513 - t395 * t628 + t396 * t446;
t663 = t204 / 0.2e1 + t166 / 0.2e1;
t661 = t693 / 0.2e1;
t660 = -t593 * pkin(1) + pkin(8) * t738;
t449 = pkin(3) * t668 - pkin(9) * t487;
t659 = 2 * m(3);
t658 = 2 * m(4);
t656 = 2 * m(5);
t654 = 2 * m(6);
t652 = 0.2e1 * m(7);
t651 = t596 * t700;
t646 = pkin(9) * t667;
t644 = -t517 + t684;
t643 = -t397 + t673;
t642 = t518 * t738 + t519 * t740 + t707;
t639 = -t546 + t660;
t638 = -pkin(1) * t697 + pkin(8) * t668;
t637 = t674 * t596;
t26 = t732 + t748 + t749 - t755 + t756;
t150 = -t315 * t513 + t317 * t433 + t319 * t434;
t151 = -t316 * t513 + t318 * t433 + t320 * t434;
t197 = -t359 * t513 + t360 * t433 + t361 * t434;
t43 = -t171 * t513 + t173 * t433 + t175 * t434 + t266 * t317 + t267 * t319 + t315 * t390;
t44 = -t170 * t513 + t172 * t433 + t174 * t434 + t266 * t318 + t267 * t320 + t316 * t390;
t88 = -t246 * t513 + t247 * t433 + t248 * t434 + t266 * t360 + t267 * t361 + t359 * t390;
t7 = t150 * t392 + t151 * t390 + t197 * t507 + t43 * t511 - t44 * t513 - t619 * t88;
t148 = t315 * t511 + t317 * t431 + t319 * t432;
t149 = t316 * t511 + t318 * t431 + t320 * t432;
t196 = t359 * t511 + t360 * t431 + t361 * t432;
t75 = t148 * t511 - t149 * t513 - t196 * t619;
t76 = t150 * t511 - t151 * t513 - t197 * t619;
t45 = t171 * t511 + t173 * t431 + t175 * t432 + t268 * t317 + t269 * t319 + t315 * t392;
t46 = t170 * t511 + t172 * t431 + t174 * t432 + t268 * t318 + t269 * t320 + t316 * t392;
t89 = t246 * t511 + t247 * t431 + t248 * t432 + t268 * t360 + t269 * t361 + t359 * t392;
t8 = t148 * t392 + t149 * t390 + t196 * t507 + t45 * t511 - t46 * t513 - t619 * t89;
t93 = t162 * t511 - t163 * t513 - t211 * t619;
t635 = -t26 * t619 + t390 * t76 + t392 * t75 + t507 * t93 + t511 * t8 - t513 * t7;
t634 = -rSges(3,1) * t489 + rSges(3,2) * t488;
t633 = -rSges(5,1) * t391 + rSges(5,2) * t390;
t632 = -rSges(5,1) * t514 - rSges(5,2) * t513;
t629 = t673 - t711;
t627 = t643 * t596;
t626 = t492 + t698;
t625 = rSges(4,1) * t738 - rSges(4,3) * t556;
t624 = t448 * t740 + t449 * t738 + t519 * t668 + t679;
t623 = t427 * t738 + t428 * t740 + t642;
t621 = t639 - t519;
t620 = t629 * t596;
t351 = t589 * t353;
t618 = t589 * t449 - t593 * t646 + t351;
t522 = t560 * t669;
t617 = t561 * t669 - t596 * t646 + t522;
t470 = rSges(3,1) * t557 - rSges(3,2) * t556 - rSges(3,3) * t738;
t616 = t589 * t303 + t618;
t615 = t490 * t669 + t617;
t612 = t518 + t626;
t253 = Icges(5,5) * t393 - Icges(5,6) * t392 - Icges(5,3) * t487;
t255 = Icges(5,4) * t393 - Icges(5,2) * t392 - Icges(5,6) * t487;
t257 = Icges(5,1) * t393 - Icges(5,4) * t392 - Icges(5,5) * t487;
t398 = Icges(5,5) * t512 - Icges(5,6) * t511 + Icges(5,3) * t559;
t400 = Icges(5,4) * t512 - Icges(5,2) * t511 + Icges(5,6) * t559;
t402 = Icges(5,1) * t512 - Icges(5,4) * t511 + Icges(5,5) * t559;
t126 = t255 * t619 + t257 * t555 - t400 * t507 + t402 * t506 + (t253 * t592 + t398 * t693) * t588;
t145 = -t392 * t451 + t393 * t452 + t416 * t559 - t417 * t511 + t418 * t512 - t450 * t487;
t611 = t126 / 0.2e1 + t56 / 0.2e1 + t89 / 0.2e1 + t145 / 0.2e1 + t685;
t252 = Icges(5,5) * t391 - Icges(5,6) * t390 + Icges(5,3) * t489;
t254 = Icges(5,4) * t391 - Icges(5,2) * t390 + Icges(5,6) * t489;
t256 = Icges(5,1) * t391 - Icges(5,4) * t390 + Icges(5,5) * t489;
t399 = Icges(5,5) * t514 + Icges(5,6) * t513 + Icges(5,3) * t557;
t401 = Icges(5,4) * t514 + Icges(5,2) * t513 + Icges(5,6) * t557;
t403 = Icges(5,1) * t514 + Icges(5,4) * t513 + Icges(5,5) * t557;
t127 = t254 * t619 + t256 * t555 - t401 * t507 + t403 * t506 + (t252 * t592 + t399 * t693) * t588;
t144 = -t390 * t451 + t391 * t452 + t416 * t557 + t417 * t513 + t418 * t514 + t450 * t489;
t610 = t127 / 0.2e1 + t57 / 0.2e1 + t88 / 0.2e1 + t144 / 0.2e1 - t686;
t609 = t302 * t740 + t303 * t738 + t428 * t668 + t624;
t19 = t589 * t88 + (t43 * t593 - t44 * t596 + (t150 * t596 + t151 * t593) * qJD(1)) * t588;
t20 = t589 * t89 + (t45 * t593 - t46 * t596 + (t148 * t596 + t149 * t593) * qJD(1)) * t588;
t111 = t112 * t589;
t31 = t111 + (t56 * t593 - t57 * t596 + (t162 * t596 + t163 * t593) * qJD(1)) * t588;
t79 = t196 * t589 + (t148 * t593 - t149 * t596) * t588;
t80 = t197 * t589 + (t150 * t593 - t151 * t596) * t588;
t98 = t211 * t589 + (t162 * t593 - t163 * t596) * t588;
t608 = t20 * t773 + t19 * t772 + t26 * t767 + t31 * t771 + t8 * t740 / 0.2e1 - t7 * t738 / 0.2e1 + t80 * t778 + t79 * t777 + t98 * t774 + (t593 * t76 + t596 * t75) * qJD(1) * t768;
t228 = t399 * t741 + t401 * t619 + t403 * t555;
t251 = t450 * t557 + t451 * t513 + t452 * t514;
t607 = t228 / 0.2e1 + t163 / 0.2e1 + t251 / 0.2e1 + t197 / 0.2e1 + t663;
t227 = t398 * t741 + t400 * t619 + t402 * t555;
t250 = t450 * t559 - t451 * t511 + t452 * t512;
t606 = -t250 / 0.2e1 - t196 / 0.2e1 - t227 / 0.2e1 - t162 / 0.2e1 - t664;
t605 = t749 / 0.2e1 + t748 / 0.2e1 + t197 * t778 + t196 * t777 + t756 / 0.2e1 - t755 / 0.2e1 + t89 * t773 + t88 * t772 + t732;
t11 = -t150 * t487 + t151 * t489 + t43 * t559 + t44 * t557 + (t197 * t693 + t592 * t88) * t588;
t12 = -t148 * t487 + t149 * t489 + t45 * t559 + t46 * t557 + (t196 * t693 + t592 * t89) * t588;
t28 = -t162 * t487 + t163 * t489 + t57 * t557 + t56 * t559 + t731;
t77 = t148 * t559 + t149 * t557 + t196 * t741;
t78 = t150 * t559 + t151 * t557 + t197 * t741;
t97 = t162 * t559 + t163 * t557 + t211 * t741;
t604 = t11 * t772 + t12 * t773 + t7 * t770 + t8 * t769 + t26 * t741 / 0.2e1 + t28 * t771 + t78 * t778 + t77 * t777 + t75 * t776 + t76 * t775 + t588 * t93 * t661 + t97 * t774;
t602 = t353 + t638;
t601 = (-t584 + (-pkin(3) - pkin(8)) * t740) * qJD(1) - t485 + t705;
t599 = t449 + t602;
t598 = t601 - t762;
t544 = (rSges(3,1) * t595 - rSges(3,2) * t592) * t695;
t532 = rSges(3,3) * t589 + (rSges(3,1) * t592 + rSges(3,2) * t595) * t588;
t531 = Icges(4,1) * t589 + (-Icges(4,4) * t592 - Icges(4,5) * t595) * t588;
t526 = Icges(3,3) * t589 + (Icges(3,5) * t592 + Icges(3,6) * t595) * t588;
t469 = -rSges(4,2) * t557 - t625;
t441 = t557 * t490;
t436 = t471 + t698;
t435 = -t470 + t660;
t423 = t427 * t741;
t414 = -t470 * t589 - t532 * t738;
t413 = t471 * t589 - t532 * t740;
t408 = t559 * t428;
t405 = rSges(5,3) * t557 - t632;
t382 = rSges(3,3) * t669 - t634;
t379 = rSges(4,1) * t669 - rSges(4,2) * t489 + t758;
t356 = t626 + t468;
t355 = t557 * t765 + t625 + t639;
t350 = (-t584 + (-rSges(3,3) - pkin(8)) * t740) * qJD(1) + t634;
t349 = t638 + t381;
t346 = t526 * t740 - t527 * t558 + t528 * t559;
t345 = -t526 * t738 - t527 * t556 + t528 * t557;
t344 = t529 * t556 - t530 * t557 - t531 * t738;
t343 = t529 * t558 - t530 * t559 + t531 * t740;
t340 = t513 * t363;
t336 = (-t469 - t491) * t589 + t588 * t651;
t335 = t468 * t589 + t700 * t740 + t477;
t326 = t381 * t589 + (-t532 * t696 - t544 * t593) * t588;
t325 = -t382 * t589 + (t532 * t697 - t544 * t596) * t588;
t322 = t404 * t741 - t453 * t559;
t321 = -t405 * t741 + t453 * t557;
t314 = t460 * t589 + (-t456 * t595 - t458 * t592) * t588;
t313 = t459 * t589 + (-t455 * t595 - t457 * t592) * t588;
t312 = t462 * t589 + (t464 * t595 + t466 * t592) * t588;
t311 = t461 * t589 + (t463 * t595 + t465 * t592) * t588;
t308 = t612 + t404;
t307 = t557 * t764 + t621 + t632;
t304 = t619 * t323;
t286 = t511 * t324;
t275 = (t468 * t596 + t469 * t593) * t588 + t707;
t265 = -t404 * t557 + t405 * t559;
t262 = -t758 + t765 * t489 + (-t584 + (-rSges(4,1) - pkin(8)) * t740) * qJD(1) + t705;
t261 = t602 + t380;
t258 = rSges(5,3) * t489 - t633;
t240 = -t488 * t527 + t489 * t528 - t541 * t556 + t542 * t557 + (t526 * t697 - t539 * t596) * t588;
t239 = t486 * t527 - t487 * t528 - t541 * t558 + t542 * t559 + (t526 * t696 + t539 * t593) * t588;
t238 = -t486 * t529 + t487 * t530 + t537 * t558 - t538 * t559 + (t531 * t696 + t540 * t593) * t588;
t237 = t488 * t529 - t489 * t530 + t537 * t556 - t538 * t557 + (t531 * t697 - t540 * t596) * t588;
t236 = (-t405 + t704) * t589 + t588 * t637;
t235 = t404 * t589 + t674 * t740 + t706;
t232 = t334 * t619 - t397 * t513;
t231 = -t333 * t619 - t397 * t511;
t230 = t612 - t714;
t229 = t621 + t713 - t761;
t226 = t324 * t619 - t340;
t225 = -t363 * t511 - t304;
t224 = t380 * t589 + t351 + (qJD(1) * t651 + t593 * t701) * t588;
t223 = t522 + (-t354 - t379) * t589 + (t533 * t697 + t596 * t701) * t588;
t222 = t612 + t672 + t323;
t221 = -t514 * t581 + (-pkin(5) * t590 - pkin(2)) * t557 + t754 * t513 + t621 + t630;
t220 = t399 * t557 + t401 * t513 + t403 * t514;
t219 = t398 * t557 + t400 * t513 + t402 * t514;
t218 = t399 * t559 - t401 * t511 + t403 * t512;
t217 = t398 * t559 - t400 * t511 + t402 * t512;
t216 = (t404 * t596 + t405 * t593) * t588 + t642;
t213 = t333 * t513 + t334 * t511;
t210 = t323 * t513 + t286;
t208 = t489 * t764 + t601 + t633;
t207 = t599 + t259;
t206 = t333 * t741 - t559 * t710 + t423;
t205 = t397 * t557 + t713 * t741 + t441;
t202 = t372 * t589 + (-t368 * t595 - t370 * t592 + (t455 * t592 - t457 * t595) * qJD(2)) * t588;
t201 = t371 * t589 + (-t367 * t595 - t369 * t592 + (t456 * t592 - t458 * t595) * qJD(2)) * t588;
t200 = t373 * t589 + (t375 * t595 + t377 * t592 + (-t464 * t592 + t466 * t595) * qJD(2)) * t588;
t199 = t374 * t589 + (t376 * t595 + t378 * t592 + (-t463 * t592 + t465 * t595) * qJD(2)) * t588;
t193 = (-t334 + t675) * t589 + t588 * t627;
t192 = t333 * t589 + t643 * t740 + t677;
t182 = t183 * t589;
t180 = t334 * t559 + t557 * t714 + t408;
t179 = t419 * t557 + t453 * t489 + (-t258 * t592 - t405 * t693) * t588;
t178 = -t419 * t559 + t453 * t487 + (t259 * t592 + t404 * t693) * t588;
t161 = -t328 * t513 - t330 * t628 + t332 * t446;
t160 = -t327 * t513 - t329 * t628 + t331 * t446;
t159 = t328 * t511 + t330 * t443 + t332 * t444;
t158 = t327 * t511 + t329 * t443 + t331 * t444;
t157 = -t364 * t513 + t619 * t715 - t340;
t156 = -t309 * t619 - t511 * t711 - t304;
t155 = (t333 * t596 + t334 * t593) * t588 + t623;
t154 = (t379 * t593 + t380 * t596 + (t469 * t596 + (-t468 - t492) * t593) * qJD(1)) * t588 + t679;
t153 = t259 * t589 + (qJD(1) * t637 + t593 * t708) * t588 + t618;
t152 = (-t258 + t712) * t589 + (t453 * t697 + t596 * t708) * t588 + t617;
t143 = -t559 * t678 + t716 * t741 + t423;
t142 = t557 * t711 + t680 * t741 + t441;
t141 = t258 * t559 - t259 * t557 - t404 * t489 - t405 * t487;
t140 = (t675 - t715) * t589 + t588 * t620;
t139 = t589 * t716 + t629 * t740 + t677;
t138 = t310 * t511 + t513 * t716 + t286;
t137 = t598 + t723;
t136 = t599 - t722;
t135 = t557 * t681 + t559 * t715 + t408;
t134 = t251 * t589 + (t219 * t593 - t220 * t596) * t588;
t133 = t250 * t589 + (t217 * t593 - t218 * t596) * t588;
t132 = t219 * t559 + t220 * t557 + t251 * t741;
t131 = t217 * t559 + t218 * t557 + t250 * t741;
t130 = -t754 * t390 - t391 * t581 + t598 + t631 + t800;
t129 = t599 + t614 + t177;
t128 = (t593 * t715 + t596 * t716) * t588 + t623;
t125 = -t191 * t619 - t273 * t511 + t333 * t507 - t392 * t397;
t124 = t190 * t619 - t273 * t513 - t334 * t507 + t390 * t397;
t122 = (t258 * t593 + t259 * t596 + (t405 * t596 + (-t404 + t703) * t593) * qJD(1)) * t588 + t624;
t121 = t123 * t589;
t119 = -t249 * t511 - t363 * t392 + t727;
t118 = t176 * t619 - t324 * t507 + t717;
t116 = t252 * t559 - t254 * t511 + t256 * t512 - t392 * t401 + t393 * t403 - t399 * t487;
t115 = t253 * t559 - t255 * t511 + t257 * t512 - t392 * t400 + t393 * t402 - t398 * t487;
t114 = t252 * t557 + t254 * t513 + t256 * t514 - t390 * t401 + t391 * t403 + t399 * t489;
t113 = t253 * t557 + t255 * t513 + t257 * t514 - t390 * t400 + t391 * t402 + t398 * t489;
t108 = t191 * t589 + (qJD(1) * t627 + t593 * t683) * t588 + t616;
t107 = (-t190 + t682) * t589 + (t397 * t697 + t596 * t683) * t588 + t615;
t106 = t273 * t557 + t397 * t489 + (t592 * t723 + t693 * t713) * t588 + t709;
t105 = (t191 * t592 + t333 * t693) * t588 + t720 * t559 + t710 * t487 + t718;
t104 = t215 * t589 + (t165 * t593 - t166 * t596) * t588;
t103 = t165 * t559 + t166 * t557 + t215 * t741;
t100 = t165 * t511 - t166 * t513 - t215 * t619;
t99 = t190 * t511 + t191 * t513 - t333 * t390 + t334 * t392;
t95 = t177 * t513 - t323 * t390 + t728;
t91 = t204 * t589 + (t160 * t593 - t161 * t596) * t588;
t90 = t203 * t589 + (t158 * t593 - t159 * t596) * t588;
t86 = t160 * t559 + t161 * t557 + t204 * t741;
t85 = t158 * t559 + t159 * t557 + t203 * t741;
t82 = t160 * t511 - t161 * t513 - t204 * t619;
t81 = t158 * t511 - t159 * t513 - t203 * t619;
t64 = t190 * t559 - t334 * t487 + t489 * t714 + t557 * t722 + t719;
t61 = t725 * t589 + (qJD(1) * t620 + t593 * t644) * t588 + t616;
t60 = (t682 - t726) * t589 + (t596 * t644 + t697 * t711) * t588 + t615;
t59 = -t195 * t619 + t309 * t507 - t392 * t711 - t511 * t721 + t727;
t58 = -t287 * t513 + t364 * t390 - t507 * t715 + t619 * t726 + t717;
t55 = (t190 * t593 + t191 * t596 + (t334 * t596 + (-t333 + t676) * t593) * qJD(1)) * t588 + t609;
t52 = t721 * t557 + t711 * t489 + ((-t302 - t726) * t592 + t680 * t693) * t588 + t709;
t51 = t684 * t559 + t678 * t487 + (t592 * t725 + t693 * t716) * t588 + t718;
t50 = t184 * t511 + t186 * t443 + t188 * t444 + t281 * t330 + t282 * t332 + t328 * t392;
t49 = t185 * t511 + t187 * t443 + t189 * t444 + t281 * t329 + t282 * t331 + t327 * t392;
t48 = -t184 * t513 - t186 * t628 + t188 * t446 + t279 * t330 + t280 * t332 + t328 * t390;
t47 = -t185 * t513 - t187 * t628 + t189 * t446 + t279 * t329 + t280 * t331 + t327 * t390;
t42 = t194 * t511 + t310 * t392 - t390 * t716 + t513 * t725 + t728;
t41 = t182 + (t126 * t593 - t127 * t596 + (t227 * t596 + t228 * t593) * qJD(1)) * t588;
t40 = t726 * t559 - t715 * t487 + (-t303 - t725) * t557 + t681 * t489 + t719;
t39 = t126 * t559 + t127 * t557 - t227 * t487 + t228 * t489 + t724;
t38 = (t725 * t596 + t726 * t593 + (t715 * t596 + (t676 - t716) * t593) * qJD(1)) * t588 + t609;
t37 = t145 * t589 + (t115 * t593 - t116 * t596 + (t217 * t596 + t218 * t593) * qJD(1)) * t588;
t36 = t144 * t589 + (t113 * t593 - t114 * t596 + (t219 * t596 + t220 * t593) * qJD(1)) * t588;
t35 = t115 * t559 + t116 * t557 - t217 * t487 + t218 * t489 + (t145 * t592 + t250 * t693) * t588;
t34 = t113 * t559 + t114 * t557 - t219 * t487 + t220 * t489 + (t144 * t592 + t251 * t693) * t588;
t33 = t121 + (t62 * t593 - t63 * t596 + (t165 * t596 + t166 * t593) * qJD(1)) * t588;
t32 = -t165 * t487 + t166 * t489 + t63 * t557 + t62 * t559 + t729;
t29 = t165 * t392 + t166 * t390 + t62 * t511 - t63 * t513 + t730;
t22 = t102 * t589 + (t49 * t593 - t50 * t596 + (t158 * t596 + t159 * t593) * qJD(1)) * t588;
t21 = t101 * t589 + (t47 * t593 - t48 * t596 + (t160 * t596 + t161 * t593) * qJD(1)) * t588;
t18 = -t158 * t487 + t159 * t489 + t49 * t559 + t50 * t557 + (t102 * t592 + t203 * t693) * t588;
t17 = -t160 * t487 + t161 * t489 + t47 * t559 + t48 * t557 + (t101 * t592 + t204 * t693) * t588;
t14 = -t102 * t619 + t158 * t392 + t159 * t390 + t203 * t507 + t49 * t511 - t50 * t513;
t13 = -t101 * t619 + t160 * t392 + t161 * t390 + t204 * t507 + t47 * t511 - t48 * t513;
t1 = [t123 + t112 - t530 * t665 + t183 + (t129 * t222 + t130 * t221) * t652 + (t136 * t230 + t137 * t229) * t654 + (t207 * t308 + t208 * t307) * t656 + (t261 * t356 + t262 * t355) * t658 + (t349 * t436 + t350 * t435) * t659 - t537 * t739 - t588 * t737 + t799; t182 + t121 + t111 + (t152 * t307 + t153 * t308 + t207 * t235 + t208 * t236) * m(5) + (t107 * t229 + t108 * t230 + t136 * t192 + t137 * t193) * m(6) + (t129 * t139 + t130 * t140 + t221 * t60 + t222 * t61) * m(7) + ((-t240 / 0.2e1 - t199 / 0.2e1 - t201 / 0.2e1 - t237 / 0.2e1 - t610) * t596 + (t238 / 0.2e1 + t239 / 0.2e1 + t200 / 0.2e1 + t202 / 0.2e1 + t611) * t593 + ((t312 / 0.2e1 + t313 / 0.2e1 + t343 / 0.2e1 + t346 / 0.2e1 - t606) * t596 + (t344 / 0.2e1 + t345 / 0.2e1 + t311 / 0.2e1 + t314 / 0.2e1 + t607) * t593) * qJD(1)) * t588 + m(3) * (t325 * t435 + t326 * t436 + t349 * t413 + t350 * t414) + m(4) * (t223 * t355 + t224 * t356 + t261 * t335 + t262 * t336) + t793; (t128 * t38 + t139 * t61 + t140 * t60) * t652 + (t414 * t325 + t413 * t326 + (t470 * t593 + t471 * t596) * (t381 * t596 + t382 * t593 + (t470 * t596 - t471 * t593) * qJD(1)) * t585) * t659 + (t107 * t193 + t108 * t192 + t155 * t55) * t654 + (t122 * t216 + t152 * t236 + t153 * t235) * t656 + (t154 * t275 + t223 * t336 + t224 * t335) * t658 + (t37 + t22 + t20) * t740 + (-t36 - t21 - t19) * t738 + (t134 + t91 + t80) * t669 + (t133 + t90 + t79) * t668 + ((-t593 * t795 + t596 * t796) * t669 + (t593 * t797 - t596 * t794) * t668 + (t794 * t697 + t797 * t696 + (t486 * t785 + t487 * t783 + t558 * t790 - t559 * t788 - t668 * t786) * t596 + ((t593 * t792 + t696 * t787 + t781) * t588 + t789 * t559 - t791 * t558 + t782 * t487 + t784 * t486) * t593) * t740 + (t796 * t697 + t795 * t696 + ((t697 * t786 + t781) * t588 + t788 * t557 - t790 * t556 + t783 * t489 + t785 * t488) * t596 + ((t596 * t792 - t697 * t787) * t588 - t789 * t557 + t791 * t556 + t782 * t489 + t784 * t488) * t593) * t738) * t588 + (t31 + t33 + t41 + (t239 + t238) * t740 + (-t240 - t237) * t738 + (t345 + t344) * t669 + (t346 + t343) * t668 + ((-t199 - t201) * t596 + (t200 + t202) * t593 + ((t312 + t313) * t596 + (t311 + t314) * t593) * qJD(1)) * t588 + t793) * t589; (m(4) * t262 + m(5) * t208 + m(6) * t137 + m(7) * t130) * t558 + (m(4) * t261 + m(5) * t207 + m(6) * t136 + m(7) * t129) * t556 + (m(4) * t356 + m(5) * t308 + m(6) * t230 + m(7) * t222) * t488 + (-m(4) * t355 - m(5) * t307 - m(6) * t229 - m(7) * t221) * t486; (-m(4) * t154 - m(5) * t122 - m(6) * t55 - m(7) * t38) * t739 + (m(4) * t223 + m(5) * t152 + m(6) * t107 + m(7) * t60) * t558 + (m(4) * t224 + m(5) * t153 + m(6) * t108 + m(7) * t61) * t556 + (m(4) * t335 + m(5) * t235 + m(6) * t192 + m(7) * t139) * t488 + (-m(4) * t336 - m(5) * t236 - m(6) * t193 - m(7) * t140) * t486 + (m(4) * t275 + m(5) * t216 + m(6) * t155 + m(7) * t128) * t666; 0.4e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + t779) * (-t486 * t558 + t488 * t556 - t592 * t667); (t178 * t308 + t179 * t307 + t207 * t322 + t208 * t321) * m(5) + (t105 * t230 + t106 * t229 + t136 * t206 + t137 * t205) * m(6) + (t129 * t143 + t130 * t142 + t221 * t52 + t222 * t51) * m(7) + t611 * t559 + t610 * t557 + t607 * t489 + t606 * t487 + t724 + t729 + t731; (t122 * t265 + t141 * t216 + t152 * t321 + t153 * t322 + t178 * t235 + t179 * t236) * m(5) + (t105 * t192 + t106 * t193 + t107 * t205 + t108 * t206 + t155 * t64 + t180 * t55) * m(6) + (t128 * t40 + t135 * t38 + t139 * t51 + t140 * t52 + t142 * t60 + t143 * t61) * m(7) + (t39 / 0.2e1 + t32 / 0.2e1 + t28 / 0.2e1) * t589 + (t37 / 0.2e1 + t22 / 0.2e1 + t20 / 0.2e1) * t559 + (t36 / 0.2e1 + t21 / 0.2e1 + t19 / 0.2e1) * t557 + (t91 / 0.2e1 + t80 / 0.2e1 + t134 / 0.2e1) * t489 + (-t90 / 0.2e1 - t79 / 0.2e1 - t133 / 0.2e1) * t487 + ((-t11 / 0.2e1 - t17 / 0.2e1 - t34 / 0.2e1) * t596 + (t12 / 0.2e1 + t18 / 0.2e1 + t35 / 0.2e1) * t593 + (t31 / 0.2e1 + t33 / 0.2e1 + t41 / 0.2e1) * t592 + (t291 * t767 + t104 / 0.2e1 + t98 / 0.2e1 + (t227 * t593 - t228 * t596) * t768) * t693 + ((t131 / 0.2e1 + t77 / 0.2e1 + t85 / 0.2e1) * t596 + (t132 / 0.2e1 + t78 / 0.2e1 + t86 / 0.2e1) * t593) * qJD(1)) * t588; (-m(5) * t141 - m(6) * t64 - m(7) * t40) * t739 + (m(5) * t179 + m(6) * t106 + m(7) * t52) * t558 + (m(5) * t178 + m(6) * t105 + m(7) * t51) * t556 + (m(5) * t322 + m(6) * t206 + m(7) * t143) * t488 + (-m(5) * t321 - m(6) * t205 - m(7) * t142) * t486 + (m(5) * t265 + m(6) * t180 + m(7) * t135) * t666; (t135 * t40 + t142 * t52 + t143 * t51) * t652 + (t105 * t206 + t106 * t205 + t180 * t64) * t654 + (t141 * t265 + t178 * t322 + t179 * t321) * t656 + (t12 + t18 + t35) * t559 + (t11 + t17 + t34) * t557 + (t78 + t86 + t132) * t489 + (-t77 - t85 - t131) * t487 + ((t28 + t32 + t39) * t592 + (t227 * t559 + t228 * t557 + t291 * t741 + t103 + t97) * t693) * t588; (t124 * t229 + t125 * t230 + t136 * t231 + t137 * t232) * m(6) + (t129 * t156 + t130 * t157 + t221 * t58 + t222 * t59) * m(7) + t686 * t513 + t685 * t511 + t663 * t390 + t664 * t392 + t605 + t730; t608 + (-t596 * t13 / 0.2e1 + t14 * t766 + (t82 * t766 + t596 * t81 / 0.2e1) * qJD(1)) * t588 + (t107 * t232 + t108 * t231 + t124 * t193 + t125 * t192 + t155 * t99 + t213 * t55) * m(6) + (t128 * t42 + t138 * t38 + t139 * t59 + t140 * t58 + t156 * t61 + t157 * t60) * m(7) + t29 * t767 + t33 * t771 + t21 * t772 + t22 * t773 + t104 * t774 + t90 * t777 + t91 * t778; (-m(6) * t99 - m(7) * t42) * t739 + (m(6) * t124 + m(7) * t58) * t558 + (m(6) * t125 + m(7) * t59) * t556 + (m(6) * t231 + m(7) * t156) * t488 + (-m(6) * t232 - m(7) * t157) * t486 + (m(6) * t213 + m(7) * t138) * t666; (t100 * t661 + t592 * t29 / 0.2e1) * t588 + t604 + (t105 * t231 + t106 * t232 + t124 * t205 + t125 * t206 + t180 * t99 + t213 * t64) * m(6) + (t135 * t42 + t138 * t40 + t142 * t58 + t143 * t59 + t156 * t51 + t157 * t52) * m(7) + t14 * t769 + t13 * t770 + t32 * t771 + t17 * t772 + t18 * t773 + t103 * t774 + t82 * t775 + t81 * t776 + t85 * t777 + t86 * t778; (t138 * t42 + t156 * t59 + t157 * t58) * t652 + (t124 * t232 + t125 * t231 + t213 * t99) * t654 + t392 * t81 + t511 * t14 + t390 * t82 - t513 * t13 + t507 * t100 - t619 * t29 + t635; (t118 * t221 + t119 * t222 + t129 * t225 + t130 * t226) * m(7) + t605; (t118 * t140 + t119 * t139 + t128 * t95 + t210 * t38 + t225 * t61 + t226 * t60) * m(7) + t608; 0.2e1 * (t118 * t558 + t119 * t556 + t225 * t488 - t226 * t486 + (t210 * t694 - t595 * t95) * t588) * t779; t604 + (t118 * t142 + t119 * t143 + t135 * t95 + t210 * t40 + t225 * t51 + t226 * t52) * m(7); (t118 * t157 + t119 * t156 + t138 * t95 + t210 * t42 + t225 * t59 + t226 * t58) * m(7) + t635; (t118 * t226 + t119 * t225 + t210 * t95) * t652 + t635;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;