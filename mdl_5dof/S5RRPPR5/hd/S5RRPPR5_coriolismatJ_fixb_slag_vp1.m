% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:55
% EndTime: 2019-12-31 19:29:21
% DurationCPUTime: 17.03s
% Computational Cost: add. (30203->604), mult. (39395->791), div. (0->0), fcn. (42591->8), ass. (0->369)
t429 = sin(qJ(1));
t666 = -t429 / 0.4e1;
t425 = qJ(2) + pkin(8);
t407 = sin(t425);
t408 = cos(t425);
t657 = sin(qJ(5));
t658 = cos(qJ(5));
t344 = t407 * t658 - t408 * t657;
t447 = t407 * t657 + t408 * t658;
t240 = -Icges(6,5) * t447 - Icges(6,6) * t344;
t431 = cos(qJ(1));
t315 = t344 * t431;
t316 = t447 * t431;
t342 = Icges(6,4) * t344;
t244 = -Icges(6,2) * t447 + t342;
t248 = Icges(6,1) * t447 + t342;
t565 = t244 + t248;
t341 = Icges(6,4) * t447;
t245 = Icges(6,2) * t344 + t341;
t247 = Icges(6,1) * t344 - t341;
t567 = t245 - t247;
t82 = t240 * t429 + t315 * t567 + t316 * t565;
t286 = Icges(6,4) * t316;
t205 = Icges(6,2) * t315 - Icges(6,6) * t429 + t286;
t573 = -Icges(6,1) * t315 + t205 + t286;
t285 = Icges(6,4) * t315;
t207 = Icges(6,1) * t316 - Icges(6,5) * t429 + t285;
t701 = -Icges(6,2) * t316 + t207 + t285;
t91 = t344 * t573 + t447 * t701;
t784 = (t82 + t91) * t666;
t313 = t344 * t429;
t314 = t447 * t429;
t80 = -t240 * t431 + t313 * t567 + t314 * t565;
t783 = -t80 / 0.2e1;
t782 = -t80 / 0.4e1;
t780 = Icges(3,3) + Icges(4,3);
t779 = t91 / 0.2e1 + t82 / 0.2e1;
t284 = Icges(6,4) * t314;
t607 = Icges(6,2) * t313;
t478 = t607 + t284;
t604 = Icges(6,6) * t431;
t443 = t604 + t478;
t283 = Icges(6,4) * t313;
t616 = Icges(6,1) * t314;
t482 = t283 + t616;
t610 = Icges(6,5) * t431;
t444 = t610 + t482;
t466 = -t313 * t205 - t314 * t207;
t473 = Icges(6,5) * t314 + Icges(6,6) * t313;
t441 = Icges(6,3) * t431 + t473;
t707 = t441 - t473;
t427 = t431 ^ 2;
t540 = 0.2e1 * t284;
t477 = t540 + 0.2e1 * t604;
t539 = 0.2e1 * t610;
t693 = -t313 * (t477 + t607) - t314 * (t539 + t616) - Icges(6,3) * t427;
t763 = t693 * t431;
t15 = -t763 + (t707 * t429 + (-t444 + t482) * t316 - (t443 - t478) * t315 + t466) * t429;
t440 = Icges(6,5) * t316 + Icges(6,6) * t315 - Icges(6,3) * t429;
t439 = t431 * t440;
t108 = t439 - t466;
t109 = t315 * t205 + t316 * t207 - t429 * t440;
t16 = (-t315 * t478 - t316 * t482 + t108 - 0.2e1 * t439 + t466) * t431 + (-t313 * t443 - t314 * t444 - t707 * t431 + t109 - t693) * t429;
t54 = t108 * t429 + t763;
t55 = t109 * t429 - (t315 * t443 + t316 * t444 - t429 * t441) * t431;
t660 = t431 / 0.4e1;
t661 = -t431 / 0.4e1;
t778 = 0.2e1 * t16 * t660 + 0.2e1 * t55 * t661 + 0.2e1 * (t15 + t54) * t666 + (t82 / 0.4e1 + t91 / 0.4e1) * t429;
t250 = rSges(6,1) * t344 - rSges(6,2) * t447;
t600 = qJ(4) * t408;
t365 = pkin(3) * t407 - t600;
t428 = sin(qJ(2));
t635 = pkin(2) * t428;
t489 = pkin(4) * t407 + t635;
t454 = t250 + t365 + t489;
t181 = t454 * t429;
t183 = t454 * t431;
t628 = rSges(5,1) * t407;
t366 = -rSges(5,3) * t408 + t628;
t491 = t365 + t366 + t635;
t233 = t491 * t429;
t235 = t491 * t431;
t422 = t431 * rSges(5,2);
t430 = cos(qJ(2));
t634 = pkin(2) * t430;
t406 = pkin(1) + t634;
t633 = -qJ(3) - pkin(6);
t547 = -t429 * t406 - t431 * t633;
t621 = rSges(5,3) + qJ(4);
t659 = rSges(5,1) + pkin(3);
t696 = -t407 * t621 - t408 * t659;
t199 = t696 * t429 + t422 + t547;
t404 = t429 * t633;
t626 = t429 * rSges(5,2);
t200 = t626 - t404 + (t406 - t696) * t431;
t587 = t408 * t431;
t588 = t408 * t429;
t575 = t199 * t587 + t200 * t588;
t601 = qJ(4) * t407;
t686 = pkin(3) + pkin(4);
t697 = -t408 * t686 - t601;
t210 = t316 * rSges(6,1) + t315 * rSges(6,2) - t429 * rSges(6,3);
t706 = -t429 * pkin(7) + t210;
t155 = -t404 + (t406 - t697) * t431 + t706;
t705 = -t314 * rSges(6,1) - t313 * rSges(6,2);
t741 = (-rSges(6,3) - pkin(7)) * t431 + t697 * t429 + t547 + t705;
t578 = t155 * t588 + t587 * t741;
t591 = t407 * t431;
t592 = t407 * t429;
t687 = m(6) / 0.2e1;
t688 = m(5) / 0.2e1;
t620 = (-t233 * t591 + t235 * t592 + t575) * t688 + (-t181 * t591 + t183 * t592 + t578) * t687;
t392 = pkin(3) * t592;
t231 = t392 + (-t408 * t621 + t628 + t635) * t429;
t384 = qJ(4) * t587;
t391 = rSges(5,3) * t587;
t232 = t384 + t391 + (-t407 * t659 - t635) * t431;
t221 = -rSges(6,1) * t313 + rSges(6,2) * t314;
t178 = t392 + (t489 - t600) * t429 - t221;
t224 = t315 * rSges(6,1) - t316 * rSges(6,2);
t179 = t384 + (-t407 * t686 - t635) * t431 - t224;
t469 = t178 * t431 + t179 * t429;
t631 = (t407 * t469 + t578) * t687 + ((t231 * t431 + t232 * t429) * t407 + t575) * t688;
t736 = t620 - t631;
t777 = t736 * qJD(1);
t356 = Icges(4,5) * t408 - Icges(4,6) * t407;
t375 = Icges(3,5) * t430 - Icges(3,6) * t428;
t774 = (t356 + t375) * t431 + t780 * t429;
t402 = Icges(5,5) * t407;
t617 = Icges(5,1) * t408;
t483 = t402 + t617;
t296 = Icges(5,4) * t429 + t431 * t483;
t612 = Icges(4,4) * t407;
t363 = Icges(4,1) * t408 - t612;
t298 = Icges(4,5) * t429 + t363 * t431;
t773 = t296 + t298;
t584 = t429 * t430;
t586 = t428 * t429;
t762 = -Icges(3,5) * t584 - Icges(4,5) * t588 + Icges(3,6) * t586 + Icges(4,6) * t592 + t780 * t431;
t579 = Icges(6,1) - Icges(6,2);
t517 = t579 * t314;
t451 = t517 + t610;
t727 = -0.2e1 * t283;
t761 = t451 - t727;
t772 = t761 * t315 + t429 * (-Icges(6,5) * t313 + Icges(6,6) * t314);
t749 = -m(6) / 0.2e1;
t771 = -t248 / 0.2e1;
t208 = rSges(6,3) * t431 - t705;
t140 = -t429 * t208 - t431 * t210;
t753 = -t429 * t221 + t224 * t431;
t770 = t140 * t753;
t229 = t429 * t233;
t180 = t429 * t181;
t468 = t183 * t431 + t180;
t367 = rSges(4,1) * t407 + rSges(4,2) * t408;
t448 = t367 + t635;
t714 = t448 * t431;
t715 = t448 * t429;
t728 = t429 * t715 + t431 * t714;
t527 = t468 * t749 + (-t235 * t431 - t229) * t688 - m(4) * t728 / 0.2e1;
t689 = m(4) / 0.2e1;
t529 = (t429 * t178 - t179 * t431) * t687 + (t429 * t231 - t232 * t431) * t688 + t728 * t689;
t18 = t529 - t527;
t769 = t18 * qJD(1);
t215 = -Icges(6,5) * t315 + Icges(6,6) * t316;
t768 = (-t215 * t429 - t315 * t701 + t316 * t573) * t429;
t767 = (t215 * t431 - t313 * t701 + t314 * t573) * t429;
t766 = t753 * t408;
t765 = (t247 / 0.2e1 - t245 / 0.2e1) * t447;
t613 = Icges(3,4) * t428;
t379 = Icges(3,1) * t430 - t613;
t326 = Icges(3,5) * t429 + t379 * t431;
t764 = -t298 * t588 - t326 * t584;
t426 = t429 ^ 2;
t545 = t426 + t427;
t249 = -rSges(6,1) * t447 - rSges(6,2) * t344;
t710 = t545 * t407;
t639 = m(6) * (t249 * t710 + t766);
t760 = t774 * t431 + t764;
t358 = Icges(4,2) * t408 + t612;
t603 = Icges(5,3) * t408;
t474 = t603 - t402;
t759 = (-t358 - t474) * t431 + t773;
t388 = Icges(4,4) * t592;
t297 = Icges(4,1) * t588 - Icges(4,5) * t431 - t388;
t399 = Icges(3,4) * t586;
t325 = Icges(3,1) * t584 - Icges(3,5) * t431 - t399;
t583 = t430 * t431;
t758 = -t297 * t587 - t325 * t583 + t762 * t429;
t293 = Icges(4,4) * t588 - Icges(4,2) * t592 - Icges(4,6) * t431;
t323 = Icges(3,4) * t584 - Icges(3,2) * t586 - Icges(3,6) * t431;
t757 = t293 * t407 + t323 * t428;
t387 = Icges(5,5) * t587;
t288 = Icges(5,6) * t429 + Icges(5,3) * t591 + t387;
t357 = Icges(5,4) * t408 + Icges(5,6) * t407;
t292 = Icges(5,2) * t429 + t357 * t431;
t756 = t288 * t591 + t326 * t583 + t773 * t587 + (t292 + t774) * t429;
t755 = -Icges(3,5) * t428 - Icges(3,6) * t430 + (-Icges(4,6) + Icges(5,6)) * t408 + (-Icges(5,4) - Icges(4,5)) * t407;
t708 = t155 * t429 + t431 * t741;
t670 = t244 / 0.2e1;
t750 = t567 * t447 / 0.2e1 + (-t670 + t771) * t344;
t665 = t429 / 0.2e1;
t720 = -t431 / 0.2e1;
t744 = t249 * t708;
t403 = Icges(4,4) * t408;
t608 = Icges(4,2) * t407;
t294 = Icges(4,6) * t429 + (t403 - t608) * t431;
t418 = Icges(3,4) * t430;
t377 = -Icges(3,2) * t428 + t418;
t324 = Icges(3,6) * t429 + t377 * t431;
t739 = t294 * t407 + t324 * t428 + t762;
t516 = t579 * t313;
t738 = t516 - t604;
t737 = t545 * t408;
t585 = t428 * t431;
t735 = -t294 * t591 - t324 * t585 + t756;
t596 = (-Icges(5,2) * t431 + t429 * t357) * t431;
t734 = t596 + t756;
t732 = -t293 * t591 - t294 * t592 - t323 * t585 - t324 * t586 - t758 - t760;
t611 = Icges(5,5) * t408;
t355 = Icges(5,3) * t407 + t611;
t360 = Icges(5,1) * t407 - t611;
t376 = Icges(3,2) * t430 + t613;
t618 = Icges(4,1) * t407;
t731 = -(t379 / 0.2e1 - t376 / 0.2e1) * t428 - (t403 + t618 / 0.2e1 - t608 / 0.2e1 + t360 / 0.2e1 - t355 / 0.2e1) * t408 - (t363 / 0.2e1 - t358 / 0.2e1 + t402 + t617 / 0.2e1 - t603 / 0.2e1) * t407;
t726 = -0.2e1 * t710;
t721 = -t429 / 0.2e1;
t287 = -Icges(5,6) * t431 + t355 * t429;
t295 = -Icges(5,4) * t431 + t429 * t483;
t717 = (t287 * t407 + t295 * t408) * t429;
t538 = m(5) / 0.4e1 + m(6) / 0.4e1;
t548 = t737 * t407;
t711 = t538 * (-t407 * t408 + t548);
t504 = t545 * t250;
t378 = Icges(3,1) * t428 + t418;
t704 = t759 * t429;
t484 = -t403 - t618;
t512 = (t484 * t431 - t294) * t429;
t514 = (-Icges(5,1) * t591 + t288 + t387) * t429;
t703 = t512 + t514;
t702 = qJD(2) * t431;
t700 = t755 * t431;
t699 = t755 * t429;
t501 = -0.2e1 * t284;
t436 = t501 + t738;
t459 = (-t768 - (t436 * t316 + t772) * t431) * t665 + (-t767 - ((t501 - 0.2e1 * t604) * t314 - (-0.2e1 * t517 + t727 - 0.2e1 * t610) * t313) * t431) * t720;
t438 = t540 - t738;
t460 = (t768 - (t438 * t316 - t772) * t431) * t665 + (t767 - (-(-t727 + t539) * t313 + (t477 - 0.2e1 * t516) * t314) * t431) * t720;
t349 = t376 * t431;
t351 = t378 * t431;
t695 = (-t326 + t349) * t586 + (-t324 - t351) * t584;
t348 = -Icges(3,2) * t584 - t399;
t350 = t378 * t429;
t694 = (t325 + t348) * t428 + (t323 + t350) * t430;
t692 = 0.4e1 * qJD(1);
t691 = 0.2e1 * qJD(2);
t369 = rSges(5,1) * t408 + rSges(5,3) * t407;
t368 = pkin(3) * t408 + t601;
t424 = t431 * pkin(6);
t560 = -t429 * (pkin(1) * t429 - t424 + t547) + t431 * (-t429 * pkin(6) - t404 + (-pkin(1) + t406) * t431);
t496 = t545 * t368 + t560;
t131 = t429 * (t369 * t429 - t422) + (t369 * t431 + t626) * t431 + t496;
t569 = -t408 * t229 - t235 * t587;
t683 = m(5) * (t131 * t710 + t569);
t679 = m(6) * (-t181 * t224 - t183 * t221 - t744);
t678 = m(6) * (t250 * t469 - t744);
t576 = -t408 * t180 - t183 * t587;
t95 = (pkin(4) * t587 + t706) * t431 + (pkin(4) * t588 + pkin(7) * t431 + t208) * t429 + t496;
t677 = m(6) * (t710 * t95 + t576);
t503 = t545 * t249;
t568 = t737 * t250;
t675 = m(6) * (-t766 + (t140 - t503) * t407 + t568);
t672 = m(6) * (t155 * t179 + t178 * t741);
t671 = m(6) * (t155 * t224 + t221 * t741);
t630 = rSges(3,1) * t430;
t523 = pkin(1) + t630;
t546 = rSges(3,2) * t586 + t431 * rSges(3,3);
t267 = -t429 * t523 + t424 + t546;
t401 = rSges(3,2) * t585;
t268 = -t401 + t523 * t431 + (rSges(3,3) + pkin(6)) * t429;
t380 = rSges(3,1) * t428 + rSges(3,2) * t430;
t352 = t380 * t429;
t353 = t380 * t431;
t656 = m(3) * (t267 * t352 - t268 * t353);
t462 = rSges(4,1) * t588 - rSges(4,2) * t592 - t431 * rSges(4,3);
t238 = -t462 + t547;
t518 = -rSges(4,2) * t591 + t429 * rSges(4,3);
t629 = rSges(4,1) * t408;
t239 = -t404 + (t406 + t629) * t431 + t518;
t655 = m(4) * (t238 * t715 - t239 * t714);
t654 = m(4) * (t238 * t431 + t239 * t429);
t648 = m(5) * (t199 * t231 + t200 * t232);
t647 = m(5) * (-t199 * t592 + t200 * t591);
t646 = m(5) * (t199 * t431 + t200 * t429);
t642 = m(6) * (t140 * t710 + t568);
t641 = m(6) * (t155 * t591 - t592 * t741);
t640 = m(6) * t708;
t452 = (t221 * t431 + t224 * t429) * t407;
t98 = t452 * t749;
t580 = t98 * qJD(4);
t557 = -t360 * t429 + t287;
t556 = -t484 * t429 + t293;
t555 = -t474 * t429 + t295;
t553 = -Icges(4,2) * t588 + t297 - t388;
t446 = t753 * t749;
t495 = m(6) * t504;
t115 = t446 - t495 / 0.2e1;
t543 = t115 * qJD(1);
t537 = t687 + t688;
t172 = t537 * t726;
t542 = t172 * qJD(1);
t535 = t15 / 0.2e1 + t54 / 0.2e1;
t534 = t55 / 0.2e1 - t16 / 0.2e1;
t521 = -t368 - t634;
t520 = rSges(4,2) * t407 - t629 - t634;
t519 = -t460 - t459;
t515 = t557 * t431;
t513 = t556 * t431;
t511 = t555 * t431;
t509 = t553 * t431;
t494 = t545 * t635;
t88 = t436 * t344 - (t451 + 0.2e1 * t283) * t447;
t493 = t679 / 0.2e1 + t784 + (-t80 + t88) * t661;
t90 = t438 * t344 + t447 * t761;
t492 = t678 / 0.2e1 + t784 + (t80 + t90) * t660;
t490 = -t369 + t521;
t488 = -t288 * t592 + t292 * t431 - t296 * t588;
t487 = t356 / 0.2e1 + t357 / 0.2e1 + t375 / 0.2e1;
t453 = -pkin(4) * t408 + t249 + t521;
t182 = t453 * t429;
t184 = t453 * t431;
t467 = t182 * t429 + t184 * t431;
t449 = t429 * (qJ(4) * t588 - t392) + t431 * (-pkin(3) * t591 + t384) - t494;
t437 = t429 * t535 + t431 * t534;
t159 = -t596 + t717;
t435 = t534 + (t159 - t717 + t734) * t721 + t735 * t665 + ((t757 + t774) * t431 + t732 + t758 + t764) * t720;
t434 = t488 * t721 + t535 + (t159 + t762 * t431 + (t297 * t408 + t325 * t430 - t757) * t429) * t720 + (t739 * t431 - t734 + t735) * t431 / 0.2e1 + (t287 * t591 + t295 * t587 + t739 * t429 + t488 + t732 + t760) * t665;
t382 = -rSges(3,2) * t428 + t630;
t281 = t520 * t431;
t279 = t520 * t429;
t236 = t490 * t431;
t234 = t490 * t429;
t174 = 0.4e1 * t711;
t171 = t538 * t726 + (m(5) + m(6)) * t710 / 0.2e1;
t149 = t431 * (-rSges(5,1) * t591 + t391) - t366 * t426 + t449;
t120 = -pkin(4) * t710 + t449 - t753;
t114 = t446 + t495 / 0.2e1;
t110 = -t639 / 0.2e1;
t103 = t642 / 0.2e1;
t99 = t452 * t687;
t69 = t641 + t647;
t67 = -t250 * t503 + t770;
t61 = t675 / 0.2e1;
t56 = t640 + t646 + t654;
t45 = t671 + (t771 - t244 / 0.2e1) * t344 - t765;
t27 = t677 + t683;
t25 = t103 + t61 + t639 / 0.2e1;
t24 = t110 + t103 - t675 / 0.2e1;
t23 = t110 + t61 - t642 / 0.2e1;
t19 = t527 + t529;
t17 = (t378 / 0.2e1 + t377 / 0.2e1) * t430 + (t248 / 0.2e1 + t670) * t344 + t765 + t656 + t655 + t648 + t672 - t731;
t10 = m(6) * t67 + t460;
t8 = t620 + t631;
t6 = m(6) * (t249 * t468 + t753 * t95) + t459;
t4 = t429 * t434 + t431 * t435;
t3 = (t782 + t88 / 0.4e1) * t431 - t679 / 0.2e1 + t492 + t778;
t2 = t437 + t492 + t493;
t1 = (t782 - t90 / 0.4e1) * t431 - t678 / 0.2e1 + t493 + t778;
t5 = [t17 * qJD(2) + t56 * qJD(3) + t69 * qJD(4) + t45 * qJD(5), t17 * qJD(1) + t19 * qJD(3) + t8 * qJD(4) + t2 * qJD(5) + ((t238 * t281 + t239 * t279) * t689 + (t199 * t236 + t200 * t234 - t231 * t235 - t232 * t233) * t688 + (t155 * t182 - t178 * t183 - t179 * t181 + t184 * t741) * t687) * t691 + (m(3) * (-t267 * t382 - t352 * t380) - t90 / 0.2e1 + t783 + t487 * t431 + (-t325 / 0.2e1 - t348 / 0.2e1) * t430 + (t323 / 0.2e1 + t350 / 0.2e1) * t428 - t435) * t702 + ((m(3) * (-t268 * t382 + t353 * t380) + t487 * t429 + (t326 / 0.2e1 - t349 / 0.2e1) * t430 + (-t324 / 0.2e1 - t351 / 0.2e1) * t428 - t434 + t779) * t429 + (-t515 / 0.2e1 + t513 / 0.2e1 + t514 / 0.2e1 + t512 / 0.2e1) * t407 + (-t511 / 0.2e1 - t509 / 0.2e1 + t759 * t665) * t408) * qJD(2), qJD(1) * t56 + qJD(2) * t19 + qJD(4) * t171 + qJD(5) * t114, qJD(1) * t69 + qJD(2) * t8 + qJD(3) * t171 - qJD(5) * t98, t45 * qJD(1) + t2 * qJD(2) + t114 * qJD(3) - t580 + ((m(6) * (t221 * t250 + t249 * t741) + t783 + t88 / 0.2e1 - t534) * t431 + (m(6) * (t155 * t249 + t224 * t250) - t535 + t779) * t429) * qJD(5); t4 * qJD(2) - t18 * qJD(3) + t736 * qJD(4) + t1 * qJD(5) + (-t656 / 0.4e1 - t655 / 0.4e1 - t648 / 0.4e1 - t672 / 0.4e1) * t692 + (-(t378 + t377) * t430 / 0.2e1 + t731 + t750) * qJD(1), t4 * qJD(1) + t27 * qJD(4) + t6 * qJD(5) - (((t513 - t515 + t703) * t408 - t704 * t407 + ((t553 + t555) * t407 + t694 - t700) * t431 + t695) * t429 + t699 * t427) * t702 / 0.2e1 + (m(3) * ((t429 * (rSges(3,1) * t584 - t546) + t431 * (rSges(3,1) * t583 + t429 * rSges(3,3) - t401)) * (-t429 * t352 - t353 * t431) + t545 * t382 * t380) + m(6) * (t120 * t95 - t181 * t182 - t183 * t184) + m(5) * (t131 * t149 - t233 * t234 - t235 * t236) + m(4) * (-t714 * t281 - t715 * t279 + (t429 * t462 + t431 * (rSges(4,1) * t587 + t518) + t560) * (-t545 * t367 - t494)) + t460 + ((t694 * t431 + (t509 + t511 - t704) * t407 - t699 * t429 + ((t556 - t557) * t431 + t703) * t408 + t695) * t431 + t700 * t426) * t665) * qJD(2), -t769, t777 + t27 * qJD(2) + t24 * qJD(5) + (-0.4e1 * t711 + 0.2e1 * t537 * (-t408 * t710 + t548)) * qJD(4), t1 * qJD(1) + t6 * qJD(2) + t24 * qJD(4) + ((-t67 + (t140 - t95) * t753 + (-t468 - t504) * t249) * m(6) - t459 + t519) * qJD(5); t18 * qJD(2) + t172 * qJD(4) + t115 * qJD(5) + (-t640 / 0.4e1 - t646 / 0.4e1 - t654 / 0.4e1) * t692, t769 + ((-t182 * t431 + t429 * t184) * t687 + (-t234 * t431 + t429 * t236) * t688 + (-t279 * t431 + t429 * t281) * t689) * t691, 0, t542, t543; -t736 * qJD(2) - t172 * qJD(3) + t99 * qJD(5) + (-t641 / 0.4e1 - t647 / 0.4e1) * t692, -t777 + t174 * qJD(4) + t23 * qJD(5) + 0.4e1 * (-t677 / 0.4e1 - t683 / 0.4e1) * qJD(2) + ((-t408 * t120 + t576) * t687 + (-t408 * t149 + t569) * t688 + ((t467 + t95) * t687 + (t234 * t429 + t236 * t431 + t131) * t688) * t407) * t691, -t542, t174 * qJD(2), t99 * qJD(1) + t23 * qJD(2) + qJD(5) * t639; (-t671 - t750) * qJD(1) + t3 * qJD(2) - t115 * qJD(3) + t580 + t437 * qJD(5), t3 * qJD(1) + ((t140 * t120 + t250 * t467) * m(6) - t460 + t519) * qJD(2) + t25 * qJD(4) + t10 * qJD(5), -t543, qJD(1) * t98 + qJD(2) * t25, t437 * qJD(1) + t10 * qJD(2) + (m(6) * (t249 * t504 - t770) + t459) * qJD(5);];
Cq = t5;
