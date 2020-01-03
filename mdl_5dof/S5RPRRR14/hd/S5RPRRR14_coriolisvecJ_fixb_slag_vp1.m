% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR14_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:54
% EndTime: 2019-12-31 19:18:49
% DurationCPUTime: 92.08s
% Computational Cost: add. (105874->1396), mult. (310956->1886), div. (0->0), fcn. (380248->14), ass. (0->504)
t673 = sin(pkin(11));
t678 = cos(pkin(5));
t602 = t678 * t673;
t676 = cos(pkin(11));
t694 = sin(qJ(1));
t697 = cos(qJ(1));
t477 = t602 * t697 + t676 * t694;
t509 = sin(qJ(3));
t696 = cos(qJ(3));
t603 = t678 * t676;
t538 = -t603 * t697 + t673 * t694;
t674 = sin(pkin(6));
t675 = sin(pkin(5));
t600 = t675 * t674;
t677 = cos(pkin(6));
t784 = t538 * t677 + t697 * t600;
t746 = t784 * t696;
t423 = t477 * t509 + t746;
t601 = t677 * t675;
t454 = t538 * t674 - t697 * t601;
t756 = qJD(3) * t454;
t381 = qJD(4) * t423 + t756;
t734 = t784 * t509;
t426 = -t477 * t696 + t734;
t508 = sin(qJ(4));
t695 = cos(qJ(4));
t749 = t426 * t508 + t454 * t695;
t257 = -qJD(5) * t749 + t381;
t478 = -t602 * t694 + t676 * t697;
t539 = t603 * t694 + t673 * t697;
t747 = -t539 * t677 + t694 * t600;
t733 = t747 * t696;
t427 = t478 * t509 - t733;
t535 = t539 * t674;
t568 = t694 * t601;
t455 = t535 + t568;
t449 = qJD(3) * t455;
t382 = qJD(4) * t427 + t449;
t730 = t747 * t509;
t428 = t478 * t696 + t730;
t576 = -t428 * t508 + t455 * t695;
t258 = -qJD(5) * t576 + t382;
t476 = -t600 * t676 + t677 * t678;
t464 = qJD(3) * t476 + qJD(1);
t599 = t675 * t673;
t732 = t676 * t601 + t674 * t678;
t520 = -t509 * t599 + t696 * t732;
t429 = -qJD(4) * t520 + t464;
t452 = t509 * t732 + t696 * t599;
t575 = -t452 * t508 + t476 * t695;
t350 = -qJD(5) * t575 + t429;
t377 = t426 * t695 - t454 * t508;
t507 = sin(qJ(5));
t510 = cos(qJ(5));
t280 = t377 * t507 + t423 * t510;
t283 = t377 * t510 - t423 * t507;
t133 = Icges(6,5) * t283 - Icges(6,6) * t280 + Icges(6,3) * t749;
t666 = Icges(6,4) * t283;
t136 = -Icges(6,2) * t280 + Icges(6,6) * t749 + t666;
t276 = Icges(6,4) * t280;
t139 = Icges(6,1) * t283 + Icges(6,5) * t749 - t276;
t48 = t133 * t749 - t136 * t280 + t139 * t283;
t379 = t428 * t695 + t455 * t508;
t284 = -t379 * t507 + t427 * t510;
t285 = t379 * t510 + t427 * t507;
t134 = Icges(6,5) * t285 + Icges(6,6) * t284 - Icges(6,3) * t576;
t665 = Icges(6,4) * t285;
t137 = Icges(6,2) * t284 - Icges(6,6) * t576 + t665;
t277 = Icges(6,4) * t284;
t140 = Icges(6,1) * t285 - Icges(6,5) * t576 + t277;
t49 = -t134 * t749 + t280 * t137 - t140 * t283;
t422 = t452 * t695 + t476 * t508;
t366 = -t422 * t507 - t510 * t520;
t367 = t422 * t510 - t507 * t520;
t204 = Icges(6,5) * t367 + Icges(6,6) * t366 - Icges(6,3) * t575;
t664 = Icges(6,4) * t367;
t205 = Icges(6,2) * t366 - Icges(6,6) * t575 + t664;
t365 = Icges(6,4) * t366;
t206 = Icges(6,1) * t367 - Icges(6,5) * t575 + t365;
t67 = -t204 * t749 + t205 * t280 - t206 * t283;
t16 = t257 * t48 + t258 * t49 + t350 * t67;
t50 = t133 * t576 - t136 * t284 - t139 * t285;
t51 = -t134 * t576 + t284 * t137 + t285 * t140;
t68 = -t204 * t576 + t205 * t284 + t206 * t285;
t17 = t257 * t50 + t258 * t51 + t68 * t350;
t54 = t133 * t575 - t136 * t366 - t139 * t367;
t596 = rSges(6,1) * t283 - t280 * rSges(6,2);
t141 = -rSges(6,3) * t749 - t596;
t207 = rSges(6,1) * t367 + rSges(6,2) * t366 - rSges(6,3) * t575;
t692 = pkin(4) * t377;
t249 = pkin(10) * t749 + t692;
t336 = pkin(4) * t422 - pkin(10) * t575;
t806 = -t141 * t350 + t207 * t257 + t249 * t429 + t336 * t381;
t385 = Icges(4,5) * t452 + Icges(4,6) * t520 + Icges(4,3) * t476;
t670 = Icges(4,4) * t452;
t386 = Icges(4,2) * t520 + Icges(4,6) * t476 + t670;
t443 = Icges(4,4) * t520;
t387 = Icges(4,1) * t452 + Icges(4,5) * t476 + t443;
t145 = t385 * t455 - t386 * t427 + t387 * t428;
t293 = Icges(4,5) * t426 + Icges(4,6) * t423 - Icges(4,3) * t454;
t672 = Icges(4,4) * t426;
t296 = Icges(4,2) * t423 - Icges(4,6) * t454 + t672;
t413 = Icges(4,4) * t423;
t299 = Icges(4,1) * t426 - Icges(4,5) * t454 + t413;
t102 = -t293 * t455 + t296 * t427 - t299 * t428;
t294 = Icges(4,5) * t428 - Icges(4,6) * t427 + Icges(4,3) * t455;
t671 = Icges(4,4) * t428;
t297 = -Icges(4,2) * t427 + Icges(4,6) * t455 + t671;
t414 = Icges(4,4) * t427;
t300 = Icges(4,1) * t428 + Icges(4,5) * t455 - t414;
t103 = t455 * t294 - t427 * t297 + t428 * t300;
t590 = t102 * t454 + t103 * t455;
t47 = qJD(3) * t590 + t145 * t464;
t209 = Icges(5,5) * t377 - Icges(5,6) * t749 - Icges(5,3) * t423;
t669 = Icges(5,4) * t377;
t212 = -Icges(5,2) * t749 - Icges(5,6) * t423 + t669;
t368 = Icges(5,4) * t749;
t215 = Icges(5,1) * t377 - Icges(5,5) * t423 - t368;
t72 = -t209 * t423 - t212 * t749 + t215 * t377;
t210 = Icges(5,5) * t379 + Icges(5,6) * t576 + Icges(5,3) * t427;
t668 = Icges(5,4) * t379;
t213 = Icges(5,2) * t576 + Icges(5,6) * t427 + t668;
t369 = Icges(5,4) * t576;
t216 = Icges(5,1) * t379 + Icges(5,5) * t427 + t369;
t73 = t423 * t210 + t213 * t749 - t216 * t377;
t286 = Icges(5,5) * t422 + Icges(5,6) * t575 - Icges(5,3) * t520;
t667 = Icges(5,4) * t422;
t287 = Icges(5,2) * t575 - Icges(5,6) * t520 + t667;
t409 = Icges(5,4) * t575;
t288 = Icges(5,1) * t422 - Icges(5,5) * t520 + t409;
t783 = -t286 * t423 - t287 * t749 + t288 * t377;
t30 = t381 * t72 + t382 * t73 - t783 * t429;
t74 = -t209 * t427 - t212 * t576 - t215 * t379;
t75 = t427 * t210 + t213 * t576 + t379 * t216;
t87 = t286 * t427 + t287 * t576 + t288 * t379;
t31 = t381 * t74 + t382 * t75 + t87 * t429;
t79 = t209 * t520 - t212 * t575 - t215 * t422;
t101 = t454 * t294 - t423 * t297 - t300 * t426;
t767 = t293 * t454 - t296 * t423 - t299 * t426;
t591 = t101 * t455 - t767 * t454;
t218 = rSges(5,1) * t377 - rSges(5,2) * t749 - t423 * rSges(5,3);
t290 = rSges(5,1) * t422 + rSges(5,2) * t575 - rSges(5,3) * t520;
t798 = t218 * t429 + t290 * t381;
t122 = -t293 * t476 - t296 * t520 - t299 * t452;
t302 = rSges(4,1) * t426 + t423 * rSges(4,2) - t454 * rSges(4,3);
t388 = rSges(4,1) * t452 + rSges(4,2) * t520 + rSges(4,3) * t476;
t789 = t302 * t464 + t388 * t756;
t346 = -pkin(3) * t426 + t423 * pkin(9);
t407 = pkin(3) * t452 - pkin(9) * t520;
t788 = -t346 * t464 + t407 * t756;
t622 = t756 / 0.2e1;
t766 = -t385 * t454 + t386 * t423 + t387 * t426;
t757 = t454 * pkin(8);
t748 = -t477 * pkin(2) - t757;
t470 = t477 * qJD(1);
t643 = qJD(3) * t509;
t320 = qJD(1) * t734 + qJD(3) * t733 - t470 * t696 - t478 * t643;
t440 = t454 * qJD(1);
t200 = qJD(4) * t379 + t320 * t508 + t440 * t695;
t319 = -qJD(1) * t746 + qJD(3) * t428 - t470 * t509;
t435 = qJD(3) * t440;
t274 = qJD(4) * t319 - t435;
t130 = qJD(5) * t200 + t274;
t471 = t478 * qJD(1);
t322 = qJD(1) * t730 - qJD(3) * t746 + t471 * t696 - t477 * t643;
t441 = t455 * qJD(1);
t202 = -qJD(4) * t377 + t322 * t508 - t441 * t695;
t321 = -qJD(1) * t733 - qJD(3) * t426 + t471 * t509;
t436 = qJD(3) * t441;
t275 = qJD(4) * t321 + t436;
t131 = qJD(5) * t202 + t275;
t438 = t520 * qJD(3);
t363 = qJD(4) * t575 + t438 * t695;
t439 = t452 * qJD(3);
t222 = -qJD(5) * t367 - t363 * t507 + t439 * t510;
t223 = qJD(5) * t366 + t363 * t510 + t439 * t507;
t362 = qJD(4) * t422 + t438 * t508;
t104 = Icges(6,5) * t223 + Icges(6,6) * t222 + Icges(6,3) * t362;
t105 = Icges(6,4) * t223 + Icges(6,2) * t222 + Icges(6,6) * t362;
t106 = Icges(6,1) * t223 + Icges(6,4) * t222 + Icges(6,5) * t362;
t201 = qJD(4) * t576 + t320 * t695 - t440 * t508;
t118 = -qJD(5) * t285 - t201 * t507 + t319 * t510;
t119 = qJD(5) * t284 + t201 * t510 + t319 * t507;
t19 = -t104 * t576 + t105 * t284 + t106 * t285 + t118 * t205 + t119 * t206 + t200 * t204;
t640 = qJD(4) * t439;
t308 = qJD(5) * t362 + t640;
t203 = qJD(4) * t749 + t322 * t695 + t441 * t508;
t120 = qJD(5) * t283 - t203 * t507 + t321 * t510;
t121 = qJD(5) * t280 + t203 * t510 + t321 * t507;
t57 = Icges(6,5) * t121 + Icges(6,6) * t120 + Icges(6,3) * t202;
t59 = Icges(6,4) * t121 + Icges(6,2) * t120 + Icges(6,6) * t202;
t61 = Icges(6,1) * t121 + Icges(6,4) * t120 + Icges(6,5) * t202;
t8 = -t118 * t136 - t119 * t139 - t133 * t200 + t284 * t59 + t285 * t61 - t57 * t576;
t56 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t200;
t58 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t200;
t60 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t200;
t9 = t118 * t137 + t119 * t140 + t134 * t200 + t284 * t58 + t285 * t60 - t56 * t576;
t1 = t130 * t51 + t131 * t50 + t19 * t350 + t257 * t8 + t258 * t9 + t308 * t68;
t89 = Icges(5,5) * t203 - Icges(5,6) * t202 + Icges(5,3) * t321;
t91 = Icges(5,4) * t203 - Icges(5,2) * t202 + Icges(5,6) * t321;
t93 = Icges(5,1) * t203 - Icges(5,4) * t202 + Icges(5,5) * t321;
t21 = t200 * t212 - t201 * t215 - t209 * t319 + t379 * t93 + t427 * t89 + t576 * t91;
t88 = Icges(5,5) * t201 - Icges(5,6) * t200 + Icges(5,3) * t319;
t90 = Icges(5,4) * t201 - Icges(5,2) * t200 + Icges(5,6) * t319;
t92 = Icges(5,1) * t201 - Icges(5,4) * t200 + Icges(5,5) * t319;
t22 = -t200 * t213 + t201 * t216 + t210 * t319 + t379 * t92 + t427 * t88 + t576 * t90;
t224 = Icges(5,5) * t363 - Icges(5,6) * t362 + Icges(5,3) * t439;
t225 = Icges(5,4) * t363 - Icges(5,2) * t362 + Icges(5,6) * t439;
t226 = Icges(5,1) * t363 - Icges(5,4) * t362 + Icges(5,5) * t439;
t35 = -t200 * t287 + t201 * t288 + t224 * t427 + t225 * t576 + t226 * t379 + t286 * t319;
t745 = t21 * t381 + t382 * t22 + t274 * t75 + t275 * t74 + t35 * t429 + t640 * t87 + t1;
t10 = -t120 * t136 - t121 * t139 - t133 * t202 + t280 * t59 - t283 * t61 - t57 * t749;
t11 = t120 * t137 + t121 * t140 + t134 * t202 + t280 * t58 - t283 * t60 - t56 * t749;
t20 = -t104 * t749 + t105 * t280 - t106 * t283 + t120 * t205 + t121 * t206 + t202 * t204;
t2 = t257 * t10 + t11 * t258 + t130 * t49 + t131 * t48 + t20 * t350 + t308 * t67;
t23 = t202 * t212 - t203 * t215 - t209 * t321 - t377 * t93 + t423 * t89 + t749 * t91;
t24 = -t202 * t213 + t203 * t216 + t210 * t321 - t377 * t92 + t423 * t88 + t749 * t90;
t36 = -t202 * t287 + t203 * t288 + t224 * t423 + t225 * t749 - t226 * t377 + t286 * t321;
t744 = t381 * t23 + t24 * t382 + t274 * t73 + t275 * t72 + t36 * t429 - t640 * t783 + t2;
t55 = -t134 * t575 + t137 * t366 + t140 * t367;
t683 = t55 * t130;
t684 = t54 * t131;
t13 = t134 * t362 + t137 * t222 + t140 * t223 + t366 * t58 + t367 * t60 - t56 * t575;
t687 = t13 * t258;
t12 = -t133 * t362 - t136 * t222 - t139 * t223 + t366 * t59 + t367 * t61 - t57 * t575;
t688 = t12 * t257;
t26 = -t104 * t575 + t105 * t366 + t106 * t367 + t204 * t362 + t205 * t222 + t206 * t223;
t70 = -t204 * t575 + t205 * t366 + t206 * t367;
t690 = t26 * t350 + t70 * t308;
t3 = t683 + t684 + t687 + t688 + t690;
t80 = -t210 * t520 + t213 * t575 + t216 * t422;
t681 = t80 * t274;
t682 = t79 * t275;
t29 = t210 * t439 - t213 * t362 + t216 * t363 + t422 * t92 - t520 * t88 + t575 * t90;
t685 = t29 * t382;
t28 = -t209 * t439 + t212 * t362 - t215 * t363 + t422 * t93 - t520 * t89 + t575 * t91;
t686 = t28 * t381;
t53 = -t224 * t520 + t225 * t575 + t226 * t422 + t286 * t439 - t287 * t362 + t288 * t363;
t99 = -t286 * t520 + t287 * t575 + t288 * t422;
t689 = t53 * t429 + t99 * t640;
t743 = t681 + t682 + t685 + t686 + t689 + t3;
t742 = t16 + t30;
t741 = t31 + t17;
t490 = pkin(8) * t568;
t506 = t697 * pkin(1);
t607 = t675 * t694;
t646 = qJ(2) * t607 + t506;
t739 = -t490 - t646;
t537 = t538 * rSges(3,2);
t608 = t697 * t675;
t738 = -t477 * rSges(3,1) + rSges(3,3) * t608 + t537;
t657 = t141 - t249;
t143 = t285 * rSges(6,1) + t284 * rSges(6,2) - rSges(6,3) * t576;
t251 = t379 * pkin(4) - pkin(10) * t576;
t656 = t143 + t251;
t654 = t207 + t336;
t530 = pkin(8) * t535;
t521 = t478 * pkin(2) + t530 - t739;
t633 = t694 * pkin(1);
t483 = -qJ(2) * t608 + t633;
t735 = -t483 + t738;
t473 = t539 * rSges(3,2);
t731 = t478 * rSges(3,1) + rSges(3,3) * t607 - t473 + t646;
t729 = -t17 / 0.2e1;
t728 = -t31 / 0.2e1;
t727 = t130 / 0.2e1;
t726 = t131 / 0.2e1;
t725 = -t257 / 0.2e1;
t724 = t257 / 0.2e1;
t723 = -t258 / 0.2e1;
t722 = t258 / 0.2e1;
t721 = t274 / 0.2e1;
t720 = t275 / 0.2e1;
t719 = t308 / 0.2e1;
t716 = -t350 / 0.2e1;
t715 = t350 / 0.2e1;
t714 = -t381 / 0.2e1;
t713 = t381 / 0.2e1;
t712 = -t382 / 0.2e1;
t711 = t382 / 0.2e1;
t708 = -t429 / 0.2e1;
t707 = t429 / 0.2e1;
t698 = -rSges(6,3) - pkin(10);
t693 = t203 * pkin(4);
t124 = t201 * pkin(4) + t200 * pkin(10);
t62 = t119 * rSges(6,1) + t118 * rSges(6,2) + t200 * rSges(6,3);
t680 = t124 + t62;
t125 = t202 * pkin(10) + t693;
t597 = -t121 * rSges(6,1) - t120 * rSges(6,2);
t63 = t202 * rSges(6,3) - t597;
t679 = t125 + t63;
t349 = t428 * pkin(3) + t427 * pkin(9);
t663 = t349 * t441;
t662 = t407 * t440;
t661 = t423 * t508;
t660 = t427 * t508;
t659 = t520 * t508;
t107 = rSges(6,1) * t223 + rSges(6,2) * t222 + rSges(6,3) * t362;
t234 = pkin(4) * t363 + pkin(10) * t362;
t658 = t107 + t234;
t221 = t322 * pkin(3) + t321 * pkin(9);
t655 = t455 * t221 - t440 * t346;
t219 = t379 * rSges(5,1) + rSges(5,2) * t576 + t427 * rSges(5,3);
t653 = -t219 - t349;
t652 = -t290 - t407;
t395 = pkin(3) * t438 + pkin(9) * t439;
t651 = t454 * t395 + t441 * t407;
t650 = Icges(4,1) * t520 - t386 - t670;
t649 = Icges(4,2) * t452 - t387 - t443;
t559 = -t483 + t748;
t499 = qJD(2) * t607;
t582 = qJD(1) * t608;
t571 = qJ(2) * t582 - qJD(1) * t633 + t499;
t648 = (t499 + t571) * qJD(1);
t647 = -qJD(1) * t483 + t499;
t642 = qJD(4) * t426;
t641 = qJD(4) * t428;
t639 = qJD(4) * t452;
t638 = qJD(5) * t508;
t637 = 2 * m(3);
t636 = 2 * m(4);
t635 = 2 * m(5);
t634 = 2 * m(6);
t94 = t201 * rSges(5,1) - t200 * rSges(5,2) + t319 * rSges(5,3);
t190 = t320 * rSges(4,1) - t319 * rSges(4,2) - t440 * rSges(4,3);
t523 = -t470 * pkin(2) - qJD(1) * t757;
t632 = qJD(1) * t523 + t648;
t303 = t428 * rSges(4,1) - t427 * rSges(4,2) + t455 * rSges(4,3);
t631 = qJD(1) * t748 + t647;
t630 = -t470 * rSges(3,1) + rSges(3,3) * t582 + qJD(1) * t537;
t627 = t507 * t695;
t626 = t510 * t695;
t625 = -t435 / 0.2e1;
t624 = t436 / 0.2e1;
t623 = -t756 / 0.2e1;
t620 = t449 / 0.2e1;
t619 = t640 / 0.2e1;
t618 = t471 * pkin(2) + qJD(1) * t530;
t220 = t320 * pkin(3) + t319 * pkin(9);
t616 = -t636 / 0.2e1;
t615 = t636 / 0.2e1;
t614 = -t635 / 0.2e1;
t613 = t635 / 0.2e1;
t612 = -t634 / 0.2e1;
t611 = t634 / 0.2e1;
t345 = -pkin(3) * t423 - pkin(9) * t426;
t406 = pkin(3) * t520 + pkin(9) * t452;
t609 = -t345 * t464 + t406 * t756;
t348 = -pkin(3) * t427 + pkin(9) * t428;
t605 = t345 * t449 - t348 * t756;
t604 = t464 * t348 - t406 * t449;
t595 = -rSges(6,1) * t510 + rSges(6,2) * t507;
t594 = -Icges(6,1) * t510 + Icges(6,4) * t507;
t593 = -Icges(6,4) * t510 + Icges(6,2) * t507;
t592 = -Icges(6,5) * t510 + Icges(6,6) * t507;
t589 = -t302 * t455 - t303 * t454;
t588 = (-Icges(4,5) * t423 + Icges(4,6) * t426) * t454 + (-Icges(4,5) * t427 - Icges(4,6) * t428) * t455;
t394 = rSges(4,1) * t438 - rSges(4,2) * t439;
t587 = t388 * t440 - t394 * t455;
t586 = t388 * t441 + t394 * t454;
t500 = qJD(2) * t608;
t460 = qJD(1) * t646 - t500;
t496 = qJD(1) * t500;
t585 = t496 + (-qJD(1) * t490 - t460 - t618) * qJD(1);
t584 = qJD(1) * t559 + t499;
t583 = qJD(1) * t521 - t500;
t581 = qJD(1) * t607;
t580 = -pkin(4) * t695 - pkin(10) * t508;
t505 = qJD(2) * t678;
t579 = t346 * t449 - t349 * t756 + t505;
t578 = -rSges(5,1) * t695 + rSges(5,2) * t508;
t574 = -Icges(5,1) * t695 + Icges(5,4) * t508;
t573 = -Icges(5,4) * t695 + Icges(5,2) * t508;
t572 = -Icges(5,5) * t695 + Icges(5,6) * t508;
t191 = t322 * rSges(4,1) - t321 * rSges(4,2) + t441 * rSges(4,3);
t95 = t203 * rSges(5,1) - t202 * rSges(5,2) + t321 * rSges(5,3);
t564 = (Icges(6,5) * t280 + Icges(6,6) * t283) * t257 + (Icges(6,5) * t284 - Icges(6,6) * t285) * t258 + (Icges(6,5) * t366 - Icges(6,6) * t367) * t350;
t563 = (Icges(5,5) * t749 + Icges(5,6) * t377) * t381 + (Icges(5,5) * t576 - Icges(5,6) * t379) * t382 + (Icges(5,5) * t575 - Icges(5,6) * t422) * t429;
t562 = t631 + t788;
t561 = (-Icges(4,1) * t427 - t297 - t671) * t455 + (-Icges(4,1) * t423 + t296 + t672) * t454;
t560 = (Icges(4,2) * t428 - t300 + t414) * t455 + (-Icges(4,2) * t426 + t299 + t413) * t454;
t557 = t584 + t788;
t556 = t349 * t464 - t407 * t449 + t583;
t553 = t221 * t449 - t346 * t435 + (-t220 * t454 - t663) * qJD(3);
t551 = -t221 * t464 + t395 * t756 + t407 * t436 + t585;
t184 = Icges(4,5) * t320 - Icges(4,6) * t319 - Icges(4,3) * t440;
t185 = Icges(4,5) * t322 - Icges(4,6) * t321 + Icges(4,3) * t441;
t186 = Icges(4,4) * t320 - Icges(4,2) * t319 - Icges(4,6) * t440;
t187 = Icges(4,4) * t322 - Icges(4,2) * t321 + Icges(4,6) * t441;
t188 = Icges(4,1) * t320 - Icges(4,4) * t319 - Icges(4,5) * t440;
t189 = Icges(4,1) * t322 - Icges(4,4) * t321 + Icges(4,5) * t441;
t550 = -t767 * t441 - t101 * t440 + (t185 * t454 - t187 * t423 - t189 * t426 - t293 * t441 + t296 * t321 - t299 * t322) * t454 + (t184 * t454 - t186 * t423 - t188 * t426 + t294 * t441 - t297 * t321 + t300 * t322) * t455;
t549 = t102 * t441 - t103 * t440 + (t185 * t455 - t187 * t427 + t189 * t428 + t293 * t440 + t296 * t319 - t299 * t320) * t454 + (t184 * t455 - t186 * t427 + t188 * t428 - t294 * t440 - t297 * t319 + t300 * t320) * t455;
t123 = t294 * t476 + t297 * t520 + t300 * t452;
t64 = t185 * t476 + t187 * t520 + t189 * t452 + t296 * t439 - t299 * t438;
t65 = t184 * t476 + t186 * t520 + t188 * t452 - t297 * t439 + t300 * t438;
t548 = t122 * t441 - t123 * t440 + t64 * t454 + t65 * t455;
t547 = -t190 * t454 + t191 * t455 + t302 * t440 - t303 * t441;
t546 = t464 * t220 + (-t395 * t455 + t662) * qJD(3) + t632;
t545 = -t346 + t559;
t543 = (Icges(6,1) * t284 - t137 - t665) * t258 + (Icges(6,1) * t280 + t136 + t666) * t257 + (Icges(6,1) * t366 - t205 - t664) * t350;
t542 = (-Icges(6,2) * t285 + t140 + t277) * t258 + (Icges(6,2) * t283 - t139 + t276) * t257 + (-Icges(6,2) * t367 + t206 + t365) * t350;
t541 = (Icges(5,1) * t576 - t213 - t668) * t382 + (Icges(5,1) * t749 + t212 + t669) * t381 + (Icges(5,1) * t575 - t287 - t667) * t429;
t540 = (Icges(5,2) * t379 - t216 - t369) * t382 + (-Icges(5,2) * t377 + t215 - t368) * t381 + (Icges(5,2) * t422 - t288 - t409) * t429;
t519 = t349 + t521;
t518 = qJD(1) * t739 + t500 - t618;
t517 = (Icges(6,3) * t379 + t137 * t507 - t140 * t510 - t576 * t592) * t258 + (-Icges(6,3) * t377 - t136 * t507 + t139 * t510 - t592 * t749) * t257 + (Icges(6,3) * t422 + t205 * t507 - t206 * t510 - t575 * t592) * t350;
t516 = t471 * rSges(3,1) + rSges(3,3) * t581 - qJD(1) * t473;
t515 = t523 + t571;
t514 = -t221 + t518;
t513 = (Icges(5,3) * t428 + t213 * t508 - t216 * t695 + t427 * t572) * t382 + (-Icges(5,3) * t426 - t212 * t508 + t215 * t695 + t423 * t572) * t381 + (Icges(5,3) * t452 + t287 * t508 - t288 * t695 - t520 * t572) * t429;
t512 = t220 + t515;
t405 = rSges(4,1) * t520 - rSges(4,2) * t452;
t402 = Icges(4,5) * t520 - Icges(4,6) * t452;
t401 = qJD(1) * t731 - t500;
t400 = qJD(1) * t735 + t499;
t399 = t520 * t638 + t639;
t398 = t580 * t520;
t397 = t452 * t507 + t520 * t626;
t396 = t452 * t510 - t520 * t627;
t393 = Icges(4,1) * t438 - Icges(4,4) * t439;
t392 = Icges(4,4) * t438 - Icges(4,2) * t439;
t391 = Icges(4,5) * t438 - Icges(4,6) * t439;
t384 = t496 + (-t460 - t516) * qJD(1);
t383 = qJD(1) * t630 + t648;
t364 = t454 * t407;
t356 = t452 * rSges(5,3) - t520 * t578;
t354 = Icges(5,5) * t452 - t520 * t574;
t353 = Icges(5,6) * t452 - t520 * t573;
t344 = -rSges(4,1) * t427 - rSges(4,2) * t428;
t343 = -rSges(4,1) * t423 + rSges(4,2) * t426;
t335 = pkin(4) * t575 + pkin(10) * t422;
t334 = rSges(5,1) * t575 - rSges(5,2) * t422;
t330 = -t427 * t638 + t641;
t329 = -t423 * t638 - t642;
t328 = t580 * t427;
t327 = t580 * t423;
t326 = -t427 * t626 + t428 * t507;
t325 = t427 * t627 + t428 * t510;
t324 = -t423 * t626 - t426 * t507;
t323 = t423 * t627 - t426 * t510;
t311 = t476 * t349;
t307 = t455 * t346;
t273 = t428 * rSges(5,3) + t427 * t578;
t272 = -rSges(5,3) * t426 + t423 * t578;
t271 = Icges(5,5) * t428 + t427 * t574;
t270 = -Icges(5,5) * t426 + t423 * t574;
t269 = Icges(5,6) * t428 + t427 * t573;
t268 = -Icges(5,6) * t426 + t423 * t573;
t262 = rSges(6,3) * t422 - t575 * t595;
t261 = Icges(6,5) * t422 - t575 * t594;
t260 = Icges(6,6) * t422 - t575 * t593;
t255 = rSges(6,1) * t397 + rSges(6,2) * t396 + rSges(6,3) * t659;
t254 = Icges(6,1) * t397 + Icges(6,4) * t396 + Icges(6,5) * t659;
t253 = Icges(6,4) * t397 + Icges(6,2) * t396 + Icges(6,6) * t659;
t252 = Icges(6,5) * t397 + Icges(6,6) * t396 + Icges(6,3) * t659;
t250 = pkin(4) * t576 + pkin(10) * t379;
t247 = pkin(4) * t749 - pkin(10) * t377;
t246 = rSges(5,1) * t576 - rSges(5,2) * t379;
t245 = rSges(5,1) * t749 + rSges(5,2) * t377;
t238 = rSges(6,1) * t366 - rSges(6,2) * t367;
t227 = rSges(5,1) * t363 - rSges(5,2) * t362 + rSges(5,3) * t439;
t195 = t476 * t220;
t183 = rSges(6,3) * t379 - t576 * t595;
t182 = -rSges(6,3) * t377 - t595 * t749;
t181 = Icges(6,5) * t379 - t576 * t594;
t180 = -Icges(6,5) * t377 - t594 * t749;
t179 = Icges(6,6) * t379 - t576 * t593;
t178 = -Icges(6,6) * t377 - t593 * t749;
t172 = rSges(6,1) * t326 + rSges(6,2) * t325 - rSges(6,3) * t660;
t171 = rSges(6,1) * t324 + rSges(6,2) * t323 - rSges(6,3) * t661;
t170 = Icges(6,1) * t326 + Icges(6,4) * t325 - Icges(6,5) * t660;
t169 = Icges(6,1) * t324 + Icges(6,4) * t323 - Icges(6,5) * t661;
t168 = Icges(6,4) * t326 + Icges(6,2) * t325 - Icges(6,6) * t660;
t167 = Icges(6,4) * t324 + Icges(6,2) * t323 - Icges(6,6) * t661;
t166 = Icges(6,5) * t326 + Icges(6,6) * t325 - Icges(6,3) * t660;
t165 = Icges(6,5) * t324 + Icges(6,6) * t323 - Icges(6,3) * t661;
t164 = rSges(6,1) * t284 - rSges(6,2) * t285;
t163 = rSges(6,1) * t280 + rSges(6,2) * t283;
t150 = t303 * t464 - t388 * t449 + t583;
t149 = t584 + t789;
t148 = qJD(3) * t589 + t505;
t115 = -t386 * t439 + t387 * t438 + t391 * t476 + t392 * t520 + t393 * t452;
t108 = t115 * t464;
t98 = qJD(3) * t586 - t191 * t464 + t585;
t97 = qJD(3) * t587 + t190 * t464 + t632;
t84 = t219 * t429 - t290 * t382 + t556;
t83 = t557 + t798;
t78 = -t321 * t386 + t322 * t387 + t385 * t441 + t391 * t454 - t392 * t423 - t393 * t426;
t77 = -t319 * t386 + t320 * t387 - t385 * t440 + t391 * t455 - t392 * t427 + t393 * t428;
t76 = -t218 * t382 - t219 * t381 + t579;
t71 = t547 * qJD(3);
t41 = t143 * t350 - t207 * t258 + t251 * t429 - t336 * t382 + t556;
t40 = t557 + t806;
t37 = t141 * t258 - t143 * t257 - t249 * t382 - t251 * t381 + t579;
t34 = t218 * t640 + t227 * t381 + t275 * t290 - t429 * t95 + t551;
t33 = t219 * t640 - t227 * t382 - t274 * t290 + t429 * t94 + t546;
t32 = t381 * t79 + t382 * t80 + t429 * t99;
t27 = -t218 * t274 - t219 * t275 - t381 * t94 + t382 * t95 + t553;
t18 = t257 * t54 + t258 * t55 + t350 * t70;
t15 = t107 * t257 - t125 * t429 + t131 * t207 - t141 * t308 + t234 * t381 + t249 * t640 + t275 * t336 - t350 * t63 + t551;
t14 = -t107 * t258 + t124 * t429 - t130 * t207 + t143 * t308 - t234 * t382 + t251 * t640 - t274 * t336 + t350 * t62 + t546;
t7 = -t124 * t381 + t125 * t382 + t130 * t141 - t131 * t143 - t249 * t274 - t251 * t275 - t257 * t62 + t258 * t63 + t553;
t4 = [t47 * t623 + (t65 + t77) * t620 + t17 * t725 + (-t40 + t562 + t806) * t41 * t612 + (t122 - t766) * t624 + (t98 * (t302 + t559) + t149 * (-t191 + t518) + t97 * (t521 + t303) + t150 * (t515 + t190)) * t615 - t381 * t728 + (t34 * (t218 + t545) + t83 * (t514 - t95) + t33 * (t519 + t219) + t84 * (t512 + t94)) * t613 + (t123 + t145) * t625 + (t562 + t798 - t83) * t84 * t614 + t87 * t721 + t19 * t722 + t20 * t724 + t67 * t726 + t35 * t711 + t36 * t713 + t681 / 0.2e1 + t682 / 0.2e1 - t257 * t729 + (-t149 + t631 + t789) * t150 * t616 - t783 * t720 + t686 / 0.2e1 + t108 + t683 / 0.2e1 + t684 / 0.2e1 + t31 * t714 + t68 * t727 + t689 + t690 + t687 / 0.2e1 + t685 / 0.2e1 + (t15 * (-t698 * t749 + t545 + t596 + t692) + t40 * (t202 * t698 + t514 + t597 - t693) + t14 * (t519 + t656) + t41 * (t512 + t680)) * t611 + (t47 + t64 + t78) * t622 + t688 / 0.2e1 - (qJD(1) * t738 - t400 + t647) * t401 * t637 / 0.2e1 + (t384 * t735 + t400 * (-qJ(2) * t581 - qJD(1) * t506 + t500 - t516) + t383 * t731 + t401 * (t571 + t630)) * t637 / 0.2e1; (-t383 * t608 + t384 * t607) * m(3) + (-t14 * t608 + t15 * t607 + t678 * t7) * m(6) + (t27 * t678 - t33 * t608 + t34 * t607) * m(5) + (t607 * t98 - t608 * t97 + t678 * t71) * m(4); -((t455 * t402 + t427 * t649 + t428 * t650) * t464 + (t427 * t560 + t428 * t561 + t455 * t588) * qJD(3)) * t449 / 0.2e1 + (t145 * t476 + t590) * t625 + (-t476 * t766 + t591) * t624 + (t98 * (t302 * t476 + t388 * t454) + t149 * (-t191 * t476 + t586) + t97 * (t303 * t476 - t388 * t455) + t150 * (t190 * t476 + t587) + t71 * t589 + t148 * t547) * t615 - t329 * t16 / 0.2e1 + (qJD(3) * t548 + t108 + t743) * t476 / 0.2e1 + (qJD(3) * t550 + t464 * t78 + t744) * t454 / 0.2e1 + (qJD(3) * t549 + t464 * t77 + t745) * t455 / 0.2e1 + (t476 * t77 + t549) * t620 + (t476 * t78 + t550) * t622 + (t454 * t54 + t455 * t55 + t476 * t70) * t719 + (t454 * t74 + t455 * t75 + t476 * t87) * t721 + (t454 * t48 + t455 * t49 + t476 * t67) * t726 + (t454 * t50 + t455 * t51 + t476 * t68) * t727 - ((t476 * t402 + t452 * t650 - t520 * t649) * t464 + (t452 * t561 + t476 * t588 - t520 * t560) * qJD(3)) * t464 / 0.2e1 + (t40 * (-t141 * t399 - t171 * t350 + t207 * t329 + t255 * t257 - t327 * t429 - t381 * t398 + (t249 * t452 - t336 * t426) * qJD(4) + t609) + t41 * (t143 * t399 + t172 * t350 - t207 * t330 - t255 * t258 + t328 * t429 + t382 * t398 + (t251 * t452 - t336 * t428) * qJD(4) + t604) + t37 * (t141 * t330 - t143 * t329 + t171 * t258 - t172 * t257 + t327 * t382 - t328 * t381 + (-t249 * t428 + t251 * t426) * qJD(4) + t605)) * t612 + (t149 * (-t343 * t464 + t405 * t756) + t150 * (t344 * t464 - t405 * t449) + (t343 * t455 - t344 * t454) * qJD(3) * t148) * t616 + (t454 * t72 + t455 * t73 - t476 * t783) * t720 + ((t428 * t210 + t269 * t576 + t379 * t271) * t382 + (-t209 * t428 + t268 * t576 + t379 * t270) * t381 + (t428 * t286 + t353 * t576 + t379 * t354) * t429 + (-t426 * t74 + t428 * t75 + t452 * t87) * qJD(4) + t513 * t427) * t712 + ((t452 * t210 + t269 * t575 + t422 * t271) * t382 + (-t209 * t452 + t268 * t575 + t422 * t270) * t381 + (t452 * t286 + t353 * t575 + t422 * t354) * t429 + (-t426 * t79 + t428 * t80 + t452 * t99) * qJD(4) - t513 * t520) * t708 + ((-t210 * t426 + t269 * t749 - t271 * t377) * t382 + (t209 * t426 + t268 * t749 - t270 * t377) * t381 + (-t286 * t426 + t353 * t749 - t354 * t377) * t429 + (-t426 * t72 + t428 * t73 - t452 * t783) * qJD(4) + t513 * t423) * t714 + (t83 * (-t272 * t429 + t356 * t381 + (t218 * t452 - t290 * t426) * qJD(4) + t609) + t84 * (t273 * t429 - t356 * t382 + (t219 * t452 - t290 * t428) * qJD(4) + t604) + t76 * (t272 * t382 - t273 * t381 + (-t218 * t428 + t219 * t426) * qJD(4) + t605)) * t614 + (t34 * (t290 * t454 + t364 + (t218 - t346) * t476) + t83 * (t227 * t454 + t290 * t441 + (-t221 - t95) * t476 + t651) + t33 * (t219 * t476 + t455 * t652 + t311) + t84 * (t476 * t94 + t195 + (-t227 - t395) * t455 - t652 * t440) + t27 * (-t218 * t455 + t454 * t653 + t307) + t76 * (t218 * t440 + t455 * t95 + (-t220 - t94) * t454 + t653 * t441 + t655)) * t613 + (qJD(3) * t591 - t766 * t464 + t742) * t441 / 0.2e1 + t464 * (t115 * t476 + t548) / 0.2e1 + ((t454 * t402 + t423 * t649 - t426 * t650) * t464 + (t423 * t560 - t426 * t561 + t454 * t588) * qJD(3)) * t623 + (t19 * t476 - t440 * t51 + t441 * t50 + t454 * t8 + t455 * t9) * t722 + (t10 * t454 + t11 * t455 + t20 * t476 - t440 * t49 + t441 * t48) * t724 + (t28 * t454 + t29 * t455 - t440 * t80 + t441 * t79 + t476 * t53) * t707 + (t21 * t454 + t22 * t455 + t35 * t476 - t440 * t75 + t441 * t74) * t711 + (t23 * t454 + t24 * t455 + t36 * t476 - t440 * t73 + t441 * t72) * t713 + (t12 * t454 + t13 * t455 + t26 * t476 - t440 * t55 + t441 * t54) * t715 + m(6) * (t15 * t364 + t40 * (t441 * t654 + t651) + t14 * t311 + t41 * (t440 * t654 + t195 + t662) + t7 * t307 + t37 * (-t440 * t657 - t441 * t656 + t655 - t663) + (t15 * (-t346 - t657) + t40 * (-t221 - t679) + t14 * t656 + t41 * t680) * t476 + (t14 * (-t407 - t654) + t41 * (-t395 - t658) + t7 * t657 + t37 * t679) * t455 + (t15 * t654 + t40 * t658 + t7 * (-t349 - t656) + t37 * (-t220 - t680)) * t454) + t641 * t728 + t330 * t729 + (t454 * t79 + t455 * t80 + t476 * t99) * t619 - t32 * t639 / 0.2e1 + t30 * t642 / 0.2e1 + ((-t134 * t661 + t137 * t323 + t140 * t324 - t166 * t749 + t168 * t280 - t170 * t283) * t258 + t49 * t330 + (t133 * t661 - t136 * t323 - t139 * t324 - t165 * t749 + t167 * t280 - t169 * t283) * t257 + t48 * t329 + (-t204 * t661 + t205 * t323 + t206 * t324 - t252 * t749 + t253 * t280 - t254 * t283) * t350 + t67 * t399) * t725 + ((-t134 * t660 + t137 * t325 + t140 * t326 - t166 * t576 + t168 * t284 + t170 * t285) * t258 + t51 * t330 + (t133 * t660 - t136 * t325 - t139 * t326 - t165 * t576 + t167 * t284 + t169 * t285) * t257 + t50 * t329 + (-t204 * t660 + t205 * t325 + t206 * t326 - t252 * t576 + t253 * t284 + t254 * t285) * t350 + t68 * t399) * t723 + ((t134 * t659 + t137 * t396 + t140 * t397 - t166 * t575 + t168 * t366 + t170 * t367) * t258 + t55 * t330 + (-t133 * t659 - t136 * t396 - t139 * t397 - t165 * t575 + t167 * t366 + t169 * t367) * t257 + t54 * t329 + (t204 * t659 + t205 * t396 + t206 * t397 - t252 * t575 + t253 * t366 + t254 * t367) * t350 + t70 * t399) * t716 - t399 * t18 / 0.2e1 - (t47 + t741) * t440 / 0.2e1; -t743 * t520 / 0.2e1 + t744 * t423 / 0.2e1 + t745 * t427 / 0.2e1 + (t423 * t72 + t427 * t73 + t520 * t783) * t720 + (t23 * t423 + t24 * t427 + t319 * t73 + t321 * t72 - t36 * t520 - t439 * t783) * t713 + (t34 * (-t218 * t520 + t290 * t423) + t83 * (t218 * t439 + t227 * t423 + t290 * t321 + t520 * t95) + t33 * (-t219 * t520 - t290 * t427) + t84 * (t219 * t439 - t227 * t427 - t290 * t319 - t520 * t94) + t27 * (-t218 * t427 - t219 * t423) + t76 * (-t218 * t319 - t219 * t321 - t423 * t94 + t427 * t95)) * t613 + ((t134 * t422 + t179 * t366 + t181 * t367) * t258 + (-t133 * t422 + t178 * t366 + t180 * t367) * t257 + (t204 * t422 + t260 * t366 + t261 * t367) * t350 + (-t377 * t54 + t379 * t55 + t422 * t70) * qJD(5) - t517 * t575) * t716 + ((-t134 * t377 + t179 * t280 - t181 * t283) * t258 + (t133 * t377 + t178 * t280 - t180 * t283) * t257 + (-t204 * t377 + t260 * t280 - t261 * t283) * t350 + (-t377 * t48 + t379 * t49 + t422 * t67) * qJD(5) - t517 * t749) * t725 + ((t134 * t379 + t179 * t284 + t181 * t285) * t258 + (-t133 * t379 + t178 * t284 + t180 * t285) * t257 + (t204 * t379 + t260 * t284 + t261 * t285) * t350 + (-t377 * t50 + t379 * t51 + t422 * t68) * qJD(5) - t517 * t576) * t723 + (t32 + t18) * t439 / 0.2e1 + (-t377 * t541 + t423 * t563 - t540 * t749) * t714 - (-t16 * t377 + t17 * t379 + t18 * t422) * qJD(5) / 0.2e1 + (t40 * (-t182 * t350 - t247 * t429 + t257 * t262 + t335 * t381 + (-t141 * t422 - t207 * t377) * qJD(5)) + t41 * (t183 * t350 + t250 * t429 - t258 * t262 - t335 * t382 + (t143 * t422 - t207 * t379) * qJD(5)) + t37 * (t182 * t258 - t183 * t257 + t247 * t382 - t250 * t381 + (t141 * t379 + t143 * t377) * qJD(5))) * t612 + (t83 * (-t245 * t429 + t334 * t381) + t84 * (t246 * t429 - t334 * t382) + t76 * (t245 * t382 - t246 * t381)) * t614 + (t379 * t541 + t427 * t563 - t540 * t576) * t712 + (t15 * (t423 * t654 + t520 * t657) + t40 * (t321 * t654 + t423 * t658 - t439 * t657 + t520 * t679) + t14 * (-t427 * t654 - t520 * t656) + t41 * (-t319 * t654 - t427 * t658 + t439 * t656 - t520 * t680) + t7 * (-t423 * t656 + t427 * t657) + t37 * (t319 * t657 - t321 * t656 - t423 * t680 + t427 * t679)) * t611 + (t423 * t79 + t427 * t80 - t520 * t99) * t619 + (t423 * t54 + t427 * t55 - t520 * t70) * t719 + (t423 * t74 + t427 * t75 - t520 * t87) * t721 + (-t19 * t520 + t319 * t51 + t321 * t50 + t423 * t8 + t427 * t9 + t439 * t68) * t722 + (t10 * t423 + t11 * t427 - t20 * t520 + t319 * t49 + t321 * t48 + t439 * t67) * t724 + (t423 * t48 + t427 * t49 - t520 * t67) * t726 + (t423 * t50 + t427 * t51 - t520 * t68) * t727 + (t28 * t423 + t29 * t427 + t319 * t80 + t321 * t79 + t439 * t99 - t520 * t53) * t707 + (t21 * t423 + t22 * t427 + t319 * t75 + t321 * t74 - t35 * t520 + t439 * t87) * t711 + (t12 * t423 + t13 * t427 - t26 * t520 + t319 * t55 + t321 * t54 + t439 * t70) * t715 + (t422 * t541 - t520 * t563 - t540 * t575) * t708 + t741 * t319 / 0.2e1 + t742 * t321 / 0.2e1; (t15 * (t141 * t575 - t207 * t749) + t40 * (-t107 * t749 - t141 * t362 + t202 * t207 + t575 * t63) + t14 * (-t143 * t575 + t207 * t576) + t41 * (t107 * t576 + t143 * t362 - t200 * t207 - t575 * t62) + t7 * (-t141 * t576 + t143 * t749) + t37 * (t141 * t200 - t143 * t202 - t576 * t63 + t62 * t749)) * t611 + t200 * t17 / 0.2e1 - t576 * t1 / 0.2e1 + (-t50 * t749 - t51 * t576 - t575 * t68) * t727 + (-t19 * t575 + t200 * t51 + t202 * t50 + t362 * t68 - t576 * t9 - t749 * t8) * t722 + t202 * t16 / 0.2e1 - t749 * t2 / 0.2e1 + (-t48 * t749 - t49 * t576 - t575 * t67) * t726 + (-t10 * t749 - t11 * t576 - t20 * t575 + t200 * t49 + t202 * t48 + t362 * t67) * t724 + t362 * t18 / 0.2e1 - t575 * t3 / 0.2e1 + (-t54 * t749 - t55 * t576 - t575 * t70) * t719 + (-t12 * t749 - t13 * t576 + t200 * t55 + t202 * t54 - t26 * t575 + t362 * t70) * t715 + (t40 * (-t163 * t350 + t238 * t257) + t41 * (t164 * t350 - t238 * t258) + t37 * (t163 * t258 - t164 * t257)) * t612 + (t284 * t542 + t285 * t543 - t564 * t576) * t723 + (t280 * t542 - t283 * t543 - t564 * t749) * t725 + (t366 * t542 + t367 * t543 - t564 * t575) * t716;];
tauc = t4(:);