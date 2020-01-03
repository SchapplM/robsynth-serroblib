% Calculate vector of inverse dynamics joint torques for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR7_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:54
% EndTime: 2019-12-31 18:19:50
% DurationCPUTime: 48.45s
% Computational Cost: add. (23762->936), mult. (22784->1223), div. (0->0), fcn. (20600->10), ass. (0->452)
t745 = Icges(4,3) + Icges(5,3);
t379 = qJ(3) + pkin(9);
t372 = sin(t379);
t374 = cos(t379);
t275 = Icges(5,5) * t374 - Icges(5,6) * t372;
t383 = sin(qJ(3));
t386 = cos(qJ(3));
t319 = Icges(4,5) * t386 - Icges(4,6) * t383;
t738 = t275 + t319;
t380 = qJ(1) + pkin(8);
t375 = cos(t380);
t744 = t745 * t375;
t373 = sin(t380);
t591 = t373 * t386;
t593 = t373 * t383;
t595 = t373 * t374;
t597 = t372 * t373;
t732 = -Icges(4,5) * t591 - Icges(5,5) * t595 + Icges(4,6) * t593 + Icges(5,6) * t597 + t744;
t735 = t745 * t373 + t738 * t375;
t612 = Icges(5,6) * t375;
t198 = Icges(5,4) * t595 - Icges(5,2) * t597 - t612;
t613 = Icges(4,6) * t375;
t210 = Icges(4,4) * t591 - Icges(4,2) * t593 - t613;
t743 = t198 * t372 + t210 * t383;
t622 = Icges(5,4) * t372;
t279 = Icges(5,1) * t374 - t622;
t201 = Icges(5,5) * t373 + t279 * t375;
t623 = Icges(4,4) * t383;
t323 = Icges(4,1) * t386 - t623;
t215 = Icges(4,5) * t373 + t323 * t375;
t742 = -t201 * t595 - t215 * t591;
t303 = Icges(5,4) * t597;
t617 = Icges(5,5) * t375;
t200 = Icges(5,1) * t595 - t303 - t617;
t331 = Icges(4,4) * t593;
t618 = Icges(4,5) * t375;
t214 = Icges(4,1) * t591 - t331 - t618;
t727 = -t200 * t374 - t214 * t386 + t743;
t276 = Icges(5,2) * t374 + t622;
t355 = Icges(5,4) * t374;
t278 = Icges(5,1) * t372 + t355;
t320 = Icges(4,2) * t386 + t623;
t376 = Icges(4,4) * t386;
t322 = Icges(4,1) * t383 + t376;
t737 = t276 * t372 - t278 * t374 + t320 * t383 - t322 * t386;
t741 = t735 * t375 + t742;
t584 = t375 * t386;
t590 = t374 * t375;
t740 = -t200 * t590 - t214 * t584 + t732 * t373;
t707 = t201 * t590 + t215 * t584 + t735 * t373;
t700 = -t727 * t373 + t732 * t375;
t456 = -Icges(5,2) * t372 + t355;
t199 = Icges(5,6) * t373 + t375 * t456;
t457 = -Icges(4,2) * t383 + t376;
t211 = Icges(4,6) * t373 + t375 * t457;
t699 = -t199 * t597 - t211 * t593 - t741;
t586 = t375 * t383;
t596 = t372 * t375;
t698 = -t198 * t596 - t210 * t586 - t740;
t697 = -t199 * t596 - t211 * t586 + t707;
t274 = Icges(5,5) * t372 + Icges(5,6) * t374;
t318 = Icges(4,5) * t383 + Icges(4,6) * t386;
t739 = t274 + t318;
t736 = -t276 * t374 - t278 * t372 - t320 * t386 - t322 * t383;
t598 = t318 * t375;
t600 = t274 * t375;
t712 = -t737 * t373 - t598 - t600;
t599 = t318 * t373;
t601 = t274 * t373;
t711 = -t737 * t375 + t599 + t601;
t381 = -qJ(4) - pkin(6);
t340 = t375 * t381;
t377 = t386 * pkin(3);
t368 = t377 + pkin(2);
t554 = -t373 * t368 - t340;
t384 = sin(qJ(1));
t655 = pkin(1) * t384;
t734 = t554 - t655;
t733 = t199 * t372 + t211 * t383;
t695 = t198 * t374 + t200 * t372 + t210 * t386 + t214 * t383;
t694 = t199 * t374 + t201 * t372 + t211 * t386 + t215 * t383;
t265 = t456 * qJD(3);
t266 = t279 * qJD(3);
t295 = t457 * qJD(3);
t296 = t323 * qJD(3);
t731 = t739 * qJD(1) + t736 * qJD(3) - t265 * t372 + t266 * t374 - t295 * t383 + t296 * t386;
t730 = t739 * qJD(3);
t527 = rSges(5,1) * t595;
t729 = -t527 + t734;
t728 = t201 * t374 + t215 * t386 - t733;
t726 = t697 * t373 - t698 * t375;
t725 = t699 * t373 - t700 * t375;
t724 = t737 * qJD(1) + t738 * qJD(3);
t387 = cos(qJ(1));
t378 = t387 * pkin(1);
t723 = t711 * qJD(1);
t283 = rSges(3,1) * t373 + rSges(3,2) * t375;
t261 = -t283 - t655;
t722 = t712 * qJD(1);
t721 = t735 * qJD(1);
t385 = cos(qJ(5));
t585 = t375 * t385;
t382 = sin(qJ(5));
t589 = t374 * t382;
t241 = t373 * t589 + t585;
t587 = t375 * t382;
t588 = t374 * t385;
t242 = t373 * t588 - t587;
t477 = t242 * rSges(6,1) - t241 * rSges(6,2);
t130 = -rSges(6,3) * t597 - t477;
t640 = rSges(6,2) * t382;
t643 = rSges(6,1) * t385;
t476 = -t640 + t643;
t216 = -rSges(6,3) * t374 + t372 * t476;
t537 = qJD(5) * t372;
t539 = qJD(3) * t375;
t256 = -t373 * t537 + t539;
t536 = qJD(5) * t374;
t335 = qJD(1) - t536;
t347 = qJD(4) * t373;
t653 = pkin(4) * t372;
t286 = -pkin(7) * t374 + t653;
t654 = pkin(3) * t383;
t499 = -t286 - t654;
t720 = t335 * t130 - t256 * t216 + t499 * t539 + t347;
t118 = Icges(6,5) * t242 - Icges(6,6) * t241 + Icges(6,3) * t597;
t227 = Icges(6,4) * t242;
t121 = -Icges(6,2) * t241 + Icges(6,6) * t597 + t227;
t226 = Icges(6,4) * t241;
t125 = -Icges(6,1) * t242 - Icges(6,5) * t597 + t226;
t701 = t121 * t382 + t125 * t385;
t48 = -t118 * t374 - t372 * t701;
t719 = t375 ^ 2;
t389 = qJD(1) ^ 2;
t528 = t389 * t378;
t718 = qJD(3) * t726 + t723;
t717 = qJD(3) * t725 + t722;
t423 = qJD(3) * t276;
t113 = qJD(1) * t199 - t373 * t423;
t425 = qJD(3) * t278;
t115 = qJD(1) * t201 - t373 * t425;
t424 = qJD(3) * t320;
t145 = qJD(1) * t211 - t373 * t424;
t426 = qJD(3) * t322;
t148 = qJD(1) * t215 - t373 * t426;
t716 = qJD(3) * t727 - t113 * t374 - t115 * t372 - t145 * t386 - t148 * t383;
t112 = -t375 * t423 + (-t373 * t456 + t612) * qJD(1);
t114 = -t375 * t425 + (-t279 * t373 + t617) * qJD(1);
t144 = -t375 * t424 + (-t373 * t457 + t613) * qJD(1);
t147 = -t375 * t426 + (-t323 * t373 + t618) * qJD(1);
t715 = qJD(3) * t728 + t112 * t374 + t114 * t372 + t144 * t386 + t147 * t383;
t714 = t373 * t724 + t375 * t731;
t713 = t373 * t731 - t375 * t724;
t710 = t732 * qJD(1) + qJD(3) * t695 + t113 * t372 - t115 * t374 + t145 * t383 - t148 * t386;
t709 = -qJD(3) * t694 - t112 * t372 + t114 * t374 - t144 * t383 + t147 * t386 + t721;
t708 = t732 + t733;
t706 = -t730 * t375 + (-t738 * t373 - t728 + t744) * qJD(1);
t705 = qJD(1) * t727 - t373 * t730 + t721;
t540 = qJD(3) * t373;
t704 = rSges(4,2) * t383;
t365 = t374 * pkin(4);
t683 = t372 * pkin(7) + t365;
t236 = t683 * t373;
t366 = t375 * pkin(6);
t287 = pkin(2) * t373 - t366;
t192 = t287 + t554;
t502 = -t287 - t655;
t488 = t192 + t502;
t438 = -t236 + t488;
t41 = qJD(1) * t438 + t720;
t592 = t373 * t385;
t243 = -t374 * t587 + t592;
t594 = t373 * t382;
t244 = t374 * t585 + t594;
t131 = t244 * rSges(6,1) + t243 * rSges(6,2) + rSges(6,3) * t596;
t238 = pkin(4) * t590 + pkin(7) * t596;
t255 = t375 * t537 + t540;
t364 = t373 * pkin(6);
t289 = t375 * pkin(2) + t364;
t493 = t375 * t368 - t373 * t381;
t193 = t493 - t289;
t501 = t289 + t378;
t487 = t193 + t501;
t531 = pkin(3) * t593;
t552 = qJD(3) * t531 + qJD(4) * t375;
t42 = -t286 * t540 + t131 * t335 - t216 * t255 + (t238 + t487) * qJD(1) - t552;
t633 = t373 * t42;
t702 = t375 * t41 + t633;
t37 = t118 * t597 - t121 * t241 - t125 * t242;
t120 = Icges(6,5) * t244 + Icges(6,6) * t243 + Icges(6,3) * t596;
t621 = Icges(6,4) * t244;
t123 = Icges(6,2) * t243 + Icges(6,6) * t596 + t621;
t228 = Icges(6,4) * t243;
t126 = Icges(6,1) * t244 + Icges(6,5) * t596 + t228;
t38 = t120 * t597 - t241 * t123 + t242 * t126;
t454 = Icges(6,5) * t385 - Icges(6,6) * t382;
t204 = -Icges(6,3) * t374 + t372 * t454;
t619 = Icges(6,4) * t385;
t455 = -Icges(6,2) * t382 + t619;
t208 = -Icges(6,6) * t374 + t372 * t455;
t620 = Icges(6,4) * t382;
t458 = Icges(6,1) * t385 - t620;
t212 = -Icges(6,5) * t374 + t372 * t458;
t66 = t204 * t597 - t208 * t241 + t212 * t242;
t12 = t255 * t38 - t256 * t37 + t335 * t66;
t39 = t118 * t596 + t243 * t121 - t125 * t244;
t40 = t120 * t596 + t243 * t123 + t244 * t126;
t67 = t204 * t596 + t208 * t243 + t212 * t244;
t13 = t255 * t40 - t256 * t39 + t67 * t335;
t690 = t710 * t719 + (t706 * t373 + (-t705 + t709) * t375) * t373;
t689 = t705 * t719 + (t709 * t373 + (-t706 + t710) * t375) * t373;
t273 = qJD(1) * t287;
t687 = qJD(1) * t192 - t273;
t358 = t373 * rSges(4,3);
t219 = rSges(4,1) * t584 - rSges(4,2) * t586 + t358;
t166 = t219 + t501;
t684 = t374 * rSges(5,1) - rSges(5,2) * t372;
t686 = t684 + t377;
t285 = t375 * rSges(3,1) - rSges(3,2) * t373;
t262 = t285 + t378;
t685 = -rSges(5,2) * t597 - t375 * rSges(5,3);
t674 = t373 * (-t276 * t375 + t201) - t375 * (-Icges(5,2) * t595 + t200 - t303);
t562 = -Icges(4,2) * t591 + t214 - t331;
t565 = t322 * t373 + t210;
t673 = -t383 * t562 - t386 * t565;
t205 = Icges(6,3) * t372 + t374 * t454;
t447 = -t208 * t382 + t212 * t385;
t452 = -t123 * t382 + t126 * t385;
t672 = t255 * (-t204 * t375 - t452) - t256 * (-t204 * t373 + t701) + t335 * (t205 - t447);
t249 = (-Icges(6,2) * t385 - t620) * t372;
t671 = t255 * (-Icges(6,2) * t244 + t126 + t228) - t256 * (-Icges(6,2) * t242 - t125 - t226) + t335 * (t212 + t249);
t670 = -m(5) - m(6);
t534 = qJD(1) * qJD(3);
t270 = qJDD(3) * t373 + t375 * t534;
t516 = t374 * t539;
t544 = qJD(1) * t373;
t519 = t372 * t544;
t417 = t516 - t519;
t532 = qJDD(5) * t372;
t135 = qJD(5) * t417 + t375 * t532 + t270;
t669 = t135 / 0.2e1;
t271 = -qJDD(3) * t375 + t373 * t534;
t542 = qJD(1) * t375;
t418 = t372 * t542 + t374 * t540;
t136 = qJD(5) * t418 + t373 * t532 + t271;
t668 = t136 / 0.2e1;
t667 = -t255 / 0.2e1;
t666 = t255 / 0.2e1;
t665 = -t256 / 0.2e1;
t664 = t256 / 0.2e1;
t263 = qJD(3) * t537 - qJDD(5) * t374 + qJDD(1);
t663 = t263 / 0.2e1;
t662 = t270 / 0.2e1;
t661 = t271 / 0.2e1;
t660 = -t335 / 0.2e1;
t659 = t335 / 0.2e1;
t658 = t373 / 0.2e1;
t657 = -t375 / 0.2e1;
t656 = -rSges(6,3) - pkin(7);
t652 = g(1) * t373;
t651 = g(2) * t373;
t541 = qJD(3) * t372;
t409 = t335 * t385 + t382 * t541;
t543 = qJD(1) * t374;
t492 = -qJD(5) + t543;
t104 = t373 * t409 - t492 * t587;
t408 = t335 * t382 - t385 * t541;
t105 = t373 * t408 + t492 * t585;
t56 = Icges(6,5) * t105 + Icges(6,6) * t104 + Icges(6,3) * t418;
t58 = Icges(6,4) * t105 + Icges(6,2) * t104 + Icges(6,6) * t418;
t60 = Icges(6,1) * t105 + Icges(6,4) * t104 + Icges(6,5) * t418;
t8 = (-qJD(3) * t701 - t56) * t374 + (qJD(3) * t118 - t382 * t58 + t385 * t60 + (-t121 * t385 + t125 * t382) * qJD(5)) * t372;
t650 = t8 * t256;
t102 = t375 * t409 + t492 * t594;
t103 = t375 * t408 - t492 * t592;
t55 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t417;
t57 = Icges(6,4) * t103 + Icges(6,2) * t102 + Icges(6,6) * t417;
t59 = Icges(6,1) * t103 + Icges(6,4) * t102 + Icges(6,5) * t417;
t9 = (qJD(3) * t452 - t55) * t374 + (qJD(3) * t120 - t382 * t57 + t385 * t59 + (-t123 * t385 - t126 * t382) * qJD(5)) * t372;
t649 = t9 * t255;
t646 = pkin(2) - t368;
t246 = (-Icges(6,5) * t382 - Icges(6,6) * t385) * t372;
t140 = qJD(3) * t205 + qJD(5) * t246;
t209 = Icges(6,6) * t372 + t374 * t455;
t143 = qJD(3) * t209 + qJD(5) * t249;
t213 = Icges(6,5) * t372 + t374 * t458;
t252 = (-Icges(6,1) * t382 - t619) * t372;
t146 = qJD(3) * t213 + qJD(5) * t252;
t23 = (qJD(3) * t447 - t140) * t374 + (qJD(3) * t204 - t143 * t382 + t146 * t385 + (-t208 * t385 - t212 * t382) * qJD(5)) * t372;
t78 = -t204 * t374 + t372 * t447;
t645 = t23 * t335 + t78 * t263;
t644 = rSges(4,1) * t386;
t152 = t418 * pkin(7) + (-t372 * t540 + t374 * t542) * pkin(4);
t257 = (-rSges(6,1) * t382 - rSges(6,2) * t385) * t372;
t356 = t372 * rSges(6,3);
t153 = qJD(5) * t257 + (t374 * t476 + t356) * qJD(3);
t268 = t683 * qJD(3);
t533 = qJD(1) * qJD(4);
t428 = qJDD(4) * t373 + t271 * t654 + t375 * t533 - t528;
t583 = t386 * qJD(3) ^ 2;
t529 = pkin(3) * t583;
t520 = t381 * t544 + t552;
t134 = (-t375 * t646 - t364) * qJD(1) - t520;
t269 = t289 * qJD(1);
t576 = -t134 - t269;
t478 = rSges(6,1) * t105 + rSges(6,2) * t104;
t62 = rSges(6,3) * t418 + t478;
t10 = t263 * t130 + t136 * t216 - t256 * t153 + t271 * t286 - t335 * t62 + (-qJD(3) * t268 - t529) * t375 + (-t152 + t576) * qJD(1) + t438 * qJDD(1) + t428;
t639 = t10 * t373;
t293 = pkin(7) * t516;
t517 = t372 * t539;
t419 = -t373 * t543 - t517;
t151 = pkin(4) * t419 - pkin(7) * t519 + t293;
t344 = pkin(6) * t542;
t513 = t383 * t539;
t437 = -pkin(3) * t513 + t347;
t133 = -t344 + (t373 * t646 - t340) * qJD(1) + t437;
t489 = qJDD(1) * t378 - t389 * t655;
t427 = qJD(1) * (-pkin(2) * t544 + t344) + qJDD(1) * t289 + t489;
t392 = qJD(1) * t133 + qJDD(1) * t193 + t373 * t533 + (-t270 * t383 - t373 * t583) * pkin(3) - qJDD(4) * t375 + t427;
t522 = t103 * rSges(6,1) + t102 * rSges(6,2) + rSges(6,3) * t516;
t61 = -rSges(6,3) * t519 + t522;
t11 = qJD(1) * t151 + qJDD(1) * t238 + t263 * t131 - t135 * t216 - t255 * t153 - t268 * t540 - t270 * t286 + t335 * t61 + t392;
t638 = t11 * t375;
t325 = rSges(4,1) * t383 + rSges(4,2) * t386;
t259 = t325 * t375;
t99 = qJD(1) * t166 - t325 * t540;
t635 = t259 * t99;
t357 = t373 * rSges(5,3);
t549 = rSges(4,2) * t593 + t375 * rSges(4,3);
t217 = rSges(4,1) * t591 - t549;
t486 = -t217 + t502;
t514 = t325 * t539;
t98 = qJD(1) * t486 - t514;
t632 = t373 * t98;
t631 = t375 * t98;
t630 = t48 * t136;
t49 = -t120 * t374 + t372 * t452;
t629 = t49 * t135;
t282 = rSges(5,1) * t372 + rSges(5,2) * t374;
t420 = -t282 - t654;
t416 = t420 * t539 + t347;
t202 = t527 + t685;
t439 = -t202 + t488;
t70 = qJD(1) * t439 + t416;
t628 = t70 * t282;
t578 = -t130 + t236;
t577 = t131 + t238;
t572 = -t373 * t192 + t375 * t193;
t203 = rSges(5,1) * t590 - rSges(5,2) * t596 + t357;
t569 = -t193 - t203;
t568 = -t193 - t238;
t564 = -t322 * t375 - t211;
t561 = -t320 * t375 + t215;
t559 = -t276 + t279;
t558 = t278 + t456;
t525 = t372 * t640;
t557 = rSges(6,3) * t595 + t373 * t525;
t556 = rSges(6,3) * t590 + t375 * t525;
t555 = rSges(5,2) * t519 + rSges(5,3) * t542;
t553 = rSges(4,3) * t542 + t544 * t704;
t551 = -t320 + t323;
t550 = t322 + t457;
t546 = qJD(1) * t275;
t545 = qJD(1) * t319;
t538 = qJD(3) * t386;
t530 = pkin(3) * t586;
t526 = t372 * t643;
t523 = pkin(3) * t538;
t521 = t375 * t133 + t373 * t134 - t192 * t542;
t512 = -t192 * t540 + t193 * t539 + qJD(2);
t511 = -pkin(2) - t644;
t508 = t542 / 0.2e1;
t507 = -t540 / 0.2e1;
t506 = t540 / 0.2e1;
t505 = -t539 / 0.2e1;
t504 = t539 / 0.2e1;
t491 = (-t373 ^ 2 - t719) * t654;
t490 = t133 * t539 + t134 * t540 - t270 * t192 + qJDD(2);
t484 = -t216 + t499;
t267 = t684 * qJD(3);
t483 = -t267 - t523;
t481 = t378 + t493;
t328 = rSges(2,1) * t387 - rSges(2,2) * t384;
t326 = rSges(2,1) * t384 + rSges(2,2) * t387;
t327 = t644 - t704;
t471 = t37 * t375 - t373 * t38;
t470 = t37 * t373 + t375 * t38;
t469 = t373 * t40 - t375 * t39;
t468 = t373 * t39 + t375 * t40;
t467 = t373 * t49 - t375 * t48;
t466 = t373 * t48 + t375 * t49;
t461 = -t373 * t99 - t631;
t451 = -t130 * t375 - t131 * t373;
t154 = -rSges(4,2) * t375 * t538 + (-t386 * t544 - t513) * rSges(4,1) + t553;
t258 = t325 * t373;
t155 = -qJD(3) * t258 + (t327 * t375 + t358) * qJD(1);
t450 = t154 * t375 + t155 * t373;
t444 = t217 * t373 + t219 * t375;
t218 = rSges(6,1) * t588 - rSges(6,2) * t589 + t356;
t436 = -t153 - t268 - t523;
t434 = t372 * t656 - t365;
t229 = t282 * t373;
t415 = -t118 * t256 + t120 * t255 + t204 * t335;
t414 = (-Icges(6,5) * t241 - Icges(6,6) * t242) * t256 - (Icges(6,5) * t243 - Icges(6,6) * t244) * t255 - t246 * t335;
t413 = t198 * t375 - t199 * t373;
t412 = -t383 * t561 + t386 * t564;
t411 = t375 * t420;
t410 = t372 * t414;
t404 = (-t372 * t558 + t374 * t559) * qJD(1);
t403 = (-t383 * t550 + t386 * t551) * qJD(1);
t401 = (Icges(6,1) * t243 - t123 - t621) * t255 - (-Icges(6,1) * t241 - t121 - t227) * t256 + (-t208 + t252) * t335;
t32 = -t130 * t255 + t131 * t256 + (t236 * t373 + t238 * t375) * qJD(3) + t512;
t399 = t32 * t451 + (t373 * t41 - t375 * t42) * t216;
t391 = -t372 * t674 + t413 * t374;
t390 = t672 * t372;
t311 = pkin(7) * t590;
t309 = pkin(7) * t595;
t297 = t327 * qJD(3);
t237 = -pkin(4) * t596 + t311;
t235 = -pkin(4) * t597 + t309;
t230 = t282 * t375;
t177 = -t375 * t526 + t556;
t176 = -t373 * t526 + t557;
t175 = t212 * t375;
t174 = t212 * t373;
t173 = t208 * t375;
t172 = t208 * t373;
t164 = rSges(6,1) * t243 - rSges(6,2) * t244;
t163 = -rSges(6,1) * t241 - rSges(6,2) * t242;
t117 = -qJD(3) * t229 + (t375 * t684 + t357) * qJD(1);
t116 = rSges(5,1) * t419 - rSges(5,2) * t516 + t555;
t97 = qJD(3) * t444 + qJD(2);
t71 = -t282 * t540 + (t203 + t487) * qJD(1) - t552;
t54 = (t202 * t373 + t203 * t375) * qJD(3) + t512;
t53 = qJD(1) * t154 + qJDD(1) * t219 - t270 * t325 - t297 * t540 + t427;
t52 = -t528 - t297 * t539 + t271 * t325 + (-t155 - t269) * qJD(1) + t486 * qJDD(1);
t47 = qJD(3) * t450 + t217 * t270 - t219 * t271 + qJDD(2);
t29 = qJD(1) * t116 + qJDD(1) * t203 - t267 * t540 - t270 * t282 + t392;
t28 = t271 * t282 + (-qJD(3) * t267 - t529) * t375 + (-t117 + t576) * qJD(1) + t439 * qJDD(1) + t428;
t17 = t202 * t270 + t569 * t271 + (t116 * t375 + t117 * t373) * qJD(3) + t490;
t16 = t104 * t208 + t105 * t212 + t140 * t597 - t143 * t241 + t146 * t242 + t204 * t418;
t15 = t102 * t208 + t103 * t212 + t140 * t596 + t143 * t243 + t146 * t244 + t204 * t417;
t14 = t255 * t49 - t256 * t48 + t335 * t78;
t7 = t104 * t123 + t105 * t126 + t120 * t418 - t241 * t57 + t242 * t59 + t55 * t597;
t6 = t104 * t121 - t105 * t125 + t118 * t418 - t241 * t58 + t242 * t60 + t56 * t597;
t5 = t102 * t123 + t103 * t126 + t120 * t417 + t243 * t57 + t244 * t59 + t55 * t596;
t4 = t102 * t121 - t103 * t125 + t118 * t417 + t243 * t58 + t244 * t60 + t56 * t596;
t3 = -t130 * t135 - t131 * t136 + t236 * t270 + t255 * t62 + t256 * t61 + t568 * t271 + (t151 * t375 + t152 * t373) * qJD(3) + t490;
t2 = t135 * t38 + t136 * t37 + t16 * t335 + t255 * t7 - t256 * t6 + t263 * t66;
t1 = t135 * t40 + t136 * t39 + t15 * t335 + t255 * t5 - t256 * t4 + t263 * t67;
t18 = [t66 * t668 + t67 * t669 + t16 * t665 + t15 * t666 - t650 / 0.2e1 + t649 / 0.2e1 + t645 + t629 / 0.2e1 + t630 / 0.2e1 - m(2) * (-g(1) * t326 + g(2) * t328) + ((t707 * t373 + ((t735 + t743) * t375 + t699 + t740 + t742) * t375) * qJD(3) + t723) * t504 + (t665 + t664) * t13 + (-t737 * qJD(3) + t265 * t374 + t266 * t372 + t295 * t386 + t296 * t383) * qJD(1) + ((-t283 * t389 - g(2) + t489) * t262 + (-t528 + (-0.2e1 * t285 - t378 + t262) * t389 - g(1)) * t261) * m(3) + (t41 * (-t478 + t520) + t42 * (-pkin(4) * t517 + t293 + t437 + t522) + (t10 * t434 + t41 * (t374 * t656 + t653) * qJD(3)) * t373 + ((-t384 * t42 - t387 * t41) * pkin(1) + (t41 * (-t368 - t683 - t356) - t42 * t381) * t375 + (-t368 + t434) * t633) * qJD(1) - t434 * t652 - (-t41 + (-t236 - t655) * qJD(1) + t687 + t720) * t42 + (t11 - g(2)) * (t481 + t577) + (t10 - g(1)) * (-t477 + t734)) * m(6) + (t628 * t540 + (t520 + (-t357 - t378 + (-t368 - t684) * t375) * qJD(1)) * t70 + (t411 * qJD(3) + t347 - t416 + t555 - t687 + t70 + (t202 + t655 + t729) * qJD(1)) * t71 + (t29 - g(2)) * (t203 + t481) + (t28 - g(1)) * (-t685 + t729)) * m(5) + (-(-t514 - t273 - t98 + (-t217 - t655) * qJD(1)) * t99 + t99 * (t344 + t553) + (t325 * t632 - t635) * qJD(3) + ((-t384 * t99 - t387 * t98) * pkin(1) + (-pkin(2) - t327) * t631 + (t98 * (-rSges(4,3) - pkin(6)) + t99 * t511) * t373) * qJD(1) + (t53 - g(2)) * t166 + (t52 - g(1)) * (t511 * t373 + t366 + t549 - t655)) * m(4) + (t711 + t694) * t662 + (t712 + t695) * t661 + (((t375 * t708 + t697 - t707) * t375 + (t373 * t708 + t698 + t741) * t373) * qJD(3) + t717 - t722) * t507 + (t714 + t715) * t506 + (m(3) * (t261 ^ 2 + t285 * t262) + m(2) * (t326 ^ 2 + t328 ^ 2) + Icges(2,3) + Icges(3,3) - t736) * qJDD(1) + (t713 - t716 + t718) * t505; m(3) * qJDD(2) + m(4) * t47 + m(5) * t17 + m(6) * t3 + (-m(3) - m(4) + t670) * g(3); ((t373 * t698 + t375 * t697) * qJD(1) + t690) * t506 + ((t373 * t700 + t375 * t699) * qJD(1) + t689) * t505 + t469 * t669 + (qJD(1) * t466 + t373 * t9 - t375 * t8) * t659 + t467 * t663 + (qJD(1) * t470 + t373 * t7 - t375 * t6) * t665 + (qJD(1) * t468 + t373 * t5 - t375 * t4) * t666 + (qJD(1) * t714 + t690 * qJD(3) + qJDD(1) * t711 + t697 * t270 + t698 * t271 + t1) * t658 + (t716 * t375 + t715 * t373 + (t695 * t373 + t694 * t375) * qJD(1)) * qJD(1) / 0.2e1 + (t12 + t717) * t544 / 0.2e1 + (t13 + t718) * t508 + (qJD(1) * t713 + t689 * qJD(3) + qJDD(1) * t712 + t699 * t270 + t700 * t271 + t2) * t657 + (-g(1) * (t311 - t530 + t556) - g(2) * (t309 - t531 + t557) - g(3) * t218 - (g(1) * t375 + t651) * t372 * (-pkin(4) - t643) - t41 * (-qJD(1) * t235 - t176 * t335 - t218 * t256) - t42 * (t177 * t335 - t218 * t255 + (t237 - t530) * qJD(1)) - ((t130 * t41 + t131 * t42) * t372 + t399 * t374) * qJD(5) + t3 * t572 + (t3 * t577 + t41 * t436 + (qJD(1) * t42 + t10) * t484) * t375 + (t11 * t484 + t42 * t436 + t3 * t578 + t41 * (t216 + t286) * qJD(1)) * t373 + (-qJD(3) * t702 + g(3)) * (-t683 - t377) + (-t176 * t255 - t177 * t256 - (t235 * t373 + t237 * t375 + t491) * qJD(3) + t521 + (qJD(1) * t578 + t151 + t61) * t375 + (t152 + t62 + (-t131 + t568) * qJD(1)) * t373) * t32) * m(6) - (t12 * t373 + t13 * t375) * t536 / 0.2e1 - (((t373 * t561 - t375 * t562) * t386 + (t373 * t564 + t375 * t565) * t383 + t413 * t372 + t374 * t674) * qJD(3) + (t372 * t559 + t374 * t558 + t383 * t551 + t386 * t550) * qJD(1)) * qJD(1) / 0.2e1 + ((-t540 * t600 + t546) * t373 + (t404 + (t373 * t601 + t391) * qJD(3)) * t375 + (-t540 * t598 + t545) * t373 + (t403 + (-t673 * t375 + (t599 + t412) * t373) * qJD(3)) * t375) * t507 + ((-t539 * t601 - t546) * t375 + (t404 + (t375 * t600 + t391) * qJD(3)) * t373 + (-t539 * t599 - t545) * t375 + (t403 + (t412 * t373 + (t598 - t673) * t375) * qJD(3)) * t373) * t504 + (-g(1) * t411 - g(3) * t686 - t420 * t651 + t17 * t572 + t54 * t521 + (t28 * t420 + t70 * t483 + t17 * t203 + t54 * t116 + (t54 * t202 + t420 * t71) * qJD(1)) * t375 + (t29 * t420 + t71 * t483 + t17 * t202 + t54 * t117 + (t54 * t569 + t628) * qJD(1)) * t373 - (t70 * t229 + t71 * (-t230 - t530)) * qJD(1) - (t54 * t491 + (-t54 * t230 - t686 * t70) * t375 + (-t54 * t229 - t686 * t71) * t373) * qJD(3)) * m(5) + (t373 * t694 - t375 * t695) * qJDD(1) / 0.2e1 - t136 * t471 / 0.2e1 + (t47 * t444 + t97 * ((t217 * t375 - t219 * t373) * qJD(1) + t450) + t461 * t297 + (-t53 * t373 - t52 * t375 + (-t375 * t99 + t632) * qJD(1)) * t325 + g(1) * t259 + g(2) * t258 - g(3) * t327 - (t258 * t98 - t635) * qJD(1) - (t97 * (-t258 * t373 - t259 * t375) + t461 * t327) * qJD(3)) * m(4) + ((t173 * t241 - t175 * t242) * t255 - (t172 * t241 - t174 * t242) * t256 + (-t209 * t241 + t213 * t242) * t335 + (t372 * t66 + t38 * t590) * qJD(5) + ((qJD(5) * t37 + t415) * t374 + t390) * t373) * t664 + ((-t173 * t243 - t175 * t244) * t255 - (-t172 * t243 - t174 * t244) * t256 + (t209 * t243 + t213 * t244) * t335 + (t372 * t67 + t39 * t595) * qJD(5) + ((qJD(5) * t40 + t415) * t374 + t390) * t375) * t667 + (((t173 * t382 - t175 * t385 + t120) * t255 - (t172 * t382 - t174 * t385 + t118) * t256 + (-t209 * t382 + t213 * t385 + t204) * t335 + t78 * qJD(5)) * t372 + (qJD(5) * t466 - t672) * t374) * t660 + t725 * t661 + t726 * t662 - t14 * t537 / 0.2e1; t670 * (-g(2) * t375 + t652) + 0.2e1 * (t639 / 0.2e1 - t638 / 0.2e1) * m(6) + 0.2e1 * (t28 * t658 + t29 * t657) * m(5); t1 * t596 / 0.2e1 + (t372 * t468 - t374 * t67) * t669 + ((qJD(3) * t468 - t15) * t374 + (-qJD(1) * t469 + qJD(3) * t67 + t373 * t4 + t375 * t5) * t372) * t666 + t2 * t597 / 0.2e1 + (t372 * t470 - t374 * t66) * t668 + ((qJD(3) * t470 - t16) * t374 + (qJD(1) * t471 + qJD(3) * t66 + t373 * t6 + t375 * t7) * t372) * t665 + t14 * t541 / 0.2e1 - t374 * (t629 + t630 + t645 + t649 - t650) / 0.2e1 + (t372 * t466 - t374 * t78) * t663 + ((qJD(3) * t466 - t23) * t374 + (-qJD(1) * t467 + qJD(3) * t78 + t373 * t8 + t375 * t9) * t372) * t659 + (t243 * t671 + t401 * t244 - t375 * t410) * t667 + (-t241 * t671 + t242 * t401 - t373 * t410) * t664 + (t414 * t374 + (-t382 * t671 + t385 * t401) * t372) * t660 + (-t519 / 0.2e1 + t374 * t504) * t13 + (t372 * t508 + t374 * t506) * t12 + ((qJD(3) * t399 - t10 * t130 - t11 * t131 + t41 * t62 - t42 * t61) * t374 + (t41 * (qJD(3) * t130 + t153 * t373) + t42 * (qJD(3) * t131 - t153 * t375) + t3 * t451 + t32 * (t130 * t544 - t131 * t542 - t373 * t61 + t375 * t62) + (qJD(1) * t702 - t638 + t639) * t216) * t372 - t41 * (-t163 * t335 - t256 * t257) - t42 * (t164 * t335 - t255 * t257) - t32 * (t163 * t255 + t164 * t256) - g(1) * t164 - g(2) * t163 - g(3) * t257) * m(6);];
tau = t18;
