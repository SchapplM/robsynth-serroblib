% Calculate vector of inverse dynamics joint torques for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR8_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:48
% EndTime: 2019-12-31 17:08:22
% DurationCPUTime: 28.88s
% Computational Cost: add. (7805->759), mult. (20649->959), div. (0->0), fcn. (19135->6), ass. (0->371)
t385 = cos(qJ(2));
t383 = sin(qJ(2));
t577 = Icges(3,4) * t383;
t302 = Icges(3,2) * t385 + t577;
t366 = Icges(4,5) * t383;
t720 = Icges(4,3) * t385 + t302 - t366;
t572 = Icges(4,5) * t385;
t304 = Icges(4,1) * t383 - t572;
t369 = Icges(3,4) * t385;
t718 = Icges(3,1) * t383 + t304 + t369;
t299 = Icges(3,5) * t385 - Icges(3,6) * t383;
t384 = sin(qJ(1));
t386 = cos(qJ(1));
t206 = Icges(3,3) * t384 + t299 * t386;
t301 = Icges(4,4) * t385 + Icges(4,6) * t383;
t208 = Icges(4,2) * t384 + t301 * t386;
t719 = t206 + t208;
t443 = Icges(4,1) * t385 + t366;
t211 = -Icges(4,4) * t386 + t384 * t443;
t551 = t383 * t384;
t340 = Icges(3,4) * t551;
t549 = t384 * t385;
t573 = Icges(3,5) * t386;
t213 = Icges(3,1) * t549 - t340 - t573;
t715 = -t211 - t213;
t212 = Icges(4,4) * t384 + t386 * t443;
t307 = Icges(3,1) * t385 - t577;
t214 = Icges(3,5) * t384 + t307 * t386;
t714 = t212 + t214;
t297 = Icges(4,3) * t383 + t572;
t203 = -Icges(4,6) * t386 + t297 * t384;
t566 = Icges(3,6) * t386;
t209 = Icges(3,4) * t549 - Icges(3,2) * t551 - t566;
t717 = t203 - t209;
t547 = t385 * t386;
t339 = Icges(4,5) * t547;
t550 = t383 * t386;
t565 = Icges(4,6) * t384;
t204 = Icges(4,3) * t550 + t339 + t565;
t441 = -Icges(3,2) * t383 + t369;
t210 = Icges(3,6) * t384 + t386 * t441;
t716 = -t204 + t210;
t713 = t297 - t441;
t712 = (Icges(3,6) - Icges(4,6)) * t385 + (Icges(4,4) + Icges(3,5)) * t383;
t709 = t307 + t443;
t708 = t720 * qJD(2);
t707 = t718 * qJD(2);
t706 = t204 * t550 + t719 * t384 + t714 * t547;
t207 = -Icges(4,2) * t386 + t301 * t384;
t190 = t384 * t207;
t563 = Icges(3,3) * t386;
t205 = Icges(3,5) * t549 - Icges(3,6) * t551 - t563;
t705 = -t203 * t550 - t384 * t205 + t715 * t547 - t190;
t696 = -t383 * t720 + t718 * t385;
t557 = t209 * t383;
t433 = -t213 * t385 + t557;
t558 = t207 * t386;
t436 = t203 * t383 + t211 * t385;
t637 = t384 * t436;
t71 = -t558 + t637;
t660 = -t205 * t386 - t384 * t433 + t71;
t659 = -t209 * t550 - t705;
t658 = -t210 * t550 + t706;
t461 = -t204 * t551 + t208 * t386 - t212 * t549;
t181 = t214 * t549;
t469 = t206 * t386 - t181;
t74 = -t210 * t551 - t469;
t657 = -t461 + t74;
t704 = t708 * t386 + (t384 * t441 - t203 - t566) * qJD(1);
t703 = t708 * t384 + (t297 * t386 - t210 + t565) * qJD(1);
t702 = -t707 * t386 + (-t307 * t384 - t211 + t573) * qJD(1);
t701 = -t714 * qJD(1) + t707 * t384;
t700 = t713 * qJD(2);
t699 = t709 * qJD(2);
t698 = -t299 - t301;
t697 = t712 * qJD(2);
t695 = -t383 * t718 - t385 * t720;
t556 = t210 * t383;
t694 = t204 * t383 + t714 * t385 - t556;
t693 = t433 - t436;
t648 = t714 * t383 + t716 * t385;
t647 = t715 * t383 + t717 * t385;
t642 = t712 * t386;
t641 = t712 * t384;
t653 = t696 * t384 - t642;
t652 = t696 * t386 + t641;
t692 = t719 * qJD(1);
t633 = qJD(2) - qJD(4);
t294 = t633 * t384;
t501 = qJD(2) * t386;
t295 = -qJD(4) * t386 + t501;
t382 = sin(qJ(4));
t601 = cos(qJ(4));
t486 = t383 * t601;
t282 = -t385 * t382 + t486;
t237 = t282 * t384;
t421 = t383 * t382 + t385 * t601;
t238 = t421 * t384;
t113 = Icges(5,5) * t238 + Icges(5,6) * t237 + Icges(5,3) * t386;
t574 = Icges(5,4) * t238;
t116 = Icges(5,2) * t237 + Icges(5,6) * t386 + t574;
t575 = Icges(5,4) * t237;
t120 = -Icges(5,1) * t238 - Icges(5,5) * t386 - t575;
t239 = t382 * t547 - t386 * t486;
t240 = t421 * t386;
t40 = -t113 * t384 - t239 * t116 - t120 * t240;
t115 = Icges(5,5) * t240 - Icges(5,6) * t239 - Icges(5,3) * t384;
t202 = Icges(5,4) * t240;
t118 = -Icges(5,2) * t239 - Icges(5,6) * t384 + t202;
t201 = Icges(5,4) * t239;
t121 = Icges(5,1) * t240 - Icges(5,5) * t384 - t201;
t463 = t115 * t384 + t239 * t118 - t240 * t121;
t165 = Icges(5,5) * t282 - Icges(5,6) * t421;
t269 = Icges(5,4) * t282;
t168 = -Icges(5,2) * t421 + t269;
t268 = Icges(5,4) * t421;
t171 = Icges(5,1) * t282 - t268;
t62 = -t165 * t384 - t168 * t239 + t171 * t240;
t14 = t62 * qJD(1) - t294 * t463 - t295 * t40;
t505 = qJD(1) * t386;
t163 = t633 * t421;
t654 = t633 * t282;
t79 = Icges(5,5) * t163 + Icges(5,6) * t654;
t80 = Icges(5,4) * t163 + Icges(5,2) * t654;
t81 = Icges(5,1) * t163 + Icges(5,4) * t654;
t506 = qJD(1) * t384;
t91 = t163 * t386 - t282 * t506;
t92 = -t386 * t654 - t421 * t506;
t15 = -t165 * t505 + t168 * t91 + t171 * t92 - t239 * t80 + t240 * t81 - t384 * t79;
t93 = t163 * t384 + t282 * t505;
t94 = qJD(1) * t240 - t384 * t654;
t16 = -t165 * t506 + t168 * t93 + t171 * t94 + t237 * t80 + t238 * t81 + t386 * t79;
t56 = Icges(5,4) * t94 + Icges(5,2) * t93 - Icges(5,6) * t506;
t58 = Icges(5,1) * t94 + Icges(5,4) * t93 - Icges(5,5) * t506;
t17 = t116 * t654 - t120 * t163 + t282 * t58 - t421 * t56;
t55 = Icges(5,4) * t92 + Icges(5,2) * t91 - Icges(5,6) * t505;
t57 = Icges(5,1) * t92 + Icges(5,4) * t91 - Icges(5,5) * t505;
t18 = t118 * t654 + t121 * t163 + t282 * t57 - t421 * t55;
t496 = qJD(1) * qJD(2);
t290 = qJDD(2) * t384 + t386 * t496;
t495 = qJD(1) * qJD(4);
t192 = -qJDD(4) * t384 - t386 * t495 + t290;
t350 = t384 * t496;
t193 = -t384 * t495 + t350 + (-qJDD(2) + qJDD(4)) * t386;
t38 = t113 * t386 + t116 * t237 - t120 * t238;
t39 = t386 * t115 + t237 * t118 + t238 * t121;
t401 = qJD(1) * (Icges(5,1) * t421 + t168 + t269) + t294 * (Icges(5,1) * t239 + t118 + t202) - t295 * (-Icges(5,1) * t237 + t116 + t574);
t410 = qJD(1) * (-Icges(5,5) * t421 - Icges(5,6) * t282) + (-Icges(5,5) * t237 + Icges(5,6) * t238) * t295 - (Icges(5,5) * t239 + Icges(5,6) * t240) * t294;
t595 = -qJD(1) / 0.2e1;
t54 = Icges(5,5) * t94 + Icges(5,6) * t93 - Icges(5,3) * t506;
t6 = -t113 * t505 + t116 * t91 - t120 * t92 - t239 * t56 + t240 * t58 - t384 * t54;
t606 = t295 / 0.2e1;
t61 = t165 * t386 + t168 * t237 + t171 * t238;
t618 = qJD(1) * (Icges(5,2) * t282 - t171 + t268) - t294 * (-Icges(5,2) * t240 + t121 - t201) + t295 * (-Icges(5,2) * t238 - t120 + t575);
t640 = t40 / 0.2e1 + t39 / 0.2e1;
t65 = -t116 * t421 - t120 * t282;
t66 = -t118 * t421 + t121 * t282;
t53 = Icges(5,5) * t92 + Icges(5,6) * t91 - Icges(5,3) * t505;
t7 = -t115 * t505 + t118 * t91 + t121 * t92 - t239 * t55 + t240 * t57 - t384 * t53;
t8 = -t113 * t506 + t116 * t93 - t120 * t94 + t237 * t56 + t238 * t58 + t386 * t54;
t9 = -t115 * t506 + t118 * t93 + t121 * t94 + t237 * t55 + t238 * t57 + t386 * t53;
t691 = (-t17 * t386 + t18 * t384 + (t384 * t65 + t386 * t66) * qJD(1)) * t595 - qJDD(1) * (t384 * t66 - t386 * t65) / 0.2e1 - t384 * (qJD(1) * t15 + qJDD(1) * t62 + t294 * t7 - t295 * t6) / 0.2e1 + t386 * (qJD(1) * t16 + qJDD(1) * t61 + t294 * t9 - t295 * t8) / 0.2e1 - (qJD(1) * t61 + t294 * t39 - t295 * t38) * t506 / 0.2e1 - t14 * t505 / 0.2e1 + (t38 * t386 - t384 * t640) * t193 + (t384 * t463 + t386 * t640) * t192 - (t239 * t618 - t240 * t401 + (-t463 * qJD(1) - t6) * t386 + (t40 * qJD(1) - t410 + t7) * t384) * t294 / 0.2e1 + (-t237 * t618 - t238 * t401 + (t38 * qJD(1) + t9) * t384 + (t39 * qJD(1) + t410 - t8) * t386) * t606;
t690 = t386 ^ 2;
t457 = -t238 * rSges(5,1) - t237 * rSges(5,2);
t123 = -rSges(5,3) * t386 + t457;
t598 = pkin(6) * t386;
t288 = pkin(3) * t549 + t598;
t689 = t123 - t288;
t688 = t712 * qJD(1) + t695 * qJD(2) + t700 * t383 + t699 * t385;
t687 = -t648 * qJD(2) + t383 * t704 + t702 * t385 + t692;
t634 = qJD(1) * t207;
t686 = -qJD(1) * t205 - qJD(2) * t647 - t383 * t703 + t385 * t701 - t634;
t685 = -t717 * t386 + (-Icges(4,1) * t550 + t304 * t386 + t339 - t716) * t384;
t684 = t658 * t384 - t659 * t386;
t683 = t657 * t384 - t660 * t386;
t682 = t718 - t713;
t681 = t709 - t720;
t680 = (Icges(3,2) * t549 + t340 + t715) * t386 + (-t302 * t386 + t714) * t384;
t679 = qJD(1) * t696 + qJD(2) * t698;
t678 = qJD(1) * t693 - t384 * t697 + t692;
t677 = -t634 - t697 * t386 + (-t299 * t384 + t563 - t694) * qJD(1);
t676 = -t401 * t282 + t421 * t618;
t673 = t652 * qJD(1);
t175 = rSges(5,1) * t421 + t282 * rSges(5,2);
t672 = t175 * t294;
t671 = t175 * t295;
t668 = t653 * qJD(1);
t159 = t237 * rSges(5,1) - t238 * rSges(5,2);
t161 = -t239 * rSges(5,1) - t240 * rSges(5,2);
t667 = t159 * t294 + t161 * t295;
t666 = t683 * qJD(2) + t668;
t665 = t684 * qJD(2) + t673;
t664 = qJD(2) * t693 + t383 * t701 + t385 * t703;
t663 = t694 * qJD(2) + t702 * t383 - t385 * t704;
t662 = -t679 * t384 + t688 * t386;
t661 = t688 * t384 + t679 * t386;
t651 = -t680 * t383 + t685 * t385;
t650 = (-t682 * t383 + t681 * t385) * qJD(1);
t649 = t558 + t706;
t646 = t678 * t690 + (t687 * t384 + (-t677 + t686) * t386) * t384;
t645 = t686 * t690 + (t677 * t384 + (-t678 + t687) * t386) * t384;
t494 = qJD(2) * qJD(3);
t644 = qJDD(3) * t383 + t385 * t494;
t643 = t698 * qJD(1);
t361 = t383 * qJ(3);
t635 = t385 * pkin(2) + t361;
t263 = t635 * t384;
t378 = t386 * pkin(5);
t316 = pkin(1) * t384 - t378;
t293 = qJD(1) * t316;
t636 = -qJD(1) * t263 - t293;
t347 = pkin(3) * t547;
t289 = -pkin(6) * t384 + t347;
t458 = t385 * rSges(4,1) + t383 * rSges(4,3);
t317 = t386 * pkin(1) + t384 * pkin(5);
t372 = t384 * rSges(4,2);
t232 = rSges(4,1) * t547 + rSges(4,3) * t550 + t372;
t348 = pkin(2) * t547;
t267 = qJ(3) * t550 + t348;
t465 = t267 + t317;
t129 = t465 + t232;
t596 = g(2) * t384;
t627 = (g(1) * t386 + t596) * t383;
t616 = m(4) / 0.2e1;
t615 = m(5) / 0.2e1;
t614 = -pkin(2) - pkin(3);
t611 = t290 / 0.2e1;
t291 = -qJDD(2) * t386 + t350;
t610 = t291 / 0.2e1;
t603 = -rSges(4,1) - pkin(2);
t602 = -rSges(5,3) - pkin(6);
t600 = pkin(3) * t385;
t593 = t92 * rSges(5,1) + t91 * rSges(5,2);
t592 = rSges(3,1) * t385;
t591 = rSges(4,1) * t383;
t371 = t384 * rSges(3,3);
t310 = rSges(3,1) * t383 + rSges(3,2) * t385;
t483 = t310 * t501;
t512 = rSges(3,2) * t551 + t386 * rSges(3,3);
t231 = rSges(3,1) * t549 - t512;
t529 = -t231 - t316;
t99 = qJD(1) * t529 - t483;
t588 = t384 * t99;
t587 = t386 * t99;
t585 = -rSges(4,3) - qJ(3);
t412 = -t383 * t501 - t385 * t506;
t194 = pkin(3) * t412 - pkin(6) * t505;
t59 = -rSges(5,3) * t505 + t593;
t584 = t194 + t59;
t503 = qJD(2) * t384;
t482 = t383 * t503;
t326 = pkin(3) * t482;
t195 = qJD(1) * t289 - t326;
t462 = rSges(5,1) * t94 + rSges(5,2) * t93;
t60 = -rSges(5,3) * t506 + t462;
t583 = t195 + t60;
t233 = rSges(3,1) * t547 - rSges(3,2) * t550 + t371;
t177 = t233 + t317;
t100 = qJD(1) * t177 - t310 * t503;
t266 = t310 * t386;
t560 = t100 * t266;
t546 = t385 * qJD(2) ^ 2;
t531 = t240 * rSges(5,1) - t239 * rSges(5,2);
t124 = -rSges(5,3) * t384 + t531;
t540 = t124 + t289;
t528 = -t232 - t267;
t527 = t384 * t263 + t386 * t267;
t500 = qJD(3) * t385;
t241 = qJD(2) * t635 - t500;
t526 = -t458 * qJD(2) - t241;
t336 = qJ(3) * t547;
t264 = -pkin(2) * t550 + t336;
t525 = qJD(1) * t264 + t384 * t500;
t524 = -t263 - t316;
t523 = -t267 - t289;
t357 = pkin(5) * t505;
t522 = qJD(1) * (-pkin(1) * t506 + t357) + qJDD(1) * t317;
t308 = pkin(2) * t383 - qJ(3) * t385;
t309 = -rSges(4,3) * t385 + t591;
t517 = -t308 - t309;
t516 = -t635 - t458;
t499 = qJD(3) * t386;
t332 = t383 * t499;
t481 = t385 * t501;
t515 = qJ(3) * t481 + t332;
t514 = rSges(4,2) * t505 + rSges(4,3) * t481;
t484 = t383 * t506;
t513 = rSges(3,2) * t484 + rSges(3,3) * t505;
t511 = t384 ^ 2 + t690;
t504 = qJD(2) * t383;
t502 = qJD(2) * t385;
t360 = qJD(3) * t383;
t126 = pkin(2) * t412 - qJ(3) * t484 + t515;
t246 = t383 * t505 + t384 * t502;
t327 = pkin(2) * t482;
t479 = t384 * t360;
t460 = t327 - t479;
t127 = qJ(3) * t246 + qJD(1) * t348 - t460;
t492 = t386 * t126 + t384 * t127 + t263 * t505;
t491 = -t124 + t523;
t490 = t291 * t308 + t386 * t644;
t334 = qJ(3) * t549;
t260 = -pkin(2) * t551 + t334;
t489 = t260 * t503 + t264 * t501 + t360;
t375 = t386 * rSges(4,2);
t230 = t384 * t458 - t375;
t488 = -t230 + t524;
t487 = t357 + t515;
t485 = t386 * t603;
t478 = -pkin(1) - t592;
t474 = -t503 / 0.2e1;
t473 = t503 / 0.2e1;
t472 = -t501 / 0.2e1;
t471 = t501 / 0.2e1;
t470 = -pkin(3) * t383 - t308;
t468 = -t205 + t556;
t467 = qJD(2) * t526;
t466 = t524 + t689;
t174 = rSges(5,1) * t282 - rSges(5,2) * t421;
t464 = -t174 + t470;
t459 = t263 * t503 + t267 * t501 - t500;
t315 = rSges(2,1) * t386 - rSges(2,2) * t384;
t311 = rSges(2,1) * t384 + rSges(2,2) * t386;
t314 = -rSges(3,2) * t383 + t592;
t287 = t317 * qJD(1);
t423 = -t127 - t287 - t479;
t82 = rSges(5,1) * t163 + rSges(5,2) * t654;
t11 = -t241 * t501 + t174 * t193 - t295 * t82 + (t291 * t383 - t386 * t546) * pkin(3) + t466 * qJDD(1) + (t423 - t583) * qJD(1) + t490;
t413 = qJDD(1) * t267 + t522 + t644 * t384 + (t126 + t332) * qJD(1);
t12 = -t241 * t503 - t174 * t192 - t290 * t308 - t294 * t82 + t540 * qJDD(1) + t584 * qJD(1) + (-t290 * t383 - t384 * t546) * pkin(3) + t413;
t455 = t11 * t386 + t12 * t384;
t403 = -t295 * t174 + t470 * t501 + t332;
t42 = qJD(1) * t466 + t403;
t272 = t308 * t503;
t43 = t479 - t174 * t294 - t272 - t326 + (t317 - t491) * qJD(1);
t450 = -t384 * t43 - t386 * t42;
t445 = -t100 * t384 - t587;
t155 = rSges(3,1) * t412 - rSges(3,2) * t481 + t513;
t262 = t310 * t384;
t157 = -qJD(2) * t262 + (t314 * t386 + t371) * qJD(1);
t438 = t155 * t386 + t157 * t384;
t431 = t231 * t384 + t233 * t386;
t426 = -pkin(3) * t502 - t241 - t82;
t425 = -pkin(1) - t635;
t422 = t501 * t517 + t332;
t414 = -qJDD(3) * t385 + t126 * t501 + t127 * t503 + t290 * t263 + t383 * t494;
t409 = t425 - t600;
t406 = t383 * t585 + t385 * t603 - pkin(1);
t344 = rSges(4,3) * t547;
t342 = rSges(4,3) * t549;
t333 = t385 * t499;
t284 = t314 * qJD(2);
t273 = t308 * t506;
t265 = -rSges(4,1) * t550 + t344;
t261 = -rSges(4,1) * t551 + t342;
t247 = t481 - t484;
t245 = t511 * t504;
t156 = -t309 * t503 + (t386 * t458 + t372) * qJD(1);
t154 = rSges(4,1) * t412 - rSges(4,3) * t484 + t514;
t95 = t431 * qJD(2);
t69 = -t272 + (-qJD(2) * t309 + t360) * t384 + t129 * qJD(1);
t68 = qJD(1) * t488 + t422;
t67 = (t230 * t384 + t232 * t386) * qJD(2) + t459;
t64 = qJD(1) * t155 + qJDD(1) * t233 - t284 * t503 - t290 * t310 + t522;
t63 = -t284 * t501 + t291 * t310 + t529 * qJDD(1) + (-t157 - t287) * qJD(1);
t37 = -t123 * t294 + t124 * t295 + (t288 * t384 + t289 * t386) * qJD(2) + t459;
t31 = qJD(1) * t154 + qJDD(1) * t232 + t290 * t517 + t384 * t467 + t413;
t30 = t291 * t309 + t386 * t467 + t488 * qJDD(1) + (-t156 + t423) * qJD(1) + t490;
t21 = t230 * t290 + t528 * t291 + (t154 * t386 + t156 * t384) * qJD(2) + t414;
t10 = -t123 * t192 - t124 * t193 + t288 * t290 + t294 * t60 + t295 * t59 + t523 * t291 + (t194 * t386 + t195 * t384) * qJD(2) + t414;
t1 = [t14 * t606 - m(2) * (-g(1) * t311 + g(2) * t315) + (t66 + t62) * t192 / 0.2e1 + (t65 + t61) * t193 / 0.2e1 + (t15 + t18) * t294 / 0.2e1 + (((t74 - t181 + (t206 + t557) * t386 + t705) * t386 + (t71 - t637 + t649) * t384) * qJD(2) + t673) * t471 - (t16 + t14 + t17) * t295 / 0.2e1 + (qJD(2) * t696 + t163 * t171 + t654 * t168 + t282 * t81 + t383 * t699 - t385 * t700 - t421 * t80) * qJD(1) + (-(t689 * qJD(1) + t403 - t42 + t636) * t43 + t11 * (t378 + t457) + t42 * (t326 + t327 - t462) + t43 * (t487 + t593) + (t43 * t504 * t614 + t11 * t602) * t386 + (t11 * t409 + t42 * (-qJ(3) * t502 - t360)) * t384 + ((t409 * t42 + t43 * t602) * t386 + (t42 * (-pkin(5) - t602) + t43 * t409) * t384) * qJD(1) - (t123 + t378 - t598 + (t385 * t614 - pkin(1) - t361) * t384) * g(1) + (t12 - g(2)) * (t384 * t602 + t347 + t465 + t531)) * m(5) + (-(-qJD(1) * t230 + t422 + t636 - t68) * t69 + t68 * t460 + t69 * (t487 + t514) + (t69 * t383 * t485 + t68 * (t385 * t585 + t591) * t384) * qJD(2) + (t68 * t406 * t386 + (t68 * (-rSges(4,2) - pkin(5)) + t69 * (t425 - t458)) * t384) * qJD(1) + (t31 - g(2)) * t129 + (t30 - g(1)) * (t384 * t406 + t375 + t378)) * m(4) + (t100 * (t357 + t513) + (t310 * t588 - t560) * qJD(2) + ((-pkin(1) - t314) * t587 + (t99 * (-rSges(3,3) - pkin(5)) + t100 * t478) * t384) * qJD(1) - (-qJD(1) * t231 - t293 - t483 - t99) * t100 + (t64 - g(2)) * t177 + (t63 - g(1)) * (t478 * t384 + t378 + t512)) * m(3) + (t648 + t652) * t611 + (-t647 + t653) * t610 + (t662 + t663) * t473 + (((t386 * t468 - t649 + t658) * t386 + (t384 * t468 - t190 + t461 + t469 + t659) * t384) * qJD(2) + t666 - t668) * t474 + (-t168 * t421 + t171 * t282 + m(2) * (t311 ^ 2 + t315 ^ 2) + Icges(2,3) - t695) * qJDD(1) + (t661 - t664 + t665) * t472; t684 * t611 + t683 * t610 + (t662 * qJD(1) + t645 * qJD(2) + t652 * qJDD(1) + t658 * t290 + t659 * t291) * t384 / 0.2e1 - (t661 * qJD(1) + t646 * qJD(2) + t653 * qJDD(1) + t657 * t290 + t660 * t291) * t386 / 0.2e1 + (t664 * t386 + t663 * t384 + (-t647 * t384 + t648 * t386) * qJD(1)) * qJD(1) / 0.2e1 + (t648 * t384 + t647 * t386) * qJDD(1) / 0.2e1 + t666 * t506 / 0.2e1 + t665 * t505 / 0.2e1 + ((-t503 * t642 - t643) * t384 + ((t384 * t641 + t651) * qJD(2) + t650) * t386) * t474 + ((t659 * t384 + t386 * t658) * qJD(1) + t645) * t473 + ((t384 * t660 + t657 * t386) * qJD(1) + t646) * t472 + ((-t501 * t641 + t643) * t386 + ((t386 * t642 + t651) * qJD(2) + t650) * t384) * t471 + ((t685 * t383 + t680 * t385) * qJD(2) + (t681 * t383 + t682 * t385) * qJD(1) - t676) * t595 + (t42 * t273 + t10 * t527 + t37 * t492 + (t12 * t464 + t43 * t426 - t10 * t689 + t37 * t583 + (t42 * t174 + t37 * t491) * qJD(1)) * t384 + (t11 * t464 + t42 * t426 + t10 * t540 + t37 * t584 + (-t37 * t689 + t43 * t464) * qJD(1)) * t386 - g(1) * (t336 - t161) - g(2) * (t334 - t159) - g(3) * (t175 + t635 + t600) - t614 * t627 - t42 * (t333 - t671) - t43 * (t525 - t672) - t37 * (t489 - t667) - (t42 * (t159 - t260) + t43 * (-pkin(3) * t550 - t161)) * qJD(1) - (t450 * t635 + (-t37 * t383 * t511 + t385 * t450) * pkin(3)) * qJD(2)) * m(5) + (-t68 * (t333 + (-t260 - t261) * qJD(1)) - t69 * (qJD(1) * t265 + t525) - t67 * t489 - ((t67 * t265 + t516 * t68) * t386 + (t67 * t261 + t516 * t69) * t384) * qJD(2) - g(1) * (t336 + t344) - g(2) * (t334 + t342) + g(3) * t516 - (g(1) * t485 + t596 * t603) * t383 + t68 * t273 + t21 * t527 + t67 * t492 + (t30 * t517 + t68 * t526 + t21 * t232 + t67 * t154 + (t67 * t230 + t517 * t69) * qJD(1)) * t386 + (t31 * t517 + t69 * t526 + t21 * t230 + t67 * t156 + (t68 * t309 + t528 * t67) * qJD(1)) * t384) * m(4) + (g(1) * t266 + g(2) * t262 - g(3) * t314 - (t262 * t99 - t560) * qJD(1) - (t95 * (-t262 * t384 - t266 * t386) + t445 * t314) * qJD(2) + (qJD(2) * t438 + t231 * t290 - t233 * t291) * t431 + t95 * ((t231 * t386 - t233 * t384) * qJD(1) + t438) + t445 * t284 + (-t64 * t384 - t63 * t386 + (-t100 * t386 + t588) * qJD(1)) * t310) * m(3) - t691; (-m(4) - m(5)) * (-g(3) * t385 + t627) - m(4) * (t245 * t67 + t246 * t69 + t247 * t68) - m(5) * (t245 * t37 + t246 * t43 + t247 * t42) + 0.2e1 * ((t501 * t68 + t503 * t69 - t21) * t616 + (t42 * t501 + t43 * t503 - t10) * t615) * t385 + 0.2e1 * ((qJD(2) * t67 + t30 * t386 + t31 * t384 + t505 * t69 - t506 * t68) * t616 + (qJD(2) * t37 - t42 * t506 + t43 * t505 + t455) * t615) * t383; t676 * t595 + ((t42 * t82 - t10 * t124 + t37 * (qJD(1) * t123 - t59)) * t386 + (t43 * t82 + t10 * t123 + t37 * (qJD(1) * t124 - t60)) * t384 + ((-t384 * t42 + t386 * t43) * qJD(1) + t455) * t174 - t42 * (-qJD(1) * t159 + t671) - t43 * (qJD(1) * t161 + t672) - t37 * t667 - g(1) * t161 - g(2) * t159 + g(3) * t175) * m(5) + t691;];
tau = t1;
