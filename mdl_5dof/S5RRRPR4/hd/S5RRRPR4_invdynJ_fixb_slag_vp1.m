% Calculate vector of inverse dynamics joint torques for
% S5RRRPR4
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:00
% EndTime: 2019-12-31 21:11:38
% DurationCPUTime: 32.36s
% Computational Cost: add. (21597->869), mult. (25940->1040), div. (0->0), fcn. (23785->8), ass. (0->439)
t471 = cos(qJ(3));
t469 = sin(qJ(3));
t684 = Icges(4,4) * t469;
t391 = Icges(4,2) * t471 + t684;
t459 = Icges(5,5) * t469;
t552 = Icges(5,3) * t471 - t459;
t832 = t391 + t552;
t681 = Icges(5,5) * t471;
t393 = Icges(5,1) * t469 - t681;
t460 = Icges(4,4) * t471;
t395 = Icges(4,1) * t469 + t460;
t828 = t393 + t395;
t467 = qJ(1) + qJ(2);
t456 = sin(t467);
t388 = Icges(4,5) * t471 - Icges(4,6) * t469;
t457 = cos(t467);
t520 = t388 * t457;
t260 = Icges(4,3) * t456 + t520;
t390 = Icges(5,4) * t471 + Icges(5,6) * t469;
t521 = t390 * t457;
t262 = Icges(5,2) * t456 + t521;
t831 = t260 + t262;
t555 = Icges(5,1) * t471 + t459;
t265 = -Icges(5,4) * t457 + t456 * t555;
t662 = t456 * t469;
t422 = Icges(4,4) * t662;
t661 = t456 * t471;
t267 = Icges(4,1) * t661 - Icges(4,5) * t457 - t422;
t830 = -t265 - t267;
t523 = t555 * t457;
t266 = Icges(5,4) * t456 + t523;
t396 = Icges(4,1) * t471 - t684;
t524 = t396 * t457;
t268 = Icges(4,5) * t456 + t524;
t827 = t266 + t268;
t386 = Icges(5,3) * t469 + t681;
t553 = -Icges(4,2) * t469 + t460;
t826 = t386 - t553;
t387 = Icges(4,5) * t469 + Icges(4,6) * t471;
t389 = Icges(5,4) * t469 - Icges(5,6) * t471;
t825 = t387 + t389;
t824 = t396 + t555;
t466 = qJD(1) + qJD(2);
t823 = (-Icges(4,6) + Icges(5,6)) * t466 + t832 * qJD(3);
t822 = (Icges(5,4) + Icges(4,5)) * t466 - t828 * qJD(3);
t656 = t457 * t471;
t421 = Icges(5,5) * t656;
t657 = t457 * t469;
t258 = Icges(5,6) * t456 + Icges(5,3) * t657 + t421;
t821 = t258 * t657 + t831 * t456 + t827 * t656;
t811 = -t469 * t832 + t828 * t471;
t261 = -Icges(5,2) * t457 + t390 * t456;
t237 = t456 * t261;
t257 = -Icges(5,6) * t457 + t386 * t456;
t259 = Icges(4,5) * t661 - Icges(4,6) * t662 - Icges(4,3) * t457;
t820 = -t257 * t657 - t456 * t259 + t830 * t656 - t237;
t659 = t457 * t261;
t544 = t257 * t469 + t265 * t471;
t748 = t456 * t544;
t107 = -t659 + t748;
t263 = Icges(4,4) * t661 - Icges(4,2) * t662 - Icges(4,6) * t457;
t672 = t263 * t469;
t541 = -t267 * t471 + t672;
t766 = -t457 * t259 - t456 * t541 + t107;
t765 = -t263 * t657 - t820;
t522 = t553 * t457;
t264 = Icges(4,6) * t456 + t522;
t764 = -t264 * t657 + t821;
t663 = t456 * t466;
t819 = t823 * t457 - t826 * t663;
t658 = t457 * t466;
t818 = t386 * t658 + t823 * t456 - t466 * t522;
t817 = t822 * t457 - t824 * t663;
t816 = (-t523 - t524) * t466 - t822 * t456;
t815 = t826 * qJD(3);
t814 = t824 * qJD(3);
t813 = t388 + t390;
t812 = (-Icges(5,2) - Icges(4,3)) * t466 + t825 * qJD(3);
t810 = -t828 * t469 - t471 * t832;
t671 = t264 * t469;
t809 = t258 * t469 + t827 * t471 - t671;
t808 = t541 - t544;
t807 = (t258 - t264) * t471 - t827 * t469;
t758 = (t257 - t263) * t471 + t830 * t469;
t224 = t268 * t661;
t573 = t457 * t260 - t224;
t110 = -t264 * t662 - t573;
t564 = -t258 * t662 + t262 * t457 - t266 * t661;
t757 = -t564 + t110;
t665 = t389 * t457;
t668 = t387 * t457;
t763 = t456 * t811 - t665 - t668;
t666 = t389 * t456;
t669 = t387 * t456;
t762 = t457 * t811 + t666 + t669;
t468 = sin(qJ(5));
t704 = cos(qJ(5));
t528 = t469 * t468 + t471 * t704;
t741 = qJD(3) - qJD(5);
t240 = t741 * t528;
t364 = -t471 * t468 + t469 * t704;
t128 = t240 * t457 - t364 * t663;
t767 = t741 * t364;
t129 = -t457 * t767 - t528 * t663;
t304 = t364 * t457;
t305 = t528 * t457;
t184 = Icges(6,5) * t305 + Icges(6,6) * t304 - Icges(6,3) * t456;
t289 = Icges(6,4) * t305;
t187 = Icges(6,2) * t304 - Icges(6,6) * t456 + t289;
t288 = Icges(6,4) * t304;
t190 = Icges(6,1) * t305 - Icges(6,5) * t456 + t288;
t78 = Icges(6,5) * t129 + Icges(6,6) * t128 - Icges(6,3) * t658;
t80 = Icges(6,4) * t129 + Icges(6,2) * t128 - Icges(6,6) * t658;
t82 = Icges(6,1) * t129 + Icges(6,4) * t128 - Icges(6,5) * t658;
t10 = t128 * t187 + t129 * t190 - t184 * t658 + t304 * t80 + t305 * t82 - t456 * t78;
t130 = t240 * t456 + t364 * t658;
t131 = t305 * t466 - t456 * t767;
t302 = t364 * t456;
t303 = t528 * t456;
t182 = Icges(6,5) * t303 + Icges(6,6) * t302 + Icges(6,3) * t457;
t682 = Icges(6,4) * t303;
t185 = Icges(6,2) * t302 + Icges(6,6) * t457 + t682;
t683 = Icges(6,4) * t302;
t189 = -Icges(6,1) * t303 - Icges(6,5) * t457 - t683;
t79 = Icges(6,5) * t131 + Icges(6,6) * t130 - Icges(6,3) * t663;
t81 = Icges(6,4) * t131 + Icges(6,2) * t130 - Icges(6,6) * t663;
t83 = Icges(6,1) * t131 + Icges(6,4) * t130 - Icges(6,5) * t663;
t11 = t130 * t185 - t131 * t189 - t182 * t663 + t302 * t81 + t303 * t83 + t457 * t79;
t12 = t130 * t187 + t131 * t190 - t184 * t663 + t302 * t80 + t303 * t82 + t457 * t78;
t610 = qJD(3) * t466;
t334 = qJDD(3) * t456 + t457 * t610;
t606 = qJD(5) * t466;
t220 = -qJDD(5) * t456 - t457 * t606 + t334;
t383 = t456 * t610;
t221 = -t456 * t606 + t383 + (-qJDD(3) + qJDD(5)) * t457;
t350 = t741 * t456;
t611 = qJD(3) * t457;
t351 = -qJD(5) * t457 + t611;
t565 = t184 * t456 - t304 * t187 - t305 * t190;
t60 = -t182 * t456 + t304 * t185 - t189 * t305;
t242 = Icges(6,5) * t364 - Icges(6,6) * t528;
t347 = Icges(6,4) * t364;
t245 = -Icges(6,2) * t528 + t347;
t346 = Icges(6,4) * t528;
t248 = Icges(6,1) * t364 - t346;
t96 = -t242 * t456 + t245 * t304 + t248 * t305;
t23 = -t350 * t565 - t351 * t60 + t96 * t466;
t26 = t185 * t767 - t189 * t240 + t364 * t83 - t528 * t81;
t27 = t187 * t767 + t190 * t240 + t364 * t82 - t528 * t80;
t140 = Icges(6,5) * t240 + Icges(6,6) * t767;
t141 = Icges(6,4) * t240 + Icges(6,2) * t767;
t142 = Icges(6,1) * t240 + Icges(6,4) * t767;
t30 = t128 * t245 + t129 * t248 - t140 * t456 + t141 * t304 + t142 * t305 - t242 * t658;
t31 = t130 * t245 + t131 * t248 + t140 * t457 + t141 * t302 + t142 * t303 - t242 * t663;
t465 = qJDD(1) + qJDD(2);
t490 = t350 * (-Icges(6,1) * t304 + t187 + t289) - t351 * (-Icges(6,1) * t302 + t185 + t682) + t466 * (Icges(6,1) * t528 + t245 + t347);
t514 = (-Icges(6,5) * t302 + Icges(6,6) * t303) * t351 - (-Icges(6,5) * t304 + Icges(6,6) * t305) * t350 + (-Icges(6,5) * t528 - Icges(6,6) * t364) * t466;
t58 = t182 * t457 + t185 * t302 - t189 * t303;
t59 = t457 * t184 + t302 * t187 + t303 * t190;
t708 = -t466 / 0.2e1;
t712 = t351 / 0.2e1;
t750 = t350 * (-Icges(6,2) * t305 + t190 + t288) - t351 * (-Icges(6,2) * t303 - t189 + t683) - t466 * (Icges(6,2) * t364 - t248 + t346);
t753 = t60 / 0.2e1 + t59 / 0.2e1;
t9 = t128 * t185 - t129 * t189 - t182 * t658 + t304 * t81 + t305 * t83 - t456 * t79;
t95 = t242 * t457 + t245 * t302 + t248 * t303;
t97 = -t185 * t528 - t189 * t364;
t98 = -t187 * t528 + t190 * t364;
t806 = -t456 * (t10 * t350 + t30 * t466 - t351 * t9 + t465 * t96) / 0.2e1 + t457 * (-t11 * t351 + t12 * t350 + t31 * t466 + t465 * t95) / 0.2e1 - t465 * (t456 * t98 - t457 * t97) / 0.2e1 + ((t466 * t98 - t26) * t457 + (t466 * t97 + t27) * t456) * t708 - (t350 * t59 - t351 * t58 + t466 * t95) * t663 / 0.2e1 - t23 * t658 / 0.2e1 + (-t456 * t753 + t457 * t58) * t221 + (t456 * t565 + t457 * t753) * t220 - (t304 * t750 - t305 * t490 + (-t466 * t565 - t9) * t457 + (t466 * t60 + t10 - t514) * t456) * t350 / 0.2e1 + (t302 * t750 - t303 * t490 + (t466 * t58 + t12) * t456 + (t466 * t59 - t11 + t514) * t457) * t712;
t559 = t303 * rSges(6,1) + t302 * rSges(6,2);
t191 = rSges(6,3) * t457 + t559;
t701 = pkin(8) * t457;
t344 = pkin(4) * t661 + t701;
t646 = t191 + t344;
t805 = t810 * qJD(3) + t825 * t466 + t815 * t469 + t814 * t471;
t804 = t816 * t471 - t818 * t469 + (-t259 - t261) * t466 - t758 * qJD(3);
t803 = t807 * qJD(3) + t831 * t466 + t819 * t469 + t817 * t471;
t802 = t764 * t456 - t765 * t457;
t801 = t757 * t456 - t766 * t457;
t800 = -t813 * qJD(3) + t811 * t466;
t799 = (-t521 - t520 - t808) * t466 + t812 * t456;
t798 = t812 * t457 + t809 * t466 + t813 * t663;
t795 = -t490 * t364 - t750 * t528;
t252 = rSges(6,1) * t528 + t364 * rSges(6,2);
t794 = t252 * t350;
t793 = t252 * t351;
t791 = t762 * t466;
t790 = t763 * t466;
t211 = t302 * rSges(6,1) - t303 * rSges(6,2);
t213 = t304 * rSges(6,1) - t305 * rSges(6,2);
t789 = -t211 * t350 - t213 * t351;
t143 = rSges(6,1) * t240 + rSges(6,2) * t767;
t251 = rSges(6,1) * t364 - rSges(6,2) * t528;
t607 = qJD(4) * t471;
t458 = t469 * qJ(4);
t743 = t471 * pkin(3) + t458;
t340 = qJD(3) * t743 - t607;
t403 = pkin(3) * t469 - qJ(4) * t471;
t609 = qJD(3) * t469;
t585 = t457 * t609;
t654 = t466 * t471;
t515 = -t456 * t654 - t585;
t655 = t466 * t469;
t601 = t456 * t655;
t455 = qJD(4) * t469;
t414 = t457 * t455;
t608 = qJD(3) * t471;
t584 = t457 * t608;
t623 = qJ(4) * t584 + t414;
t174 = pkin(3) * t515 - qJ(4) * t601 + t623;
t332 = pkin(3) * t656 + t457 * t458;
t355 = t457 * pkin(2) + t456 * pkin(7);
t402 = pkin(7) * t658;
t472 = cos(qJ(1));
t464 = t472 * pkin(1);
t474 = qJD(1) ^ 2;
t470 = sin(qJ(1));
t703 = pkin(1) * t470;
t567 = qJDD(1) * t464 - t474 * t703;
t527 = t466 * (-pkin(2) * t663 + t402) + t465 * t355 + t567;
t605 = qJD(3) * qJD(4);
t756 = qJDD(4) * t469 + t471 * t605;
t493 = t465 * t332 + t527 + (t174 + t414) * t466 + t756 * t456;
t612 = qJD(3) * t456;
t628 = t305 * rSges(6,1) + t304 * rSges(6,2);
t193 = -rSges(6,3) * t456 + t628;
t429 = pkin(4) * t656;
t345 = -pkin(8) * t456 + t429;
t645 = t193 + t345;
t652 = t471 * qJD(3) ^ 2;
t231 = pkin(4) * t515 - pkin(8) * t658;
t651 = t129 * rSges(6,1) + t128 * rSges(6,2);
t84 = -rSges(6,3) * t658 + t651;
t691 = t231 + t84;
t14 = t493 - t340 * t612 - t334 * t403 + (-t334 * t469 - t456 * t652) * pkin(4) - t143 * t350 - t220 * t251 + t691 * t466 + t645 * t465;
t786 = -g(2) + t14;
t588 = t456 * t609;
t377 = rSges(5,1) * t588;
t440 = t456 * rSges(5,2);
t560 = t471 * rSges(5,1) + t469 * rSges(5,3);
t587 = t456 * t608;
t172 = rSges(5,3) * t587 - t377 + (t457 * t560 + t440) * t466;
t335 = -qJDD(3) * t457 + t383;
t404 = rSges(5,1) * t469 - rSges(5,3) * t471;
t519 = (-qJDD(1) * t470 - t472 * t474) * pkin(1);
t494 = t335 * t403 + t457 * t756 + t519;
t306 = t457 * t655 + t587;
t382 = pkin(3) * t588;
t583 = t456 * t455;
t562 = -t382 + t583;
t599 = t457 * t654;
t175 = pkin(3) * t599 + qJ(4) * t306 + t562;
t310 = t355 * t466;
t531 = -t175 - t310 - t583;
t624 = -t560 * qJD(3) - t340;
t571 = qJD(3) * t624;
t442 = t457 * rSges(5,2);
t275 = t456 * t560 - t442;
t328 = t743 * t456;
t446 = t457 * pkin(7);
t354 = pkin(2) * t456 - t446;
t626 = -t328 - t354;
t595 = -t275 + t626;
t44 = t335 * t404 + t457 * t571 + t595 * t465 + (-t172 + t531) * t466 + t494;
t785 = t44 - g(1);
t622 = rSges(5,2) * t658 + rSges(5,3) * t584;
t170 = rSges(5,1) * t515 - rSges(5,3) * t601 + t622;
t277 = rSges(5,1) * t656 + rSges(5,3) * t657 + t440;
t616 = -t403 - t404;
t45 = t170 * t466 + t277 * t465 + t334 * t616 + t456 * t571 + t493;
t784 = t45 - g(2);
t593 = rSges(4,1) * t588 + rSges(4,2) * t306;
t744 = rSges(4,1) * t656 + t456 * rSges(4,3);
t173 = t466 * t744 - t593;
t698 = rSges(4,1) * t471;
t409 = -rSges(4,2) * t469 + t698;
t368 = t409 * qJD(3);
t405 = rSges(4,1) * t469 + rSges(4,2) * t471;
t614 = rSges(4,2) * t662 + t457 * rSges(4,3);
t276 = rSges(4,1) * t661 - t614;
t631 = -t276 - t354;
t64 = -t368 * t611 + t335 * t405 + (-t173 - t310) * t466 + t631 * t465 + t519;
t783 = t64 - g(1);
t307 = t584 - t601;
t525 = -rSges(4,2) * t307 + rSges(4,3) * t658;
t171 = rSges(4,1) * t515 + t525;
t278 = -rSges(4,2) * t657 + t744;
t65 = t171 * t466 + t278 * t465 - t334 * t405 - t368 * t612 + t527;
t782 = t65 - g(2);
t309 = rSges(3,1) * t658 - rSges(3,2) * t663;
t352 = rSges(3,1) * t456 + rSges(3,2) * t457;
t781 = -t309 * t466 - t352 * t465 - g(1) + t519;
t353 = t457 * rSges(3,1) - rSges(3,2) * t456;
t670 = t352 * t466;
t780 = t353 * t465 - t466 * t670 - g(2) + t567;
t779 = t799 * t456 + t804 * t457;
t778 = -t798 * t456 + t803 * t457;
t777 = t804 * t456 - t799 * t457;
t776 = t803 * t456 + t798 * t457;
t775 = t801 * qJD(3) + t790;
t774 = t802 * qJD(3) + t791;
t773 = t808 * qJD(3) + t816 * t469 + t818 * t471;
t772 = t809 * qJD(3) + t817 * t469 - t819 * t471;
t771 = -t800 * t456 + t805 * t457;
t770 = t805 * t456 + t800 * t457;
t254 = t466 * t276;
t339 = t466 * t354;
t761 = -rSges(4,1) * t585 + t254 + t339 + t402 + t525;
t760 = t659 + t821;
t696 = pkin(1) * qJD(1);
t602 = t470 * t696;
t280 = -t602 - t670;
t299 = t466 * t328;
t745 = t299 + t339;
t755 = t466 * t275 + t622 + t745;
t574 = -pkin(4) * t469 - t403;
t495 = -t351 * t251 + t574 * t611 + t414;
t754 = t646 * t466 - t495 + t651 + t745;
t569 = t626 - t646;
t381 = pkin(4) * t588;
t621 = -pkin(8) * t663 - t381;
t232 = pkin(4) * t599 + t621;
t85 = rSges(6,1) * t131 + rSges(6,2) * t130 - rSges(6,3) * t663;
t690 = t232 + t85;
t13 = -t340 * t611 - t351 * t143 + t221 * t251 + (t335 * t469 - t457 * t652) * pkin(4) + t569 * t465 + (t531 - t690) * t466 + t494;
t533 = -pkin(2) - t743;
t702 = pkin(4) * t471;
t511 = t533 - t702;
t62 = t466 * t569 + t495 - t602;
t341 = t403 * t612;
t603 = t472 * t696;
t516 = t583 + t603;
t510 = -t341 + t516;
t625 = -t332 - t345;
t597 = -t193 + t625;
t570 = t355 - t597;
t747 = t350 * t251 + t381;
t63 = t466 * t570 + t510 - t747;
t705 = -rSges(6,3) - pkin(8);
t720 = -pkin(3) - pkin(4);
t752 = (t609 * t63 * t720 + t13 * t705 + (t511 * t62 + t63 * t705) * t466) * t457 + (t13 * t511 + t62 * (-qJ(4) * t608 - t455) + (-t62 * pkin(7) + t511 * t63) * t466 - g(1) * (t471 * t720 - pkin(2) - t458)) * t456;
t568 = t332 + t355;
t742 = t568 + t277;
t699 = g(2) * t456;
t734 = (g(1) * t457 + t699) * t469;
t218 = t278 + t355;
t729 = -t218 * t466 + t405 * t612;
t728 = -t404 * t612 + t466 * t742;
t634 = -Icges(4,2) * t661 + t267 - t422;
t638 = t395 * t456 + t263;
t725 = -t469 * t634 - t471 * t638;
t636 = -t552 * t456 + t265;
t640 = -t393 * t456 + t257;
t724 = -t469 * t636 + t471 * t640;
t722 = m(5) / 0.2e1;
t721 = m(6) / 0.2e1;
t717 = t334 / 0.2e1;
t716 = t335 / 0.2e1;
t706 = -rSges(5,1) - pkin(3);
t694 = t466 * t62;
t692 = -rSges(5,3) - qJ(4);
t589 = t405 * t611;
t517 = -t589 - t602;
t136 = t466 * t631 + t517;
t675 = t136 * t457;
t667 = t388 * t466;
t664 = t390 * t466;
t639 = -Icges(5,1) * t657 + t258 + t421;
t637 = -t395 * t457 - t264;
t635 = -t552 * t457 + t266;
t633 = -t391 * t457 + t268;
t630 = -t277 - t332;
t629 = t456 * t328 + t457 * t332;
t418 = qJ(4) * t656;
t329 = -pkin(3) * t657 + t418;
t627 = t466 * t329 + t456 * t607;
t620 = -t552 + t555;
t619 = t386 - t393;
t618 = -t391 + t396;
t617 = t395 + t553;
t615 = -t743 - t560;
t613 = t456 ^ 2 + t457 ^ 2;
t598 = t456 * t175 + (t174 + t299) * t457;
t416 = qJ(4) * t661;
t325 = -pkin(3) * t662 + t416;
t596 = t325 * t612 + t329 * t611 + t455;
t592 = t402 + t623;
t591 = t457 * t706;
t580 = -pkin(2) - t698;
t578 = -t612 / 0.2e1;
t577 = t612 / 0.2e1;
t576 = -t611 / 0.2e1;
t575 = t611 / 0.2e1;
t572 = -t259 + t671;
t566 = -t251 + t574;
t563 = -t341 + t583;
t561 = t328 * t612 + t332 * t611 - t607;
t410 = rSges(2,1) * t472 - rSges(2,2) * t470;
t406 = rSges(2,1) * t470 + rSges(2,2) * t472;
t557 = t456 * t63 + t457 * t62;
t137 = t603 - t729;
t546 = -t137 * t456 - t675;
t539 = t276 * t456 + t278 * t457;
t534 = -pkin(4) * t608 - t143 - t340;
t530 = t446 - t559;
t529 = t611 * t616 + t414;
t526 = t592 - t602;
t518 = -qJDD(4) * t471 + t174 * t611 + t175 * t612 + t334 * t328 + t469 * t605;
t512 = -t191 + t446 - t701;
t509 = -t469 * t635 + t471 * t639;
t508 = -t469 * t633 + t471 * t637;
t217 = t456 * t580 + t446 + t614;
t507 = t469 * t692 + t471 * t706 - pkin(2);
t506 = t382 - t85 - t621;
t504 = t529 - t602;
t503 = (t469 * t619 + t471 * t620) * t466;
t502 = (-t469 * t617 + t471 * t618) * t466;
t106 = t456 * t705 + t429 + t568 + t628;
t178 = t456 * t507 + t442 + t446;
t479 = (t580 * t675 + (t136 * (-rSges(4,3) - pkin(7)) + t137 * t580) * t456) * t466;
t477 = t23 * t712 + (t98 + t96) * t220 / 0.2e1 + (t97 + t95) * t221 / 0.2e1 + (t30 + t27) * t350 / 0.2e1 + (((t110 - t224 + (t260 + t672) * t457 + t820) * t457 + (t107 - t748 + t760) * t456) * qJD(3) + t791) * t575 - (t31 + t23 + t26) * t351 / 0.2e1 + (t811 * qJD(3) - t141 * t528 + t142 * t364 + t240 * t248 + t767 * t245 + t814 * t469 - t815 * t471) * t466 + (-t807 + t762) * t717 + (-t758 + t763) * t716 + (t771 + t772) * t577 + (-t245 * t528 + t248 * t364 + Icges(3,3) - t810) * t465 + (((t457 * t572 - t760 + t764) * t457 + (t456 * t572 - t237 + t564 + t573 + t765) * t456) * qJD(3) + t775 - t790) * t578 + (t770 - t773 + t774) * t576;
t103 = t466 * t595 + t504;
t104 = t510 + t728;
t476 = (t103 * t507 * t457 + (t103 * (-rSges(5,2) - pkin(7)) + t104 * (t533 - t560)) * t456) * t466 + (t103 * t661 * t692 + t104 * t469 * t591) * qJD(3);
t426 = rSges(5,3) * t656;
t424 = rSges(5,3) * t661;
t415 = t457 * t607;
t333 = t403 * t663;
t331 = t405 * t457;
t330 = -rSges(5,1) * t657 + t426;
t327 = t405 * t456;
t326 = -rSges(5,1) * t662 + t424;
t312 = t613 * t609;
t281 = t353 * t466 + t603;
t144 = t539 * qJD(3);
t102 = (t275 * t456 + t277 * t457) * qJD(3) + t561;
t56 = t191 * t350 + t193 * t351 + (t344 * t456 + t345 * t457) * qJD(3) + t561;
t35 = t275 * t334 + t630 * t335 + (t170 * t457 + t172 * t456) * qJD(3) + t518;
t8 = t191 * t220 - t193 * t221 + t334 * t344 + t350 * t85 + t351 * t84 + t625 * t335 + (t231 * t457 + t232 * t456) * qJD(3) + t518;
t1 = [Icges(2,3) * qJDD(1) + t477 + (t780 * (t353 + t464) + t781 * (-t352 - t703) + (-t309 - t603 + t281) * t280) * m(3) + ((t406 ^ 2 + t410 ^ 2) * qJDD(1) + g(1) * t406 - g(2) * t410) * m(2) + (-g(1) * (t512 - t703) + t13 * (t530 - t703) + t62 * (t506 - t603) + (t62 + t602 + t526 + t754) * t63 + t786 * (t464 + t106) + t752) * m(6) + (t103 * (t377 + t382 - t516) + t476 + t784 * (t464 + t742) + t785 * (t178 - t703) + (t526 + t103 - t504 + t755) * t104) * m(5) + (t136 * (t593 - t603) + t479 + t782 * (t218 + t464) + t783 * (t217 - t703) + (-t602 + t136 - t517 + t761) * t137) * m(4); t477 + (-g(1) * t512 + t13 * t530 + t570 * t694 + (t592 + t754) * t63 + (t563 - t747 + t506) * t62 + t786 * t106 + t752) * m(6) + (t476 + t784 * t742 + t785 * t178 + (t592 - t529 + t755) * t104 + (t377 - t562 + t563 + t728) * t103) * m(5) + (t479 + t782 * t218 + t783 * t217 + (t589 + t761) * t137 + (-t729 + t593) * t136) * m(4) + (-t280 * t309 - t281 * t670 + (t280 * t466 + t780) * t353 + (t281 * t466 - t781) * t352) * m(3); t802 * t717 + t801 * t716 + (t771 * t466 + t762 * t465 + t765 * t335 + t764 * t334 + (t778 * t456 + t779 * t457) * qJD(3)) * t456 / 0.2e1 - (t770 * t466 + t763 * t465 + t766 * t335 + t757 * t334 + (t776 * t456 + t777 * t457) * qJD(3)) * t457 / 0.2e1 + (-t456 * t807 + t457 * t758) * t465 / 0.2e1 + ((-t466 * t807 + t773) * t457 + (-t466 * t758 + t772) * t456) * t466 / 0.2e1 + t775 * t663 / 0.2e1 + t774 * t658 / 0.2e1 + ((-t612 * t665 + t664) * t456 + (t503 + (-t724 * t457 + (t666 + t509) * t456) * qJD(3)) * t457 + (-t612 * t668 + t667) * t456 + (t502 + (-t725 * t457 + (t669 + t508) * t456) * qJD(3)) * t457) * t578 + ((t466 * t764 + t779) * t457 + (t466 * t765 + t778) * t456) * t577 + ((t466 * t757 + t777) * t457 + (t466 * t766 + t776) * t456) * t576 + ((-t611 * t666 - t664) * t457 + (t503 + (t509 * t456 + (t665 - t724) * t457) * qJD(3)) * t456 + (-t611 * t669 - t667) * t457 + (t502 + (t508 * t456 + (t668 - t725) * t457) * qJD(3)) * t456) * t575 + (((t617 - t619) * t471 + (t618 + t620) * t469) * t466 + (((-t634 - t636) * t457 + (t633 + t635) * t456) * t471 + ((t638 - t640) * t457 + (t637 + t639) * t456) * t469) * qJD(3) - t795) * t708 + (-t62 * (t415 - t793) - t63 * (t627 - t794) - t56 * (t596 + t789) - (t62 * (t211 - t325) + t63 * (-pkin(4) * t657 - t213)) * t466 - (-t557 * t743 + (-t469 * t56 * t613 - t471 * t557) * pkin(4)) * qJD(3) + t62 * t333 + t8 * t629 + t56 * t598 + (t14 * t566 + t63 * t534 + t8 * t646 + t56 * t690 + (t62 * t251 + t56 * t597) * t466) * t456 + (t13 * t566 + t62 * t534 + t8 * t645 + t56 * t691 + (t56 * t646 + t566 * t63) * t466) * t457 - g(1) * (t418 - t213) - g(2) * (t416 - t211) - g(3) * (t252 + t743 + t702) - t720 * t734) * m(6) + (-t103 * (t415 + (-t325 - t326) * t466) - t104 * (t330 * t466 + t627) - t102 * t596 - ((t102 * t330 + t103 * t615) * t457 + (t102 * t326 + t104 * t615) * t456) * qJD(3) - g(1) * (t418 + t426) - g(2) * (t416 + t424) + g(3) * t615 - (g(1) * t591 + t699 * t706) * t469 + t103 * t333 + t35 * t629 + t102 * t598 + (t44 * t616 + t103 * t624 + t35 * t277 + t102 * t170 + (t102 * t275 + t104 * t616) * t466) * t457 + (t45 * t616 + t104 * t624 + t35 * t275 + t102 * t172 + (t102 * t630 + t103 * t404) * t466) * t456) * m(5) + (-(t136 * t327 - t137 * t331) * t466 - (t144 * (-t327 * t456 - t331 * t457) + t546 * t409) * qJD(3) + g(1) * t331 + g(2) * t327 - g(3) * t409 + (t276 * t334 - t278 * t335 + (t171 * t457 + t173 * t456) * qJD(3)) * t539 + t144 * ((t171 + t254) * t457 + (-t278 * t466 + t173) * t456) + t546 * t368 + ((-t137 * t466 - t64) * t457 + (t136 * t466 - t65) * t456) * t405) * m(4) - t806; (-m(5) - m(6)) * (-g(3) * t471 + t734) - m(5) * (t102 * t312 + t103 * t307 + t104 * t306) - m(6) * (t306 * t63 + t307 * t62 + t312 * t56) + 0.2e1 * ((t103 * t611 + t104 * t612 - t35) * t722 + (t611 * t62 + t612 * t63 - t8) * t721) * t471 + 0.2e1 * ((qJD(3) * t102 - t103 * t663 + t104 * t658 + t44 * t457 + t45 * t456) * t722 + (qJD(3) * t56 + t13 * t457 + t14 * t456 - t62 * t663 + t63 * t658) * t721) * t469; t795 * t708 + (t8 * (-t191 * t456 - t193 * t457) + t557 * t143 + ((t466 * t63 + t13) * t457 + (t14 - t694) * t456) * t251 - t62 * (-t211 * t466 + t793) - t63 * (t213 * t466 + t794) - g(1) * t213 - g(2) * t211 + g(3) * t252 + ((-t191 * t466 - t84) * t457 + (t193 * t466 - t85) * t456 + t789) * t56) * m(6) + t806;];
tau = t1;
