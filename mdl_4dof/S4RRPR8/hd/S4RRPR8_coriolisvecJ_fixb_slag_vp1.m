% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR8_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:48
% EndTime: 2019-12-31 17:08:15
% DurationCPUTime: 23.41s
% Computational Cost: add. (7217->663), mult. (19560->858), div. (0->0), fcn. (18067->6), ass. (0->336)
t672 = Icges(4,4) + Icges(3,5);
t351 = sin(qJ(1));
t353 = cos(qJ(1));
t350 = sin(qJ(2));
t336 = Icges(4,5) * t350;
t352 = cos(qJ(2));
t415 = Icges(4,1) * t352 + t336;
t671 = -t351 * t415 + t672 * t353;
t515 = t350 * t351;
t316 = Icges(3,4) * t515;
t513 = t351 * t352;
t656 = Icges(3,1) * t513 - t316 - t671;
t412 = Icges(4,3) * t352 - t336;
t540 = Icges(3,4) * t350;
t668 = Icges(3,2) * t352 + t412 + t540;
t535 = Icges(4,5) * t352;
t284 = Icges(4,1) * t350 - t535;
t339 = Icges(3,4) * t352;
t667 = Icges(3,1) * t350 + t284 + t339;
t279 = Icges(3,5) * t352 - Icges(3,6) * t350;
t197 = Icges(3,3) * t351 + t279 * t353;
t281 = Icges(4,4) * t352 + Icges(4,6) * t350;
t199 = Icges(4,2) * t351 + t281 * t353;
t670 = t197 + t199;
t203 = Icges(4,4) * t351 + t353 * t415;
t287 = Icges(3,1) * t352 - t540;
t205 = Icges(3,5) * t351 + t287 * t353;
t655 = t203 + t205;
t277 = Icges(4,3) * t350 + t535;
t413 = -Icges(3,2) * t350 + t339;
t669 = t277 - t413;
t662 = (Icges(3,6) - Icges(4,6)) * t352 + t672 * t350;
t666 = t287 + t415;
t526 = Icges(3,3) * t353;
t196 = Icges(3,5) * t513 - Icges(3,6) * t515 - t526;
t194 = -Icges(4,6) * t353 + t277 * t351;
t529 = Icges(3,6) * t353;
t200 = Icges(3,4) * t513 - Icges(3,2) * t515 - t529;
t521 = t200 * t350;
t633 = -t194 * t350 - t656 * t352 + t521;
t665 = -t196 * t353 - t633 * t351;
t469 = qJD(2) * t350;
t467 = qJD(2) * t352;
t663 = t655 * t350;
t647 = -t350 * t668 + t667 * t352;
t661 = t467 * t668 + t469 * t667;
t511 = t352 * t353;
t315 = Icges(4,5) * t511;
t514 = t350 * t353;
t528 = Icges(4,6) * t351;
t195 = Icges(4,3) * t514 + t315 + t528;
t660 = t195 * t514 + t351 * t670 + t655 * t511;
t198 = -Icges(4,2) * t353 + t281 * t351;
t183 = t351 * t198;
t659 = -t194 * t514 - t351 * t196 - t511 * t656 - t183;
t658 = -t194 + t200;
t201 = Icges(3,6) * t351 + t353 * t413;
t657 = -t195 + t201;
t654 = t669 * qJD(2);
t653 = t666 * qJD(2);
t652 = -t279 - t281;
t592 = t662 * t353;
t593 = t662 * t351;
t460 = qJD(2) - qJD(4);
t274 = t460 * t351;
t466 = qJD(2) * t353;
t275 = -qJD(4) * t353 + t466;
t349 = sin(qJ(4));
t560 = cos(qJ(4));
t455 = t350 * t560;
t266 = -t352 * t349 + t455;
t221 = t266 * t351;
t395 = t350 * t349 + t352 * t560;
t222 = t395 * t351;
t111 = Icges(5,5) * t222 + Icges(5,6) * t221 + Icges(5,3) * t353;
t537 = Icges(5,4) * t222;
t114 = Icges(5,2) * t221 + Icges(5,6) * t353 + t537;
t538 = Icges(5,4) * t221;
t118 = -Icges(5,1) * t222 - Icges(5,5) * t353 - t538;
t223 = t349 * t511 - t353 * t455;
t224 = t395 * t353;
t39 = -t111 * t351 - t223 * t114 - t118 * t224;
t113 = Icges(5,5) * t224 - Icges(5,6) * t223 - Icges(5,3) * t351;
t193 = Icges(5,4) * t224;
t116 = -Icges(5,2) * t223 - Icges(5,6) * t351 + t193;
t192 = Icges(5,4) * t223;
t119 = Icges(5,1) * t224 - Icges(5,5) * t351 - t192;
t429 = t113 * t351 + t223 * t116 - t224 * t119;
t161 = Icges(5,5) * t266 - Icges(5,6) * t395;
t253 = Icges(5,4) * t266;
t164 = -Icges(5,2) * t395 + t253;
t252 = Icges(5,4) * t395;
t167 = Icges(5,1) * t266 - t252;
t61 = -t161 * t351 - t164 * t223 + t167 * t224;
t12 = t61 * qJD(1) - t274 * t429 - t275 * t39;
t470 = qJD(1) * t353;
t159 = t460 * t395;
t617 = t460 * t266;
t77 = Icges(5,5) * t159 + Icges(5,6) * t617;
t78 = Icges(5,4) * t159 + Icges(5,2) * t617;
t79 = Icges(5,1) * t159 + Icges(5,4) * t617;
t471 = qJD(1) * t351;
t89 = t159 * t353 - t266 * t471;
t90 = -t353 * t617 - t395 * t471;
t13 = -t161 * t470 + t164 * t89 + t167 * t90 - t223 * t78 + t224 * t79 - t351 * t77;
t91 = t159 * t351 + t266 * t470;
t92 = qJD(1) * t224 - t351 * t617;
t14 = -t161 * t471 + t164 * t91 + t167 * t92 + t221 * t78 + t222 * t79 + t353 * t77;
t55 = Icges(5,4) * t92 + Icges(5,2) * t91 - Icges(5,6) * t471;
t57 = Icges(5,1) * t92 + Icges(5,4) * t91 - Icges(5,5) * t471;
t15 = t114 * t617 - t118 * t159 + t266 * t57 - t395 * t55;
t54 = Icges(5,4) * t90 + Icges(5,2) * t89 - Icges(5,6) * t470;
t56 = Icges(5,1) * t90 + Icges(5,4) * t89 - Icges(5,5) * t470;
t16 = t116 * t617 + t119 * t159 + t266 * t56 - t395 * t54;
t434 = qJD(1) * t460;
t257 = t351 * t434;
t258 = t353 * t434;
t37 = t111 * t353 + t114 * t221 - t118 * t222;
t372 = qJD(1) * (Icges(5,1) * t395 + t164 + t253) + t274 * (Icges(5,1) * t223 + t116 + t193) - t275 * (-Icges(5,1) * t221 + t114 + t537);
t38 = t353 * t113 + t221 * t116 + t222 * t119;
t381 = qJD(1) * (-Icges(5,5) * t395 - Icges(5,6) * t266) + (-Icges(5,5) * t221 + Icges(5,6) * t222) * t275 - (Icges(5,5) * t223 + Icges(5,6) * t224) * t274;
t557 = -qJD(1) / 0.2e1;
t565 = t275 / 0.2e1;
t53 = Icges(5,5) * t92 + Icges(5,6) * t91 - Icges(5,3) * t471;
t6 = -t111 * t470 + t114 * t89 - t118 * t90 - t223 * t55 + t224 * t57 - t351 * t53;
t60 = t161 * t353 + t164 * t221 + t167 * t222;
t503 = Icges(5,2) * t266 - t167 + t252;
t574 = -t274 * (-Icges(5,2) * t224 + t119 - t192) + t275 * (-Icges(5,2) * t222 - t118 + t538);
t607 = -qJD(1) * t503 - t574;
t609 = t39 / 0.2e1 + t38 / 0.2e1;
t62 = -t114 * t395 - t118 * t266;
t63 = -t116 * t395 + t119 * t266;
t52 = Icges(5,5) * t90 + Icges(5,6) * t89 - Icges(5,3) * t470;
t7 = -t113 * t470 + t116 * t89 + t119 * t90 - t223 * t54 + t224 * t56 - t351 * t52;
t8 = -t111 * t471 + t114 * t91 - t118 * t92 + t221 * t55 + t222 * t57 + t353 * t53;
t9 = -t113 * t471 + t116 * t91 + t119 * t92 + t221 * t54 + t222 * t56 + t353 * t52;
t646 = (-t15 * t353 + t16 * t351 + (t351 * t62 + t353 * t63) * qJD(1)) * t557 - t351 * (qJD(1) * t13 + t274 * t7 - t275 * t6) / 0.2e1 + t353 * (qJD(1) * t14 + t274 * t9 - t275 * t8) / 0.2e1 - (qJD(1) * t60 + t274 * t38 - t275 * t37) * t471 / 0.2e1 - t12 * t470 / 0.2e1 + (t351 * t429 + t353 * t609) * t258 + (-t351 * t609 + t353 * t37) * t257 - (-t223 * t607 - t224 * t372 + (-qJD(1) * t429 - t6) * t353 + (qJD(1) * t39 - t381 + t7) * t351) * t274 / 0.2e1 + (t221 * t607 - t222 * t372 + (qJD(1) * t37 + t9) * t351 + (qJD(1) * t38 + t381 - t8) * t353) * t565;
t522 = t198 * t353;
t625 = -t522 + t665;
t624 = -t200 * t514 - t659;
t623 = -t201 * t514 + t660;
t424 = -t195 * t515 + t199 * t353 - t203 * t513;
t174 = t205 * t513;
t436 = t197 * t353 - t174;
t72 = -t201 * t515 - t436;
t645 = -t424 + t72;
t644 = t372 * t266;
t643 = t351 * t647 - t592;
t642 = t353 * t647 + t593;
t641 = t662 * qJD(2);
t520 = t201 * t350;
t640 = t195 * t350 + t352 * t655 - t520;
t637 = (-t353 * t668 + t655) * t351 - (-Icges(3,2) * t513 - t412 * t351 - t316 + t656) * t353;
t636 = (-t641 * t351 + (t633 + t670) * qJD(1)) * t353;
t635 = t653 * t352 + t654 * t350 + (-t350 * t667 - t352 * t668) * qJD(2) + t662 * qJD(1);
t634 = t658 * t353 + (-Icges(4,1) * t514 + t284 * t353 + t315 - t657) * t351;
t632 = t667 - t669;
t631 = -t668 + t666;
t628 = t647 * qJD(1) + qJD(2) * t652;
t622 = t642 * qJD(1);
t169 = -rSges(5,1) * t395 - rSges(5,2) * t266;
t619 = t169 * t274;
t618 = t169 * t275;
t420 = -rSges(5,1) * t222 - rSges(5,2) * t221;
t120 = rSges(5,3) * t353 - t420;
t271 = pkin(3) * t513 + pkin(6) * t353;
t506 = t120 + t271;
t614 = (t623 * t351 - t624 * t353) * qJD(2);
t613 = (t645 * t351 - t625 * t353) * qJD(2);
t612 = t643 * qJD(1);
t610 = -t641 * t353 + (-t279 * t351 - t198 + t526 - t640) * qJD(1);
t154 = -rSges(5,1) * t221 + rSges(5,2) * t222;
t156 = rSges(5,1) * t223 + rSges(5,2) * t224;
t606 = -t154 * t274 - t156 * t275;
t605 = 0.2e1 * qJD(2);
t604 = t612 + t613;
t603 = t614 + t622;
t602 = t633 * qJD(2) + t661 * t351 + ((t277 * t353 - t201 + t528) * t352 - t663) * qJD(1);
t601 = t640 * qJD(2) - t661 * t353 + ((-t351 * t413 + t194 + t529) * t352 + (-t287 * t351 + t671) * t350) * qJD(1);
t600 = -t628 * t351 + t635 * t353;
t599 = t635 * t351 + t628 * t353;
t598 = t350 * t656 + t352 * t658;
t551 = t350 * rSges(4,1);
t289 = -rSges(4,3) * t352 + t551;
t331 = qJD(3) * t350;
t596 = (qJD(2) * t289 - t331) * t351;
t595 = t352 * t657 + t663;
t292 = pkin(2) * t352 + qJ(3) * t350;
t247 = t292 * t351;
t226 = qJD(1) * t247;
t345 = t353 * pkin(5);
t296 = pkin(1) * t351 - t345;
t273 = qJD(1) * t296;
t594 = -t226 - t273;
t321 = pkin(3) * t511;
t272 = -pkin(6) * t351 + t321;
t297 = t353 * pkin(1) + t351 * pkin(5);
t591 = -t637 * t350 + t634 * t352;
t590 = (-t632 * t350 + t631 * t352) * qJD(1);
t589 = t522 + t660;
t341 = t351 * rSges(4,2);
t217 = rSges(4,1) * t511 + rSges(4,3) * t514 + t341;
t322 = pkin(2) * t511;
t251 = qJ(3) * t514 + t322;
t431 = t251 + t297;
t588 = t431 + t217;
t587 = t652 * qJD(1);
t572 = m(4) / 0.2e1;
t571 = m(5) / 0.2e1;
t562 = -rSges(4,1) - pkin(2);
t561 = -rSges(5,3) - pkin(6);
t559 = pkin(3) * t352;
t556 = qJD(1) / 0.2e1;
t555 = t90 * rSges(5,1) + t89 * rSges(5,2);
t554 = rSges(3,1) * t352;
t290 = rSges(3,1) * t350 + rSges(3,2) * t352;
t250 = t290 * t353;
t468 = qJD(2) * t351;
t453 = t290 * t468;
t340 = t351 * rSges(3,3);
t218 = rSges(3,1) * t511 - rSges(3,2) * t514 + t340;
t493 = t218 + t297;
t98 = qJD(1) * t493 - t453;
t553 = t250 * t98;
t477 = rSges(3,2) * t515 + t353 * rSges(3,3);
t216 = rSges(3,1) * t513 - t477;
t452 = t290 * t466;
t97 = -t452 + (-t216 - t296) * qJD(1);
t550 = t351 * t97;
t549 = t353 * t97;
t548 = -rSges(4,3) - qJ(3);
t383 = -t350 * t466 - t352 * t471;
t185 = pkin(3) * t383 - pkin(6) * t470;
t58 = -rSges(5,3) * t470 + t555;
t547 = t185 + t58;
t451 = t350 * t468;
t306 = pkin(3) * t451;
t186 = qJD(1) * t272 - t306;
t425 = rSges(5,1) * t92 + rSges(5,2) * t91;
t59 = -rSges(5,3) * t471 + t425;
t546 = t186 + t59;
t230 = t350 * t470 + t351 * t467;
t307 = pkin(2) * t451;
t448 = t351 * t331;
t125 = qJ(3) * t230 + qJD(1) * t322 - t307 + t448;
t270 = t297 * qJD(1);
t505 = -t125 - t270;
t496 = t224 * rSges(5,1) - t223 * rSges(5,2);
t494 = -t217 - t251;
t492 = t351 * t247 + t353 * t251;
t465 = qJD(3) * t352;
t225 = qJD(2) * t292 - t465;
t293 = rSges(4,1) * t352 + rSges(4,3) * t350;
t491 = -t293 * qJD(2) - t225;
t288 = pkin(2) * t350 - qJ(3) * t352;
t248 = t288 * t353;
t490 = -qJD(1) * t248 + t351 * t465;
t254 = t288 * t468;
t461 = qJD(2) * qJD(3);
t446 = t352 * t461;
t489 = qJD(1) * t254 + t353 * t446;
t488 = -t247 - t296;
t487 = -t251 - t272;
t482 = -t288 - t289;
t481 = -t292 - t293;
t464 = qJD(3) * t353;
t310 = t350 * t464;
t450 = t352 * t466;
t480 = qJ(3) * t450 + t310;
t479 = rSges(4,2) * t470 + rSges(4,3) * t450;
t454 = t350 * t471;
t478 = rSges(3,2) * t454 + rSges(3,3) * t470;
t476 = t351 ^ 2 + t353 ^ 2;
t124 = pkin(2) * t383 - qJ(3) * t454 + t480;
t459 = t353 * t124 + t351 * t125 + t247 * t470;
t122 = -rSges(5,3) * t351 + t496;
t458 = -t122 + t487;
t244 = t288 * t351;
t457 = -t244 * t468 - t248 * t466 + t331;
t329 = pkin(5) * t470;
t456 = t329 + t480;
t447 = -pkin(1) - t554;
t443 = -t468 / 0.2e1;
t440 = t466 / 0.2e1;
t439 = -pkin(3) * t350 - t288;
t437 = t353 * t482;
t435 = -t196 + t520;
t433 = t125 * t468 + t350 * t461 + (t124 + t226) * t466;
t256 = qJD(1) * (-pkin(1) * t471 + t329);
t432 = t351 * t446 + t256 + (t124 + t310) * qJD(1);
t170 = rSges(5,1) * t266 - rSges(5,2) * t395;
t430 = -t170 + t439;
t426 = t352 * t562 - pkin(1);
t423 = t247 * t468 + t251 * t466 - t465;
t421 = -rSges(3,2) * t350 + t554;
t397 = (-t559 * qJD(2) - t225) * qJD(2);
t398 = t439 * t466;
t80 = rSges(5,1) * t159 + rSges(5,2) * t617;
t17 = -t170 * t258 - t274 * t80 + t397 * t351 + (t398 + t547) * qJD(1) + t432;
t18 = t170 * t257 - t275 * t80 + t397 * t353 + ((pkin(3) * qJD(2) - qJD(3)) * t515 + t505 - t546) * qJD(1) + t489;
t419 = t17 * t351 + t18 * t353;
t375 = -t275 * t170 + t310 + t398;
t41 = (t488 - t506) * qJD(1) + t375;
t42 = t448 - t170 * t274 - t254 - t306 + (t297 - t458) * qJD(1);
t418 = -t351 * t42 - t353 * t41;
t417 = -t351 * t98 - t549;
t402 = -pkin(3) * t467 - t225 - t80;
t400 = -pkin(1) - t292;
t396 = qJD(2) * t437 + t310;
t246 = t290 * t351;
t245 = t289 * t351;
t93 = (t216 * t351 + t218 * t353) * qJD(2);
t380 = t400 - t559;
t343 = t353 * rSges(4,2);
t311 = t352 * t464;
t268 = t421 * qJD(2);
t255 = t288 * t471;
t249 = t289 * t353;
t231 = t450 - t454;
t229 = t476 * t469;
t215 = t293 * t351 - t343;
t153 = -qJD(2) * t246 + (t353 * t421 + t340) * qJD(1);
t152 = -qJD(2) * t245 + (t293 * t353 + t341) * qJD(1);
t151 = rSges(3,1) * t383 - rSges(3,2) * t450 + t478;
t150 = rSges(4,1) * t383 - rSges(4,3) * t454 + t479;
t68 = qJD(1) * t588 - t254 - t596;
t67 = (-t215 + t488) * qJD(1) + t396;
t66 = -t268 * t466 + (-t153 - t270 + t453) * qJD(1);
t65 = -t268 * t468 + t256 + (t151 - t452) * qJD(1);
t64 = (t215 * t351 + t217 * t353) * qJD(2) + t423;
t36 = t120 * t274 + t122 * t275 + (t271 * t351 + t272 * t353) * qJD(2) + t423;
t35 = t491 * t466 + (-t152 + t505 + t596) * qJD(1) + t489;
t34 = qJD(1) * t150 + (qJD(1) * t437 + t351 * t491) * qJD(2) + t432;
t21 = (t150 * t353 + t152 * t351 + (t215 * t353 + t351 * t494) * qJD(1)) * qJD(2) + t433;
t10 = t120 * t258 - t122 * t257 + t274 * t59 + t275 * t58 + (t185 * t353 + t186 * t351 + (t271 * t353 + t351 * t487) * qJD(1)) * qJD(2) + t433;
t1 = [t12 * t565 + (t62 + t60) * t257 / 0.2e1 + (t63 + t61) * t258 / 0.2e1 + (t13 + t16) * t274 / 0.2e1 + (((t72 - t174 + (t197 + t521) * t353 + t659) * t353 + (t589 + t625 - t665) * t351) * qJD(2) + t622) * t440 + (t647 * qJD(2) + t159 * t167 + t617 * t164 + t266 * t79 + t653 * t350 - t654 * t352 - t395 * t78) * qJD(1) + (-(-t506 * qJD(1) + t375 - t41 + t594) * t42 + t18 * (t345 + t420) + t41 * (t306 + t307 - t425) + t17 * (t321 + t431 + t496) + t42 * (t456 + t555) + (t18 * t561 + t42 * (-pkin(2) - pkin(3)) * t469) * t353 + (t18 * t380 + t41 * (-qJ(3) * t467 - t331) + t17 * t561) * t351 + ((t380 * t41 + t42 * t561) * t353 + (t41 * (-pkin(5) - t561) + t42 * t380) * t351) * qJD(1)) * m(5) + (-(-qJD(1) * t215 + t396 + t594 - t67) * t68 + t35 * (t343 + t345) + t67 * t307 + t34 * t588 + t68 * (t456 + t479) + (t68 * t562 * t469 + t67 * (t350 * t548 + t426) * qJD(1)) * t353 + (t35 * t426 + (-t67 * qJD(3) + t35 * t548) * t350 + t67 * (t352 * t548 + t551) * qJD(2) + (t67 * (-rSges(4,2) - pkin(5)) + t68 * (-t293 + t400)) * qJD(1)) * t351) * m(4) + (t66 * (t351 * t447 + t345 + t477) + t65 * t493 + t98 * (t329 + t478) + (t290 * t550 - t553) * qJD(2) + ((-pkin(1) - t421) * t549 + (t97 * (-rSges(3,3) - pkin(5)) + t98 * t447) * t351) * qJD(1) - (-qJD(1) * t216 - t273 - t452 - t97) * t98) * m(3) - (t14 + t12 + t15) * t275 / 0.2e1 + (t600 + t601) * t468 / 0.2e1 + (((t353 * t435 - t589 + t623) * t353 + (t351 * t435 - t183 + t424 + t436 + t624) * t351) * qJD(2) + t604 - t612) * t443 - (t599 - t602 + t603) * t466 / 0.2e1 + ((t598 + t643) * t351 + (t595 + t642) * t353) * qJD(2) * t556; (t602 * t353 + t601 * t351 + (t351 * t598 + t353 * t595) * qJD(1)) * t556 + ((-t468 * t592 - t587) * t351 + ((t351 * t593 + t591) * qJD(2) + t590) * t353) * t443 + ((-t466 * t593 + t587) * t353 + ((t353 * t592 + t591) * qJD(2) + t590) * t351) * t440 + (t41 * t255 + t10 * t492 + t36 * t459 + (t17 * t430 + t42 * t402 + t10 * t506 + t36 * t546 + (t41 * t170 + t36 * t458) * qJD(1)) * t351 + (t18 * t430 + t41 * t402 + t10 * (t122 + t272) + t36 * t547 + (t36 * t506 + t42 * t430) * qJD(1)) * t353 - t41 * (t311 + t618) - t42 * (t490 + t619) - t36 * (t457 - t606) - (t41 * (-t154 + t244) + t42 * (-pkin(3) * t514 + t156)) * qJD(1) - (t418 * t292 + (-t350 * t36 * t476 + t352 * t418) * pkin(3)) * qJD(2)) * m(5) + (t67 * t255 + t21 * t492 + t64 * t459 + (t35 * t482 + t67 * t491 + t21 * t217 + t64 * t150 + (t64 * t215 + t482 * t68) * qJD(1)) * t353 + (t34 * t482 + t68 * t491 + t21 * t215 + t64 * t152 + (t67 * t289 + t494 * t64) * qJD(1)) * t351 - t67 * (t311 + (t244 + t245) * qJD(1)) - t68 * (-qJD(1) * t249 + t490) - t64 * t457 - ((-t64 * t249 + t481 * t67) * t353 + (-t64 * t245 + t481 * t68) * t351) * qJD(2)) * m(4) + (-(t246 * t97 - t553) * qJD(1) - (t93 * (-t246 * t351 - t250 * t353) + t417 * t421) * qJD(2) + 0.2e1 * t93 * (t151 * t353 + t153 * t351 + (t216 * t353 - t218 * t351) * qJD(1)) + t417 * t268 + (-t65 * t351 - t66 * t353 + (-t353 * t98 + t550) * qJD(1)) * t290) * m(3) + (t600 * qJD(1) + (t623 * t470 + (t624 * qJD(1) + t610 * t351 - t636) * t351) * t605) * t351 / 0.2e1 - (t599 * qJD(1) + ((t645 * qJD(1) + t636) * t353 + (t625 * qJD(1) - t353 * t610) * t351) * t605) * t353 / 0.2e1 + (t644 - t574 * t395 + (t634 * t350 + t637 * t352) * qJD(2) + (t631 * t350 + t632 * t352 - t503 * t395) * qJD(1)) * t557 + (t604 + t613) * t471 / 0.2e1 + (t603 + t614) * t470 / 0.2e1 - t646; -m(4) * (t229 * t64 + t230 * t68 + t231 * t67) - m(5) * (t229 * t36 + t230 * t42 + t231 * t41) + 0.2e1 * ((t466 * t67 + t468 * t68 - t21) * t572 + (t41 * t466 + t42 * t468 - t10) * t571) * t352 + 0.2e1 * ((qJD(2) * t64 + t34 * t351 + t35 * t353 + t470 * t68 - t471 * t67) * t572 + (qJD(2) * t36 - t41 * t471 + t42 * t470 + t419) * t571) * t350; (-t395 * t607 - t644) * t557 + ((t41 * t80 - t10 * t122 + t36 * (-qJD(1) * t120 - t58)) * t353 + (t42 * t80 - t10 * t120 + t36 * (qJD(1) * t122 - t59)) * t351 + ((-t351 * t41 + t353 * t42) * qJD(1) + t419) * t170 - t41 * (qJD(1) * t154 - t618) - t42 * (-qJD(1) * t156 - t619) - t36 * t606) * m(5) + t646;];
tauc = t1(:);
