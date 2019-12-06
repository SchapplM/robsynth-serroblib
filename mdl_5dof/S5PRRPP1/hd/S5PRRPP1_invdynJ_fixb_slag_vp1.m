% Calculate vector of inverse dynamics joint torques for
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:12
% EndTime: 2019-12-05 16:06:51
% DurationCPUTime: 27.30s
% Computational Cost: add. (13108->607), mult. (12058->745), div. (0->0), fcn. (9365->6), ass. (0->325)
t685 = Icges(4,3) + Icges(5,3);
t328 = qJ(3) + pkin(8);
t322 = sin(t328);
t324 = cos(t328);
t330 = sin(qJ(3));
t331 = cos(qJ(3));
t684 = Icges(4,5) * t331 + Icges(5,5) * t324 - Icges(4,6) * t330 - Icges(5,6) * t322;
t224 = Icges(6,4) * t324 + Icges(6,6) * t322;
t327 = pkin(7) + qJ(2);
t321 = sin(t327);
t323 = cos(t327);
t154 = Icges(6,2) * t321 + t224 * t323;
t668 = t685 * t321 + t684 * t323;
t683 = t154 + t668;
t300 = Icges(6,5) * t322;
t388 = Icges(6,3) * t324 - t300;
t545 = Icges(5,4) * t322;
t681 = Icges(5,2) * t324 + t388 + t545;
t304 = Icges(5,4) * t324;
t229 = Icges(5,1) * t322 + t304;
t541 = Icges(6,5) * t324;
t654 = Icges(6,1) * t322 + t229 - t541;
t682 = t685 * t323;
t517 = t321 * t331;
t518 = t321 * t330;
t519 = t321 * t324;
t520 = t321 * t322;
t659 = -Icges(4,5) * t517 - Icges(5,5) * t519 + Icges(4,6) * t518 + Icges(5,6) * t520 + t682;
t391 = Icges(6,1) * t324 + t300;
t157 = -Icges(6,4) * t323 + t321 * t391;
t256 = Icges(5,4) * t520;
t542 = Icges(5,5) * t323;
t159 = Icges(5,1) * t519 - t256 - t542;
t676 = -t157 - t159;
t158 = Icges(6,4) * t321 + t323 * t391;
t230 = Icges(5,1) * t324 - t545;
t160 = Icges(5,5) * t321 + t230 * t323;
t675 = t158 + t160;
t536 = Icges(5,6) * t323;
t155 = Icges(5,4) * t519 - Icges(5,2) * t520 - t536;
t537 = Icges(4,6) * t323;
t167 = Icges(4,4) * t517 - Icges(4,2) * t518 - t537;
t680 = t155 * t322 + t167 * t330;
t546 = Icges(4,4) * t330;
t274 = Icges(4,1) * t331 - t546;
t170 = Icges(4,5) * t321 + t274 * t323;
t679 = -t160 * t519 - t170 * t517;
t220 = Icges(6,3) * t322 + t541;
t149 = -Icges(6,6) * t323 + t220 * t321;
t678 = t149 - t155;
t512 = t323 * t324;
t255 = Icges(6,5) * t512;
t516 = t322 * t323;
t535 = Icges(6,6) * t321;
t150 = Icges(6,3) * t516 + t255 + t535;
t389 = -Icges(5,2) * t322 + t304;
t156 = Icges(5,6) * t321 + t323 * t389;
t677 = -t150 + t156;
t674 = t220 - t389;
t673 = t230 + t391;
t672 = t681 * qJD(3);
t671 = t654 * qJD(3);
t280 = Icges(4,4) * t518;
t543 = Icges(4,5) * t323;
t169 = Icges(4,1) * t517 - t280 - t543;
t670 = t159 * t324 + t169 * t331 - t680;
t669 = Icges(4,5) * t330 + Icges(4,6) * t331 + (Icges(5,6) - Icges(6,6)) * t324 + (Icges(6,4) + Icges(5,5)) * t322;
t409 = -t150 * t520 + t154 * t323 - t158 * t519;
t325 = Icges(4,4) * t331;
t390 = -Icges(4,2) * t330 + t325;
t168 = Icges(4,6) * t321 + t323 * t390;
t667 = -t668 * t323 - t679;
t627 = -t156 * t520 - t168 * t518 + t667;
t609 = -t409 + t627;
t666 = t156 * t322 + t168 * t330;
t510 = t323 * t331;
t665 = t150 * t516 + t170 * t510 + t683 * t321 + t675 * t512;
t153 = -Icges(6,2) * t323 + t224 * t321;
t138 = t321 * t153;
t664 = -t149 * t516 - t169 * t510 + t659 * t321 + t676 * t512 - t138;
t271 = Icges(4,2) * t331 + t546;
t273 = Icges(4,1) * t330 + t325;
t651 = t271 * t330 - t273 * t331 + t681 * t322 - t654 * t324;
t663 = t672 * t323 + (t321 * t389 - t149 - t536) * qJD(2);
t662 = t672 * t321 + (t220 * t323 - t156 + t535) * qJD(2);
t661 = -t671 * t323 + (-t230 * t321 - t157 + t542) * qJD(2);
t660 = -t675 * qJD(2) + t671 * t321;
t658 = t674 * qJD(3);
t657 = t673 * qJD(3);
t514 = t323 * t153;
t385 = t149 * t322 + t157 * t324;
t600 = t321 * t385;
t50 = -t514 + t600;
t615 = t670 * t321 + t659 * t323 + t50;
t511 = t323 * t330;
t614 = -t155 * t516 - t167 * t511 - t664;
t613 = -t156 * t516 - t168 * t511 + t665;
t610 = -t167 * t331 - t169 * t330 + t676 * t322 + t678 * t324;
t608 = t168 * t331 + t170 * t330 + t675 * t322 + t677 * t324;
t653 = t224 + t684;
t652 = t669 * qJD(3);
t650 = -t271 * t331 - t273 * t330 - t654 * t322 - t324 * t681;
t649 = t150 * t322 + t170 * t331 + t675 * t324 - t666;
t648 = -t385 - t670;
t602 = t669 * t323;
t601 = t669 * t321;
t612 = -t321 * t651 - t602;
t611 = -t323 * t651 + t601;
t647 = t683 * qJD(2);
t646 = t323 ^ 2;
t645 = t654 - t674;
t644 = -t681 + t673;
t643 = -t323 * t681 + t675;
t642 = -Icges(5,2) * t519 - t388 * t321 - t256 - t676;
t641 = -Icges(6,1) * t516 - t229 * t323 + t255 - t677;
t640 = t654 * t321 - t678;
t246 = t390 * qJD(3);
t247 = t274 * qJD(3);
t639 = t669 * qJD(2) + t650 * qJD(3) - t246 * t330 + t247 * t331 + t658 * t322 + t657 * t324;
t363 = qJD(3) * t271;
t107 = -t323 * t363 + (-t321 * t390 + t537) * qJD(2);
t366 = qJD(3) * t273;
t109 = -t323 * t366 + (-t274 * t321 + t543) * qJD(2);
t638 = -t608 * qJD(3) - t107 * t330 + t109 * t331 + t663 * t322 + t661 * t324 + t647;
t108 = qJD(2) * t168 - t321 * t363;
t110 = qJD(2) * t170 - t321 * t366;
t589 = qJD(2) * t153;
t637 = t659 * qJD(2) - t610 * qJD(3) + t108 * t330 - t110 * t331 - t662 * t322 + t660 * t324 - t589;
t636 = t321 * t609 - t615 * t323;
t635 = t613 * t321 - t614 * t323;
t634 = t651 * qJD(2) + t653 * qJD(3);
t633 = t648 * qJD(2) - t652 * t321 + t647;
t632 = -t589 - t652 * t323 + (-t321 * t684 - t649 + t682) * qJD(2);
t631 = rSges(4,2) * t330;
t626 = rSges(6,3) + qJ(5);
t630 = t611 * qJD(2);
t629 = t612 * qJD(2);
t448 = qJD(2) * qJD(3);
t216 = -qJDD(3) * t323 + t321 * t448;
t509 = t331 * qJD(3) ^ 2;
t445 = pkin(3) * t509;
t450 = qJD(5) * t324;
t591 = t324 * pkin(4) + t322 * qJ(5);
t593 = t324 * rSges(6,1) + t322 * rSges(6,3);
t622 = t591 + t593;
t482 = -t622 * qJD(3) + t450;
t343 = -t445 + qJDD(5) * t322 + (t450 + t482) * qJD(3);
t311 = t323 * rSges(6,2);
t489 = t622 * t321 - t311;
t315 = t323 * pkin(6);
t239 = pkin(2) * t321 - t315;
t329 = -qJ(4) - pkin(6);
t285 = t323 * t329;
t326 = t331 * pkin(3);
t318 = t326 + pkin(2);
t470 = -t321 * t318 - t285;
t145 = t239 + t470;
t499 = t145 - t239;
t416 = -t489 + t499;
t451 = qJD(5) * t322;
t432 = t321 * t451;
t447 = qJD(2) * qJD(4);
t568 = pkin(3) * t330;
t440 = qJDD(4) * t321 + t216 * t568 + t323 * t447;
t232 = pkin(4) * t322 - qJ(5) * t324;
t233 = rSges(6,1) * t322 - rSges(6,3) * t324;
t474 = t232 + t233;
t314 = t321 * pkin(6);
t457 = qJD(2) * t321;
t453 = qJD(3) * t330;
t468 = t321 * pkin(3) * t453 + qJD(4) * t323;
t439 = t329 * t457 + t468;
t564 = pkin(2) - t318;
t103 = (-t323 * t564 - t314) * qJD(2) - t439;
t240 = t323 * pkin(2) + t314;
t214 = t240 * qJD(2);
t506 = -t103 - t214;
t307 = t321 * rSges(6,2);
t455 = qJD(3) * t321;
t456 = qJD(2) * t323;
t562 = t591 * t456 + (-qJD(3) * t232 + t451) * t321 - t233 * t455 + (t323 * t593 + t307) * qJD(2);
t2 = t474 * t216 + t416 * qJDD(2) + t343 * t323 + (-t432 + t506 - t562) * qJD(2) + t440;
t628 = -g(1) + t2;
t569 = rSges(6,1) + pkin(4);
t625 = t641 * t321 + t640 * t323;
t624 = t659 + t666;
t623 = t321 * t643 - t323 * t642;
t621 = qJD(3) * t636 + t629;
t620 = qJD(3) * t635 + t630;
t619 = t648 * qJD(3) - t108 * t331 - t110 * t330 + t660 * t322 + t662 * t324;
t618 = t649 * qJD(3) + t107 * t331 + t109 * t330 + t661 * t322 - t663 * t324;
t617 = t321 * t634 + t323 * t639;
t616 = t321 * t639 - t323 * t634;
t466 = t273 + t390;
t467 = -t271 + t274;
t607 = (-t322 * t645 + t324 * t644 - t330 * t466 + t331 * t467) * qJD(2);
t606 = t514 + t665;
t605 = t633 * t646 + (t638 * t321 + (-t632 + t637) * t323) * t321;
t604 = t637 * t646 + (t632 * t321 + (-t633 + t638) * t323) * t321;
t603 = t653 * qJD(2);
t599 = t322 * t569;
t218 = qJD(2) * t239;
t598 = qJD(2) * t145 - t218;
t313 = t324 * rSges(5,1);
t592 = -rSges(5,2) * t322 + t313;
t597 = t592 + t326;
t454 = qJD(3) * t323;
t434 = t324 * t454;
t596 = rSges(6,2) * t456 + t626 * t434;
t595 = t626 * t519;
t594 = t626 * t512;
t414 = -t474 - t568;
t249 = t323 * t451;
t291 = qJD(4) * t321;
t471 = t249 + t291;
t590 = t414 * t454 + t471;
t588 = t326 + t622;
t485 = -Icges(4,2) * t517 + t169 - t280;
t487 = t273 * t321 + t167;
t576 = -t330 * t485 - t331 * t487;
t575 = m(2) + m(3);
t574 = -m(5) - m(6);
t215 = qJDD(3) * t321 + t323 * t448;
t573 = t215 / 0.2e1;
t572 = t216 / 0.2e1;
t571 = t321 / 0.2e1;
t570 = -t323 / 0.2e1;
t567 = g(2) * t321;
t355 = -t322 * t454 - t324 * t457;
t438 = t322 * t457;
t563 = t569 * t355 - t626 * t438 + t249 + t596;
t561 = rSges(4,1) * t331;
t276 = rSges(4,1) * t330 + rSges(4,2) * t331;
t202 = t276 * t323;
t306 = t321 * rSges(4,3);
t172 = rSges(4,1) * t510 - rSges(4,2) * t511 + t306;
t135 = t172 + t240;
t94 = qJD(2) * t135 - t276 * t455;
t559 = t202 * t94;
t418 = t323 * t318 - t321 * t329;
t146 = t418 - t240;
t431 = -t145 * t455 + t146 * t454 + qJD(1);
t488 = t569 * t512 + t626 * t516 + t307;
t25 = -t450 + (t321 * t489 + t323 * t488) * qJD(3) + t431;
t558 = t25 * t322;
t305 = t321 * rSges(5,3);
t436 = t276 * t454;
t465 = rSges(4,2) * t518 + t323 * rSges(4,3);
t171 = rSges(4,1) * t517 - t465;
t483 = -t171 - t239;
t93 = qJD(2) * t483 - t436;
t556 = t321 * t93;
t555 = t323 * t93;
t234 = rSges(5,1) * t322 + rSges(5,2) * t324;
t356 = -t234 - t568;
t353 = t356 * t454 + t291;
t162 = rSges(5,1) * t519 - rSges(5,2) * t520 - t323 * rSges(5,3);
t442 = -t162 + t499;
t48 = qJD(2) * t442 + t353;
t553 = t48 * t234;
t502 = -t321 * t145 + t323 * t146;
t164 = rSges(5,1) * t512 - rSges(5,2) * t516 + t305;
t498 = -t146 - t164;
t486 = -t273 * t323 - t168;
t484 = -t271 * t323 + t170;
t481 = -t569 * t520 + t595;
t480 = -t569 * t516 + t594;
t289 = pkin(6) * t456;
t479 = qJD(2) * (-pkin(2) * t457 + t289) + qJDD(2) * t240;
t472 = rSges(5,2) * t438 + rSges(5,3) * t456;
t469 = rSges(4,3) * t456 + t457 * t631;
t464 = t321 ^ 2 + t646;
t452 = qJD(3) * t331;
t446 = pkin(3) * t511;
t433 = t323 * t453;
t102 = -pkin(3) * t433 - t289 + t291 + (t321 * t564 - t285) * qJD(2);
t444 = t323 * t102 + t321 * t103 - t145 * t456;
t443 = pkin(3) * t452;
t441 = -t146 - t488;
t430 = -pkin(2) - t561;
t427 = -t455 / 0.2e1;
t426 = t455 / 0.2e1;
t425 = -t454 / 0.2e1;
t424 = t454 / 0.2e1;
t417 = t102 * t454 + t103 * t455 - t215 * t145 + qJDD(1);
t415 = t464 * t568;
t211 = t592 * qJD(3);
t412 = -t211 - t443;
t235 = rSges(3,1) * t323 - rSges(3,2) * t321;
t231 = rSges(3,1) * t321 + rSges(3,2) * t323;
t277 = t561 - t631;
t37 = qJD(2) * t416 + t590;
t38 = (-qJD(3) * t474 + t451) * t321 + (t240 - t441) * qJD(2) - t468;
t401 = t321 * t38 + t323 * t37;
t394 = -t321 * t94 - t555;
t113 = -rSges(4,2) * t323 * t452 + (-t331 * t457 - t433) * rSges(4,1) + t469;
t201 = t276 * t321;
t114 = -qJD(3) * t201 + (t277 * t323 + t306) * qJD(2);
t387 = t113 * t323 + t114 * t321;
t378 = t171 * t321 + t172 * t323;
t371 = -t443 + t482;
t369 = -t568 - t599;
t188 = t234 * t321;
t368 = t401 * t324;
t354 = t323 * t369;
t352 = qJD(2) * t102 + qJDD(2) * t146 - qJDD(4) * t323 + t321 * t447 + t479;
t349 = -t330 * t484 + t331 * t486;
t348 = t323 * t356;
t344 = -t318 - t622;
t248 = t277 * qJD(3);
t192 = t234 * t323;
t98 = -qJD(3) * t188 + (t323 * t592 + t305) * qJD(2);
t96 = rSges(5,1) * t355 - rSges(5,2) * t434 + t472;
t74 = qJD(3) * t378 + qJD(1);
t49 = -t234 * t455 + (t240 - t498) * qJD(2) - t468;
t43 = (t162 * t321 + t164 * t323) * qJD(3) + t431;
t42 = qJD(2) * t113 + qJDD(2) * t172 - t215 * t276 - t248 * t455 + t479;
t41 = -t248 * t454 + t216 * t276 + t483 * qJDD(2) + (-t114 - t214) * qJD(2);
t34 = qJD(3) * t387 + t171 * t215 - t172 * t216 + qJDD(1);
t22 = -t211 * t455 + qJD(2) * t96 + qJDD(2) * t164 - t215 * t234 + (-t215 * t330 - t321 * t509) * pkin(3) + t352;
t21 = t216 * t234 + (-qJD(3) * t211 - t445) * t323 + t442 * qJDD(2) + (-t98 + t506) * qJD(2) + t440;
t4 = t162 * t215 + t498 * t216 + (t321 * t98 + t323 * t96) * qJD(3) + t417;
t3 = t488 * qJDD(2) + t414 * t215 + (t249 + t563) * qJD(2) + t343 * t321 + t352;
t1 = -qJDD(5) * t324 + t489 * t215 + t441 * t216 + (t321 * t562 + t323 * t563 + t451) * qJD(3) + t417;
t5 = [t575 * qJDD(1) + m(4) * t34 + m(5) * t4 + m(6) * t1 + (-m(4) + t574 - t575) * g(3); -m(3) * (-g(1) * t231 + g(2) * t235) + (((t50 - t600 + t606) * t321 + ((t668 + t680) * t323 + t627 + t664 + t679) * t323) * qJD(3) + t630) * t424 + (-t651 * qJD(3) + t246 * t331 + t247 * t330 + t657 * t322 - t658 * t324) * qJD(2) + (t37 * (-t432 + t439) + t38 * (t471 + t596) + (t38 * t354 + (-t324 * t626 + t599) * t321 * t37) * qJD(3) + ((-t38 * t329 + t344 * t37) * t323 + (-t37 * rSges(6,2) + t344 * t38) * t321) * qJD(2) - (-qJD(2) * t489 - t37 + t590 + t598) * t38 + (-g(2) + t3) * (t418 + t488) + t628 * (t311 + (-t322 * t626 - t324 * t569) * t321 + t470)) * m(6) + (t48 * t439 + t49 * (t291 + t472) + (t321 * t553 + t348 * t49) * qJD(3) + ((-t48 * rSges(5,3) + t49 * (-t318 - t313)) * t321 + (t48 * (-t318 - t592) - t49 * t329) * t323) * qJD(2) - (-qJD(2) * t162 + t353 - t48 + t598) * t49 + (-g(2) + t22) * (t164 + t418) + (-g(1) + t21) * (-t162 + t470)) * m(5) + (t94 * (t289 + t469) + (t276 * t556 - t559) * qJD(3) + ((-pkin(2) - t277) * t555 + (t93 * (-rSges(4,3) - pkin(6)) + t94 * t430) * t321) * qJD(2) - (-qJD(2) * t171 - t218 - t436 - t93) * t94 + (-g(2) + t42) * t135 + (t41 - g(1)) * (t430 * t321 + t315 + t465)) * m(4) + (m(3) * (t231 ^ 2 + t235 ^ 2) + Icges(3,3) - t650) * qJDD(2) + (t608 + t611) * t573 + (-t610 + t612) * t572 + (t617 + t618) * t426 + (((t624 * t323 - t606 + t613) * t323 + (t624 * t321 - t138 + t409 + t614 - t667) * t321) * qJD(3) + t621 - t629) * t427 + (t616 - t619 + t620) * t425; t635 * t573 + t636 * t572 + (t617 * qJD(2) + t604 * qJD(3) + t611 * qJDD(2) + t613 * t215 + t614 * t216) * t571 + (t616 * qJD(2) + t605 * qJD(3) + t612 * qJDD(2) + t609 * t215 + t615 * t216) * t570 - (((t321 * t484 - t323 * t485) * t331 + (t321 * t486 + t323 * t487) * t330 + t623 * t324 + t625 * t322) * qJD(3) + (t644 * t322 + t645 * t324 + t330 * t467 + t331 * t466) * qJD(2)) * qJD(2) / 0.2e1 + (t619 * t323 + t618 * t321 + (-t610 * t321 + t608 * t323) * qJD(2)) * qJD(2) / 0.2e1 + (t608 * t321 + t610 * t323) * qJDD(2) / 0.2e1 + t621 * t457 / 0.2e1 + t620 * t456 / 0.2e1 + ((-t602 * t455 + t603) * t321 + (((t322 * t642 + t324 * t640 - t576) * t323 + (-t322 * t643 + t324 * t641 + t349 + t601) * t321) * qJD(3) + t607) * t323) * t427 + ((t614 * t321 + t613 * t323) * qJD(2) + t604) * t426 + ((t615 * t321 + t609 * t323) * qJD(2) + t605) * t425 + ((-t601 * t454 - t603) * t323 + ((t349 * t321 + t625 * t324 - t623 * t322 + (-t576 + t602) * t323) * qJD(3) + t607) * t321) * t424 + (-(t201 * t93 - t559) * qJD(2) - (t74 * (-t201 * t321 - t202 * t323) + t394 * t277) * qJD(3) + g(1) * t202 + g(2) * t201 - g(3) * t277 + t34 * t378 + t74 * ((t171 * t323 - t172 * t321) * qJD(2) + t387) + t394 * t248 + (-t42 * t321 - t41 * t323 + (-t323 * t94 + t556) * qJD(2)) * t276) * m(4) + (-(t368 + t558) * qJD(5) - (-t37 * t481 + t38 * (-t446 + t480)) * qJD(2) - (-t25 * t415 + (t25 * t480 - t37 * t588) * t323 + (t25 * t481 - t38 * t588) * t321) * qJD(3) + t1 * t502 + t25 * t444 + (t2 * t414 + t37 * t371 + t1 * t488 + t25 * t563 + (t25 * t489 + t38 * t414) * qJD(2)) * t323 + (t3 * t414 + t38 * t371 + t1 * t489 + t25 * t562 + (t25 * t441 + t37 * t474) * qJD(2)) * t321 - g(2) * t595 - g(3) * t588 - t369 * t567 - (t354 + t594) * g(1)) * m(6) + (-(t48 * t188 + t49 * (-t192 - t446)) * qJD(2) - (-t43 * t415 + (-t43 * t192 - t48 * t597) * t323 + (-t43 * t188 - t49 * t597) * t321) * qJD(3) + t4 * t502 + t43 * t444 + (t21 * t356 + t48 * t412 + t4 * t164 + t43 * t96 + (t43 * t162 + t356 * t49) * qJD(2)) * t323 + (t22 * t356 + t49 * t412 + t4 * t162 + t43 * t98 + (t43 * t498 + t553) * qJD(2)) * t321 - g(1) * t348 - g(3) * t597 - t356 * t567) * m(5); t574 * (g(1) * t321 - g(2) * t323) + 0.2e1 * (t2 * t571 + t3 * t570) * m(6) + 0.2e1 * (t21 * t571 + t22 * t570) * m(5); (-(t464 * t558 + t368) * qJD(3) + (qJD(3) * t401 + g(3) - t1) * t324 + (qJD(3) * t25 + t3 * t321 + t628 * t323 - t567) * t322) * m(6);];
tau = t5;
