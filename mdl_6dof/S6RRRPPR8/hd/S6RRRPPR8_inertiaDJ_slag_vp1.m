% Calculate time derivative of joint inertia matrix for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR8_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:39
% EndTime: 2019-03-09 16:04:38
% DurationCPUTime: 35.32s
% Computational Cost: add. (64250->1470), mult. (179222->1940), div. (0->0), fcn. (197846->10), ass. (0->562)
t532 = cos(pkin(6));
t534 = sin(qJ(2));
t531 = sin(pkin(6));
t678 = cos(qJ(3));
t612 = t531 * t678;
t677 = sin(qJ(3));
t508 = t532 * t677 + t534 * t612;
t679 = cos(qJ(2));
t599 = qJD(2) * t679;
t580 = t531 * t599;
t473 = qJD(3) * t508 + t580 * t677;
t597 = qJD(3) * t677;
t579 = t531 * t597;
t598 = qJD(3) * t678;
t474 = t532 * t598 - t534 * t579 + t580 * t678;
t632 = qJD(2) * t534;
t602 = t531 * t632;
t372 = Icges(6,5) * t473 - Icges(6,6) * t474 - Icges(6,3) * t602;
t375 = Icges(6,4) * t473 - Icges(6,2) * t474 - Icges(6,6) * t602;
t378 = Icges(6,1) * t473 - Icges(6,4) * t474 - Icges(6,5) * t602;
t611 = t531 * t677;
t507 = -t532 * t678 + t534 * t611;
t613 = t531 * t679;
t418 = Icges(6,5) * t507 - Icges(6,6) * t508 + Icges(6,3) * t613;
t421 = Icges(6,4) * t507 - Icges(6,2) * t508 + Icges(6,6) * t613;
t424 = Icges(6,1) * t507 - Icges(6,4) * t508 + Icges(6,5) * t613;
t134 = t372 * t613 - t508 * t375 + t507 * t378 - t418 * t602 - t474 * t421 + t473 * t424;
t373 = Icges(5,5) * t474 + Icges(5,6) * t602 + Icges(5,3) * t473;
t376 = Icges(5,4) * t474 + Icges(5,2) * t602 + Icges(5,6) * t473;
t379 = Icges(5,1) * t474 + Icges(5,4) * t602 + Icges(5,5) * t473;
t419 = Icges(5,5) * t508 - Icges(5,6) * t613 + Icges(5,3) * t507;
t422 = Icges(5,4) * t508 - Icges(5,2) * t613 + Icges(5,6) * t507;
t425 = Icges(5,1) * t508 - Icges(5,4) * t613 + Icges(5,5) * t507;
t135 = t507 * t373 - t376 * t613 + t508 * t379 + t473 * t419 + t422 * t602 + t474 * t425;
t374 = Icges(4,5) * t474 - Icges(4,6) * t473 + Icges(4,3) * t602;
t377 = Icges(4,4) * t474 - Icges(4,2) * t473 + Icges(4,6) * t602;
t380 = Icges(4,1) * t474 - Icges(4,4) * t473 + Icges(4,5) * t602;
t420 = Icges(4,5) * t508 - Icges(4,6) * t507 - Icges(4,3) * t613;
t423 = Icges(4,4) * t508 - Icges(4,2) * t507 - Icges(4,6) * t613;
t426 = Icges(4,1) * t508 - Icges(4,4) * t507 - Icges(4,5) * t613;
t136 = -t374 * t613 - t507 * t377 + t508 * t380 + t420 * t602 - t473 * t423 + t474 * t426;
t701 = t134 + t135 + t136;
t535 = sin(qJ(1));
t610 = t535 * t679;
t537 = cos(qJ(1));
t666 = t537 * t534;
t555 = -t532 * t666 - t610;
t583 = t537 * t612;
t477 = -t555 * t677 + t583;
t608 = t537 * t677;
t478 = -t531 * t608 - t555 * t678;
t609 = t537 * t679;
t582 = t532 * t609;
t667 = t535 * t534;
t554 = t582 - t667;
t333 = Icges(6,5) * t477 - Icges(6,6) * t478 + Icges(6,3) * t554;
t339 = Icges(6,4) * t477 - Icges(6,2) * t478 + Icges(6,6) * t554;
t345 = Icges(6,1) * t477 - Icges(6,4) * t478 + Icges(6,5) * t554;
t178 = t333 * t613 - t508 * t339 + t507 * t345;
t335 = Icges(5,5) * t478 - Icges(5,6) * t554 + Icges(5,3) * t477;
t341 = Icges(5,4) * t478 - Icges(5,2) * t554 + Icges(5,6) * t477;
t347 = Icges(5,1) * t478 - Icges(5,4) * t554 + Icges(5,5) * t477;
t180 = t507 * t335 - t341 * t613 + t508 * t347;
t337 = Icges(4,5) * t478 - Icges(4,6) * t477 - Icges(4,3) * t554;
t343 = Icges(4,4) * t478 - Icges(4,2) * t477 - Icges(4,6) * t554;
t349 = Icges(4,1) * t478 - Icges(4,4) * t477 - Icges(4,5) * t554;
t182 = -t337 * t613 - t507 * t343 + t508 * t349;
t700 = t178 + t180 + t182;
t512 = -t532 * t667 + t609;
t479 = t512 * t677 - t535 * t612;
t584 = t535 * t611;
t480 = t512 * t678 + t584;
t556 = -t532 * t610 - t666;
t334 = Icges(6,5) * t479 - Icges(6,6) * t480 + Icges(6,3) * t556;
t340 = Icges(6,4) * t479 - Icges(6,2) * t480 + Icges(6,6) * t556;
t346 = Icges(6,1) * t479 - Icges(6,4) * t480 + Icges(6,5) * t556;
t179 = t334 * t613 - t508 * t340 + t507 * t346;
t336 = Icges(5,5) * t480 - Icges(5,6) * t556 + Icges(5,3) * t479;
t342 = Icges(5,4) * t480 - Icges(5,2) * t556 + Icges(5,6) * t479;
t348 = Icges(5,1) * t480 - Icges(5,4) * t556 + Icges(5,5) * t479;
t181 = t507 * t336 - t342 * t613 + t508 * t348;
t338 = Icges(4,5) * t480 - Icges(4,6) * t479 - Icges(4,3) * t556;
t344 = Icges(4,4) * t480 - Icges(4,2) * t479 - Icges(4,6) * t556;
t350 = Icges(4,1) * t480 - Icges(4,4) * t479 - Icges(4,5) * t556;
t183 = -t338 * t613 - t507 * t344 + t508 * t350;
t699 = t179 + t181 + t183;
t698 = t701 * t532;
t533 = sin(qJ(6));
t536 = cos(qJ(6));
t411 = -t477 * t533 + t536 * t554;
t412 = t477 * t536 + t533 * t554;
t276 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t478;
t278 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t478;
t280 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t478;
t117 = t276 * t478 + t278 * t411 + t280 * t412;
t413 = -t479 * t533 + t536 * t556;
t414 = t479 * t536 + t533 * t556;
t277 = Icges(7,5) * t414 + Icges(7,6) * t413 + Icges(7,3) * t480;
t279 = Icges(7,4) * t414 + Icges(7,2) * t413 + Icges(7,6) * t480;
t281 = Icges(7,1) * t414 + Icges(7,4) * t413 + Icges(7,5) * t480;
t118 = t277 * t478 + t279 * t411 + t281 * t412;
t475 = -t507 * t533 + t536 * t613;
t557 = -t507 * t536 - t533 * t613;
t328 = -Icges(7,5) * t557 + Icges(7,6) * t475 + Icges(7,3) * t508;
t329 = -Icges(7,4) * t557 + Icges(7,2) * t475 + Icges(7,6) * t508;
t330 = -Icges(7,1) * t557 + Icges(7,4) * t475 + Icges(7,5) * t508;
t151 = t328 * t478 + t329 * t411 + t330 * t412;
t455 = qJD(1) * t512 + qJD(2) * t554;
t581 = qJD(1) * t612;
t325 = t455 * t677 - t535 * t581 - t537 * t579 - t555 * t598;
t454 = -qJD(1) * t556 - qJD(2) * t555;
t250 = -qJD(6) * t412 - t325 * t533 - t454 * t536;
t251 = qJD(6) * t411 + t325 * t536 - t454 * t533;
t326 = qJD(1) * t584 - qJD(3) * t583 + t455 * t678 + t555 * t597;
t138 = Icges(7,5) * t251 + Icges(7,6) * t250 + Icges(7,3) * t326;
t140 = Icges(7,4) * t251 + Icges(7,2) * t250 + Icges(7,6) * t326;
t142 = Icges(7,1) * t251 + Icges(7,4) * t250 + Icges(7,5) * t326;
t32 = t138 * t478 + t140 * t411 + t142 * t412 + t250 * t278 + t251 * t280 + t276 * t326;
t453 = qJD(1) * t555 + qJD(2) * t556;
t324 = t453 * t678 - t512 * t597 + (qJD(1) * t608 + t535 * t598) * t531;
t323 = qJD(3) * t480 + t453 * t677 - t537 * t581;
t452 = -qJD(1) * t582 - t537 * t599 + (qJD(2) * t532 + qJD(1)) * t667;
t248 = -qJD(6) * t414 - t323 * t533 + t452 * t536;
t249 = qJD(6) * t413 + t323 * t536 + t452 * t533;
t137 = Icges(7,5) * t249 + Icges(7,6) * t248 + Icges(7,3) * t324;
t139 = Icges(7,4) * t249 + Icges(7,2) * t248 + Icges(7,6) * t324;
t141 = Icges(7,1) * t249 + Icges(7,4) * t248 + Icges(7,5) * t324;
t33 = t137 * t478 + t139 * t411 + t141 * t412 + t250 * t279 + t251 * t281 + t277 * t326;
t361 = qJD(6) * t557 - t473 * t533 - t536 * t602;
t362 = qJD(6) * t475 + t473 * t536 - t533 * t602;
t239 = Icges(7,5) * t362 + Icges(7,6) * t361 + Icges(7,3) * t474;
t240 = Icges(7,4) * t362 + Icges(7,2) * t361 + Icges(7,6) * t474;
t241 = Icges(7,1) * t362 + Icges(7,4) * t361 + Icges(7,5) * t474;
t52 = t239 * t478 + t240 * t411 + t241 * t412 + t250 * t329 + t251 * t330 + t326 * t328;
t2 = t117 * t326 + t118 * t324 + t151 * t474 + t32 * t478 + t33 * t480 + t508 * t52;
t697 = -t2 / 0.2e1;
t143 = t249 * rSges(7,1) + t248 * rSges(7,2) + t324 * rSges(7,3);
t665 = t323 * pkin(5) + pkin(10) * t324 + t143;
t571 = -t251 * rSges(7,1) - t250 * rSges(7,2);
t144 = t326 * rSges(7,3) - t571;
t673 = t325 * pkin(5);
t664 = t326 * pkin(10) + t144 + t673;
t570 = -rSges(7,1) * t412 - rSges(7,2) * t411;
t282 = rSges(7,3) * t478 - t570;
t676 = pkin(5) * t477;
t658 = pkin(10) * t478 + t282 + t676;
t283 = t414 * rSges(7,1) + t413 * rSges(7,2) + t480 * rSges(7,3);
t657 = t479 * pkin(5) + pkin(10) * t480 + t283;
t530 = t537 * pkin(1);
t669 = t531 * t535;
t636 = pkin(8) * t669 + t530;
t204 = Icges(6,5) * t325 - Icges(6,6) * t326 - Icges(6,3) * t454;
t210 = Icges(6,4) * t325 - Icges(6,2) * t326 - Icges(6,6) * t454;
t216 = Icges(6,1) * t325 - Icges(6,4) * t326 - Icges(6,5) * t454;
t79 = -t508 * t210 + t507 * t216 - t474 * t339 + t473 * t345 + (t204 * t679 - t333 * t632) * t531;
t203 = Icges(6,5) * t323 - Icges(6,6) * t324 + Icges(6,3) * t452;
t209 = Icges(6,4) * t323 - Icges(6,2) * t324 + Icges(6,6) * t452;
t215 = Icges(6,1) * t323 - Icges(6,4) * t324 + Icges(6,5) * t452;
t80 = -t508 * t209 + t507 * t215 - t474 * t340 + t473 * t346 + (t203 * t679 - t334 * t632) * t531;
t206 = Icges(5,5) * t326 + Icges(5,6) * t454 + Icges(5,3) * t325;
t212 = Icges(5,4) * t326 + Icges(5,2) * t454 + Icges(5,6) * t325;
t218 = Icges(5,1) * t326 + Icges(5,4) * t454 + Icges(5,5) * t325;
t81 = t507 * t206 + t508 * t218 + t473 * t335 + t474 * t347 + (-t212 * t679 + t341 * t632) * t531;
t205 = Icges(5,5) * t324 - Icges(5,6) * t452 + Icges(5,3) * t323;
t211 = Icges(5,4) * t324 - Icges(5,2) * t452 + Icges(5,6) * t323;
t217 = Icges(5,1) * t324 - Icges(5,4) * t452 + Icges(5,5) * t323;
t82 = t507 * t205 + t508 * t217 + t473 * t336 + t474 * t348 + (-t211 * t679 + t342 * t632) * t531;
t208 = Icges(4,5) * t326 - Icges(4,6) * t325 + Icges(4,3) * t454;
t214 = Icges(4,4) * t326 - Icges(4,2) * t325 + Icges(4,6) * t454;
t220 = Icges(4,1) * t326 - Icges(4,4) * t325 + Icges(4,5) * t454;
t83 = -t507 * t214 + t508 * t220 - t473 * t343 + t474 * t349 + (-t208 * t679 + t337 * t632) * t531;
t207 = Icges(4,5) * t324 - Icges(4,6) * t323 - Icges(4,3) * t452;
t213 = Icges(4,4) * t324 - Icges(4,2) * t323 - Icges(4,6) * t452;
t219 = Icges(4,1) * t324 - Icges(4,4) * t323 - Icges(4,5) * t452;
t84 = -t507 * t213 + t508 * t219 - t473 * t344 + t474 * t350 + (-t207 * t679 + t338 * t632) * t531;
t122 = t276 * t508 + t278 * t475 - t280 * t557;
t123 = t277 * t508 + t279 * t475 - t281 * t557;
t36 = t138 * t508 + t140 * t475 - t142 * t557 + t276 * t474 + t278 * t361 + t280 * t362;
t37 = t137 * t508 + t139 * t475 - t141 * t557 + t277 * t474 + t279 * t361 + t281 * t362;
t72 = t508 * t239 + t475 * t240 - t241 * t557 + t474 * t328 + t361 * t329 + t362 * t330;
t71 = t72 * t532;
t9 = t71 + (-t36 * t537 + t37 * t535 + (t122 * t535 + t123 * t537) * qJD(1)) * t531;
t696 = t9 + t698 + ((-t79 - t81 - t83) * t537 + (t80 + t82 + t84) * t535 + (t535 * t700 + t699 * t537) * qJD(1)) * t531;
t695 = -t72 - t701;
t694 = 2 * m(3);
t693 = 2 * m(4);
t692 = 2 * m(5);
t691 = 2 * m(6);
t690 = 2 * m(7);
t689 = t531 ^ 2;
t688 = -pkin(3) - pkin(4);
t687 = t324 / 0.2e1;
t686 = t326 / 0.2e1;
t685 = t474 / 0.2e1;
t684 = t478 / 0.2e1;
t683 = t480 / 0.2e1;
t682 = t508 / 0.2e1;
t681 = t535 / 0.2e1;
t680 = -rSges(5,1) - pkin(3);
t675 = pkin(9) * t452;
t674 = pkin(9) * t556;
t672 = pkin(9) - qJ(5);
t161 = t328 * t508 + t329 * t475 - t330 * t557;
t671 = t161 * t474 + t72 * t508;
t670 = Icges(3,4) * t534;
t668 = t531 * t537;
t656 = -t325 * qJ(4) - t477 * qJD(4);
t233 = t326 * pkin(3) - t656;
t462 = t477 * qJ(4);
t394 = pkin(3) * t478 + t462;
t663 = -t233 * t556 - t452 * t394;
t232 = t324 * pkin(3) + qJ(4) * t323 + qJD(4) * t479;
t450 = t453 * pkin(2);
t359 = t450 - t675;
t332 = t532 * t359;
t662 = t532 * t232 + t332;
t594 = t324 * pkin(4) + qJD(5) * t556;
t274 = qJ(5) * t452 + t594;
t661 = -t232 - t274;
t360 = t455 * pkin(2) + t454 * pkin(9);
t660 = -t233 - t360;
t242 = rSges(7,1) * t362 + rSges(7,2) * t361 + rSges(7,3) * t474;
t659 = pkin(5) * t473 + pkin(10) * t474 + t242;
t327 = pkin(3) * t474 + qJ(4) * t473 + qJD(4) * t507;
t382 = rSges(5,1) * t474 + rSges(5,2) * t602 + rSges(5,3) * t473;
t655 = -t327 - t382;
t417 = t474 * pkin(4) + (-qJ(5) * t632 + qJD(5) * t679) * t531;
t654 = -t327 - t417;
t331 = -rSges(7,1) * t557 + rSges(7,2) * t475 + rSges(7,3) * t508;
t653 = pkin(5) * t507 + pkin(10) * t508 + t331;
t568 = rSges(5,2) * t554 - rSges(5,3) * t477;
t352 = rSges(5,1) * t478 - t568;
t652 = -t352 - t394;
t355 = t480 * rSges(5,1) - rSges(5,2) * t556 + t479 * rSges(5,3);
t395 = t480 * pkin(3) + qJ(4) * t479;
t651 = -t355 - t395;
t356 = t480 * rSges(4,1) - t479 * rSges(4,2) - rSges(4,3) * t556;
t506 = t512 * pkin(2);
t459 = t506 - t674;
t650 = -t356 - t459;
t358 = t556 * t394;
t499 = t554 * qJ(5);
t415 = pkin(4) * t478 + t499;
t649 = -t415 * t556 - t358;
t383 = rSges(4,1) * t474 - rSges(4,2) * t473 + rSges(4,3) * t602;
t633 = qJD(2) * t531;
t497 = (pkin(2) * t679 + pkin(9) * t534) * t633;
t648 = -t383 - t497;
t384 = t395 * t602;
t471 = t480 * pkin(4);
t416 = qJ(5) * t556 + t471;
t647 = t416 * t602 + t384;
t456 = pkin(3) * t508 + qJ(4) * t507;
t646 = t394 * t613 - t456 * t554;
t443 = t532 * t459;
t645 = t532 * t395 + t443;
t644 = -t394 - t415;
t643 = -t395 - t416;
t513 = (pkin(2) * t534 - pkin(9) * t679) * t531;
t635 = qJD(1) * t535;
t604 = t531 * t635;
t485 = t513 * t604;
t642 = t456 * t604 + t485;
t429 = t508 * rSges(5,1) - rSges(5,2) * t613 + t507 * rSges(5,3);
t641 = t429 + t456;
t430 = t508 * rSges(4,1) - t507 * rSges(4,2) - rSges(4,3) * t613;
t640 = -t430 - t513;
t458 = -pkin(2) * t555 - pkin(9) * t554;
t639 = t458 * t669 + t459 * t668;
t638 = t454 * qJ(5) - qJD(5) * t554;
t482 = t508 * pkin(4) + qJ(5) * t613;
t637 = t456 + t482;
t634 = qJD(1) * t537;
t631 = rSges(6,2) + t688;
t630 = m(6) / 0.2e1 + m(7) / 0.2e1;
t629 = -t679 / 0.2e1;
t51 = t239 * t480 + t240 * t413 + t241 * t414 + t248 * t329 + t249 * t330 + t324 * t328;
t628 = t51 / 0.2e1 + t37 / 0.2e1;
t627 = t52 / 0.2e1 + t36 / 0.2e1;
t626 = t233 * t613 - t327 * t554 + t454 * t456;
t625 = t532 * t274 + t662;
t275 = t326 * pkin(4) - t638;
t624 = -t275 + t660;
t221 = t323 * rSges(6,1) - t324 * rSges(6,2) + t452 * rSges(6,3);
t222 = t324 * rSges(5,1) - t452 * rSges(5,2) + t323 * rSges(5,3);
t223 = t324 * rSges(4,1) - t323 * rSges(4,2) - t452 * rSges(4,3);
t603 = t531 * t634;
t623 = t359 * t668 + t360 * t669 + t458 * t603;
t381 = rSges(6,1) * t473 - rSges(6,2) * t474 - rSges(6,3) * t602;
t622 = -t381 + t654;
t621 = -t497 + t655;
t572 = -rSges(6,1) * t477 - rSges(6,3) * t554;
t351 = -rSges(6,2) * t478 - t572;
t620 = -t351 + t644;
t354 = t479 * rSges(6,1) - t480 * rSges(6,2) + rSges(6,3) * t556;
t619 = -t354 + t643;
t618 = t532 * t416 + t645;
t617 = t482 * t604 + t642;
t428 = t507 * rSges(6,1) - t508 * rSges(6,2) + rSges(6,3) * t613;
t616 = t428 + t637;
t615 = -t513 - t641;
t307 = t453 * rSges(3,1) + t452 * rSges(3,2) + rSges(3,3) * t603;
t438 = t512 * rSges(3,1) + rSges(3,2) * t556 + rSges(3,3) * t669;
t614 = -rSges(7,3) - pkin(10) + t688;
t607 = t679 * Icges(3,4);
t606 = t679 * t232;
t605 = t679 * t395;
t601 = t151 / 0.2e1 + t122 / 0.2e1;
t152 = t328 * t480 + t329 * t413 + t330 * t414;
t600 = t152 / 0.2e1 + t123 / 0.2e1;
t596 = -t535 * pkin(1) + pkin(8) * t668;
t595 = t640 * t537;
t593 = -t275 * t556 - t452 * t415 + t663;
t592 = t654 - t659;
t591 = t644 - t658;
t590 = t643 - t657;
t589 = -t497 + t622;
t588 = t637 + t653;
t587 = t394 * t669 + t395 * t668 + t639;
t586 = t415 * t613 - t482 * t554 + t646;
t585 = -t513 - t616;
t578 = -pkin(1) * t635 + pkin(8) * t603;
t577 = t615 * t537;
t576 = -t178 / 0.2e1 - t180 / 0.2e1 - t182 / 0.2e1;
t575 = t179 / 0.2e1 + t181 / 0.2e1 + t183 / 0.2e1;
t574 = -t455 * rSges(3,1) + t454 * rSges(3,2);
t573 = -t325 * rSges(6,1) + t454 * rSges(6,3);
t569 = -t454 * rSges(5,2) - t325 * rSges(5,3);
t567 = -t497 + t592;
t566 = -t513 - t588;
t565 = t450 + t578;
t564 = t585 * t537;
t563 = t232 * t668 + t233 * t669 + t394 * t603 + t623;
t562 = t275 * t613 - t417 * t554 + t454 * t482 + t626;
t561 = t415 * t669 + t416 * t668 + t587;
t560 = t566 * t537;
t559 = -t458 + t596;
t558 = t395 + t506 + t636;
t226 = t326 * rSges(4,1) - t325 * rSges(4,2) + t454 * rSges(4,3);
t353 = rSges(4,1) * t478 - rSges(4,2) * t477 - rSges(4,3) * t554;
t553 = -t274 * t679 - t606;
t552 = -t416 * t679 - t605;
t551 = -t462 + t559;
t437 = -rSges(3,1) * t555 + rSges(3,2) * t554 - rSges(3,3) * t668;
t550 = -t499 + t551;
t488 = Icges(3,6) * t532 + (Icges(3,2) * t679 + t670) * t531;
t489 = Icges(3,5) * t532 + (Icges(3,1) * t534 + t607) * t531;
t493 = (Icges(3,5) * t679 - Icges(3,6) * t534) * t633;
t494 = (-Icges(3,2) * t534 + t607) * t633;
t495 = (Icges(3,1) * t679 - t670) * t633;
t549 = t531 * t534 * t495 - t488 * t602 + t489 * t580 + t532 * t493 + t494 * t613;
t548 = t274 * t668 + t275 * t669 + t415 * t603 + t563;
t547 = -qJD(1) * t636 - t360;
t546 = -t556 * t672 + t471 + t558;
t545 = t232 + t565;
t110 = t325 * t424 - t326 * t421 + t372 * t554 - t375 * t478 + t378 * t477 - t418 * t454;
t111 = t325 * t419 + t326 * t425 + t373 * t477 - t376 * t554 + t379 * t478 + t422 * t454;
t112 = -t325 * t423 + t326 * t426 - t374 * t554 - t377 * t477 + t380 * t478 + t420 * t454;
t544 = t83 / 0.2e1 + t81 / 0.2e1 + t79 / 0.2e1 + t110 / 0.2e1 + t111 / 0.2e1 + t112 / 0.2e1 + t627;
t107 = t323 * t424 - t324 * t421 + t372 * t556 - t375 * t480 + t378 * t479 + t418 * t452;
t108 = t323 * t419 + t324 * t425 + t373 * t479 - t376 * t556 + t379 * t480 - t422 * t452;
t109 = -t323 * t423 + t324 * t426 - t374 * t556 - t377 * t479 + t380 * t480 - t420 * t452;
t543 = t84 / 0.2e1 + t82 / 0.2e1 + t80 / 0.2e1 + t107 / 0.2e1 + t108 / 0.2e1 + t109 / 0.2e1 + t628;
t197 = t418 * t554 - t421 * t478 + t424 * t477;
t198 = t419 * t477 - t422 * t554 + t425 * t478;
t199 = -t420 * t554 - t423 * t477 + t426 * t478;
t542 = t198 / 0.2e1 + t199 / 0.2e1 + t197 / 0.2e1 - t576 + t601;
t200 = t418 * t556 - t421 * t480 + t424 * t479;
t201 = t419 * t479 - t422 * t556 + t425 * t480;
t202 = -t420 * t556 - t423 * t479 + t426 * t480;
t541 = -t200 / 0.2e1 - t201 / 0.2e1 - t202 / 0.2e1 - t575 - t600;
t540 = t547 + t656;
t539 = t540 + t638;
t538 = -t452 * t672 + t545 + t594;
t496 = (rSges(3,1) * t679 - rSges(3,2) * t534) * t633;
t492 = t532 * rSges(3,3) + (rSges(3,1) * t534 + rSges(3,2) * t679) * t531;
t487 = Icges(3,3) * t532 + (Icges(3,5) * t534 + Icges(3,6) * t679) * t531;
t436 = Icges(3,1) * t512 + Icges(3,4) * t556 + Icges(3,5) * t669;
t435 = -Icges(3,1) * t555 + Icges(3,4) * t554 - Icges(3,5) * t668;
t434 = Icges(3,4) * t512 + Icges(3,2) * t556 + Icges(3,6) * t669;
t433 = -Icges(3,4) * t555 + Icges(3,2) * t554 - Icges(3,6) * t668;
t432 = Icges(3,5) * t512 + Icges(3,6) * t556 + Icges(3,3) * t669;
t431 = -Icges(3,5) * t555 + Icges(3,6) * t554 - Icges(3,3) * t668;
t404 = t438 + t636;
t403 = -t437 + t596;
t370 = -t532 * t437 - t492 * t668;
t369 = t438 * t532 - t492 * t669;
t308 = rSges(3,3) * t604 - t574;
t306 = Icges(3,1) * t455 - Icges(3,4) * t454 + Icges(3,5) * t604;
t305 = Icges(3,1) * t453 + Icges(3,4) * t452 + Icges(3,5) * t603;
t304 = Icges(3,4) * t455 - Icges(3,2) * t454 + Icges(3,6) * t604;
t303 = Icges(3,4) * t453 + Icges(3,2) * t452 + Icges(3,6) * t603;
t302 = Icges(3,5) * t455 - Icges(3,6) * t454 + Icges(3,3) * t604;
t301 = Icges(3,5) * t453 + Icges(3,6) * t452 + Icges(3,3) * t603;
t292 = (-t530 + (-rSges(3,3) - pkin(8)) * t669) * qJD(1) + t574;
t291 = t578 + t307;
t290 = t487 * t669 + t488 * t556 + t489 * t512;
t289 = -t487 * t668 + t488 * t554 - t489 * t555;
t286 = t636 - t650;
t285 = -t353 + t559;
t284 = t549 * t532;
t273 = t532 * t307 + (-t492 * t634 - t496 * t535) * t531;
t272 = -t532 * t308 + (t492 * t635 - t496 * t537) * t531;
t267 = -t356 * t613 + t430 * t556;
t266 = t353 * t613 - t430 * t554;
t265 = t532 * t432 + (t434 * t679 + t436 * t534) * t531;
t264 = t532 * t431 + (t433 * t679 + t435 * t534) * t531;
t258 = t432 * t669 + t434 * t556 + t436 * t512;
t257 = t431 * t669 + t433 * t556 + t435 * t512;
t256 = -t432 * t668 + t434 * t554 - t436 * t555;
t255 = -t431 * t668 + t433 * t554 - t435 * t555;
t254 = -t420 * t613 - t507 * t423 + t508 * t426;
t253 = t507 * t419 - t422 * t613 + t508 * t425;
t252 = t418 * t613 - t508 * t421 + t507 * t424;
t245 = t254 * t602;
t244 = t253 * t602;
t243 = t252 * t602;
t238 = t355 + t558 - t674;
t237 = t478 * t680 + t551 + t568;
t236 = -t353 * t556 + t356 * t554;
t235 = (-t353 - t458) * t532 + t531 * t595;
t234 = t356 * t532 + t640 * t669 + t443;
t225 = t326 * rSges(5,1) - t569;
t224 = -t326 * rSges(6,2) - t573;
t193 = t546 + t354;
t192 = t478 * t631 + t550 + t572;
t191 = -t454 * t488 + t455 * t489 + t554 * t494 - t555 * t495 + (t487 * t635 - t493 * t537) * t531;
t190 = t452 * t488 + t453 * t489 + t556 * t494 + t512 * t495 + (t487 * t634 + t493 * t535) * t531;
t189 = (t353 * t535 + t356 * t537) * t531 + t639;
t188 = t283 * t508 - t331 * t480;
t187 = -t282 * t508 + t331 * t478;
t186 = (-t355 * t679 - t605) * t531 + t641 * t556;
t185 = t352 * t613 - t429 * t554 + t646;
t177 = -t226 + t547;
t176 = t223 + t565 - t675;
t175 = (-t458 + t652) * t532 + t531 * t577;
t174 = t355 * t532 + t615 * t669 + t645;
t173 = -t338 * t556 - t344 * t479 + t350 * t480;
t172 = -t337 * t556 - t343 * t479 + t349 * t480;
t171 = t336 * t479 - t342 * t556 + t348 * t480;
t170 = t335 * t479 - t341 * t556 + t347 * t480;
t169 = t334 * t556 - t340 * t480 + t346 * t479;
t168 = t333 * t556 - t339 * t480 + t345 * t479;
t167 = -t338 * t554 - t344 * t477 + t350 * t478;
t166 = -t337 * t554 - t343 * t477 + t349 * t478;
t165 = t336 * t477 - t342 * t554 + t348 * t478;
t164 = t335 * t477 - t341 * t554 + t347 * t478;
t163 = t334 * t554 - t340 * t478 + t346 * t477;
t162 = t333 * t554 - t339 * t478 + t345 * t477;
t160 = t161 * t602;
t159 = t282 * t480 - t283 * t478;
t158 = -t352 * t556 - t554 * t651 - t358;
t156 = t546 + t657;
t155 = t478 * t614 + t550 + t570 - t676;
t154 = (-t354 * t679 + t552) * t531 + t616 * t556;
t153 = t351 * t613 - t428 * t554 + t586;
t150 = t532 * t301 + (t679 * t303 + t305 * t534 + (-t434 * t534 + t436 * t679) * qJD(2)) * t531;
t149 = t532 * t302 + (t679 * t304 + t306 * t534 + (-t433 * t534 + t435 * t679) * qJD(2)) * t531;
t148 = (-t458 + t620) * t532 + t531 * t564;
t147 = t354 * t532 + t585 * t669 + t618;
t145 = (t352 * t535 + t355 * t537) * t531 + t587;
t130 = t532 * t223 + t332 + (qJD(1) * t595 + t535 * t648) * t531;
t129 = t485 + (-t226 - t360) * t532 + (t430 * t635 + t537 * t648) * t531;
t128 = -t554 * t383 + t454 * t430 + (t226 * t679 - t353 * t632) * t531;
t127 = t556 * t383 + t452 * t430 + (-t223 * t679 + t356 * t632) * t531;
t126 = -t351 * t556 - t554 * t619 + t649;
t125 = t326 * t680 + t540 + t569;
t124 = t222 + t545 - t675;
t121 = (t351 * t535 + t354 * t537) * t531 + t561;
t120 = t277 * t480 + t279 * t413 + t281 * t414;
t119 = t276 * t480 + t278 * t413 + t280 * t414;
t116 = (-t657 * t679 + t552) * t531 + t588 * t556;
t115 = -t554 * t653 + t613 * t658 + t586;
t114 = (-t458 + t591) * t532 + t531 * t560;
t113 = t532 * t657 + t566 * t669 + t618;
t106 = t326 * t631 + t539 + t573;
t105 = t538 + t221;
t104 = t223 * t554 - t226 * t556 - t353 * t452 - t356 * t454;
t103 = -t554 * t590 - t556 * t658 + t649;
t102 = (t535 * t658 + t537 * t657) * t531 + t561;
t101 = t532 * t222 + (qJD(1) * t577 + t535 * t621) * t531 + t662;
t100 = (-t225 + t660) * t532 + (t429 * t635 + t537 * t621) * t531 + t642;
t99 = (t223 * t537 + t226 * t535 + (t353 * t537 + t535 * t650) * qJD(1)) * t531 + t623;
t98 = t202 * t532 + (-t172 * t537 + t173 * t535) * t531;
t97 = t201 * t532 + (-t170 * t537 + t171 * t535) * t531;
t96 = t200 * t532 + (-t168 * t537 + t169 * t535) * t531;
t95 = t199 * t532 + (-t166 * t537 + t167 * t535) * t531;
t94 = t198 * t532 + (-t164 * t537 + t165 * t535) * t531;
t93 = t197 * t532 + (-t162 * t537 + t163 * t535) * t531;
t92 = -t172 * t554 - t173 * t556 - t202 * t613;
t91 = -t170 * t554 - t171 * t556 - t201 * t613;
t90 = -t168 * t554 - t169 * t556 - t200 * t613;
t89 = -t166 * t554 - t167 * t556 - t199 * t613;
t88 = -t164 * t554 - t165 * t556 - t198 * t613;
t87 = -t162 * t554 - t163 * t556 - t197 * t613;
t86 = -t554 * t382 + t454 * t429 + (t225 * t679 + t632 * t652) * t531 + t626;
t85 = t384 - t655 * t556 + t641 * t452 + (-t222 * t679 + t355 * t632 - t606) * t531;
t78 = -t144 * t508 + t242 * t478 - t282 * t474 + t326 * t331;
t77 = t143 * t508 - t242 * t480 + t283 * t474 - t324 * t331;
t76 = t532 * t221 + (qJD(1) * t564 + t535 * t589) * t531 + t625;
t75 = (-t224 + t624) * t532 + (t428 * t635 + t537 * t589) * t531 + t617;
t74 = t326 * t614 + t539 + t571 - t673;
t73 = t538 + t665;
t69 = -t207 * t554 - t213 * t477 + t219 * t478 - t325 * t344 + t326 * t350 + t338 * t454;
t68 = -t208 * t554 - t214 * t477 + t220 * t478 - t325 * t343 + t326 * t349 + t337 * t454;
t67 = t205 * t477 - t211 * t554 + t217 * t478 + t325 * t336 + t326 * t348 + t342 * t454;
t66 = t206 * t477 - t212 * t554 + t218 * t478 + t325 * t335 + t326 * t347 + t341 * t454;
t65 = t203 * t554 - t209 * t478 + t215 * t477 + t325 * t346 - t326 * t340 - t334 * t454;
t64 = t204 * t554 - t210 * t478 + t216 * t477 + t325 * t345 - t326 * t339 - t333 * t454;
t63 = -t207 * t556 - t213 * t479 + t219 * t480 - t323 * t344 + t324 * t350 - t338 * t452;
t62 = -t208 * t556 - t214 * t479 + t220 * t480 - t323 * t343 + t324 * t349 - t337 * t452;
t61 = t205 * t479 - t211 * t556 + t217 * t480 + t323 * t336 + t324 * t348 - t342 * t452;
t60 = t206 * t479 - t212 * t556 + t218 * t480 + t323 * t335 + t324 * t347 - t341 * t452;
t59 = t203 * t556 - t209 * t480 + t215 * t479 + t323 * t346 - t324 * t340 + t334 * t452;
t58 = t204 * t556 - t210 * t480 + t216 * t479 + t323 * t345 - t324 * t339 + t333 * t452;
t57 = -t554 * t381 + t454 * t428 + (t224 * t679 + t620 * t632) * t531 + t562;
t56 = -t622 * t556 + t616 * t452 + (-t221 * t679 + t354 * t632 + t553) * t531 + t647;
t55 = -t225 * t556 - t352 * t452 - (-t222 - t232) * t554 + t651 * t454 + t663;
t54 = t161 * t532 + (-t122 * t537 + t123 * t535) * t531;
t53 = -t122 * t554 - t123 * t556 - t161 * t613;
t50 = t122 * t478 + t123 * t480 + t161 * t508;
t49 = -t143 * t478 + t144 * t480 + t282 * t324 - t283 * t326;
t48 = (t222 * t537 + t225 * t535 + (t352 * t537 + (-t459 + t651) * t535) * qJD(1)) * t531 + t563;
t47 = t152 * t532 + (-t119 * t537 + t120 * t535) * t531;
t46 = t151 * t532 + (-t117 * t537 + t118 * t535) * t531;
t45 = -t119 * t554 - t120 * t556 - t152 * t613;
t44 = -t117 * t554 - t118 * t556 - t151 * t613;
t43 = t119 * t478 + t120 * t480 + t152 * t508;
t42 = t117 * t478 + t118 * t480 + t151 * t508;
t41 = t665 * t532 + (qJD(1) * t560 + t535 * t567) * t531 + t625;
t40 = (t624 - t664) * t532 + (t537 * t567 + t635 * t653) * t531 + t617;
t39 = -t224 * t556 - t351 * t452 - (-t221 + t661) * t554 + t619 * t454 + t593;
t38 = (t221 * t537 + t224 * t535 + (t351 * t537 + (-t459 + t619) * t535) * qJD(1)) * t531 + t548;
t35 = -t659 * t554 + t653 * t454 + (t591 * t632 + t664 * t679) * t531 + t562;
t34 = -t592 * t556 + t588 * t452 + (t632 * t657 - t665 * t679 + t553) * t531 + t647;
t31 = t137 * t480 + t139 * t413 + t141 * t414 + t248 * t279 + t249 * t281 + t277 * t324;
t30 = t138 * t480 + t140 * t413 + t142 * t414 + t248 * t278 + t249 * t280 + t276 * t324;
t29 = (t665 * t537 + t664 * t535 + (t658 * t537 + (-t459 + t590) * t535) * qJD(1)) * t531 + t548;
t28 = -t664 * t556 - t658 * t452 - (t661 - t665) * t554 + t590 * t454 + t593;
t24 = -t136 * t613 + t182 * t454 - t183 * t452 - t554 * t83 - t556 * t84 + t245;
t23 = -t135 * t613 + t180 * t454 - t181 * t452 - t554 * t81 - t556 * t82 + t244;
t22 = -t134 * t613 + t178 * t454 - t179 * t452 - t554 * t79 - t556 * t80 + t243;
t21 = t112 * t532 + (t535 * t69 - t537 * t68 + (t166 * t535 + t167 * t537) * qJD(1)) * t531;
t20 = t111 * t532 + (t535 * t67 - t537 * t66 + (t164 * t535 + t165 * t537) * qJD(1)) * t531;
t19 = t110 * t532 + (t535 * t65 - t537 * t64 + (t162 * t535 + t163 * t537) * qJD(1)) * t531;
t18 = t109 * t532 + (t535 * t63 - t537 * t62 + (t172 * t535 + t173 * t537) * qJD(1)) * t531;
t17 = t108 * t532 + (t535 * t61 - t537 * t60 + (t170 * t535 + t171 * t537) * qJD(1)) * t531;
t16 = t107 * t532 + (t535 * t59 - t537 * t58 + (t168 * t535 + t169 * t537) * qJD(1)) * t531;
t15 = t166 * t454 - t167 * t452 - t68 * t554 - t69 * t556 + (-t112 * t679 + t199 * t632) * t531;
t14 = t164 * t454 - t165 * t452 - t66 * t554 - t67 * t556 + (-t111 * t679 + t198 * t632) * t531;
t13 = t162 * t454 - t163 * t452 - t64 * t554 - t65 * t556 + (-t110 * t679 + t197 * t632) * t531;
t12 = t172 * t454 - t173 * t452 - t62 * t554 - t63 * t556 + (-t109 * t679 + t202 * t632) * t531;
t11 = t170 * t454 - t171 * t452 - t60 * t554 - t61 * t556 + (-t108 * t679 + t201 * t632) * t531;
t10 = t168 * t454 - t169 * t452 - t58 * t554 - t59 * t556 + (-t107 * t679 + t200 * t632) * t531;
t8 = t122 * t454 - t123 * t452 - t36 * t554 - t37 * t556 - t613 * t72 + t160;
t7 = t122 * t326 + t123 * t324 + t36 * t478 + t37 * t480 + t671;
t6 = t52 * t532 + (-t32 * t537 + t33 * t535 + (t117 * t535 + t118 * t537) * qJD(1)) * t531;
t5 = t51 * t532 + (-t30 * t537 + t31 * t535 + (t119 * t535 + t120 * t537) * qJD(1)) * t531;
t4 = t117 * t454 - t118 * t452 - t32 * t554 - t33 * t556 + (t151 * t632 - t52 * t679) * t531;
t3 = t119 * t454 - t120 * t452 - t30 * t554 - t31 * t556 + (t152 * t632 - t51 * t679) * t531;
t1 = t119 * t326 + t120 * t324 + t152 * t474 + t30 * t478 + t31 * t480 + t508 * t51;
t25 = [(t176 * t286 + t177 * t285) * t693 + (t124 * t238 + t125 * t237) * t692 + (t105 * t193 + t106 * t192) * t691 + (t155 * t74 + t156 * t73) * t690 + (t291 * t404 + t292 * t403) * t694 + t549 - t695; t71 + t284 + m(4) * (t129 * t285 + t130 * t286 + t176 * t234 + t177 * t235) + m(5) * (t100 * t237 + t101 * t238 + t124 * t174 + t125 * t175) + m(6) * (t105 * t147 + t106 * t148 + t192 * t75 + t193 * t76) + m(7) * (t113 * t73 + t114 * t74 + t155 * t40 + t156 * t41) + m(3) * (t272 * t403 + t273 * t404 + t291 * t369 + t292 * t370) + ((-t149 / 0.2e1 - t191 / 0.2e1 - t544) * t537 + (t150 / 0.2e1 + t190 / 0.2e1 + t543) * t535 + ((t265 / 0.2e1 + t290 / 0.2e1 - t541) * t537 + (t264 / 0.2e1 + t289 / 0.2e1 + t542) * t535) * qJD(1)) * t531 + t698; (t102 * t29 + t113 * t41 + t114 * t40) * t690 + (t121 * t38 + t147 * t76 + t148 * t75) * t691 + (t100 * t175 + t101 * t174 + t145 * t48) * t692 + (t129 * t235 + t130 * t234 + t189 * t99) * t693 + (t370 * t272 + t369 * t273 + (t437 * t535 + t438 * t537) * (t307 * t537 + t308 * t535 + (t437 * t537 - t438 * t535) * qJD(1)) * t689) * t694 + (((t303 * t556 + t512 * t305 + t452 * t434 + t453 * t436) * t535 + t258 * t634 - (t304 * t556 + t512 * t306 + t452 * t433 + t453 * t435) * t537 + t257 * t635 + ((t301 * t535 + t432 * t634) * t535 - (t302 * t535 + t431 * t634) * t537) * t531) * t531 + t18 + t17 + t16 + t5) * t669 + (-t21 - t20 - t19 - t6 - ((t303 * t554 - t305 * t555 - t454 * t434 + t455 * t436) * t535 + t256 * t634 - (t304 * t554 - t306 * t555 - t454 * t433 + t455 * t435) * t537 + t255 * t635 + ((-t301 * t537 + t432 * t635) * t535 - (-t302 * t537 + t431 * t635) * t537) * t531) * t531) * t668 + ((-t255 * t537 + t256 * t535) * t531 + t95 + t94 + t93 + t46) * t604 + ((-t257 * t537 + t258 * t535) * t531 + t96 + t47 + t98 + t97) * t603 + (-t191 * t668 + t289 * t604 + t284 + (-t149 * t537 + t150 * t535 + (t264 * t535 + t265 * t537) * qJD(1)) * t531 + t290 * t603 + t190 * t669 + t696) * t532; t160 + t243 + t244 + t245 + m(4) * (t127 * t286 + t128 * t285 + t176 * t267 + t177 * t266) + m(5) * (t124 * t186 + t125 * t185 + t237 * t86 + t238 * t85) + m(6) * (t105 * t154 + t106 * t153 + t192 * t57 + t193 * t56) + m(7) * (t115 * t74 + t116 * t73 + t155 * t35 + t156 * t34) - t543 * t556 - t544 * t554 + t542 * t454 + t541 * t452 + t695 * t613; (t24 / 0.2e1 + t23 / 0.2e1 + t22 / 0.2e1 + t8 / 0.2e1) * t532 - (t18 / 0.2e1 + t17 / 0.2e1 + t16 / 0.2e1 + t5 / 0.2e1) * t556 - (t21 / 0.2e1 + t20 / 0.2e1 + t19 / 0.2e1 + t6 / 0.2e1) * t554 + (t95 / 0.2e1 + t93 / 0.2e1 + t94 / 0.2e1 + t46 / 0.2e1) * t454 + (-t98 / 0.2e1 - t96 / 0.2e1 - t97 / 0.2e1 - t47 / 0.2e1) * t452 + m(4) * (t104 * t189 + t127 * t234 + t128 * t235 + t129 * t266 + t130 * t267 + t236 * t99) + m(5) * (t100 * t185 + t101 * t186 + t145 * t55 + t158 * t48 + t174 * t85 + t175 * t86) + m(6) * (t121 * t39 + t126 * t38 + t147 * t56 + t148 * t57 + t153 * t75 + t154 * t76) + m(7) * (t102 * t28 + t103 * t29 + t113 * t34 + t114 * t35 + t115 * t40 + t116 * t41) + ((-t13 / 0.2e1 - t14 / 0.2e1 - t15 / 0.2e1 - t4 / 0.2e1) * t537 + (t10 / 0.2e1 + t11 / 0.2e1 + t12 / 0.2e1 + t3 / 0.2e1) * t535 + (t54 / 0.2e1 + t576 * t668 + t575 * t669 + (t252 / 0.2e1 + t253 / 0.2e1 + t254 / 0.2e1) * t532) * t632 + ((t90 / 0.2e1 + t91 / 0.2e1 + t92 / 0.2e1 + t45 / 0.2e1) * t537 + (t87 / 0.2e1 + t88 / 0.2e1 + t89 / 0.2e1 + t44 / 0.2e1) * t535) * qJD(1) + t696 * t629) * t531; (-t22 - t23 - t24 - t8) * t613 - (t10 + t11 + t12 + t3) * t556 - (t13 + t14 + t15 + t4) * t554 + (t89 + t87 + t88 + t44) * t454 + (-t92 - t90 - t91 - t45) * t452 + (t53 - t699 * t556 - t700 * t554 + (-t252 - t253 - t254) * t613) * t602 + (t104 * t236 + t127 * t267 + t128 * t266) * t693 + (t158 * t55 + t185 * t86 + t186 * t85) * t692 + (t126 * t39 + t153 * t57 + t154 * t56) * t691 + (t103 * t28 + t115 * t35 + t116 * t34) * t690; m(7) * (t155 * t323 + t156 * t325 + t477 * t73 + t479 * t74) + m(5) * (t124 * t477 + t125 * t479 + t237 * t323 + t238 * t325) + m(6) * (t105 * t477 + t106 * t479 + t192 * t323 + t193 * t325); m(7) * (t102 * t473 + t113 * t325 + t114 * t323 + t29 * t507 + t40 * t479 + t41 * t477) + m(6) * (t121 * t473 + t147 * t325 + t148 * t323 + t38 * t507 + t477 * t76 + t479 * t75) + m(5) * (t100 * t479 + t101 * t477 + t145 * t473 + t174 * t325 + t175 * t323 + t48 * t507); m(7) * (t103 * t473 + t115 * t323 + t116 * t325 + t28 * t507 + t34 * t477 + t35 * t479) + m(6) * (t126 * t473 + t153 * t323 + t154 * t325 + t39 * t507 + t477 * t56 + t479 * t57) + m(5) * (t158 * t473 + t185 * t323 + t186 * t325 + t477 * t85 + t479 * t86 + t507 * t55); 0.4e1 * (m(5) / 0.2e1 + t630) * (t323 * t479 + t325 * t477 + t473 * t507); m(7) * (t155 * t452 - t156 * t454 + t554 * t73 + t556 * t74) + m(6) * (t105 * t554 + t106 * t556 + t192 * t452 - t193 * t454); m(7) * (-t454 * t113 + t452 * t114 + t556 * t40 + t554 * t41 + (-t102 * t632 + t29 * t679) * t531) + m(6) * (-t454 * t147 + t452 * t148 + t554 * t76 + t556 * t75 + (-t121 * t632 + t38 * t679) * t531); m(7) * (t452 * t115 - t454 * t116 + t554 * t34 + t556 * t35 + (-t103 * t632 + t28 * t679) * t531) + m(6) * (t452 * t153 - t454 * t154 + t554 * t56 + t556 * t57 + (-t126 * t632 + t39 * t679) * t531); 0.2e1 * t630 * (t556 * t323 + t554 * t325 + t452 * t479 - t454 * t477 + (t473 * t679 - t507 * t632) * t531); 0.4e1 * t630 * (-t534 * t599 * t689 + t452 * t556 - t454 * t554); m(7) * (t155 * t78 + t156 * t77 + t187 * t74 + t188 * t73) + t628 * t480 + t627 * t478 + t601 * t326 + t600 * t324 + t671; m(7) * (t102 * t49 + t113 * t77 + t114 * t78 + t159 * t29 + t187 * t40 + t188 * t41) + t47 * t687 + t5 * t683 + t532 * t7 / 0.2e1 + t46 * t686 + t6 * t684 + t54 * t685 + t9 * t682 + (t537 * t697 + t1 * t681 + (t42 * t681 + t537 * t43 / 0.2e1) * qJD(1)) * t531; m(7) * (t103 * t49 + t115 * t78 + t116 * t77 + t159 * t28 + t187 * t35 + t188 * t34) + t45 * t687 + t3 * t683 + t454 * t42 / 0.2e1 + t554 * t697 + t44 * t686 + t4 * t684 - t452 * t43 / 0.2e1 - t556 * t1 / 0.2e1 + t53 * t685 + t8 * t682 + (t50 * t632 / 0.2e1 + t7 * t629) * t531; m(7) * (t159 * t473 + t187 * t323 + t188 * t325 + t477 * t77 + t479 * t78 + t49 * t507); m(7) * (t187 * t452 - t188 * t454 + t77 * t554 + t78 * t556 + (-t159 * t632 + t49 * t679) * t531); t324 * t43 + t480 * t1 + t326 * t42 + t478 * t2 + t474 * t50 + t508 * t7 + (t159 * t49 + t187 * t78 + t188 * t77) * t690;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
