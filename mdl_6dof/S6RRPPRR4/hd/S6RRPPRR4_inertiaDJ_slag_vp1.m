% Calculate time derivative of joint inertia matrix for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:16
% EndTime: 2019-03-09 09:02:03
% DurationCPUTime: 29.05s
% Computational Cost: add. (82209->1267), mult. (228214->1684), div. (0->0), fcn. (269514->12), ass. (0->517)
t489 = sin(qJ(1));
t491 = cos(qJ(1));
t486 = cos(pkin(6));
t616 = sin(pkin(11));
t617 = cos(pkin(11));
t633 = sin(qJ(2));
t635 = cos(qJ(2));
t503 = t633 * t616 - t635 * t617;
t501 = t486 * t503;
t504 = t616 * t635 + t617 * t633;
t413 = -t489 * t504 - t491 * t501;
t500 = t486 * t504;
t414 = -t489 * t503 + t491 * t500;
t485 = sin(pkin(6));
t613 = t485 * t491;
t285 = -Icges(5,5) * t613 - Icges(5,6) * t414 - Icges(5,3) * t413;
t292 = Icges(4,4) * t414 + Icges(4,2) * t413 - Icges(4,6) * t613;
t666 = t292 - t285;
t415 = t489 * t501 - t491 * t504;
t416 = -t489 * t500 - t491 * t503;
t614 = t485 * t489;
t284 = Icges(5,5) * t614 - Icges(5,6) * t416 - Icges(5,3) * t415;
t293 = Icges(4,4) * t416 + Icges(4,2) * t415 + Icges(4,6) * t614;
t665 = -t293 + t284;
t287 = -Icges(5,4) * t613 - Icges(5,2) * t414 - Icges(5,6) * t413;
t294 = Icges(4,1) * t414 + Icges(4,4) * t413 - Icges(4,5) * t613;
t664 = t294 - t287;
t286 = Icges(5,4) * t614 - Icges(5,2) * t416 - Icges(5,6) * t415;
t295 = Icges(4,1) * t416 + Icges(4,4) * t415 + Icges(4,5) * t614;
t663 = -t295 + t286;
t288 = Icges(5,1) * t614 - Icges(5,4) * t416 - Icges(5,5) * t415;
t291 = Icges(4,5) * t416 + Icges(4,6) * t415 + Icges(4,3) * t614;
t569 = t491 * t633;
t572 = t489 * t635;
t462 = -t486 * t572 - t569;
t570 = t491 * t635;
t571 = t489 * t633;
t509 = t486 * t571 - t570;
t393 = -Icges(3,5) * t509 + Icges(3,6) * t462 + Icges(3,3) * t614;
t654 = t288 + t291 + t393;
t289 = -Icges(5,1) * t613 - Icges(5,4) * t414 - Icges(5,5) * t413;
t290 = Icges(4,5) * t414 + Icges(4,6) * t413 - Icges(4,3) * t613;
t507 = -t486 * t569 - t572;
t508 = -t486 * t570 + t571;
t392 = -Icges(3,5) * t507 - Icges(3,6) * t508 - Icges(3,3) * t613;
t653 = t289 + t290 + t392;
t452 = t504 * t485;
t439 = qJD(2) * t452;
t451 = t503 * t485;
t440 = t451 * qJD(2);
t373 = Icges(5,6) * t440 + Icges(5,3) * t439;
t374 = Icges(5,2) * t440 + Icges(5,6) * t439;
t375 = -Icges(4,5) * t440 - Icges(4,6) * t439;
t376 = Icges(5,4) * t440 + Icges(5,5) * t439;
t377 = -Icges(4,4) * t440 - Icges(4,2) * t439;
t378 = -Icges(4,1) * t440 - Icges(4,4) * t439;
t383 = Icges(5,5) * t486 - Icges(5,6) * t452 + Icges(5,3) * t451;
t385 = Icges(5,4) * t486 - Icges(5,2) * t452 + Icges(5,6) * t451;
t386 = Icges(4,4) * t452 - Icges(4,2) * t451 + Icges(4,6) * t486;
t388 = Icges(4,1) * t452 - Icges(4,4) * t451 + Icges(4,5) * t486;
t567 = t633 * Icges(3,4);
t445 = Icges(3,6) * t486 + (Icges(3,2) * t635 + t567) * t485;
t568 = t635 * Icges(3,4);
t446 = Icges(3,5) * t486 + (Icges(3,1) * t633 + t568) * t485;
t593 = qJD(2) * t485;
t453 = (Icges(3,5) * t635 - Icges(3,6) * t633) * t593;
t454 = (-Icges(3,2) * t633 + t568) * t593;
t455 = (Icges(3,1) * t635 - t567) * t593;
t564 = qJD(2) * t635;
t541 = t485 * t564;
t563 = qJD(2) * t633;
t573 = t485 * t633;
t674 = t446 * t541 + t455 * t573 + (-t445 * t563 + t454 * t635) * t485 + (-t374 + t378) * t452 + (t373 - t377) * t451 + (t385 - t388) * t440 + (t383 - t386) * t439 + (t453 + t376 + t375) * t486;
t488 = sin(qJ(5));
t634 = cos(qJ(5));
t363 = t415 * t634 + t488 * t614;
t574 = t485 * t634;
t364 = -t415 * t488 + t489 * t574;
t238 = Icges(6,5) * t364 - Icges(6,6) * t363 + Icges(6,3) * t416;
t240 = Icges(6,4) * t364 - Icges(6,2) * t363 + Icges(6,6) * t416;
t242 = Icges(6,1) * t364 - Icges(6,4) * t363 + Icges(6,5) * t416;
t365 = -t413 * t634 + t488 * t613;
t543 = t491 * t574;
t366 = -t413 * t488 - t543;
t107 = t238 * t414 + t240 * t365 + t242 * t366;
t239 = Icges(6,5) * t366 + Icges(6,6) * t365 + Icges(6,3) * t414;
t241 = Icges(6,4) * t366 + Icges(6,2) * t365 + Icges(6,6) * t414;
t243 = Icges(6,1) * t366 + Icges(6,4) * t365 + Icges(6,5) * t414;
t108 = t239 * t414 + t241 * t365 + t243 * t366;
t424 = -t451 * t634 + t486 * t488;
t425 = t451 * t488 + t486 * t634;
t326 = Icges(6,5) * t425 - Icges(6,6) * t424 + Icges(6,3) * t452;
t327 = Icges(6,4) * t425 - Icges(6,2) * t424 + Icges(6,6) * t452;
t328 = Icges(6,1) * t425 - Icges(6,4) * t424 + Icges(6,5) * t452;
t138 = t326 * t414 + t327 * t365 + t328 * t366;
t487 = sin(qJ(6));
t490 = cos(qJ(6));
t352 = -t425 * t487 + t452 * t490;
t353 = t425 * t490 + t452 * t487;
t246 = Icges(7,5) * t353 + Icges(7,6) * t352 + Icges(7,3) * t424;
t247 = Icges(7,4) * t353 + Icges(7,2) * t352 + Icges(7,6) * t424;
t248 = Icges(7,1) * t353 + Icges(7,4) * t352 + Icges(7,5) * t424;
t278 = -t366 * t487 + t414 * t490;
t279 = t366 * t490 + t414 * t487;
t101 = -t246 * t365 + t247 * t278 + t248 * t279;
t499 = qJD(2) * t504;
t441 = t486 * t499;
t594 = qJD(1) * t491;
t660 = -qJD(1) * t501 - qJD(2) * t503;
t316 = -t491 * t441 - t489 * t660 - t504 * t594;
t542 = qJD(1) * t574;
t235 = qJD(5) * t365 - t316 * t488 + t489 * t542;
t442 = qJD(2) * t501;
t661 = qJD(1) * t500 + t499;
t317 = -t491 * t442 - t489 * t661 - t503 * t594;
t158 = -qJD(6) * t279 - t235 * t487 + t317 * t490;
t159 = qJD(6) * t278 + t235 * t490 + t317 * t487;
t276 = -t364 * t487 + t416 * t490;
t277 = t364 * t490 + t416 * t487;
t181 = Icges(7,5) * t277 + Icges(7,6) * t276 + Icges(7,3) * t363;
t183 = Icges(7,4) * t277 + Icges(7,2) * t276 + Icges(7,6) * t363;
t185 = Icges(7,1) * t277 + Icges(7,4) * t276 + Icges(7,5) * t363;
t595 = qJD(1) * t489;
t566 = t485 * t595;
t234 = t316 * t634 - qJD(5) * t543 + (-qJD(5) * t413 + t566) * t488;
t314 = -t489 * t441 + t491 * t660 - t504 * t595;
t237 = -qJD(5) * t363 + t314 * t488 + t491 * t542;
t315 = -t489 * t442 + t491 * t661 - t503 * t595;
t160 = -qJD(6) * t277 - t237 * t487 - t315 * t490;
t161 = qJD(6) * t276 + t237 * t490 - t315 * t487;
t565 = t485 * t594;
t236 = qJD(5) * t364 - t314 * t634 + t488 * t565;
t81 = Icges(7,5) * t161 + Icges(7,6) * t160 + Icges(7,3) * t236;
t83 = Icges(7,4) * t161 + Icges(7,2) * t160 + Icges(7,6) * t236;
t85 = Icges(7,1) * t161 + Icges(7,4) * t160 + Icges(7,5) * t236;
t16 = t158 * t183 + t159 * t185 + t181 * t234 + t278 * t83 + t279 * t85 - t365 * t81;
t182 = Icges(7,5) * t279 + Icges(7,6) * t278 - Icges(7,3) * t365;
t184 = Icges(7,4) * t279 + Icges(7,2) * t278 - Icges(7,6) * t365;
t186 = Icges(7,1) * t279 + Icges(7,4) * t278 - Icges(7,5) * t365;
t80 = Icges(7,5) * t159 + Icges(7,6) * t158 + Icges(7,3) * t234;
t82 = Icges(7,4) * t159 + Icges(7,2) * t158 + Icges(7,6) * t234;
t84 = Icges(7,1) * t159 + Icges(7,4) * t158 + Icges(7,5) * t234;
t17 = t158 * t184 + t159 * t186 + t182 * t234 + t278 * t82 + t279 * t84 - t365 * t80;
t349 = -qJD(5) * t424 + t439 * t488;
t252 = -qJD(6) * t353 - t349 * t487 - t440 * t490;
t253 = qJD(6) * t352 + t349 * t490 - t440 * t487;
t350 = qJD(5) * t425 - t439 * t634;
t154 = Icges(7,5) * t253 + Icges(7,6) * t252 + Icges(7,3) * t350;
t155 = Icges(7,4) * t253 + Icges(7,2) * t252 + Icges(7,6) * t350;
t156 = Icges(7,1) * t253 + Icges(7,4) * t252 + Icges(7,5) * t350;
t31 = -t154 * t365 + t155 * t278 + t156 * t279 + t158 * t247 + t159 * t248 + t234 * t246;
t68 = -t181 * t365 + t183 * t278 + t185 * t279;
t69 = -t182 * t365 + t184 * t278 + t186 * t279;
t3 = -t101 * t440 + t16 * t416 + t17 * t414 + t31 * t452 - t315 * t68 + t317 * t69;
t130 = Icges(6,5) * t237 - Icges(6,6) * t236 - Icges(6,3) * t315;
t132 = Icges(6,4) * t237 - Icges(6,2) * t236 - Icges(6,6) * t315;
t134 = Icges(6,1) * t237 - Icges(6,4) * t236 - Icges(6,5) * t315;
t41 = t130 * t414 + t132 * t365 + t134 * t366 - t234 * t240 + t235 * t242 + t238 * t317;
t129 = Icges(6,5) * t235 - Icges(6,6) * t234 + Icges(6,3) * t317;
t131 = Icges(6,4) * t235 - Icges(6,2) * t234 + Icges(6,6) * t317;
t133 = Icges(6,1) * t235 - Icges(6,4) * t234 + Icges(6,5) * t317;
t42 = t129 * t414 + t131 * t365 + t133 * t366 - t234 * t241 + t235 * t243 + t239 * t317;
t254 = Icges(6,5) * t349 - Icges(6,6) * t350 - Icges(6,3) * t440;
t255 = Icges(6,4) * t349 - Icges(6,2) * t350 - Icges(6,6) * t440;
t256 = Icges(6,1) * t349 - Icges(6,4) * t350 - Icges(6,5) * t440;
t63 = -t234 * t327 + t235 * t328 + t254 * t414 + t255 * t365 + t256 * t366 + t317 * t326;
t673 = -t107 * t315 + t108 * t317 - t138 * t440 + t41 * t416 + t414 * t42 + t452 * t63 + t3;
t105 = t238 * t416 - t240 * t363 + t242 * t364;
t106 = t239 * t416 - t241 * t363 + t243 * t364;
t137 = t326 * t416 - t327 * t363 + t328 * t364;
t100 = t246 * t363 + t247 * t276 + t248 * t277;
t18 = t160 * t183 + t161 * t185 + t181 * t236 + t276 * t83 + t277 * t85 + t363 * t81;
t19 = t160 * t184 + t161 * t186 + t182 * t236 + t276 * t82 + t277 * t84 + t363 * t80;
t32 = t154 * t363 + t155 * t276 + t156 * t277 + t160 * t247 + t161 * t248 + t236 * t246;
t66 = t181 * t363 + t183 * t276 + t185 * t277;
t67 = t182 * t363 + t184 * t276 + t186 * t277;
t4 = -t100 * t440 + t18 * t416 + t19 * t414 - t315 * t66 + t317 * t67 + t32 * t452;
t43 = t130 * t416 - t132 * t363 + t134 * t364 - t236 * t240 + t237 * t242 - t238 * t315;
t44 = t129 * t416 - t131 * t363 + t133 * t364 - t236 * t241 + t237 * t243 - t239 * t315;
t64 = -t236 * t327 + t237 * t328 + t254 * t416 - t255 * t363 + t256 * t364 - t315 * t326;
t672 = -t105 * t315 + t106 * t317 - t137 * t440 + t414 * t44 + t416 * t43 + t452 * t64 + t4;
t484 = pkin(2) * t635 + pkin(1);
t625 = pkin(1) - t484;
t671 = t489 * t625;
t208 = Icges(5,5) * t565 + Icges(5,6) * t315 + Icges(5,3) * t314;
t215 = -Icges(4,4) * t315 - Icges(4,2) * t314 + Icges(4,6) * t565;
t670 = -t215 + t208;
t207 = Icges(5,5) * t566 - Icges(5,6) * t317 - Icges(5,3) * t316;
t216 = Icges(4,4) * t317 + Icges(4,2) * t316 + Icges(4,6) * t566;
t669 = t216 - t207;
t210 = Icges(5,4) * t565 + Icges(5,2) * t315 + Icges(5,6) * t314;
t217 = -Icges(4,1) * t315 - Icges(4,4) * t314 + Icges(4,5) * t565;
t668 = t217 - t210;
t209 = Icges(5,4) * t566 - Icges(5,2) * t317 - Icges(5,6) * t316;
t218 = Icges(4,1) * t317 + Icges(4,4) * t316 + Icges(4,5) * t566;
t667 = t218 - t209;
t482 = pkin(8) * t613;
t662 = -t489 * pkin(1) + t482;
t395 = -Icges(3,4) * t509 + Icges(3,2) * t462 + Icges(3,6) * t614;
t397 = -Icges(3,1) * t509 + Icges(3,4) * t462 + Icges(3,5) * t614;
t659 = t395 * t462 - t397 * t509 - t665 * t415 - t663 * t416 + t654 * t614;
t394 = -Icges(3,4) * t507 - Icges(3,2) * t508 - Icges(3,6) * t613;
t396 = -Icges(3,1) * t507 - Icges(3,4) * t508 - Icges(3,5) * t613;
t658 = t394 * t508 + t396 * t507 - t666 * t413 - t664 * t414 + t653 * t613;
t657 = t674 * t486;
t656 = t395 * t508 + t397 * t507 + t665 * t413 + t663 * t414 + t654 * t613;
t655 = t394 * t462 - t396 * t509 + t666 * t415 + t664 * t416 + t653 * t614;
t211 = Icges(5,1) * t566 - Icges(5,4) * t317 - Icges(5,5) * t316;
t212 = Icges(5,1) * t565 + Icges(5,4) * t315 + Icges(5,5) * t314;
t213 = -Icges(4,5) * t315 - Icges(4,6) * t314 + Icges(4,3) * t565;
t214 = Icges(4,5) * t317 + Icges(4,6) * t316 + Icges(4,3) * t566;
t419 = qJD(1) * t508 + qJD(2) * t509;
t420 = qJD(1) * t507 + qJD(2) * t462;
t318 = Icges(3,5) * t420 + Icges(3,6) * t419 + Icges(3,3) * t565;
t421 = qJD(1) * t462 + qJD(2) * t507;
t422 = -qJD(1) * t509 - qJD(2) * t508;
t319 = Icges(3,5) * t422 + Icges(3,6) * t421 + Icges(3,3) * t566;
t652 = (t212 + t318 + t213) * t489 + (-t319 - t214 - t211) * t491;
t651 = 0.2e1 * t486;
t650 = 0.2e1 * t489;
t649 = m(7) / 0.2e1;
t21 = t181 * t350 + t183 * t252 + t185 * t253 + t352 * t83 + t353 * t85 + t424 * t81;
t22 = t182 * t350 + t184 * t252 + t186 * t253 + t352 * t82 + t353 * t84 + t424 * t80;
t109 = t246 * t424 + t247 * t352 + t248 * t353;
t48 = t424 * t154 + t352 * t155 + t353 * t156 + t350 * t246 + t252 * t247 + t253 * t248;
t622 = t109 * t350 + t48 * t424;
t73 = t181 * t424 + t183 * t352 + t185 * t353;
t74 = t182 * t424 + t184 * t352 + t186 * t353;
t7 = t21 * t363 - t22 * t365 + t74 * t234 + t73 * t236 + t622;
t648 = t7 / 0.2e1;
t647 = t234 / 0.2e1;
t646 = t236 / 0.2e1;
t645 = t350 / 0.2e1;
t644 = t363 / 0.2e1;
t643 = -t365 / 0.2e1;
t642 = t424 / 0.2e1;
t641 = -t440 / 0.2e1;
t640 = t489 / 0.2e1;
t639 = -t491 / 0.2e1;
t638 = rSges(5,2) - pkin(3);
t637 = -rSges(6,3) - pkin(3);
t636 = rSges(7,3) + pkin(10);
t632 = pkin(1) * t491;
t631 = pkin(2) * t486;
t630 = t235 * pkin(5);
t629 = t317 * pkin(3);
t628 = t366 * pkin(5);
t627 = t414 * pkin(3);
t624 = t316 * rSges(5,3);
t623 = t413 * rSges(5,3);
t621 = -t109 * t440 + t48 * t452;
t151 = (t326 * t452 - t327 * t424 + t328 * t425) * t440;
t76 = t452 * t254 - t424 * t255 + t425 * t256 - t440 * t326 - t350 * t327 + t349 * t328;
t620 = t76 * t452 - t151;
t530 = -t159 * rSges(7,1) - t158 * rSges(7,2);
t86 = t234 * rSges(7,3) - t530;
t619 = t234 * pkin(10) + t630 + t86;
t87 = t161 * rSges(7,1) + t160 * rSges(7,2) + t236 * rSges(7,3);
t618 = t237 * pkin(5) + t236 * pkin(10) + t87;
t457 = t633 * t631 + (-pkin(8) - qJ(3)) * t485;
t615 = t457 * t491;
t612 = t489 * t484;
t472 = t491 * t484;
t157 = rSges(7,1) * t253 + rSges(7,2) * t252 + rSges(7,3) * t350;
t611 = pkin(5) * t349 + pkin(10) * t350 + t157;
t187 = t277 * rSges(7,1) + t276 * rSges(7,2) + t363 * rSges(7,3);
t610 = t364 * pkin(5) + pkin(10) * t363 + t187;
t529 = -t279 * rSges(7,1) - t278 * rSges(7,2);
t188 = -t365 * rSges(7,3) - t529;
t609 = -t365 * pkin(10) + t188 + t628;
t200 = -t315 * pkin(3) + t314 * qJ(4) - t415 * qJD(4);
t465 = -t485 * qJD(3) + t564 * t631;
t549 = pkin(2) * t563;
t511 = -t489 * t465 - t491 * t549;
t562 = -t485 * pkin(8) - t457;
t337 = (t491 * t562 + t671) * qJD(1) + t511;
t334 = t486 * t337;
t608 = t486 * t200 + t334;
t605 = t316 * qJ(4) + t413 * qJD(4);
t201 = -t605 + t629;
t536 = -t457 * t595 + t465 * t491 - t489 * t549;
t588 = pkin(8) * t614;
t338 = (-t491 * t625 - t588) * qJD(1) + t536;
t607 = -t201 - t338;
t249 = rSges(7,1) * t353 + rSges(7,2) * t352 + rSges(7,3) * t424;
t606 = pkin(5) * t425 + pkin(10) * t424 + t249;
t331 = t416 * pkin(3) - qJ(4) * t415;
t404 = t489 * t562 + t472 - t632;
t391 = t486 * t404;
t604 = t486 * t331 + t391;
t405 = t413 * qJ(4);
t330 = -t405 + t627;
t403 = t482 + t615 - t671;
t603 = -t330 - t403;
t602 = -t331 - t404;
t464 = pkin(2) * t541 + qJD(3) * t486;
t601 = pkin(3) * t440 - qJ(4) * t439 - qJD(4) * t451 - t464;
t401 = pkin(3) * t452 + qJ(4) * t451;
t468 = pkin(2) * t573 + t486 * qJ(3);
t447 = t468 * t566;
t600 = t401 * t566 + t447;
t599 = rSges(4,1) * t440 + rSges(4,2) * t439 - t464;
t598 = t403 * t614 + t404 * t613;
t390 = rSges(4,1) * t452 - rSges(4,2) * t451 + rSges(4,3) * t486;
t597 = -t390 - t468;
t596 = -t401 - t468;
t372 = -pkin(4) * t613 + t414 * pkin(9);
t590 = 0.2e1 * qJD(1);
t589 = pkin(9) * t440 * t485;
t481 = pkin(4) * t614;
t480 = rSges(5,1) * t614;
t478 = rSges(4,3) * t614;
t587 = t21 / 0.2e1 + t32 / 0.2e1;
t586 = -t22 / 0.2e1 - t31 / 0.2e1;
t585 = t73 / 0.2e1 + t100 / 0.2e1;
t584 = t74 / 0.2e1 + t101 / 0.2e1;
t312 = t317 * pkin(9);
t282 = pkin(4) * t566 + t312;
t583 = -t282 + t607;
t136 = t237 * rSges(6,1) - t236 * rSges(6,2) - t315 * rSges(6,3);
t257 = rSges(6,1) * t349 - rSges(6,2) * t350 - rSges(6,3) * t440;
t582 = -t257 + t601;
t221 = -t315 * rSges(4,1) - t314 * rSges(4,2) + rSges(4,3) * t565;
t371 = pkin(9) * t416 + t481;
t581 = t486 * t371 + t604;
t580 = -t372 + t603;
t579 = -t371 + t602;
t578 = t337 * t613 + t338 * t614 + t403 * t565;
t577 = -rSges(5,2) * t440 - rSges(5,3) * t439 + t601;
t244 = t364 * rSges(6,1) - t363 * rSges(6,2) + t416 * rSges(6,3);
t389 = rSges(5,1) * t486 - rSges(5,2) * t452 + rSges(5,3) * t451;
t576 = -t389 + t596;
t429 = pkin(4) * t486 + pkin(9) * t452;
t575 = -t429 + t596;
t300 = t416 * rSges(4,1) + t415 * rSges(4,2) + t478;
t324 = t420 * rSges(3,1) + t419 * rSges(3,2) + rSges(3,3) * t565;
t400 = -rSges(3,1) * t509 + t462 * rSges(3,2) + rSges(3,3) * t614;
t220 = rSges(5,1) * t565 + t315 * rSges(5,2) + t314 * rSges(5,3);
t297 = -t416 * rSges(5,2) - t415 * rSges(5,3) + t480;
t283 = pkin(4) * t565 - t315 * pkin(9);
t561 = 2 * m(3);
t559 = 2 * m(4);
t557 = 2 * m(5);
t555 = 2 * m(6);
t553 = 0.2e1 * m(7);
t552 = t597 * t491;
t551 = -t489 * t457 + t472;
t550 = m(5) / 0.2e1 + m(6) / 0.2e1 + t649;
t548 = t601 - t611;
t547 = t486 * t283 + t489 * t589 + t608;
t546 = t330 * t614 + t331 * t613 + t598;
t329 = rSges(6,1) * t425 - rSges(6,2) * t424 + rSges(6,3) * t452;
t545 = -t329 + t575;
t544 = t429 * t566 + t491 * t589 + t600;
t23 = -t315 * t609 - t317 * t610 - t414 * t618 + t416 * t619;
t532 = -t235 * rSges(6,1) + t234 * rSges(6,2);
t135 = t317 * rSges(6,3) - t532;
t531 = -t366 * rSges(6,1) - t365 * rSges(6,2);
t245 = t414 * rSges(6,3) - t531;
t60 = t135 * t416 - t136 * t414 - t244 * t317 - t245 * t315;
t540 = m(6) * t60 + m(7) * t23;
t539 = t491 * t576;
t535 = -t422 * rSges(3,1) - t421 * rSges(3,2);
t534 = -t317 * rSges(4,1) - t316 * rSges(4,2);
t533 = -t414 * rSges(4,1) - t413 * rSges(4,2);
t528 = t575 - t606;
t525 = -t612 - t615;
t524 = t545 * t491;
t51 = t130 * t452 - t132 * t424 + t134 * t425 - t238 * t440 - t240 * t350 + t242 * t349;
t521 = t64 / 0.2e1 + t51 / 0.2e1 + t587;
t52 = t129 * t452 - t131 * t424 + t133 * t425 - t239 * t440 - t241 * t350 + t243 * t349;
t520 = t63 / 0.2e1 + t52 / 0.2e1 - t586;
t519 = t200 * t613 + t201 * t614 + t330 * t565 + t578;
t518 = t371 * t613 + t372 * t614 + t546;
t114 = t238 * t452 - t240 * t424 + t242 * t425;
t517 = -t137 / 0.2e1 - t114 / 0.2e1 - t585;
t115 = t239 * t452 - t241 * t424 + t243 * t425;
t516 = t138 / 0.2e1 + t115 / 0.2e1 + t584;
t506 = t282 * t614 + t283 * t613 + t372 * t565 + t519;
t20 = (t618 * t491 + t619 * t489 + (t609 * t491 + (t579 - t610) * t489) * qJD(1)) * t485 + t506;
t40 = (t135 * t489 + t136 * t491 + (t245 * t491 + (-t244 + t579) * t489) * qJD(1)) * t485 + t506;
t219 = rSges(5,1) * t566 - t317 * rSges(5,2) - t624;
t298 = -rSges(5,1) * t613 - t414 * rSges(5,2) - t623;
t59 = (t219 * t489 + t220 * t491 + (t298 * t491 + (-t297 + t602) * t489) * qJD(1)) * t485 + t519;
t515 = m(5) * t59 + m(6) * t40 + m(7) * t20;
t514 = t528 * t491;
t513 = -t536 + t605;
t512 = t331 + t551;
t510 = t405 + t525 - t372;
t399 = -rSges(3,1) * t507 - rSges(3,2) * t508 - rSges(3,3) * t613;
t502 = t371 + t512;
t495 = qJD(1) * t525 + t511;
t494 = (-t481 - t472) * qJD(1) - t312 + t513;
t493 = t200 + t495;
t492 = t283 + t493;
t456 = (rSges(3,1) * t635 - rSges(3,2) * t633) * t593;
t448 = t486 * rSges(3,3) + (rSges(3,1) * t633 + rSges(3,2) * t635) * t485;
t444 = Icges(3,3) * t486 + (Icges(3,5) * t633 + Icges(3,6) * t635) * t485;
t387 = Icges(5,1) * t486 - Icges(5,4) * t452 + Icges(5,5) * t451;
t384 = Icges(4,5) * t452 - Icges(4,6) * t451 + Icges(4,3) * t486;
t361 = t400 + t588 + t632;
t360 = -t399 + t662;
t336 = -t486 * t399 - t448 * t613;
t335 = t400 * t486 - t448 * t614;
t325 = rSges(3,3) * t566 - t535;
t323 = Icges(3,1) * t422 + Icges(3,4) * t421 + Icges(3,5) * t566;
t322 = Icges(3,1) * t420 + Icges(3,4) * t419 + Icges(3,5) * t565;
t321 = Icges(3,4) * t422 + Icges(3,2) * t421 + Icges(3,6) * t566;
t320 = Icges(3,4) * t420 + Icges(3,2) * t419 + Icges(3,6) * t565;
t299 = -rSges(4,3) * t613 - t533;
t281 = (-t632 + (-rSges(3,3) - pkin(8)) * t614) * qJD(1) + t535;
t280 = qJD(1) * t662 + t324;
t275 = t444 * t614 + t445 * t462 - t446 * t509;
t274 = -t444 * t613 - t445 * t508 - t446 * t507;
t264 = t551 + t300;
t263 = -t612 + (rSges(4,3) * t485 - t457) * t491 + t533;
t261 = t486 * t324 + (-t448 * t594 - t456 * t489) * t485;
t260 = -t486 * t325 + (t448 * t595 - t456 * t491) * t485;
t251 = t486 * t393 + (t395 * t635 + t397 * t633) * t485;
t250 = t486 * t392 + (t394 * t635 + t396 * t633) * t485;
t222 = rSges(4,3) * t566 - t534;
t206 = t512 + t297;
t205 = t623 - t612 + t405 + (rSges(5,1) * t485 - t457) * t491 + t638 * t414;
t203 = (-t299 - t403) * t486 + t485 * t552;
t202 = t300 * t486 + t597 * t614 + t391;
t196 = t384 * t614 + t386 * t415 + t388 * t416;
t195 = -t384 * t613 + t413 * t386 + t414 * t388;
t194 = -t413 * t383 - t414 * t385 - t387 * t613;
t193 = -t383 * t415 - t385 * t416 + t387 * t614;
t192 = t421 * t445 + t422 * t446 - t508 * t454 - t507 * t455 + (t444 * t595 - t453 * t491) * t485;
t191 = t419 * t445 + t420 * t446 + t462 * t454 - t509 * t455 + (t444 * t594 + t453 * t489) * t485;
t190 = (-t478 - t472) * qJD(1) + t534 - t536;
t189 = t495 + t221;
t178 = t244 * t452 - t329 * t416;
t177 = -t245 * t452 + t329 * t414;
t174 = t502 + t244;
t173 = t414 * t637 + t510 + t531;
t170 = t291 * t486 - t293 * t451 + t295 * t452;
t169 = t290 * t486 - t292 * t451 + t294 * t452;
t168 = t285 * t451 - t287 * t452 + t289 * t486;
t167 = t284 * t451 - t286 * t452 + t288 * t486;
t163 = (-t298 + t603) * t486 + t485 * t539;
t162 = t297 * t486 + t576 * t614 + t604;
t142 = t486 * t318 + (t633 * t322 + t635 * t320 + (-t395 * t633 + t397 * t635) * qJD(2)) * t485;
t141 = t486 * t319 + (t633 * t323 + t635 * t321 + (-t394 * t633 + t396 * t635) * qJD(2)) * t485;
t140 = -t244 * t414 + t245 * t416;
t126 = t486 * t221 + t334 + (qJD(1) * t552 + t489 * t599) * t485;
t125 = t447 + (-t222 - t338) * t486 + (t390 * t595 + t491 * t599) * t485;
t124 = t624 + t638 * t317 + (-t480 - t472) * qJD(1) + t513;
t123 = t493 + t220;
t122 = (t297 * t491 + t298 * t489) * t485 + t546;
t121 = -t188 * t424 - t249 * t365;
t120 = t187 * t424 - t249 * t363;
t119 = t502 + t610;
t118 = t365 * t636 + t510 + t529 - t627 - t628;
t117 = (-t245 + t580) * t486 + t485 * t524;
t116 = t244 * t486 + t545 * t614 + t581;
t113 = t316 * t386 + t317 * t388 + t413 * t377 + t414 * t378 + (-t375 * t491 + t384 * t595) * t485;
t112 = -t314 * t386 - t315 * t388 + t415 * t377 + t416 * t378 + (t375 * t489 + t384 * t594) * t485;
t111 = t314 * t383 + t315 * t385 - t415 * t373 - t416 * t374 + (t376 * t489 + t387 * t594) * t485;
t110 = -t316 * t383 - t317 * t385 - t413 * t373 - t414 * t374 + (-t376 * t491 + t387 * t595) * t485;
t103 = t187 * t365 + t188 * t363;
t99 = -t416 * t606 + t452 * t610;
t98 = t414 * t606 - t452 * t609;
t97 = (t244 * t491 + t245 * t489) * t485 + t518;
t96 = t486 * t220 + (qJD(1) * t539 + t489 * t577) * t485 + t608;
t95 = (-t219 + t607) * t486 + (t389 * t595 + t491 * t577) * t485 + t600;
t94 = t317 * t637 + t494 + t532;
t93 = t492 + t136;
t92 = t213 * t486 - t215 * t451 + t217 * t452 - t293 * t439 - t295 * t440;
t91 = t214 * t486 - t216 * t451 + t218 * t452 - t292 * t439 - t294 * t440;
t90 = t208 * t451 - t210 * t452 + t212 * t486 + t284 * t439 + t286 * t440;
t89 = t207 * t451 - t209 * t452 + t211 * t486 + t285 * t439 + t287 * t440;
t88 = (t221 * t491 + t222 * t489 + (t299 * t491 + (-t300 - t404) * t489) * qJD(1)) * t485 + t578;
t79 = (t580 - t609) * t486 + t485 * t514;
t78 = t486 * t610 + t528 * t614 + t581;
t77 = -t414 * t610 + t416 * t609;
t75 = t76 * t486;
t71 = -t135 * t452 + t245 * t440 + t257 * t414 + t317 * t329;
t70 = t136 * t452 - t244 * t440 - t257 * t416 + t315 * t329;
t65 = (t489 * t609 + t491 * t610) * t485 + t518;
t62 = t486 * t136 + (qJD(1) * t524 + t489 * t582) * t485 + t547;
t61 = (-t135 + t583) * t486 + (t329 * t595 + t491 * t582) * t485 + t544;
t58 = -t234 * t636 + t494 + t530 - t629 - t630;
t57 = t492 + t618;
t56 = t138 * t486 + (t107 * t489 - t108 * t491) * t485;
t55 = t137 * t486 + (t105 * t489 - t106 * t491) * t485;
t54 = t107 * t416 + t108 * t414 + t138 * t452;
t53 = t105 * t416 + t106 * t414 + t137 * t452;
t50 = -t157 * t363 + t187 * t350 - t236 * t249 + t424 * t87;
t49 = -t157 * t365 - t188 * t350 + t234 * t249 - t424 * t86;
t47 = t48 * t486;
t39 = t109 * t486 + (t489 * t73 - t491 * t74) * t485;
t38 = t109 * t452 + t414 * t74 + t416 * t73;
t37 = t618 * t486 + (qJD(1) * t514 + t489 * t548) * t485 + t547;
t36 = (t583 - t619) * t486 + (t491 * t548 + t595 * t606) * t485 + t544;
t35 = t109 * t424 + t363 * t73 - t365 * t74;
t34 = t317 * t606 + t414 * t611 + t440 * t609 - t452 * t619;
t33 = t315 * t606 - t416 * t611 - t440 * t610 + t452 * t618;
t30 = -t187 * t234 + t188 * t236 + t363 * t86 + t365 * t87;
t29 = t101 * t486 + (t489 * t68 - t491 * t69) * t485;
t28 = t100 * t486 + (t489 * t66 - t491 * t67) * t485;
t27 = t101 * t452 + t414 * t69 + t416 * t68;
t26 = t100 * t452 + t414 * t67 + t416 * t66;
t25 = t101 * t424 + t363 * t68 - t365 * t69;
t24 = t100 * t424 + t363 * t66 - t365 * t67;
t15 = t75 + (t51 * t489 - t52 * t491 + (t114 * t491 + t115 * t489) * qJD(1)) * t485;
t14 = -t114 * t315 + t115 * t317 + t52 * t414 + t51 * t416 + t620;
t13 = t64 * t486 + (t43 * t489 - t44 * t491 + (t105 * t491 + t106 * t489) * qJD(1)) * t485;
t12 = t63 * t486 + (t41 * t489 - t42 * t491 + (t107 * t491 + t108 * t489) * qJD(1)) * t485;
t9 = t47 + (t21 * t489 - t22 * t491 + (t74 * t489 + t73 * t491) * qJD(1)) * t485;
t8 = t21 * t416 + t22 * t414 - t73 * t315 + t74 * t317 + t621;
t6 = t32 * t486 + (t18 * t489 - t19 * t491 + (t489 * t67 + t491 * t66) * qJD(1)) * t485;
t5 = t31 * t486 + (t16 * t489 - t17 * t491 + (t489 * t69 + t491 * t68) * qJD(1)) * t485;
t2 = t100 * t350 + t18 * t363 - t19 * t365 + t234 * t67 + t236 * t66 + t32 * t424;
t1 = t101 * t350 + t16 * t363 - t17 * t365 + t234 * t69 + t236 * t68 + t31 * t424;
t10 = [t48 + t76 + (t118 * t58 + t119 * t57) * t553 + (t123 * t206 + t124 * t205) * t557 + (t173 * t94 + t174 * t93) * t555 + (t189 * t264 + t190 * t263) * t559 + (t280 * t361 + t281 * t360) * t561 + t674; t75 + m(3) * (t260 * t360 + t261 * t361 + t280 * t335 + t281 * t336) + (t125 * t263 + t126 * t264 + t189 * t202 + t190 * t203) * m(4) + (t123 * t162 + t124 * t163 + t205 * t95 + t206 * t96) * m(5) + (t116 * t93 + t117 * t94 + t173 * t61 + t174 * t62) * m(6) + (t118 * t36 + t119 * t37 + t57 * t78 + t58 * t79) * m(7) + t47 + ((-t141 / 0.2e1 - t89 / 0.2e1 - t91 / 0.2e1 - t110 / 0.2e1 - t113 / 0.2e1 - t192 / 0.2e1 - t520) * t491 + (t111 / 0.2e1 + t112 / 0.2e1 + t191 / 0.2e1 + t142 / 0.2e1 + t90 / 0.2e1 + t92 / 0.2e1 + t521) * t489 + ((t251 / 0.2e1 + t167 / 0.2e1 + t170 / 0.2e1 + t193 / 0.2e1 + t275 / 0.2e1 + t196 / 0.2e1 - t517) * t491 + (t194 / 0.2e1 + t195 / 0.2e1 + t274 / 0.2e1 + t250 / 0.2e1 + t168 / 0.2e1 + t169 / 0.2e1 + t516) * t489) * qJD(1)) * t485 + t657; (t203 * t125 + t202 * t126 + t598 * t88) * t559 + (t20 * t65 + t36 * t79 + t37 * t78) * t553 + (t116 * t62 + t117 * t61 + t40 * t97) * t555 + (t122 * t59 + t162 * t96 + t163 * t95) * t557 + (t336 * t260 + t335 * t261) * t561 + (t13 + t6) * t614 + (-t12 - t5) * t613 + (t56 + t29) * t566 + (t55 + t28) * t565 + ((t299 * t489 + t300 * t491) * t88 * t559 + (t399 * t489 + t400 * t491) * (t324 * t491 + t325 * t489 + (t399 * t491 - t400 * t489) * qJD(1)) * t561 * t485 + (-t489 * t656 + t491 * t658) * t566 + (t489 * t659 - t491 * t655) * t565 + (t655 * t595 + t659 * t594 + (t314 * t666 + t315 * t664 - t462 * t321 + t323 * t509 - t419 * t394 - t420 * t396 - t415 * t669 - t416 * t667 - t565 * t653) * t491 + (t462 * t320 - t509 * t322 + t419 * t395 + t420 * t397 + (t594 * t654 + t652) * t485 + t668 * t416 - t670 * t415 + t663 * t315 + t665 * t314) * t489) * t614 + (t658 * t595 + t656 * t594 + (t316 * t665 + t317 * t663 + t320 * t508 + t322 * t507 - t421 * t395 - t422 * t397 + t413 * t670 - t414 * t668 - t566 * t654) * t489 + (-t508 * t321 - t507 * t323 + t421 * t394 + t422 * t396 + (t595 * t653 + t652) * t485 + t667 * t414 + t669 * t413 + t664 * t317 + t666 * t316) * t491) * t613) * t485 + (t15 + t9 + (t191 + t112 + t111) * t614 + (-t192 - t113 - t110) * t613 + (t274 + t195 + t194) * t566 + (t275 + t196 + t193) * t565 + ((-t141 - t89 - t91) * t491 + (t142 + t90 + t92) * t489 + ((t167 + t170 + t251) * t491 + (t168 + t169 + t250) * t489) * qJD(1)) * t485 + t657) * t486; ((t118 * t594 + t119 * t595 + t489 * t58 - t491 * t57) * m(7) + (t173 * t594 + t174 * t595 + t489 * t94 - t491 * t93) * m(6) + (-t123 * t491 + t124 * t489 + t205 * t594 + t206 * t595) * m(5) + (-t189 * t491 + t190 * t489 + t263 * t594 + t264 * t595) * m(4)) * t485; (m(4) * t88 + t515) * t486 + ((t36 * t489 - t37 * t491 + t594 * t79 + t595 * t78) * m(7) + (t116 * t595 + t117 * t594 + t489 * t61 - t491 * t62) * m(6) + (t162 * t595 + t163 * t594 + t489 * t95 - t491 * t96) * m(5) + (t125 * t489 - t126 * t491 + t202 * t595 + t203 * t594) * m(4)) * t485; 0; (-m(5) * t124 - m(6) * t94 - m(7) * t58) * t415 + (-m(5) * t123 - m(6) * t93 - m(7) * t57) * t413 + (-m(5) * t206 - m(6) * t174 - m(7) * t119) * t316 + (m(5) * t205 + m(6) * t173 + m(7) * t118) * t314; t515 * t451 + (m(5) * t122 + m(6) * t97 + m(7) * t65) * t439 + (-m(5) * t95 - m(6) * t61 - m(7) * t36) * t415 + (-m(5) * t96 - m(6) * t62 - m(7) * t37) * t413 + (-m(5) * t162 - m(6) * t116 - m(7) * t78) * t316 + (m(5) * t163 + m(6) * t117 + m(7) * t79) * t314; t550 * (t439 * t651 + (t314 * t650 + 0.2e1 * t316 * t491 + (-t413 * t489 - t415 * t491) * t590) * t485); 0.4e1 * t550 * (-t314 * t415 + t316 * t413 + t439 * t451); (t118 * t34 + t119 * t33 + t57 * t99 + t58 * t98) * m(7) + (t173 * t71 + t174 * t70 + t177 * t94 + t178 * t93) * m(6) + t521 * t416 + t520 * t414 + t516 * t317 + t517 * t315 + t620 + t621; (t20 * t77 + t23 * t65 + t33 * t78 + t34 * t79 + t36 * t98 + t37 * t99) * m(7) + t39 * t641 + (t116 * t70 + t117 * t71 + t140 * t40 + t177 * t61 + t178 * t62 + t60 * t97) * m(6) + (t9 / 0.2e1 + t15 / 0.2e1) * t452 + (t6 / 0.2e1 + t13 / 0.2e1) * t416 + (t5 / 0.2e1 + t12 / 0.2e1) * t414 + (t29 / 0.2e1 + t56 / 0.2e1) * t317 + (-t28 / 0.2e1 - t55 / 0.2e1) * t315 + (t8 / 0.2e1 + t14 / 0.2e1 - t151 / 0.2e1) * t486 + ((t114 * t489 - t115 * t491) * t641 + ((t26 / 0.2e1 + t53 / 0.2e1) * t491 + (t27 / 0.2e1 + t54 / 0.2e1) * t489) * qJD(1) + t672 * t640 + t673 * t639) * t485; t540 * t486 + ((t177 * t594 + t178 * t595 + t489 * t71 - t491 * t70) * m(6) + (-t33 * t491 + t34 * t489 + t594 * t98 + t595 * t99) * m(7)) * t485; t540 * t451 + (m(6) * t140 + m(7) * t77) * t439 + (-m(6) * t71 - m(7) * t34) * t415 + (-m(6) * t70 - m(7) * t33) * t413 + (-m(6) * t178 - m(7) * t99) * t316 + (m(6) * t177 + m(7) * t98) * t314; (t23 * t77 + t33 * t99 + t34 * t98) * t553 - t440 * t38 + (t140 * t60 + t177 * t71 + t178 * t70) * t555 + (t14 + t8 - t151) * t452 + (-t114 * t440 + t672) * t416 + (-t115 * t440 + t673) * t414 + (t27 + t54) * t317 + (-t26 - t53) * t315; (t118 * t49 + t119 * t50 + t120 * t57 + t121 * t58) * m(7) + t586 * t365 + t587 * t363 + t585 * t236 + t584 * t234 + t622; t29 * t647 + t5 * t643 + t28 * t646 + t6 * t644 + (t103 * t20 + t120 * t37 + t121 * t36 + t30 * t65 + t49 * t79 + t50 * t78) * m(7) + t486 * t648 + t39 * t645 + t9 * t642 + (t2 * t640 + t1 * t639 + (t491 * t24 / 0.2e1 + t25 * t640) * qJD(1)) * t485; (t30 * t651 + (t49 * t650 - 0.2e1 * t491 * t50 + (t120 * t489 + t121 * t491) * t590) * t485) * t649; (t103 * t439 - t120 * t316 + t121 * t314 + t30 * t451 - t413 * t50 - t415 * t49) * m(7); (t103 * t23 + t120 * t33 + t121 * t34 + t30 * t77 + t49 * t98 + t50 * t99) * m(7) + t27 * t647 + t3 * t643 + t26 * t646 + t4 * t644 + t35 * t641 + t452 * t648 - t315 * t24 / 0.2e1 + t416 * t2 / 0.2e1 + t38 * t645 + t8 * t642 + t317 * t25 / 0.2e1 + t414 * t1 / 0.2e1; t236 * t24 + t363 * t2 + t234 * t25 - t365 * t1 + t350 * t35 + t424 * t7 + (t103 * t30 + t120 * t50 + t121 * t49) * t553;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;