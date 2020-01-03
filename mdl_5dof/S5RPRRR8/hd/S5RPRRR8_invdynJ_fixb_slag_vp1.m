% Calculate vector of inverse dynamics joint torques for
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR8_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR8_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:43
% EndTime: 2019-12-31 19:06:05
% DurationCPUTime: 17.17s
% Computational Cost: add. (15262->674), mult. (27292->857), div. (0->0), fcn. (28803->8), ass. (0->338)
t564 = sin(qJ(3));
t565 = sin(qJ(1));
t566 = cos(qJ(3));
t567 = cos(qJ(1));
t286 = -t565 * t564 - t567 * t566;
t337 = qJD(1) - qJD(3);
t214 = t337 * t286;
t287 = t567 * t564 - t565 * t566;
t215 = t337 * t287;
t339 = sin(qJ(4));
t465 = qJD(4) * t339;
t446 = t286 * t465;
t246 = pkin(4) * t446;
t340 = cos(qJ(4));
t562 = pkin(4) * t340;
t320 = pkin(3) + t562;
t341 = -pkin(8) - pkin(7);
t448 = -t214 * t341 + t215 * t320 + t246;
t336 = qJD(4) + qJD(5);
t216 = t286 * t336;
t338 = qJ(4) + qJ(5);
t326 = sin(t338);
t327 = cos(t338);
t422 = rSges(6,1) * t326 + rSges(6,2) * t327;
t435 = t287 * pkin(3) + t286 * pkin(7);
t479 = t286 * t341 - t287 * t320;
t126 = t435 + t479;
t512 = t287 * t327;
t513 = t287 * t326;
t151 = -rSges(6,1) * t512 + rSges(6,2) * t513 - t286 * rSges(6,3);
t498 = t126 + t151;
t584 = t216 * t422 - t337 * t498 + t246;
t507 = t326 * t336;
t389 = t215 * t327 + t286 * t507;
t506 = t327 * t336;
t390 = -t215 * t326 + t286 * t506;
t66 = rSges(6,1) * t389 + t390 * rSges(6,2) + t214 * rSges(6,3);
t651 = t584 - t448 - t66;
t123 = -t214 * pkin(3) + t215 * pkin(7);
t170 = qJD(4) * t215 - qJDD(4) * t286;
t282 = -rSges(6,1) * t327 + t326 * rSges(6,2);
t233 = t282 * t336;
t335 = qJDD(1) - qJDD(3);
t343 = qJD(1) ^ 2;
t329 = t567 * qJ(2);
t461 = t565 * pkin(1);
t301 = t461 - t329;
t325 = qJD(2) * t567;
t470 = t567 * pkin(1) + t565 * qJ(2);
t370 = -qJDD(1) * t301 + qJDD(2) * t565 + (-qJD(1) * t470 + 0.2e1 * t325) * qJD(1);
t346 = (-qJDD(1) * t565 - t343 * t567) * pkin(2) + t370;
t505 = t340 * qJD(4) ^ 2;
t444 = t287 * t465;
t428 = pkin(4) * t444;
t382 = -t214 * t320 - t215 * t341 + t428;
t57 = t382 - t123;
t391 = -t214 * t327 + t287 * t507;
t392 = t214 * t326 + t287 * t506;
t65 = rSges(6,1) * t391 + t392 * rSges(6,2) + t215 * rSges(6,3);
t91 = qJD(5) * t215 - qJDD(5) * t286 + t170;
t17 = -t216 * t233 - t91 * t422 + (-t170 * t339 + t286 * t505) * pkin(4) + (-t123 - t57 - t65) * t337 + (t435 - t498) * t335 + t346;
t650 = t17 - g(1);
t423 = rSges(5,1) * t339 + rSges(5,2) * t340;
t304 = -rSges(5,1) * t340 + rSges(5,2) * t339;
t289 = t304 * qJD(4);
t467 = qJD(4) * t289;
t510 = t287 * t340;
t511 = t287 * t339;
t167 = -rSges(5,1) * t510 + rSges(5,2) * t511 - t286 * rSges(5,3);
t488 = -t167 + t435;
t387 = -t214 * t340 + t444;
t380 = rSges(5,1) * t387 + t215 * rSges(5,3);
t464 = qJD(4) * t340;
t528 = t214 * t339;
t388 = t287 * t464 + t528;
t595 = t388 * rSges(5,2);
t75 = t380 + t595;
t30 = -t286 * t467 - t170 * t423 + (-t123 - t75) * t337 + t488 * t335 + t346;
t649 = t30 - g(1);
t217 = t287 * t336;
t477 = t286 * pkin(3) - t287 * pkin(7);
t590 = t286 * t320 + t287 * t341;
t127 = t477 - t590;
t519 = t286 * t327;
t520 = t286 * t326;
t152 = -rSges(6,1) * t519 + rSges(6,2) * t520 + t287 * rSges(6,3);
t497 = t127 + t152;
t648 = -t217 * t422 - t428 - (-t477 + t497) * t337;
t169 = qJD(4) * t214 + qJDD(4) * t287;
t333 = t567 * pkin(2);
t324 = qJD(2) * t565;
t441 = qJD(1) * t565;
t442 = qJD(1) * t567;
t473 = qJ(2) * t442 + t324;
t375 = qJDD(1) * t470 - qJDD(2) * t567 + (-pkin(1) * t441 + t324 + t473) * qJD(1);
t460 = t565 * pkin(2);
t350 = qJDD(1) * t333 - t343 * t460 + t375;
t487 = t215 * pkin(3) + t214 * pkin(7);
t347 = -t335 * t477 + t337 * t487 + t350;
t58 = t448 - t487;
t558 = t58 + t66;
t90 = qJD(5) * t214 + qJDD(5) * t287 + t169;
t18 = -t217 * t233 + t90 * t422 + t558 * t337 + t497 * t335 + (t169 * t339 + t287 * t505) * pkin(4) + t347;
t647 = -g(2) + t18;
t405 = Icges(6,5) * t327 - Icges(6,6) * t326;
t137 = -Icges(6,3) * t287 + t286 * t405;
t547 = Icges(6,4) * t326;
t413 = Icges(6,1) * t327 - t547;
t146 = Icges(6,5) * t286 + t287 * t413;
t546 = Icges(6,4) * t327;
t409 = -Icges(6,2) * t326 + t546;
t142 = Icges(6,6) * t286 + t287 * t409;
t533 = t142 * t326;
t635 = -t146 * t327 + t533;
t645 = -t137 + t635;
t138 = Icges(6,3) * t286 + t287 * t405;
t145 = -Icges(6,5) * t287 + t286 * t413;
t141 = -Icges(6,6) * t287 + t286 * t409;
t532 = t141 * t326;
t636 = -t145 * t327 + t532;
t644 = t138 + t636;
t592 = t337 * t167;
t385 = t215 * t340 + t446;
t445 = t286 * t464;
t386 = -t215 * t339 + t445;
t76 = rSges(5,1) * t385 + t386 * rSges(5,2) + t214 * rSges(5,3);
t641 = -t487 - t76 - t592;
t504 = t138 * t286 + t146 * t512;
t640 = t286 * t636 + t504 + (t137 - t533) * t287;
t502 = -t137 * t287 + t145 * t519;
t639 = t287 * t635 - (t138 + t532) * t286 + t502;
t121 = -t214 * rSges(4,1) - t215 * rSges(4,2);
t478 = t287 * rSges(4,1) - t286 * rSges(4,2);
t55 = -t337 * t121 + t335 * t478 + t346;
t638 = t55 - g(1);
t193 = t423 * t287;
t203 = t337 * t435;
t369 = t324 + (-t460 - t301) * qJD(1);
t362 = t203 + t369;
t466 = qJD(4) * t423;
t79 = t286 * t466 + t362 - t592;
t637 = t193 * t79;
t605 = t151 + t479;
t517 = t286 * t340;
t480 = rSges(5,1) * t517 - t287 * rSges(5,3);
t518 = t286 * t339;
t168 = rSges(5,2) * t518 - t480;
t603 = t168 - t477;
t549 = Icges(5,4) * t339;
t415 = Icges(5,1) * t340 - t549;
t161 = -Icges(5,5) * t287 + t286 * t415;
t548 = Icges(5,4) * t340;
t411 = -Icges(5,2) * t339 + t548;
t157 = -Icges(5,6) * t287 + t286 * t411;
t530 = t157 * t339;
t396 = t161 * t340 - t530;
t614 = t590 - t152;
t162 = Icges(5,5) * t286 + t287 * t415;
t158 = Icges(5,6) * t286 + t287 * t411;
t531 = t158 * t339;
t633 = -t162 * t340 + t531;
t407 = Icges(5,5) * t340 - Icges(5,6) * t339;
t153 = -Icges(5,3) * t287 + t286 * t407;
t154 = Icges(5,3) * t286 + t287 * t407;
t500 = -t154 * t286 - t162 * t510;
t631 = -(t153 - t531) * t287 + t500;
t397 = -t157 * t340 - t161 * t339;
t399 = -t158 * t340 - t162 * t339;
t468 = qJD(4) * t287;
t469 = qJD(4) * t286;
t569 = -t337 / 0.2e1;
t630 = ((-t397 * t286 + t399 * t287) * qJD(4) + t397 * t469 - t399 * t468) * t569;
t122 = t215 * rSges(4,1) - t214 * rSges(4,2);
t220 = rSges(4,1) * t286 + rSges(4,2) * t287;
t56 = t337 * t122 - t220 * t335 + t350;
t629 = t56 - g(2);
t626 = t153 * t287;
t625 = t154 * t287;
t622 = t337 * t478;
t617 = t154 + t530;
t615 = qJD(1) * t301 - t324 + t473;
t85 = t141 * t327 + t145 * t326;
t84 = t142 * t327 + t146 * t326;
t612 = -t216 / 0.2e1;
t611 = t217 / 0.2e1;
t503 = t137 * t286 + t145 * t512;
t47 = -t141 * t513 + t503;
t610 = t47 * t216;
t501 = -t138 * t287 + t146 * t519;
t48 = -t142 * t520 + t501;
t609 = t48 * t217;
t606 = t220 * t337;
t406 = Icges(5,5) * t339 + Icges(5,6) * t340;
t188 = t406 * t286;
t187 = t406 * t287;
t412 = Icges(6,1) * t326 + t546;
t601 = -t409 - t412;
t437 = pkin(4) * t339 + t422;
t377 = -t461 - t460;
t600 = pkin(2) * t441 + t377 * qJD(1) + t615;
t180 = t422 * t287;
t181 = t422 * t286;
t36 = t151 * t217 + t152 * t216 + (t126 * t287 + t127 * t286) * qJD(4);
t44 = t362 + t584;
t447 = t333 + t470;
t368 = t447 * qJD(1) - t325;
t45 = t368 - t648;
t599 = -(t337 * t181 - t217 * t282) * t45 - t36 * (t217 * t180 + t181 * t216) - t44 * (-t180 * t337 - t216 * t282);
t410 = Icges(5,2) * t340 + t549;
t414 = Icges(5,1) * t339 + t548;
t393 = -t339 * t410 + t340 * t414;
t102 = t287 * t393 + t188;
t99 = t102 * t337;
t103 = t286 * t393 - t187;
t100 = t103 * t337;
t471 = t567 * rSges(3,1) + t565 * rSges(3,3);
t591 = t470 + t471;
t587 = t382 + t65;
t583 = t600 - t203;
t490 = t410 * t287 - t162;
t492 = -t414 * t287 - t158;
t582 = t339 * t490 + t340 * t492;
t408 = Icges(6,2) * t327 + t547;
t581 = t217 * (t408 * t286 - t145) - t216 * (t408 * t287 - t146) + t337 * t601;
t576 = t214 / 0.2e1;
t575 = t215 / 0.2e1;
t574 = -t217 / 0.2e1;
t573 = t216 / 0.2e1;
t572 = -t286 / 0.2e1;
t571 = t287 / 0.2e1;
t570 = t335 / 0.2e1;
t568 = t337 / 0.2e1;
t559 = t214 * t151 + t287 * t65;
t80 = t287 * t466 + t337 * t603 + t368;
t550 = t80 * t423;
t404 = Icges(6,5) * t326 + Icges(6,6) * t327;
t523 = t404 * t337;
t508 = t407 * t337;
t499 = t153 * t286 + t161 * t510;
t496 = t162 * t517 - t625;
t495 = t161 * t517 - t626;
t491 = -t414 * t286 - t157;
t489 = t410 * t286 - t161;
t475 = -t410 + t415;
t474 = -t411 - t414;
t463 = -t567 / 0.2e1;
t462 = t565 / 0.2e1;
t451 = t565 * rSges(3,1);
t440 = t469 / 0.2e1;
t439 = -t468 / 0.2e1;
t438 = t468 / 0.2e1;
t434 = -Icges(6,1) * t391 - Icges(6,4) * t392 - Icges(6,5) * t215 - t142 * t336;
t433 = -Icges(6,1) * t389 - Icges(6,4) * t390 - Icges(6,5) * t214 - t141 * t336;
t432 = Icges(6,4) * t391 + Icges(6,2) * t392 + Icges(6,6) * t215 - t146 * t336;
t431 = Icges(6,4) * t389 + Icges(6,2) * t390 + Icges(6,6) * t214 - t145 * t336;
t429 = t601 * t336;
t72 = Icges(5,4) * t385 + Icges(5,2) * t386 + Icges(5,6) * t214;
t74 = Icges(5,1) * t385 + Icges(5,4) * t386 + Icges(5,5) * t214;
t352 = qJD(4) * t397 + t339 * t72 - t340 * t74;
t71 = Icges(5,4) * t387 + Icges(5,2) * t388 + Icges(5,6) * t215;
t73 = Icges(5,1) * t387 + Icges(5,4) * t388 + Icges(5,5) * t215;
t353 = qJD(4) * t399 + t339 * t71 - t340 * t73;
t69 = Icges(5,5) * t387 + Icges(5,6) * t388 + Icges(5,3) * t215;
t70 = Icges(5,5) * t385 + Icges(5,6) * t386 + Icges(5,3) * t214;
t421 = -(-t154 * t215 - t214 * t633 - t286 * t69 + t287 * t353) * t286 + (-t153 * t215 + t214 * t396 - t286 * t70 + t287 * t352) * t287;
t420 = -(-t154 * t214 + t215 * t633 + t286 * t353 + t287 * t69) * t286 + (-t153 * t214 - t215 * t396 + t286 * t352 + t287 * t70) * t287;
t50 = -t158 * t511 - t500;
t51 = -t157 * t511 + t499;
t419 = -t286 * t50 + t287 * t51;
t52 = -t158 * t518 + t496;
t53 = -t157 * t518 + t495;
t418 = -t286 * t52 + t287 * t53;
t417 = t286 * t76 + t287 * t75;
t416 = -t286 * t79 - t287 * t80;
t395 = t167 * t287 + t168 * t286;
t394 = -t326 * t408 + t327 * t412;
t228 = -t339 * t414 - t340 * t410;
t378 = t394 * t337;
t376 = -t451 - t461;
t307 = rSges(2,1) * t567 - rSges(2,2) * t565;
t303 = rSges(2,1) * t565 + rSges(2,2) * t567;
t367 = -t337 * t405 + (-t216 * t287 + t286 * t217) * t404;
t366 = t329 + t377;
t365 = t339 * t489 + t340 * t491;
t364 = t380 + t123;
t360 = (t339 * t474 + t340 * t475) * t337;
t355 = t326 * t431 + t327 * t433;
t60 = Icges(6,5) * t389 + Icges(6,6) * t390 + Icges(6,3) * t214;
t10 = -t137 * t214 + t215 * t636 + t286 * t355 + t287 * t60;
t46 = -t142 * t513 + t504;
t94 = t286 * t404 + t287 * t394;
t92 = t94 * t337;
t21 = -t216 * t46 + t217 * t47 + t92;
t49 = -t141 * t520 + t502;
t95 = t286 * t394 - t287 * t404;
t93 = t95 * t337;
t22 = -t216 * t48 + t217 * t49 + t93;
t28 = t326 * t434 - t327 * t432;
t29 = t326 * t433 - t327 * t431;
t351 = (-t412 * t286 - t141) * t217 - (-t412 * t287 - t142) * t216 + (-t408 + t413) * t337;
t345 = t581 * t326 + t351 * t327;
t230 = t405 * t336;
t232 = t413 * t336;
t354 = (-t336 * t408 + t232) * t327 + t429 * t326;
t37 = t214 * t394 - t215 * t404 + t230 * t286 + t287 * t354;
t38 = -t214 * t404 - t215 * t394 - t230 * t287 + t286 * t354;
t356 = t326 * t432 + t327 * t434;
t59 = Icges(6,5) * t391 + Icges(6,6) * t392 + Icges(6,3) * t215;
t7 = -t138 * t215 - t214 * t635 - t286 * t59 + t287 * t356;
t8 = -t137 * t215 - t214 * t636 - t286 * t60 + t287 * t355;
t9 = -t138 * t214 + t215 * t635 + t286 * t356 + t287 * t59;
t359 = (t10 * t217 - t216 * t9 + t335 * t95 + t337 * t38 + t48 * t91 + t49 * t90) * t571 + (t286 * t345 + t287 * t367) * t574 + (-t286 * t367 + t287 * t345) * t573 + t22 * t576 + (-t216 * t7 + t217 * t8 + t335 * t94 + t337 * t37 + t46 * t91 + t47 * t90) * t572 + t21 * t575 + (t351 * t326 - t581 * t327) * t569 + t90 * (-t286 * t48 + t287 * t49) / 0.2e1 + t91 * (-t286 * t46 + t287 * t47) / 0.2e1 + (t10 * t287 + t214 * t49 + t215 * t48 - t286 * t9) * t611 + (-t286 * t84 + t287 * t85) * t570 + (t214 * t47 + t215 * t46 - t286 * t7 + t287 * t8) * t612 + (t214 * t85 + t215 * t84 - t28 * t286 + t287 * t29) * t568;
t284 = t411 * qJD(4);
t285 = t415 * qJD(4);
t349 = qJD(4) * t228 - t284 * t339 + t285 * t340;
t26 = qJD(4) * t419 + t99;
t27 = qJD(4) * t418 + t100;
t33 = -t633 * qJD(4) - t339 * t73 - t340 * t71;
t34 = qJD(4) * t396 - t339 * t74 - t340 * t72;
t283 = t407 * qJD(4);
t40 = t214 * t393 - t215 * t406 + t283 * t286 + t287 * t349;
t41 = -t214 * t406 - t215 * t393 - t283 * t287 + t286 * t349;
t344 = t26 * t438 - (t85 + t95) * t90 / 0.2e1 - (t84 + t94) * t91 / 0.2e1 - (t103 - t397) * t169 / 0.2e1 - (t102 - t399) * t170 / 0.2e1 + (t38 + t29) * t574 + (t37 + t28) * t573 + (t41 + t34) * t439 + (-t232 * t326 + t408 * t507 + t327 * t429 - t285 * t339 + t410 * t465 + (-qJD(4) * t414 - t284) * t340) * t337 + (t27 + t40 + t33) * t440 + (-t326 * t412 - t327 * t408 - Icges(4,3) + t228) * t335;
t331 = t567 * rSges(3,3);
t322 = rSges(3,3) * t442;
t302 = t451 - t331;
t194 = t423 * t286;
t136 = t368 - t606;
t135 = t369 + t622;
t116 = t287 * t151;
t110 = qJDD(1) * t471 + qJD(1) * (-rSges(3,1) * t441 + t322) + t375;
t109 = -qJDD(1) * t302 - t343 * t471 + t370;
t81 = t395 * qJD(4);
t35 = (-t157 * t287 - t158 * t286) * t339 + t496 + t499;
t32 = (-t141 * t287 - t142 * t286) * t326 + t501 + t503;
t31 = t335 * t168 + t169 * t423 - t287 * t467 + t337 * t76 + t347;
t5 = t126 * t169 - t127 * t170 + t151 * t90 - t152 * t91 + t217 * t65 + t216 * t66 + (t286 * t58 + t287 * t57) * qJD(4);
t1 = [(t610 - t32 * t216 + t93 + (t46 + t639) * t217) * t573 - t344 + t22 * t612 - m(2) * (-g(1) * t303 + g(2) * t307) + (-t99 + ((t52 + (t154 - t396) * t287) * t287 + (t617 * t286 - t495 + t53 - t626) * t286) * qJD(4)) * t439 + (t100 + ((t633 * t287 + t495 + t50) * t287 + (-t617 * t287 - t35 + t51) * t286) * qJD(4)) * t440 + (t609 + (t645 * t217 - t523) * t286 + (t644 * t217 - t378) * t287 + t21 + (t49 - t639) * t216) * t574 + t630 + (m(2) * (t303 ^ 2 + t307 ^ 2) + Icges(2,3) + Icges(3,2)) * qJDD(1) + (t647 * (t447 - t614) + (-t368 - t587) * t44 + (t44 + t583 - t651) * t45 + t650 * (t366 - t605)) * m(6) + (-t550 * t469 + (-g(2) + t31) * (t447 + t603) + (-t364 - t368 - t595) * t79 + (t583 + t79 - t641) * t80 + t649 * (t366 + t488)) * m(5) + (t629 * (-t220 + t447) + (-t121 - t368) * t135 + (t122 + t135 + t600 - t622) * t136 + t638 * (t366 + t478)) * m(4) + ((-g(2) + t110) * t591 + (-g(1) + t109) * (t329 + t331 + t376) + (t322 + (t302 + t376) * qJD(1) + t615) * (qJD(1) * t591 - t325)) * m(3); (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t565 - g(2) * t567) + 0.2e1 * (t17 * t462 + t18 * t463) * m(6) + 0.2e1 * (t30 * t462 + t31 * t463) * m(5) + 0.2e1 * (t462 * t55 + t463 * t56) * m(4) + 0.2e1 * (t109 * t462 + t110 * t463) * m(3); t344 + t21 * t611 + (t99 + ((t35 - t52) * t287 + (t396 * t286 - t53 + t631) * t286) * qJD(4)) * t439 + (-t100 + ((-t50 - t631) * t287 + (-t51 + (t153 - t633) * t286 - t625) * t286) * qJD(4)) * t440 + (t32 * t217 - t609 + t92 - (t49 + t640) * t216) * t574 + (-t610 + (-t644 * t216 + t523) * t287 + (-t645 * t216 - t378) * t286 + t22 + (-t46 + t640) * t217) * t573 + t630 + (t650 * t605 + (t203 + t651) * t45 + (t587 + t648) * t44 + t647 * t614) * m(6) + (-(-t286 * t550 + t637) * qJD(4) + t79 * t364 + t31 * (t477 + t480) + (-t31 * t518 + t388 * t79) * rSges(5,2) + (t203 + t641) * t80 - t649 * t488 - (t79 * t337 - g(2)) * t603) * m(5) + (-t122 * t136 + (t121 + t606) * t135 - (-t337 * t136 + t638) * t478 + t629 * t220) * m(4); -(t214 * t51 + t215 * t50 + t421) * t469 / 0.2e1 + (t214 * t53 + t215 * t52 + t420) * t438 + t359 + (qJD(4) * t420 + t103 * t335 + t169 * t53 + t170 * t52 + t337 * t41) * t571 + (qJD(4) * t421 + t102 * t335 + t169 * t51 + t170 * t50 + t337 * t40) * t572 + t27 * t576 + t26 * t575 + ((t187 * t469 + t508) * t286 + (t360 + (t365 * t287 + (-t188 - t582) * t286) * qJD(4)) * t287) * t440 + (-t214 * t397 - t215 * t399 - t286 * t33 + t287 * t34) * t568 + (t286 * t399 - t287 * t397) * t570 + ((t188 * t468 - t508) * t287 + (t360 + (-t582 * t286 + (-t187 + t365) * t287) * qJD(4)) * t286) * t439 + t169 * t418 / 0.2e1 + t170 * t419 / 0.2e1 + ((t339 * t475 - t340 * t474) * t337 + ((t286 * t490 - t287 * t489) * t340 + (-t286 * t492 + t287 * t491) * t339) * qJD(4)) * t569 + (-g(3) * (t282 - t562) - (g(1) * t286 + g(2) * t287) * t437 - (t44 * t445 + ((t286 * t45 - t287 * t44) * t337 + t36 * (t286 ^ 2 + t287 ^ 2) * qJD(4)) * t339) * pkin(4) + t45 * (pkin(4) * t528 + t214 * t422) + t5 * t116 + t36 * (t126 * t214 + t559) + (t5 * t126 + t18 * t437 - t45 * t233 + t36 * t57) * t287 + (-t36 * t497 - t437 * t44) * t215 + (t17 * t437 + t44 * (pkin(4) * t464 - t233) + t5 * t497 + t36 * t558) * t286 + t599) * m(6) + (-g(1) * t194 - g(2) * t193 - g(3) * t304 + (qJD(4) * t417 + t167 * t169 - t168 * t170) * t395 + t81 * (t167 * t214 - t168 * t215 + t417) + t416 * t289 - (-t214 * t80 + t215 * t79 - t286 * t30 - t287 * t31) * t423 - (t194 * t80 - t637) * t337 - (t81 * (t193 * t287 + t194 * t286) + t416 * t304) * qJD(4)) * m(5); t359 + (t5 * (t152 * t286 + t116) + t36 * (-t152 * t215 + t286 * t66 + t559) + (-t286 * t44 - t287 * t45) * t233 - (-t17 * t286 - t18 * t287 - t214 * t45 + t215 * t44) * t422 - g(1) * t181 - g(2) * t180 - g(3) * t282 + t599) * m(6);];
tau = t1;
