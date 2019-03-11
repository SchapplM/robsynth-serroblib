% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:13
% EndTime: 2019-03-09 01:37:41
% DurationCPUTime: 24.47s
% Computational Cost: add. (15327->866), mult. (35915->1108), div. (0->0), fcn. (39890->8), ass. (0->390)
t339 = sin(qJ(5));
t341 = cos(qJ(5));
t337 = cos(pkin(9));
t342 = cos(qJ(1));
t522 = sin(pkin(9));
t534 = sin(qJ(1));
t426 = t534 * t522;
t269 = t337 * t342 - t426;
t270 = t337 * t534 + t342 * t522;
t340 = cos(qJ(6));
t338 = sin(qJ(6));
t499 = t338 * t341;
t179 = t269 * t499 - t270 * t340;
t169 = Icges(7,4) * t179;
t498 = t340 * t341;
t176 = t269 * t498 + t270 * t338;
t504 = t269 * t339;
t100 = Icges(7,1) * t176 + Icges(7,5) * t504 - t169;
t518 = Icges(7,4) * t176;
t97 = -Icges(7,2) * t179 + Icges(7,6) * t504 + t518;
t593 = t100 * t340 - t338 * t97;
t94 = Icges(7,5) * t176 - Icges(7,6) * t179 + Icges(7,3) * t504;
t34 = -t593 * t339 + t341 * t94;
t461 = qJD(6) * t339;
t469 = qJD(5) * t270;
t197 = t269 * t461 + t469;
t470 = qJD(5) * t269;
t199 = t270 * t461 - t470;
t177 = -t269 * t340 - t270 * t499;
t168 = Icges(7,4) * t177;
t178 = -t269 * t338 + t270 * t498;
t501 = t270 * t339;
t101 = Icges(7,1) * t178 + Icges(7,5) * t501 + t168;
t95 = Icges(7,5) * t178 + Icges(7,6) * t177 + Icges(7,3) * t501;
t519 = Icges(7,4) * t178;
t98 = Icges(7,2) * t177 + Icges(7,6) * t501 + t519;
t24 = t178 * t101 + t177 * t98 + t95 * t501;
t594 = t100 * t178 + t177 * t97;
t25 = -t501 * t94 - t594;
t460 = qJD(6) * t341;
t302 = qJD(1) - t460;
t397 = Icges(7,5) * t340 - Icges(7,6) * t338;
t218 = -Icges(7,3) * t341 + t339 * t397;
t516 = Icges(7,4) * t340;
t398 = -Icges(7,2) * t338 + t516;
t220 = -Icges(7,6) * t341 + t339 * t398;
t517 = Icges(7,4) * t338;
t400 = Icges(7,1) * t340 - t517;
t222 = -Icges(7,5) * t341 + t339 * t400;
t58 = t177 * t220 + t178 * t222 + t218 * t501;
t10 = -t197 * t25 + t199 * t24 + t58 * t302;
t27 = t100 * t176 - t179 * t97 + t504 * t94;
t417 = t176 * rSges(7,1) - t179 * rSges(7,2);
t105 = -rSges(7,3) * t504 - t417;
t529 = rSges(7,2) * t338;
t530 = rSges(7,1) * t340;
t416 = -t529 + t530;
t224 = -rSges(7,3) * t341 + t339 * t416;
t471 = qJD(1) * t342;
t316 = pkin(3) * t471;
t320 = qJD(3) * t342;
t321 = qJD(2) * t534;
t474 = t320 + t321;
t446 = t316 + t474;
t533 = pkin(5) * t339;
t296 = -pkin(8) * t341 + t533;
t465 = qJD(5) * t296;
t590 = t302 * t105 + t197 * t224 + t270 * t465 - t446;
t568 = t176 * t222 - t179 * t220 + t218 * t504;
t589 = t197 * t27 + t568 * t302;
t281 = Icges(6,5) * t341 - Icges(6,6) * t339;
t147 = -Icges(6,3) * t269 + t270 * t281;
t587 = t270 * t147;
t219 = Icges(7,3) * t339 + t341 * t397;
t392 = -t220 * t338 + t222 * t340;
t403 = t101 * t340 - t338 * t98;
t555 = t197 * (-t218 * t269 - t593) - t199 * (t218 * t270 + t403) + t302 * (t219 - t392);
t586 = t555 * t339;
t251 = t270 * qJD(1);
t252 = -qJD(1) * t426 + t337 * t471;
t585 = t252 * pkin(4) + t251 * pkin(7);
t324 = t342 * qJ(3);
t333 = t534 * pkin(3);
t382 = -t333 - t324;
t323 = t534 * qJ(2);
t561 = t342 * pkin(1) + t323;
t367 = t561 - t382;
t322 = qJD(2) * t342;
t436 = qJD(3) * t534;
t419 = t436 - t322;
t584 = -t367 * qJD(1) - t419;
t146 = Icges(6,3) * t270 + t269 * t281;
t326 = Icges(6,4) * t341;
t399 = -Icges(6,2) * t339 + t326;
t149 = Icges(6,6) * t270 + t269 * t399;
t520 = Icges(6,4) * t339;
t285 = Icges(6,1) * t341 - t520;
t152 = Icges(6,5) * t270 + t269 * t285;
t394 = -t149 * t339 + t152 * t341;
t55 = t146 * t270 + t269 * t394;
t78 = -t149 * t341 - t152 * t339;
t286 = rSges(6,1) * t339 + rSges(6,2) * t341;
t466 = qJD(5) * t286;
t581 = t269 * t466;
t443 = t534 * qJ(3);
t455 = t534 * pkin(1);
t366 = -t455 - t443;
t580 = t366 + t443;
t463 = qJD(5) * t339;
t372 = t252 * t341 - t270 * t463;
t462 = qJD(5) * t341;
t506 = t252 * t339;
t373 = t270 * t462 + t506;
t93 = rSges(6,1) * t372 - rSges(6,2) * t373 + t251 * rSges(6,3);
t343 = qJD(1) ^ 2;
t579 = qJDD(1) * t534 + t342 * t343;
t330 = t534 * rSges(4,1);
t289 = t342 * rSges(4,3) + t330;
t444 = t324 + t561;
t578 = t444 + t289;
t563 = -t342 * rSges(3,2) + t534 * rSges(3,3);
t230 = t561 + t563;
t334 = t341 * pkin(5);
t562 = t339 * pkin(8) + t334;
t577 = qJD(5) * t562;
t564 = Icges(6,1) * t339 + t326;
t478 = t564 + t399;
t282 = Icges(6,2) * t341 + t520;
t479 = t282 - t285;
t575 = (t339 * t478 + t341 * t479) * qJD(1);
t242 = t252 * pkin(7);
t571 = t242 + t584;
t309 = qJ(2) * t471;
t429 = t309 + t446;
t265 = t270 * pkin(7);
t431 = t269 * pkin(4) + t265;
t570 = -qJD(1) * t431 + t429 + t585;
t569 = t101 * t176 - t179 * t98;
t205 = -t269 * rSges(5,1) + t270 * rSges(5,2);
t565 = -0.2e1 * qJD(1) * t436 + qJDD(3) * t342;
t389 = t282 * t339 - t341 * t564;
t280 = Icges(6,5) * t339 + Icges(6,6) * t341;
t502 = t270 * t280;
t121 = t269 * t389 - t502;
t560 = qJD(1) * t121;
t150 = -Icges(6,6) * t269 + t270 * t399;
t153 = -Icges(6,5) * t269 + t270 * t285;
t77 = t150 * t341 + t153 * t339;
t88 = Icges(6,4) * t372 - Icges(6,2) * t373 + Icges(6,6) * t251;
t90 = Icges(6,1) * t372 - Icges(6,4) * t373 + Icges(6,5) * t251;
t559 = qJD(5) * t77 + t339 * t88 - t341 * t90;
t374 = t251 * t341 + t269 * t463;
t507 = t251 * t339;
t375 = t269 * t462 - t507;
t87 = Icges(6,4) * t374 + Icges(6,2) * t375 - Icges(6,6) * t252;
t89 = Icges(6,1) * t374 + Icges(6,4) * t375 - Icges(6,5) * t252;
t558 = qJD(5) * t78 + t339 * t87 - t341 * t89;
t272 = t399 * qJD(5);
t273 = t285 * qJD(5);
t390 = t282 * t341 + t339 * t564;
t557 = qJD(5) * t390 + t272 * t339 - t273 * t341;
t256 = (-Icges(7,2) * t340 - t517) * t339;
t554 = t197 * (Icges(7,2) * t176 - t100 + t169) - t199 * (-Icges(7,2) * t178 + t101 + t168) - t302 * (t222 + t256);
t553 = -t10 / 0.2e1;
t26 = -t95 * t504 - t569;
t11 = t199 * t26 - t589;
t552 = t11 / 0.2e1;
t195 = -qJD(5) * t252 - qJDD(5) * t270;
t457 = qJDD(6) * t339;
t106 = -qJD(6) * t375 - t269 * t457 + t195;
t551 = t106 / 0.2e1;
t196 = qJD(5) * t251 - qJDD(5) * t269;
t107 = qJD(6) * t373 + t270 * t457 + t196;
t550 = t107 / 0.2e1;
t549 = t195 / 0.2e1;
t548 = t196 / 0.2e1;
t547 = t197 / 0.2e1;
t546 = -t197 / 0.2e1;
t545 = -t199 / 0.2e1;
t544 = t199 / 0.2e1;
t268 = qJD(5) * t461 - qJDD(6) * t341 + qJDD(1);
t541 = t268 / 0.2e1;
t538 = -t302 / 0.2e1;
t537 = t302 / 0.2e1;
t536 = -t342 / 0.2e1;
t535 = rSges(7,3) + pkin(8);
t356 = qJD(6) * t269 - t372;
t439 = t270 * t460;
t421 = t251 - t439;
t83 = t338 * t356 + t340 * t421;
t84 = t338 * t421 - t340 * t356;
t43 = Icges(7,5) * t84 + Icges(7,6) * t83 + Icges(7,3) * t373;
t45 = Icges(7,4) * t84 + Icges(7,2) * t83 + Icges(7,6) * t373;
t47 = Icges(7,1) * t84 + Icges(7,4) * t83 + Icges(7,5) * t373;
t8 = (qJD(5) * t403 - t43) * t341 + (qJD(5) * t95 - t338 * t45 + t340 * t47 + (-t101 * t338 - t340 * t98) * qJD(6)) * t339;
t532 = t8 * t199;
t357 = -qJD(6) * t270 + t374;
t440 = t269 * t460;
t420 = -t252 + t440;
t81 = -t338 * t357 + t340 * t420;
t82 = t338 * t420 + t340 * t357;
t42 = Icges(7,5) * t82 + Icges(7,6) * t81 - Icges(7,3) * t375;
t44 = Icges(7,4) * t82 + Icges(7,2) * t81 - Icges(7,6) * t375;
t46 = Icges(7,1) * t82 + Icges(7,4) * t81 - Icges(7,5) * t375;
t9 = (-qJD(5) * t593 - t42) * t341 + (-qJD(5) * t94 - t338 * t44 + t340 * t46 + (t100 * t338 + t340 * t97) * qJD(6)) * t339;
t531 = t9 * t197;
t528 = t252 * rSges(6,3);
t260 = t270 * rSges(6,3);
t33 = t339 * t403 - t341 * t95;
t527 = t33 * t107;
t327 = t339 * rSges(7,3);
t526 = t34 * t106;
t124 = -t218 * t341 + t339 * t392;
t255 = (-Icges(7,5) * t338 - Icges(7,6) * t340) * t339;
t162 = qJD(5) * t219 + qJD(6) * t255;
t221 = Icges(7,6) * t339 + t341 * t398;
t163 = qJD(5) * t221 + qJD(6) * t256;
t223 = Icges(7,5) * t339 + t341 * t400;
t257 = (-Icges(7,1) * t338 - t516) * t339;
t164 = qJD(5) * t223 + qJD(6) * t257;
t41 = (qJD(5) * t392 - t162) * t341 + (qJD(5) * t218 - t163 * t338 + t164 * t340 + (-t220 * t340 - t222 * t338) * qJD(6)) * t339;
t524 = t124 * t268 + t41 * t302;
t113 = pkin(5) * t372 + pkin(8) * t373;
t49 = t84 * rSges(7,1) + t83 * rSges(7,2) + rSges(7,3) * t373;
t523 = t113 + t49;
t508 = t150 * t339;
t182 = t269 * t280;
t503 = t269 * t341;
t500 = t270 * t341;
t104 = t178 * rSges(7,1) + t177 * rSges(7,2) + rSges(7,3) * t501;
t192 = pkin(5) * t500 + pkin(8) * t501;
t494 = t104 + t192;
t194 = t562 * t269;
t493 = t105 - t194;
t492 = -t153 * t503 - t587;
t491 = t269 * t147 - t153 * t500;
t490 = t270 * t564 + t150;
t489 = t269 * t564 + t149;
t488 = -t282 * t270 + t153;
t487 = t282 * t269 - t152;
t267 = (-rSges(7,1) * t338 - rSges(7,2) * t340) * t339;
t165 = qJD(6) * t267 + (t341 * t416 + t327) * qJD(5);
t486 = t165 + t577;
t250 = qJD(1) * t561 - t322;
t485 = -t251 * pkin(4) + t242 - t250;
t453 = t339 * t529;
t484 = rSges(7,3) * t500 + t270 * t453;
t483 = -rSges(7,3) * t503 - t269 * t453;
t481 = t224 + t296;
t480 = t252 * rSges(5,1) - t251 * rSges(5,2);
t190 = pkin(5) * t503 + pkin(8) * t504;
t204 = t270 * rSges(5,1) + t269 * rSges(5,2);
t477 = qJD(1) * t322 + qJDD(2) * t534;
t476 = t309 + t321;
t437 = qJD(1) * t534;
t475 = rSges(3,2) * t437 + rSges(3,3) * t471;
t325 = t342 * qJ(2);
t287 = t455 - t325;
t279 = qJD(1) * t287;
t473 = t321 - t279;
t472 = qJD(1) * t281;
t291 = rSges(6,1) * t341 - rSges(6,2) * t339;
t274 = t291 * qJD(5);
t468 = qJD(5) * t274;
t467 = qJD(5) * t577;
t459 = -m(5) - m(6) - m(7);
t458 = qJDD(1) * t342;
t456 = t534 / 0.2e1;
t454 = t339 * t530;
t451 = -m(4) + t459;
t328 = t534 * rSges(4,3);
t445 = t320 + t473;
t435 = -t470 / 0.2e1;
t434 = t470 / 0.2e1;
t433 = -t469 / 0.2e1;
t432 = t469 / 0.2e1;
t206 = t270 * pkin(4) - pkin(7) * t269;
t294 = t342 * rSges(4,1) - t328;
t428 = t316 + t445;
t427 = t333 + t444;
t425 = t82 * rSges(7,1) + t81 * rSges(7,2);
t423 = -t287 - t443;
t422 = -t443 + t294;
t418 = -t251 * rSges(5,1) - t252 * rSges(5,2);
t395 = t153 * t341 - t508;
t85 = Icges(6,5) * t374 + Icges(6,6) * t375 - Icges(6,3) * t252;
t86 = Icges(6,5) * t372 - Icges(6,6) * t373 + Icges(6,3) * t251;
t415 = -(-t147 * t252 + t395 * t251 + t269 * t559 - t270 * t86) * t269 - (t146 * t252 - t394 * t251 + t269 * t558 - t270 * t85) * t270;
t414 = -(t147 * t251 + t395 * t252 - t269 * t86 - t270 * t559) * t269 - (-t146 * t251 - t394 * t252 - t269 * t85 - t270 * t558) * t270;
t413 = t24 * t270 - t25 * t269;
t412 = t26 * t270 - t269 * t27;
t411 = -t269 * t34 + t270 * t33;
t52 = -t150 * t501 - t491;
t53 = t146 * t269 + t149 * t501 - t152 * t500;
t410 = -t269 * t52 - t270 * t53;
t54 = t150 * t504 + t492;
t409 = -t269 * t54 - t270 * t55;
t157 = -t269 * t291 - t260;
t387 = t423 + t431;
t370 = -t157 + t387;
t65 = qJD(1) * t370 - t270 * t466 + t446;
t156 = rSges(6,1) * t500 - rSges(6,2) * t501 - t269 * rSges(6,3);
t361 = t206 + t367;
t66 = t581 + (t156 + t361) * qJD(1) + t419;
t408 = t269 * t66 - t270 * t65;
t92 = rSges(6,1) * t374 + rSges(6,2) * t375 - t528;
t407 = -t269 * t92 + t270 * t93;
t396 = t104 * t269 + t105 * t270;
t393 = t156 * t270 - t157 * t269;
t155 = rSges(6,1) * t503 - rSges(6,2) * t504 + t260;
t225 = rSges(7,1) * t498 - rSges(7,2) * t499 + t327;
t388 = pkin(4) + t291;
t386 = -t205 + t423;
t384 = qJDD(1) * t561 - qJDD(2) * t342 + (-pkin(1) * t437 + t321 + t476) * qJD(1);
t295 = t342 * rSges(2,1) - rSges(2,2) * t534;
t290 = rSges(2,1) * t534 + t342 * rSges(2,2);
t288 = rSges(3,2) * t534 + t342 * rSges(3,3);
t380 = t339 * t42 - t462 * t94;
t379 = t339 * t43 + t462 * t95;
t369 = t194 + t387;
t365 = t206 + t427;
t364 = t197 * t94 + t199 * t95 + t218 * t302;
t363 = -qJD(1) * t250 - qJDD(1) * t287 + t477;
t362 = (Icges(7,5) * t177 - Icges(7,6) * t178) * t199 - (Icges(7,5) * t179 + Icges(7,6) * t176) * t197 + t255 * t302;
t360 = t325 + t366;
t359 = t339 * t488 + t341 * t490;
t358 = -t339 * t487 + t341 * t489;
t355 = t342 * pkin(3) + t360;
t354 = t265 + t355;
t352 = (Icges(7,1) * t177 - t519 - t98) * t199 - (Icges(7,1) * t179 + t518 + t97) * t197 + (-t220 + t257) * t302;
t350 = t355 + t431;
t349 = pkin(3) * t458 + t343 * t382 + t477 + t565;
t348 = qJ(3) * t458 + 0.2e1 * qJD(1) * t320 + qJDD(3) * t534 - t343 * t443 + t384;
t30 = t104 * t197 + t105 * t199 + qJD(4) + (t192 * t270 + t194 * t269) * qJD(5);
t38 = qJD(1) * t369 - t590;
t39 = t269 * t465 + t302 * t104 - t199 * t224 + (t192 + t361) * qJD(1) + t419;
t347 = t30 * t396 + (-t269 * t38 - t270 * t39) * t224;
t346 = pkin(3) * t579 + t348;
t345 = qJD(1) * t585 + qJDD(1) * t206 + t346;
t314 = rSges(4,1) * t471;
t271 = t281 * qJD(5);
t249 = pkin(8) * t503;
t247 = pkin(8) * t500;
t193 = pkin(5) * t504 - t249;
t191 = -pkin(5) * t501 + t247;
t188 = t286 * t269;
t187 = t286 * t270;
t166 = (-t287 + t422) * qJD(1) + t474;
t145 = t269 * t454 + t483;
t144 = -t270 * t454 + t484;
t143 = t222 * t269;
t142 = t222 * t270;
t141 = t220 * t269;
t140 = t220 * t270;
t133 = qJD(1) * t386 + t446;
t129 = qJD(1) * t475 + qJDD(1) * t563 + t384;
t128 = qJDD(1) * t288 - t343 * t563 + t363;
t123 = rSges(7,1) * t179 + rSges(7,2) * t176;
t122 = rSges(7,1) * t177 - rSges(7,2) * t178;
t120 = -t270 * t389 - t182;
t112 = pkin(5) * t374 - pkin(8) * t375;
t111 = t120 * qJD(1);
t109 = qJDD(1) * t289 + qJD(1) * (-rSges(4,3) * t437 + t314) + t348;
t108 = -qJ(3) * t579 + qJDD(1) * t294 - t343 * t289 + t363 + t565;
t67 = qJD(5) * t393 + qJD(4);
t61 = qJD(1) * t480 + qJDD(1) * t204 + t346;
t60 = (-t250 + t418) * qJD(1) + t386 * qJDD(1) + t349;
t51 = t251 * t280 - t389 * t252 - t269 * t271 - t270 * t557;
t50 = -t389 * t251 - t252 * t280 + t269 * t557 - t270 * t271;
t48 = -rSges(7,3) * t375 + t425;
t32 = -qJD(5) * t394 + t339 * t89 + t341 * t87;
t31 = qJD(5) * t395 + t339 * t90 + t341 * t88;
t29 = qJD(1) * t93 + qJDD(1) * t156 - t196 * t286 + t269 * t468 + t345;
t28 = -t270 * t468 + t195 * t286 + (-t92 + t485) * qJD(1) + t370 * qJDD(1) + t349;
t23 = qJD(5) * t407 - t156 * t195 + t157 * t196 + qJDD(4);
t22 = qJD(5) * t409 + t560;
t21 = qJD(5) * t410 + t111;
t20 = t162 * t501 + t163 * t177 + t164 * t178 + t218 * t373 + t220 * t83 + t222 * t84;
t19 = -t162 * t504 + t163 * t179 - t164 * t176 - t218 * t375 + t220 * t81 + t222 * t82;
t14 = t124 * t302 - t197 * t34 + t199 * t33;
t13 = -t270 * t467 - t268 * t105 + t106 * t224 - t197 * t165 + t195 * t296 - t302 * t48 + (-t112 + t485) * qJD(1) + t369 * qJDD(1) + t349;
t12 = qJD(1) * t113 + qJDD(1) * t192 + t268 * t104 - t107 * t224 - t199 * t165 - t196 * t296 + t269 * t467 + t302 * t49 + t345;
t7 = -t104 * t106 + t105 * t107 - t192 * t195 - t194 * t196 + t197 * t49 + t199 * t48 + qJDD(4) + (-t112 * t269 + t113 * t270) * qJD(5);
t6 = -t100 * t84 + t177 * t44 + t178 * t46 + t270 * t380 - t506 * t94 - t83 * t97;
t5 = t101 * t84 + t177 * t45 + t178 * t47 + t270 * t379 + t506 * t95 + t83 * t98;
t4 = -t100 * t82 - t176 * t46 + t179 * t44 - t269 * t380 - t507 * t94 - t81 * t97;
t3 = t101 * t82 - t176 * t47 + t179 * t45 - t269 * t379 + t507 * t95 + t81 * t98;
t2 = t106 * t25 + t107 * t24 - t197 * t6 + t199 * t5 + t20 * t302 + t268 * t58;
t1 = t106 * t27 + t107 * t26 + t19 * t302 - t197 * t4 + t199 * t3 - t268 * t568;
t15 = [(-t560 + ((t52 + (t146 + t508) * t270 + t491) * t270 + (-t53 + (t146 - t395) * t269 - t587) * t269) * qJD(5) + t22) * t434 + (t32 + t50 + t21) * t433 + (t120 + t77) * t548 + (-qJD(5) * t389 + t272 * t341 + t273 * t339) * qJD(1) + ((-g(2) + t129) * t230 + (-g(1) + t128) * (-t287 + t288) + (-t473 + t475 + t476 + (-t288 - t455) * qJD(1)) * (qJD(1) * t230 - t322)) * m(3) + t58 * t550 + t197 * t553 + t20 * t544 + t19 * t546 + t527 / 0.2e1 + t526 / 0.2e1 + ((-g(2) + t109) * t578 + (-g(1) + t108) * (t294 + t360) + (-t419 + (-t330 - t323 + (-rSges(4,3) - pkin(1) - qJ(3)) * t342) * qJD(1)) * t166 + (t166 - t445 + t309 + t314 + t474 + (-t422 - t328 + t366) * qJD(1)) * (t578 * qJD(1) + t419)) * m(4) + t532 / 0.2e1 - t531 / 0.2e1 + t524 + t10 * t547 - m(2) * (-g(1) * t290 + g(2) * t295) + (t31 + t51) * t435 + (-t30 * (t190 - t194) * t469 - g(1) * (t350 - t105 + t190) + (t12 - g(2)) * (t365 + t494) + (t354 + t417 + (t339 * t535 + pkin(4) + t334) * t269) * t13 + (-t425 + (-pkin(4) - t562 - t327) * t251 + (t341 * t535 - t533) * t470 + t571) * t38 + (t279 + t38 + t523 + (-t190 + t580) * qJD(1) + t570 + t590) * t39) * m(7) + (-t67 * (t155 + t157) * t469 - g(1) * (t155 + t350) + (-g(2) + t29) * (t156 + t365) + (t269 * t388 + t260 + t354) * t28 + (-t251 * t388 + t528 + t571 - t581) * t65 + (t286 * t469 - t428 + t65 + (-t155 + t580) * qJD(1) + t570 + t93) * t66) * m(6) + ((t25 + (t269 * t95 + t270 * t94) * t339 + t569 + t594) * t199 + t11 + t589) * t545 + (t121 + t78) * t549 + (m(2) * (t290 ^ 2 + t295 ^ 2) + t390 + Icges(2,3) + Icges(3,1) + Icges(4,2) + Icges(5,3)) * qJDD(1) + ((-g(1) + t60) * (-t205 + t355) + (-g(2) + t61) * (t427 + t204) + (t418 + t584) * t133 + (t133 - t428 + t429 + t480 + (t205 + t580) * qJD(1)) * ((t204 + t367) * qJD(1) + t419)) * m(5) - t568 * t551 + (t111 + ((-t53 + t54 - t492) * t270 + t491 * t269) * qJD(5)) * t432; (-m(3) + t451) * (g(1) * t534 - g(2) * t342) + 0.2e1 * (t12 * t536 + t13 * t456) * m(7) + 0.2e1 * (t28 * t456 + t29 * t536) * m(6) + 0.2e1 * (t456 * t60 + t536 * t61) * m(5) + 0.2e1 * (t108 * t456 + t109 * t536) * m(4) + 0.2e1 * (t128 * t456 + t129 * t536) * m(3); t451 * (g(1) * t342 + g(2) * t534) + m(4) * (t108 * t342 + t109 * t534) + m(5) * (t60 * t342 + t534 * t61) + m(6) * (t28 * t342 + t29 * t534) + m(7) * (t12 * t534 + t13 * t342); m(5) * qJDD(4) + m(6) * t23 + m(7) * t7 + g(3) * t459; -(qJD(1) * t50 + qJD(5) * t415 + qJDD(1) * t121 + t195 * t55 + t196 * t54 + t1) * t270 / 0.2e1 - (qJD(1) * t51 + qJD(5) * t414 + qJDD(1) * t120 + t195 * t53 + t196 * t52 + t2) * t269 / 0.2e1 - (t22 + t11) * t252 / 0.2e1 + (t21 + t10) * t251 / 0.2e1 + ((t30 * t494 - t38 * t481) * t252 + (t30 * t493 - t39 * t481) * t251 + (-t13 * t481 + t30 * t523 - t38 * t486 + t494 * t7) * t270 + (t12 * t481 + t39 * t486 - t7 * t493 + t30 * (-t112 - t48)) * t269 - t38 * (-qJD(1) * t193 - t145 * t302 - t197 * t225 - t270 * t577) - t39 * (qJD(1) * t191 + t144 * t302 - t199 * t225 + t269 * t577) - t30 * (t144 * t197 + t145 * t199 + t191 * t469 - t193 * t470) - ((t104 * t39 - t105 * t38) * t339 + t347 * t341) * qJD(6) - g(1) * (t247 + t484) - g(2) * (-t249 + t483) - g(3) * (t225 + t562) - (-g(1) * t270 + g(2) * t269) * t339 * (pkin(5) + t530)) * m(7) + ((-t140 * t177 - t142 * t178) * t199 - (t141 * t177 + t143 * t178) * t197 + (t177 * t221 + t178 * t223) * t302 + (-t25 * t503 + t339 * t58) * qJD(6) + ((qJD(6) * t24 + t364) * t341 + t586) * t270) * t545 - t14 * t461 / 0.2e1 + (-t24 * t269 - t25 * t270) * t550 + (-t26 * t269 - t27 * t270) * t551 + t440 * t552 + t439 * t553 + (-t269 * t33 - t270 * t34) * t541 + (t24 * t251 - t25 * t252 - t269 * t5 - t270 * t6) * t544 + (t251 * t26 - t252 * t27 - t269 * t3 - t270 * t4) * t546 + t410 * t548 + t409 * t549 + (t251 * t33 - t252 * t34 - t269 * t8 - t270 * t9) * t537 + (t251 * t54 - t252 * t55 + t415) * t433 - qJD(1) * ((-t479 * t339 + t478 * t341) * qJD(1) + ((-t269 * t488 - t270 * t487) * t341 + (t269 * t490 - t270 * t489) * t339) * qJD(5)) / 0.2e1 + ((-t470 * t502 - t472) * t269 + (-t575 + (-t358 * t270 + (t182 + t359) * t269) * qJD(5)) * t270) * t434 + ((t182 * t469 - t472) * t270 + (t575 + (-t359 * t269 + (-t502 + t358) * t270) * qJD(5)) * t269) * t432 + (t251 * t52 - t252 * t53 + t414) * t435 + qJD(1) * (t251 * t77 - t252 * t78 - t269 * t31 - t270 * t32) / 0.2e1 + qJDD(1) * (-t269 * t77 - t270 * t78) / 0.2e1 + (g(1) * t187 - g(2) * t188 - g(3) * t291 - (-t187 * t66 - t188 * t65) * qJD(1) - (t67 * (-t187 * t270 - t188 * t269) + t408 * t291) * qJD(5) + t23 * t393 + t67 * (t156 * t252 + t157 * t251 + t407) + t408 * t274 + (-t251 * t66 - t252 * t65 + t269 * t29 - t270 * t28) * t286) * m(6) + ((-t140 * t179 + t142 * t176) * t199 - (t141 * t179 - t143 * t176) * t197 + (-t176 * t223 + t179 * t221) * t302 + (t26 * t500 - t339 * t568) * qJD(6) + ((-qJD(6) * t27 - t364) * t341 - t586) * t269) * t547 + (((t140 * t338 - t142 * t340 + t95) * t199 - (-t141 * t338 + t143 * t340 - t94) * t197 + (-t221 * t338 + t223 * t340 + t218) * t302 + t124 * qJD(6)) * t339 + (qJD(6) * t411 - t555) * t341) * t538; t2 * t501 / 0.2e1 + (t339 * t413 - t341 * t58) * t550 + ((qJD(5) * t413 - t20) * t341 + (qJD(5) * t58 + t24 * t252 + t25 * t251 - t269 * t6 + t270 * t5) * t339) * t544 + t507 * t552 + t341 * t11 * t435 - t1 * t504 / 0.2e1 + (t339 * t412 + t341 * t568) * t551 + ((qJD(5) * t412 - t19) * t341 + (-qJD(5) * t568 + t251 * t27 + t252 * t26 - t269 * t4 + t270 * t3) * t339) * t546 + t14 * t463 / 0.2e1 - t341 * (t524 + t526 + t527 - t531 + t532) / 0.2e1 + (-t124 * t341 + t339 * t411) * t541 + ((qJD(5) * t411 - t41) * t341 + (qJD(5) * t124 + t251 * t34 + t252 * t33 - t269 * t9 + t270 * t8) * t339) * t537 + (-t177 * t554 + t178 * t352 + t362 * t501) * t545 + (-t176 * t352 - t179 * t554 - t362 * t504) * t547 + (-t362 * t341 + (t338 * t554 + t352 * t340) * t339) * t538 + (t506 / 0.2e1 + t341 * t432) * t10 + ((qJD(5) * t347 - t12 * t104 + t13 * t105 + t38 * t48 - t39 * t49) * t341 + (t38 * (-qJD(5) * t105 - t165 * t269) + t39 * (qJD(5) * t104 - t165 * t270) + t7 * t396 + t30 * (-t104 * t251 + t105 * t252 + t269 * t49 + t270 * t48) + (-t12 * t270 - t13 * t269 + t251 * t38 - t252 * t39) * t224) * t339 - t38 * (-t123 * t302 - t197 * t267) - t39 * (t122 * t302 - t199 * t267) - t30 * (t122 * t197 + t123 * t199) - g(1) * t122 - g(2) * t123 - g(3) * t267) * m(7);];
tau  = t15;
