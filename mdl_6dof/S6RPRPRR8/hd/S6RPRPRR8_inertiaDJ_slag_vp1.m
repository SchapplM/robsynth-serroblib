% Calculate time derivative of joint inertia matrix for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR8_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:46
% EndTime: 2019-03-09 03:58:15
% DurationCPUTime: 17.18s
% Computational Cost: add. (36124->969), mult. (40544->1336), div. (0->0), fcn. (38570->10), ass. (0->478)
t356 = qJ(3) + pkin(10);
t344 = cos(t356);
t366 = cos(qJ(1));
t343 = sin(t356);
t510 = qJD(3) * t343;
t462 = -t510 / 0.2e1;
t363 = sin(qJ(1));
t512 = qJD(1) * t363;
t463 = -t512 / 0.2e1;
t628 = t344 * t463 + t366 * t462;
t461 = t510 / 0.2e1;
t511 = qJD(1) * t366;
t476 = t344 * t511;
t627 = t363 * t461 - t476 / 0.2e1;
t362 = sin(qJ(3));
t365 = cos(qJ(3));
t505 = qJD(3) * t365;
t626 = t362 * t511 + t363 * t505;
t507 = qJD(3) * t363;
t471 = t344 * t507;
t625 = t343 * t511 + t471;
t504 = qJD(3) * t366;
t472 = t343 * t504;
t477 = t344 * t512;
t624 = t472 + t477;
t473 = t343 * t507;
t376 = t473 - t476;
t359 = qJ(5) + qJ(6);
t348 = sin(t359);
t349 = cos(t359);
t574 = Icges(7,4) * t349;
t414 = -Icges(7,2) * t348 + t574;
t233 = Icges(7,6) * t343 + t344 * t414;
t575 = Icges(7,4) * t348;
t420 = Icges(7,1) * t349 - t575;
t234 = Icges(7,5) * t343 + t344 * t420;
t623 = -t233 * t348 + t234 * t349;
t364 = cos(qJ(5));
t342 = pkin(5) * t364 + pkin(4);
t593 = pkin(4) - t342;
t622 = t343 * t593;
t581 = Icges(4,4) * t362;
t418 = Icges(4,2) * t365 + t581;
t263 = Icges(4,6) * t366 + t363 * t418;
t580 = Icges(4,4) * t365;
t424 = Icges(4,1) * t362 + t580;
t265 = Icges(4,5) * t366 + t363 * t424;
t396 = t263 * t365 + t265 * t362;
t621 = t366 * t396;
t579 = Icges(5,4) * t343;
t416 = Icges(5,2) * t344 + t579;
t246 = Icges(5,6) * t366 + t363 * t416;
t578 = Icges(5,4) * t344;
t422 = Icges(5,1) * t343 + t578;
t248 = Icges(5,5) * t366 + t363 * t422;
t398 = t246 * t344 + t248 * t343;
t620 = t366 * t398;
t438 = rSges(5,1) * t343 + rSges(5,2) * t344;
t384 = t366 * t438;
t439 = rSges(4,1) * t362 + rSges(4,2) * t365;
t385 = t366 * t439;
t367 = -pkin(9) - pkin(8);
t592 = -pkin(8) - t367;
t229 = t343 * t592 - t344 * t593;
t432 = rSges(7,1) * t349 - rSges(7,2) * t348;
t235 = rSges(7,3) * t343 + t344 * t432;
t527 = t229 + t235;
t456 = t527 * t363;
t474 = t365 * t511;
t484 = rSges(4,1) * t626 + rSges(4,2) * t474;
t500 = -rSges(4,3) - pkin(1) - pkin(7);
t508 = qJD(3) * t362;
t517 = qJ(2) * t511 + qJD(2) * t363;
t166 = (-rSges(4,2) * t508 + qJD(1) * t500) * t363 + t484 + t517;
t591 = rSges(4,2) * t362;
t317 = rSges(4,1) * t365 - t591;
t347 = qJD(2) * t366;
t167 = t347 + t317 * t504 + (t500 * t366 + (-qJ(2) - t439) * t363) * qJD(1);
t619 = -t166 * t366 + t167 * t363;
t553 = t342 * t343;
t584 = rSges(7,3) - t367;
t618 = -t344 * t584 + t553;
t513 = qJD(1) * t343;
t449 = qJD(5) + t513;
t470 = t344 * t504;
t617 = t363 * t449 - t470;
t355 = qJD(5) + qJD(6);
t453 = t355 + t513;
t616 = t363 * t453 - t470;
t410 = Icges(5,5) * t343 + Icges(5,6) * t344;
t615 = -Icges(5,3) * t363 + t366 * t410;
t412 = Icges(4,5) * t362 + Icges(4,6) * t365;
t614 = -Icges(4,3) * t363 + t366 * t412;
t613 = -Icges(5,6) * t363 + t366 * t416;
t612 = -Icges(4,6) * t363 + t366 * t418;
t611 = -Icges(5,5) * t363 + t366 * t422;
t610 = -Icges(4,5) * t363 + t366 * t424;
t609 = t366 * t453 + t471;
t608 = 2 * m(4);
t607 = 2 * m(5);
t606 = 2 * m(6);
t605 = 2 * m(7);
t357 = t363 ^ 2;
t358 = t366 ^ 2;
t604 = t343 / 0.2e1;
t603 = -t363 / 0.2e1;
t602 = t363 / 0.2e1;
t601 = t366 / 0.2e1;
t600 = rSges(3,2) - pkin(1);
t599 = -rSges(5,3) - pkin(1);
t598 = m(4) * t317;
t597 = pkin(3) * t362;
t596 = pkin(3) * t365;
t595 = pkin(4) * t343;
t594 = pkin(7) * t366;
t354 = t366 * pkin(1);
t590 = rSges(5,2) * t343;
t589 = rSges(6,3) * t344;
t542 = t349 * t366;
t545 = t348 * t363;
t255 = -t343 * t545 + t542;
t543 = t349 * t363;
t544 = t348 * t366;
t256 = t343 * t543 + t544;
t548 = t344 * t363;
t171 = Icges(7,5) * t256 + Icges(7,6) * t255 - Icges(7,3) * t548;
t173 = Icges(7,4) * t256 + Icges(7,2) * t255 - Icges(7,6) * t548;
t175 = Icges(7,1) * t256 + Icges(7,4) * t255 - Icges(7,5) * t548;
t407 = t173 * t348 - t175 * t349;
t454 = t343 * t355 + qJD(1);
t393 = t454 * t363;
t156 = -t348 * t609 - t349 * t393;
t157 = -t348 * t393 + t349 * t609;
t95 = Icges(7,5) * t157 + Icges(7,6) * t156 + Icges(7,3) * t376;
t97 = Icges(7,4) * t157 + Icges(7,2) * t156 + Icges(7,6) * t376;
t99 = Icges(7,1) * t157 + Icges(7,4) * t156 + Icges(7,5) * t376;
t24 = (qJD(3) * t407 + t95) * t343 + (qJD(3) * t171 + (-t173 * t355 + t99) * t349 + (-t175 * t355 - t97) * t348) * t344;
t588 = t24 * t366;
t257 = t343 * t544 + t543;
t258 = -t343 * t542 + t545;
t547 = t344 * t366;
t172 = Icges(7,5) * t258 + Icges(7,6) * t257 + Icges(7,3) * t547;
t174 = Icges(7,4) * t258 + Icges(7,2) * t257 + Icges(7,6) * t547;
t176 = Icges(7,1) * t258 + Icges(7,4) * t257 + Icges(7,5) * t547;
t406 = t174 * t348 - t176 * t349;
t394 = t366 * t454;
t154 = -t348 * t616 + t349 * t394;
t155 = t348 * t394 + t349 * t616;
t94 = Icges(7,5) * t155 + Icges(7,6) * t154 - Icges(7,3) * t624;
t96 = Icges(7,4) * t155 + Icges(7,2) * t154 - Icges(7,6) * t624;
t98 = Icges(7,1) * t155 + Icges(7,4) * t154 - Icges(7,5) * t624;
t25 = (qJD(3) * t406 + t94) * t343 + (qJD(3) * t172 + (-t174 * t355 + t98) * t349 + (-t176 * t355 - t96) * t348) * t344;
t587 = t25 * t363;
t586 = t363 * rSges(4,3);
t585 = t363 * rSges(5,3);
t353 = t366 * rSges(4,3);
t352 = t366 * rSges(5,3);
t408 = Icges(7,5) * t349 - Icges(7,6) * t348;
t232 = Icges(7,3) * t343 + t344 * t408;
t124 = t232 * t343 + t344 * t623;
t549 = t344 * t355;
t162 = (-Icges(7,2) * t349 - t575) * t549 + (Icges(7,6) * t344 - t343 * t414) * qJD(3);
t161 = (-Icges(7,5) * t348 - Icges(7,6) * t349) * t549 + (Icges(7,3) * t344 - t343 * t408) * qJD(3);
t163 = (-Icges(7,1) * t348 - t574) * t549 + (Icges(7,5) * t344 - t343 * t420) * qJD(3);
t509 = qJD(3) * t344;
t373 = t344 * t349 * t163 + t343 * t161 + t232 * t509 - t510 * t623;
t496 = t355 * t349 * t233;
t583 = t124 * t509 + ((-t496 + (-t234 * t355 - t162) * t348) * t344 + t373) * t343;
t361 = sin(qJ(5));
t409 = Icges(6,5) * t364 - Icges(6,6) * t361;
t238 = Icges(6,3) * t343 + t344 * t409;
t577 = Icges(6,4) * t361;
t421 = Icges(6,1) * t364 - t577;
t240 = Icges(6,5) * t343 + t344 * t421;
t576 = Icges(6,4) * t364;
t415 = -Icges(6,2) * t361 + t576;
t239 = Icges(6,6) * t343 + t344 * t415;
t562 = t239 * t361;
t129 = t238 * t343 + (t240 * t364 - t562) * t344;
t503 = qJD(5) * t344;
t185 = (-Icges(6,5) * t361 - Icges(6,6) * t364) * t503 + (Icges(6,3) * t344 - t343 * t409) * qJD(3);
t187 = (-Icges(6,1) * t361 - t576) * t503 + (Icges(6,5) * t344 - t343 * t421) * qJD(3);
t506 = qJD(3) * t364;
t372 = t344 * t364 * t187 + t238 * t509 + t510 * t562 + (-t240 * t506 + t185) * t343;
t186 = (-Icges(6,2) * t364 - t577) * t503 + (Icges(6,6) * t344 - t343 * t415) * qJD(3);
t541 = t361 * t186;
t582 = t129 * t509 + ((-t541 + (-t239 * t364 - t240 * t361) * qJD(5)) * t344 + t372) * t343;
t537 = t363 * t364;
t539 = t361 * t366;
t281 = t343 * t539 + t537;
t535 = t364 * t366;
t540 = t361 * t363;
t282 = -t343 * t535 + t540;
t436 = -rSges(6,1) * t282 - rSges(6,2) * t281;
t209 = rSges(6,3) * t547 - t436;
t565 = t209 * t366;
t561 = t246 * t343;
t560 = t613 * t343;
t559 = t248 * t344;
t558 = t611 * t344;
t557 = t263 * t362;
t556 = t612 * t362;
t555 = t265 * t365;
t554 = t610 * t365;
t552 = t343 * t363;
t551 = t343 * t366;
t550 = t344 * t342;
t546 = t344 * t367;
t538 = t362 * t363;
t536 = t363 * t365;
t434 = t155 * rSges(7,1) + t154 * rSges(7,2);
t100 = -rSges(7,3) * t624 + t434;
t338 = pkin(5) * t539;
t502 = qJD(5) * t361;
t452 = pkin(5) * t343 * t502;
t485 = pkin(4) * t470 + pkin(8) * t624;
t501 = qJD(5) * t364;
t497 = pkin(5) * t501;
t534 = t100 + t363 * t497 + (t452 + (t343 * t367 - t550) * qJD(3)) * t366 + (t338 + (t546 - t622) * t363) * qJD(1) + t485;
t101 = t157 * rSges(7,1) + t156 * rSges(7,2) + rSges(7,3) * t376;
t450 = qJD(5) * t343 + qJD(1);
t391 = t450 * t361;
t370 = -pkin(5) * t391 - t367 * t510;
t486 = pkin(4) * t625 + pkin(8) * t473;
t374 = pkin(8) * t476 - t486;
t447 = t342 * t625 + t366 * t497 + t367 * t476;
t120 = t363 * t370 + t374 + t447;
t533 = -t101 - t120;
t164 = (-rSges(7,1) * t348 - rSges(7,2) * t349) * t549 + (rSges(7,3) * t344 - t343 * t432) * qJD(3);
t464 = t592 * t344;
t468 = t344 * t502;
t201 = -pkin(5) * t468 + (t464 + t622) * qJD(3);
t532 = -t164 - t201;
t177 = t256 * rSges(7,1) + t255 * rSges(7,2) - rSges(7,3) * t548;
t320 = pkin(4) * t552;
t271 = -pkin(8) * t548 + t320;
t488 = t342 * t552 + t363 * t546 + t338;
t193 = -t271 + t488;
t531 = -t177 - t193;
t433 = -t258 * rSges(7,1) - t257 * rSges(7,2);
t178 = rSges(7,3) * t547 - t433;
t321 = pkin(4) * t551;
t194 = pkin(5) * t540 + t321 + (t464 - t553) * t366;
t530 = t178 + t194;
t336 = pkin(3) * t538;
t360 = -qJ(4) - pkin(7);
t498 = pkin(3) * t504;
t519 = -t360 * t511 - t365 * t498;
t440 = qJD(4) * t363 + t519;
t222 = t366 * ((t336 - t594) * qJD(1) + t440);
t529 = t366 * (qJD(1) * t320 - t485) + t222;
t446 = pkin(3) * t626 + qJD(4) * t366 + t360 * t512;
t227 = pkin(7) * t512 + t446;
t528 = t374 - t227;
t130 = t343 * t177 + t235 * t548;
t250 = rSges(5,1) * t552 + rSges(5,2) * t548 + t352;
t278 = t336 + (-pkin(7) - t360) * t366;
t526 = -t250 - t278;
t518 = t363 * t360 + t366 * t597;
t277 = -pkin(7) * t363 - t518;
t259 = t366 * t277;
t272 = pkin(8) * t547 - t321;
t525 = t366 * t272 + t259;
t279 = -t343 * t540 + t535;
t280 = t343 * t537 + t539;
t524 = t280 * rSges(6,1) + t279 * rSges(6,2);
t523 = -t271 - t278;
t522 = -t272 - t277;
t296 = t344 * pkin(4) + t343 * pkin(8);
t337 = pkin(3) * t536;
t521 = t363 * t296 + t337;
t520 = qJD(1) * t337 + t362 * t498;
t516 = t363 * qJ(2) + t354;
t244 = Icges(5,3) * t366 + t363 * t410;
t515 = qJD(1) * t244;
t261 = Icges(4,3) * t366 + t363 * t412;
t514 = qJD(1) * t261;
t499 = pkin(3) * t508;
t121 = -t238 * t548 + t239 * t279 + t240 * t280;
t202 = Icges(6,5) * t280 + Icges(6,6) * t279 - Icges(6,3) * t548;
t204 = Icges(6,4) * t280 + Icges(6,2) * t279 - Icges(6,6) * t548;
t206 = Icges(6,1) * t280 + Icges(6,4) * t279 - Icges(6,5) * t548;
t403 = t204 * t361 - t206 * t364;
t91 = t202 * t343 - t344 * t403;
t495 = -t91 / 0.2e1 - t121 / 0.2e1;
t122 = t238 * t547 + t239 * t281 + t240 * t282;
t203 = Icges(6,5) * t282 + Icges(6,6) * t281 + Icges(6,3) * t547;
t205 = Icges(6,4) * t282 + Icges(6,2) * t281 + Icges(6,6) * t547;
t207 = Icges(6,1) * t282 + Icges(6,4) * t281 + Icges(6,5) * t547;
t402 = t205 * t361 - t207 * t364;
t92 = t203 * t343 - t344 * t402;
t494 = -t92 / 0.2e1 - t122 / 0.2e1;
t493 = t177 * t624 + t178 * t473;
t190 = -t450 * t537 + (-t366 * t449 - t471) * t361;
t191 = t449 * t535 + (t344 * t506 - t391) * t363;
t492 = t191 * rSges(6,1) + t190 * rSges(6,2) + rSges(6,3) * t473;
t208 = -rSges(6,3) * t548 + t524;
t491 = -t208 + t523;
t288 = (pkin(8) * t344 - t595) * qJD(3);
t331 = pkin(3) * t474;
t490 = t363 * t288 + t296 * t511 + t331;
t489 = t296 * t512 + t520;
t487 = -rSges(5,1) * t625 - rSges(5,2) * t476;
t269 = rSges(4,1) * t538 + rSges(4,2) * t536 + t353;
t351 = t366 * qJ(2);
t483 = t351 + t518;
t482 = -pkin(5) * t361 - pkin(1);
t481 = t344 * (-rSges(6,3) - pkin(8));
t435 = rSges(6,1) * t364 - rSges(6,2) * t361;
t241 = rSges(6,3) * t343 + t344 * t435;
t480 = t241 * t512;
t479 = t241 * t511;
t467 = -t548 / 0.2e1;
t466 = t547 / 0.2e1;
t106 = -t232 * t548 + t233 * t255 + t234 * t256;
t72 = -t171 * t548 + t173 * t255 + t175 * t256;
t73 = -t172 * t548 + t174 * t255 + t176 * t256;
t431 = t363 * t72 - t366 * t73;
t38 = t106 * t343 - t344 * t431;
t107 = t232 * t547 + t233 * t257 + t234 * t258;
t18 = t154 * t173 + t155 * t175 - t171 * t624 + t257 * t97 + t258 * t99 + t547 * t95;
t19 = t154 * t174 + t155 * t176 - t172 * t624 + t257 * t96 + t258 * t98 + t547 * t94;
t74 = t171 * t547 + t173 * t257 + t175 * t258;
t75 = t172 * t547 + t174 * t257 + t176 * t258;
t430 = t363 * t74 - t366 * t75;
t46 = t154 * t233 + t155 * t234 + t161 * t547 + t162 * t257 + t163 * t258 - t232 * t624;
t54 = t363 * t75 + t366 * t74;
t4 = (qJD(3) * t430 + t46) * t343 + (-qJD(1) * t54 + qJD(3) * t107 - t18 * t363 + t19 * t366) * t344;
t84 = t171 * t343 - t344 * t407;
t85 = t172 * t343 - t344 * t406;
t428 = t363 * t84 - t366 * t85;
t429 = t363 * t85 + t366 * t84;
t465 = t4 * t547 + t38 * t473 + (t124 * t343 - t344 * t428) * t509 + t343 * (t428 * t510 + (-qJD(1) * t429 - t24 * t363 + t25 * t366) * t344 + t583);
t460 = -qJ(2) - t597;
t459 = -t296 - t596;
t86 = -t202 * t548 + t204 * t279 + t206 * t280;
t87 = -t203 * t548 + t205 * t279 + t207 * t280;
t427 = t363 * t86 - t366 * t87;
t42 = t343 * t121 - t344 * t427;
t458 = t343 * t91 + t42;
t300 = t439 * qJD(3);
t457 = t300 * (t357 + t358);
t455 = qJD(1) * t527;
t451 = t343 * t101 + t164 * t548 + t177 * t509 + t235 * t476;
t448 = t523 + t531;
t39 = t107 * t343 - t344 * t430;
t88 = t202 * t547 + t204 * t281 + t206 * t282;
t89 = t203 * t547 + t205 * t281 + t207 * t282;
t426 = t363 * t88 - t366 * t89;
t43 = t343 * t122 - t344 * t426;
t445 = -t343 * t92 - t39 - t43;
t294 = rSges(5,1) * t344 - t590;
t392 = t366 * t450;
t188 = -t361 * t617 + t364 * t392;
t189 = t361 * t392 + t364 * t617;
t437 = rSges(6,1) * t189 + rSges(6,2) * t188;
t53 = t363 * t73 + t366 * t72;
t59 = t363 * t87 + t366 * t86;
t60 = t363 * t89 + t366 * t88;
t425 = Icges(4,1) * t365 - t581;
t423 = Icges(5,1) * t344 - t579;
t419 = -Icges(4,2) * t362 + t580;
t417 = -Icges(5,2) * t343 + t578;
t413 = Icges(4,5) * t365 - Icges(4,6) * t362;
t411 = Icges(5,5) * t344 - Icges(5,6) * t343;
t401 = t208 * t366 + t209 * t363;
t397 = -t343 * t611 - t344 * t613;
t395 = -t362 * t610 - t365 * t612;
t390 = -t360 * t366 + t336 + t516;
t389 = t347 - t440;
t387 = t446 + t517;
t386 = rSges(3,3) * t366 + t363 * t600;
t383 = t397 * t363;
t382 = t395 * t363;
t381 = qJD(3) * t425;
t380 = qJD(3) * t423;
t379 = qJD(3) * t419;
t378 = qJD(3) * t417;
t377 = -t363 * pkin(1) + t366 * t481;
t11 = -qJD(1) * t430 + t18 * t366 + t19 * t363;
t20 = t156 * t173 + t157 * t175 + t171 * t376 + t255 * t97 + t256 * t99 - t548 * t95;
t21 = t156 * t174 + t157 * t176 + t172 * t376 + t255 * t96 + t256 * t98 - t548 * t94;
t12 = -qJD(1) * t431 + t20 * t366 + t21 * t363;
t47 = t156 * t233 + t157 * t234 - t161 * t548 + t162 * t255 + t163 * t256 + t232 * t376;
t5 = (qJD(3) * t431 + t47) * t343 + (-qJD(1) * t53 + qJD(3) * t106 - t20 * t363 + t21 * t366) * t344;
t371 = t11 * t466 + t12 * t467 + t4 * t602 + t5 * t601 + (-qJD(1) * t428 + t587 + t588) * t604 + t38 * t463 + t39 * t511 / 0.2e1 + t429 * t509 / 0.2e1 + t628 * t54 + t627 * t53;
t369 = t583 + (t24 + t47) * t467 + (t25 + t46) * t466 + (t107 + t85) * t628 + (t106 + t84) * t627;
t368 = (-t363 * t5 + (-t363 * t39 - t366 * t38) * qJD(1)) * t344 - t39 * t472 + t465;
t287 = t438 * qJD(3);
t276 = -rSges(3,2) * t366 + rSges(3,3) * t363 + t516;
t275 = t351 + t386;
t270 = t586 - t385;
t251 = t585 - t384;
t243 = (-t294 - t596) * t366;
t242 = t294 * t363 + t337;
t237 = t347 + (t600 * t366 + (-rSges(3,3) - qJ(2)) * t363) * qJD(1);
t236 = qJD(1) * t386 + t517;
t231 = t269 + t516 + t594;
t230 = t363 * t500 + t351 + t385;
t224 = t235 * t547;
t216 = qJD(1) * t614 + t413 * t507;
t215 = -t413 * t504 + t514;
t212 = t390 + t250;
t211 = t363 * t599 + t384 + t483;
t196 = qJD(1) * t615 + t411 * t507;
t195 = -t411 * t504 + t515;
t192 = (-rSges(6,1) * t361 - rSges(6,2) * t364) * t503 + (-t343 * t435 + t589) * qJD(3);
t184 = (-t241 + t459) * t366;
t183 = t241 * t363 + t521;
t180 = t294 * t511 + t331 + (-t287 - t499) * t363;
t179 = t287 * t366 + t294 * t512 + t520;
t150 = t164 * t547;
t147 = -t363 * t614 - t366 * t395;
t146 = t261 * t363 - t621;
t145 = -t366 * t614 + t382;
t144 = t261 * t366 + t363 * t396;
t143 = t363 * t481 + t320 + t390 + t524;
t142 = t321 + t377 + t436 + t483;
t141 = -t209 * t343 + t241 * t547;
t140 = t208 * t343 + t241 * t548;
t139 = -t363 * t615 - t366 * t397;
t138 = t244 * t363 - t620;
t137 = -t366 * t615 + t383;
t136 = t244 * t366 + t363 * t398;
t135 = (t459 - t527) * t366;
t134 = t521 + t456;
t133 = t294 * t504 + (t599 * t366 + (-t438 + t460) * t363) * qJD(1) + t389;
t132 = (-rSges(5,2) * t510 + qJD(1) * t599) * t363 + t387 - t487;
t131 = -t178 * t343 + t224;
t127 = t177 + t390 + t488;
t126 = t482 * t363 + t366 * t618 + t433 + t483;
t125 = t401 * t344;
t118 = (-t177 * t366 - t178 * t363) * t344;
t117 = -rSges(6,3) * t476 + t492;
t116 = -rSges(6,3) * t624 + t437;
t115 = Icges(6,1) * t191 + Icges(6,4) * t190 + Icges(6,5) * t376;
t114 = Icges(6,1) * t189 + Icges(6,4) * t188 - Icges(6,5) * t624;
t113 = Icges(6,4) * t191 + Icges(6,2) * t190 + Icges(6,6) * t376;
t112 = Icges(6,4) * t189 + Icges(6,2) * t188 - Icges(6,6) * t624;
t111 = Icges(6,5) * t191 + Icges(6,6) * t190 + Icges(6,3) * t376;
t110 = Icges(6,5) * t189 + Icges(6,6) * t188 - Icges(6,3) * t624;
t109 = t479 + (t192 - t499) * t363 + t490;
t108 = t480 + (-t192 - t288) * t366 + t489;
t90 = t363 * t491 + t525 + t565;
t83 = t229 * t547 - t343 * t530 + t224;
t82 = t193 * t343 + t229 * t548 + t130;
t77 = rSges(6,3) * t472 + (-t354 + (t460 + t589 - t595) * t363) * qJD(1) + t389 - t437 + t485;
t76 = qJD(1) * t377 + t387 + t486 + t492;
t71 = (-t363 * t530 + t366 * t531) * t344;
t70 = t366 * t455 + (-t499 - t532) * t363 + t490;
t69 = t363 * t455 + (-t288 + t532) * t366 + t489;
t68 = t347 + (-qJD(4) - t497) * t363 + (-t452 + (t343 * t584 + t550) * qJD(3)) * t366 + (t482 * t366 + (t460 - t618) * t363) * qJD(1) - t434 - t519;
t67 = (-pkin(1) * qJD(1) + t370) * t363 + t101 + t387 + t447;
t66 = t363 * t448 + t366 * t530 + t525;
t65 = (-t241 * t507 + t117) * t343 + (qJD(3) * t208 + t192 * t363 + t479) * t344;
t64 = (-t241 * t504 - t116) * t343 + (-qJD(3) * t209 + t192 * t366 - t480) * t344;
t62 = -t235 * t473 + t451;
t61 = -t235 * t477 - t100 * t343 + t150 + (-t178 * t344 - t235 * t551) * qJD(3);
t57 = -t185 * t548 + t186 * t279 + t187 * t280 + t190 * t239 + t191 * t240 + t238 * t376;
t56 = t185 * t547 + t186 * t281 + t187 * t282 + t188 * t239 + t189 * t240 - t238 * t624;
t48 = t401 * t510 + (-t116 * t363 - t117 * t366 + (t208 * t363 - t565) * qJD(1)) * t344;
t40 = (-t100 * t363 + (-qJD(1) * t178 - t101) * t366) * t344 + t493;
t34 = t116 * t366 + (-t117 + t528) * t363 + (t491 * t366 + (-t209 + t522) * t363) * qJD(1) + t529;
t33 = (qJD(3) * t402 + t110) * t343 + (qJD(3) * t203 - t112 * t361 + t114 * t364 + (-t205 * t364 - t207 * t361) * qJD(5)) * t344;
t32 = (qJD(3) * t403 + t111) * t343 + (qJD(3) * t202 - t113 * t361 + t115 * t364 + (-t204 * t364 - t206 * t361) * qJD(5)) * t344;
t31 = t120 * t343 + (t201 * t363 + t229 * t511) * t344 + (t193 * t344 - t343 * t456) * qJD(3) + t451;
t30 = t150 + (-t504 * t527 - t534) * t343 + (-qJD(1) * t456 - qJD(3) * t530 + t201 * t366) * t344;
t29 = -t110 * t548 + t112 * t279 + t114 * t280 + t190 * t205 + t191 * t207 + t203 * t376;
t28 = -t111 * t548 + t113 * t279 + t115 * t280 + t190 * t204 + t191 * t206 + t202 * t376;
t27 = t110 * t547 + t112 * t281 + t114 * t282 + t188 * t205 + t189 * t207 - t203 * t624;
t26 = t111 * t547 + t113 * t281 + t115 * t282 + t188 * t204 + t189 * t206 - t202 * t624;
t17 = t534 * t366 + (t528 + t533) * t363 + (t448 * t366 + (t522 - t530) * t363) * qJD(1) + t529;
t16 = (t193 * t366 + t194 * t363) * t510 + ((qJD(1) * t193 - t534) * t363 + (-qJD(1) * t530 + t533) * t366) * t344 + t493;
t15 = -qJD(1) * t427 + t28 * t366 + t29 * t363;
t14 = -qJD(1) * t426 + t26 * t366 + t27 * t363;
t8 = (qJD(3) * t427 + t57) * t343 + (-qJD(1) * t59 + qJD(3) * t121 - t28 * t363 + t29 * t366) * t344;
t7 = (qJD(3) * t426 + t56) * t343 + (-qJD(1) * t60 + qJD(3) * t122 - t26 * t363 + t27 * t366) * t344;
t1 = [-t362 * t381 - t365 * t379 - t424 * t505 + t418 * t508 - t343 * t380 - t422 * t509 + t416 * t510 - t240 * t468 - t348 * t234 * t549 + 0.2e1 * m(3) * (t236 * t276 + t237 * t275) + (t166 * t231 + t167 * t230) * t608 + (t132 * t212 + t133 * t211) * t607 + (t142 * t77 + t143 * t76) * t606 + (t126 * t68 + t127 * t67) * t605 + t372 + t373 + (-t162 * t348 - t239 * t501 - t378 - t496 - t541) * t344; m(7) * (t363 * t68 - t366 * t67 + (t126 * t366 + t127 * t363) * qJD(1)) + m(6) * (t363 * t77 - t366 * t76 + (t142 * t366 + t143 * t363) * qJD(1)) + m(5) * (-t132 * t366 + t133 * t363 + (t211 * t366 + t212 * t363) * qJD(1)) + m(4) * ((t230 * t366 + t231 * t363) * qJD(1) + t619) + m(3) * (-t236 * t366 + t237 * t363 + (t275 * t366 + t276 * t363) * qJD(1)); 0; t588 / 0.2e1 + t587 / 0.2e1 + m(4) * (t619 * t317 - (t230 * t363 - t231 * t366) * t300) + m(5) * (t132 * t243 + t133 * t242 + t179 * t212 + t180 * t211) + m(6) * (t108 * t143 + t109 * t142 + t183 * t77 + t184 * t76) + m(7) * (t126 * t70 + t127 * t69 + t134 * t68 + t135 * t67) + ((t231 * t598 + t557 / 0.2e1 - t555 / 0.2e1 - t106 / 0.2e1 + t561 / 0.2e1 - t559 / 0.2e1 - t84 / 0.2e1 + t495) * t363 + (t107 / 0.2e1 + t556 / 0.2e1 - t554 / 0.2e1 + t230 * t598 + t560 / 0.2e1 - t558 / 0.2e1 + t85 / 0.2e1 - t494) * t366) * qJD(1) + (-t412 - t410) * qJD(3) * (t357 / 0.2e1 + t358 / 0.2e1) + (-(qJD(1) * t246 - t417 * t504) * t343 + (qJD(1) * t248 - t423 * t504) * t344 - (qJD(1) * t263 - t419 * t504) * t362 + (qJD(1) * t265 - t425 * t504) * t365 + t33 + t46 + t56 + (-t395 - t397) * qJD(3)) * t602 + (-(qJD(1) * t613 + t363 * t378) * t343 + (qJD(1) * t611 + t363 * t380) * t344 - (qJD(1) * t612 + t363 * t379) * t362 + (qJD(1) * t610 + t363 * t381) * t365 + t32 + t47 + t57 + (-t396 - t398) * qJD(3)) * t601; m(5) * (-t179 * t366 + t180 * t363 + (t242 * t366 + t243 * t363) * qJD(1)) + m(6) * (-t108 * t366 + t109 * t363 + (t183 * t366 + t184 * t363) * qJD(1)) + m(7) * (t363 * t70 - t366 * t69 + (t134 * t366 + t135 * t363) * qJD(1)) - m(4) * t457; (t134 * t70 + t135 * t69 + t66 * t17) * t605 + t363 * t11 + t366 * t12 + (t108 * t184 + t109 * t183 + t34 * t90) * t606 + t363 * t14 + t366 * t15 + (t242 * t180 + t243 * t179 + (t251 * t366 + t363 * t526 + t259) * (t222 + (-t227 + t487) * t363 + (-t294 * t358 + t357 * t590) * qJD(3) + ((t526 + t352) * t366 + (-t251 - t277 + t384 + t585) * t363) * qJD(1))) * t607 + t363 * ((t363 * t195 + (-t138 + t383) * qJD(1)) * t363 + (t139 * qJD(1) + (t246 * t510 - t248 * t509 + t515) * t366 + (t196 + (t558 - t560) * qJD(3) + t398 * qJD(1)) * t363) * t366) + t363 * ((t363 * t215 + (-t146 + t382) * qJD(1)) * t363 + (t147 * qJD(1) + (t263 * t508 - t265 * t505 + t514) * t366 + (t216 + (t554 - t556) * qJD(3) + t396 * qJD(1)) * t363) * t366) + t366 * ((t366 * t196 + (t137 + t620) * qJD(1)) * t366 + (-t136 * qJD(1) + (-t509 * t611 + t510 * t613) * t363 + (t195 + (t559 - t561) * qJD(3) + (-t244 + t397) * qJD(1)) * t366) * t363) + ((-t269 * t363 + t270 * t366) * (-t363 * t484 + (-t317 * t358 + t357 * t591) * qJD(3) + ((-t269 + t353) * t366 + (-t270 + t385 + t586) * t363) * qJD(1)) - t317 * t457) * t608 + t366 * ((t366 * t216 + (t145 + t621) * qJD(1)) * t366 + (-t144 * qJD(1) + (-t505 * t610 + t508 * t612) * t363 + (t215 + (t555 - t557) * qJD(3) + (-t261 + t395) * qJD(1)) * t366) * t363) + (-t53 - t59 + (-t136 - t144) * t366 + (-t137 - t145) * t363) * t512 + (t54 + t60 + (t138 + t146) * t366 + (t139 + t147) * t363) * t511; m(7) * (t363 * t67 + t366 * t68 + (-t126 * t363 + t127 * t366) * qJD(1)) + m(6) * (t363 * t76 + t366 * t77 + (-t142 * t363 + t143 * t366) * qJD(1)) + m(5) * (t132 * t363 + t133 * t366 + (-t211 * t363 + t212 * t366) * qJD(1)); 0; m(7) * (t363 * t69 + t366 * t70 + (-t134 * t363 + t135 * t366) * qJD(1)) + m(6) * (t108 * t363 + t109 * t366 + (-t183 * t363 + t184 * t366) * qJD(1)) + m(5) * (t179 * t363 + t180 * t366 + (-t242 * t363 + t243 * t366) * qJD(1)); 0; ((t56 / 0.2e1 + t33 / 0.2e1) * t366 + (-t57 / 0.2e1 - t32 / 0.2e1) * t363 + (t363 * t494 + t366 * t495) * qJD(1)) * t344 + (-t363 * t495 + t366 * t494) * t510 + m(6) * (t140 * t76 + t141 * t77 + t142 * t64 + t143 * t65) + m(7) * (t126 * t30 + t127 * t31 + t67 * t82 + t68 * t83) + t369 + t582; m(6) * (t363 * t64 - t366 * t65 + (t140 * t363 + t141 * t366) * qJD(1)) + m(7) * (t30 * t363 - t31 * t366 + (t363 * t82 + t366 * t83) * qJD(1)); (-qJD(1) * t42 / 0.2e1 + (-qJD(1) * t91 + t33) * t604 + t7 / 0.2e1 + t59 * t461) * t363 + m(6) * (t108 * t140 + t109 * t141 - t125 * t34 + t183 * t64 + t184 * t65 + t48 * t90) + m(7) * (t134 * t30 + t135 * t31 + t16 * t66 + t17 * t71 + t69 * t82 + t70 * t83) + (qJD(3) * (t363 * t92 + t366 * t91) / 0.2e1 + t14 * t601 + t15 * t603 + (-t366 * t59 / 0.2e1 + t60 * t603) * qJD(1)) * t344 + t371 + (qJD(1) * t43 / 0.2e1 + t60 * t462 + (qJD(1) * t92 + t32) * t604 + t8 / 0.2e1) * t366; m(6) * (t363 * t65 + t366 * t64 + (t140 * t366 - t141 * t363) * qJD(1)) + m(7) * (t30 * t366 + t31 * t363 + (-t363 * t83 + t366 * t82) * qJD(1)); (t16 * t71 + t30 * t83 + t31 * t82) * t605 + (-t125 * t48 + t140 * t65 + t141 * t64) * t606 + ((t363 * t458 + t366 * t445) * qJD(3) + t582) * t343 + ((t33 * t343 + t7) * t366 + (-t32 * t343 - t5 - t8) * t363 + (t129 * t343 + (-t363 * t91 + t366 * t92) * t344) * qJD(3) + ((-t38 - t458) * t366 + t445 * t363) * qJD(1)) * t344 + t465; m(7) * (t126 * t61 + t127 * t62 + t130 * t67 + t131 * t68) + t369; m(7) * (t363 * t61 - t366 * t62 + (t130 * t363 + t131 * t366) * qJD(1)); t371 + m(7) * (t118 * t17 + t130 * t69 + t131 * t70 + t134 * t61 + t135 * t62 + t40 * t66); m(7) * (t363 * t62 + t366 * t61 + (t130 * t366 - t131 * t363) * qJD(1)); m(7) * (t118 * t16 + t130 * t31 + t131 * t30 + t40 * t71 + t61 * t83 + t62 * t82) + t368; (t118 * t40 + t130 * t62 + t131 * t61) * t605 + t368;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;