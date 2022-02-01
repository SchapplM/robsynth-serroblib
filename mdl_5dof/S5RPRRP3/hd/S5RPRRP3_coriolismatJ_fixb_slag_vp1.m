% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:35
% EndTime: 2022-01-23 09:29:50
% DurationCPUTime: 9.32s
% Computational Cost: add. (29351->368), mult. (22022->496), div. (0->0), fcn. (20324->8), ass. (0->246)
t608 = Icges(5,4) + Icges(6,4);
t606 = Icges(5,1) + Icges(6,1);
t605 = Icges(5,5) + Icges(6,5);
t600 = Icges(5,6) + Icges(6,6);
t360 = qJ(3) + qJ(4);
t355 = sin(t360);
t609 = t608 * t355;
t604 = Icges(5,2) + Icges(6,2);
t356 = cos(t360);
t607 = -t600 * t355 + t605 * t356;
t586 = t606 * t356 - t609;
t603 = Icges(5,3) + Icges(6,3);
t602 = t608 * t356;
t359 = qJ(1) + pkin(8);
t353 = sin(t359);
t354 = cos(t359);
t595 = t605 * t353 + t586 * t354;
t599 = -t604 * t355 + t602;
t472 = t353 * t355;
t598 = t608 * t472;
t597 = t607 * t354;
t471 = t353 * t356;
t588 = -t603 * t354 + t605 * t471 - t600 * t472;
t596 = t603 * t353 + t597;
t574 = -t605 * t354 + t606 * t471 - t598;
t587 = t604 * t356 + t609;
t594 = t606 * t355 + t602;
t593 = t595 * t471;
t592 = -t600 * t354 + t608 * t471 - t472 * t604;
t591 = t600 * t353 + t599 * t354;
t351 = t353 ^ 2;
t352 = t354 ^ 2;
t415 = t351 + t352;
t498 = rSges(6,2) * t356;
t300 = rSges(6,1) * t355 + t498;
t417 = rSges(6,1) * t472 + rSges(6,2) * t471;
t424 = -t352 * t300 - t353 * t417;
t509 = pkin(4) * t355;
t126 = -t415 * t509 + t424;
t361 = sin(qJ(3));
t510 = pkin(3) * t361;
t317 = -t509 - t510;
t382 = t300 - t317;
t178 = t382 * t353;
t180 = t382 * t354;
t499 = rSges(6,1) * t356;
t508 = pkin(4) * t356;
t409 = rSges(6,2) * t355 - t499 - t508;
t216 = t409 * t353;
t218 = t409 * t354;
t442 = -t178 * t216 - t180 * t218;
t365 = -pkin(7) - pkin(6);
t334 = t353 * t365;
t347 = t354 * pkin(6);
t458 = t354 * t365;
t363 = cos(qJ(3));
t358 = t363 * pkin(3);
t350 = t358 + pkin(2);
t506 = pkin(2) - t350;
t435 = -t353 * (t353 * t506 - t347 - t458) + t354 * (-pkin(6) * t353 - t354 * t506 - t334);
t304 = t350 + t508;
t460 = t354 * t356;
t589 = rSges(6,3) + qJ(5) - t365;
t556 = -rSges(6,1) * t471 + rSges(6,2) * t472 - t353 * t304 + t589 * t354;
t461 = t354 * t355;
t557 = -rSges(6,2) * t461 + t589 * t353;
t83 = (rSges(6,1) * t460 + t334 + (t304 - t350) * t354 + t557) * t354 + (-t353 * t350 - t458 - t556) * t353;
t54 = t83 + t435;
t32 = t54 * t126 + t442;
t301 = rSges(5,1) * t355 + rSges(5,2) * t356;
t260 = t301 * t353;
t261 = t301 * t354;
t159 = -t353 * t260 - t354 * t261;
t500 = rSges(5,1) * t356;
t303 = -rSges(5,2) * t355 + t500;
t378 = t301 + t510;
t549 = t378 * t354;
t550 = t378 * t353;
t406 = -rSges(5,2) * t461 + rSges(5,3) * t353;
t418 = rSges(5,2) * t472 + t354 * rSges(5,3);
t144 = t353 * (rSges(5,1) * t471 - t418) + t354 * (rSges(5,1) * t460 + t406);
t88 = t144 + t435;
t44 = t88 * t159 + (t353 * t550 + t354 * t549) * t303;
t590 = -m(5) * t44 - m(6) * t32;
t585 = -t605 * t355 - t600 * t356;
t584 = -t596 * t354 + t593;
t583 = -t587 * t354 + t595;
t582 = -t471 * t604 + t574 - t598;
t581 = -t594 * t354 - t591;
t580 = t594 * t353 + t592;
t579 = -t588 * t353 - t574 * t460;
t570 = t596 * t353 + t595 * t460;
t578 = -t594 - t599;
t577 = -t591 * t472 + t584;
t576 = t592 * t461 + t579;
t575 = -t591 * t461 + t570;
t571 = t591 * t355 - t588;
t569 = t592 * t355;
t541 = m(5) / 0.2e1;
t540 = m(6) / 0.2e1;
t507 = sin(qJ(1)) * pkin(1);
t157 = -t507 + t556;
t511 = cos(qJ(1)) * pkin(1);
t158 = t511 + (t304 + t499) * t354 + t557;
t568 = m(6) * (t157 * t354 + t158 * t353);
t522 = t353 / 0.2e1;
t521 = -t354 / 0.2e1;
t567 = t354 / 0.2e1;
t486 = Icges(4,4) * t361;
t319 = Icges(4,2) * t363 + t486;
t322 = Icges(4,1) * t363 - t486;
t563 = (t322 / 0.2e1 - t319 / 0.2e1) * t361;
t562 = t585 * t354;
t561 = t585 * t353;
t560 = (t586 - t587) * t356 + t578 * t355;
t559 = -t583 * t355 + t581 * t356;
t558 = t582 * t355 + t580 * t356;
t182 = -t317 * t353 + t417;
t289 = t354 * t317;
t183 = -t300 * t354 + t289;
t410 = t300 + t509;
t215 = t410 * t353;
t217 = t410 * t354;
t407 = t350 + t500;
t160 = -t353 * t407 + t418 - t458 - t507;
t161 = t354 * t407 - t334 + t406 + t511;
t380 = (-t160 * t354 - t161 * t353) * t303;
t443 = t218 * t157 + t216 * t158;
t504 = (-t182 * t217 - t183 * t215 + t443) * t540 + (t380 + (t353 * t549 - t354 * t550) * t301) * t541;
t204 = pkin(4) * t472 + t417;
t205 = (-t498 + (-rSges(6,1) - pkin(4)) * t355) * t354;
t505 = (-t178 * t205 - t180 * t204 + t443) * t540 + (-t260 * t549 + t261 * t550 + t380) * t541;
t555 = t504 - t505;
t551 = qJD(1) * t568;
t357 = Icges(4,4) * t363;
t320 = -Icges(4,2) * t361 + t357;
t321 = Icges(4,1) * t361 + t357;
t377 = -t578 * t356 / 0.2e1 + (-t587 / 0.2e1 + t586 / 0.2e1) * t355;
t379 = (t575 * t353 + t576 * t354) * t567 + (t570 * t353 + ((t569 + t596) * t354 + t577 + t579 - t593) * t354) * t521 + ((t571 * t353 - t576 + t577 - t584) * t353 + ((t571 + t588) * t354 + (-t574 * t356 + t569) * t353 - t570 + t575) * t354) * t522;
t546 = 0.4e1 * qJD(1);
t545 = 2 * qJD(3);
t543 = 2 * qJD(4);
t542 = m(4) / 0.2e1;
t501 = rSges(4,1) * t363;
t412 = pkin(2) + t501;
t470 = t353 * t361;
t416 = rSges(4,2) * t470 + t354 * rSges(4,3);
t170 = -t353 * t412 + t347 + t416 - t507;
t459 = t354 * t361;
t331 = rSges(4,2) * t459;
t171 = t511 - t331 + t412 * t354 + (rSges(4,3) + pkin(6)) * t353;
t323 = rSges(4,1) * t361 + rSges(4,2) * t363;
t284 = t323 * t353;
t285 = t323 * t354;
t538 = m(4) * (t170 * t284 - t171 * t285);
t80 = t415 * t301 * t303 + t144 * t159;
t79 = m(5) * t80;
t532 = m(5) * (t160 * t550 - t161 * t549);
t531 = m(5) * (t160 * t260 - t161 * t261);
t441 = -t215 * t216 - t217 * t218;
t99 = t351 * (t317 + t510) + t354 * (pkin(3) * t459 + t289) + t424;
t528 = m(6) * (t83 * t99 + t441);
t525 = m(6) * (t157 * t182 + t158 * t183);
t524 = m(6) * (t157 * t204 + t158 * t205);
t523 = -t353 / 0.2e1;
t517 = m(6) * (t182 * t353 - t183 * t354);
t516 = m(6) * (-t353 * t178 - t180 * t354);
t515 = m(6) * t126;
t514 = m(6) * (t204 * t353 - t205 * t354);
t512 = m(6) * (-t353 * t215 - t217 * t354);
t502 = m(6) * qJD(3);
t469 = t353 * t363;
t243 = Icges(4,4) * t469 - Icges(4,2) * t470 - Icges(4,6) * t354;
t449 = t361 * t243;
t244 = Icges(4,6) * t353 + t320 * t354;
t448 = t361 * t244;
t329 = Icges(4,4) * t470;
t245 = Icges(4,1) * t469 - Icges(4,5) * t354 - t329;
t447 = t363 * t245;
t246 = Icges(4,5) * t353 + t322 * t354;
t446 = t363 * t246;
t76 = 0.2e1 * (t126 / 0.4e1 - t99 / 0.4e1) * m(6);
t444 = t76 * qJD(2);
t241 = Icges(4,5) * t469 - Icges(4,6) * t470 - Icges(4,3) * t354;
t434 = -t353 * t241 - t354 * t447;
t393 = Icges(4,5) * t363 - Icges(4,6) * t361;
t242 = Icges(4,3) * t353 + t354 * t393;
t433 = t353 * t242 + t354 * t446;
t423 = t321 * t353 + t243;
t422 = -t321 * t354 - t244;
t421 = -Icges(4,2) * t469 + t245 - t329;
t420 = -t319 * t354 + t246;
t413 = t83 * t126 + t441;
t411 = -t303 - t358;
t408 = (t562 * t351 + (t558 * t354 + (t559 - t561) * t353) * t354) * t522 + (t561 * t352 + (t559 * t353 + (t558 - t562) * t354) * t353) * t521;
t201 = t353 * t446;
t402 = t354 * t242 - t201;
t399 = -t241 + t448;
t398 = t415 * t510;
t397 = t79 + t408;
t392 = -Icges(4,5) * t361 - Icges(4,6) * t363;
t381 = t409 - t358;
t372 = t361 * t421 + t363 * t423;
t371 = -t361 * t420 + t363 * t422;
t367 = -t379 + (t607 * t353 + t560 * t354 + t581 * t355 + t583 * t356) * t522 + (t560 * t353 - t580 * t355 + t582 * t356 - t597) * t521;
t366 = -t377 + (t574 * t355 + t592 * t356) * (t522 + t523);
t325 = -rSges(4,2) * t361 + t501;
t279 = t392 * t354;
t278 = t392 * t353;
t236 = t411 * t354;
t234 = t411 * t353;
t181 = t381 * t354;
t179 = t381 * t353;
t167 = -t284 * t353 - t285 * t354;
t154 = m(5) * t159;
t142 = -t216 * t354 + t218 * t353;
t135 = -t398 + t159;
t134 = t512 / 0.2e1;
t128 = m(6) * t142 * qJD(4);
t125 = t514 / 0.2e1;
t118 = t516 / 0.2e1;
t117 = t517 / 0.2e1;
t115 = -t354 * t448 + t433;
t114 = -t354 * t449 - t434;
t113 = -t353 * t448 - t402;
t90 = -t398 + t99;
t75 = -t114 * t354 + t115 * t353;
t74 = -(-t353 * (-t447 + t449) - t354 * t241) * t354 + t113 * t353;
t73 = t134 - t514 / 0.2e1;
t72 = t134 + t125;
t71 = t125 - t512 / 0.2e1;
t48 = t154 + t515 / 0.2e1 + t99 * t540;
t47 = t118 + t117;
t46 = t118 - t517 / 0.2e1;
t45 = t117 - t516 / 0.2e1;
t28 = t377 + t524 + t531;
t25 = (t113 - t201 + (t242 + t449) * t354 + t434) * t354 + t433 * t353;
t24 = (t354 * t399 + t115 - t433) * t354 + (t353 * t399 + t114 + t402) * t353;
t18 = (t321 / 0.2e1 + t320 / 0.2e1) * t363 + t563 + t538 + t532 + t525 + t377;
t7 = t397 + t528;
t6 = t408 - t590;
t4 = t379 + t555;
t3 = t379 - t555;
t2 = (t75 / 0.2e1 - t25 / 0.2e1) * t354 + (t24 / 0.2e1 + t74 / 0.2e1) * t353 + t379;
t1 = t367 + t504 + t505;
t5 = [t18 * qJD(3) + t28 * qJD(4) + qJD(5) * t568, 0, t18 * qJD(1) + t1 * qJD(4) + t47 * qJD(5) + (((-t170 * t354 - t171 * t353) * t325 + (-t284 * t354 + t285 * t353) * t323) * t542 + (t157 * t181 + t158 * t179 - t178 * t183 - t180 * t182) * t540 + (t160 * t236 + t161 * t234) * t541) * t545 + (t367 + (t361 * t422 + t363 * t420) * t522 + t25 * t567 + (t24 + t74) * t523 + (-t361 * t423 + t363 * t421 + t75) * t521 + (t352 / 0.2e1 + t351 / 0.2e1) * t393) * qJD(3), t28 * qJD(1) + t1 * qJD(3) + t367 * qJD(4) + t72 * qJD(5) + ((-t204 * t217 - t205 * t215 + t443) * t540 + (t380 + (-t260 * t354 + t261 * t353) * t301) * t541) * t543, t47 * qJD(3) + t72 * qJD(4) + t551; 0, 0, t48 * qJD(4) + (t135 * t541 + t167 * t542 + t540 * t90) * t545, t48 * qJD(3) + (t154 + t515) * qJD(4), 0; (t366 - (t321 + t320) * t363 / 0.2e1 - t563) * qJD(1) + t2 * qJD(3) + t3 * qJD(4) + t46 * qJD(5) + (-t525 / 0.4e1 - t532 / 0.4e1 - t538 / 0.4e1) * t546, qJD(4) * t76, t2 * qJD(1) + (m(6) * (-t178 * t179 - t180 * t181 + t54 * t90) + m(5) * (t135 * t88 - t234 * t550 - t236 * t549) + (t351 * t279 + (t372 * t354 + (-t278 + t371) * t353) * t354) * t522 + (t352 * t278 + (t371 * t353 + (-t279 + t372) * t354) * t353) * t521 + m(4) * (t323 * t325 * t415 + (t353 * (rSges(4,1) * t469 - t416) + t354 * (rSges(4,3) * t353 + t354 * t501 - t331)) * t167) + t408) * qJD(3) + t6 * qJD(4), t3 * qJD(1) + t444 + t6 * qJD(3) + ((t413 + t32) * t540 + (t44 + t80) * t541) * t543 + (t408 - t79 - t528) * qJD(4), t46 * qJD(1); t366 * qJD(1) + t4 * qJD(3) + t379 * qJD(4) + t73 * qJD(5) + (-t524 / 0.4e1 - t531 / 0.4e1) * t546, -qJD(3) * t76, t4 * qJD(1) - t444 + t7 * qJD(4) + ((-t179 * t215 - t181 * t217 + t54 * t99 + t83 * t90 + t442) * t540 + (t135 * t144 + (-t234 * t353 - t236 * t354) * t301 + t44) * t541) * t545 + (t408 + t590) * qJD(3), t379 * qJD(1) + t7 * qJD(3) + (m(6) * t413 + t397) * qJD(4), t73 * qJD(1); t45 * qJD(3) + t71 * qJD(4) - t551, 0, t45 * qJD(1) + (-t179 * t354 + t181 * t353) * t502 + t128, t71 * qJD(1) + t142 * t502 + t128, 0;];
Cq = t5;
