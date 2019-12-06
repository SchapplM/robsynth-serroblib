% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:32
% EndTime: 2019-12-05 16:43:52
% DurationCPUTime: 9.34s
% Computational Cost: add. (29217->362), mult. (21880->494), div. (0->0), fcn. (20188->6), ass. (0->241)
t599 = Icges(5,4) + Icges(6,4);
t597 = Icges(5,1) + Icges(6,1);
t596 = Icges(5,5) + Icges(6,5);
t591 = Icges(5,6) + Icges(6,6);
t358 = qJ(3) + qJ(4);
t353 = sin(t358);
t600 = t599 * t353;
t595 = Icges(5,2) + Icges(6,2);
t354 = cos(t358);
t598 = -t591 * t353 + t596 * t354;
t577 = t597 * t354 - t600;
t594 = Icges(5,3) + Icges(6,3);
t593 = t599 * t354;
t357 = pkin(8) + qJ(2);
t351 = sin(t357);
t352 = cos(t357);
t586 = t596 * t351 + t577 * t352;
t590 = -t595 * t353 + t593;
t462 = t351 * t353;
t589 = t599 * t462;
t588 = t598 * t352;
t461 = t351 * t354;
t579 = -t594 * t352 + t596 * t461 - t591 * t462;
t587 = t594 * t351 + t588;
t565 = -t596 * t352 + t597 * t461 - t589;
t578 = t595 * t354 + t600;
t585 = t597 * t353 + t593;
t584 = t586 * t461;
t583 = -t591 * t352 + t599 * t461 - t462 * t595;
t582 = t591 * t351 + t590 * t352;
t349 = t351 ^ 2;
t350 = t352 ^ 2;
t411 = t349 + t350;
t492 = rSges(6,2) * t354;
t300 = rSges(6,1) * t353 + t492;
t413 = rSges(6,1) * t462 + rSges(6,2) * t461;
t420 = -t350 * t300 - t351 * t413;
t502 = pkin(4) * t353;
t126 = -t411 * t502 + t420;
t359 = sin(qJ(3));
t503 = pkin(3) * t359;
t317 = -t502 - t503;
t378 = t300 - t317;
t176 = t378 * t351;
t178 = t378 * t352;
t493 = rSges(6,1) * t354;
t501 = pkin(4) * t354;
t405 = rSges(6,2) * t353 - t493 - t501;
t216 = t405 * t351;
t218 = t405 * t352;
t438 = -t176 * t216 - t178 * t218;
t361 = -pkin(7) - pkin(6);
t332 = t351 * t361;
t345 = t352 * pkin(6);
t450 = t352 * t361;
t360 = cos(qJ(3));
t356 = t360 * pkin(3);
t348 = t356 + pkin(2);
t500 = pkin(2) - t348;
t431 = -t351 * (t500 * t351 - t345 - t450) + t352 * (-pkin(6) * t351 - t500 * t352 - t332);
t304 = t348 + t501;
t580 = rSges(6,3) + qJ(5) - t361;
t158 = -rSges(6,1) * t461 + rSges(6,2) * t462 - t351 * t304 + t580 * t352;
t453 = t352 * t354;
t454 = t352 * t353;
t548 = -rSges(6,2) * t454 + t580 * t351;
t83 = (rSges(6,1) * t453 + t332 + (t304 - t348) * t352 + t548) * t352 + (-t351 * t348 - t158 - t450) * t351;
t53 = t83 + t431;
t32 = t53 * t126 + t438;
t301 = rSges(5,1) * t353 + rSges(5,2) * t354;
t260 = t301 * t351;
t261 = t301 * t352;
t157 = -t351 * t260 - t352 * t261;
t494 = rSges(5,1) * t354;
t303 = -rSges(5,2) * t353 + t494;
t374 = t301 + t503;
t541 = t374 * t352;
t542 = t374 * t351;
t402 = -rSges(5,2) * t454 + rSges(5,3) * t351;
t414 = rSges(5,2) * t462 + t352 * rSges(5,3);
t144 = t351 * (rSges(5,1) * t461 - t414) + t352 * (rSges(5,1) * t453 + t402);
t87 = t144 + t431;
t44 = t87 * t157 + (t351 * t542 + t352 * t541) * t303;
t581 = -m(5) * t44 - m(6) * t32;
t576 = -t596 * t353 - t591 * t354;
t575 = -t587 * t352 + t584;
t574 = -t578 * t352 + t586;
t573 = -t461 * t595 + t565 - t589;
t572 = -t585 * t352 - t582;
t571 = t585 * t351 + t583;
t570 = -t579 * t351 - t565 * t453;
t561 = t587 * t351 + t586 * t453;
t569 = -t585 - t590;
t568 = -t582 * t462 + t575;
t567 = t583 * t454 + t570;
t566 = -t582 * t454 + t561;
t562 = t582 * t353 - t579;
t560 = t583 * t353;
t533 = m(5) / 0.2e1;
t532 = m(6) / 0.2e1;
t159 = (t304 + t493) * t352 + t548;
t559 = m(6) * (t158 * t352 + t159 * t351);
t515 = t351 / 0.2e1;
t514 = -t352 / 0.2e1;
t558 = t352 / 0.2e1;
t480 = Icges(4,4) * t359;
t319 = Icges(4,2) * t360 + t480;
t322 = Icges(4,1) * t360 - t480;
t554 = (t322 / 0.2e1 - t319 / 0.2e1) * t359;
t553 = t576 * t351;
t552 = t576 * t352;
t551 = (t577 - t578) * t354 + t569 * t353;
t550 = -t574 * t353 + t572 * t354;
t549 = t573 * t353 + t571 * t354;
t180 = -t317 * t351 + t413;
t287 = t352 * t317;
t181 = -t300 * t352 + t287;
t406 = t300 + t502;
t215 = t406 * t351;
t217 = t406 * t352;
t403 = t348 + t494;
t166 = -t403 * t351 + t414 - t450;
t167 = t403 * t352 - t332 + t402;
t376 = (-t166 * t352 - t167 * t351) * t303;
t439 = t218 * t158 + t216 * t159;
t498 = (-t180 * t217 - t181 * t215 + t439) * t532 + (t376 + (t351 * t541 - t352 * t542) * t301) * t533;
t204 = pkin(4) * t462 + t413;
t205 = (-t492 + (-rSges(6,1) - pkin(4)) * t353) * t352;
t499 = (-t176 * t205 - t178 * t204 + t439) * t532 + (-t260 * t541 + t261 * t542 + t376) * t533;
t547 = t498 - t499;
t543 = qJD(2) * t559;
t355 = Icges(4,4) * t360;
t320 = -Icges(4,2) * t359 + t355;
t321 = Icges(4,1) * t359 + t355;
t373 = -t569 * t354 / 0.2e1 + (-t578 / 0.2e1 + t577 / 0.2e1) * t353;
t375 = (t566 * t351 + t567 * t352) * t558 + (t561 * t351 + ((t560 + t587) * t352 + t568 + t570 - t584) * t352) * t514 + ((t562 * t351 - t567 + t568 - t575) * t351 + ((t562 + t579) * t352 + (-t565 * t354 + t560) * t351 - t561 + t566) * t352) * t515;
t538 = 0.4e1 * qJD(2);
t537 = 2 * qJD(3);
t535 = 2 * qJD(4);
t534 = m(4) / 0.2e1;
t79 = t411 * t301 * t303 + t144 * t157;
t78 = m(5) * t79;
t525 = m(5) * (t166 * t542 - t167 * t541);
t524 = m(5) * (t166 * t260 - t167 * t261);
t437 = -t215 * t216 - t217 * t218;
t452 = t352 * t359;
t98 = t349 * (t317 + t503) + t352 * (pkin(3) * t452 + t287) + t420;
t521 = m(6) * (t83 * t98 + t437);
t518 = m(6) * (t158 * t180 + t159 * t181);
t517 = m(6) * (t158 * t204 + t159 * t205);
t516 = -t351 / 0.2e1;
t495 = rSges(4,1) * t360;
t408 = pkin(2) + t495;
t460 = t351 * t359;
t412 = rSges(4,2) * t460 + t352 * rSges(4,3);
t182 = -t408 * t351 + t345 + t412;
t329 = rSges(4,2) * t452;
t183 = -t329 + t408 * t352 + (rSges(4,3) + pkin(6)) * t351;
t323 = rSges(4,1) * t359 + rSges(4,2) * t360;
t284 = t323 * t351;
t285 = t323 * t352;
t512 = m(4) * (t182 * t284 - t183 * t285);
t509 = m(6) * (t180 * t351 - t181 * t352);
t508 = m(6) * (-t351 * t176 - t178 * t352);
t507 = m(6) * t126;
t506 = m(6) * (t204 * t351 - t205 * t352);
t504 = m(6) * (-t351 * t215 - t217 * t352);
t496 = m(6) * qJD(3);
t459 = t351 * t360;
t243 = Icges(4,4) * t459 - Icges(4,2) * t460 - Icges(4,6) * t352;
t469 = t243 * t359;
t451 = t352 * t360;
t76 = 0.2e1 * (t126 / 0.4e1 - t98 / 0.4e1) * m(6);
t440 = t76 * qJD(1);
t241 = Icges(4,5) * t459 - Icges(4,6) * t460 - Icges(4,3) * t352;
t327 = Icges(4,4) * t460;
t245 = Icges(4,1) * t459 - Icges(4,5) * t352 - t327;
t430 = -t351 * t241 - t245 * t451;
t389 = Icges(4,5) * t360 - Icges(4,6) * t359;
t242 = Icges(4,3) * t351 + t389 * t352;
t246 = Icges(4,5) * t351 + t322 * t352;
t429 = t351 * t242 + t246 * t451;
t419 = t321 * t351 + t243;
t244 = Icges(4,6) * t351 + t320 * t352;
t418 = -t321 * t352 - t244;
t417 = -Icges(4,2) * t459 + t245 - t327;
t416 = -t319 * t352 + t246;
t409 = t83 * t126 + t437;
t407 = -t303 - t356;
t404 = (t552 * t349 + (t549 * t352 + (t550 - t553) * t351) * t352) * t515 + (t553 * t350 + (t550 * t351 + (t549 - t552) * t352) * t351) * t514;
t201 = t246 * t459;
t398 = t242 * t352 - t201;
t395 = t244 * t359 - t241;
t394 = t411 * t503;
t393 = t78 + t404;
t388 = -Icges(4,5) * t359 - Icges(4,6) * t360;
t377 = t405 - t356;
t368 = t417 * t359 + t419 * t360;
t367 = -t416 * t359 + t418 * t360;
t363 = -t375 + (t598 * t351 + t551 * t352 + t572 * t353 + t574 * t354) * t515 + (t551 * t351 - t571 * t353 + t573 * t354 - t588) * t514;
t362 = -t373 + (t565 * t353 + t583 * t354) * (t515 + t516);
t324 = -rSges(4,2) * t359 + t495;
t279 = t388 * t352;
t278 = t388 * t351;
t236 = t407 * t352;
t234 = t407 * t351;
t179 = t377 * t352;
t177 = t377 * t351;
t165 = -t284 * t351 - t285 * t352;
t153 = m(5) * t157;
t142 = -t216 * t352 + t218 * t351;
t135 = -t394 + t157;
t134 = t504 / 0.2e1;
t128 = m(6) * t142 * qJD(4);
t125 = t506 / 0.2e1;
t118 = t508 / 0.2e1;
t117 = t509 / 0.2e1;
t115 = -t244 * t452 + t429;
t114 = -t243 * t452 - t430;
t113 = -t244 * t460 - t398;
t90 = -t394 + t98;
t75 = -t114 * t352 + t115 * t351;
t74 = -(-(-t245 * t360 + t469) * t351 - t241 * t352) * t352 + t113 * t351;
t72 = t134 - t506 / 0.2e1;
t71 = t134 + t125;
t70 = t125 - t504 / 0.2e1;
t48 = t153 + t507 / 0.2e1 + t98 * t532;
t47 = t118 + t117;
t46 = t118 - t509 / 0.2e1;
t45 = t117 - t508 / 0.2e1;
t28 = t373 + t517 + t524;
t25 = (t113 - t201 + (t242 + t469) * t352 + t430) * t352 + t429 * t351;
t24 = (t395 * t352 + t115 - t429) * t352 + (t395 * t351 + t114 + t398) * t351;
t22 = (t321 / 0.2e1 + t320 / 0.2e1) * t360 + t554 + t512 + t525 + t518 + t373;
t7 = t393 + t521;
t6 = t404 - t581;
t4 = t375 - t547;
t3 = t375 + t547;
t2 = (-t25 / 0.2e1 + t75 / 0.2e1) * t352 + (t74 / 0.2e1 + t24 / 0.2e1) * t351 + t375;
t1 = t363 + t498 + t499;
t5 = [0, 0, t48 * qJD(4) + (t135 * t533 + t165 * t534 + t90 * t532) * t537, t48 * qJD(3) + (t153 + t507) * qJD(4), 0; 0, t22 * qJD(3) + t28 * qJD(4) + qJD(5) * t559, t22 * qJD(2) + t1 * qJD(4) + t47 * qJD(5) + (((-t182 * t352 - t183 * t351) * t324 + (-t284 * t352 + t285 * t351) * t323) * t534 + (t166 * t236 + t167 * t234) * t533 + (t158 * t179 + t159 * t177 - t176 * t181 - t178 * t180) * t532) * t537 + (t363 + t25 * t558 + (t359 * t418 + t360 * t416) * t515 + (t74 + t24) * t516 + (-t359 * t419 + t360 * t417 + t75) * t514 + (t349 / 0.2e1 + t350 / 0.2e1) * t389) * qJD(3), t28 * qJD(2) + t1 * qJD(3) + t363 * qJD(4) + t71 * qJD(5) + ((-t204 * t217 - t205 * t215 + t439) * t532 + (t376 + (-t260 * t352 + t261 * t351) * t301) * t533) * t535, t47 * qJD(3) + t71 * qJD(4) + t543; qJD(4) * t76, (t362 - (t321 + t320) * t360 / 0.2e1 - t554) * qJD(2) + t2 * qJD(3) + t4 * qJD(4) + t46 * qJD(5) + (-t512 / 0.4e1 - t525 / 0.4e1 - t518 / 0.4e1) * t538, t2 * qJD(2) + (m(6) * (-t176 * t177 - t178 * t179 + t53 * t90) + m(5) * (t135 * t87 - t234 * t542 - t236 * t541) + (t350 * t278 + (t367 * t351 + (-t279 + t368) * t352) * t351) * t514 + m(4) * (t323 * t324 * t411 + (t351 * (rSges(4,1) * t459 - t412) + t352 * (rSges(4,1) * t451 + rSges(4,3) * t351 - t329)) * t165) + (t349 * t279 + (t368 * t352 + (-t278 + t367) * t351) * t352) * t515 + t404) * qJD(3) + t6 * qJD(4), t440 + t4 * qJD(2) + t6 * qJD(3) + ((t409 + t32) * t532 + (t44 + t79) * t533) * t535 + (t404 - t78 - t521) * qJD(4), t46 * qJD(2); -qJD(3) * t76, t362 * qJD(2) + t3 * qJD(3) + t375 * qJD(4) + t72 * qJD(5) + (-t517 / 0.4e1 - t524 / 0.4e1) * t538, -t440 + t3 * qJD(2) + t7 * qJD(4) + ((-t177 * t215 - t179 * t217 + t53 * t98 + t83 * t90 + t438) * t532 + (t135 * t144 + (-t234 * t351 - t236 * t352) * t301 + t44) * t533) * t537 + (t404 + t581) * qJD(3), t375 * qJD(2) + t7 * qJD(3) + (m(6) * t409 + t393) * qJD(4), t72 * qJD(2); 0, t45 * qJD(3) + t70 * qJD(4) - t543, t45 * qJD(2) + (-t177 * t352 + t179 * t351) * t496 + t128, t70 * qJD(2) + t142 * t496 + t128, 0;];
Cq = t5;
