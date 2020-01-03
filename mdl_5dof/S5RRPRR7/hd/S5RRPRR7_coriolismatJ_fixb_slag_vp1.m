% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:29
% EndTime: 2019-12-31 20:15:38
% DurationCPUTime: 5.67s
% Computational Cost: add. (33791->395), mult. (26712->512), div. (0->0), fcn. (24056->8), ass. (0->265)
t374 = qJ(1) + qJ(2);
t372 = cos(t374);
t375 = sin(qJ(4));
t483 = pkin(4) * t375;
t352 = t372 * t483;
t354 = t372 * qJ(3);
t370 = sin(t374);
t373 = qJ(4) + qJ(5);
t369 = sin(t373);
t371 = cos(t373);
t405 = rSges(6,1) * t369 + rSges(6,2) * t371;
t453 = t372 * t405;
t387 = -t370 * rSges(6,3) + t453;
t533 = -pkin(8) - pkin(7);
t211 = t352 + t354 + (-pkin(2) + t533) * t370 + t387;
t481 = sin(qJ(1)) * pkin(1);
t204 = t211 - t481;
t353 = t372 * t533;
t560 = -t405 - t483;
t212 = -t353 + (rSges(6,3) + pkin(2)) * t372 + (qJ(3) - t560) * t370;
t484 = cos(qJ(1)) * pkin(1);
t205 = t212 + t484;
t111 = t204 * t372 + t205 * t370;
t119 = t211 * t372 + t212 * t370;
t377 = cos(qJ(4));
t406 = rSges(5,1) * t375 + rSges(5,2) * t377;
t386 = -t370 * rSges(5,3) + t372 * t406;
t534 = pkin(2) + pkin(7);
t237 = -t534 * t370 + t354 + t386;
t233 = t237 - t481;
t238 = (rSges(5,3) + t534) * t372 + (qJ(3) + t406) * t370;
t234 = t238 + t484;
t139 = t233 * t372 + t234 * t370;
t143 = t237 * t372 + t238 * t370;
t565 = pkin(2) - rSges(4,2);
t286 = t372 * rSges(4,3) - t565 * t370 + t354;
t269 = t286 - t481;
t287 = (rSges(4,3) + qJ(3)) * t370 + t565 * t372;
t270 = t287 + t484;
t179 = t269 * t372 + t270 * t370;
t194 = t286 * t372 + t287 * t370;
t535 = m(6) / 0.2e1;
t536 = m(5) / 0.2e1;
t569 = m(4) / 0.2e1;
t416 = (t194 + t179) * t569 + (t119 + t111) * t535 + (t143 + t139) * t536;
t417 = (t179 - t194) * t569 + (t111 - t119) * t535 + (t139 - t143) * t536;
t14 = t417 - t416;
t574 = t14 * qJD(1);
t469 = Icges(6,4) * t369;
t401 = Icges(6,2) * t371 + t469;
t273 = Icges(6,6) * t372 + t370 * t401;
t447 = t370 * t371;
t339 = Icges(6,4) * t447;
t448 = t369 * t370;
t275 = Icges(6,1) * t448 + Icges(6,5) * t372 + t339;
t396 = -t273 * t371 - t275 * t369;
t573 = t372 * t396;
t572 = t205 - t212;
t468 = Icges(6,4) * t371;
t332 = -Icges(6,2) * t369 + t468;
t403 = Icges(6,1) * t369 + t468;
t571 = t332 + t403;
t470 = Icges(5,4) * t377;
t343 = -Icges(5,2) * t375 + t470;
t404 = Icges(5,1) * t375 + t470;
t570 = t343 + t404;
t568 = -t370 / 0.2e1;
t513 = t370 / 0.2e1;
t511 = t372 / 0.2e1;
t509 = m(3) * (t484 * (-rSges(3,1) * t370 - rSges(3,2) * t372) + (t372 * rSges(3,1) - t370 * rSges(3,2)) * t481);
t505 = m(4) * (-t287 * t269 + t270 * t286);
t397 = Icges(6,5) * t369 + Icges(6,6) * t371;
t564 = t370 * t397;
t399 = Icges(5,5) * t375 + Icges(5,6) * t377;
t563 = t370 * t399;
t454 = t370 * t405;
t562 = t372 * t397;
t561 = t372 * t399;
t367 = t370 ^ 2;
t368 = t372 ^ 2;
t419 = t367 + t368;
t244 = t419 * t405;
t512 = -t372 / 0.2e1;
t559 = t511 + t512;
t276 = -Icges(6,5) * t370 + t372 * t403;
t428 = t332 * t372 + t276;
t429 = -Icges(6,2) * t448 + t275 + t339;
t274 = -Icges(6,6) * t370 + t372 * t401;
t334 = Icges(6,1) * t371 - t469;
t430 = -t334 * t372 + t274;
t431 = -t334 * t370 + t273;
t558 = -(t370 * t430 - t372 * t431) * t369 + (t370 * t428 - t372 * t429) * t371;
t302 = rSges(6,1) * t447 - rSges(6,2) * t448;
t337 = rSges(6,1) * t371 - rSges(6,2) * t369;
t303 = t337 * t372;
t106 = t211 * t303 + t212 * t302;
t557 = m(6) * t106;
t250 = t372 * t387;
t277 = t372 * rSges(6,3) + t454;
t446 = t370 * t375;
t122 = -t372 * (-pkin(8) * t370 + t352) - t250 + (-pkin(4) * t446 + t372 * pkin(7) - t277 + t353) * t370;
t218 = -t302 * t370 - t372 * t303;
t445 = t370 * t377;
t351 = pkin(4) * t445;
t452 = t337 * t370;
t280 = t351 + t452;
t482 = pkin(4) * t377;
t408 = (t337 + t482) * t372;
t398 = Icges(6,5) * t371 - Icges(6,6) * t369;
t296 = t370 * t398;
t297 = t398 * t372;
t480 = (t368 * t296 + (-t372 * t297 - t558) * t370) * t511 + (-t367 * t297 + (t370 * t296 + t558) * t372) * t513;
t18 = t480 + m(6) * (t122 * t218 - t280 * t454 - t408 * t453);
t556 = t18 * qJD(5);
t555 = t406 * t536;
t554 = t370 * t408;
t471 = Icges(5,4) * t375;
t402 = Icges(5,2) * t377 + t471;
t291 = -Icges(5,6) * t370 + t372 * t402;
t293 = -Icges(5,5) * t370 + t372 * t404;
t553 = (t291 * t377 + t293 * t375) * t372;
t552 = (t274 * t371 + t276 * t369) * t372;
t155 = t408 * t204;
t161 = t408 * t211;
t87 = -t212 * t204 + t205 * t211;
t101 = -t238 * t233 + t234 * t237;
t418 = qJD(1) + qJD(2);
t261 = t302 + t351;
t348 = rSges(5,1) * t377 - rSges(5,2) * t375;
t318 = t348 * t370;
t319 = t348 * t372;
t434 = (-t261 * t372 + t554) * t535 + (-t318 * t372 + t319 * t370) * t536;
t461 = t238 * t318;
t463 = t234 * t318;
t479 = (t280 * t572 + t155 - t161) * t535 + (t463 - t461 + (t233 - t237) * t319) * t536;
t100 = t205 * t261 + t155;
t103 = t212 * t261 + t161;
t112 = t233 * t319 + t463;
t120 = t237 * t319 + t461;
t550 = (t103 + t100) * t535 + (t120 + t112) * t536;
t549 = t369 * t571 - t371 * (-t401 + t334);
t345 = Icges(5,1) * t377 - t471;
t546 = t375 * t570 - t377 * (-t402 + t345);
t424 = t343 * t372 + t293;
t426 = -t345 * t372 + t291;
t545 = t375 * t426 - t377 * t424;
t350 = Icges(5,4) * t445;
t292 = Icges(5,1) * t446 + Icges(5,5) * t372 + t350;
t425 = -Icges(5,2) * t446 + t292 + t350;
t290 = Icges(5,6) * t372 + t370 * t402;
t427 = -t345 * t370 + t290;
t544 = t375 * t427 - t377 * t425;
t543 = -t570 * t377 / 0.2e1 + (-t345 / 0.2e1 + t402 / 0.2e1) * t375;
t411 = -t571 * t371 / 0.2e1 + (t401 / 0.2e1 - t334 / 0.2e1) * t369;
t133 = t372 * (Icges(6,3) * t372 + t564) + t273 * t447 + t275 * t448;
t272 = -Icges(6,3) * t370 + t562;
t134 = -t372 * t272 - t274 * t447 - t276 * t448;
t136 = -t370 * t272 + t552;
t11 = ((t134 - t573) * t370 + (t136 - t552 + (t272 + t396) * t370 + t133) * t372) * t513 + (t133 * t372 + t134 * t370) * t568 + (t136 * t370 + t367 * t272 + (t134 + (t272 - t396) * t372 + t573) * t372) * t511;
t542 = 0.4e1 * qJD(1);
t540 = 0.4e1 * qJD(2);
t539 = 2 * qJD(4);
t105 = t204 * t303 + t205 * t302;
t524 = m(6) * (t106 + t105);
t523 = m(6) * ((t204 - t211) * t303 + t572 * t452);
t166 = t204 * t454;
t432 = t280 * t303 - t302 * t408;
t521 = m(6) * (t205 * t453 - t166 + t432);
t228 = t408 * t452;
t459 = t261 * t337;
t467 = t205 * t405;
t520 = m(6) * (-t166 + t228 + (-t459 + t467) * t372);
t171 = t211 * t454;
t517 = m(6) * (t212 * t453 - t171 + t432);
t465 = t212 * t405;
t516 = m(6) * (-t171 + t228 + (-t459 + t465) * t372);
t515 = m(6) * t87;
t288 = Icges(5,3) * t372 + t563;
t144 = t372 * t288 + t290 * t445 + t292 * t446;
t289 = -Icges(5,3) * t370 + t561;
t145 = -t372 * t289 - t291 * t445 - t293 * t446;
t266 = t370 * t288;
t392 = -t290 * t377 - t292 * t375;
t146 = t372 * t392 + t266;
t147 = -t370 * t289 + t553;
t32 = (-t146 + t266 + t145) * t370 + (t147 - t553 + (t289 + t392) * t370 + t144) * t372;
t33 = t367 * t289 + (t145 - t266 + (t289 - t392) * t372) * t372;
t95 = t144 * t372 + t145 * t370;
t96 = t146 * t372 + t147 * t370;
t2 = (t33 / 0.2e1 + t96 / 0.2e1) * t372 + (-t95 / 0.2e1 + t32 / 0.2e1) * t370 + t11;
t487 = m(6) * (t280 * t372 - t554);
t177 = t487 / 0.2e1;
t64 = t177 - t434;
t510 = t64 * qJD(3) + t2 * qJD(4);
t503 = m(4) * t179;
t502 = m(4) * t194;
t501 = m(5) * t101;
t499 = m(5) * t112;
t498 = m(5) * t120;
t497 = m(5) * t139;
t496 = m(5) * t143;
t494 = m(6) * t100;
t493 = m(6) * t103;
t492 = m(6) * t105;
t490 = m(6) * t111;
t489 = m(6) * t119;
t217 = -t302 * t372 + t303 * t370;
t485 = m(6) * t217;
t474 = -qJD(3) * t485 / 0.2e1 + t11 * qJD(5);
t415 = t485 / 0.2e1;
t201 = qJD(5) * t415;
t62 = t177 + t434;
t473 = t62 * qJD(4) + t201;
t63 = -t487 / 0.2e1 + t434;
t472 = t63 * qJD(4) + t201;
t456 = t302 * t337;
t433 = (-t261 + t280) * t408;
t413 = qJD(1) / 0.4e1 + qJD(2) / 0.4e1;
t412 = t419 * t406;
t407 = t524 / 0.2e1 + t411;
t400 = Icges(5,5) * t377 - Icges(5,6) * t375;
t279 = t560 * t370;
t281 = t352 + t453;
t393 = t279 * t370 - t281 * t372;
t312 = t370 * t400;
t388 = -t11 + (t369 * t428 + t371 * t430 + t372 * t549 - t564) * t513 + (-t369 * t429 - t370 * t549 - t371 * t431 - t562) * t511;
t385 = -t411 + t559 * (t274 * t369 - t276 * t371);
t384 = t411 + t543;
t381 = t384 + t550;
t380 = t385 - t543 + t559 * (t375 * t291 - t377 * t293);
t379 = t62 * qJD(3) + (t388 + t32 * t568 + (t372 * t546 + t375 * t424 + t377 * t426 - t563 + t95) * t513 + (t33 + t96) * t512 + (-t370 * t546 - t375 * t425 - t377 * t427 - t561) * t511) * qJD(4);
t313 = t400 * t372;
t239 = t303 * t452;
t236 = m(6) * t244 * qJD(5);
t202 = qJD(3) * t415;
t181 = -t370 * t277 - t250;
t173 = -t419 * t482 + t218;
t81 = t411 + t557;
t78 = t411 + t492;
t76 = t516 / 0.2e1;
t75 = t517 / 0.2e1;
t66 = t520 / 0.2e1;
t65 = t521 / 0.2e1;
t56 = t489 + t496 + t502;
t54 = t523 / 0.2e1;
t53 = t490 + t497 + t503;
t40 = t384 + t493 + t498;
t39 = t384 + t494 + t499;
t29 = t501 + t505 + t509 + t515;
t22 = -t523 / 0.2e1 + t407;
t21 = t54 + t407;
t20 = m(6) * (t181 * t218 - t337 * t244) + t480;
t19 = t20 * qJD(5);
t17 = t54 - t524 / 0.2e1 + t385;
t15 = t416 + t417;
t13 = t381 + t479;
t12 = t381 - t479;
t9 = t380 + t479 - t550;
t8 = t75 - t516 / 0.2e1 + t11;
t7 = t76 - t517 / 0.2e1 + t11;
t6 = t65 - t520 / 0.2e1 + t11;
t5 = t66 - t521 / 0.2e1 + t11;
t4 = t75 + t76 + t388;
t3 = t65 + t66 + t388;
t1 = [qJD(2) * t29 + qJD(3) * t53 + qJD(4) * t39 + qJD(5) * t78, t29 * qJD(1) + t15 * qJD(3) + t13 * qJD(4) + t21 * qJD(5) + 0.2e1 * (t505 / 0.2e1 + t101 * t536 + t509 / 0.2e1 + t87 * t535) * qJD(2), qJD(1) * t53 + qJD(2) * t15 + t473, t39 * qJD(1) + t13 * qJD(2) + t3 * qJD(5) + ((t204 * t279 + t205 * t281 + t433) * t535 - (t233 * t370 - t234 * t372) * t555) * t539 + t379, t78 * qJD(1) + t21 * qJD(2) + t202 + t3 * qJD(4) + ((-t166 + t239 + (-t456 + t467) * t372) * m(6) + t388) * qJD(5); -t14 * qJD(3) + t12 * qJD(4) + t22 * qJD(5) + (-t515 / 0.4e1 - t501 / 0.4e1 - t505 / 0.4e1 - t509 / 0.4e1) * t542, qJD(3) * t56 + qJD(4) * t40 + qJD(5) * t81, qJD(2) * t56 + t473 - t574, t12 * qJD(1) + t40 * qJD(2) + t4 * qJD(5) + ((t211 * t279 + t212 * t281 + t433) * t535 - (t237 * t370 - t238 * t372) * t555) * t539 + t379, t22 * qJD(1) + t81 * qJD(2) + t202 + t4 * qJD(4) + ((-t171 + t239 + (-t456 + t465) * t372) * m(6) + t388) * qJD(5); t14 * qJD(2) + (-t490 / 0.4e1 - t497 / 0.4e1 - t503 / 0.4e1) * t542 + t472, t574 + (-t489 / 0.4e1 - t496 / 0.4e1 - t502 / 0.4e1) * t540 + t472, 0, (t393 * t535 - t412 * t536) * t539 - t236 + t418 * t63, -t236 + 0.2e1 * (-t244 * qJD(4) / 0.2e1 + t413 * t217) * m(6); t9 * qJD(2) + t6 * qJD(5) + (-t494 / 0.4e1 - t499 / 0.4e1) * t542 + t510 + t380 * qJD(1), t9 * qJD(1) + t380 * qJD(2) + t8 * qJD(5) + (-t493 / 0.4e1 - t498 / 0.4e1) * t540 + t510, t418 * t64, (m(5) * ((-t372 * t386 + (-t372 * rSges(5,3) - t370 * t406) * t370) * (-t318 * t370 - t319 * t372) - t348 * t412) + (t368 * t312 + (t545 * t370 + (-t313 - t544) * t372) * t370) * t511 + (-t367 * t313 + (t544 * t372 + (t312 - t545) * t370) * t372) * t513 + m(6) * (t122 * t173 + t279 * t280 - t281 * t408) + t480) * qJD(4) + t556 + t418 * t2, t6 * qJD(1) + t8 * qJD(2) + t18 * qJD(4) + t556; (t385 - t492) * qJD(1) + t17 * qJD(2) + t5 * qJD(4) + t474, t17 * qJD(1) + (t385 - t557) * qJD(2) + t7 * qJD(4) + t474, -0.2e1 * t413 * t485, t5 * qJD(1) + t7 * qJD(2) + ((t173 * t181 + t337 * t393) * m(6) + t480) * qJD(4) + t19, qJD(4) * t20 + t11 * t418 + t19;];
Cq = t1;
