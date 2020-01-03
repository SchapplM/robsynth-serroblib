% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR10
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR10_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:11
% EndTime: 2019-12-31 17:11:30
% DurationCPUTime: 14.37s
% Computational Cost: add. (16906->609), mult. (42060->862), div. (0->0), fcn. (45466->6), ass. (0->342)
t371 = sin(qJ(1));
t618 = -t371 / 0.2e1;
t374 = cos(qJ(1));
t550 = t374 / 0.2e1;
t370 = sin(qJ(2));
t373 = cos(qJ(2));
t334 = pkin(2) * t370 - qJ(3) * t373;
t369 = sin(qJ(4));
t372 = cos(qJ(4));
t414 = rSges(5,1) * t369 + rSges(5,2) * t372;
t580 = -rSges(5,3) * t370 + t373 * t414;
t418 = pkin(6) * t370 + t334 - t580;
t207 = t418 * t371;
t209 = t418 * t374;
t533 = rSges(4,2) * t370;
t413 = rSges(4,3) * t373 + t533;
t456 = t334 - t413;
t248 = t456 * t371;
t250 = t456 * t374;
t366 = t374 * pkin(5);
t496 = t370 * t371;
t429 = rSges(4,1) * t374 - rSges(4,3) * t496;
t507 = qJ(3) * t370;
t431 = -pkin(1) - t507;
t222 = t366 + ((rSges(4,2) - pkin(2)) * t373 + t431) * t371 + t429;
t487 = t373 * t374;
t358 = rSges(4,2) * t487;
t526 = rSges(4,3) + qJ(3);
t537 = pkin(2) * t373;
t223 = -t358 + (rSges(4,1) + pkin(5)) * t371 + (t370 * t526 + pkin(1) + t537) * t374;
t491 = t371 * t373;
t478 = t222 * t487 + t223 * t491;
t360 = pkin(6) * t491;
t430 = pkin(3) * t374 - t360;
t492 = t371 * t372;
t299 = t369 * t374 + t370 * t492;
t488 = t372 * t374;
t493 = t371 * t369;
t300 = -t370 * t493 + t488;
t584 = t300 * rSges(5,1) - t299 * rSges(5,2);
t166 = t366 + ((-rSges(5,3) - pkin(2)) * t373 + t431) * t371 + t430 + t584;
t297 = t370 * t488 - t493;
t494 = t370 * t374;
t298 = t369 * t494 + t492;
t416 = t298 * rSges(5,1) + t297 * rSges(5,2);
t549 = rSges(5,3) + pkin(6);
t446 = pkin(2) + t549;
t168 = (pkin(3) + pkin(5)) * t371 + (t373 * t446 - t431) * t374 + t416;
t483 = t166 * t487 + t168 * t491;
t574 = m(5) / 0.2e1;
t575 = m(4) / 0.2e1;
t520 = (-t248 * t494 + t250 * t496 + t478) * t575 + (-t207 * t494 + t209 * t496 + t483) * t574;
t359 = pkin(2) * t496;
t197 = t359 + (t549 * t370 + (-qJ(3) - t414) * t373) * t371;
t345 = qJ(3) * t487;
t457 = t414 * t487;
t198 = -t446 * t494 + t345 + t457;
t236 = t359 + (-t373 * t526 - t533) * t371;
t419 = -pkin(2) * t494 + t345;
t452 = rSges(4,2) * t494 + rSges(4,3) * t487;
t237 = t419 + t452;
t525 = ((t236 * t374 + t237 * t371) * t370 + t478) * t575 + ((t197 * t374 + t198 * t371) * t370 + t483) * t574;
t595 = t520 - t525;
t617 = t595 * qJD(1);
t616 = Icges(4,1) + Icges(3,3);
t406 = Icges(5,5) * t369 + Icges(5,6) * t372;
t382 = -Icges(5,3) * t370 + t406 * t373;
t515 = Icges(5,4) * t369;
t408 = Icges(5,2) * t372 + t515;
t383 = -Icges(5,6) * t370 + t408 * t373;
t514 = Icges(5,4) * t372;
t410 = Icges(5,1) * t369 + t514;
t384 = -Icges(5,5) * t370 + t410 * t373;
t148 = -t299 * t383 + t300 * t384 - t382 * t491;
t211 = Icges(5,5) * t298 + Icges(5,6) * t297 + Icges(5,3) * t487;
t516 = Icges(5,4) * t298;
t214 = Icges(5,2) * t297 + Icges(5,6) * t487 + t516;
t292 = Icges(5,4) * t297;
t217 = Icges(5,1) * t298 + Icges(5,5) * t487 + t292;
t107 = t211 * t491 + t299 * t214 - t300 * t217;
t213 = -Icges(5,5) * t300 + Icges(5,6) * t299 + Icges(5,3) * t491;
t294 = Icges(5,4) * t300;
t216 = Icges(5,2) * t299 + Icges(5,6) * t491 - t294;
t293 = Icges(5,4) * t299;
t218 = Icges(5,1) * t300 - Icges(5,5) * t491 - t293;
t108 = t213 * t491 + t216 * t299 + t218 * t300;
t403 = t107 * t374 + t108 * t371;
t51 = t148 * t370 + t373 * t403;
t146 = -t297 * t383 - t298 * t384 - t382 * t487;
t105 = t211 * t487 + t297 * t214 + t298 * t217;
t106 = t213 * t487 + t297 * t216 - t298 * t218;
t404 = t105 * t374 + t371 * t106;
t614 = t146 * t370 + t373 * t404;
t326 = Icges(3,5) * t373 - Icges(3,6) * t370;
t327 = -Icges(4,4) * t373 + Icges(4,5) * t370;
t613 = (t326 + t327) * t374 + t616 * t371;
t609 = (-Icges(4,5) + Icges(3,6)) * t496 + (Icges(4,4) - Icges(3,5)) * t491 + t616 * t374;
t505 = t213 * t370;
t605 = t216 * t372 - t218 * t369;
t134 = t373 * t605 - t505;
t517 = Icges(3,4) * t370;
t331 = Icges(3,1) * t373 - t517;
t275 = Icges(3,5) * t371 + t331 * t374;
t510 = Icges(4,6) * t373;
t323 = Icges(4,3) * t370 - t510;
t276 = Icges(4,5) * t371 + t323 * t374;
t611 = -t275 * t491 - t276 * t496;
t610 = t105 * t371 - t106 * t374;
t221 = rSges(5,3) * t491 - t584;
t180 = -t221 * t370 - t580 * t491;
t246 = t580 * t371;
t367 = t371 ^ 2;
t368 = t374 ^ 2;
t451 = t367 + t368;
t608 = t613 * t374 + t611;
t596 = t275 * t487 + t276 * t494 + t613 * t371;
t352 = Icges(3,4) * t496;
t274 = Icges(3,1) * t491 - Icges(3,5) * t374 - t352;
t277 = Icges(4,5) * t374 + Icges(4,6) * t491 - Icges(4,3) * t496;
t607 = -t274 * t487 + t277 * t494 + t609 * t371;
t270 = Icges(3,4) * t491 - Icges(3,2) * t496 - Icges(3,6) * t374;
t347 = Icges(4,6) * t496;
t279 = Icges(4,4) * t374 + Icges(4,2) * t491 - t347;
t606 = t270 * t370 - t279 * t373;
t364 = Icges(3,4) * t373;
t512 = Icges(3,2) * t370;
t271 = Icges(3,6) * t371 + (t364 - t512) * t374;
t348 = Icges(4,6) * t494;
t278 = Icges(4,4) * t371 - Icges(4,2) * t487 + t348;
t598 = -t271 * t494 - t278 * t487 + t596;
t597 = t271 * t370 + t278 * t373 + t609;
t594 = -t270 * t494 - t271 * t496 - t278 * t491 + t279 * t487 - t607 - t608;
t554 = t371 / 0.2e1;
t592 = t373 / 0.2e1;
t552 = -t374 / 0.2e1;
t264 = Icges(5,3) * t373 + t370 * t406;
t328 = Icges(3,2) * t373 + t517;
t489 = t372 * t383;
t498 = t369 * t384;
t588 = t370 * (t331 / 0.2e1 - t328 / 0.2e1 + Icges(4,2) * t592 - Icges(4,6) * t370 - Icges(4,3) * t373 / 0.2e1 - t498 / 0.2e1 - t489 / 0.2e1 + t264 / 0.2e1);
t495 = t370 * t373;
t454 = t451 * t495;
t587 = (m(4) / 0.4e1 + m(5) / 0.4e1) * (t454 - t495);
t586 = t451 * t370;
t407 = -Icges(3,5) * t370 - Icges(3,6) * t373;
t409 = Icges(4,4) * t370 + Icges(4,5) * t373;
t583 = (t407 + t409) * t374;
t388 = rSges(5,3) * t487 + t416;
t581 = -t374 * t221 + t371 * t388;
t308 = (Icges(5,2) * t369 - t514) * t373;
t313 = (-Icges(5,1) * t372 + t515) * t373;
t579 = -t369 * (t313 / 0.2e1 + t383 / 0.2e1) - t372 * (-t384 / 0.2e1 + t308 / 0.2e1);
t578 = 0.4e1 * qJD(1);
t577 = 2 * qJD(2);
t573 = t51 / 0.2e1;
t243 = t383 * t374;
t245 = t384 * t374;
t401 = -t214 * t372 - t217 * t369;
t392 = t382 * t374 - t401;
t93 = (-t243 * t372 - t245 * t369 + t211) * t373 + t392 * t370;
t572 = t93 / 0.2e1;
t242 = t383 * t371;
t244 = t384 * t371;
t391 = t382 * t371 + t605;
t94 = (-t242 * t372 - t244 * t369 + t213) * t373 + t391 * t370;
t571 = -t94 / 0.2e1;
t224 = Icges(5,5) * t297 - Icges(5,6) * t298;
t475 = -Icges(5,2) * t298 + t217 + t292;
t477 = -Icges(5,1) * t297 + t214 + t516;
t98 = t224 * t370 + (t369 * t477 - t372 * t475) * t373;
t570 = t98 / 0.2e1;
t247 = -rSges(5,3) * t494 + t457;
t112 = (t246 * t374 - t247 * t371) * t373 + t581 * t370;
t530 = rSges(5,3) * t373;
t282 = t370 * t414 + t530;
t138 = (t282 * t371 - t221) * t373;
t139 = ((-t282 + t530) * t374 + t416) * t373 + (-t374 * t580 + t247) * t370;
t150 = t581 * t373;
t181 = -t370 * t388 - t487 * t580;
t482 = t180 * t487 - t181 * t491;
t567 = m(5) * (-t112 * t373 + (t138 * t374 + t139 * t371 - t150) * t370 + t482);
t566 = m(5) * (-t112 * t150 + t138 * t180 - t139 * t181);
t565 = m(5) * (t138 * t166 + t139 * t168 + t180 * t197 - t181 * t198);
t230 = rSges(5,1) * t297 - rSges(5,2) * t298;
t231 = rSges(5,1) * t299 + rSges(5,2) * t300;
t318 = (-rSges(5,1) * t372 + rSges(5,2) * t369) * t373;
t562 = m(5) * (-t207 * t230 + t209 * t231 + (-t166 * t374 - t168 * t371) * t318);
t338 = t507 + t537;
t459 = t451 * t338;
t118 = (t221 - t430) * t371 + (t371 * pkin(3) + t487 * t549 + t416) * t374 + t459;
t479 = -t207 * t491 - t209 * t487;
t560 = m(5) * (t118 * t586 + t479);
t559 = m(5) * (-t150 * t586 + t482);
t556 = m(5) * (t166 * t197 + t168 * t198);
t555 = t370 / 0.2e1;
t553 = t371 / 0.4e1;
t551 = -t374 / 0.4e1;
t535 = rSges(3,1) * t373;
t432 = pkin(1) + t535;
t453 = rSges(3,2) * t496 + t374 * rSges(3,3);
t238 = -t371 * t432 + t366 + t453;
t356 = rSges(3,2) * t494;
t239 = -t356 + t432 * t374 + (rSges(3,3) + pkin(5)) * t371;
t336 = rSges(3,1) * t370 + rSges(3,2) * t373;
t316 = t336 * t371;
t319 = t336 * t374;
t548 = m(3) * (t238 * t316 - t239 * t319);
t156 = -t371 * (rSges(4,2) * t491 + t429) + t374 * (t371 * rSges(4,1) + rSges(4,3) * t494 - t358) + t459;
t472 = -t248 * t491 - t250 * t487;
t545 = m(4) * (t156 * t586 + t472);
t543 = m(4) * (t222 * t236 + t223 * t237);
t542 = m(4) * (-t222 * t496 + t223 * t494);
t541 = m(5) * (-t166 * t496 + t168 * t494);
t540 = m(5) * (-t180 * t496 - t181 * t494);
t162 = t230 * t374 + t371 * t231;
t539 = m(5) * (-t162 * t373 - t318 * t586);
t399 = t230 * t371 - t231 * t374;
t538 = m(5) * t399 * t370;
t528 = t371 * t51;
t527 = t374 * t614;
t518 = Icges(3,1) * t370;
t506 = t211 * t370;
t504 = t382 * t370;
t272 = Icges(5,5) * t373 + t370 * t410;
t499 = t369 * t272;
t305 = (-Icges(5,5) * t372 + Icges(5,6) * t369) * t373;
t497 = t370 * t305;
t268 = Icges(5,6) * t373 + t370 * t408;
t490 = t372 * t268;
t476 = -Icges(5,1) * t299 + t216 - t294;
t474 = Icges(5,2) * t300 - t218 + t293;
t467 = -t383 - t313;
t411 = -t364 - t518;
t314 = t411 * t371;
t466 = t270 - t314;
t465 = -t384 + t308;
t311 = -Icges(3,2) * t491 - t352;
t464 = t274 + t311;
t405 = Icges(4,2) * t370 + t510;
t303 = t405 * t371;
t463 = t277 + t303;
t301 = Icges(4,3) * t491 + t347;
t462 = t279 - t301;
t461 = -t282 - t338;
t460 = t371 * (qJ(3) * t491 - t359) + t374 * t419;
t455 = rSges(4,2) * t373 - rSges(4,3) * t370 - t338;
t450 = qJD(1) * t373;
t449 = qJD(2) * t370;
t448 = qJD(2) * t373;
t447 = qJD(4) * t373;
t443 = -t51 / 0.2e1 + t573;
t398 = t489 + t498;
t390 = t264 - t398;
t378 = t373 * t390 + t504;
t103 = t297 * t268 + t298 * t272 + t374 * t378;
t380 = t373 * t392 - t506;
t75 = t297 * t243 + t298 * t245 + t374 * t380;
t379 = t373 * t391 - t505;
t76 = t297 * t242 + t298 * t244 + t374 * t379;
t12 = (t371 * t76 + t374 * t75 + t146) * t373 + (t103 - t404) * t370;
t85 = t224 * t487 + t297 * t475 - t298 * t477;
t225 = Icges(5,5) * t299 + Icges(5,6) * t300;
t86 = t225 * t487 + t297 * t474 - t298 * t476;
t40 = t85 * t371 - t374 * t86;
t442 = t40 / 0.2e1 - t12 / 0.2e1;
t102 = t268 * t299 - t272 * t300 + t371 * t378;
t73 = t243 * t299 - t245 * t300 + t371 * t380;
t74 = t242 * t299 - t244 * t300 + t371 * t379;
t11 = (t371 * t74 + t374 * t73 + t148) * t373 + (t102 - t403) * t370;
t87 = t224 * t491 + t299 * t475 + t300 * t477;
t88 = t225 * t491 + t299 * t474 + t300 * t476;
t41 = t87 * t371 - t374 * t88;
t441 = t41 / 0.2e1 - t11 / 0.2e1;
t437 = t491 / 0.4e1;
t433 = t327 / 0.2e1 + t326 / 0.2e1;
t315 = t411 * t374;
t428 = (-t271 + t315) * t371;
t312 = t328 * t374;
t427 = (-t275 + t312) * t371;
t304 = t405 * t374;
t426 = (t276 - t304) * t371;
t302 = Icges(4,3) * t487 + t348;
t425 = (t278 + t302) * t371;
t125 = t297 * t465 - t298 * t467 + t305 * t487;
t126 = t299 * t465 + t300 * t467 + t305 * t491;
t99 = t225 * t370 + (t369 * t476 - t372 * t474) * t373;
t417 = t562 / 0.2e1 + (t125 + t98) * t553 + (t126 + t99) * t551;
t133 = t373 * t401 + t506;
t402 = t133 * t374 - t134 * t371;
t395 = -m(5) * (-t166 * t231 + t168 * t230) - t497 / 0.2e1;
t387 = t614 * t551 - t51 * t553 + t528 / 0.4e1 + t527 / 0.4e1 + (t437 - t491 / 0.4e1) * t610;
t386 = (t609 * t374 + (t274 * t373 - t277 * t370 - t606) * t371) * t552 + (t597 * t374 - t596 + t598) * t550 + (t597 * t371 + t594 + t608) * t554;
t385 = t596 * t618 + t598 * t554 + ((t606 + t613) * t374 + t594 + t607 + t611) * t552;
t132 = (-t382 - t490 - t499) * t373 + t390 * t370;
t171 = t373 * t398 - t504;
t381 = t132 * t555 + t171 * t592 + t565 / 0.2e1 - (-t134 + t148) * t496 / 0.4e1 - (t133 + t146) * t494 / 0.4e1 + (t102 + t94) * t437 + (t103 + t93) * t487 / 0.4e1;
t375 = t499 / 0.2e1 + t490 / 0.2e1 + t382 / 0.2e1 - t364 - t518 / 0.2e1 + t512 / 0.2e1 - t405 / 0.2e1 + t323 / 0.2e1;
t340 = -rSges(3,2) * t370 + t535;
t309 = t409 * t371;
t306 = t407 * t371;
t251 = t455 * t374;
t249 = t455 * t371;
t210 = (-pkin(6) * t373 + t461) * t374;
t208 = t371 * t461 - t360;
t190 = 0.4e1 * t587;
t183 = t370 * t230 - t318 * t487;
t182 = -t231 * t370 + t318 * t491;
t178 = t367 * t413 + t374 * t452 + t460;
t155 = t399 * t373;
t152 = t538 / 0.2e1;
t142 = -pkin(6) * t586 + t371 * t246 + t247 * t374 + t460;
t141 = (t497 + (t369 * t467 - t372 * t465) * t373) * t370;
t137 = t539 / 0.2e1;
t117 = t540 / 0.2e1;
t80 = t541 + t542;
t79 = t559 / 0.2e1;
t61 = t373 * t579 - t395;
t56 = t118 * t162 + (t207 * t371 + t209 * t374) * t318;
t55 = t171 * t370 + t373 * t402;
t52 = t545 + t560;
t39 = t75 * t371 - t374 * t76;
t38 = t73 * t371 - t374 * t74;
t34 = t117 - t538 / 0.2e1;
t33 = t152 + t117;
t32 = t152 - t540 / 0.2e1;
t30 = t567 / 0.2e1;
t29 = t126 * t370 + (t371 * t88 + t374 * t87) * t373;
t28 = t125 * t370 + (t371 * t86 + t374 * t85) * t373;
t27 = -t373 * t375 + t543 + t548 + t556 + t588;
t20 = (t94 * t371 + t93 * t374 + t171) * t373 + (t132 - t402) * t370;
t19 = t79 + t30 - t539 / 0.2e1;
t18 = t137 + t79 - t567 / 0.2e1;
t17 = t137 + t30 - t559 / 0.2e1;
t10 = t520 + t525;
t7 = m(5) * t56 + t40 * t554 + t41 * t552;
t6 = t443 * t487;
t5 = t566 + (t12 * t550 + t11 * t554 + t55 / 0.2e1) * t373 + (-t527 / 0.2e1 - t528 / 0.2e1 + t20 / 0.2e1) * t370;
t4 = t371 * t386 + t374 * t385;
t3 = t381 + t417;
t2 = (t126 / 0.4e1 + t99 / 0.4e1) * t374 - t562 / 0.2e1 + t387 + (-t125 / 0.4e1 - t98 / 0.4e1) * t371 + t381;
t1 = (-t171 / 0.2e1 + (-t103 / 0.4e1 - t93 / 0.4e1) * t374 + (-t94 / 0.4e1 - t102 / 0.4e1) * t371) * t373 - t565 / 0.2e1 + t387 + (-t132 / 0.2e1 + (t146 / 0.4e1 + t133 / 0.4e1) * t374 + (t148 / 0.4e1 - t134 / 0.4e1) * t371) * t370 + t417;
t8 = [t27 * qJD(2) + t80 * qJD(3) + t61 * qJD(4), t27 * qJD(1) + t10 * qJD(3) + t3 * qJD(4) + ((t166 * t210 + t168 * t208 - t197 * t209 - t198 * t207) * t574 + (t222 * t251 + t223 * t249 - t236 * t250 - t237 * t248) * t575) * t577 + ((m(3) * (-t238 * t340 - t316 * t336) + t571 - t102 / 0.2e1 + t433 * t374 - t385) * qJD(2) + (-t279 / 0.2e1 + t301 / 0.2e1 - t274 / 0.2e1 - t311 / 0.2e1) * t448 + (t277 / 0.2e1 + t303 / 0.2e1 + t270 / 0.2e1 - t314 / 0.2e1) * t449) * t374 + ((m(3) * (-t239 * t340 + t319 * t336) + t572 + t103 / 0.2e1 + t433 * t371 - t386) * qJD(2) + (-t278 / 0.2e1 - t302 / 0.2e1 + t275 / 0.2e1 - t312 / 0.2e1) * t448 + (t276 / 0.2e1 - t304 / 0.2e1 - t271 / 0.2e1 + t315 / 0.2e1) * t449) * t371, qJD(1) * t80 + qJD(2) * t10 + qJD(4) * t33, t61 * qJD(1) + t3 * qJD(2) + t33 * qJD(3) + (t141 + m(5) * (t166 * t182 + t168 * t183 - t180 * t231 - t181 * t230)) * qJD(4) + ((t570 + t125 / 0.2e1 - t443) * t374 + (t99 / 0.2e1 + t126 / 0.2e1) * t371) * t447; t4 * qJD(2) + t595 * qJD(3) + t1 * qJD(4) + (-t556 / 0.4e1 - t543 / 0.4e1 - t548 / 0.4e1) * t578 + t375 * t450 - t588 * qJD(1), t4 * qJD(1) + t52 * qJD(3) + t7 * qJD(4) + (m(3) * ((t371 * (rSges(3,1) * t491 - t453) + t374 * (rSges(3,1) * t487 + t371 * rSges(3,3) - t356)) * (-t371 * t316 - t319 * t374) + t451 * t340 * t336) + m(5) * (t118 * t142 - t207 * t208 - t209 * t210) + m(4) * (t156 * t178 - t248 * t249 - t250 * t251) + ((-t371 * t309 + (t374 * t463 + t426) * t373 + (t374 * t462 + t425) * t370) * t374 + (-t371 * t306 + (t374 * t466 + t428) * t373 + (t374 * t464 + t427) * t370) * t374 + t39 + t583 * t367) * t554 + (t38 + (-t583 * t374 + (t426 + t428 + (t463 + t466) * t374) * t373 + (t425 + t427 + (t462 + t464) * t374) * t370) * t371 + (t306 + t309) * t368) * t552) * qJD(2), t617 + t52 * qJD(2) + t18 * qJD(4) + (-0.4e1 * t587 + 0.2e1 * (t574 + t575) * (-t373 * t586 + t454)) * qJD(3), t1 * qJD(1) + t7 * qJD(2) + t18 * qJD(3) + (-t55 / 0.2e1 + t442 * t374 + t441 * t371) * t447 + (m(5) * (-t155 * t118 - t150 * t162 - t182 * t209 - t183 * t207 + (-t180 * t374 + t181 * t371) * t318) + t28 * t554 + t29 * t552 - t566 + (-t20 / 0.2e1 + (-t99 / 0.2e1 + t614 / 0.2e1) * t374 + (t570 + t573) * t371) * t370) * qJD(4); -t595 * qJD(2) + t32 * qJD(4) + (-t542 / 0.4e1 - t541 / 0.4e1) * t578, -t617 + t190 * qJD(3) + t17 * qJD(4) + 0.4e1 * (-t560 / 0.4e1 - t545 / 0.4e1) * qJD(2) + ((-t373 * t142 + t479) * t574 + (-t373 * t178 + t472) * t575 + ((t208 * t371 + t210 * t374 + t118) * t574 + (t249 * t371 + t251 * t374 + t156) * t575) * t370) * t577, t190 * qJD(2), t32 * qJD(1) + t17 * qJD(2) + m(5) * (t155 * t373 + (t182 * t374 + t183 * t371) * t370) * qJD(4); t395 * qJD(1) + t2 * qJD(2) + t34 * qJD(3) + t6 * qJD(4) - t579 * t450, t2 * qJD(1) + (((t39 / 0.2e1 + t134 / 0.2e1) * t373 + t441) * t374 + ((t38 / 0.2e1 + t133 / 0.2e1) * t373 - t442) * t371 + ((-t610 / 0.2e1 + t571) * t374 + (t107 * t618 + t108 * t550 + t572) * t371) * t370 + (t112 * t118 - t138 * t209 - t139 * t207 - t142 * t150 + t180 * t210 - t181 * t208 - t56) * m(5)) * qJD(2) + t19 * qJD(3) + t5 * qJD(4), qJD(1) * t34 + qJD(2) * t19, t6 * qJD(1) + t5 * qJD(2) + (m(5) * (t150 * t155 + t180 * t182 - t181 * t183) + t141 * t555 + (t28 * t550 + t29 * t554 + (t99 * t371 + t98 * t374) * t555) * t373) * qJD(4);];
Cq = t8;
