% Calculate kinetic energy for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:02:09
% EndTime: 2019-03-09 18:02:11
% DurationCPUTime: 2.35s
% Computational Cost: add. (1975->285), mult. (1719->444), div. (0->0), fcn. (1564->12), ass. (0->166)
t534 = Icges(4,3) + Icges(5,3);
t454 = qJ(2) + qJ(3);
t444 = pkin(11) + t454;
t439 = sin(t444);
t440 = cos(t444);
t449 = sin(t454);
t450 = cos(t454);
t533 = Icges(4,5) * t450 + Icges(5,5) * t440 - Icges(4,6) * t449 - Icges(5,6) * t439;
t457 = sin(qJ(1));
t460 = cos(qJ(1));
t517 = Icges(5,4) * t440;
t483 = -Icges(5,2) * t439 + t517;
t376 = -Icges(5,6) * t460 + t457 * t483;
t377 = Icges(5,6) * t457 + t460 * t483;
t518 = Icges(5,4) * t439;
t487 = Icges(5,1) * t440 - t518;
t378 = -Icges(5,5) * t460 + t457 * t487;
t379 = Icges(5,5) * t457 + t460 * t487;
t519 = Icges(4,4) * t450;
t484 = -Icges(4,2) * t449 + t519;
t389 = -Icges(4,6) * t460 + t457 * t484;
t390 = Icges(4,6) * t457 + t460 * t484;
t520 = Icges(4,4) * t449;
t488 = Icges(4,1) * t450 - t520;
t391 = -Icges(4,5) * t460 + t457 * t488;
t392 = Icges(4,5) * t457 + t460 * t488;
t415 = Icges(5,2) * t440 + t518;
t416 = Icges(5,1) * t439 + t517;
t422 = Icges(4,2) * t450 + t520;
t423 = Icges(4,1) * t449 + t519;
t448 = qJD(2) * t457;
t428 = qJD(3) * t457 + t448;
t500 = -qJD(2) - qJD(3);
t429 = t500 * t460;
t532 = (-t376 * t439 + t378 * t440 - t389 * t449 + t391 * t450) * t429 + (-t377 * t439 + t379 * t440 - t390 * t449 + t392 * t450) * t428 + (-t415 * t439 + t416 * t440 - t422 * t449 + t423 * t450) * qJD(1);
t531 = (t533 * t457 - t460 * t534) * t429 + (t457 * t534 + t533 * t460) * t428 + (Icges(4,5) * t449 + Icges(5,5) * t439 + Icges(4,6) * t450 + Icges(5,6) * t440) * qJD(1);
t527 = pkin(3) * t449;
t526 = pkin(4) * t439;
t459 = cos(qJ(2));
t524 = t459 * pkin(2);
t456 = sin(qJ(2));
t522 = Icges(3,4) * t456;
t521 = Icges(3,4) * t459;
t442 = qJ(5) + t444;
t437 = sin(t442);
t516 = Icges(6,4) * t437;
t438 = cos(t442);
t515 = Icges(6,4) * t438;
t514 = t437 * t457;
t513 = t437 * t460;
t455 = sin(qJ(6));
t512 = t455 * t457;
t511 = t455 * t460;
t458 = cos(qJ(6));
t510 = t457 * t458;
t509 = t458 * t460;
t383 = -pkin(8) * t460 + t457 * t524;
t384 = pkin(8) * t457 + t460 * t524;
t502 = qJD(2) * t460;
t508 = t383 * t448 + t384 * t502;
t436 = pkin(1) * t457 - pkin(7) * t460;
t507 = -t383 - t436;
t506 = pkin(4) * t440;
t505 = pkin(3) * t450;
t501 = qJD(6) * t437;
t499 = pkin(2) * qJD(2) * t456;
t356 = -qJ(4) * t460 + t457 * t505;
t498 = t428 * t356 + t508;
t497 = -t356 + t507;
t419 = qJD(5) * t457 + t428;
t496 = t460 * t499;
t349 = -pkin(9) * t460 + t457 * t506;
t495 = -t349 + t497;
t494 = pkin(5) * t438 + pkin(10) * t437;
t493 = rSges(3,1) * t459 - rSges(3,2) * t456;
t492 = rSges(4,1) * t450 - rSges(4,2) * t449;
t491 = rSges(5,1) * t440 - rSges(5,2) * t439;
t490 = rSges(6,1) * t438 - rSges(6,2) * t437;
t420 = (-qJD(5) + t500) * t460;
t489 = Icges(3,1) * t459 - t522;
t486 = Icges(6,1) * t438 - t516;
t485 = -Icges(3,2) * t456 + t521;
t482 = -Icges(6,2) * t437 + t515;
t481 = Icges(3,5) * t459 - Icges(3,6) * t456;
t478 = Icges(6,5) * t438 - Icges(6,6) * t437;
t402 = -Icges(3,6) * t460 + t457 * t485;
t404 = -Icges(3,5) * t460 + t457 * t489;
t477 = t402 * t456 - t404 * t459;
t403 = Icges(3,6) * t457 + t460 * t485;
t405 = Icges(3,5) * t457 + t460 * t489;
t476 = -t403 * t456 + t405 * t459;
t431 = Icges(3,2) * t459 + t522;
t432 = Icges(3,1) * t456 + t521;
t475 = -t431 * t456 + t432 * t459;
t427 = qJD(1) * (pkin(1) * t460 + pkin(7) * t457);
t474 = qJD(1) * t384 - t457 * t499 + t427;
t473 = qJD(4) * t457 + t429 * t527 - t496;
t472 = t429 * t526 + t473;
t471 = (Icges(6,5) * t437 + Icges(6,6) * t438) * qJD(1) + (-Icges(6,3) * t460 + t457 * t478) * t420 + (Icges(6,3) * t457 + t460 * t478) * t419;
t350 = pkin(9) * t457 + t460 * t506;
t357 = qJ(4) * t457 + t460 * t505;
t468 = t428 * t349 + (-t350 - t357) * t429 + t498;
t467 = qJD(1) * t357 - qJD(4) * t460 + t474;
t466 = qJD(1) * t350 + (-t526 - t527) * t428 + t467;
t364 = -Icges(6,6) * t460 + t457 * t482;
t365 = Icges(6,6) * t457 + t460 * t482;
t366 = -Icges(6,5) * t460 + t457 * t486;
t367 = Icges(6,5) * t457 + t460 * t486;
t407 = Icges(6,2) * t438 + t516;
t408 = Icges(6,1) * t437 + t515;
t465 = (-t365 * t437 + t367 * t438) * t419 + (-t364 * t437 + t366 * t438) * t420 + (-t407 * t437 + t408 * t438) * qJD(1);
t435 = rSges(2,1) * t460 - rSges(2,2) * t457;
t434 = rSges(2,1) * t457 + rSges(2,2) * t460;
t433 = rSges(3,1) * t456 + rSges(3,2) * t459;
t430 = Icges(3,5) * t456 + Icges(3,6) * t459;
t426 = -qJD(6) * t438 + qJD(1);
t424 = rSges(4,1) * t449 + rSges(4,2) * t450;
t417 = rSges(5,1) * t439 + rSges(5,2) * t440;
t412 = pkin(5) * t437 - pkin(10) * t438;
t411 = rSges(6,1) * t437 + rSges(6,2) * t438;
t410 = rSges(3,3) * t457 + t460 * t493;
t409 = -rSges(3,3) * t460 + t457 * t493;
t401 = Icges(3,3) * t457 + t460 * t481;
t400 = -Icges(3,3) * t460 + t457 * t481;
t398 = t438 * t509 + t512;
t397 = -t438 * t511 + t510;
t396 = t438 * t510 - t511;
t395 = -t438 * t512 - t509;
t394 = rSges(4,3) * t457 + t460 * t492;
t393 = -rSges(4,3) * t460 + t457 * t492;
t386 = t494 * t460;
t385 = t494 * t457;
t381 = rSges(5,3) * t457 + t460 * t491;
t380 = -rSges(5,3) * t460 + t457 * t491;
t373 = t457 * t501 + t420;
t372 = t460 * t501 + t419;
t371 = rSges(6,3) * t457 + t460 * t490;
t370 = -rSges(6,3) * t460 + t457 * t490;
t361 = -rSges(7,3) * t438 + (rSges(7,1) * t458 - rSges(7,2) * t455) * t437;
t360 = -Icges(7,5) * t438 + (Icges(7,1) * t458 - Icges(7,4) * t455) * t437;
t359 = -Icges(7,6) * t438 + (Icges(7,4) * t458 - Icges(7,2) * t455) * t437;
t358 = -Icges(7,3) * t438 + (Icges(7,5) * t458 - Icges(7,6) * t455) * t437;
t354 = qJD(1) * t410 - t433 * t448 + t427;
t353 = -t433 * t502 + (-t409 - t436) * qJD(1);
t352 = (t409 * t457 + t410 * t460) * qJD(2);
t347 = rSges(7,1) * t398 + rSges(7,2) * t397 + rSges(7,3) * t513;
t346 = rSges(7,1) * t396 + rSges(7,2) * t395 + rSges(7,3) * t514;
t345 = Icges(7,1) * t398 + Icges(7,4) * t397 + Icges(7,5) * t513;
t344 = Icges(7,1) * t396 + Icges(7,4) * t395 + Icges(7,5) * t514;
t343 = Icges(7,4) * t398 + Icges(7,2) * t397 + Icges(7,6) * t513;
t342 = Icges(7,4) * t396 + Icges(7,2) * t395 + Icges(7,6) * t514;
t341 = Icges(7,5) * t398 + Icges(7,6) * t397 + Icges(7,3) * t513;
t340 = Icges(7,5) * t396 + Icges(7,6) * t395 + Icges(7,3) * t514;
t338 = qJD(1) * t394 - t424 * t428 + t474;
t337 = -t496 + t424 * t429 + (-t393 + t507) * qJD(1);
t336 = t393 * t428 - t394 * t429 + t508;
t335 = qJD(1) * t381 + (-t417 - t527) * t428 + t467;
t334 = t417 * t429 + (-t380 + t497) * qJD(1) + t473;
t333 = t380 * t428 + (-t357 - t381) * t429 + t498;
t332 = qJD(1) * t371 - t411 * t419 + t466;
t331 = t411 * t420 + (-t370 + t495) * qJD(1) + t472;
t330 = qJD(1) * t386 + t347 * t426 - t361 * t372 - t412 * t419 + t466;
t329 = -t346 * t426 + t361 * t373 + t412 * t420 + (-t385 + t495) * qJD(1) + t472;
t328 = t370 * t419 - t371 * t420 + t468;
t327 = t346 * t372 - t347 * t373 + t385 * t419 - t386 * t420 + t468;
t1 = ((t457 * t430 + t460 * t475) * qJD(1) + (t457 ^ 2 * t401 + (t477 * t460 + (-t400 + t476) * t457) * t460) * qJD(2)) * t448 / 0.2e1 - ((-t460 * t430 + t457 * t475) * qJD(1) + (t460 ^ 2 * t400 + (t476 * t457 + (-t401 + t477) * t460) * t457) * qJD(2)) * t502 / 0.2e1 + t373 * ((t341 * t514 + t343 * t395 + t345 * t396) * t372 + (t340 * t514 + t395 * t342 + t396 * t344) * t373 + (t358 * t514 + t359 * t395 + t360 * t396) * t426) / 0.2e1 + t372 * ((t341 * t513 + t397 * t343 + t398 * t345) * t372 + (t340 * t513 + t342 * t397 + t344 * t398) * t373 + (t358 * t513 + t359 * t397 + t360 * t398) * t426) / 0.2e1 + t426 * ((-t340 * t373 - t341 * t372 - t358 * t426) * t438 + ((-t343 * t455 + t345 * t458) * t372 + (-t342 * t455 + t344 * t458) * t373 + (-t359 * t455 + t360 * t458) * t426) * t437) / 0.2e1 + t419 * (t471 * t457 + t465 * t460) / 0.2e1 + t420 * (t465 * t457 - t471 * t460) / 0.2e1 + m(7) * (t327 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(6) * (t328 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(5) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(4) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(3) * (t352 ^ 2 + t353 ^ 2 + t354 ^ 2) / 0.2e1 + (t531 * t457 + t532 * t460) * t428 / 0.2e1 + (t532 * t457 - t531 * t460) * t429 / 0.2e1 + (Icges(2,3) + m(2) * (t434 ^ 2 + t435 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t365 * t438 + t367 * t437) * t419 + (t364 * t438 + t366 * t437) * t420 + ((t403 * t459 + t405 * t456) * t457 - (t402 * t459 + t404 * t456) * t460) * qJD(2) + (t376 * t440 + t378 * t439 + t389 * t450 + t391 * t449) * t429 + (t377 * t440 + t379 * t439 + t390 * t450 + t392 * t449) * t428 + (t438 * t407 + t437 * t408 + t440 * t415 + t439 * t416 + t450 * t422 + t449 * t423 + t459 * t431 + t456 * t432) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
