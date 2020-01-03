% Calculate kinetic energy for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR12_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR12_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:33
% EndTime: 2019-12-31 22:46:35
% DurationCPUTime: 2.06s
% Computational Cost: add. (3388->310), mult. (9199->499), div. (0->0), fcn. (11754->14), ass. (0->144)
t480 = sin(pkin(5));
t481 = cos(pkin(5));
t488 = cos(qJ(2));
t489 = cos(qJ(1));
t506 = t488 * t489;
t485 = sin(qJ(2));
t486 = sin(qJ(1));
t509 = t485 * t486;
t466 = t481 * t506 - t509;
t507 = t486 * t488;
t508 = t485 * t489;
t467 = t481 * t508 + t507;
t468 = -t481 * t507 - t508;
t469 = -t481 * t509 + t506;
t510 = t480 * t489;
t511 = t480 * t486;
t497 = (Icges(3,5) * t467 + Icges(3,6) * t466 - Icges(3,3) * t510) * t489 - (Icges(3,5) * t469 + Icges(3,6) * t468 + Icges(3,3) * t511) * t486;
t518 = t480 * t497;
t516 = cos(qJ(3));
t515 = cos(qJ(4));
t514 = cos(pkin(6));
t513 = sin(pkin(6));
t512 = t480 * t485;
t501 = t480 * t514;
t453 = -t466 * t513 - t489 * t501;
t434 = t467 * pkin(2) + pkin(9) * t453;
t454 = -t468 * t513 + t486 * t501;
t435 = t469 * pkin(2) + pkin(9) * t454;
t503 = qJD(2) * t480;
t477 = t486 * t503;
t502 = t489 * t503;
t505 = t434 * t477 + t435 * t502;
t444 = qJD(3) * t454 + t477;
t504 = qJD(1) * (pkin(1) * t486 - pkin(8) * t510);
t478 = qJD(2) * t481 + qJD(1);
t484 = sin(qJ(3));
t498 = t516 * t513;
t496 = t480 * t498;
t499 = t514 * t516;
t432 = -t468 * t499 + t469 * t484 - t486 * t496;
t407 = qJD(4) * t432 + t444;
t500 = t480 * t513;
t465 = t481 * t514 - t488 * t500;
t455 = qJD(3) * t465 + t478;
t451 = -t480 * t488 * t499 - t481 * t498 + t484 * t512;
t422 = qJD(4) * t451 + t455;
t445 = qJD(3) * t453 - t502;
t430 = -t466 * t499 + t467 * t484 + t489 * t496;
t431 = t467 * t516 + (t466 * t514 - t489 * t500) * t484;
t403 = pkin(3) * t431 + pkin(10) * t430;
t433 = t469 * t516 + (t468 * t514 + t486 * t500) * t484;
t404 = pkin(3) * t433 + pkin(10) * t432;
t495 = t444 * t403 - t404 * t445 + t505;
t408 = qJD(4) * t430 + t445;
t456 = pkin(2) * t512 + pkin(9) * t465;
t470 = qJD(1) * (pkin(1) * t489 + pkin(8) * t511);
t494 = t478 * t435 - t456 * t477 + t470;
t493 = -t434 * t478 - t456 * t502 - t504;
t452 = t481 * t513 * t484 + (t484 * t488 * t514 + t485 * t516) * t480;
t421 = pkin(3) * t452 + pkin(10) * t451;
t492 = t455 * t404 - t421 * t444 + t494;
t491 = -t403 * t455 + t445 * t421 + t493;
t487 = cos(qJ(5));
t483 = sin(qJ(4));
t482 = sin(qJ(5));
t475 = rSges(2,1) * t489 - rSges(2,2) * t486;
t474 = rSges(2,1) * t486 + rSges(2,2) * t489;
t462 = rSges(3,3) * t481 + (rSges(3,1) * t485 + rSges(3,2) * t488) * t480;
t461 = Icges(3,5) * t481 + (Icges(3,1) * t485 + Icges(3,4) * t488) * t480;
t460 = Icges(3,6) * t481 + (Icges(3,4) * t485 + Icges(3,2) * t488) * t480;
t459 = Icges(3,3) * t481 + (Icges(3,5) * t485 + Icges(3,6) * t488) * t480;
t443 = rSges(3,1) * t469 + rSges(3,2) * t468 + rSges(3,3) * t511;
t442 = rSges(3,1) * t467 + rSges(3,2) * t466 - rSges(3,3) * t510;
t441 = Icges(3,1) * t469 + Icges(3,4) * t468 + Icges(3,5) * t511;
t440 = Icges(3,1) * t467 + Icges(3,4) * t466 - Icges(3,5) * t510;
t439 = Icges(3,4) * t469 + Icges(3,2) * t468 + Icges(3,6) * t511;
t438 = Icges(3,4) * t467 + Icges(3,2) * t466 - Icges(3,6) * t510;
t429 = t452 * t515 + t465 * t483;
t428 = t452 * t483 - t465 * t515;
t420 = rSges(4,1) * t452 - rSges(4,2) * t451 + rSges(4,3) * t465;
t419 = Icges(4,1) * t452 - Icges(4,4) * t451 + Icges(4,5) * t465;
t418 = Icges(4,4) * t452 - Icges(4,2) * t451 + Icges(4,6) * t465;
t417 = Icges(4,5) * t452 - Icges(4,6) * t451 + Icges(4,3) * t465;
t416 = t433 * t515 + t454 * t483;
t415 = t433 * t483 - t454 * t515;
t414 = t431 * t515 + t453 * t483;
t413 = t431 * t483 - t453 * t515;
t412 = t429 * t487 + t451 * t482;
t411 = -t429 * t482 + t451 * t487;
t410 = t443 * t478 - t462 * t477 + t470;
t409 = -t442 * t478 - t462 * t502 - t504;
t406 = (t442 * t486 + t443 * t489) * t503;
t402 = pkin(4) * t429 + pkin(11) * t428;
t401 = qJD(5) * t428 + t422;
t399 = rSges(4,1) * t433 - rSges(4,2) * t432 + rSges(4,3) * t454;
t398 = rSges(4,1) * t431 - rSges(4,2) * t430 + rSges(4,3) * t453;
t397 = Icges(4,1) * t433 - Icges(4,4) * t432 + Icges(4,5) * t454;
t396 = Icges(4,1) * t431 - Icges(4,4) * t430 + Icges(4,5) * t453;
t395 = Icges(4,4) * t433 - Icges(4,2) * t432 + Icges(4,6) * t454;
t394 = Icges(4,4) * t431 - Icges(4,2) * t430 + Icges(4,6) * t453;
t393 = Icges(4,5) * t433 - Icges(4,6) * t432 + Icges(4,3) * t454;
t392 = Icges(4,5) * t431 - Icges(4,6) * t430 + Icges(4,3) * t453;
t391 = rSges(5,1) * t429 - rSges(5,2) * t428 + rSges(5,3) * t451;
t390 = Icges(5,1) * t429 - Icges(5,4) * t428 + Icges(5,5) * t451;
t389 = Icges(5,4) * t429 - Icges(5,2) * t428 + Icges(5,6) * t451;
t388 = Icges(5,5) * t429 - Icges(5,6) * t428 + Icges(5,3) * t451;
t387 = t416 * t487 + t432 * t482;
t386 = -t416 * t482 + t432 * t487;
t385 = t414 * t487 + t430 * t482;
t384 = -t414 * t482 + t430 * t487;
t382 = pkin(4) * t416 + pkin(11) * t415;
t381 = pkin(4) * t414 + pkin(11) * t413;
t380 = qJD(5) * t413 + t408;
t379 = qJD(5) * t415 + t407;
t378 = rSges(5,1) * t416 - rSges(5,2) * t415 + rSges(5,3) * t432;
t377 = rSges(5,1) * t414 - rSges(5,2) * t413 + rSges(5,3) * t430;
t376 = Icges(5,1) * t416 - Icges(5,4) * t415 + Icges(5,5) * t432;
t375 = Icges(5,1) * t414 - Icges(5,4) * t413 + Icges(5,5) * t430;
t374 = Icges(5,4) * t416 - Icges(5,2) * t415 + Icges(5,6) * t432;
t373 = Icges(5,4) * t414 - Icges(5,2) * t413 + Icges(5,6) * t430;
t372 = Icges(5,5) * t416 - Icges(5,6) * t415 + Icges(5,3) * t432;
t371 = Icges(5,5) * t414 - Icges(5,6) * t413 + Icges(5,3) * t430;
t370 = rSges(6,1) * t412 + rSges(6,2) * t411 + rSges(6,3) * t428;
t369 = Icges(6,1) * t412 + Icges(6,4) * t411 + Icges(6,5) * t428;
t368 = Icges(6,4) * t412 + Icges(6,2) * t411 + Icges(6,6) * t428;
t367 = Icges(6,5) * t412 + Icges(6,6) * t411 + Icges(6,3) * t428;
t366 = rSges(6,1) * t387 + rSges(6,2) * t386 + rSges(6,3) * t415;
t365 = rSges(6,1) * t385 + rSges(6,2) * t384 + rSges(6,3) * t413;
t364 = t399 * t455 - t420 * t444 + t494;
t363 = -t398 * t455 + t420 * t445 + t493;
t362 = Icges(6,1) * t387 + Icges(6,4) * t386 + Icges(6,5) * t415;
t361 = Icges(6,1) * t385 + Icges(6,4) * t384 + Icges(6,5) * t413;
t360 = Icges(6,4) * t387 + Icges(6,2) * t386 + Icges(6,6) * t415;
t359 = Icges(6,4) * t385 + Icges(6,2) * t384 + Icges(6,6) * t413;
t358 = Icges(6,5) * t387 + Icges(6,6) * t386 + Icges(6,3) * t415;
t357 = Icges(6,5) * t385 + Icges(6,6) * t384 + Icges(6,3) * t413;
t356 = t398 * t444 - t399 * t445 + t505;
t355 = t378 * t422 - t391 * t407 + t492;
t354 = -t377 * t422 + t391 * t408 + t491;
t353 = t377 * t407 - t378 * t408 + t495;
t352 = t366 * t401 - t370 * t379 + t382 * t422 - t402 * t407 + t492;
t351 = -t365 * t401 + t370 * t380 - t381 * t422 + t402 * t408 + t491;
t350 = t365 * t379 - t366 * t380 + t381 * t407 - t382 * t408 + t495;
t1 = m(3) * (t406 ^ 2 + t409 ^ 2 + t410 ^ 2) / 0.2e1 + ((t459 * t511 + t460 * t468 + t461 * t469) * t478 + (-(t438 * t468 + t440 * t469) * t489 + (t468 * t439 + t469 * t441 - t518) * t486) * t503) * t477 / 0.2e1 - ((-t459 * t510 + t460 * t466 + t461 * t467) * t478 + ((t439 * t466 + t441 * t467) * t486 + (-t466 * t438 - t467 * t440 + t518) * t489) * t503) * t502 / 0.2e1 + t478 * ((t481 * t459 + (t460 * t488 + t461 * t485) * t480) * t478 + (((t439 * t488 + t441 * t485) * t486 - (t438 * t488 + t440 * t485) * t489) * t480 - t497 * t481) * t503) / 0.2e1 + m(4) * (t356 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + t444 * ((t454 * t393 - t432 * t395 + t433 * t397) * t444 + (t392 * t454 - t394 * t432 + t396 * t433) * t445 + (t417 * t454 - t418 * t432 + t419 * t433) * t455) / 0.2e1 + t445 * ((t393 * t453 - t395 * t430 + t397 * t431) * t444 + (t453 * t392 - t430 * t394 + t431 * t396) * t445 + (t417 * t453 - t418 * t430 + t419 * t431) * t455) / 0.2e1 + t455 * ((t393 * t465 - t395 * t451 + t397 * t452) * t444 + (t392 * t465 - t394 * t451 + t396 * t452) * t445 + (t465 * t417 - t451 * t418 + t452 * t419) * t455) / 0.2e1 + m(5) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + t407 * ((t432 * t372 - t415 * t374 + t416 * t376) * t407 + (t371 * t432 - t373 * t415 + t375 * t416) * t408 + (t388 * t432 - t389 * t415 + t390 * t416) * t422) / 0.2e1 + t408 * ((t372 * t430 - t374 * t413 + t376 * t414) * t407 + (t430 * t371 - t413 * t373 + t414 * t375) * t408 + (t388 * t430 - t389 * t413 + t390 * t414) * t422) / 0.2e1 + t422 * ((t372 * t451 - t374 * t428 + t376 * t429) * t407 + (t371 * t451 - t373 * t428 + t375 * t429) * t408 + (t451 * t388 - t428 * t389 + t429 * t390) * t422) / 0.2e1 + m(6) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + t379 * ((t415 * t358 + t386 * t360 + t387 * t362) * t379 + (t357 * t415 + t359 * t386 + t361 * t387) * t380 + (t367 * t415 + t368 * t386 + t369 * t387) * t401) / 0.2e1 + t380 * ((t358 * t413 + t360 * t384 + t362 * t385) * t379 + (t413 * t357 + t384 * t359 + t385 * t361) * t380 + (t367 * t413 + t368 * t384 + t369 * t385) * t401) / 0.2e1 + t401 * ((t358 * t428 + t360 * t411 + t362 * t412) * t379 + (t357 * t428 + t359 * t411 + t361 * t412) * t380 + (t428 * t367 + t411 * t368 + t412 * t369) * t401) / 0.2e1 + (m(2) * (t474 ^ 2 + t475 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
