% Calculate kinetic energy for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR10_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR10_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:13
% EndTime: 2019-12-05 17:23:15
% DurationCPUTime: 1.75s
% Computational Cost: add. (3336->303), mult. (9155->493), div. (0->0), fcn. (11726->14), ass. (0->139)
t507 = qJD(2) ^ 2;
t506 = cos(qJ(3));
t505 = cos(qJ(4));
t504 = cos(pkin(6));
t503 = sin(pkin(6));
t474 = sin(pkin(11));
t475 = sin(pkin(5));
t502 = t474 * t475;
t476 = cos(pkin(11));
t501 = t475 * t476;
t481 = sin(qJ(2));
t500 = t475 * t481;
t477 = cos(pkin(5));
t499 = t477 * t481;
t483 = cos(qJ(2));
t498 = t477 * t483;
t467 = -t474 * t498 - t476 * t481;
t493 = t475 * t504;
t453 = -t467 * t503 + t474 * t493;
t497 = qJD(2) * t475;
t472 = t474 * t497;
t443 = qJD(3) * t453 + t472;
t492 = t475 * t503;
t464 = t477 * t504 - t483 * t492;
t473 = qJD(2) * t477;
t455 = qJD(3) * t464 + t473;
t468 = -t474 * t499 + t476 * t483;
t480 = sin(qJ(3));
t490 = t506 * t503;
t489 = t475 * t490;
t491 = t504 * t506;
t428 = -t467 * t491 + t468 * t480 - t474 * t489;
t406 = qJD(4) * t428 + t443;
t450 = -t475 * t483 * t491 - t477 * t490 + t480 * t500;
t423 = qJD(4) * t450 + t455;
t495 = t476 * t497;
t465 = -t474 * t481 + t476 * t498;
t452 = -t465 * t503 - t476 * t493;
t466 = t474 * t483 + t476 * t499;
t433 = t466 * pkin(2) + t452 * pkin(8);
t434 = t468 * pkin(2) + t453 * pkin(8);
t494 = t433 * t472 + t434 * t495 + qJD(1);
t444 = qJD(3) * t452 - t495;
t454 = pkin(2) * t500 + t464 * pkin(8);
t488 = t434 * t473 - t454 * t472;
t426 = -t465 * t491 + t466 * t480 + t476 * t489;
t407 = qJD(4) * t426 + t444;
t427 = t466 * t506 + (t504 * t465 - t476 * t492) * t480;
t401 = pkin(3) * t427 + pkin(9) * t426;
t429 = t468 * t506 + (t504 * t467 + t474 * t492) * t480;
t402 = pkin(3) * t429 + pkin(9) * t428;
t487 = t443 * t401 - t402 * t444 + t494;
t486 = (-t433 * t477 - t454 * t501) * qJD(2);
t451 = t477 * t503 * t480 + (t504 * t480 * t483 + t506 * t481) * t475;
t420 = pkin(3) * t451 + pkin(9) * t450;
t485 = t455 * t402 - t420 * t443 + t488;
t484 = -t401 * t455 + t444 * t420 + t486;
t482 = cos(qJ(5));
t479 = sin(qJ(4));
t478 = sin(qJ(5));
t461 = t477 * rSges(3,3) + (rSges(3,1) * t481 + rSges(3,2) * t483) * t475;
t460 = Icges(3,5) * t477 + (Icges(3,1) * t481 + Icges(3,4) * t483) * t475;
t459 = Icges(3,6) * t477 + (Icges(3,4) * t481 + Icges(3,2) * t483) * t475;
t458 = Icges(3,3) * t477 + (Icges(3,5) * t481 + Icges(3,6) * t483) * t475;
t442 = rSges(3,1) * t468 + rSges(3,2) * t467 + rSges(3,3) * t502;
t441 = rSges(3,1) * t466 + rSges(3,2) * t465 - rSges(3,3) * t501;
t440 = Icges(3,1) * t468 + Icges(3,4) * t467 + Icges(3,5) * t502;
t439 = Icges(3,1) * t466 + Icges(3,4) * t465 - Icges(3,5) * t501;
t438 = Icges(3,4) * t468 + Icges(3,2) * t467 + Icges(3,6) * t502;
t437 = Icges(3,4) * t466 + Icges(3,2) * t465 - Icges(3,6) * t501;
t436 = Icges(3,5) * t468 + Icges(3,6) * t467 + Icges(3,3) * t502;
t435 = Icges(3,5) * t466 + Icges(3,6) * t465 - Icges(3,3) * t501;
t432 = t451 * t505 + t464 * t479;
t431 = t451 * t479 - t464 * t505;
t419 = (-t441 * t477 - t461 * t501) * qJD(2);
t418 = (t442 * t477 - t461 * t502) * qJD(2);
t417 = rSges(4,1) * t451 - rSges(4,2) * t450 + rSges(4,3) * t464;
t416 = Icges(4,1) * t451 - Icges(4,4) * t450 + Icges(4,5) * t464;
t415 = Icges(4,4) * t451 - Icges(4,2) * t450 + Icges(4,6) * t464;
t414 = Icges(4,5) * t451 - Icges(4,6) * t450 + Icges(4,3) * t464;
t413 = t429 * t505 + t453 * t479;
t412 = t429 * t479 - t453 * t505;
t411 = t427 * t505 + t452 * t479;
t410 = t427 * t479 - t452 * t505;
t409 = t432 * t482 + t450 * t478;
t408 = -t432 * t478 + t450 * t482;
t405 = qJD(1) + (t441 * t474 + t442 * t476) * t497;
t403 = pkin(4) * t432 + pkin(10) * t431;
t400 = qJD(5) * t431 + t423;
t398 = rSges(5,1) * t432 - rSges(5,2) * t431 + rSges(5,3) * t450;
t397 = rSges(4,1) * t429 - rSges(4,2) * t428 + rSges(4,3) * t453;
t396 = rSges(4,1) * t427 - rSges(4,2) * t426 + rSges(4,3) * t452;
t395 = Icges(5,1) * t432 - Icges(5,4) * t431 + Icges(5,5) * t450;
t394 = Icges(5,4) * t432 - Icges(5,2) * t431 + Icges(5,6) * t450;
t393 = Icges(5,5) * t432 - Icges(5,6) * t431 + Icges(5,3) * t450;
t392 = Icges(4,1) * t429 - Icges(4,4) * t428 + Icges(4,5) * t453;
t391 = Icges(4,1) * t427 - Icges(4,4) * t426 + Icges(4,5) * t452;
t390 = Icges(4,4) * t429 - Icges(4,2) * t428 + Icges(4,6) * t453;
t389 = Icges(4,4) * t427 - Icges(4,2) * t426 + Icges(4,6) * t452;
t388 = Icges(4,5) * t429 - Icges(4,6) * t428 + Icges(4,3) * t453;
t387 = Icges(4,5) * t427 - Icges(4,6) * t426 + Icges(4,3) * t452;
t386 = t413 * t482 + t428 * t478;
t385 = -t413 * t478 + t428 * t482;
t384 = t411 * t482 + t426 * t478;
t383 = -t411 * t478 + t426 * t482;
t381 = pkin(4) * t413 + pkin(10) * t412;
t380 = pkin(4) * t411 + pkin(10) * t410;
t379 = qJD(5) * t410 + t407;
t378 = qJD(5) * t412 + t406;
t377 = rSges(6,1) * t409 + rSges(6,2) * t408 + rSges(6,3) * t431;
t376 = rSges(5,1) * t413 - rSges(5,2) * t412 + rSges(5,3) * t428;
t375 = rSges(5,1) * t411 - rSges(5,2) * t410 + rSges(5,3) * t426;
t374 = Icges(5,1) * t413 - Icges(5,4) * t412 + Icges(5,5) * t428;
t373 = Icges(5,1) * t411 - Icges(5,4) * t410 + Icges(5,5) * t426;
t372 = Icges(6,1) * t409 + Icges(6,4) * t408 + Icges(6,5) * t431;
t371 = Icges(5,4) * t413 - Icges(5,2) * t412 + Icges(5,6) * t428;
t370 = Icges(5,4) * t411 - Icges(5,2) * t410 + Icges(5,6) * t426;
t369 = Icges(6,4) * t409 + Icges(6,2) * t408 + Icges(6,6) * t431;
t368 = Icges(5,5) * t413 - Icges(5,6) * t412 + Icges(5,3) * t428;
t367 = Icges(5,5) * t411 - Icges(5,6) * t410 + Icges(5,3) * t426;
t366 = Icges(6,5) * t409 + Icges(6,6) * t408 + Icges(6,3) * t431;
t365 = -t396 * t455 + t417 * t444 + t486;
t364 = t397 * t455 - t417 * t443 + t488;
t363 = rSges(6,1) * t386 + rSges(6,2) * t385 + rSges(6,3) * t412;
t362 = rSges(6,1) * t384 + rSges(6,2) * t383 + rSges(6,3) * t410;
t361 = Icges(6,1) * t386 + Icges(6,4) * t385 + Icges(6,5) * t412;
t360 = Icges(6,1) * t384 + Icges(6,4) * t383 + Icges(6,5) * t410;
t359 = Icges(6,4) * t386 + Icges(6,2) * t385 + Icges(6,6) * t412;
t358 = Icges(6,4) * t384 + Icges(6,2) * t383 + Icges(6,6) * t410;
t357 = Icges(6,5) * t386 + Icges(6,6) * t385 + Icges(6,3) * t412;
t356 = Icges(6,5) * t384 + Icges(6,6) * t383 + Icges(6,3) * t410;
t355 = t396 * t443 - t397 * t444 + t494;
t354 = -t375 * t423 + t398 * t407 + t484;
t353 = t376 * t423 - t398 * t406 + t485;
t352 = t375 * t406 - t376 * t407 + t487;
t351 = -t362 * t400 + t377 * t379 - t380 * t423 + t403 * t407 + t484;
t350 = t363 * t400 - t377 * t378 + t381 * t423 - t403 * t406 + t485;
t349 = t362 * t378 - t363 * t379 + t380 * t406 - t381 * t407 + t487;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t405 ^ 2 + t418 ^ 2 + t419 ^ 2) / 0.2e1 - t507 * ((-t436 * t501 + t438 * t465 + t440 * t466) * t502 - (-t435 * t501 + t437 * t465 + t439 * t466) * t501 + (-t458 * t501 + t459 * t465 + t460 * t466) * t477) * t501 / 0.2e1 + m(4) * (t355 ^ 2 + t364 ^ 2 + t365 ^ 2) / 0.2e1 + t443 * ((t453 * t388 - t428 * t390 + t429 * t392) * t443 + (t387 * t453 - t389 * t428 + t391 * t429) * t444 + (t414 * t453 - t415 * t428 + t416 * t429) * t455) / 0.2e1 + t444 * ((t388 * t452 - t390 * t426 + t392 * t427) * t443 + (t452 * t387 - t426 * t389 + t427 * t391) * t444 + (t414 * t452 - t415 * t426 + t416 * t427) * t455) / 0.2e1 + t455 * ((t388 * t464 - t390 * t450 + t392 * t451) * t443 + (t387 * t464 - t389 * t450 + t391 * t451) * t444 + (t464 * t414 - t450 * t415 + t451 * t416) * t455) / 0.2e1 + m(5) * (t352 ^ 2 + t353 ^ 2 + t354 ^ 2) / 0.2e1 + t406 * ((t428 * t368 - t412 * t371 + t413 * t374) * t406 + (t367 * t428 - t370 * t412 + t373 * t413) * t407 + (t393 * t428 - t394 * t412 + t395 * t413) * t423) / 0.2e1 + t407 * ((t368 * t426 - t371 * t410 + t374 * t411) * t406 + (t426 * t367 - t410 * t370 + t411 * t373) * t407 + (t393 * t426 - t394 * t410 + t395 * t411) * t423) / 0.2e1 + t423 * ((t368 * t450 - t371 * t431 + t374 * t432) * t406 + (t367 * t450 - t370 * t431 + t373 * t432) * t407 + (t450 * t393 - t431 * t394 + t432 * t395) * t423) / 0.2e1 + m(6) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + t378 * ((t412 * t357 + t385 * t359 + t386 * t361) * t378 + (t356 * t412 + t358 * t385 + t360 * t386) * t379 + (t366 * t412 + t369 * t385 + t372 * t386) * t400) / 0.2e1 + t379 * ((t357 * t410 + t359 * t383 + t361 * t384) * t378 + (t410 * t356 + t383 * t358 + t384 * t360) * t379 + (t366 * t410 + t369 * t383 + t372 * t384) * t400) / 0.2e1 + t400 * ((t357 * t431 + t359 * t408 + t361 * t409) * t378 + (t356 * t431 + t358 * t408 + t360 * t409) * t379 + (t431 * t366 + t408 * t369 + t409 * t372) * t400) / 0.2e1 + (((t436 * t502 + t438 * t467 + t440 * t468) * t502 - (t435 * t502 + t437 * t467 + t439 * t468) * t501 + (t458 * t502 + t459 * t467 + t460 * t468) * t477) * t502 + t477 * (t477 ^ 2 * t458 + (((t438 * t483 + t440 * t481) * t474 - (t437 * t483 + t439 * t481) * t476) * t475 + (-t435 * t476 + t436 * t474 + t459 * t483 + t460 * t481) * t477) * t475)) * t507 / 0.2e1;
T = t1;
