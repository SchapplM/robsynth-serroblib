% Calculate kinetic energy for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR14_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR14_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR14_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:52
% EndTime: 2019-12-31 19:16:53
% DurationCPUTime: 1.56s
% Computational Cost: add. (3132->266), mult. (8588->426), div. (0->0), fcn. (11106->14), ass. (0->126)
t458 = cos(pkin(5));
t488 = cos(pkin(11));
t490 = sin(qJ(1));
t477 = t490 * t488;
t456 = sin(pkin(11));
t463 = cos(qJ(1));
t484 = t463 * t456;
t470 = t458 * t477 + t484;
t457 = sin(pkin(5));
t487 = sin(pkin(6));
t480 = t457 * t487;
t489 = cos(pkin(6));
t494 = t470 * t489 - t490 * t480;
t492 = cos(qJ(3));
t491 = cos(qJ(4));
t486 = t457 * t456;
t485 = t457 * t463;
t445 = t458 * t484 + t477;
t461 = sin(qJ(3));
t479 = t463 * t488;
t482 = t490 * t456;
t469 = -t458 * t479 + t482;
t465 = t469 * t489;
t478 = t492 * t487;
t423 = t445 * t461 + t465 * t492 + t478 * t485;
t481 = t457 * t489;
t437 = -t463 * t481 + t469 * t487;
t433 = qJD(3) * t437;
t410 = qJD(4) * t423 + t433;
t446 = -t458 * t482 + t479;
t425 = t446 * t461 + t494 * t492;
t438 = t470 * t487 + t481 * t490;
t434 = qJD(3) * t438;
t411 = qJD(4) * t425 + t434;
t444 = t458 * t489 - t480 * t488;
t441 = qJD(3) * t444 + qJD(1);
t483 = t457 * t490;
t475 = t489 * t488;
t435 = -t457 * t475 * t492 - t458 * t478 + t461 * t486;
t427 = qJD(4) * t435 + t441;
t476 = -qJD(2) * t485 + qJD(1) * (t463 * pkin(1) + qJ(2) * t483);
t448 = pkin(1) * t490 - qJ(2) * t485;
t454 = qJD(2) * t483;
t474 = t454 + (-t445 * pkin(2) - pkin(8) * t437 - t448) * qJD(1);
t472 = qJD(1) * (t446 * pkin(2) + pkin(8) * t438) + t476;
t424 = t445 * t492 + (-t463 * t480 - t465) * t461;
t400 = pkin(3) * t424 + pkin(9) * t423;
t426 = t446 * t492 - t494 * t461;
t401 = pkin(3) * t426 + pkin(9) * t425;
t455 = qJD(2) * t458;
t471 = t400 * t434 - t401 * t433 + t455;
t436 = t458 * t487 * t461 + (t456 * t492 + t461 * t475) * t457;
t418 = pkin(3) * t436 + pkin(9) * t435;
t468 = -t400 * t441 + t418 * t433 + t474;
t467 = t441 * t401 - t418 * t434 + t472;
t462 = cos(qJ(5));
t460 = sin(qJ(4));
t459 = sin(qJ(5));
t452 = t463 * rSges(2,1) - rSges(2,2) * t490;
t451 = rSges(2,1) * t490 + t463 * rSges(2,2);
t422 = t436 * t491 + t444 * t460;
t421 = t436 * t460 - t444 * t491;
t417 = qJD(1) * (t446 * rSges(3,1) - rSges(3,2) * t470 + rSges(3,3) * t483) + t476;
t416 = t454 + (-t445 * rSges(3,1) + rSges(3,2) * t469 + rSges(3,3) * t485 - t448) * qJD(1);
t415 = rSges(4,1) * t436 - rSges(4,2) * t435 + rSges(4,3) * t444;
t414 = Icges(4,1) * t436 - Icges(4,4) * t435 + Icges(4,5) * t444;
t413 = Icges(4,4) * t436 - Icges(4,2) * t435 + Icges(4,6) * t444;
t412 = Icges(4,5) * t436 - Icges(4,6) * t435 + Icges(4,3) * t444;
t409 = t426 * t491 + t438 * t460;
t408 = t426 * t460 - t438 * t491;
t407 = t424 * t491 + t437 * t460;
t406 = t424 * t460 - t437 * t491;
t405 = t422 * t462 + t435 * t459;
t404 = -t422 * t459 + t435 * t462;
t402 = qJD(5) * t421 + t427;
t399 = pkin(4) * t422 + pkin(10) * t421;
t396 = rSges(4,1) * t426 - rSges(4,2) * t425 + rSges(4,3) * t438;
t395 = rSges(4,1) * t424 - rSges(4,2) * t423 + rSges(4,3) * t437;
t394 = Icges(4,1) * t426 - Icges(4,4) * t425 + Icges(4,5) * t438;
t393 = Icges(4,1) * t424 - Icges(4,4) * t423 + Icges(4,5) * t437;
t392 = Icges(4,4) * t426 - Icges(4,2) * t425 + Icges(4,6) * t438;
t391 = Icges(4,4) * t424 - Icges(4,2) * t423 + Icges(4,6) * t437;
t390 = Icges(4,5) * t426 - Icges(4,6) * t425 + Icges(4,3) * t438;
t389 = Icges(4,5) * t424 - Icges(4,6) * t423 + Icges(4,3) * t437;
t388 = rSges(5,1) * t422 - rSges(5,2) * t421 + rSges(5,3) * t435;
t387 = Icges(5,1) * t422 - Icges(5,4) * t421 + Icges(5,5) * t435;
t386 = Icges(5,4) * t422 - Icges(5,2) * t421 + Icges(5,6) * t435;
t385 = Icges(5,5) * t422 - Icges(5,6) * t421 + Icges(5,3) * t435;
t384 = t409 * t462 + t425 * t459;
t383 = -t409 * t459 + t425 * t462;
t382 = t407 * t462 + t423 * t459;
t381 = -t407 * t459 + t423 * t462;
t380 = qJD(5) * t408 + t411;
t379 = qJD(5) * t406 + t410;
t378 = pkin(4) * t409 + pkin(10) * t408;
t377 = pkin(4) * t407 + pkin(10) * t406;
t376 = rSges(5,1) * t409 - rSges(5,2) * t408 + rSges(5,3) * t425;
t375 = rSges(5,1) * t407 - rSges(5,2) * t406 + rSges(5,3) * t423;
t374 = Icges(5,1) * t409 - Icges(5,4) * t408 + Icges(5,5) * t425;
t373 = Icges(5,1) * t407 - Icges(5,4) * t406 + Icges(5,5) * t423;
t372 = Icges(5,4) * t409 - Icges(5,2) * t408 + Icges(5,6) * t425;
t371 = Icges(5,4) * t407 - Icges(5,2) * t406 + Icges(5,6) * t423;
t370 = Icges(5,5) * t409 - Icges(5,6) * t408 + Icges(5,3) * t425;
t369 = Icges(5,5) * t407 - Icges(5,6) * t406 + Icges(5,3) * t423;
t368 = rSges(6,1) * t405 + rSges(6,2) * t404 + rSges(6,3) * t421;
t367 = Icges(6,1) * t405 + Icges(6,4) * t404 + Icges(6,5) * t421;
t366 = Icges(6,4) * t405 + Icges(6,2) * t404 + Icges(6,6) * t421;
t365 = Icges(6,5) * t405 + Icges(6,6) * t404 + Icges(6,3) * t421;
t364 = t396 * t441 - t415 * t434 + t472;
t363 = -t395 * t441 + t415 * t433 + t474;
t362 = t455 + (t395 * t438 - t396 * t437) * qJD(3);
t361 = rSges(6,1) * t384 + rSges(6,2) * t383 + rSges(6,3) * t408;
t360 = rSges(6,1) * t382 + rSges(6,2) * t381 + rSges(6,3) * t406;
t359 = Icges(6,1) * t384 + Icges(6,4) * t383 + Icges(6,5) * t408;
t358 = Icges(6,1) * t382 + Icges(6,4) * t381 + Icges(6,5) * t406;
t357 = Icges(6,4) * t384 + Icges(6,2) * t383 + Icges(6,6) * t408;
t356 = Icges(6,4) * t382 + Icges(6,2) * t381 + Icges(6,6) * t406;
t355 = Icges(6,5) * t384 + Icges(6,6) * t383 + Icges(6,3) * t408;
t354 = Icges(6,5) * t382 + Icges(6,6) * t381 + Icges(6,3) * t406;
t353 = t376 * t427 - t388 * t411 + t467;
t352 = -t375 * t427 + t388 * t410 + t468;
t351 = t375 * t411 - t376 * t410 + t471;
t350 = t361 * t402 - t368 * t380 + t378 * t427 - t399 * t411 + t467;
t349 = -t360 * t402 + t368 * t379 - t377 * t427 + t399 * t410 + t468;
t348 = t360 * t380 - t361 * t379 + t377 * t411 - t378 * t410 + t471;
t1 = m(3) * (qJD(2) ^ 2 * t458 ^ 2 + t416 ^ 2 + t417 ^ 2) / 0.2e1 + m(4) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + ((t412 * t438 - t413 * t425 + t414 * t426) * t441 + ((t390 * t438 - t392 * t425 + t394 * t426) * t438 + (t389 * t438 - t391 * t425 + t393 * t426) * t437) * qJD(3)) * t434 / 0.2e1 + ((t412 * t437 - t413 * t423 + t414 * t424) * t441 + ((t390 * t437 - t392 * t423 + t394 * t424) * t438 + (t389 * t437 - t391 * t423 + t393 * t424) * t437) * qJD(3)) * t433 / 0.2e1 + t441 * ((t412 * t444 - t435 * t413 + t414 * t436) * t441 + ((t390 * t444 - t392 * t435 + t394 * t436) * t438 + (t389 * t444 - t391 * t435 + t393 * t436) * t437) * qJD(3)) / 0.2e1 + m(5) * (t351 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + t411 * ((t425 * t370 - t372 * t408 + t409 * t374) * t411 + (t369 * t425 - t371 * t408 + t373 * t409) * t410 + (t385 * t425 - t386 * t408 + t387 * t409) * t427) / 0.2e1 + t410 * ((t370 * t423 - t372 * t406 + t374 * t407) * t411 + (t369 * t423 - t371 * t406 + t373 * t407) * t410 + (t385 * t423 - t386 * t406 + t387 * t407) * t427) / 0.2e1 + t427 * ((t370 * t435 - t372 * t421 + t374 * t422) * t411 + (t369 * t435 - t371 * t421 + t373 * t422) * t410 + (t385 * t435 - t386 * t421 + t387 * t422) * t427) / 0.2e1 + m(6) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + t380 * ((t355 * t408 + t357 * t383 + t359 * t384) * t380 + (t354 * t408 + t356 * t383 + t358 * t384) * t379 + (t365 * t408 + t366 * t383 + t367 * t384) * t402) / 0.2e1 + t379 * ((t355 * t406 + t357 * t381 + t359 * t382) * t380 + (t354 * t406 + t356 * t381 + t358 * t382) * t379 + (t365 * t406 + t366 * t381 + t367 * t382) * t402) / 0.2e1 + t402 * ((t355 * t421 + t357 * t404 + t359 * t405) * t380 + (t354 * t421 + t356 * t404 + t358 * t405) * t379 + (t365 * t421 + t366 * t404 + t367 * t405) * t402) / 0.2e1 + (m(2) * (t451 ^ 2 + t452 ^ 2) + Icges(2,3) + (Icges(3,5) * t458 + (Icges(3,1) * t456 + Icges(3,4) * t488) * t457) * t486 + t457 * t488 * (Icges(3,6) * t458 + (Icges(3,4) * t456 + Icges(3,2) * t488) * t457) + t458 * (Icges(3,3) * t458 + (Icges(3,5) * t456 + Icges(3,6) * t488) * t457)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
