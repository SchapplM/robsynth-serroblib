% Calculate kinetic energy for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:45
% EndTime: 2019-03-09 07:04:47
% DurationCPUTime: 1.56s
% Computational Cost: add. (1829->253), mult. (1450->411), div. (0->0), fcn. (1328->12), ass. (0->147)
t433 = pkin(11) + qJ(3);
t425 = qJ(4) + t433;
t419 = sin(t425);
t496 = pkin(4) * t419;
t435 = cos(pkin(11));
t494 = t435 * pkin(2);
t423 = sin(t433);
t493 = Icges(4,4) * t423;
t424 = cos(t433);
t492 = Icges(4,4) * t424;
t491 = Icges(5,4) * t419;
t420 = cos(t425);
t490 = Icges(5,4) * t420;
t422 = qJ(5) + t425;
t416 = sin(t422);
t489 = Icges(6,4) * t416;
t417 = cos(t422);
t488 = Icges(6,4) * t417;
t438 = sin(qJ(1));
t487 = t416 * t438;
t440 = cos(qJ(1));
t486 = t416 * t440;
t437 = sin(qJ(6));
t485 = t437 * t438;
t484 = t437 * t440;
t439 = cos(qJ(6));
t483 = t438 * t439;
t482 = t439 * t440;
t477 = pkin(3) * t424;
t349 = -pkin(8) * t440 + t438 * t477;
t350 = pkin(8) * t438 + t440 * t477;
t428 = qJD(3) * t438;
t474 = qJD(3) * t440;
t480 = t349 * t428 + t350 * t474;
t413 = pkin(1) * t438 - qJ(2) * t440;
t479 = pkin(7) * t440 - t438 * t494 - t413;
t478 = pkin(4) * t420;
t411 = qJD(4) * t438 + t428;
t473 = qJD(6) * t416;
t472 = -qJD(3) - qJD(4);
t471 = pkin(3) * qJD(3) * t423;
t470 = -t349 + t479;
t402 = qJD(5) * t438 + t411;
t344 = -pkin(9) * t440 + t438 * t478;
t469 = -t344 + t470;
t468 = pkin(5) * t417 + pkin(10) * t416;
t410 = qJD(1) * (pkin(1) * t440 + qJ(2) * t438);
t467 = -qJD(2) * t440 + qJD(1) * (pkin(7) * t438 + t440 * t494) + t410;
t434 = sin(pkin(11));
t466 = rSges(3,1) * t435 - rSges(3,2) * t434;
t465 = rSges(4,1) * t424 - rSges(4,2) * t423;
t464 = rSges(5,1) * t420 - rSges(5,2) * t419;
t463 = rSges(6,1) * t417 - rSges(6,2) * t416;
t403 = (-qJD(5) + t472) * t440;
t462 = Icges(4,1) * t424 - t493;
t461 = Icges(5,1) * t420 - t491;
t460 = Icges(6,1) * t417 - t489;
t459 = -Icges(4,2) * t423 + t492;
t458 = -Icges(5,2) * t419 + t490;
t457 = -Icges(6,2) * t416 + t488;
t456 = Icges(4,5) * t424 - Icges(4,6) * t423;
t455 = Icges(5,5) * t420 - Icges(5,6) * t419;
t454 = Icges(6,5) * t417 - Icges(6,6) * t416;
t379 = -Icges(4,6) * t440 + t438 * t459;
t381 = -Icges(4,5) * t440 + t438 * t462;
t453 = t379 * t423 - t381 * t424;
t380 = Icges(4,6) * t438 + t440 * t459;
t382 = Icges(4,5) * t438 + t440 * t462;
t452 = -t380 * t423 + t382 * t424;
t405 = Icges(4,2) * t424 + t493;
t406 = Icges(4,1) * t423 + t492;
t451 = -t405 * t423 + t406 * t424;
t429 = qJD(2) * t438;
t450 = -t440 * t471 + t429;
t345 = pkin(9) * t438 + t440 * t478;
t412 = t472 * t440;
t449 = t411 * t344 - t345 * t412 + t480;
t448 = t412 * t496 + t450;
t447 = (Icges(6,5) * t416 + Icges(6,6) * t417) * qJD(1) + (-Icges(6,3) * t440 + t438 * t454) * t403 + (Icges(6,3) * t438 + t440 * t454) * t402;
t446 = (Icges(5,5) * t419 + Icges(5,6) * t420) * qJD(1) + (-Icges(5,3) * t440 + t438 * t455) * t412 + (Icges(5,3) * t438 + t440 * t455) * t411;
t445 = qJD(1) * t350 - t438 * t471 + t467;
t444 = qJD(1) * t345 - t411 * t496 + t445;
t359 = -Icges(6,6) * t440 + t438 * t457;
t360 = Icges(6,6) * t438 + t440 * t457;
t361 = -Icges(6,5) * t440 + t438 * t460;
t362 = Icges(6,5) * t438 + t440 * t460;
t393 = Icges(6,2) * t417 + t489;
t394 = Icges(6,1) * t416 + t488;
t443 = (-t360 * t416 + t362 * t417) * t402 + (-t359 * t416 + t361 * t417) * t403 + (-t393 * t416 + t394 * t417) * qJD(1);
t370 = -Icges(5,6) * t440 + t438 * t458;
t371 = Icges(5,6) * t438 + t440 * t458;
t372 = -Icges(5,5) * t440 + t438 * t461;
t373 = Icges(5,5) * t438 + t440 * t461;
t399 = Icges(5,2) * t420 + t491;
t400 = Icges(5,1) * t419 + t490;
t442 = (-t371 * t419 + t373 * t420) * t411 + (-t370 * t419 + t372 * t420) * t412 + (-t399 * t419 + t400 * t420) * qJD(1);
t415 = rSges(2,1) * t440 - rSges(2,2) * t438;
t414 = rSges(2,1) * t438 + rSges(2,2) * t440;
t409 = -qJD(6) * t417 + qJD(1);
t407 = rSges(4,1) * t423 + rSges(4,2) * t424;
t404 = Icges(4,5) * t423 + Icges(4,6) * t424;
t401 = rSges(5,1) * t419 + rSges(5,2) * t420;
t396 = pkin(5) * t416 - pkin(10) * t417;
t395 = rSges(6,1) * t416 + rSges(6,2) * t417;
t391 = t417 * t482 + t485;
t390 = -t417 * t484 + t483;
t389 = t417 * t483 - t484;
t388 = -t417 * t485 - t482;
t386 = rSges(4,3) * t438 + t440 * t465;
t385 = -rSges(4,3) * t440 + t438 * t465;
t384 = t468 * t440;
t383 = t468 * t438;
t378 = Icges(4,3) * t438 + t440 * t456;
t377 = -Icges(4,3) * t440 + t438 * t456;
t375 = rSges(5,3) * t438 + t440 * t464;
t374 = -rSges(5,3) * t440 + t438 * t464;
t367 = t438 * t473 + t403;
t366 = t440 * t473 + t402;
t364 = rSges(6,3) * t438 + t440 * t463;
t363 = -rSges(6,3) * t440 + t438 * t463;
t356 = -rSges(7,3) * t417 + (rSges(7,1) * t439 - rSges(7,2) * t437) * t416;
t355 = -Icges(7,5) * t417 + (Icges(7,1) * t439 - Icges(7,4) * t437) * t416;
t354 = -Icges(7,6) * t417 + (Icges(7,4) * t439 - Icges(7,2) * t437) * t416;
t353 = -Icges(7,3) * t417 + (Icges(7,5) * t439 - Icges(7,6) * t437) * t416;
t352 = qJD(1) * t438 * rSges(3,3) + t410 + (qJD(1) * t466 - qJD(2)) * t440;
t351 = t429 + (t440 * rSges(3,3) - t438 * t466 - t413) * qJD(1);
t343 = rSges(7,1) * t391 + rSges(7,2) * t390 + rSges(7,3) * t486;
t342 = rSges(7,1) * t389 + rSges(7,2) * t388 + rSges(7,3) * t487;
t341 = Icges(7,1) * t391 + Icges(7,4) * t390 + Icges(7,5) * t486;
t340 = Icges(7,1) * t389 + Icges(7,4) * t388 + Icges(7,5) * t487;
t339 = Icges(7,4) * t391 + Icges(7,2) * t390 + Icges(7,6) * t486;
t338 = Icges(7,4) * t389 + Icges(7,2) * t388 + Icges(7,6) * t487;
t337 = Icges(7,5) * t391 + Icges(7,6) * t390 + Icges(7,3) * t486;
t336 = Icges(7,5) * t389 + Icges(7,6) * t388 + Icges(7,3) * t487;
t334 = (t385 * t438 + t386 * t440) * qJD(3);
t332 = qJD(1) * t386 - t407 * t428 + t467;
t331 = -t407 * t474 + t429 + (-t385 + t479) * qJD(1);
t330 = qJD(1) * t375 - t401 * t411 + t445;
t329 = t401 * t412 + (-t374 + t470) * qJD(1) + t450;
t328 = t374 * t411 - t375 * t412 + t480;
t327 = qJD(1) * t364 - t395 * t402 + t444;
t326 = t395 * t403 + (-t363 + t469) * qJD(1) + t448;
t325 = t363 * t402 - t364 * t403 + t449;
t324 = qJD(1) * t384 + t343 * t409 - t356 * t366 - t396 * t402 + t444;
t323 = -t342 * t409 + t356 * t367 + t396 * t403 + (-t383 + t469) * qJD(1) + t448;
t322 = t342 * t366 - t343 * t367 + t383 * t402 - t384 * t403 + t449;
t1 = t366 * ((t337 * t486 + t390 * t339 + t391 * t341) * t366 + (t336 * t486 + t338 * t390 + t340 * t391) * t367 + (t353 * t486 + t354 * t390 + t355 * t391) * t409) / 0.2e1 + t409 * ((-t336 * t367 - t337 * t366 - t353 * t409) * t417 + ((-t339 * t437 + t341 * t439) * t366 + (-t338 * t437 + t340 * t439) * t367 + (-t354 * t437 + t355 * t439) * t409) * t416) / 0.2e1 + t367 * ((t337 * t487 + t339 * t388 + t341 * t389) * t366 + (t336 * t487 + t388 * t338 + t389 * t340) * t367 + (t353 * t487 + t354 * t388 + t355 * t389) * t409) / 0.2e1 + t402 * (t438 * t447 + t440 * t443) / 0.2e1 + t403 * (t438 * t443 - t440 * t447) / 0.2e1 + t412 * (t438 * t442 - t440 * t446) / 0.2e1 + t411 * (t438 * t446 + t440 * t442) / 0.2e1 - ((-t440 * t404 + t438 * t451) * qJD(1) + (t440 ^ 2 * t377 + (t452 * t438 + (-t378 + t453) * t440) * t438) * qJD(3)) * t474 / 0.2e1 + m(7) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + m(6) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(5) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(4) * (t331 ^ 2 + t332 ^ 2 + t334 ^ 2) / 0.2e1 + m(3) * (t351 ^ 2 + t352 ^ 2) / 0.2e1 + ((t438 * t404 + t440 * t451) * qJD(1) + (t438 ^ 2 * t378 + (t453 * t440 + (-t377 + t452) * t438) * t440) * qJD(3)) * t428 / 0.2e1 + (m(2) * (t414 ^ 2 + t415 ^ 2) + Icges(3,2) * t435 ^ 2 + (Icges(3,1) * t434 + 0.2e1 * Icges(3,4) * t435) * t434 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t360 * t417 + t362 * t416) * t402 + (t359 * t417 + t361 * t416) * t403 + (t371 * t420 + t373 * t419) * t411 + (t370 * t420 + t372 * t419) * t412 + ((t380 * t424 + t382 * t423) * t438 - (t379 * t424 + t381 * t423) * t440) * qJD(3) + (t417 * t393 + t416 * t394 + t420 * t399 + t419 * t400 + t424 * t405 + t423 * t406) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
