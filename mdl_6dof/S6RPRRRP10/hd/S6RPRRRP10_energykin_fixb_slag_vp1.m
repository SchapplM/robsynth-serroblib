% Calculate kinetic energy for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:16
% EndTime: 2019-03-09 06:30:18
% DurationCPUTime: 2.09s
% Computational Cost: add. (1152->245), mult. (1910->385), div. (0->0), fcn. (1910->8), ass. (0->128)
t494 = Icges(6,1) + Icges(7,1);
t493 = Icges(6,4) - Icges(7,5);
t492 = Icges(7,4) + Icges(6,5);
t491 = Icges(6,2) + Icges(7,3);
t490 = Icges(7,6) - Icges(6,6);
t489 = Icges(6,3) + Icges(7,2);
t488 = rSges(7,1) + pkin(5);
t487 = rSges(7,3) + qJ(6);
t431 = qJ(4) + qJ(5);
t429 = sin(t431);
t430 = cos(t431);
t437 = cos(qJ(1));
t460 = t437 * t430;
t433 = sin(qJ(3));
t434 = sin(qJ(1));
t466 = t433 * t434;
t394 = t429 * t466 - t460;
t395 = t429 * t437 + t430 * t466;
t436 = cos(qJ(3));
t463 = t434 * t436;
t486 = t491 * t394 - t493 * t395 - t490 * t463;
t465 = t433 * t437;
t396 = t429 * t465 + t430 * t434;
t397 = t429 * t434 - t433 * t460;
t461 = t436 * t437;
t485 = -t491 * t396 - t493 * t397 + t490 * t461;
t484 = t490 * t394 + t492 * t395 - t489 * t463;
t483 = -t490 * t396 + t492 * t397 + t489 * t461;
t482 = -t493 * t394 + t494 * t395 - t492 * t463;
t481 = t493 * t396 + t494 * t397 + t492 * t461;
t480 = (t491 * t429 - t493 * t430) * t436 + t490 * t433;
t479 = (t490 * t429 + t492 * t430) * t436 + t489 * t433;
t478 = (-t493 * t429 + t494 * t430) * t436 + t492 * t433;
t435 = cos(qJ(4));
t472 = pkin(4) * t435;
t477 = -pkin(9) * t436 + t433 * t472;
t470 = Icges(4,4) * t433;
t469 = Icges(4,4) * t436;
t432 = sin(qJ(4));
t468 = t432 * t434;
t467 = t432 * t437;
t464 = t434 * t435;
t462 = t435 * t437;
t459 = -rSges(7,2) * t463 + t487 * t394 + t488 * t395;
t458 = rSges(7,2) * t461 - t487 * t396 + t488 * t397;
t457 = rSges(7,2) * t433 + (t487 * t429 + t488 * t430) * t436;
t411 = qJD(1) * (pkin(1) * t437 + qJ(2) * t434);
t456 = qJD(1) * t437 * pkin(7) + t411;
t426 = qJD(3) * t434;
t455 = qJD(4) * t436;
t408 = t437 * t455 + t426;
t427 = qJD(3) * t437;
t422 = qJD(4) * t433 + qJD(1);
t415 = pkin(1) * t434 - qJ(2) * t437;
t454 = -pkin(7) * t434 - t415;
t453 = pkin(3) * t433 - pkin(8) * t436;
t405 = t453 * t434;
t406 = t453 * t437;
t452 = -t405 * t426 - t406 * t427;
t451 = rSges(4,1) * t433 + rSges(4,2) * t436;
t450 = Icges(4,1) * t433 + t469;
t449 = Icges(4,2) * t436 + t470;
t448 = Icges(4,5) * t433 + Icges(4,6) * t436;
t386 = Icges(4,6) * t437 + t434 * t449;
t389 = Icges(4,5) * t437 + t434 * t450;
t447 = -t386 * t436 - t389 * t433;
t387 = Icges(4,6) * t434 - t437 * t449;
t390 = Icges(4,5) * t434 - t437 * t450;
t446 = t387 * t436 + t390 * t433;
t413 = -Icges(4,2) * t433 + t469;
t414 = Icges(4,1) * t436 - t470;
t445 = t413 * t436 + t414 * t433;
t366 = pkin(4) * t467 + t477 * t434;
t367 = pkin(4) * t468 - t477 * t437;
t409 = -t434 * t455 + t427;
t444 = -t366 * t408 + t409 * t367 + t452;
t419 = pkin(3) * t436 + pkin(8) * t433;
t428 = qJD(2) * t434;
t443 = t419 * t426 + t428 + (t406 + t454) * qJD(1);
t442 = qJD(1) * t405 + (-qJD(3) * t419 - qJD(2)) * t437 + t456;
t371 = pkin(9) * t433 + t436 * t472;
t441 = -t367 * t422 + t408 * t371 + t443;
t440 = t422 * t366 - t371 * t409 + t442;
t418 = rSges(2,1) * t437 - rSges(2,2) * t434;
t417 = rSges(4,1) * t436 - rSges(4,2) * t433;
t416 = rSges(2,1) * t434 + rSges(2,2) * t437;
t412 = Icges(4,5) * t436 - Icges(4,6) * t433;
t410 = qJD(5) * t433 + t422;
t404 = -t433 * t462 + t468;
t403 = t432 * t465 + t464;
t402 = t433 * t464 + t467;
t401 = -t432 * t466 + t462;
t393 = rSges(4,3) * t434 - t437 * t451;
t392 = rSges(5,3) * t433 + (rSges(5,1) * t435 - rSges(5,2) * t432) * t436;
t391 = rSges(4,3) * t437 + t434 * t451;
t388 = Icges(5,5) * t433 + (Icges(5,1) * t435 - Icges(5,4) * t432) * t436;
t385 = Icges(5,6) * t433 + (Icges(5,4) * t435 - Icges(5,2) * t432) * t436;
t384 = Icges(4,3) * t434 - t437 * t448;
t383 = Icges(4,3) * t437 + t434 * t448;
t382 = Icges(5,3) * t433 + (Icges(5,5) * t435 - Icges(5,6) * t432) * t436;
t381 = t427 + (-qJD(4) - qJD(5)) * t463;
t380 = qJD(5) * t461 + t408;
t379 = rSges(6,3) * t433 + (rSges(6,1) * t430 - rSges(6,2) * t429) * t436;
t370 = t411 - qJD(2) * t437 + qJD(1) * (-rSges(3,2) * t437 + rSges(3,3) * t434);
t369 = t428 + (rSges(3,2) * t434 + rSges(3,3) * t437 - t415) * qJD(1);
t365 = rSges(5,1) * t404 + rSges(5,2) * t403 + rSges(5,3) * t461;
t364 = rSges(5,1) * t402 + rSges(5,2) * t401 - rSges(5,3) * t463;
t363 = Icges(5,1) * t404 + Icges(5,4) * t403 + Icges(5,5) * t461;
t362 = Icges(5,1) * t402 + Icges(5,4) * t401 - Icges(5,5) * t463;
t361 = Icges(5,4) * t404 + Icges(5,2) * t403 + Icges(5,6) * t461;
t360 = Icges(5,4) * t402 + Icges(5,2) * t401 - Icges(5,6) * t463;
t359 = Icges(5,5) * t404 + Icges(5,6) * t403 + Icges(5,3) * t461;
t358 = Icges(5,5) * t402 + Icges(5,6) * t401 - Icges(5,3) * t463;
t355 = (-t391 * t434 + t393 * t437) * qJD(3);
t353 = rSges(6,1) * t397 + rSges(6,2) * t396 + rSges(6,3) * t461;
t351 = rSges(6,1) * t395 - rSges(6,2) * t394 - rSges(6,3) * t463;
t336 = qJD(1) * t391 + (-qJD(3) * t417 - qJD(2)) * t437 + t456;
t335 = t417 * t426 + t428 + (-t393 + t454) * qJD(1);
t334 = t364 * t422 - t392 * t409 + t442;
t333 = -t365 * t422 + t392 * t408 + t443;
t332 = -t364 * t408 + t365 * t409 + t452;
t331 = t351 * t410 - t379 * t381 + t440;
t330 = -t353 * t410 + t379 * t380 + t441;
t329 = -t351 * t380 + t353 * t381 + t444;
t328 = -qJD(6) * t396 - t381 * t457 + t410 * t459 + t440;
t327 = qJD(6) * t394 + t380 * t457 - t410 * t458 + t441;
t326 = qJD(6) * t429 * t436 - t380 * t459 + t381 * t458 + t444;
t1 = m(3) * (t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(4) * (t335 ^ 2 + t336 ^ 2 + t355 ^ 2) / 0.2e1 + m(5) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(6) * (t329 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(7) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + ((t434 * t412 - t437 * t445) * qJD(1) + (t434 ^ 2 * t384 + (t447 * t437 + (t383 - t446) * t434) * t437) * qJD(3)) * t426 / 0.2e1 + qJD(1) * ((-t433 * t413 + t436 * t414) * qJD(1) + ((-t386 * t433 + t389 * t436) * t437 + (-t387 * t433 + t390 * t436) * t434) * qJD(3)) / 0.2e1 + t409 * ((-t358 * t463 + t401 * t360 + t402 * t362) * t409 + (-t359 * t463 + t361 * t401 + t363 * t402) * t408 + (-t382 * t463 + t385 * t401 + t388 * t402) * t422) / 0.2e1 + t408 * ((t358 * t461 + t360 * t403 + t362 * t404) * t409 + (t359 * t461 + t403 * t361 + t404 * t363) * t408 + (t382 * t461 + t385 * t403 + t388 * t404) * t422) / 0.2e1 + t422 * ((t358 * t409 + t359 * t408 + t382 * t422) * t433 + ((-t360 * t432 + t362 * t435) * t409 + (-t361 * t432 + t363 * t435) * t408 + (-t385 * t432 + t388 * t435) * t422) * t436) / 0.2e1 + ((t437 * t412 + t434 * t445) * qJD(1) + (t437 ^ 2 * t383 + (t446 * t434 + (t384 - t447) * t437) * t434) * qJD(3)) * t427 / 0.2e1 + ((-t480 * t396 + t478 * t397 + t479 * t461) * t410 + (-t486 * t396 + t482 * t397 + t484 * t461) * t381 + (-t485 * t396 + t481 * t397 + t483 * t461) * t380) * t380 / 0.2e1 + ((t480 * t394 + t478 * t395 - t479 * t463) * t410 + (t486 * t394 + t482 * t395 - t484 * t463) * t381 + (t485 * t394 + t481 * t395 - t483 * t463) * t380) * t381 / 0.2e1 + (((t480 * t429 + t478 * t430) * t410 + (t486 * t429 + t482 * t430) * t381 + (t485 * t429 + t481 * t430) * t380) * t436 + (t483 * t380 + t484 * t381 + t479 * t410) * t433) * t410 / 0.2e1 + (m(2) * (t416 ^ 2 + t418 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
