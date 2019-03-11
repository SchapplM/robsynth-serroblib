% Calculate kinetic energy for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:32
% EndTime: 2019-03-09 03:27:34
% DurationCPUTime: 2.11s
% Computational Cost: add. (1089->249), mult. (1805->372), div. (0->0), fcn. (1805->8), ass. (0->125)
t494 = Icges(6,1) + Icges(7,1);
t493 = Icges(6,4) - Icges(7,5);
t492 = Icges(7,4) + Icges(6,5);
t491 = Icges(6,2) + Icges(7,3);
t490 = Icges(7,6) - Icges(6,6);
t489 = Icges(6,3) + Icges(7,2);
t488 = rSges(7,1) + pkin(5);
t487 = rSges(7,3) + qJ(6);
t426 = pkin(9) + qJ(5);
t422 = sin(t426);
t423 = cos(t426);
t433 = cos(qJ(1));
t462 = t433 * t423;
t430 = sin(qJ(3));
t431 = sin(qJ(1));
t466 = t430 * t431;
t384 = t422 * t466 - t462;
t385 = t422 * t433 + t423 * t466;
t432 = cos(qJ(3));
t464 = t431 * t432;
t486 = t491 * t384 - t493 * t385 - t490 * t464;
t465 = t430 * t433;
t386 = t422 * t465 + t423 * t431;
t387 = t422 * t431 - t430 * t462;
t463 = t432 * t433;
t485 = -t491 * t386 - t493 * t387 + t490 * t463;
t484 = t490 * t384 + t492 * t385 - t489 * t464;
t483 = -t490 * t386 + t492 * t387 + t489 * t463;
t482 = -t493 * t384 + t494 * t385 - t492 * t464;
t481 = t493 * t386 + t494 * t387 + t492 * t463;
t480 = (t491 * t422 - t493 * t423) * t432 + t490 * t430;
t479 = (t490 * t422 + t492 * t423) * t432 + t489 * t430;
t478 = (-t493 * t422 + t494 * t423) * t432 + t492 * t430;
t428 = cos(pkin(9));
t471 = pkin(4) * t428;
t477 = -pkin(8) * t432 + t430 * t471;
t470 = Icges(4,4) * t430;
t469 = Icges(4,4) * t432;
t427 = sin(pkin(9));
t468 = t427 * t431;
t467 = t427 * t433;
t460 = -rSges(7,2) * t464 + t487 * t384 + t488 * t385;
t459 = rSges(7,2) * t463 - t487 * t386 + t488 * t387;
t458 = rSges(7,2) * t430 + (t487 * t422 + t488 * t423) * t432;
t444 = pkin(3) * t430 - qJ(4) * t432;
t404 = t444 * t433;
t453 = qJD(3) * t433;
t457 = qJD(4) * t430 - t404 * t453;
t414 = pkin(3) * t432 + qJ(4) * t430;
t425 = qJD(2) * t431;
t454 = qJD(3) * t431;
t456 = t414 * t454 + t425;
t408 = qJD(1) * (pkin(1) * t433 + qJ(2) * t431);
t455 = qJD(1) * t433 * pkin(7) + t408;
t452 = qJD(4) * t432;
t451 = qJD(5) * t432;
t412 = pkin(1) * t431 - qJ(2) * t433;
t448 = -pkin(7) * t431 - t412;
t403 = t444 * t431;
t447 = qJD(1) * t403 + t433 * t452 + t455;
t446 = t404 + t448;
t445 = rSges(4,1) * t430 + rSges(4,2) * t432;
t443 = Icges(4,1) * t430 + t469;
t442 = Icges(4,2) * t432 + t470;
t441 = Icges(4,5) * t430 + Icges(4,6) * t432;
t390 = Icges(4,6) * t433 + t431 * t442;
t392 = Icges(4,5) * t433 + t431 * t443;
t440 = -t390 * t432 - t392 * t430;
t391 = Icges(4,6) * t431 - t433 * t442;
t393 = Icges(4,5) * t431 - t433 * t443;
t439 = t391 * t432 + t393 * t430;
t410 = -Icges(4,2) * t430 + t469;
t411 = Icges(4,1) * t432 - t470;
t438 = t410 * t432 + t411 * t430;
t366 = pkin(4) * t467 + t477 * t431;
t367 = pkin(4) * t468 - t477 * t433;
t437 = t367 * t453 + (-t366 - t403) * t454 + t457;
t371 = pkin(8) * t430 + t432 * t471;
t436 = qJD(1) * t366 + (-qJD(2) + (-t371 - t414) * qJD(3)) * t433 + t447;
t435 = t371 * t454 + (-t367 + t446) * qJD(1) - t431 * t452 + t456;
t419 = qJD(5) * t430 + qJD(1);
t416 = rSges(2,1) * t433 - rSges(2,2) * t431;
t415 = rSges(4,1) * t432 - rSges(4,2) * t430;
t413 = rSges(2,1) * t431 + rSges(2,2) * t433;
t409 = Icges(4,5) * t432 - Icges(4,6) * t430;
t407 = -t431 * t451 + t453;
t406 = t433 * t451 + t454;
t402 = -t428 * t465 + t468;
t401 = t427 * t465 + t428 * t431;
t400 = t428 * t466 + t467;
t399 = -t427 * t466 + t428 * t433;
t397 = rSges(4,3) * t431 - t433 * t445;
t396 = rSges(4,3) * t433 + t431 * t445;
t389 = Icges(4,3) * t431 - t433 * t441;
t388 = Icges(4,3) * t433 + t431 * t441;
t383 = rSges(5,3) * t430 + (rSges(5,1) * t428 - rSges(5,2) * t427) * t432;
t382 = Icges(5,5) * t430 + (Icges(5,1) * t428 - Icges(5,4) * t427) * t432;
t381 = Icges(5,6) * t430 + (Icges(5,4) * t428 - Icges(5,2) * t427) * t432;
t380 = Icges(5,3) * t430 + (Icges(5,5) * t428 - Icges(5,6) * t427) * t432;
t379 = rSges(6,3) * t430 + (rSges(6,1) * t423 - rSges(6,2) * t422) * t432;
t370 = t408 - qJD(2) * t433 + qJD(1) * (-rSges(3,2) * t433 + rSges(3,3) * t431);
t369 = t425 + (rSges(3,2) * t431 + rSges(3,3) * t433 - t412) * qJD(1);
t365 = rSges(5,1) * t402 + rSges(5,2) * t401 + rSges(5,3) * t463;
t364 = rSges(5,1) * t400 + rSges(5,2) * t399 - rSges(5,3) * t464;
t363 = Icges(5,1) * t402 + Icges(5,4) * t401 + Icges(5,5) * t463;
t362 = Icges(5,1) * t400 + Icges(5,4) * t399 - Icges(5,5) * t464;
t361 = Icges(5,4) * t402 + Icges(5,2) * t401 + Icges(5,6) * t463;
t360 = Icges(5,4) * t400 + Icges(5,2) * t399 - Icges(5,6) * t464;
t359 = Icges(5,5) * t402 + Icges(5,6) * t401 + Icges(5,3) * t463;
t358 = Icges(5,5) * t400 + Icges(5,6) * t399 - Icges(5,3) * t464;
t353 = (-t396 * t431 + t397 * t433) * qJD(3);
t352 = rSges(6,1) * t387 + rSges(6,2) * t386 + rSges(6,3) * t463;
t350 = rSges(6,1) * t385 - rSges(6,2) * t384 - rSges(6,3) * t464;
t336 = qJD(1) * t396 + (-qJD(3) * t415 - qJD(2)) * t433 + t455;
t335 = t415 * t454 + t425 + (-t397 + t448) * qJD(1);
t334 = qJD(1) * t364 + (-qJD(2) + (-t383 - t414) * qJD(3)) * t433 + t447;
t333 = (qJD(3) * t383 - t452) * t431 + (-t365 + t446) * qJD(1) + t456;
t332 = (t365 * t433 + (-t364 - t403) * t431) * qJD(3) + t457;
t331 = t350 * t419 - t379 * t407 + t436;
t330 = -t352 * t419 + t379 * t406 + t435;
t329 = -t350 * t406 + t352 * t407 + t437;
t328 = -qJD(6) * t386 - t407 * t458 + t419 * t460 + t436;
t327 = qJD(6) * t384 + t406 * t458 - t419 * t459 + t435;
t326 = qJD(6) * t422 * t432 - t406 * t460 + t407 * t459 + t437;
t1 = m(3) * (t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(4) * (t335 ^ 2 + t336 ^ 2 + t353 ^ 2) / 0.2e1 + m(5) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(6) * (t329 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(7) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + ((-t480 * t386 + t478 * t387 + t479 * t463) * t419 + (-t486 * t386 + t482 * t387 + t484 * t463) * t407 + (-t485 * t386 + t481 * t387 + t483 * t463) * t406) * t406 / 0.2e1 + ((t480 * t384 + t478 * t385 - t479 * t464) * t419 + (t486 * t384 + t482 * t385 - t484 * t464) * t407 + (t485 * t384 + t481 * t385 - t483 * t464) * t406) * t407 / 0.2e1 + (((t480 * t422 + t478 * t423) * t419 + (t486 * t422 + t482 * t423) * t407 + (t485 * t422 + t481 * t423) * t406) * t432 + (t483 * t406 + t484 * t407 + t479 * t419) * t430) * t419 / 0.2e1 + (((-t390 * t430 + t392 * t432) * t433 + (-t391 * t430 + t393 * t432) * t431 + (t358 * t433 + t359 * t431) * t430 + ((-t360 * t427 + t362 * t428) * t433 + (-t361 * t427 + t363 * t428) * t431) * t432) * qJD(3) + ((-t381 * t427 + t382 * t428 + t411) * t432 + (-t410 + t380) * t430) * qJD(1)) * qJD(1) / 0.2e1 + (((t358 * t463 + t360 * t401 + t362 * t402 + t440 * t433) * t433 + ((t388 - t439) * t433 + t359 * t463 + t361 * t401 + t363 * t402 + t389 * t431) * t431) * qJD(3) + (t380 * t463 + t381 * t401 + t382 * t402 + t431 * t409 - t433 * t438) * qJD(1)) * t454 / 0.2e1 + (((-t358 * t464 + t360 * t399 + t362 * t400 + t388 * t433) * t433 + ((t389 - t440) * t433 - t359 * t464 + t361 * t399 + t363 * t400 + t439 * t431) * t431) * qJD(3) + (-t380 * t464 + t381 * t399 + t382 * t400 + t433 * t409 + t431 * t438) * qJD(1)) * t453 / 0.2e1 + (Icges(2,3) + Icges(3,1) + m(2) * (t413 ^ 2 + t416 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
