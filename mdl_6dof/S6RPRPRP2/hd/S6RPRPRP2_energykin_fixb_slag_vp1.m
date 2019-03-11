% Calculate kinetic energy for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:04:28
% EndTime: 2019-03-09 03:04:30
% DurationCPUTime: 1.80s
% Computational Cost: add. (1742->214), mult. (1529->321), div. (0->0), fcn. (1472->10), ass. (0->123)
t499 = Icges(6,1) + Icges(7,1);
t498 = -Icges(6,4) + Icges(7,5);
t497 = Icges(7,4) + Icges(6,5);
t496 = Icges(6,2) + Icges(7,3);
t495 = -Icges(7,6) + Icges(6,6);
t494 = Icges(4,3) + Icges(5,3);
t493 = -Icges(6,3) - Icges(7,2);
t411 = qJ(3) + pkin(10);
t407 = sin(t411);
t409 = cos(t411);
t415 = sin(qJ(3));
t418 = cos(qJ(3));
t492 = Icges(4,5) * t418 + Icges(5,5) * t409 - Icges(4,6) * t415 - Icges(5,6) * t407;
t491 = rSges(7,1) + pkin(5);
t490 = rSges(7,3) + qJ(6);
t412 = qJ(1) + pkin(9);
t410 = cos(t412);
t417 = cos(qJ(5));
t454 = t410 * t417;
t408 = sin(t412);
t414 = sin(qJ(5));
t457 = t408 * t414;
t382 = t409 * t457 + t454;
t455 = t410 * t414;
t456 = t408 * t417;
t383 = t409 * t456 - t455;
t459 = t407 * t408;
t489 = t496 * t382 + t498 * t383 - t495 * t459;
t384 = t409 * t455 - t456;
t385 = t409 * t454 + t457;
t458 = t407 * t410;
t488 = t496 * t384 + t498 * t385 - t495 * t458;
t487 = -t495 * t382 + t497 * t383 - t493 * t459;
t486 = -t495 * t384 + t497 * t385 - t493 * t458;
t485 = t498 * t382 + t499 * t383 + t497 * t459;
t484 = t498 * t384 + t499 * t385 + t497 * t458;
t483 = t492 * t408 - t494 * t410;
t482 = t494 * t408 + t492 * t410;
t481 = t495 * t409 + (t496 * t414 + t498 * t417) * t407;
t480 = t493 * t409 + (-t495 * t414 + t497 * t417) * t407;
t479 = -t497 * t409 + (t498 * t414 + t499 * t417) * t407;
t478 = Icges(4,5) * t415 + Icges(5,5) * t407 + Icges(4,6) * t418 + Icges(5,6) * t409;
t461 = Icges(5,4) * t407;
t391 = Icges(5,2) * t409 + t461;
t460 = Icges(5,4) * t409;
t392 = Icges(5,1) * t407 + t460;
t463 = Icges(4,4) * t415;
t397 = Icges(4,2) * t418 + t463;
t462 = Icges(4,4) * t418;
t398 = Icges(4,1) * t415 + t462;
t477 = -t391 * t407 + t392 * t409 - t397 * t415 + t398 * t418;
t432 = -Icges(5,2) * t407 + t460;
t354 = Icges(5,6) * t408 + t410 * t432;
t434 = Icges(5,1) * t409 - t461;
t356 = Icges(5,5) * t408 + t410 * t434;
t433 = -Icges(4,2) * t415 + t462;
t368 = Icges(4,6) * t408 + t410 * t433;
t435 = Icges(4,1) * t418 - t463;
t372 = Icges(4,5) * t408 + t410 * t435;
t476 = -t354 * t407 + t356 * t409 - t368 * t415 + t372 * t418;
t353 = -Icges(5,6) * t410 + t408 * t432;
t355 = -Icges(5,5) * t410 + t408 * t434;
t367 = -Icges(4,6) * t410 + t408 * t433;
t371 = -Icges(4,5) * t410 + t408 * t435;
t475 = t353 * t407 - t355 * t409 + t367 * t415 - t371 * t418;
t416 = sin(qJ(1));
t468 = pkin(1) * t416;
t467 = pkin(3) * t415;
t465 = pkin(3) * t418;
t453 = rSges(7,2) * t459 + t490 * t382 + t383 * t491;
t452 = rSges(7,2) * t458 + t490 * t384 + t385 * t491;
t451 = -rSges(7,2) * t409 + (t490 * t414 + t417 * t491) * t407;
t419 = cos(qJ(1));
t406 = qJD(1) * t419 * pkin(1);
t450 = qJD(1) * (pkin(2) * t410 + pkin(7) * t408) + t406;
t449 = qJD(3) * t408;
t448 = qJD(3) * t410;
t447 = qJD(5) * t407;
t349 = -qJ(4) * t410 + t408 * t465;
t350 = qJ(4) * t408 + t410 * t465;
t446 = t349 * t449 + t350 * t448 + qJD(2);
t443 = -pkin(2) * t408 + pkin(7) * t410 - t468;
t442 = -t349 + t443;
t441 = pkin(4) * t409 + pkin(8) * t407;
t440 = rSges(4,1) * t418 - rSges(4,2) * t415;
t439 = rSges(5,1) * t409 - rSges(5,2) * t407;
t438 = qJD(3) * (-rSges(5,1) * t407 - rSges(5,2) * t409 - t467);
t437 = qJD(3) * (-pkin(4) * t407 + pkin(8) * t409 - t467);
t378 = t441 * t408;
t379 = t441 * t410;
t436 = t378 * t449 + t379 * t448 + t446;
t423 = qJD(1) * t350 - qJD(4) * t410 + t450;
t422 = qJD(1) * t379 + t408 * t437 + t423;
t404 = qJD(4) * t408;
t421 = t404 + (-t378 + t442) * qJD(1) + t410 * t437;
t402 = -qJD(5) * t409 + qJD(1);
t401 = rSges(2,1) * t419 - rSges(2,2) * t416;
t400 = rSges(2,1) * t416 + rSges(2,2) * t419;
t399 = rSges(4,1) * t415 + rSges(4,2) * t418;
t387 = t408 * t447 - t448;
t386 = t410 * t447 + t449;
t381 = t406 + qJD(1) * (rSges(3,1) * t410 - rSges(3,2) * t408);
t380 = (-rSges(3,1) * t408 - rSges(3,2) * t410 - t468) * qJD(1);
t376 = rSges(4,3) * t408 + t410 * t440;
t375 = -rSges(4,3) * t410 + t408 * t440;
t374 = -rSges(6,3) * t409 + (rSges(6,1) * t417 - rSges(6,2) * t414) * t407;
t358 = rSges(5,3) * t408 + t410 * t439;
t357 = -rSges(5,3) * t410 + t408 * t439;
t343 = rSges(6,1) * t385 - rSges(6,2) * t384 + rSges(6,3) * t458;
t341 = rSges(6,1) * t383 - rSges(6,2) * t382 + rSges(6,3) * t459;
t327 = qJD(1) * t376 - t399 * t449 + t450;
t326 = -t399 * t448 + (-t375 + t443) * qJD(1);
t325 = qJD(2) + (t375 * t408 + t376 * t410) * qJD(3);
t324 = qJD(1) * t358 + t408 * t438 + t423;
t323 = t404 + t410 * t438 + (-t357 + t442) * qJD(1);
t322 = (t357 * t408 + t358 * t410) * qJD(3) + t446;
t321 = t343 * t402 - t374 * t386 + t422;
t320 = -t341 * t402 + t374 * t387 + t421;
t319 = t341 * t386 - t343 * t387 + t436;
t318 = qJD(6) * t382 - t386 * t451 + t402 * t452 + t422;
t317 = qJD(6) * t384 + t387 * t451 - t402 * t453 + t421;
t316 = qJD(6) * t407 * t414 + t386 * t453 - t387 * t452 + t436;
t1 = m(7) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(6) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(5) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + m(4) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + ((t481 * t384 + t479 * t385 + t480 * t458) * t402 + (t489 * t384 + t485 * t385 + t487 * t458) * t387 + (t488 * t384 + t484 * t385 + t486 * t458) * t386) * t386 / 0.2e1 + ((t481 * t382 + t479 * t383 + t480 * t459) * t402 + (t489 * t382 + t485 * t383 + t487 * t459) * t387 + (t488 * t382 + t484 * t383 + t486 * t459) * t386) * t387 / 0.2e1 + ((-t486 * t386 - t487 * t387 - t480 * t402) * t409 + ((t481 * t414 + t479 * t417) * t402 + (t489 * t414 + t485 * t417) * t387 + (t488 * t414 + t484 * t417) * t386) * t407) * t402 / 0.2e1 + (((-t353 * t409 - t355 * t407 - t367 * t418 - t371 * t415) * t410 + (t354 * t409 + t356 * t407 + t368 * t418 + t372 * t415) * t408) * qJD(3) + (t409 * t391 + t407 * t392 + t418 * t397 + t415 * t398) * qJD(1)) * qJD(1) / 0.2e1 + ((t482 * t408 ^ 2 + (t475 * t410 + (t476 - t483) * t408) * t410) * qJD(3) + (t478 * t408 + t477 * t410) * qJD(1)) * t449 / 0.2e1 - ((t483 * t410 ^ 2 + (t476 * t408 + (t475 - t482) * t410) * t408) * qJD(3) + (t477 * t408 - t478 * t410) * qJD(1)) * t448 / 0.2e1 + (m(2) * (t400 ^ 2 + t401 ^ 2) + Icges(3,3) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
