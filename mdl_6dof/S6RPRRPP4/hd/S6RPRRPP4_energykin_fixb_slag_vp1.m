% Calculate kinetic energy for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:39
% EndTime: 2019-03-09 04:38:41
% DurationCPUTime: 2.07s
% Computational Cost: add. (1794->243), mult. (1894->361), div. (0->0), fcn. (1894->10), ass. (0->129)
t517 = Icges(6,1) + Icges(7,1);
t516 = -Icges(6,4) + Icges(7,5);
t515 = Icges(7,4) + Icges(6,5);
t514 = Icges(6,2) + Icges(7,3);
t513 = -Icges(7,6) + Icges(6,6);
t512 = -Icges(6,3) - Icges(7,2) - Icges(5,3);
t511 = rSges(7,1) + pkin(5);
t510 = rSges(7,3) + qJ(6);
t445 = pkin(9) + qJ(3);
t442 = cos(t445);
t446 = qJ(4) + pkin(10);
t443 = cos(t446);
t454 = cos(qJ(1));
t487 = t443 * t454;
t441 = sin(t446);
t452 = sin(qJ(1));
t489 = t441 * t452;
t412 = t442 * t489 + t487;
t484 = t452 * t443;
t488 = t441 * t454;
t413 = t442 * t484 - t488;
t440 = sin(t445);
t491 = t440 * t452;
t509 = t514 * t412 + t516 * t413 - t513 * t491;
t414 = t442 * t488 - t484;
t415 = t442 * t487 + t489;
t490 = t440 * t454;
t508 = t514 * t414 + t516 * t415 - t513 * t490;
t507 = t516 * t412 + t517 * t413 + t515 * t491;
t506 = t516 * t414 + t517 * t415 + t515 * t490;
t505 = t513 * t442 + (t514 * t441 + t516 * t443) * t440;
t504 = -t515 * t442 + (t516 * t441 + t517 * t443) * t440;
t453 = cos(qJ(4));
t482 = t453 * t454;
t451 = sin(qJ(4));
t486 = t451 * t452;
t419 = -t442 * t486 - t482;
t483 = t452 * t453;
t485 = t451 * t454;
t420 = t442 * t483 - t485;
t503 = Icges(5,5) * t420 + Icges(5,6) * t419 - t513 * t412 + t515 * t413 - t512 * t491;
t421 = -t442 * t485 + t483;
t422 = t442 * t482 + t486;
t502 = Icges(5,5) * t422 + Icges(5,6) * t421 - t513 * t414 + t515 * t415 - t512 * t490;
t501 = t512 * t442 + (Icges(5,5) * t453 - Icges(5,6) * t451 - t513 * t441 + t515 * t443) * t440;
t448 = cos(pkin(9));
t496 = pkin(2) * t448;
t495 = pkin(4) * t453;
t493 = Icges(4,4) * t440;
t492 = Icges(4,4) * t442;
t480 = rSges(7,2) * t491 + t510 * t412 + t511 * t413;
t479 = rSges(7,2) * t490 + t510 * t414 + t511 * t415;
t478 = -rSges(7,2) * t442 + (t510 * t441 + t511 * t443) * t440;
t433 = pkin(1) * t452 - qJ(2) * t454;
t477 = pkin(7) * t454 - t452 * t496 - t433;
t471 = pkin(3) * t442 + pkin(8) * t440;
t417 = t471 * t452;
t418 = t471 * t454;
t474 = qJD(3) * t454;
t475 = qJD(3) * t452;
t476 = t417 * t475 + t418 * t474;
t473 = qJD(4) * t440;
t472 = qJD(5) * t440;
t430 = qJD(1) * (pkin(1) * t454 + qJ(2) * t452);
t470 = -qJD(2) * t454 + qJD(1) * (pkin(7) * t452 + t454 * t496) + t430;
t447 = sin(pkin(9));
t469 = rSges(3,1) * t448 - rSges(3,2) * t447;
t468 = rSges(4,1) * t442 - rSges(4,2) * t440;
t467 = Icges(4,1) * t442 - t493;
t466 = -Icges(4,2) * t440 + t492;
t465 = Icges(4,5) * t442 - Icges(4,6) * t440;
t403 = -Icges(4,6) * t454 + t452 * t466;
t405 = -Icges(4,5) * t454 + t452 * t467;
t464 = t403 * t440 - t405 * t442;
t404 = Icges(4,6) * t452 + t454 * t466;
t406 = Icges(4,5) * t452 + t454 * t467;
t463 = -t404 * t440 + t406 * t442;
t426 = Icges(4,2) * t442 + t493;
t427 = Icges(4,1) * t440 + t492;
t462 = -t426 * t440 + t427 * t442;
t460 = qJ(5) * t440 + t442 * t495;
t371 = -pkin(4) * t485 + t452 * t460;
t423 = t454 * t473 + t475;
t461 = -qJD(5) * t442 + t423 * t371 + t476;
t429 = pkin(3) * t440 - pkin(8) * t442;
t459 = qJD(1) * t418 - t429 * t475 + t470;
t444 = qJD(2) * t452;
t458 = t444 + (-t417 + t477) * qJD(1) - t429 * t474;
t372 = pkin(4) * t486 + t454 * t460;
t436 = -qJD(4) * t442 + qJD(1);
t457 = t436 * t372 + t452 * t472 + t459;
t386 = -qJ(5) * t442 + t440 * t495;
t424 = t452 * t473 - t474;
t456 = t424 * t386 + t454 * t472 + t458;
t435 = rSges(2,1) * t454 - rSges(2,2) * t452;
t434 = rSges(2,1) * t452 + rSges(2,2) * t454;
t428 = rSges(4,1) * t440 + rSges(4,2) * t442;
t425 = Icges(4,5) * t440 + Icges(4,6) * t442;
t409 = rSges(4,3) * t452 + t454 * t468;
t408 = -rSges(4,3) * t454 + t452 * t468;
t402 = Icges(4,3) * t452 + t454 * t465;
t401 = -Icges(4,3) * t454 + t452 * t465;
t399 = -rSges(5,3) * t442 + (rSges(5,1) * t453 - rSges(5,2) * t451) * t440;
t398 = -Icges(5,5) * t442 + (Icges(5,1) * t453 - Icges(5,4) * t451) * t440;
t397 = -Icges(5,6) * t442 + (Icges(5,4) * t453 - Icges(5,2) * t451) * t440;
t394 = -rSges(6,3) * t442 + (rSges(6,1) * t443 - rSges(6,2) * t441) * t440;
t385 = qJD(1) * t452 * rSges(3,3) + t430 + (qJD(1) * t469 - qJD(2)) * t454;
t384 = t444 + (t454 * rSges(3,3) - t452 * t469 - t433) * qJD(1);
t381 = rSges(5,1) * t422 + rSges(5,2) * t421 + rSges(5,3) * t490;
t380 = rSges(5,1) * t420 + rSges(5,2) * t419 + rSges(5,3) * t491;
t379 = Icges(5,1) * t422 + Icges(5,4) * t421 + Icges(5,5) * t490;
t378 = Icges(5,1) * t420 + Icges(5,4) * t419 + Icges(5,5) * t491;
t377 = Icges(5,4) * t422 + Icges(5,2) * t421 + Icges(5,6) * t490;
t376 = Icges(5,4) * t420 + Icges(5,2) * t419 + Icges(5,6) * t491;
t370 = (t408 * t452 + t409 * t454) * qJD(3);
t369 = rSges(6,1) * t415 - rSges(6,2) * t414 + rSges(6,3) * t490;
t367 = rSges(6,1) * t413 - rSges(6,2) * t412 + rSges(6,3) * t491;
t351 = qJD(1) * t409 - t428 * t475 + t470;
t350 = -t428 * t474 + t444 + (-t408 + t477) * qJD(1);
t349 = t380 * t423 - t381 * t424 + t476;
t348 = t381 * t436 - t399 * t423 + t459;
t347 = -t380 * t436 + t399 * t424 + t458;
t346 = t369 * t436 + (-t386 - t394) * t423 + t457;
t345 = t394 * t424 + (-t367 - t371) * t436 + t456;
t344 = t367 * t423 + (-t369 - t372) * t424 + t461;
t343 = qJD(6) * t412 + t479 * t436 + (-t386 - t478) * t423 + t457;
t342 = qJD(6) * t414 + t478 * t424 + (-t371 - t480) * t436 + t456;
t341 = qJD(6) * t440 * t441 + t480 * t423 + (-t372 - t479) * t424 + t461;
t1 = -((-t454 * t425 + t452 * t462) * qJD(1) + (t454 ^ 2 * t401 + (t463 * t452 + (-t402 + t464) * t454) * t452) * qJD(3)) * t474 / 0.2e1 + qJD(1) * ((t442 * t426 + t440 * t427) * qJD(1) + ((t404 * t442 + t406 * t440) * t452 - (t403 * t442 + t405 * t440) * t454) * qJD(3)) / 0.2e1 + m(3) * (t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(4) * (t350 ^ 2 + t351 ^ 2 + t370 ^ 2) / 0.2e1 + m(5) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(6) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(7) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + ((t452 * t425 + t454 * t462) * qJD(1) + (t452 ^ 2 * t402 + (t464 * t454 + (-t401 + t463) * t452) * t454) * qJD(3)) * t475 / 0.2e1 + ((t421 * t397 + t422 * t398 + t505 * t414 + t504 * t415 + t501 * t490) * t436 + (t421 * t376 + t422 * t378 + t509 * t414 + t507 * t415 + t503 * t490) * t424 + (t421 * t377 + t422 * t379 + t508 * t414 + t506 * t415 + t502 * t490) * t423) * t423 / 0.2e1 + ((t419 * t397 + t420 * t398 + t505 * t412 + t504 * t413 + t501 * t491) * t436 + (t419 * t376 + t420 * t378 + t509 * t412 + t507 * t413 + t503 * t491) * t424 + (t419 * t377 + t420 * t379 + t508 * t412 + t506 * t413 + t502 * t491) * t423) * t424 / 0.2e1 + ((-t502 * t423 - t503 * t424 - t501 * t436) * t442 + ((-t397 * t451 + t398 * t453 + t505 * t441 + t504 * t443) * t436 + (-t376 * t451 + t378 * t453 + t509 * t441 + t507 * t443) * t424 + (-t377 * t451 + t379 * t453 + t508 * t441 + t506 * t443) * t423) * t440) * t436 / 0.2e1 + (Icges(2,3) + Icges(3,2) * t448 ^ 2 + (Icges(3,1) * t447 + 0.2e1 * Icges(3,4) * t448) * t447 + m(2) * (t434 ^ 2 + t435 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
