% Calculate kinetic energy for
% S6RPRPRP8
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:18
% EndTime: 2019-03-09 03:24:20
% DurationCPUTime: 1.87s
% Computational Cost: add. (1079->219), mult. (1528->324), div. (0->0), fcn. (1472->8), ass. (0->121)
t510 = Icges(6,1) + Icges(7,1);
t509 = Icges(6,4) - Icges(7,5);
t508 = Icges(7,4) + Icges(6,5);
t507 = Icges(6,2) + Icges(7,3);
t506 = Icges(7,6) - Icges(6,6);
t505 = Icges(4,3) + Icges(5,3);
t504 = Icges(6,3) + Icges(7,2);
t426 = qJ(3) + pkin(9);
t421 = sin(t426);
t422 = cos(t426);
t429 = sin(qJ(3));
t432 = cos(qJ(3));
t503 = Icges(4,5) * t429 + Icges(5,5) * t421 + Icges(4,6) * t432 + Icges(5,6) * t422;
t502 = rSges(7,1) + pkin(5);
t501 = rSges(7,3) + qJ(6);
t431 = cos(qJ(5));
t433 = cos(qJ(1));
t466 = t433 * t431;
t428 = sin(qJ(5));
t430 = sin(qJ(1));
t469 = t428 * t430;
t398 = t421 * t469 - t466;
t467 = t430 * t431;
t468 = t428 * t433;
t399 = t421 * t467 + t468;
t471 = t422 * t430;
t500 = t507 * t398 - t509 * t399 - t506 * t471;
t400 = t421 * t468 + t467;
t401 = -t421 * t466 + t469;
t470 = t422 * t433;
t499 = -t507 * t400 - t509 * t401 + t506 * t470;
t498 = t506 * t398 + t508 * t399 - t504 * t471;
t497 = -t506 * t400 + t508 * t401 + t504 * t470;
t496 = -t509 * t398 + t510 * t399 - t508 * t471;
t495 = t509 * t400 + t510 * t401 + t508 * t470;
t494 = (t507 * t428 - t509 * t431) * t422 + t506 * t421;
t493 = (t506 * t428 + t508 * t431) * t422 + t504 * t421;
t492 = (-t509 * t428 + t510 * t431) * t422 + t508 * t421;
t491 = t505 * t430 - t503 * t433;
t490 = t503 * t430 + t505 * t433;
t489 = Icges(4,5) * t432 + Icges(5,5) * t422 - Icges(4,6) * t429 - Icges(5,6) * t421;
t472 = Icges(5,4) * t422;
t405 = -Icges(5,2) * t421 + t472;
t473 = Icges(5,4) * t421;
t406 = Icges(5,1) * t422 - t473;
t474 = Icges(4,4) * t432;
t411 = -Icges(4,2) * t429 + t474;
t475 = Icges(4,4) * t429;
t412 = Icges(4,1) * t432 - t475;
t488 = t405 * t422 + t406 * t421 + t411 * t432 + t412 * t429;
t446 = Icges(5,2) * t422 + t473;
t375 = Icges(5,6) * t430 - t433 * t446;
t448 = Icges(5,1) * t421 + t472;
t377 = Icges(5,5) * t430 - t433 * t448;
t447 = Icges(4,2) * t432 + t475;
t387 = Icges(4,6) * t430 - t433 * t447;
t449 = Icges(4,1) * t429 + t474;
t389 = Icges(4,5) * t430 - t433 * t449;
t487 = t375 * t422 + t377 * t421 + t387 * t432 + t389 * t429;
t374 = Icges(5,6) * t433 + t430 * t446;
t376 = Icges(5,5) * t433 + t430 * t448;
t386 = Icges(4,6) * t433 + t430 * t447;
t388 = Icges(4,5) * t433 + t430 * t449;
t486 = -t374 * t422 - t376 * t421 - t386 * t432 - t388 * t429;
t479 = pkin(3) * t429;
t478 = pkin(3) * t432;
t465 = -rSges(7,2) * t471 + t501 * t398 + t502 * t399;
t464 = rSges(7,2) * t470 - t501 * t400 + t502 * t401;
t463 = rSges(7,2) * t421 + (t501 * t428 + t502 * t431) * t422;
t409 = qJD(1) * (pkin(1) * t433 + qJ(2) * t430);
t462 = qJD(1) * t433 * pkin(7) + t409;
t461 = qJD(3) * t430;
t460 = qJD(3) * t433;
t459 = qJD(5) * t422;
t425 = qJD(2) * t430;
t458 = qJD(4) * t433 + t461 * t478 + t425;
t413 = pkin(1) * t430 - qJ(2) * t433;
t455 = -pkin(7) * t430 - t413;
t397 = qJ(4) * t433 + t430 * t479;
t454 = qJD(1) * t397 + qJD(4) * t430 + t462;
t396 = qJ(4) * t430 - t433 * t479;
t453 = -t396 + t455;
t452 = pkin(4) * t421 - pkin(8) * t422;
t451 = rSges(4,1) * t429 + rSges(4,2) * t432;
t450 = rSges(5,1) * t421 + rSges(5,2) * t422;
t381 = t396 * t460;
t393 = t452 * t430;
t394 = t452 * t433;
t437 = -t394 * t460 + t381 + (-t393 - t397) * t461;
t408 = pkin(4) * t422 + pkin(8) * t421;
t436 = t408 * t461 + (t394 + t453) * qJD(1) + t458;
t435 = qJD(1) * t393 + (-qJD(2) + (-t408 - t478) * qJD(3)) * t433 + t454;
t417 = qJD(5) * t421 + qJD(1);
t416 = rSges(2,1) * t433 - rSges(2,2) * t430;
t415 = rSges(4,1) * t432 - rSges(4,2) * t429;
t414 = rSges(2,1) * t430 + rSges(2,2) * t433;
t407 = rSges(5,1) * t422 - rSges(5,2) * t421;
t403 = -t430 * t459 + t460;
t402 = t433 * t459 + t461;
t392 = rSges(4,3) * t430 - t433 * t451;
t391 = rSges(4,3) * t433 + t430 * t451;
t379 = rSges(5,3) * t430 - t433 * t450;
t378 = rSges(5,3) * t433 + t430 * t450;
t371 = rSges(6,3) * t421 + (rSges(6,1) * t431 - rSges(6,2) * t428) * t422;
t363 = t409 - qJD(2) * t433 + qJD(1) * (-rSges(3,2) * t433 + rSges(3,3) * t430);
t362 = t425 + (rSges(3,2) * t430 + rSges(3,3) * t433 - t413) * qJD(1);
t359 = (-t391 * t430 + t392 * t433) * qJD(3);
t358 = rSges(6,1) * t401 + rSges(6,2) * t400 + rSges(6,3) * t470;
t356 = rSges(6,1) * t399 - rSges(6,2) * t398 - rSges(6,3) * t471;
t342 = qJD(1) * t391 + (-qJD(3) * t415 - qJD(2)) * t433 + t462;
t341 = t415 * t461 + t425 + (-t392 + t455) * qJD(1);
t340 = qJD(1) * t378 + (-qJD(2) + (-t407 - t478) * qJD(3)) * t433 + t454;
t339 = t407 * t461 + (-t379 + t453) * qJD(1) + t458;
t338 = t381 + (t379 * t433 + (-t378 - t397) * t430) * qJD(3);
t337 = t356 * t417 - t371 * t403 + t435;
t336 = -t358 * t417 + t371 * t402 + t436;
t335 = -t356 * t402 + t358 * t403 + t437;
t334 = -qJD(6) * t400 - t403 * t463 + t417 * t465 + t435;
t333 = qJD(6) * t398 + t402 * t463 - t417 * t464 + t436;
t332 = qJD(6) * t422 * t428 - t402 * t465 + t403 * t464 + t437;
t1 = m(4) * (t341 ^ 2 + t342 ^ 2 + t359 ^ 2) / 0.2e1 + m(5) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(6) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(7) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(3) * (t362 ^ 2 + t363 ^ 2) / 0.2e1 + ((-t400 * t494 + t401 * t492 + t470 * t493) * t417 + (-t400 * t500 + t496 * t401 + t498 * t470) * t403 + (-t499 * t400 + t495 * t401 + t497 * t470) * t402) * t402 / 0.2e1 + ((t398 * t494 + t399 * t492 - t471 * t493) * t417 + (t500 * t398 + t496 * t399 - t498 * t471) * t403 + (t398 * t499 + t399 * t495 - t471 * t497) * t402) * t403 / 0.2e1 + (((t428 * t494 + t431 * t492) * t417 + (t428 * t500 + t496 * t431) * t403 + (t428 * t499 + t431 * t495) * t402) * t422 + (t402 * t497 + t403 * t498 + t417 * t493) * t421) * t417 / 0.2e1 + (((-t374 * t421 + t376 * t422 - t386 * t429 + t388 * t432) * t433 + (-t375 * t421 + t377 * t422 - t387 * t429 + t432 * t389) * t430) * qJD(3) + (-t421 * t405 + t422 * t406 - t429 * t411 + t432 * t412) * qJD(1)) * qJD(1) / 0.2e1 + ((t491 * t430 ^ 2 + (t486 * t433 + (-t487 + t490) * t430) * t433) * qJD(3) + (t430 * t489 - t433 * t488) * qJD(1)) * t461 / 0.2e1 + ((t490 * t433 ^ 2 + (t487 * t430 + (-t486 + t491) * t433) * t430) * qJD(3) + (t430 * t488 + t433 * t489) * qJD(1)) * t460 / 0.2e1 + (Icges(2,3) + Icges(3,1) + m(2) * (t414 ^ 2 + t416 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
