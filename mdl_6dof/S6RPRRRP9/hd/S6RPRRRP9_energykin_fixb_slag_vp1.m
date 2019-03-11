% Calculate kinetic energy for
% S6RPRRRP9
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:07
% EndTime: 2019-03-09 06:26:09
% DurationCPUTime: 2.11s
% Computational Cost: add. (1187->252), mult. (1943->393), div. (0->0), fcn. (1933->8), ass. (0->131)
t507 = Icges(6,1) + Icges(7,1);
t506 = Icges(6,4) + Icges(7,4);
t505 = Icges(7,5) + Icges(6,5);
t504 = Icges(6,2) + Icges(7,2);
t503 = Icges(7,6) + Icges(6,6);
t502 = Icges(7,3) + Icges(6,3);
t438 = qJ(4) + qJ(5);
t434 = sin(t438);
t435 = cos(t438);
t444 = cos(qJ(1));
t440 = sin(qJ(3));
t441 = sin(qJ(1));
t477 = t440 * t441;
t399 = -t434 * t477 + t435 * t444;
t400 = t434 * t444 + t435 * t477;
t443 = cos(qJ(3));
t474 = t441 * t443;
t501 = t503 * t399 + t505 * t400 - t502 * t474;
t476 = t440 * t444;
t401 = t434 * t476 + t435 * t441;
t402 = t434 * t441 - t435 * t476;
t472 = t443 * t444;
t500 = t503 * t401 + t505 * t402 + t502 * t472;
t499 = t504 * t399 + t506 * t400 - t503 * t474;
t498 = t504 * t401 + t506 * t402 + t503 * t472;
t497 = t506 * t399 + t507 * t400 - t505 * t474;
t496 = t506 * t401 + t507 * t402 + t505 * t472;
t439 = sin(qJ(4));
t478 = t439 * t444;
t442 = cos(qJ(4));
t484 = t442 * pkin(4);
t490 = -pkin(9) * t443 + t440 * t484;
t370 = pkin(4) * t478 + t490 * t441;
t376 = pkin(9) * t440 + t443 * t484;
t432 = qJD(3) * t444;
t465 = qJD(4) * t443;
t413 = -t441 * t465 + t432;
t427 = qJD(4) * t440 + qJD(1);
t495 = t427 * t370 - t376 * t413;
t494 = (-t503 * t434 + t505 * t435) * t443 + t502 * t440;
t493 = (-t504 * t434 + t506 * t435) * t443 + t503 * t440;
t492 = (-t506 * t434 + t507 * t435) * t443 + t505 * t440;
t468 = pkin(5) * t435;
t491 = -qJ(6) * t443 + t440 * t468;
t482 = Icges(4,4) * t440;
t481 = Icges(4,4) * t443;
t479 = t439 * t441;
t475 = t441 * t442;
t473 = t442 * t444;
t462 = pkin(5) * t434;
t471 = rSges(7,1) * t400 + rSges(7,2) * t399 - rSges(7,3) * t474 + t491 * t441 + t462 * t444;
t470 = rSges(7,1) * t402 + rSges(7,2) * t401 + rSges(7,3) * t472 + t462 * t441 - t491 * t444;
t469 = (rSges(7,1) * t435 - rSges(7,2) * t434 + t468) * t443 + (qJ(6) + rSges(7,3)) * t440;
t416 = qJD(1) * (pkin(1) * t444 + qJ(2) * t441);
t467 = qJD(1) * t444 * pkin(7) + t416;
t431 = qJD(3) * t441;
t412 = t444 * t465 + t431;
t464 = qJD(6) * t443;
t459 = pkin(3) * t440 - pkin(8) * t443;
t409 = t459 * t441;
t463 = qJD(1) * t409 + t467;
t421 = pkin(1) * t441 - qJ(2) * t444;
t461 = -pkin(7) * t441 - t421;
t425 = pkin(3) * t443 + pkin(8) * t440;
t460 = -qJD(3) * t425 - qJD(2);
t410 = t459 * t444;
t458 = -t409 * t431 - t410 * t432;
t457 = rSges(4,1) * t440 + rSges(4,2) * t443;
t456 = Icges(4,1) * t440 + t481;
t455 = Icges(4,2) * t443 + t482;
t454 = Icges(4,5) * t440 + Icges(4,6) * t443;
t391 = Icges(4,6) * t444 + t441 * t455;
t394 = Icges(4,5) * t444 + t441 * t456;
t453 = -t391 * t443 - t394 * t440;
t392 = Icges(4,6) * t441 - t444 * t455;
t395 = Icges(4,5) * t441 - t444 * t456;
t452 = t392 * t443 + t395 * t440;
t419 = -Icges(4,2) * t440 + t481;
t420 = Icges(4,1) * t443 - t482;
t451 = t419 * t443 + t420 * t440;
t371 = pkin(4) * t479 - t490 * t444;
t450 = -t370 * t412 + t413 * t371 + t458;
t433 = qJD(2) * t441;
t449 = t425 * t431 + t433 + (t410 + t461) * qJD(1);
t448 = t444 * t460 + t463;
t447 = -t371 * t427 + t412 * t376 + t449;
t424 = rSges(2,1) * t444 - rSges(2,2) * t441;
t423 = rSges(4,1) * t443 - rSges(4,2) * t440;
t422 = rSges(2,1) * t441 + rSges(2,2) * t444;
t418 = Icges(4,5) * t443 - Icges(4,6) * t440;
t415 = qJD(5) * t440 + t427;
t408 = -t440 * t473 + t479;
t407 = t439 * t476 + t475;
t406 = t440 * t475 + t478;
t405 = -t439 * t477 + t473;
t398 = rSges(4,3) * t441 - t444 * t457;
t397 = rSges(5,3) * t440 + (rSges(5,1) * t442 - rSges(5,2) * t439) * t443;
t396 = rSges(4,3) * t444 + t441 * t457;
t393 = Icges(5,5) * t440 + (Icges(5,1) * t442 - Icges(5,4) * t439) * t443;
t390 = Icges(5,6) * t440 + (Icges(5,4) * t442 - Icges(5,2) * t439) * t443;
t389 = Icges(4,3) * t441 - t444 * t454;
t388 = Icges(4,3) * t444 + t441 * t454;
t387 = Icges(5,3) * t440 + (Icges(5,5) * t442 - Icges(5,6) * t439) * t443;
t386 = t432 + (-qJD(4) - qJD(5)) * t474;
t385 = qJD(5) * t472 + t412;
t384 = rSges(6,3) * t440 + (rSges(6,1) * t435 - rSges(6,2) * t434) * t443;
t375 = t416 - qJD(2) * t444 + qJD(1) * (-rSges(3,2) * t444 + rSges(3,3) * t441);
t374 = t433 + (rSges(3,2) * t441 + rSges(3,3) * t444 - t421) * qJD(1);
t369 = rSges(5,1) * t408 + rSges(5,2) * t407 + rSges(5,3) * t472;
t368 = rSges(5,1) * t406 + rSges(5,2) * t405 - rSges(5,3) * t474;
t367 = Icges(5,1) * t408 + Icges(5,4) * t407 + Icges(5,5) * t472;
t366 = Icges(5,1) * t406 + Icges(5,4) * t405 - Icges(5,5) * t474;
t365 = Icges(5,4) * t408 + Icges(5,2) * t407 + Icges(5,6) * t472;
t364 = Icges(5,4) * t406 + Icges(5,2) * t405 - Icges(5,6) * t474;
t363 = Icges(5,5) * t408 + Icges(5,6) * t407 + Icges(5,3) * t472;
t362 = Icges(5,5) * t406 + Icges(5,6) * t405 - Icges(5,3) * t474;
t361 = (-t396 * t441 + t398 * t444) * qJD(3);
t359 = rSges(6,1) * t402 + rSges(6,2) * t401 + rSges(6,3) * t472;
t357 = rSges(6,1) * t400 + rSges(6,2) * t399 - rSges(6,3) * t474;
t342 = qJD(1) * t396 + (-qJD(3) * t423 - qJD(2)) * t444 + t467;
t341 = t423 * t431 + t433 + (-t398 + t461) * qJD(1);
t338 = t368 * t427 - t397 * t413 + t448;
t337 = -t369 * t427 + t397 * t412 + t449;
t336 = -t368 * t412 + t369 * t413 + t458;
t335 = t357 * t415 - t384 * t386 + t448 + t495;
t334 = -t359 * t415 + t384 * t385 + t447;
t333 = -t357 * t385 + t359 * t386 + t450;
t332 = t471 * t415 - t469 * t386 + (t460 + t464) * t444 + t463 + t495;
t331 = t385 * t469 - t415 * t470 - t441 * t464 + t447;
t330 = qJD(6) * t440 - t385 * t471 + t386 * t470 + t450;
t1 = ((t444 * t418 + t441 * t451) * qJD(1) + (t444 ^ 2 * t388 + (t452 * t441 + (t389 - t453) * t444) * t441) * qJD(3)) * t432 / 0.2e1 + m(3) * (t374 ^ 2 + t375 ^ 2) / 0.2e1 + t413 * ((-t362 * t474 + t405 * t364 + t406 * t366) * t413 + (-t363 * t474 + t365 * t405 + t367 * t406) * t412 + (-t387 * t474 + t390 * t405 + t393 * t406) * t427) / 0.2e1 + t412 * ((t362 * t472 + t364 * t407 + t366 * t408) * t413 + (t363 * t472 + t407 * t365 + t408 * t367) * t412 + (t387 * t472 + t390 * t407 + t393 * t408) * t427) / 0.2e1 + t427 * ((t362 * t413 + t363 * t412 + t387 * t427) * t440 + ((-t364 * t439 + t366 * t442) * t413 + (-t365 * t439 + t367 * t442) * t412 + (-t390 * t439 + t393 * t442) * t427) * t443) / 0.2e1 + ((t441 * t418 - t444 * t451) * qJD(1) + (t441 ^ 2 * t389 + (t453 * t444 + (t388 - t452) * t441) * t444) * qJD(3)) * t431 / 0.2e1 + qJD(1) * ((-t440 * t419 + t443 * t420) * qJD(1) + ((-t391 * t440 + t394 * t443) * t444 + (-t440 * t392 + t443 * t395) * t441) * qJD(3)) / 0.2e1 + m(6) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(7) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(5) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(4) * (t341 ^ 2 + t342 ^ 2 + t361 ^ 2) / 0.2e1 + ((t493 * t401 + t492 * t402 + t494 * t472) * t415 + (t499 * t401 + t497 * t402 + t501 * t472) * t386 + (t498 * t401 + t496 * t402 + t500 * t472) * t385) * t385 / 0.2e1 + ((t493 * t399 + t492 * t400 - t494 * t474) * t415 + (t499 * t399 + t497 * t400 - t501 * t474) * t386 + (t498 * t399 + t496 * t400 - t500 * t474) * t385) * t386 / 0.2e1 + (((-t493 * t434 + t492 * t435) * t415 + (-t499 * t434 + t497 * t435) * t386 + (-t498 * t434 + t496 * t435) * t385) * t443 + (t500 * t385 + t501 * t386 + t494 * t415) * t440) * t415 / 0.2e1 + (m(2) * (t422 ^ 2 + t424 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
