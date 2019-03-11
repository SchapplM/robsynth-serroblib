% Calculate kinetic energy for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:58:45
% EndTime: 2019-03-09 05:58:47
% DurationCPUTime: 2.18s
% Computational Cost: add. (1936->251), mult. (1932->393), div. (0->0), fcn. (1921->10), ass. (0->134)
t504 = Icges(6,1) + Icges(7,1);
t503 = Icges(6,4) + Icges(7,4);
t502 = -Icges(7,5) - Icges(6,5);
t501 = Icges(6,2) + Icges(7,2);
t500 = -Icges(7,6) - Icges(6,6);
t499 = -Icges(7,3) - Icges(6,3);
t430 = qJ(1) + pkin(10);
t424 = sin(t430);
t425 = cos(t430);
t431 = qJ(4) + qJ(5);
t427 = cos(t431);
t426 = sin(t431);
t436 = cos(qJ(3));
t471 = t426 * t436;
t379 = -t424 * t471 - t425 * t427;
t470 = t427 * t436;
t380 = t424 * t470 - t425 * t426;
t433 = sin(qJ(3));
t474 = t424 * t433;
t498 = -t500 * t379 - t502 * t380 - t499 * t474;
t381 = t424 * t427 - t425 * t471;
t382 = t424 * t426 + t425 * t470;
t472 = t425 * t433;
t497 = -t500 * t381 - t502 * t382 - t499 * t472;
t496 = t501 * t379 + t503 * t380 - t500 * t474;
t495 = t501 * t381 + t503 * t382 - t500 * t472;
t494 = t503 * t379 + t504 * t380 - t502 * t474;
t493 = t503 * t381 + t504 * t382 - t502 * t472;
t435 = cos(qJ(4));
t481 = t435 * pkin(4);
t442 = pkin(9) * t433 + t436 * t481;
t432 = sin(qJ(4));
t475 = t424 * t432;
t363 = pkin(4) * t475 + t425 * t442;
t378 = -pkin(9) * t436 + t433 * t481;
t420 = qJD(3) * t424;
t460 = qJD(4) * t433;
t404 = t425 * t460 + t420;
t421 = -qJD(4) * t436 + qJD(1);
t492 = t421 * t363 - t378 * t404;
t473 = t425 * t432;
t362 = -pkin(4) * t473 + t424 * t442;
t461 = qJD(3) * t425;
t405 = t424 * t460 - t461;
t491 = -t362 * t421 + t405 * t378;
t490 = t499 * t436 + (t500 * t426 - t502 * t427) * t433;
t489 = t500 * t436 + (-t501 * t426 + t503 * t427) * t433;
t488 = t502 * t436 + (-t503 * t426 + t504 * t427) * t433;
t434 = sin(qJ(1));
t483 = pkin(1) * t434;
t479 = Icges(4,4) * t433;
t478 = Icges(4,4) * t436;
t469 = t432 * t436;
t468 = t435 * t436;
t463 = pkin(5) * t427;
t441 = qJ(6) * t433 + t436 * t463;
t455 = pkin(5) * t426;
t467 = rSges(7,1) * t380 + rSges(7,2) * t379 + rSges(7,3) * t474 + t424 * t441 - t425 * t455;
t466 = rSges(7,1) * t382 + rSges(7,2) * t381 + rSges(7,3) * t472 + t424 * t455 + t425 * t441;
t465 = (-qJ(6) - rSges(7,3)) * t436 + (rSges(7,1) * t427 - rSges(7,2) * t426 + t463) * t433;
t437 = cos(qJ(1));
t423 = qJD(1) * t437 * pkin(1);
t464 = qJD(1) * (pkin(2) * t425 + pkin(7) * t424) + t423;
t459 = qJD(5) * t433;
t454 = pkin(3) * t436 + pkin(8) * t433;
t403 = t454 * t425;
t458 = qJD(1) * t403 + t464;
t402 = t454 * t424;
t457 = t402 * t420 + t403 * t461 + qJD(2);
t456 = -pkin(2) * t424 + pkin(7) * t425 - t483;
t453 = rSges(4,1) * t436 - rSges(4,2) * t433;
t452 = Icges(4,1) * t436 - t479;
t451 = -Icges(4,2) * t433 + t478;
t450 = Icges(4,5) * t436 - Icges(4,6) * t433;
t370 = -Icges(4,6) * t425 + t424 * t451;
t372 = -Icges(4,5) * t425 + t424 * t452;
t449 = t370 * t433 - t372 * t436;
t371 = Icges(4,6) * t424 + t425 * t451;
t373 = Icges(4,5) * t424 + t425 * t452;
t448 = -t371 * t433 + t373 * t436;
t412 = Icges(4,2) * t436 + t479;
t413 = Icges(4,1) * t433 + t478;
t447 = -t412 * t433 + t413 * t436;
t419 = pkin(3) * t433 - pkin(8) * t436;
t446 = -qJD(3) * t419 + qJD(6) * t433;
t445 = (-t402 + t456) * qJD(1);
t444 = -t419 * t420 + t458;
t443 = t404 * t362 - t363 * t405 + t457;
t440 = -t419 * t461 + t445;
t418 = rSges(2,1) * t437 - rSges(2,2) * t434;
t417 = rSges(2,1) * t434 + rSges(2,2) * t437;
t416 = rSges(4,1) * t433 + rSges(4,2) * t436;
t411 = Icges(4,5) * t433 + Icges(4,6) * t436;
t409 = qJD(1) + (-qJD(4) - qJD(5)) * t436;
t401 = -rSges(5,3) * t436 + (rSges(5,1) * t435 - rSges(5,2) * t432) * t433;
t400 = -Icges(5,5) * t436 + (Icges(5,1) * t435 - Icges(5,4) * t432) * t433;
t399 = -Icges(5,6) * t436 + (Icges(5,4) * t435 - Icges(5,2) * t432) * t433;
t398 = -Icges(5,3) * t436 + (Icges(5,5) * t435 - Icges(5,6) * t432) * t433;
t397 = t425 * t468 + t475;
t396 = t424 * t435 - t425 * t469;
t395 = t424 * t468 - t473;
t394 = -t424 * t469 - t425 * t435;
t392 = t423 + qJD(1) * (rSges(3,1) * t425 - rSges(3,2) * t424);
t391 = (-rSges(3,1) * t424 - rSges(3,2) * t425 - t483) * qJD(1);
t390 = -rSges(6,3) * t436 + (rSges(6,1) * t427 - rSges(6,2) * t426) * t433;
t375 = rSges(4,3) * t424 + t425 * t453;
t374 = -rSges(4,3) * t425 + t424 * t453;
t369 = Icges(4,3) * t424 + t425 * t450;
t368 = -Icges(4,3) * t425 + t424 * t450;
t367 = t424 * t459 + t405;
t366 = t425 * t459 + t404;
t361 = rSges(5,1) * t397 + rSges(5,2) * t396 + rSges(5,3) * t472;
t360 = rSges(5,1) * t395 + rSges(5,2) * t394 + rSges(5,3) * t474;
t359 = Icges(5,1) * t397 + Icges(5,4) * t396 + Icges(5,5) * t472;
t358 = Icges(5,1) * t395 + Icges(5,4) * t394 + Icges(5,5) * t474;
t357 = Icges(5,4) * t397 + Icges(5,2) * t396 + Icges(5,6) * t472;
t356 = Icges(5,4) * t395 + Icges(5,2) * t394 + Icges(5,6) * t474;
t355 = Icges(5,5) * t397 + Icges(5,6) * t396 + Icges(5,3) * t472;
t354 = Icges(5,5) * t395 + Icges(5,6) * t394 + Icges(5,3) * t474;
t352 = rSges(6,1) * t382 + rSges(6,2) * t381 + rSges(6,3) * t472;
t350 = rSges(6,1) * t380 + rSges(6,2) * t379 + rSges(6,3) * t474;
t336 = qJD(1) * t375 - t416 * t420 + t464;
t335 = -t416 * t461 + (-t374 + t456) * qJD(1);
t334 = qJD(2) + (t374 * t424 + t375 * t425) * qJD(3);
t330 = t361 * t421 - t401 * t404 + t444;
t329 = -t360 * t421 + t401 * t405 + t440;
t328 = t360 * t404 - t361 * t405 + t457;
t327 = t352 * t409 - t366 * t390 + t444 + t492;
t326 = -t350 * t409 + t367 * t390 + t440 + t491;
t325 = t350 * t366 - t352 * t367 + t443;
t324 = -t366 * t465 + t409 * t466 + t424 * t446 + t458 + t492;
t323 = t367 * t465 - t409 * t467 + t425 * t446 + t445 + t491;
t322 = -qJD(6) * t436 + t366 * t467 - t367 * t466 + t443;
t1 = t404 * ((t355 * t472 + t396 * t357 + t397 * t359) * t404 + (t354 * t472 + t356 * t396 + t358 * t397) * t405 + (t396 * t399 + t397 * t400 + t398 * t472) * t421) / 0.2e1 + t405 * ((t355 * t474 + t357 * t394 + t359 * t395) * t404 + (t354 * t474 + t394 * t356 + t395 * t358) * t405 + (t394 * t399 + t395 * t400 + t398 * t474) * t421) / 0.2e1 + t421 * ((-t354 * t405 - t355 * t404 - t398 * t421) * t436 + ((-t357 * t432 + t359 * t435) * t404 + (-t356 * t432 + t358 * t435) * t405 + (-t399 * t432 + t400 * t435) * t421) * t433) / 0.2e1 - ((-t425 * t411 + t424 * t447) * qJD(1) + (t425 ^ 2 * t368 + (t448 * t424 + (-t369 + t449) * t425) * t424) * qJD(3)) * t461 / 0.2e1 + qJD(1) * ((t436 * t412 + t433 * t413) * qJD(1) + ((t371 * t436 + t373 * t433) * t424 - (t370 * t436 + t372 * t433) * t425) * qJD(3)) / 0.2e1 + m(6) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(7) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + m(5) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(4) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + ((t424 * t411 + t425 * t447) * qJD(1) + (t424 ^ 2 * t369 + (t449 * t425 + (-t368 + t448) * t424) * t425) * qJD(3)) * t420 / 0.2e1 + ((t489 * t381 + t488 * t382 + t490 * t472) * t409 + (t496 * t381 + t494 * t382 + t498 * t472) * t367 + (t495 * t381 + t493 * t382 + t497 * t472) * t366) * t366 / 0.2e1 + ((t489 * t379 + t488 * t380 + t490 * t474) * t409 + (t496 * t379 + t494 * t380 + t498 * t474) * t367 + (t495 * t379 + t493 * t380 + t497 * t474) * t366) * t367 / 0.2e1 + ((-t497 * t366 - t498 * t367 - t490 * t409) * t436 + ((-t489 * t426 + t488 * t427) * t409 + (-t496 * t426 + t494 * t427) * t367 + (-t495 * t426 + t493 * t427) * t366) * t433) * t409 / 0.2e1 + (Icges(3,3) + Icges(2,3) + m(2) * (t417 ^ 2 + t418 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
