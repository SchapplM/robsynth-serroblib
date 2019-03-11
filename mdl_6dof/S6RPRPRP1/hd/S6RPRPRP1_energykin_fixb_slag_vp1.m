% Calculate kinetic energy for
% S6RPRPRP1
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:18
% EndTime: 2019-03-09 03:01:20
% DurationCPUTime: 1.94s
% Computational Cost: add. (1774->220), mult. (1538->328), div. (0->0), fcn. (1471->10), ass. (0->124)
t503 = Icges(6,1) + Icges(7,1);
t502 = Icges(6,4) + Icges(7,4);
t501 = -Icges(7,5) - Icges(6,5);
t500 = Icges(6,2) + Icges(7,2);
t499 = -Icges(7,6) - Icges(6,6);
t498 = Icges(4,3) + Icges(5,3);
t497 = -Icges(7,3) - Icges(6,3);
t412 = qJ(3) + pkin(10);
t408 = sin(t412);
t410 = cos(t412);
t417 = sin(qJ(3));
t420 = cos(qJ(3));
t496 = Icges(4,5) * t420 + Icges(5,5) * t410 - Icges(4,6) * t417 - Icges(5,6) * t408;
t413 = qJ(1) + pkin(9);
t411 = cos(t413);
t419 = cos(qJ(5));
t458 = t411 * t419;
t409 = sin(t413);
t416 = sin(qJ(5));
t461 = t409 * t416;
t384 = -t410 * t461 - t458;
t459 = t411 * t416;
t460 = t409 * t419;
t385 = t410 * t460 - t459;
t463 = t408 * t409;
t495 = -t499 * t384 - t501 * t385 - t497 * t463;
t386 = -t410 * t459 + t460;
t387 = t410 * t458 + t461;
t462 = t408 * t411;
t494 = -t499 * t386 - t501 * t387 - t497 * t462;
t493 = t500 * t384 + t502 * t385 - t499 * t463;
t492 = t500 * t386 + t502 * t387 - t499 * t462;
t491 = t502 * t384 + t503 * t385 - t501 * t463;
t490 = t502 * t386 + t503 * t387 - t501 * t462;
t489 = t497 * t410 + (t499 * t416 - t501 * t419) * t408;
t488 = t496 * t409 - t498 * t411;
t487 = t498 * t409 + t496 * t411;
t486 = t499 * t410 + (-t500 * t416 + t502 * t419) * t408;
t485 = t501 * t410 + (-t502 * t416 + t503 * t419) * t408;
t484 = Icges(4,5) * t417 + Icges(5,5) * t408 + Icges(4,6) * t420 + Icges(5,6) * t410;
t465 = Icges(5,4) * t408;
t392 = Icges(5,2) * t410 + t465;
t464 = Icges(5,4) * t410;
t393 = Icges(5,1) * t408 + t464;
t467 = Icges(4,4) * t417;
t398 = Icges(4,2) * t420 + t467;
t466 = Icges(4,4) * t420;
t399 = Icges(4,1) * t417 + t466;
t483 = -t392 * t408 + t393 * t410 - t398 * t417 + t399 * t420;
t436 = -Icges(5,2) * t408 + t464;
t356 = Icges(5,6) * t409 + t411 * t436;
t438 = Icges(5,1) * t410 - t465;
t358 = Icges(5,5) * t409 + t411 * t438;
t437 = -Icges(4,2) * t417 + t466;
t370 = Icges(4,6) * t409 + t411 * t437;
t439 = Icges(4,1) * t420 - t467;
t374 = Icges(4,5) * t409 + t411 * t439;
t482 = -t356 * t408 + t358 * t410 - t370 * t417 + t374 * t420;
t355 = -Icges(5,6) * t411 + t409 * t436;
t357 = -Icges(5,5) * t411 + t409 * t438;
t369 = -Icges(4,6) * t411 + t409 * t437;
t373 = -Icges(4,5) * t411 + t409 * t439;
t481 = t355 * t408 - t357 * t410 + t369 * t417 - t373 * t420;
t418 = sin(qJ(1));
t474 = pkin(1) * t418;
t473 = pkin(3) * t417;
t471 = pkin(3) * t420;
t470 = pkin(5) * t419;
t425 = qJ(6) * t408 + t410 * t470;
t457 = rSges(7,1) * t385 + rSges(7,2) * t384 + rSges(7,3) * t463 - pkin(5) * t459 + t409 * t425;
t456 = rSges(7,1) * t387 + rSges(7,2) * t386 + rSges(7,3) * t462 + pkin(5) * t461 + t411 * t425;
t455 = (-qJ(6) - rSges(7,3)) * t410 + (rSges(7,1) * t419 - rSges(7,2) * t416 + t470) * t408;
t421 = cos(qJ(1));
t407 = qJD(1) * t421 * pkin(1);
t454 = qJD(1) * (pkin(2) * t411 + pkin(7) * t409) + t407;
t453 = qJD(3) * t409;
t452 = qJD(3) * t411;
t451 = qJD(5) * t408;
t351 = -qJ(4) * t411 + t409 * t471;
t352 = qJ(4) * t409 + t411 * t471;
t450 = t351 * t453 + t352 * t452 + qJD(2);
t447 = -pkin(2) * t409 + pkin(7) * t411 - t474;
t446 = -t351 + t447;
t445 = pkin(4) * t410 + pkin(8) * t408;
t444 = rSges(4,1) * t420 - rSges(4,2) * t417;
t443 = rSges(5,1) * t410 - rSges(5,2) * t408;
t442 = qJD(3) * (-rSges(5,1) * t408 - rSges(5,2) * t410 - t473);
t441 = (-pkin(4) * t408 + pkin(8) * t410 - t473) * qJD(3);
t380 = t445 * t409;
t381 = t445 * t411;
t440 = t380 * t453 + t381 * t452 + t450;
t427 = qJD(1) * t352 - qJD(4) * t411 + t454;
t426 = qJD(1) * t381 + t427;
t424 = qJD(6) * t408 + t441;
t404 = qJD(4) * t409;
t423 = t404 + (-t380 + t446) * qJD(1);
t403 = -qJD(5) * t410 + qJD(1);
t402 = rSges(2,1) * t421 - rSges(2,2) * t418;
t401 = rSges(2,1) * t418 + rSges(2,2) * t421;
t400 = rSges(4,1) * t417 + rSges(4,2) * t420;
t389 = t409 * t451 - t452;
t388 = t411 * t451 + t453;
t383 = t407 + qJD(1) * (rSges(3,1) * t411 - rSges(3,2) * t409);
t382 = (-rSges(3,1) * t409 - rSges(3,2) * t411 - t474) * qJD(1);
t378 = rSges(4,3) * t409 + t411 * t444;
t377 = -rSges(4,3) * t411 + t409 * t444;
t376 = -rSges(6,3) * t410 + (rSges(6,1) * t419 - rSges(6,2) * t416) * t408;
t360 = rSges(5,3) * t409 + t411 * t443;
t359 = -rSges(5,3) * t411 + t409 * t443;
t346 = rSges(6,1) * t387 + rSges(6,2) * t386 + rSges(6,3) * t462;
t344 = rSges(6,1) * t385 + rSges(6,2) * t384 + rSges(6,3) * t463;
t328 = qJD(1) * t378 - t400 * t453 + t454;
t327 = -t400 * t452 + (-t377 + t447) * qJD(1);
t326 = qJD(2) + (t377 * t409 + t378 * t411) * qJD(3);
t325 = qJD(1) * t360 + t409 * t442 + t427;
t324 = t404 + t411 * t442 + (-t359 + t446) * qJD(1);
t323 = (t359 * t409 + t360 * t411) * qJD(3) + t450;
t322 = t346 * t403 - t376 * t388 + t409 * t441 + t426;
t321 = -t344 * t403 + t376 * t389 + t411 * t441 + t423;
t320 = t344 * t388 - t346 * t389 + t440;
t319 = -t388 * t455 + t403 * t456 + t409 * t424 + t426;
t318 = t389 * t455 - t403 * t457 + t411 * t424 + t423;
t317 = -qJD(6) * t410 + t388 * t457 - t389 * t456 + t440;
t1 = m(5) * (t323 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + m(7) * (t317 ^ 2 + t318 ^ 2 + t319 ^ 2) / 0.2e1 + m(6) * (t320 ^ 2 + t321 ^ 2 + t322 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(4) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + ((t386 * t486 + t387 * t485 + t462 * t489) * t403 + (t493 * t386 + t491 * t387 + t462 * t495) * t389 + (t386 * t492 + t387 * t490 + t462 * t494) * t388) * t388 / 0.2e1 + ((t384 * t486 + t385 * t485 + t463 * t489) * t403 + (t493 * t384 + t491 * t385 + t463 * t495) * t389 + (t384 * t492 + t385 * t490 + t463 * t494) * t388) * t389 / 0.2e1 + ((-t494 * t388 - t389 * t495 - t489 * t403) * t410 + ((-t416 * t486 + t419 * t485) * t403 + (-t416 * t493 + t419 * t491) * t389 + (-t416 * t492 + t419 * t490) * t388) * t408) * t403 / 0.2e1 + (((-t355 * t410 - t357 * t408 - t369 * t420 - t373 * t417) * t411 + (t356 * t410 + t358 * t408 + t370 * t420 + t374 * t417) * t409) * qJD(3) + (t410 * t392 + t408 * t393 + t420 * t398 + t417 * t399) * qJD(1)) * qJD(1) / 0.2e1 + ((t487 * t409 ^ 2 + (t481 * t411 + (t482 - t488) * t409) * t411) * qJD(3) + (t409 * t484 + t411 * t483) * qJD(1)) * t453 / 0.2e1 - ((t488 * t411 ^ 2 + (t482 * t409 + (t481 - t487) * t411) * t409) * qJD(3) + (t409 * t483 - t411 * t484) * qJD(1)) * t452 / 0.2e1 + (Icges(2,3) + Icges(3,3) + m(2) * (t401 ^ 2 + t402 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
