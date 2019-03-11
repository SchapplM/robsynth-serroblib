% Calculate kinetic energy for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:10
% EndTime: 2019-03-09 03:11:13
% DurationCPUTime: 2.34s
% Computational Cost: add. (1264->197), mult. (1520->294), div. (0->0), fcn. (1463->8), ass. (0->114)
t512 = Icges(4,4) + Icges(5,6);
t511 = Icges(4,1) + Icges(5,2);
t510 = -Icges(4,2) - Icges(5,3);
t419 = cos(qJ(3));
t509 = t512 * t419;
t416 = sin(qJ(3));
t508 = t512 * t416;
t507 = -Icges(5,4) + Icges(4,5);
t506 = Icges(5,5) - Icges(4,6);
t505 = t510 * t416 + t509;
t504 = -t511 * t419 + t508;
t503 = Icges(5,1) + Icges(4,3);
t502 = Icges(6,1) + Icges(7,1);
t501 = Icges(6,4) - Icges(7,5);
t500 = Icges(7,4) + Icges(6,5);
t499 = Icges(6,2) + Icges(7,3);
t498 = Icges(7,6) - Icges(6,6);
t497 = Icges(6,3) + Icges(7,2);
t414 = qJ(1) + pkin(9);
t412 = sin(t414);
t413 = cos(t414);
t496 = t412 * t505 + t413 * t506;
t495 = -t412 * t506 + t413 * t505;
t494 = t504 * t412 + t413 * t507;
t493 = t412 * t507 - t504 * t413;
t492 = t510 * t419 - t508;
t491 = t511 * t416 + t509;
t490 = t506 * t416 + t419 * t507;
t489 = rSges(7,1) + pkin(5);
t488 = rSges(7,3) + qJ(6);
t415 = sin(qJ(5));
t418 = cos(qJ(5));
t457 = t416 * t418;
t374 = t412 * t415 - t413 * t457;
t458 = t415 * t416;
t375 = t412 * t418 + t413 * t458;
t456 = t419 * t413;
t487 = t499 * t374 - t501 * t375 + t498 * t456;
t376 = t412 * t457 + t413 * t415;
t377 = t412 * t458 - t413 * t418;
t459 = t412 * t419;
t486 = -t499 * t376 - t501 * t377 + t498 * t459;
t485 = t498 * t374 + t500 * t375 + t497 * t456;
t484 = -t498 * t376 + t500 * t377 + t497 * t459;
t483 = -t501 * t374 + t502 * t375 + t500 * t456;
t482 = t501 * t376 + t502 * t377 + t500 * t459;
t481 = t503 * t412 + t490 * t413;
t480 = t490 * t412 - t503 * t413;
t391 = pkin(4) * t412 + pkin(8) * t456;
t392 = -pkin(4) * t413 + pkin(8) * t459;
t450 = qJD(3) * t413;
t451 = qJD(3) * t412;
t479 = t391 * t450 + t392 * t451;
t478 = (t501 * t415 + t499 * t418) * t419 + t498 * t416;
t477 = (-t500 * t415 + t498 * t418) * t419 + t497 * t416;
t476 = (-t502 * t415 - t501 * t418) * t419 + t500 * t416;
t475 = t416 * t507 - t506 * t419;
t474 = t492 * t416 + t491 * t419;
t473 = t496 * t416 + t494 * t419;
t472 = -t495 * t416 + t493 * t419;
t417 = sin(qJ(1));
t465 = pkin(1) * t417;
t455 = rSges(7,2) * t456 + t488 * t374 + t489 * t375;
t454 = rSges(7,2) * t459 - t488 * t376 + t489 * t377;
t453 = rSges(7,2) * t416 + (-t489 * t415 + t488 * t418) * t419;
t420 = cos(qJ(1));
t411 = qJD(1) * t420 * pkin(1);
t452 = qJD(1) * (pkin(2) * t413 + pkin(7) * t412) + t411;
t449 = qJD(4) * t416;
t448 = qJD(5) * t419;
t438 = pkin(3) * t419 + qJ(4) * t416;
t384 = t438 * t412;
t385 = t438 * t413;
t447 = t384 * t451 + t385 * t450 + qJD(2);
t444 = -pkin(2) * t412 + pkin(7) * t413 - t465;
t405 = pkin(3) * t416 - qJ(4) * t419;
t443 = qJD(3) * (rSges(5,2) * t416 + rSges(5,3) * t419 - t405);
t442 = qJD(1) * t385 + t412 * t449 + t452;
t441 = -t384 + t444;
t440 = rSges(4,1) * t419 - rSges(4,2) * t416;
t439 = -rSges(5,2) * t419 + rSges(5,3) * t416;
t437 = qJD(3) * (-pkin(8) * t416 - t405);
t424 = -qJD(4) * t419 + t447;
t423 = qJD(1) * t391 + t412 * t437 + t442;
t404 = t413 * t449;
t422 = t404 + (-t392 + t441) * qJD(1) + t413 * t437;
t410 = qJD(5) * t416 + qJD(1);
t409 = rSges(2,1) * t420 - rSges(2,2) * t417;
t408 = rSges(2,1) * t417 + rSges(2,2) * t420;
t407 = rSges(4,1) * t416 + rSges(4,2) * t419;
t390 = t412 * t448 - t450;
t389 = t413 * t448 + t451;
t387 = rSges(6,3) * t416 + (-rSges(6,1) * t415 - rSges(6,2) * t418) * t419;
t372 = t411 + qJD(1) * (rSges(3,1) * t413 - rSges(3,2) * t412);
t371 = (-rSges(3,1) * t412 - rSges(3,2) * t413 - t465) * qJD(1);
t368 = -rSges(5,1) * t413 + t412 * t439;
t367 = rSges(5,1) * t412 + t413 * t439;
t366 = rSges(4,3) * t412 + t413 * t440;
t365 = -rSges(4,3) * t413 + t412 * t440;
t348 = rSges(6,1) * t377 + rSges(6,2) * t376 + rSges(6,3) * t459;
t346 = rSges(6,1) * t375 - rSges(6,2) * t374 + rSges(6,3) * t456;
t332 = qJD(1) * t366 - t407 * t451 + t452;
t331 = -t407 * t450 + (-t365 + t444) * qJD(1);
t330 = qJD(2) + (t365 * t412 + t366 * t413) * qJD(3);
t329 = qJD(1) * t367 + t412 * t443 + t442;
t328 = t404 + t413 * t443 + (-t368 + t441) * qJD(1);
t327 = (t367 * t413 + t368 * t412) * qJD(3) + t424;
t326 = t346 * t410 - t387 * t389 + t423;
t325 = -t348 * t410 + t387 * t390 + t422;
t324 = -t346 * t390 + t348 * t389 + t424 + t479;
t323 = -qJD(6) * t376 - t389 * t453 + t410 * t455 + t423;
t322 = qJD(6) * t374 + t390 * t453 - t410 * t454 + t422;
t321 = (qJD(6) * t418 - qJD(4)) * t419 - t455 * t390 + t454 * t389 + t447 + t479;
t1 = m(7) * (t321 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(6) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(5) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + m(4) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + ((t478 * t374 + t476 * t375 + t477 * t456) * t410 + (t486 * t374 + t482 * t375 + t484 * t456) * t390 + (t487 * t374 + t483 * t375 + t485 * t456) * t389) * t389 / 0.2e1 + ((-t478 * t376 + t476 * t377 + t477 * t459) * t410 + (-t486 * t376 + t482 * t377 + t484 * t459) * t390 + (-t487 * t376 + t483 * t377 + t485 * t459) * t389) * t390 / 0.2e1 + (((-t476 * t415 + t478 * t418) * t410 + (-t482 * t415 + t486 * t418) * t390 + (-t483 * t415 + t487 * t418) * t389) * t419 + (t485 * t389 + t484 * t390 + t477 * t410) * t416) * t410 / 0.2e1 + (((t494 * t416 - t496 * t419) * t413 + (t493 * t416 + t495 * t419) * t412) * qJD(3) + (t491 * t416 - t492 * t419) * qJD(1)) * qJD(1) / 0.2e1 + ((t481 * t412 ^ 2 + (t473 * t413 + (t472 - t480) * t412) * t413) * qJD(3) + (t475 * t412 + t474 * t413) * qJD(1)) * t451 / 0.2e1 - ((t480 * t413 ^ 2 + (t472 * t412 + (t473 - t481) * t413) * t412) * qJD(3) + (t474 * t412 - t475 * t413) * qJD(1)) * t450 / 0.2e1 + (Icges(2,3) + Icges(3,3) + m(2) * (t408 ^ 2 + t409 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
