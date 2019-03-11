% Calculate kinetic energy for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:44
% EndTime: 2019-03-09 02:58:47
% DurationCPUTime: 2.11s
% Computational Cost: add. (541->203), mult. (1271->306), div. (0->0), fcn. (1149->6), ass. (0->112)
t487 = -Icges(4,4) - Icges(6,4) + Icges(5,5);
t486 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t485 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t405 = cos(qJ(3));
t484 = t487 * t405;
t402 = sin(qJ(3));
t483 = t487 * t402;
t482 = Icges(5,4) + Icges(4,5) + Icges(6,6);
t481 = Icges(6,5) + Icges(4,6) - Icges(5,6);
t480 = -t485 * t405 + t483;
t479 = -t486 * t402 + t484;
t478 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t403 = sin(qJ(1));
t406 = cos(qJ(1));
t477 = t480 * t403 - t481 * t406;
t476 = -t481 * t403 - t480 * t406;
t475 = t479 * t403 - t482 * t406;
t474 = t482 * t403 + t479 * t406;
t473 = t485 * t402 + t484;
t472 = t486 * t405 + t483;
t471 = t482 * t402 + t481 * t405;
t470 = t471 * t403 + t478 * t406;
t469 = t478 * t403 - t471 * t406;
t468 = -t481 * t402 + t482 * t405;
t467 = t472 * t402 - t473 * t405;
t466 = t474 * t402 - t476 * t405;
t465 = t475 * t402 + t477 * t405;
t461 = pkin(4) * t405;
t453 = t402 * t403;
t452 = t402 * t406;
t451 = t403 * t405;
t401 = sin(qJ(6));
t450 = t406 * t401;
t404 = cos(qJ(6));
t449 = t406 * t404;
t429 = pkin(3) * t402 - qJ(4) * t405;
t368 = t429 * t406;
t443 = qJD(3) * t406;
t448 = qJD(4) * t402 - t368 * t443;
t367 = t429 * t403;
t375 = pkin(4) * t453 - t406 * qJ(5);
t447 = -t367 - t375;
t390 = t405 * pkin(3) + t402 * qJ(4);
t400 = qJD(2) * t403;
t444 = qJD(3) * t403;
t446 = t390 * t444 + t400;
t377 = qJD(1) * (t406 * pkin(1) + t403 * qJ(2));
t445 = qJD(1) * t406 * pkin(7) + t377;
t442 = qJD(4) * t405;
t441 = qJD(6) * t402;
t376 = -pkin(4) * t452 - t403 * qJ(5);
t440 = t376 * t443 + t448;
t437 = -t390 - t461;
t388 = t403 * pkin(1) - t406 * qJ(2);
t436 = -pkin(7) * t403 - t388;
t435 = qJD(1) * t367 + t406 * t442 + t445;
t434 = t368 + t436;
t433 = pkin(5) * t405 - pkin(8) * t402;
t432 = rSges(4,1) * t402 + rSges(4,2) * t405;
t431 = rSges(5,1) * t402 - rSges(5,3) * t405;
t430 = rSges(6,1) * t405 + rSges(6,2) * t402;
t410 = -t376 + t434;
t409 = -qJD(5) * t406 + t444 * t461 + t446;
t408 = qJD(1) * t375 - qJD(5) * t403 + t435;
t397 = qJD(6) * t405 + qJD(1);
t394 = t402 * pkin(5) + t405 * pkin(8);
t393 = t406 * rSges(2,1) - t403 * rSges(2,2);
t392 = t405 * rSges(4,1) - t402 * rSges(4,2);
t391 = t405 * rSges(5,1) + t402 * rSges(5,3);
t389 = t403 * rSges(2,1) + t406 * rSges(2,2);
t387 = t402 * rSges(6,1) - t405 * rSges(6,2);
t374 = t403 * t441 + t443;
t373 = -t406 * t441 + t444;
t370 = t433 * t406;
t369 = t433 * t403;
t366 = -t403 * t401 + t405 * t449;
t365 = -t403 * t404 - t405 * t450;
t364 = -t404 * t451 - t450;
t363 = t401 * t451 - t449;
t360 = -t403 * rSges(6,3) + t430 * t406;
t359 = t403 * rSges(4,3) - t432 * t406;
t358 = t403 * rSges(5,2) - t431 * t406;
t357 = -t406 * rSges(6,3) - t430 * t403;
t356 = t406 * rSges(4,3) + t432 * t403;
t355 = t406 * rSges(5,2) + t431 * t403;
t354 = t405 * rSges(7,3) + (rSges(7,1) * t404 - rSges(7,2) * t401) * t402;
t346 = Icges(7,5) * t405 + (Icges(7,1) * t404 - Icges(7,4) * t401) * t402;
t339 = Icges(7,6) * t405 + (Icges(7,4) * t404 - Icges(7,2) * t401) * t402;
t332 = Icges(7,3) * t405 + (Icges(7,5) * t404 - Icges(7,6) * t401) * t402;
t331 = t377 - qJD(2) * t406 + qJD(1) * (-t406 * rSges(3,2) + t403 * rSges(3,3));
t330 = t400 + (t403 * rSges(3,2) + t406 * rSges(3,3) - t388) * qJD(1);
t329 = t366 * rSges(7,1) + t365 * rSges(7,2) - rSges(7,3) * t452;
t328 = t364 * rSges(7,1) + t363 * rSges(7,2) + rSges(7,3) * t453;
t327 = Icges(7,1) * t366 + Icges(7,4) * t365 - Icges(7,5) * t452;
t326 = Icges(7,1) * t364 + Icges(7,4) * t363 + Icges(7,5) * t453;
t325 = Icges(7,4) * t366 + Icges(7,2) * t365 - Icges(7,6) * t452;
t324 = Icges(7,4) * t364 + Icges(7,2) * t363 + Icges(7,6) * t453;
t323 = Icges(7,5) * t366 + Icges(7,6) * t365 - Icges(7,3) * t452;
t322 = Icges(7,5) * t364 + Icges(7,6) * t363 + Icges(7,3) * t453;
t321 = (-t356 * t403 + t359 * t406) * qJD(3);
t320 = qJD(1) * t356 + (-qJD(3) * t392 - qJD(2)) * t406 + t445;
t319 = t392 * t444 + t400 + (-t359 + t436) * qJD(1);
t318 = (t358 * t406 + (-t355 - t367) * t403) * qJD(3) + t448;
t317 = qJD(1) * t355 + (-qJD(2) + (-t390 - t391) * qJD(3)) * t406 + t435;
t316 = (qJD(3) * t391 - t442) * t403 + (-t358 + t434) * qJD(1) + t446;
t315 = qJD(1) * t357 + (-qJD(2) + (-t387 + t437) * qJD(3)) * t406 + t408;
t314 = (qJD(3) * t387 - t442) * t403 + (-t360 + t410) * qJD(1) + t409;
t313 = (t360 * t406 + (-t357 + t447) * t403) * qJD(3) + t440;
t312 = -qJD(1) * t369 + t397 * t328 - t374 * t354 + (-qJD(2) + (-t394 + t437) * qJD(3)) * t406 + t408;
t311 = -t397 * t329 + t373 * t354 + (qJD(3) * t394 - t442) * t403 + (-t370 + t410) * qJD(1) + t409;
t310 = -t373 * t328 + t374 * t329 + (t370 * t406 + (t369 + t447) * t403) * qJD(3) + t440;
t1 = m(7) * (t310 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + m(4) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(5) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(6) * (t313 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + m(3) * (t330 ^ 2 + t331 ^ 2) / 0.2e1 + t374 * ((t322 * t453 + t363 * t324 + t364 * t326) * t374 + (t323 * t453 + t363 * t325 + t364 * t327) * t373 + (t332 * t453 + t363 * t339 + t364 * t346) * t397) / 0.2e1 + t373 * ((-t322 * t452 + t365 * t324 + t366 * t326) * t374 + (-t323 * t452 + t365 * t325 + t366 * t327) * t373 + (-t332 * t452 + t365 * t339 + t366 * t346) * t397) / 0.2e1 + t397 * ((t322 * t374 + t323 * t373 + t332 * t397) * t405 + ((-t324 * t401 + t326 * t404) * t374 + (-t325 * t401 + t327 * t404) * t373 + (-t339 * t401 + t346 * t404) * t397) * t402) / 0.2e1 + (Icges(3,1) + m(2) * (t389 ^ 2 + t393 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t477 * t402 - t475 * t405) * t406 + (t476 * t402 + t474 * t405) * t403) * qJD(3) + (t473 * t402 + t472 * t405) * qJD(1)) * qJD(1) / 0.2e1 + ((t469 * t403 ^ 2 + (t465 * t406 + (-t466 + t470) * t403) * t406) * qJD(3) + (t403 * t468 - t467 * t406) * qJD(1)) * t444 / 0.2e1 + ((t470 * t406 ^ 2 + (t466 * t403 + (-t465 + t469) * t406) * t403) * qJD(3) + (t403 * t467 + t406 * t468) * qJD(1)) * t443 / 0.2e1;
T  = t1;
