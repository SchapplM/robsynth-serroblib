% Calculate kinetic energy for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:34
% EndTime: 2019-03-09 02:12:35
% DurationCPUTime: 1.40s
% Computational Cost: add. (1011->192), mult. (1316->294), div. (0->0), fcn. (1279->8), ass. (0->104)
t462 = Icges(6,1) + Icges(7,1);
t461 = Icges(6,4) + Icges(7,4);
t460 = Icges(7,5) + Icges(6,5);
t459 = Icges(6,2) + Icges(7,2);
t458 = Icges(7,6) + Icges(6,6);
t457 = Icges(7,3) + Icges(6,3);
t396 = pkin(9) + qJ(4);
t391 = sin(t396);
t403 = cos(qJ(5));
t404 = cos(qJ(1));
t431 = t403 * t404;
t401 = sin(qJ(5));
t402 = sin(qJ(1));
t434 = t401 * t402;
t373 = -t391 * t434 + t431;
t432 = t402 * t403;
t433 = t401 * t404;
t374 = t391 * t432 + t433;
t392 = cos(t396);
t436 = t392 * t402;
t456 = t458 * t373 + t460 * t374 - t457 * t436;
t375 = t391 * t433 + t432;
t376 = -t391 * t431 + t434;
t435 = t392 * t404;
t455 = t458 * t375 + t460 * t376 + t457 * t435;
t454 = t459 * t373 + t461 * t374 - t458 * t436;
t453 = t459 * t375 + t461 * t376 + t458 * t435;
t452 = t461 * t373 + t462 * t374 - t460 * t436;
t451 = t461 * t375 + t462 * t376 + t460 * t435;
t450 = (-t458 * t401 + t460 * t403) * t392 + t457 * t391;
t449 = (-t459 * t401 + t461 * t403) * t392 + t458 * t391;
t448 = (-t461 * t401 + t462 * t403) * t392 + t460 * t391;
t447 = qJD(1) * t404 * qJ(3) + qJD(3) * t402;
t440 = pkin(5) * t403;
t446 = -qJ(6) * t392 + t440 * t391;
t397 = sin(pkin(9));
t441 = pkin(3) * t397;
t438 = Icges(5,4) * t391;
t437 = Icges(5,4) * t392;
t429 = rSges(7,1) * t374 + rSges(7,2) * t373 - rSges(7,3) * t436 + pkin(5) * t433 + t402 * t446;
t428 = rSges(7,1) * t376 + rSges(7,2) * t375 + rSges(7,3) * t435 + pkin(5) * t434 - t404 * t446;
t427 = (rSges(7,1) * t403 - rSges(7,2) * t401 + t440) * t392 + (qJ(6) + rSges(7,3)) * t391;
t395 = qJD(2) * t402;
t426 = qJD(3) * t404 + t395;
t425 = qJD(4) * t402;
t424 = qJD(4) * t404;
t423 = qJD(5) * t392;
t422 = qJD(6) * t392;
t384 = qJD(1) * (pkin(1) * t404 + qJ(2) * t402);
t421 = -qJD(2) * t404 + t384;
t383 = pkin(4) * t392 + pkin(8) * t391;
t420 = -qJD(4) * t383 - qJD(2);
t419 = qJD(1) * (pkin(7) * t404 + t402 * t441) + t384 + t447;
t418 = pkin(4) * t391 - pkin(8) * t392;
t385 = pkin(1) * t402 - qJ(2) * t404;
t417 = t404 * t441 - t385 + (-pkin(7) - qJ(3)) * t402;
t369 = t418 * t402;
t370 = t418 * t404;
t416 = -t369 * t425 - t370 * t424;
t398 = cos(pkin(9));
t415 = rSges(4,1) * t397 + rSges(4,2) * t398;
t414 = rSges(5,1) * t391 + rSges(5,2) * t392;
t413 = qJD(1) * t369 + t419;
t412 = Icges(5,1) * t391 + t437;
t411 = Icges(5,2) * t392 + t438;
t410 = Icges(5,5) * t391 + Icges(5,6) * t392;
t360 = Icges(5,6) * t404 + t411 * t402;
t362 = Icges(5,5) * t404 + t412 * t402;
t409 = -t360 * t392 - t362 * t391;
t361 = Icges(5,6) * t402 - t411 * t404;
t363 = Icges(5,5) * t402 - t412 * t404;
t408 = t361 * t392 + t363 * t391;
t380 = -Icges(5,2) * t391 + t437;
t381 = Icges(5,1) * t392 - t438;
t407 = t380 * t392 + t381 * t391;
t406 = t383 * t425 + (t370 + t417) * qJD(1) + t426;
t388 = qJD(5) * t391 + qJD(1);
t387 = rSges(2,1) * t404 - rSges(2,2) * t402;
t386 = rSges(2,1) * t402 + rSges(2,2) * t404;
t382 = rSges(5,1) * t392 - rSges(5,2) * t391;
t379 = Icges(5,5) * t392 - Icges(5,6) * t391;
t378 = -t402 * t423 + t424;
t377 = t404 * t423 + t425;
t365 = rSges(5,3) * t402 - t414 * t404;
t364 = rSges(5,3) * t404 + t414 * t402;
t359 = Icges(5,3) * t402 - t410 * t404;
t358 = Icges(5,3) * t404 + t410 * t402;
t357 = rSges(6,3) * t391 + (rSges(6,1) * t403 - rSges(6,2) * t401) * t392;
t349 = qJD(1) * (-rSges(3,2) * t404 + rSges(3,3) * t402) + t421;
t348 = t395 + (rSges(3,2) * t402 + rSges(3,3) * t404 - t385) * qJD(1);
t346 = rSges(6,1) * t376 + rSges(6,2) * t375 + rSges(6,3) * t435;
t344 = rSges(6,1) * t374 + rSges(6,2) * t373 - rSges(6,3) * t436;
t330 = qJD(1) * (rSges(4,3) * t404 + t415 * t402) + t421 + t447;
t329 = (-t385 + t415 * t404 + (-rSges(4,3) - qJ(3)) * t402) * qJD(1) + t426;
t326 = (-t364 * t402 + t365 * t404) * qJD(4);
t325 = qJD(1) * t364 + (-qJD(4) * t382 - qJD(2)) * t404 + t419;
t324 = t382 * t425 + (-t365 + t417) * qJD(1) + t426;
t323 = -t344 * t377 + t346 * t378 + t416;
t322 = t344 * t388 - t357 * t378 + t420 * t404 + t413;
t321 = -t346 * t388 + t357 * t377 + t406;
t320 = t429 * t388 - t427 * t378 + (t420 + t422) * t404 + t413;
t319 = t427 * t377 - t428 * t388 - t402 * t422 + t406;
t318 = qJD(6) * t391 - t429 * t377 + t428 * t378 + t416;
t1 = ((t379 * t404 + t407 * t402) * qJD(1) + (t358 * t404 ^ 2 + (t408 * t402 + (t359 - t409) * t404) * t402) * qJD(4)) * t424 / 0.2e1 + m(3) * (t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(4) * (t329 ^ 2 + t330 ^ 2) / 0.2e1 + ((t379 * t402 - t407 * t404) * qJD(1) + (t359 * t402 ^ 2 + (t409 * t404 + (t358 - t408) * t402) * t404) * qJD(4)) * t425 / 0.2e1 + qJD(1) * ((-t380 * t391 + t381 * t392) * qJD(1) + ((-t360 * t391 + t362 * t392) * t404 + (-t361 * t391 + t363 * t392) * t402) * qJD(4)) / 0.2e1 + m(6) * (t321 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(7) * (t318 ^ 2 + t319 ^ 2 + t320 ^ 2) / 0.2e1 + m(5) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + ((t375 * t449 + t376 * t448 + t435 * t450) * t388 + (t454 * t375 + t452 * t376 + t435 * t456) * t378 + (t453 * t375 + t451 * t376 + t455 * t435) * t377) * t377 / 0.2e1 + ((t373 * t449 + t374 * t448 - t436 * t450) * t388 + (t454 * t373 + t452 * t374 - t456 * t436) * t378 + (t373 * t453 + t374 * t451 - t455 * t436) * t377) * t378 / 0.2e1 + (((-t449 * t401 + t448 * t403) * t388 + (-t401 * t454 + t403 * t452) * t378 + (-t401 * t453 + t403 * t451) * t377) * t392 + (t455 * t377 + t378 * t456 + t450 * t388) * t391) * t388 / 0.2e1 + (m(2) * (t386 ^ 2 + t387 ^ 2) + Icges(2,3) + Icges(3,1) + Icges(4,1) * t398 ^ 2 + (-0.2e1 * Icges(4,4) * t398 + Icges(4,2) * t397) * t397) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
