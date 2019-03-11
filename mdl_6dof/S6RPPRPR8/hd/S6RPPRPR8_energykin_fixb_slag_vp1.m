% Calculate kinetic energy for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:17
% EndTime: 2019-03-09 01:55:19
% DurationCPUTime: 1.81s
% Computational Cost: add. (798->195), mult. (1047->298), div. (0->0), fcn. (954->8), ass. (0->105)
t480 = -Icges(5,4) - Icges(6,6);
t479 = Icges(5,1) + Icges(6,2);
t478 = Icges(5,2) + Icges(6,3);
t402 = pkin(9) + qJ(4);
t398 = cos(t402);
t477 = t480 * t398;
t397 = sin(t402);
t476 = t480 * t397;
t475 = -Icges(6,4) + Icges(5,5);
t474 = Icges(6,5) - Icges(5,6);
t473 = -t478 * t398 + t476;
t472 = t479 * t397 - t477;
t471 = Icges(6,1) + Icges(5,3);
t407 = sin(qJ(1));
t409 = cos(qJ(1));
t470 = t473 * t407 + t474 * t409;
t469 = -t474 * t407 + t473 * t409;
t468 = t472 * t407 + t475 * t409;
t467 = t475 * t407 - t472 * t409;
t466 = t478 * t397 + t477;
t465 = t479 * t398 + t476;
t464 = t475 * t397 - t474 * t398;
t463 = t471 * t407 - t464 * t409;
t462 = t464 * t407 + t471 * t409;
t461 = t474 * t397 + t475 * t398;
t460 = qJD(1) * t409 * qJ(3) + qJD(3) * t407;
t459 = t465 * t397 - t466 * t398;
t458 = -t468 * t397 + t470 * t398;
t457 = t467 * t397 + t469 * t398;
t403 = sin(pkin(9));
t453 = pkin(3) * t403;
t447 = t397 * t407;
t446 = t397 * t409;
t406 = sin(qJ(6));
t445 = t406 * t407;
t444 = t406 * t409;
t408 = cos(qJ(6));
t443 = t407 * t408;
t442 = t408 * t409;
t425 = pkin(4) * t397 - qJ(5) * t398;
t369 = t425 * t409;
t437 = qJD(4) * t409;
t440 = qJD(5) * t397 - t369 * t437;
t401 = qJD(2) * t407;
t439 = qJD(3) * t409 + t401;
t438 = qJD(4) * t407;
t436 = qJD(5) * t398;
t435 = qJD(6) * t397;
t385 = pkin(4) * t398 + qJ(5) * t397;
t434 = t385 * t438 + t439;
t389 = qJD(1) * (pkin(1) * t409 + qJ(2) * t407);
t431 = -qJD(2) * t409 + t389;
t430 = qJD(1) * (pkin(7) * t409 + t407 * t453) + t389 + t460;
t391 = pkin(1) * t407 - qJ(2) * t409;
t429 = t409 * t453 - t391 + (-pkin(7) - qJ(3)) * t407;
t404 = cos(pkin(9));
t428 = rSges(4,1) * t403 + rSges(4,2) * t404;
t427 = rSges(5,1) * t397 + rSges(5,2) * t398;
t426 = rSges(6,2) * t397 + rSges(6,3) * t398;
t412 = t369 + t429;
t368 = t425 * t407;
t411 = qJD(1) * t368 + t409 * t436 + t430;
t394 = qJD(6) * t398 + qJD(1);
t393 = rSges(2,1) * t409 - rSges(2,2) * t407;
t392 = rSges(2,1) * t407 + rSges(2,2) * t409;
t388 = pkin(5) * t409 + pkin(8) * t447;
t387 = rSges(5,1) * t398 - rSges(5,2) * t397;
t386 = -rSges(6,2) * t398 + rSges(6,3) * t397;
t384 = pkin(5) * t407 - pkin(8) * t446;
t377 = t407 * t435 + t437;
t376 = -t409 * t435 + t438;
t375 = -t398 * t445 + t442;
t374 = -t398 * t443 - t444;
t373 = t398 * t444 + t443;
t372 = t398 * t442 - t445;
t365 = rSges(6,1) * t409 - t407 * t426;
t364 = rSges(6,1) * t407 + t409 * t426;
t363 = rSges(5,3) * t407 - t409 * t427;
t362 = rSges(5,3) * t409 + t407 * t427;
t348 = rSges(7,3) * t398 + (rSges(7,1) * t406 + rSges(7,2) * t408) * t397;
t347 = Icges(7,5) * t398 + (Icges(7,1) * t406 + Icges(7,4) * t408) * t397;
t346 = Icges(7,6) * t398 + (Icges(7,4) * t406 + Icges(7,2) * t408) * t397;
t345 = Icges(7,3) * t398 + (Icges(7,5) * t406 + Icges(7,6) * t408) * t397;
t344 = qJD(1) * (-rSges(3,2) * t409 + rSges(3,3) * t407) + t431;
t343 = t401 + (rSges(3,2) * t407 + rSges(3,3) * t409 - t391) * qJD(1);
t342 = rSges(7,1) * t375 + rSges(7,2) * t374 + rSges(7,3) * t447;
t341 = rSges(7,1) * t373 + rSges(7,2) * t372 - rSges(7,3) * t446;
t340 = Icges(7,1) * t375 + Icges(7,4) * t374 + Icges(7,5) * t447;
t339 = Icges(7,1) * t373 + Icges(7,4) * t372 - Icges(7,5) * t446;
t338 = Icges(7,4) * t375 + Icges(7,2) * t374 + Icges(7,6) * t447;
t337 = Icges(7,4) * t373 + Icges(7,2) * t372 - Icges(7,6) * t446;
t336 = Icges(7,5) * t375 + Icges(7,6) * t374 + Icges(7,3) * t447;
t335 = Icges(7,5) * t373 + Icges(7,6) * t372 - Icges(7,3) * t446;
t334 = qJD(1) * (rSges(4,3) * t409 + t407 * t428) + t431 + t460;
t333 = (-t391 + t428 * t409 + (-rSges(4,3) - qJ(3)) * t407) * qJD(1) + t439;
t332 = (-t362 * t407 + t363 * t409) * qJD(4);
t331 = qJD(1) * t362 + (-qJD(4) * t387 - qJD(2)) * t409 + t430;
t330 = t387 * t438 + (-t363 + t429) * qJD(1) + t439;
t329 = (t364 * t409 + (-t365 - t368) * t407) * qJD(4) + t440;
t328 = qJD(1) * t365 + (-qJD(2) + (-t385 - t386) * qJD(4)) * t409 + t411;
t327 = (qJD(4) * t386 - t436) * t407 + (-t364 + t412) * qJD(1) + t434;
t326 = t341 * t377 - t342 * t376 + (t384 * t409 + (-t368 - t388) * t407) * qJD(4) + t440;
t325 = qJD(1) * t388 + t342 * t394 - t348 * t377 + (-qJD(2) + (-pkin(8) * t398 - t385) * qJD(4)) * t409 + t411;
t324 = -t341 * t394 + t348 * t376 + (pkin(8) * qJD(4) - qJD(5)) * t407 * t398 + (-t384 + t412) * qJD(1) + t434;
t1 = m(7) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(6) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + m(5) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(4) * (t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(3) * (t343 ^ 2 + t344 ^ 2) / 0.2e1 + t394 * ((t335 * t376 + t336 * t377 + t345 * t394) * t398 + ((t338 * t408 + t340 * t406) * t377 + (t337 * t408 + t339 * t406) * t376 + (t346 * t408 + t347 * t406) * t394) * t397) / 0.2e1 + t377 * ((t336 * t447 + t374 * t338 + t375 * t340) * t377 + (t335 * t447 + t337 * t374 + t339 * t375) * t376 + (t345 * t447 + t346 * t374 + t347 * t375) * t394) / 0.2e1 + t376 * ((-t336 * t446 + t338 * t372 + t340 * t373) * t377 + (-t335 * t446 + t372 * t337 + t373 * t339) * t376 + (-t345 * t446 + t346 * t372 + t347 * t373) * t394) / 0.2e1 + (((t470 * t397 + t468 * t398) * t409 + (-t469 * t397 + t467 * t398) * t407) * qJD(4) + (t466 * t397 + t465 * t398) * qJD(1)) * qJD(1) / 0.2e1 + ((t463 * t407 ^ 2 + (t458 * t409 + (-t457 + t462) * t407) * t409) * qJD(4) + (t461 * t407 - t459 * t409) * qJD(1)) * t438 / 0.2e1 + ((t462 * t409 ^ 2 + (t457 * t407 + (-t458 + t463) * t409) * t407) * qJD(4) + (t459 * t407 + t461 * t409) * qJD(1)) * t437 / 0.2e1 + (m(2) * (t392 ^ 2 + t393 ^ 2) + Icges(3,1) + Icges(4,1) * t404 ^ 2 + (-0.2e1 * Icges(4,4) * t404 + Icges(4,2) * t403) * t403 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
