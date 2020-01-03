% Calculate kinetic energy for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP11_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP11_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:20
% EndTime: 2019-12-31 20:12:22
% DurationCPUTime: 2.21s
% Computational Cost: add. (622->184), mult. (1486->279), div. (0->0), fcn. (1449->6), ass. (0->109)
t467 = Icges(3,4) + Icges(4,6);
t466 = Icges(3,1) + Icges(4,2);
t465 = -Icges(3,2) - Icges(4,3);
t375 = cos(qJ(2));
t464 = t467 * t375;
t372 = sin(qJ(2));
t463 = t467 * t372;
t462 = -Icges(4,4) + Icges(3,5);
t461 = Icges(4,5) - Icges(3,6);
t460 = t465 * t372 + t464;
t459 = -t466 * t375 + t463;
t458 = Icges(4,1) + Icges(3,3);
t457 = Icges(5,1) + Icges(6,1);
t456 = Icges(5,4) - Icges(6,5);
t455 = Icges(6,4) + Icges(5,5);
t454 = Icges(5,2) + Icges(6,3);
t453 = Icges(6,6) - Icges(5,6);
t452 = Icges(5,3) + Icges(6,2);
t373 = sin(qJ(1));
t376 = cos(qJ(1));
t451 = t460 * t373 + t461 * t376;
t450 = -t461 * t373 + t460 * t376;
t449 = t459 * t373 + t462 * t376;
t448 = t462 * t373 - t459 * t376;
t447 = t465 * t375 - t463;
t446 = t466 * t372 + t464;
t445 = t461 * t372 + t462 * t375;
t444 = rSges(6,1) + pkin(4);
t443 = rSges(6,3) + qJ(5);
t374 = cos(qJ(4));
t411 = t374 * t376;
t371 = sin(qJ(4));
t415 = t371 * t373;
t342 = -t372 * t411 + t415;
t413 = t373 * t374;
t414 = t371 * t376;
t343 = t372 * t414 + t413;
t410 = t375 * t376;
t442 = t454 * t342 - t456 * t343 + t453 * t410;
t344 = t372 * t413 + t414;
t345 = t372 * t415 - t411;
t412 = t373 * t375;
t441 = -t454 * t344 - t456 * t345 + t453 * t412;
t440 = t453 * t342 + t455 * t343 + t452 * t410;
t439 = -t453 * t344 + t455 * t345 + t452 * t412;
t438 = -t456 * t342 + t457 * t343 + t455 * t410;
t437 = t456 * t344 + t457 * t345 + t455 * t412;
t436 = (t456 * t371 + t454 * t374) * t375 + t453 * t372;
t435 = (-t455 * t371 + t453 * t374) * t375 + t452 * t372;
t434 = t445 * t373 - t458 * t376;
t433 = t458 * t373 + t445 * t376;
t432 = (-t457 * t371 - t456 * t374) * t375 + t455 * t372;
t352 = pkin(3) * t373 + pkin(7) * t410;
t353 = -pkin(3) * t376 + pkin(7) * t412;
t403 = qJD(2) * t376;
t404 = qJD(2) * t373;
t431 = t352 * t403 + t353 * t404;
t430 = t462 * t372 - t461 * t375;
t429 = t447 * t372 + t446 * t375;
t428 = -t450 * t372 + t448 * t375;
t427 = t451 * t372 + t449 * t375;
t409 = rSges(6,2) * t410 + t443 * t342 + t444 * t343;
t408 = rSges(6,2) * t412 - t443 * t344 + t444 * t345;
t393 = pkin(2) * t375 + qJ(3) * t372;
t346 = t393 * t373;
t348 = t393 * t376;
t407 = t346 * t404 + t348 * t403;
t406 = rSges(6,2) * t372 + (-t444 * t371 + t443 * t374) * t375;
t366 = pkin(1) * t373 - pkin(6) * t376;
t405 = -t346 - t366;
t402 = qJD(3) * t372;
t401 = qJD(4) * t375;
t354 = qJD(1) * (pkin(1) * t376 + pkin(6) * t373);
t400 = qJD(1) * t348 + t373 * t402 + t354;
t361 = pkin(2) * t372 - qJ(3) * t375;
t397 = qJD(2) * (rSges(4,2) * t372 + rSges(4,3) * t375 - t361);
t396 = -qJD(3) * t375 + t407;
t395 = rSges(3,1) * t375 - rSges(3,2) * t372;
t394 = -rSges(4,2) * t375 + rSges(4,3) * t372;
t392 = qJD(2) * (-pkin(7) * t372 - t361);
t379 = qJD(1) * t352 + t373 * t392 + t400;
t369 = t376 * t402;
t378 = t369 + (-t353 + t405) * qJD(1) + t376 * t392;
t370 = qJD(4) * t372 + qJD(1);
t365 = rSges(2,1) * t376 - rSges(2,2) * t373;
t364 = rSges(2,1) * t373 + rSges(2,2) * t376;
t363 = rSges(3,1) * t372 + rSges(3,2) * t375;
t351 = t373 * t401 - t403;
t350 = t376 * t401 + t404;
t338 = -rSges(4,1) * t376 + t394 * t373;
t337 = rSges(4,1) * t373 + t394 * t376;
t336 = rSges(3,3) * t373 + t395 * t376;
t335 = rSges(5,3) * t372 + (-rSges(5,1) * t371 - rSges(5,2) * t374) * t375;
t333 = -rSges(3,3) * t376 + t395 * t373;
t310 = rSges(5,1) * t345 + rSges(5,2) * t344 + rSges(5,3) * t412;
t308 = rSges(5,1) * t343 - rSges(5,2) * t342 + rSges(5,3) * t410;
t294 = qJD(1) * t336 - t363 * t404 + t354;
t293 = -t363 * t403 + (-t333 - t366) * qJD(1);
t292 = (t333 * t373 + t336 * t376) * qJD(2);
t291 = qJD(1) * t337 + t373 * t397 + t400;
t290 = t369 + t376 * t397 + (-t338 + t405) * qJD(1);
t289 = (t337 * t376 + t338 * t373) * qJD(2) + t396;
t288 = t308 * t370 - t335 * t350 + t379;
t287 = -t310 * t370 + t335 * t351 + t378;
t286 = -t308 * t351 + t310 * t350 + t396 + t431;
t285 = -qJD(5) * t344 - t406 * t350 + t409 * t370 + t379;
t284 = qJD(5) * t342 + t406 * t351 - t408 * t370 + t378;
t283 = (qJD(5) * t374 - qJD(3)) * t375 - t409 * t351 + t408 * t350 + t407 + t431;
t1 = m(3) * (t292 ^ 2 + t293 ^ 2 + t294 ^ 2) / 0.2e1 + m(4) * (t289 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + m(5) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + m(6) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + ((t436 * t342 + t432 * t343 + t435 * t410) * t370 + (t441 * t342 + t437 * t343 + t439 * t410) * t351 + (t442 * t342 + t438 * t343 + t440 * t410) * t350) * t350 / 0.2e1 + ((-t436 * t344 + t432 * t345 + t435 * t412) * t370 + (-t441 * t344 + t437 * t345 + t439 * t412) * t351 + (-t442 * t344 + t438 * t345 + t440 * t412) * t350) * t351 / 0.2e1 + (((-t432 * t371 + t436 * t374) * t370 + (-t437 * t371 + t441 * t374) * t351 + (-t438 * t371 + t442 * t374) * t350) * t375 + (t440 * t350 + t439 * t351 + t435 * t370) * t372) * t370 / 0.2e1 + (m(2) * (t364 ^ 2 + t365 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t449 * t372 - t451 * t375) * t376 + (t448 * t372 + t450 * t375) * t373) * qJD(2) + (t446 * t372 - t447 * t375) * qJD(1)) * qJD(1) / 0.2e1 + ((t433 * t373 ^ 2 + (t427 * t376 + (t428 - t434) * t373) * t376) * qJD(2) + (t430 * t373 + t429 * t376) * qJD(1)) * t404 / 0.2e1 - ((t434 * t376 ^ 2 + (t428 * t373 + (t427 - t433) * t376) * t373) * qJD(2) + (t429 * t373 - t430 * t376) * qJD(1)) * t403 / 0.2e1;
T = t1;
