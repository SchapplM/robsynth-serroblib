% Calculate kinetic energy for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP11_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP11_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:24
% EndTime: 2019-12-31 22:14:26
% DurationCPUTime: 2.08s
% Computational Cost: add. (1797->250), mult. (4512->389), div. (0->0), fcn. (5509->10), ass. (0->119)
t478 = Icges(5,1) + Icges(6,1);
t477 = -Icges(5,4) + Icges(6,5);
t476 = Icges(6,4) + Icges(5,5);
t475 = Icges(5,2) + Icges(6,3);
t474 = Icges(6,2) + Icges(5,3);
t473 = -Icges(5,6) + Icges(6,6);
t472 = rSges(6,1) + pkin(4);
t471 = rSges(6,3) + qJ(5);
t430 = sin(qJ(2));
t431 = sin(qJ(1));
t432 = cos(qJ(2));
t433 = cos(qJ(1));
t454 = cos(pkin(5));
t441 = t433 * t454;
t411 = t430 * t431 - t432 * t441;
t412 = t430 * t441 + t431 * t432;
t427 = sin(pkin(5));
t451 = t427 * t433;
t378 = Icges(3,5) * t412 - Icges(3,6) * t411 - Icges(3,3) * t451;
t442 = t431 * t454;
t413 = t433 * t430 + t432 * t442;
t414 = -t430 * t442 + t433 * t432;
t453 = t427 * t431;
t379 = Icges(3,5) * t414 - Icges(3,6) * t413 + Icges(3,3) * t453;
t470 = (t378 * t433 - t379 * t431) * t427;
t429 = sin(qJ(3));
t456 = cos(qJ(3));
t396 = t412 * t456 - t429 * t451;
t428 = sin(qJ(4));
t455 = cos(qJ(4));
t370 = t396 * t428 - t411 * t455;
t371 = t396 * t455 + t411 * t428;
t444 = t427 * t456;
t395 = t412 * t429 + t433 * t444;
t469 = t475 * t370 + t477 * t371 + t473 * t395;
t398 = t414 * t456 + t429 * t453;
t372 = t398 * t428 - t413 * t455;
t373 = t398 * t455 + t413 * t428;
t397 = t414 * t429 - t431 * t444;
t468 = t475 * t372 + t477 * t373 + t473 * t397;
t467 = t473 * t370 + t476 * t371 + t474 * t395;
t466 = t473 * t372 + t476 * t373 + t474 * t397;
t465 = t477 * t370 + t478 * t371 + t476 * t395;
t464 = t477 * t372 + t478 * t373 + t476 * t397;
t410 = t429 * t454 + t430 * t444;
t452 = t427 * t432;
t393 = t410 * t428 + t452 * t455;
t394 = t410 * t455 - t428 * t452;
t409 = t427 * t429 * t430 - t454 * t456;
t463 = t475 * t393 + t477 * t394 + t473 * t409;
t462 = t473 * t393 + t476 * t394 + t474 * t409;
t461 = t477 * t393 + t478 * t394 + t476 * t409;
t450 = rSges(6,2) * t395 + t471 * t370 + t472 * t371;
t449 = rSges(6,2) * t397 + t471 * t372 + t472 * t373;
t448 = rSges(6,2) * t409 + t471 * t393 + t472 * t394;
t390 = pkin(2) * t412 + pkin(8) * t411;
t391 = pkin(2) * t414 + pkin(8) * t413;
t445 = qJD(2) * t427;
t423 = t431 * t445;
t443 = t433 * t445;
t447 = t390 * t423 + t391 * t443;
t399 = qJD(3) * t413 + t423;
t446 = qJD(1) * (pkin(1) * t431 - pkin(7) * t451);
t424 = qJD(2) * t454 + qJD(1);
t400 = qJD(3) * t411 - t443;
t366 = pkin(3) * t396 + pkin(9) * t395;
t367 = pkin(3) * t398 + pkin(9) * t397;
t439 = t399 * t366 - t367 * t400 + t447;
t416 = -qJD(3) * t452 + t424;
t415 = (pkin(2) * t430 - pkin(8) * t432) * t427;
t417 = qJD(1) * (pkin(1) * t433 + pkin(7) * t453);
t438 = t424 * t391 - t415 * t423 + t417;
t437 = -t390 * t424 - t415 * t443 - t446;
t389 = pkin(3) * t410 + pkin(9) * t409;
t436 = t416 * t367 - t389 * t399 + t438;
t435 = -t366 * t416 + t400 * t389 + t437;
t420 = rSges(2,1) * t433 - rSges(2,2) * t431;
t419 = rSges(2,1) * t431 + rSges(2,2) * t433;
t404 = t454 * rSges(3,3) + (rSges(3,1) * t430 + rSges(3,2) * t432) * t427;
t403 = Icges(3,5) * t454 + (Icges(3,1) * t430 + Icges(3,4) * t432) * t427;
t402 = Icges(3,6) * t454 + (Icges(3,4) * t430 + Icges(3,2) * t432) * t427;
t401 = Icges(3,3) * t454 + (Icges(3,5) * t430 + Icges(3,6) * t432) * t427;
t392 = qJD(4) * t409 + t416;
t386 = rSges(3,1) * t414 - rSges(3,2) * t413 + rSges(3,3) * t453;
t385 = rSges(3,1) * t412 - rSges(3,2) * t411 - rSges(3,3) * t451;
t383 = Icges(3,1) * t414 - Icges(3,4) * t413 + Icges(3,5) * t453;
t382 = Icges(3,1) * t412 - Icges(3,4) * t411 - Icges(3,5) * t451;
t381 = Icges(3,4) * t414 - Icges(3,2) * t413 + Icges(3,6) * t453;
t380 = Icges(3,4) * t412 - Icges(3,2) * t411 - Icges(3,6) * t451;
t377 = rSges(4,1) * t410 - rSges(4,2) * t409 - rSges(4,3) * t452;
t376 = Icges(4,1) * t410 - Icges(4,4) * t409 - Icges(4,5) * t452;
t375 = Icges(4,4) * t410 - Icges(4,2) * t409 - Icges(4,6) * t452;
t374 = Icges(4,5) * t410 - Icges(4,6) * t409 - Icges(4,3) * t452;
t369 = qJD(4) * t395 + t400;
t368 = qJD(4) * t397 + t399;
t362 = rSges(4,1) * t398 - rSges(4,2) * t397 + rSges(4,3) * t413;
t361 = rSges(4,1) * t396 - rSges(4,2) * t395 + rSges(4,3) * t411;
t360 = Icges(4,1) * t398 - Icges(4,4) * t397 + Icges(4,5) * t413;
t359 = Icges(4,1) * t396 - Icges(4,4) * t395 + Icges(4,5) * t411;
t358 = Icges(4,4) * t398 - Icges(4,2) * t397 + Icges(4,6) * t413;
t357 = Icges(4,4) * t396 - Icges(4,2) * t395 + Icges(4,6) * t411;
t356 = Icges(4,5) * t398 - Icges(4,6) * t397 + Icges(4,3) * t413;
t355 = Icges(4,5) * t396 - Icges(4,6) * t395 + Icges(4,3) * t411;
t354 = rSges(5,1) * t394 - rSges(5,2) * t393 + rSges(5,3) * t409;
t345 = t386 * t424 - t404 * t423 + t417;
t344 = -t385 * t424 - t404 * t443 - t446;
t343 = (t385 * t431 + t386 * t433) * t445;
t340 = rSges(5,1) * t373 - rSges(5,2) * t372 + rSges(5,3) * t397;
t338 = rSges(5,1) * t371 - rSges(5,2) * t370 + rSges(5,3) * t395;
t324 = t362 * t416 - t377 * t399 + t438;
t323 = -t361 * t416 + t377 * t400 + t437;
t322 = t361 * t399 - t362 * t400 + t447;
t321 = t340 * t392 - t354 * t368 + t436;
t320 = -t338 * t392 + t354 * t369 + t435;
t319 = t338 * t368 - t340 * t369 + t439;
t318 = qJD(5) * t370 - t368 * t448 + t392 * t449 + t436;
t317 = qJD(5) * t372 + t369 * t448 - t392 * t450 + t435;
t316 = qJD(5) * t393 + t368 * t450 - t369 * t449 + t439;
t1 = m(3) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + ((t401 * t453 - t402 * t413 + t403 * t414) * t424 + (-(-t380 * t413 + t382 * t414) * t433 + (-t413 * t381 + t414 * t383 - t470) * t431) * t445) * t423 / 0.2e1 - ((-t401 * t451 - t402 * t411 + t403 * t412) * t424 + ((-t381 * t411 + t383 * t412) * t431 + (t411 * t380 - t412 * t382 + t470) * t433) * t445) * t443 / 0.2e1 + t424 * ((t454 * t379 + (t381 * t432 + t383 * t430) * t427) * t423 - (t454 * t378 + (t380 * t432 + t382 * t430) * t427) * t443 + (t454 * t401 + (t402 * t432 + t403 * t430) * t427) * t424) / 0.2e1 + m(4) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + t399 * ((t413 * t356 - t397 * t358 + t398 * t360) * t399 + (t355 * t413 - t357 * t397 + t359 * t398) * t400 + (t374 * t413 - t375 * t397 + t376 * t398) * t416) / 0.2e1 + t400 * ((t356 * t411 - t358 * t395 + t360 * t396) * t399 + (t411 * t355 - t395 * t357 + t396 * t359) * t400 + (t374 * t411 - t375 * t395 + t376 * t396) * t416) / 0.2e1 + t416 * ((-t356 * t452 - t358 * t409 + t360 * t410) * t399 + (-t355 * t452 - t357 * t409 + t359 * t410) * t400 + (-t374 * t452 - t409 * t375 + t410 * t376) * t416) / 0.2e1 + m(5) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(6) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + ((t463 * t372 + t461 * t373 + t462 * t397) * t392 + (t469 * t372 + t465 * t373 + t467 * t397) * t369 + (t468 * t372 + t464 * t373 + t466 * t397) * t368) * t368 / 0.2e1 + ((t463 * t370 + t461 * t371 + t462 * t395) * t392 + (t469 * t370 + t465 * t371 + t467 * t395) * t369 + (t468 * t370 + t464 * t371 + t466 * t395) * t368) * t369 / 0.2e1 + ((t463 * t393 + t461 * t394 + t462 * t409) * t392 + (t469 * t393 + t465 * t394 + t467 * t409) * t369 + (t468 * t393 + t464 * t394 + t466 * t409) * t368) * t392 / 0.2e1 + (m(2) * (t419 ^ 2 + t420 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
