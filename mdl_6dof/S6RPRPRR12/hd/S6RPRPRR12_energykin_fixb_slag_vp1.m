% Calculate kinetic energy for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:10
% EndTime: 2019-03-09 04:18:12
% DurationCPUTime: 2.80s
% Computational Cost: add. (818->254), mult. (1559->395), div. (0->0), fcn. (1493->8), ass. (0->128)
t500 = -Icges(4,4) - Icges(5,6);
t499 = Icges(4,1) + Icges(5,2);
t498 = Icges(4,2) + Icges(5,3);
t428 = cos(qJ(3));
t497 = t500 * t428;
t425 = sin(qJ(3));
t496 = t500 * t425;
t495 = -Icges(5,4) + Icges(4,5);
t494 = Icges(5,5) - Icges(4,6);
t493 = -t498 * t428 + t496;
t492 = t499 * t425 - t497;
t491 = Icges(5,1) + Icges(4,3);
t426 = sin(qJ(1));
t429 = cos(qJ(1));
t490 = t493 * t426 + t494 * t429;
t489 = -t494 * t426 + t493 * t429;
t488 = t492 * t426 + t495 * t429;
t487 = t495 * t426 - t492 * t429;
t486 = t498 * t425 + t497;
t485 = t499 * t428 + t496;
t484 = t495 * t425 - t494 * t428;
t483 = t484 * t426 + t491 * t429;
t482 = t491 * t426 - t484 * t429;
t481 = t494 * t425 + t495 * t428;
t480 = t485 * t425 - t486 * t428;
t479 = t487 * t425 + t489 * t428;
t478 = -t488 * t425 + t490 * t428;
t424 = sin(qJ(5));
t473 = pkin(5) * t424;
t477 = pkin(9) * t425 - t428 * t473;
t472 = pkin(8) * t428;
t427 = cos(qJ(5));
t470 = pkin(5) * t427;
t464 = t425 * t426;
t463 = t425 * t429;
t462 = t426 * t428;
t461 = t428 * t429;
t447 = pkin(3) * t425 - qJ(4) * t428;
t389 = t447 * t429;
t419 = qJD(3) * t429;
t460 = qJD(4) * t425 - t389 * t419;
t406 = pkin(3) * t428 + qJ(4) * t425;
t418 = qJD(3) * t426;
t420 = qJD(2) * t426;
t459 = t406 * t418 + t420;
t397 = qJD(1) * (pkin(1) * t429 + qJ(2) * t426);
t458 = qJD(1) * t429 * pkin(7) + t397;
t456 = qJD(5) * t425;
t393 = t426 * t456 + t419;
t457 = qJD(4) * t428;
t413 = qJD(5) * t428 + qJD(1);
t404 = pkin(1) * t426 - qJ(2) * t429;
t452 = -pkin(7) * t426 - t404;
t388 = t447 * t426;
t451 = qJD(1) * t388 + t429 * t457 + t458;
t450 = t389 + t452;
t449 = rSges(4,1) * t425 + rSges(4,2) * t428;
t448 = rSges(5,2) * t425 + rSges(5,3) * t428;
t394 = pkin(4) * t426 - pkin(8) * t463;
t395 = pkin(4) * t429 + pkin(8) * t464;
t434 = t394 * t419 + (-t388 - t395) * t418 + t460;
t433 = qJD(1) * t395 + (-qJD(2) + (-t406 - t472) * qJD(3)) * t429 + t451;
t432 = t418 * t472 + (-t394 + t450) * qJD(1) - t426 * t457 + t459;
t423 = qJ(5) + qJ(6);
t422 = cos(t423);
t421 = sin(t423);
t409 = rSges(2,1) * t429 - rSges(2,2) * t426;
t408 = rSges(4,1) * t428 - rSges(4,2) * t425;
t407 = -rSges(5,2) * t428 + rSges(5,3) * t425;
t405 = rSges(2,1) * t426 + rSges(2,2) * t429;
t396 = qJD(6) * t428 + t413;
t392 = -t429 * t456 + t418;
t387 = -t424 * t462 + t427 * t429;
t386 = -t424 * t429 - t427 * t462;
t385 = t424 * t461 + t426 * t427;
t384 = -t424 * t426 + t427 * t461;
t381 = pkin(9) * t428 + t425 * t473;
t380 = -t421 * t462 + t422 * t429;
t379 = -t421 * t429 - t422 * t462;
t378 = t421 * t461 + t422 * t426;
t377 = -t421 * t426 + t422 * t461;
t376 = t429 * rSges(5,1) - t426 * t448;
t375 = t426 * rSges(5,1) + t429 * t448;
t374 = t426 * rSges(4,3) - t429 * t449;
t373 = t429 * rSges(4,3) + t426 * t449;
t372 = t428 * rSges(6,3) + (rSges(6,1) * t424 + rSges(6,2) * t427) * t425;
t362 = Icges(6,5) * t428 + (Icges(6,1) * t424 + Icges(6,4) * t427) * t425;
t359 = Icges(6,6) * t428 + (Icges(6,4) * t424 + Icges(6,2) * t427) * t425;
t356 = Icges(6,3) * t428 + (Icges(6,5) * t424 + Icges(6,6) * t427) * t425;
t355 = qJD(6) * t464 + t393;
t354 = t418 + (-qJD(5) - qJD(6)) * t463;
t353 = t428 * rSges(7,3) + (rSges(7,1) * t421 + rSges(7,2) * t422) * t425;
t352 = Icges(7,5) * t428 + (Icges(7,1) * t421 + Icges(7,4) * t422) * t425;
t351 = Icges(7,6) * t428 + (Icges(7,4) * t421 + Icges(7,2) * t422) * t425;
t350 = Icges(7,3) * t428 + (Icges(7,5) * t421 + Icges(7,6) * t422) * t425;
t349 = t397 - qJD(2) * t429 + qJD(1) * (-rSges(3,2) * t429 + rSges(3,3) * t426);
t348 = t420 + (rSges(3,2) * t426 + rSges(3,3) * t429 - t404) * qJD(1);
t347 = t426 * t477 + t470 * t429;
t346 = t470 * t426 - t429 * t477;
t345 = rSges(6,1) * t387 + rSges(6,2) * t386 + rSges(6,3) * t464;
t344 = rSges(6,1) * t385 + rSges(6,2) * t384 - rSges(6,3) * t463;
t343 = Icges(6,1) * t387 + Icges(6,4) * t386 + Icges(6,5) * t464;
t342 = Icges(6,1) * t385 + Icges(6,4) * t384 - Icges(6,5) * t463;
t341 = Icges(6,4) * t387 + Icges(6,2) * t386 + Icges(6,6) * t464;
t340 = Icges(6,4) * t385 + Icges(6,2) * t384 - Icges(6,6) * t463;
t339 = Icges(6,5) * t387 + Icges(6,6) * t386 + Icges(6,3) * t464;
t338 = Icges(6,5) * t385 + Icges(6,6) * t384 - Icges(6,3) * t463;
t337 = (-t373 * t426 + t374 * t429) * qJD(3);
t336 = rSges(7,1) * t380 + rSges(7,2) * t379 + rSges(7,3) * t464;
t335 = rSges(7,1) * t378 + rSges(7,2) * t377 - rSges(7,3) * t463;
t334 = Icges(7,1) * t380 + Icges(7,4) * t379 + Icges(7,5) * t464;
t333 = Icges(7,1) * t378 + Icges(7,4) * t377 - Icges(7,5) * t463;
t332 = Icges(7,4) * t380 + Icges(7,2) * t379 + Icges(7,6) * t464;
t331 = Icges(7,4) * t378 + Icges(7,2) * t377 - Icges(7,6) * t463;
t330 = Icges(7,5) * t380 + Icges(7,6) * t379 + Icges(7,3) * t464;
t329 = Icges(7,5) * t378 + Icges(7,6) * t377 - Icges(7,3) * t463;
t328 = qJD(1) * t373 + (-qJD(3) * t408 - qJD(2)) * t429 + t458;
t327 = t408 * t418 + t420 + (-t374 + t452) * qJD(1);
t326 = (t375 * t429 + (-t376 - t388) * t426) * qJD(3) + t460;
t325 = qJD(1) * t376 + (-qJD(2) + (-t406 - t407) * qJD(3)) * t429 + t451;
t324 = (qJD(3) * t407 - t457) * t426 + (-t375 + t450) * qJD(1) + t459;
t323 = t345 * t413 - t372 * t393 + t433;
t322 = -t344 * t413 + t372 * t392 + t432;
t321 = t344 * t393 - t345 * t392 + t434;
t320 = t336 * t396 + t347 * t413 - t353 * t355 - t381 * t393 + t433;
t319 = -t335 * t396 - t346 * t413 + t353 * t354 + t381 * t392 + t432;
t318 = t335 * t355 - t336 * t354 + t346 * t393 - t347 * t392 + t434;
t1 = m(4) * (t327 ^ 2 + t328 ^ 2 + t337 ^ 2) / 0.2e1 + m(7) * (t318 ^ 2 + t319 ^ 2 + t320 ^ 2) / 0.2e1 + m(6) * (t321 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(5) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(3) * (t348 ^ 2 + t349 ^ 2) / 0.2e1 + t355 * ((t330 * t464 + t379 * t332 + t380 * t334) * t355 + (t329 * t464 + t331 * t379 + t333 * t380) * t354 + (t350 * t464 + t351 * t379 + t352 * t380) * t396) / 0.2e1 + t354 * ((-t330 * t463 + t332 * t377 + t334 * t378) * t355 + (-t329 * t463 + t377 * t331 + t378 * t333) * t354 + (-t350 * t463 + t351 * t377 + t352 * t378) * t396) / 0.2e1 + t396 * ((t329 * t354 + t330 * t355 + t350 * t396) * t428 + ((t332 * t422 + t334 * t421) * t355 + (t331 * t422 + t333 * t421) * t354 + (t351 * t422 + t352 * t421) * t396) * t425) / 0.2e1 + t393 * ((t339 * t464 + t386 * t341 + t387 * t343) * t393 + (t338 * t464 + t340 * t386 + t342 * t387) * t392 + (t356 * t464 + t359 * t386 + t362 * t387) * t413) / 0.2e1 + t392 * ((-t339 * t463 + t341 * t384 + t343 * t385) * t393 + (-t338 * t463 + t384 * t340 + t385 * t342) * t392 + (-t356 * t463 + t359 * t384 + t362 * t385) * t413) / 0.2e1 + t413 * ((t338 * t392 + t339 * t393 + t356 * t413) * t428 + ((t341 * t427 + t343 * t424) * t393 + (t340 * t427 + t342 * t424) * t392 + (t359 * t427 + t362 * t424) * t413) * t425) / 0.2e1 + (((t490 * t425 + t488 * t428) * t429 + (-t489 * t425 + t487 * t428) * t426) * qJD(3) + (t486 * t425 + t485 * t428) * qJD(1)) * qJD(1) / 0.2e1 + ((t482 * t426 ^ 2 + (t478 * t429 + (-t479 + t483) * t426) * t429) * qJD(3) + (t426 * t481 - t429 * t480) * qJD(1)) * t418 / 0.2e1 + ((t483 * t429 ^ 2 + (t479 * t426 + (-t478 + t482) * t429) * t426) * qJD(3) + (t426 * t480 + t429 * t481) * qJD(1)) * t419 / 0.2e1 + (m(2) * (t405 ^ 2 + t409 ^ 2) + Icges(3,1) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
