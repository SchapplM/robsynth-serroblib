% Calculate kinetic energy for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:52:51
% EndTime: 2019-03-09 01:52:52
% DurationCPUTime: 1.45s
% Computational Cost: add. (1109->243), mult. (1281->371), div. (0->0), fcn. (1244->10), ass. (0->121)
t404 = sin(qJ(1));
t405 = cos(qJ(1));
t449 = qJD(1) * t405 * qJ(3) + qJD(3) * t404;
t397 = pkin(9) + qJ(4);
t390 = sin(t397);
t392 = cos(t397);
t400 = cos(pkin(10));
t444 = pkin(5) * t400;
t448 = -pkin(8) * t392 + t390 * t444;
t399 = sin(pkin(9));
t446 = pkin(3) * t399;
t443 = Icges(5,4) * t390;
t442 = Icges(5,4) * t392;
t396 = pkin(10) + qJ(6);
t389 = sin(t396);
t441 = t389 * t404;
t440 = t389 * t405;
t391 = cos(t396);
t439 = t391 * t404;
t438 = t391 * t405;
t437 = t392 * t404;
t436 = t392 * t405;
t398 = sin(pkin(10));
t435 = t398 * t404;
t434 = t398 * t405;
t433 = t400 * t404;
t432 = t400 * t405;
t415 = pkin(4) * t390 - qJ(5) * t392;
t366 = t415 * t405;
t426 = qJD(4) * t405;
t429 = qJD(5) * t390 - t366 * t426;
t395 = qJD(2) * t404;
t428 = qJD(3) * t405 + t395;
t427 = qJD(4) * t404;
t425 = qJD(5) * t392;
t424 = qJD(6) * t392;
t378 = pkin(4) * t392 + qJ(5) * t390;
t423 = t378 * t427 + t428;
t380 = qJD(1) * (pkin(1) * t405 + qJ(2) * t404);
t420 = -qJD(2) * t405 + t380;
t419 = qJD(1) * (pkin(7) * t405 + t404 * t446) + t380 + t449;
t382 = pkin(1) * t404 - qJ(2) * t405;
t418 = t405 * t446 - t382 + (-pkin(7) - qJ(3)) * t404;
t401 = cos(pkin(9));
t417 = rSges(4,1) * t399 + rSges(4,2) * t401;
t416 = rSges(5,1) * t390 + rSges(5,2) * t392;
t414 = Icges(5,1) * t390 + t442;
t413 = Icges(5,2) * t392 + t443;
t412 = Icges(5,5) * t390 + Icges(5,6) * t392;
t352 = Icges(5,6) * t405 + t404 * t413;
t354 = Icges(5,5) * t405 + t404 * t414;
t411 = -t352 * t392 - t354 * t390;
t353 = Icges(5,6) * t404 - t405 * t413;
t355 = Icges(5,5) * t404 - t405 * t414;
t410 = t353 * t392 + t355 * t390;
t376 = -Icges(5,2) * t390 + t442;
t377 = Icges(5,1) * t392 - t443;
t409 = t376 * t392 + t377 * t390;
t408 = t366 + t418;
t365 = t415 * t404;
t407 = qJD(1) * t365 + t405 * t425 + t419;
t385 = qJD(6) * t390 + qJD(1);
t384 = rSges(2,1) * t405 - rSges(2,2) * t404;
t383 = rSges(2,1) * t404 + rSges(2,2) * t405;
t379 = rSges(5,1) * t392 - rSges(5,2) * t390;
t375 = Icges(5,5) * t392 - Icges(5,6) * t390;
t374 = -t404 * t424 + t426;
t373 = t405 * t424 + t427;
t372 = -t390 * t432 + t435;
t371 = t390 * t434 + t433;
t370 = t390 * t433 + t434;
t369 = -t390 * t435 + t432;
t362 = -t390 * t438 + t441;
t361 = t390 * t440 + t439;
t360 = t390 * t439 + t440;
t359 = -t390 * t441 + t438;
t358 = rSges(5,3) * t404 - t405 * t416;
t357 = rSges(5,3) * t405 + t404 * t416;
t351 = Icges(5,3) * t404 - t405 * t412;
t350 = Icges(5,3) * t405 + t404 * t412;
t349 = rSges(6,3) * t390 + (rSges(6,1) * t400 - rSges(6,2) * t398) * t392;
t348 = Icges(6,5) * t390 + (Icges(6,1) * t400 - Icges(6,4) * t398) * t392;
t347 = Icges(6,6) * t390 + (Icges(6,4) * t400 - Icges(6,2) * t398) * t392;
t346 = Icges(6,3) * t390 + (Icges(6,5) * t400 - Icges(6,6) * t398) * t392;
t345 = qJD(1) * (-rSges(3,2) * t405 + rSges(3,3) * t404) + t420;
t344 = t395 + (rSges(3,2) * t404 + rSges(3,3) * t405 - t382) * qJD(1);
t343 = rSges(7,3) * t390 + (rSges(7,1) * t391 - rSges(7,2) * t389) * t392;
t342 = Icges(7,5) * t390 + (Icges(7,1) * t391 - Icges(7,4) * t389) * t392;
t341 = Icges(7,6) * t390 + (Icges(7,4) * t391 - Icges(7,2) * t389) * t392;
t340 = Icges(7,3) * t390 + (Icges(7,5) * t391 - Icges(7,6) * t389) * t392;
t339 = pkin(8) * t390 + t392 * t444;
t338 = qJD(1) * (rSges(4,3) * t405 + t404 * t417) + t420 + t449;
t337 = (-t382 + t417 * t405 + (-rSges(4,3) - qJ(3)) * t404) * qJD(1) + t428;
t336 = rSges(6,1) * t372 + rSges(6,2) * t371 + rSges(6,3) * t436;
t335 = rSges(6,1) * t370 + rSges(6,2) * t369 - rSges(6,3) * t437;
t334 = Icges(6,1) * t372 + Icges(6,4) * t371 + Icges(6,5) * t436;
t333 = Icges(6,1) * t370 + Icges(6,4) * t369 - Icges(6,5) * t437;
t332 = Icges(6,4) * t372 + Icges(6,2) * t371 + Icges(6,6) * t436;
t331 = Icges(6,4) * t370 + Icges(6,2) * t369 - Icges(6,6) * t437;
t330 = Icges(6,5) * t372 + Icges(6,6) * t371 + Icges(6,3) * t436;
t329 = Icges(6,5) * t370 + Icges(6,6) * t369 - Icges(6,3) * t437;
t328 = pkin(5) * t435 - t448 * t405;
t327 = pkin(5) * t434 + t448 * t404;
t326 = (-t357 * t404 + t358 * t405) * qJD(4);
t325 = rSges(7,1) * t362 + rSges(7,2) * t361 + rSges(7,3) * t436;
t324 = rSges(7,1) * t360 + rSges(7,2) * t359 - rSges(7,3) * t437;
t323 = Icges(7,1) * t362 + Icges(7,4) * t361 + Icges(7,5) * t436;
t322 = Icges(7,1) * t360 + Icges(7,4) * t359 - Icges(7,5) * t437;
t321 = Icges(7,4) * t362 + Icges(7,2) * t361 + Icges(7,6) * t436;
t320 = Icges(7,4) * t360 + Icges(7,2) * t359 - Icges(7,6) * t437;
t319 = Icges(7,5) * t362 + Icges(7,6) * t361 + Icges(7,3) * t436;
t318 = Icges(7,5) * t360 + Icges(7,6) * t359 - Icges(7,3) * t437;
t317 = qJD(1) * t357 + (-qJD(4) * t379 - qJD(2)) * t405 + t419;
t316 = t379 * t427 + (-t358 + t418) * qJD(1) + t428;
t315 = (t336 * t405 + (-t335 - t365) * t404) * qJD(4) + t429;
t314 = qJD(1) * t335 + (-qJD(2) + (-t349 - t378) * qJD(4)) * t405 + t407;
t313 = (qJD(4) * t349 - t425) * t404 + (-t336 + t408) * qJD(1) + t423;
t312 = qJD(1) * t327 + t324 * t385 - t343 * t374 + (-qJD(2) + (-t339 - t378) * qJD(4)) * t405 + t407;
t311 = -t325 * t385 + t343 * t373 + (qJD(4) * t339 - t425) * t404 + (-t328 + t408) * qJD(1) + t423;
t310 = -t324 * t373 + t325 * t374 + (t328 * t405 + (-t327 - t365) * t404) * qJD(4) + t429;
t1 = m(6) * (t313 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + m(7) * (t310 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + m(5) * (t316 ^ 2 + t317 ^ 2 + t326 ^ 2) / 0.2e1 + m(3) * (t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(4) * (t337 ^ 2 + t338 ^ 2) / 0.2e1 + t385 * ((t318 * t374 + t319 * t373 + t340 * t385) * t390 + ((-t320 * t389 + t322 * t391) * t374 + (-t321 * t389 + t323 * t391) * t373 + (-t341 * t389 + t342 * t391) * t385) * t392) / 0.2e1 + t374 * ((-t318 * t437 + t359 * t320 + t360 * t322) * t374 + (-t319 * t437 + t321 * t359 + t323 * t360) * t373 + (-t340 * t437 + t341 * t359 + t342 * t360) * t385) / 0.2e1 + t373 * ((t318 * t436 + t320 * t361 + t322 * t362) * t374 + (t319 * t436 + t361 * t321 + t362 * t323) * t373 + (t340 * t436 + t341 * t361 + t342 * t362) * t385) / 0.2e1 + (((t329 * t405 + t330 * t404) * t390 + ((-t331 * t398 + t333 * t400) * t405 + (-t332 * t398 + t334 * t400) * t404) * t392 + (-t352 * t390 + t354 * t392) * t405 + (-t353 * t390 + t355 * t392) * t404) * qJD(4) + ((-t347 * t398 + t348 * t400 + t377) * t392 + (t346 - t376) * t390) * qJD(1)) * qJD(1) / 0.2e1 + (((t329 * t436 + t331 * t371 + t333 * t372 + t411 * t405) * t405 + (t330 * t436 + t371 * t332 + t372 * t334 + (t350 - t410) * t405 + t404 * t351) * t404) * qJD(4) + (t346 * t436 + t347 * t371 + t348 * t372 + t404 * t375 - t405 * t409) * qJD(1)) * t427 / 0.2e1 + (((-t329 * t437 + t369 * t331 + t370 * t333 + t405 * t350) * t405 + (-t330 * t437 + t332 * t369 + t334 * t370 + (t351 - t411) * t405 + t410 * t404) * t404) * qJD(4) + (-t346 * t437 + t347 * t369 + t348 * t370 + t405 * t375 + t404 * t409) * qJD(1)) * t426 / 0.2e1 + (m(2) * (t383 ^ 2 + t384 ^ 2) + Icges(2,3) + Icges(3,1) + Icges(4,1) * t401 ^ 2 + (-0.2e1 * Icges(4,4) * t401 + Icges(4,2) * t399) * t399) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
