% Calculate kinetic energy for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:26
% EndTime: 2019-03-09 04:31:28
% DurationCPUTime: 1.70s
% Computational Cost: add. (1478->197), mult. (1806->296), div. (0->0), fcn. (1819->8), ass. (0->107)
t446 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t445 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t444 = Icges(7,5) - Icges(6,4) - Icges(5,5);
t443 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t442 = -Icges(6,6) + Icges(7,6) + Icges(5,6);
t441 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t440 = rSges(7,1) + pkin(5);
t439 = rSges(7,3) + qJ(6);
t386 = qJ(1) + pkin(9);
t384 = sin(t386);
t385 = cos(t386);
t390 = cos(qJ(4));
t387 = sin(qJ(4));
t391 = cos(qJ(3));
t418 = t387 * t391;
t351 = -t384 * t390 + t385 * t418;
t417 = t390 * t391;
t352 = t384 * t387 + t385 * t417;
t330 = pkin(4) * t352 + qJ(5) * t351;
t349 = t384 * t418 + t385 * t390;
t382 = -qJD(4) * t391 + qJD(1);
t438 = qJD(5) * t349 + t382 * t330;
t388 = sin(qJ(3));
t410 = qJD(4) * t388;
t411 = qJD(3) * t385;
t368 = t384 * t410 - t411;
t370 = (pkin(4) * t390 + qJ(5) * t387) * t388;
t437 = qJD(5) * t351 + t368 * t370;
t350 = t384 * t417 - t385 * t387;
t420 = t384 * t388;
t436 = t442 * t349 + t444 * t350 - t441 * t420;
t419 = t385 * t388;
t435 = t442 * t351 + t444 * t352 - t441 * t419;
t434 = t443 * t349 + t445 * t350 - t442 * t420;
t433 = t443 * t351 + t445 * t352 - t442 * t419;
t432 = t445 * t349 + t446 * t350 - t444 * t420;
t431 = t445 * t351 + t446 * t352 - t444 * t419;
t430 = t441 * t391 + (t442 * t387 + t444 * t390) * t388;
t429 = t442 * t391 + (t443 * t387 + t445 * t390) * t388;
t428 = t444 * t391 + (t445 * t387 + t446 * t390) * t388;
t389 = sin(qJ(1));
t423 = pkin(1) * t389;
t422 = Icges(4,4) * t388;
t421 = Icges(4,4) * t391;
t416 = rSges(7,2) * t349 + t440 * t350 - t439 * t420;
t415 = rSges(7,2) * t351 + t440 * t352 - t439 * t419;
t414 = t439 * t391 + (rSges(7,2) * t387 + t440 * t390) * t388;
t392 = cos(qJ(1));
t383 = qJD(1) * t392 * pkin(1);
t413 = qJD(1) * (pkin(2) * t385 + pkin(7) * t384) + t383;
t412 = qJD(3) * t384;
t406 = pkin(3) * t391 + pkin(8) * t388;
t366 = t406 * t385;
t409 = qJD(1) * t366 + t413;
t365 = t406 * t384;
t408 = t365 * t412 + t366 * t411 + qJD(2);
t407 = -pkin(2) * t384 + pkin(7) * t385 - t423;
t405 = rSges(4,1) * t391 - rSges(4,2) * t388;
t329 = pkin(4) * t350 + qJ(5) * t349;
t367 = t385 * t410 + t412;
t404 = qJD(5) * t388 * t387 + t367 * t329 + t408;
t403 = Icges(4,1) * t391 - t422;
t402 = -Icges(4,2) * t388 + t421;
t401 = Icges(4,5) * t391 - Icges(4,6) * t388;
t336 = -Icges(4,6) * t385 + t384 * t402;
t338 = -Icges(4,5) * t385 + t384 * t403;
t400 = t336 * t388 - t338 * t391;
t337 = Icges(4,6) * t384 + t385 * t402;
t339 = Icges(4,5) * t384 + t385 * t403;
t399 = -t337 * t388 + t339 * t391;
t374 = Icges(4,2) * t391 + t422;
t375 = Icges(4,1) * t388 + t421;
t398 = -t374 * t388 + t375 * t391;
t379 = pkin(3) * t388 - pkin(8) * t391;
t397 = -qJD(3) * t379 - qJD(6) * t388;
t396 = (-t365 + t407) * qJD(1);
t395 = -t379 * t412 + t409;
t394 = -t379 * t411 + t396;
t378 = rSges(2,1) * t392 - rSges(2,2) * t389;
t377 = rSges(2,1) * t389 + rSges(2,2) * t392;
t376 = rSges(4,1) * t388 + rSges(4,2) * t391;
t373 = Icges(4,5) * t388 + Icges(4,6) * t391;
t364 = -rSges(5,3) * t391 + (rSges(5,1) * t390 - rSges(5,2) * t387) * t388;
t363 = -rSges(6,2) * t391 + (rSges(6,1) * t390 + rSges(6,3) * t387) * t388;
t345 = t383 + qJD(1) * (rSges(3,1) * t385 - rSges(3,2) * t384);
t344 = (-rSges(3,1) * t384 - rSges(3,2) * t385 - t423) * qJD(1);
t341 = rSges(4,3) * t384 + t385 * t405;
t340 = -rSges(4,3) * t385 + t384 * t405;
t335 = Icges(4,3) * t384 + t385 * t401;
t334 = -Icges(4,3) * t385 + t384 * t401;
t327 = rSges(5,1) * t352 - rSges(5,2) * t351 + rSges(5,3) * t419;
t326 = rSges(6,1) * t352 + rSges(6,2) * t419 + rSges(6,3) * t351;
t324 = rSges(5,1) * t350 - rSges(5,2) * t349 + rSges(5,3) * t420;
t323 = rSges(6,1) * t350 + rSges(6,2) * t420 + rSges(6,3) * t349;
t302 = qJD(1) * t341 - t376 * t412 + t413;
t301 = -t376 * t411 + (-t340 + t407) * qJD(1);
t300 = qJD(2) + (t340 * t384 + t341 * t385) * qJD(3);
t299 = t327 * t382 - t364 * t367 + t395;
t298 = -t324 * t382 + t364 * t368 + t394;
t297 = t324 * t367 - t327 * t368 + t408;
t296 = t326 * t382 + (-t363 - t370) * t367 + t395 + t438;
t295 = t363 * t368 + (-t323 - t329) * t382 + t394 + t437;
t294 = t323 * t367 + (-t326 - t330) * t368 + t404;
t293 = t397 * t384 + t415 * t382 + (-t370 - t414) * t367 + t409 + t438;
t292 = t397 * t385 + t414 * t368 + (-t329 - t416) * t382 + t396 + t437;
t291 = qJD(6) * t391 + t416 * t367 + (-t330 - t415) * t368 + t404;
t1 = ((t384 * t373 + t385 * t398) * qJD(1) + (t384 ^ 2 * t335 + (t400 * t385 + (-t334 + t399) * t384) * t385) * qJD(3)) * t412 / 0.2e1 + m(6) * (t294 ^ 2 + t295 ^ 2 + t296 ^ 2) / 0.2e1 + m(7) * (t291 ^ 2 + t292 ^ 2 + t293 ^ 2) / 0.2e1 + m(5) * (t297 ^ 2 + t298 ^ 2 + t299 ^ 2) / 0.2e1 + m(4) * (t300 ^ 2 + t301 ^ 2 + t302 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 - ((-t385 * t373 + t384 * t398) * qJD(1) + (t385 ^ 2 * t334 + (t399 * t384 + (-t335 + t400) * t385) * t384) * qJD(3)) * t411 / 0.2e1 + qJD(1) * ((t391 * t374 + t388 * t375) * qJD(1) + ((t337 * t391 + t339 * t388) * t384 - (t336 * t391 + t338 * t388) * t385) * qJD(3)) / 0.2e1 + ((t429 * t351 + t428 * t352 - t430 * t419) * t382 + (t434 * t351 + t432 * t352 - t436 * t419) * t368 + (t433 * t351 + t431 * t352 - t435 * t419) * t367) * t367 / 0.2e1 + ((t429 * t349 + t428 * t350 - t430 * t420) * t382 + (t434 * t349 + t432 * t350 - t436 * t420) * t368 + (t433 * t349 + t431 * t350 - t435 * t420) * t367) * t368 / 0.2e1 + ((t435 * t367 + t436 * t368 + t430 * t382) * t391 + ((t429 * t387 + t428 * t390) * t382 + (t434 * t387 + t432 * t390) * t368 + (t433 * t387 + t431 * t390) * t367) * t388) * t382 / 0.2e1 + (m(2) * (t377 ^ 2 + t378 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
