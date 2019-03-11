% Calculate kinetic energy for
% S6RPPRRP4
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:08
% EndTime: 2019-03-09 02:05:09
% DurationCPUTime: 1.47s
% Computational Cost: add. (1035->181), mult. (2248->284), div. (0->0), fcn. (2706->8), ass. (0->98)
t448 = Icges(6,1) + Icges(7,1);
t447 = Icges(6,4) - Icges(7,5);
t446 = Icges(7,4) + Icges(6,5);
t445 = Icges(6,2) + Icges(7,3);
t444 = Icges(7,6) - Icges(6,6);
t443 = Icges(6,3) + Icges(7,2);
t442 = rSges(7,1) + pkin(5);
t441 = rSges(7,3) + qJ(6);
t424 = sin(pkin(9));
t425 = cos(pkin(9));
t426 = sin(qJ(1));
t427 = cos(qJ(1));
t372 = -t424 * t426 - t425 * t427;
t389 = sin(qJ(5));
t391 = cos(qJ(5));
t373 = t424 * t427 - t425 * t426;
t392 = cos(qJ(4));
t418 = t392 * t373;
t348 = t372 * t391 - t389 * t418;
t349 = -t372 * t389 - t391 * t418;
t390 = sin(qJ(4));
t420 = t373 * t390;
t440 = t445 * t348 - t447 * t349 - t444 * t420;
t419 = t392 * t372;
t350 = -t373 * t391 - t389 * t419;
t351 = t373 * t389 - t391 * t419;
t421 = t372 * t390;
t439 = t445 * t350 - t447 * t351 - t444 * t421;
t438 = t444 * t348 + t446 * t349 - t443 * t420;
t437 = t444 * t350 + t446 * t351 - t443 * t421;
t436 = -t447 * t348 + t448 * t349 - t446 * t420;
t435 = -t447 * t350 + t448 * t351 - t446 * t421;
t434 = t444 * t392 + (-t445 * t389 + t447 * t391) * t390;
t433 = t443 * t392 + (-t444 * t389 - t446 * t391) * t390;
t432 = t446 * t392 + (t447 * t389 - t448 * t391) * t390;
t423 = Icges(5,4) * t390;
t422 = Icges(5,4) * t392;
t417 = -rSges(7,2) * t420 + t441 * t348 + t442 * t349;
t416 = -rSges(7,2) * t421 + t441 * t350 + t442 * t351;
t415 = rSges(7,2) * t392 + (-t441 * t389 - t442 * t391) * t390;
t414 = qJD(4) * t372;
t413 = qJD(4) * t373;
t412 = qJD(4) * (-rSges(5,1) * t390 - rSges(5,2) * t392);
t411 = qJD(4) * (-pkin(4) * t390 + pkin(8) * t392);
t410 = qJD(5) * t390;
t379 = pkin(1) * t426 - qJ(2) * t427;
t409 = -pkin(2) * t426 - t379;
t408 = -pkin(4) * t392 - pkin(8) * t390;
t407 = -qJD(2) * t427 + qJD(1) * (pkin(1) * t427 + qJ(2) * t426);
t406 = -rSges(5,1) * t392 + rSges(5,2) * t390;
t405 = -Icges(5,1) * t392 + t423;
t404 = Icges(5,2) * t390 - t422;
t403 = -Icges(5,5) * t392 + Icges(5,6) * t390;
t340 = -Icges(5,6) * t372 + t373 * t404;
t342 = -Icges(5,5) * t372 + t373 * t405;
t402 = -t340 * t390 + t342 * t392;
t341 = Icges(5,6) * t373 + t372 * t404;
t343 = Icges(5,5) * t373 + t372 * t405;
t401 = t341 * t390 - t343 * t392;
t376 = -Icges(5,2) * t392 - t423;
t377 = -Icges(5,1) * t390 - t422;
t400 = t376 * t390 - t377 * t392;
t399 = pkin(3) * t373 + pkin(7) * t372 + t409;
t398 = qJD(1) * t427 * pkin(2) + t407;
t352 = t408 * t373;
t353 = t408 * t372;
t397 = t352 * t413 + t353 * t414 - qJD(3);
t396 = qJD(1) * (-pkin(3) * t372 + pkin(7) * t373) + t398;
t388 = qJD(2) * t426;
t395 = t388 + (-t352 + t399) * qJD(1) - t372 * t411;
t394 = qJD(1) * t353 - t373 * t411 + t396;
t383 = qJD(5) * t392 + qJD(1);
t381 = rSges(2,1) * t427 - rSges(2,2) * t426;
t380 = rSges(2,1) * t426 + rSges(2,2) * t427;
t375 = -Icges(5,5) * t390 - Icges(5,6) * t392;
t368 = rSges(6,3) * t392 + (-rSges(6,1) * t391 + rSges(6,2) * t389) * t390;
t359 = qJD(1) * (rSges(3,1) * t427 + rSges(3,3) * t426) + t407;
t358 = t388 + (-rSges(3,1) * t426 + rSges(3,3) * t427 - t379) * qJD(1);
t355 = -t373 * t410 - t414;
t354 = -t372 * t410 + t413;
t345 = rSges(5,3) * t373 + t372 * t406;
t344 = -rSges(5,3) * t372 + t373 * t406;
t339 = Icges(5,3) * t373 + t372 * t403;
t338 = -Icges(5,3) * t372 + t373 * t403;
t337 = qJD(1) * (-rSges(4,1) * t372 - rSges(4,2) * t373) + t398;
t336 = t388 + (t373 * rSges(4,1) - t372 * rSges(4,2) + t409) * qJD(1);
t333 = rSges(6,1) * t351 - rSges(6,2) * t350 - rSges(6,3) * t421;
t331 = rSges(6,1) * t349 - rSges(6,2) * t348 - rSges(6,3) * t420;
t317 = qJD(1) * t345 - t373 * t412 + t396;
t316 = -t372 * t412 + t388 + (-t344 + t399) * qJD(1);
t315 = -qJD(3) + (t344 * t373 + t345 * t372) * qJD(4);
t314 = t383 * t333 - t354 * t368 + t394;
t313 = -t383 * t331 + t355 * t368 + t395;
t312 = t331 * t354 - t333 * t355 + t397;
t311 = qJD(6) * t348 - t354 * t415 + t383 * t416 + t394;
t310 = qJD(6) * t350 + t355 * t415 - t383 * t417 + t395;
t309 = -qJD(6) * t389 * t390 + t354 * t417 - t355 * t416 + t397;
t1 = m(6) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(5) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(3) * (t358 ^ 2 + t359 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 - ((-t372 * t375 + t400 * t373) * qJD(1) + (t372 ^ 2 * t338 + (t401 * t373 + (-t339 + t402) * t372) * t373) * qJD(4)) * t414 / 0.2e1 + qJD(1) * ((-t392 * t376 - t390 * t377) * qJD(1) + ((-t341 * t392 - t343 * t390) * t373 - (-t340 * t392 - t342 * t390) * t372) * qJD(4)) / 0.2e1 + m(7) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + ((t372 * t400 + t373 * t375) * qJD(1) + (t373 ^ 2 * t339 + (t402 * t372 + (-t338 + t401) * t373) * t372) * qJD(4)) * t413 / 0.2e1 + ((t434 * t350 + t432 * t351 - t433 * t421) * t383 + (t440 * t350 + t436 * t351 - t438 * t421) * t355 + (t439 * t350 + t435 * t351 - t437 * t421) * t354) * t354 / 0.2e1 + ((t434 * t348 + t432 * t349 - t433 * t420) * t383 + (t440 * t348 + t436 * t349 - t438 * t420) * t355 + (t439 * t348 + t435 * t349 - t437 * t420) * t354) * t355 / 0.2e1 + ((t437 * t354 + t438 * t355 + t433 * t383) * t392 + ((-t434 * t389 - t432 * t391) * t383 + (-t440 * t389 - t436 * t391) * t355 + (-t439 * t389 - t435 * t391) * t354) * t390) * t383 / 0.2e1 + (m(2) * (t380 ^ 2 + t381 ^ 2) + Icges(2,3) + Icges(3,2) + Icges(4,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
