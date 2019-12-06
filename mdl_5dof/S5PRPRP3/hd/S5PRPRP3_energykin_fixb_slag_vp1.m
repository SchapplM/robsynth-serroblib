% Calculate kinetic energy for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:37
% EndTime: 2019-12-05 15:32:38
% DurationCPUTime: 1.26s
% Computational Cost: add. (918->134), mult. (1270->211), div. (0->0), fcn. (1263->8), ass. (0->79)
t423 = Icges(5,1) + Icges(6,1);
t422 = Icges(5,4) + Icges(6,4);
t421 = Icges(5,5) + Icges(6,5);
t420 = Icges(5,2) + Icges(6,2);
t419 = -Icges(6,6) - Icges(5,6);
t418 = Icges(3,3) + Icges(4,3);
t417 = -Icges(6,3) - Icges(5,3);
t344 = qJ(2) + pkin(8);
t342 = sin(t344);
t343 = cos(t344);
t350 = sin(qJ(2));
t352 = cos(qJ(2));
t416 = Icges(3,5) * t352 + Icges(4,5) * t343 - Icges(3,6) * t350 - Icges(4,6) * t342;
t346 = cos(pkin(7));
t398 = t346 ^ 2;
t345 = sin(pkin(7));
t399 = t345 ^ 2;
t401 = t398 + t399;
t400 = qJD(2) * t401;
t415 = t346 * t345;
t351 = cos(qJ(4));
t384 = t346 * t351;
t349 = sin(qJ(4));
t387 = t345 * t349;
t329 = -t343 * t387 - t384;
t385 = t346 * t349;
t386 = t345 * t351;
t330 = t343 * t386 - t385;
t389 = t342 * t345;
t414 = -t419 * t329 + t421 * t330 - t417 * t389;
t331 = -t343 * t385 + t386;
t332 = t343 * t384 + t387;
t388 = t342 * t346;
t413 = -t419 * t331 + t421 * t332 - t417 * t388;
t412 = t420 * t329 + t422 * t330 - t419 * t389;
t411 = t420 * t331 + t422 * t332 - t419 * t388;
t410 = t422 * t329 + t423 * t330 + t421 * t389;
t409 = t422 * t331 + t423 * t332 + t421 * t388;
t408 = t417 * t343 + (t419 * t349 + t421 * t351) * t342;
t407 = t419 * t343 + (-t420 * t349 + t422 * t351) * t342;
t406 = t421 * t343 + (t422 * t349 - t423 * t351) * t342;
t405 = t416 * t345 - t418 * t346;
t404 = t418 * t345 + t416 * t346;
t397 = qJD(2) ^ 2;
t394 = pkin(2) * t350;
t393 = pkin(2) * t352;
t392 = pkin(4) * t351;
t355 = qJ(5) * t342 + t343 * t392;
t383 = rSges(6,1) * t330 + rSges(6,2) * t329 + rSges(6,3) * t389 - pkin(4) * t385 + t345 * t355;
t382 = -rSges(6,1) * t332 - rSges(6,2) * t331 - rSges(6,3) * t388 - pkin(4) * t387 - t346 * t355;
t381 = (-qJ(5) - rSges(6,3)) * t343 + (rSges(6,1) * t351 - rSges(6,2) * t349 + t392) * t342;
t380 = qJD(2) * t345;
t379 = qJD(2) * t346;
t378 = qJD(3) * t346;
t377 = qJD(4) * t342;
t376 = qJD(4) * t343;
t375 = qJD(1) + (-qJ(3) * t346 + t345 * t393) * t380 + (qJ(3) * t345 + t346 * t393) * t379;
t371 = qJD(2) * (-rSges(4,1) * t342 - rSges(4,2) * t343 - t394);
t370 = (-pkin(3) * t342 + pkin(6) * t343 - t394) * qJD(2);
t369 = t375 + (pkin(3) * t343 + pkin(6) * t342) * t400;
t354 = qJD(5) * t342 + t370;
t341 = qJD(3) * t345;
t338 = rSges(3,1) * t350 + rSges(3,2) * t352;
t334 = t345 * t377 - t379;
t333 = t346 * t377 + t380;
t314 = -rSges(5,3) * t343 + (rSges(5,1) * t351 - rSges(5,2) * t349) * t342;
t304 = t346 * t371 + t341;
t303 = t345 * t371 - t378;
t301 = rSges(5,1) * t332 + rSges(5,2) * t331 + rSges(5,3) * t388;
t299 = rSges(5,1) * t330 + rSges(5,2) * t329 + rSges(5,3) * t389;
t285 = qJD(1) + (rSges(3,1) * t352 - rSges(3,2) * t350) * t400;
t282 = t375 + (rSges(4,1) * t343 - rSges(4,2) * t342) * t400;
t281 = t299 * t376 + t314 * t334 + t346 * t370 + t341;
t280 = -t301 * t376 - t314 * t333 + t345 * t370 - t378;
t279 = t299 * t333 - t301 * t334 + t369;
t278 = t334 * t381 + t346 * t354 + t376 * t383 + t341;
t277 = -t333 * t381 + t345 * t354 + t376 * t382 - t378;
t276 = -qJD(5) * t343 + t333 * t383 + t334 * t382 + t369;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t401 * t397 * t338 ^ 2 + t285 ^ 2) / 0.2e1 + m(4) * (t282 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + m(5) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(6) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + ((-t407 * t331 + t406 * t332 - t408 * t388) * t376 + (t412 * t331 + t410 * t332 + t414 * t388) * t334 + (t411 * t331 + t409 * t332 + t413 * t388) * t333) * t333 / 0.2e1 + ((-t407 * t329 + t406 * t330 - t408 * t389) * t376 + (t412 * t329 + t410 * t330 + t414 * t389) * t334 + (t411 * t329 + t409 * t330 + t413 * t389) * t333) * t334 / 0.2e1 - ((-t413 * t333 - t414 * t334 + t408 * t376) * t343 + ((t407 * t349 + t406 * t351) * t376 + (-t412 * t349 + t410 * t351) * t334 + (-t411 * t349 + t409 * t351) * t333) * t342) * t376 / 0.2e1 + (t404 * t399 - t405 * t415) * t345 * t397 / 0.2e1 - (t405 * t398 - t404 * t415) * t346 * t397 / 0.2e1;
T = t1;
