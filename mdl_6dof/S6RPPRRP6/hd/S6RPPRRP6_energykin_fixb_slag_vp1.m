% Calculate kinetic energy for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:10
% EndTime: 2019-03-09 02:10:11
% DurationCPUTime: 1.28s
% Computational Cost: add. (547->178), mult. (1277->275), div. (0->0), fcn. (1250->6), ass. (0->96)
t424 = Icges(6,1) + Icges(7,1);
t423 = -Icges(6,4) + Icges(7,5);
t422 = Icges(7,4) + Icges(6,5);
t421 = Icges(6,2) + Icges(7,3);
t420 = Icges(7,6) - Icges(6,6);
t419 = Icges(6,3) + Icges(7,2);
t418 = rSges(7,1) + pkin(5);
t417 = rSges(7,3) + qJ(6);
t368 = sin(qJ(4));
t370 = cos(qJ(5));
t372 = cos(qJ(1));
t395 = t372 * t370;
t367 = sin(qJ(5));
t369 = sin(qJ(1));
t400 = t367 * t369;
t341 = t368 * t400 - t395;
t398 = t369 * t370;
t399 = t367 * t372;
t342 = t368 * t398 + t399;
t371 = cos(qJ(4));
t397 = t369 * t371;
t416 = t421 * t341 + t423 * t342 - t420 * t397;
t343 = t368 * t399 + t398;
t344 = t368 * t395 - t400;
t396 = t371 * t372;
t415 = t421 * t343 + t423 * t344 - t420 * t396;
t414 = t420 * t341 + t422 * t342 - t419 * t397;
t413 = t420 * t343 + t422 * t344 - t419 * t396;
t412 = t423 * t341 + t424 * t342 - t422 * t397;
t411 = t423 * t343 + t424 * t344 - t422 * t396;
t410 = (t421 * t367 + t423 * t370) * t371 + t420 * t368;
t409 = (t420 * t367 + t422 * t370) * t371 + t419 * t368;
t408 = (t423 * t367 + t424 * t370) * t371 + t422 * t368;
t403 = pkin(7) * t369;
t402 = Icges(5,4) * t368;
t401 = Icges(5,4) * t371;
t394 = -rSges(7,2) * t397 + t417 * t341 + t418 * t342;
t393 = -rSges(7,2) * t396 + t417 * t343 + t418 * t344;
t392 = rSges(7,2) * t368 + (t417 * t367 + t418 * t370) * t371;
t366 = qJD(2) * t369;
t391 = qJD(3) * t372 + t366;
t390 = qJD(4) * t369;
t389 = qJD(4) * t372;
t388 = qJD(5) * t371;
t387 = -qJD(2) * t372 + qJD(1) * (pkin(1) * t372 + qJ(2) * t369);
t386 = pkin(4) * t368 - pkin(8) * t371;
t385 = rSges(5,1) * t368 + rSges(5,2) * t371;
t384 = Icges(5,1) * t368 + t401;
t383 = Icges(5,2) * t371 + t402;
t382 = Icges(5,5) * t368 + Icges(5,6) * t371;
t330 = Icges(5,6) * t372 + t369 * t383;
t334 = Icges(5,5) * t372 + t369 * t384;
t381 = t330 * t371 + t334 * t368;
t331 = -Icges(5,6) * t369 + t372 * t383;
t335 = -Icges(5,5) * t369 + t372 * t384;
t380 = -t331 * t371 - t335 * t368;
t354 = -Icges(5,2) * t368 + t401;
t355 = Icges(5,1) * t371 - t402;
t379 = t354 * t371 + t355 * t368;
t378 = qJD(1) * t372 * qJ(3) + qJD(3) * t369 + t387;
t356 = pkin(1) * t369 - qJ(2) * t372;
t377 = -pkin(7) * t372 - qJ(3) * t369 - t356;
t346 = t386 * t369;
t347 = t386 * t372;
t376 = (-t346 * t369 - t347 * t372) * qJD(4);
t360 = pkin(4) * t371 + pkin(8) * t368;
t375 = t360 * t390 + t378 + (t347 - t403) * qJD(1);
t374 = t360 * t389 + (-t346 + t377) * qJD(1) + t391;
t361 = qJD(5) * t368 + qJD(1);
t359 = rSges(2,1) * t372 - rSges(2,2) * t369;
t358 = rSges(5,1) * t371 - rSges(5,2) * t368;
t357 = rSges(2,1) * t369 + rSges(2,2) * t372;
t353 = Icges(5,5) * t371 - Icges(5,6) * t368;
t351 = -t369 * t388 + t389;
t350 = -t372 * t388 - t390;
t339 = -rSges(5,3) * t369 + t372 * t385;
t338 = rSges(6,3) * t368 + (rSges(6,1) * t370 - rSges(6,2) * t367) * t371;
t336 = rSges(5,3) * t372 + t369 * t385;
t327 = -Icges(5,3) * t369 + t372 * t382;
t326 = Icges(5,3) * t372 + t369 * t382;
t323 = qJD(1) * (-rSges(3,2) * t372 + rSges(3,3) * t369) + t387;
t322 = t366 + (rSges(3,2) * t369 + rSges(3,3) * t372 - t356) * qJD(1);
t319 = qJD(1) * (rSges(4,2) * t369 + rSges(4,3) * t372) + t378;
t318 = (t372 * rSges(4,2) - t356 + (-rSges(4,3) - qJ(3)) * t369) * qJD(1) + t391;
t317 = rSges(6,1) * t344 - rSges(6,2) * t343 - rSges(6,3) * t396;
t315 = rSges(6,1) * t342 - rSges(6,2) * t341 - rSges(6,3) * t397;
t301 = (-t336 * t369 - t339 * t372) * qJD(4);
t300 = t358 * t390 + (t339 - t403) * qJD(1) + t378;
t299 = t358 * t389 + (-t336 + t377) * qJD(1) + t391;
t298 = t317 * t361 - t338 * t350 + t375;
t297 = -t315 * t361 + t338 * t351 + t374;
t296 = t315 * t350 - t317 * t351 + t376;
t295 = qJD(6) * t341 - t350 * t392 + t361 * t393 + t375;
t294 = qJD(6) * t343 + t351 * t392 - t361 * t394 + t374;
t293 = qJD(6) * t367 * t371 + t350 * t394 - t351 * t393 + t376;
t1 = m(3) * (t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(4) * (t318 ^ 2 + t319 ^ 2) / 0.2e1 + ((t372 * t353 + t369 * t379) * qJD(1) + (t372 ^ 2 * t326 + (t380 * t369 + (-t327 + t381) * t372) * t369) * qJD(4)) * t389 / 0.2e1 + qJD(1) * ((-t368 * t354 + t371 * t355) * qJD(1) + (-(-t331 * t368 + t335 * t371) * t369 + (-t330 * t368 + t334 * t371) * t372) * qJD(4)) / 0.2e1 + m(6) * (t296 ^ 2 + t297 ^ 2 + t298 ^ 2) / 0.2e1 + m(7) * (t293 ^ 2 + t294 ^ 2 + t295 ^ 2) / 0.2e1 + m(5) * (t299 ^ 2 + t300 ^ 2 + t301 ^ 2) / 0.2e1 - ((-t369 * t353 + t372 * t379) * qJD(1) + (t369 ^ 2 * t327 + (t381 * t372 + (-t326 + t380) * t369) * t372) * qJD(4)) * t390 / 0.2e1 + ((t343 * t410 + t344 * t408 - t396 * t409) * t361 + (t343 * t416 + t412 * t344 - t414 * t396) * t351 + (t415 * t343 + t411 * t344 - t413 * t396) * t350) * t350 / 0.2e1 + ((t341 * t410 + t342 * t408 - t397 * t409) * t361 + (t416 * t341 + t412 * t342 - t414 * t397) * t351 + (t341 * t415 + t342 * t411 - t397 * t413) * t350) * t351 / 0.2e1 + (((t367 * t410 + t370 * t408) * t361 + (t367 * t416 + t412 * t370) * t351 + (t367 * t415 + t370 * t411) * t350) * t371 + (t350 * t413 + t351 * t414 + t361 * t409) * t368) * t361 / 0.2e1 + (m(2) * (t357 ^ 2 + t359 ^ 2) + Icges(3,1) + Icges(4,1) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
