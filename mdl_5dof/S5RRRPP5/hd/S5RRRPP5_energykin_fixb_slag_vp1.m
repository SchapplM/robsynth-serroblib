% Calculate kinetic energy for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:33
% EndTime: 2019-12-31 20:57:35
% DurationCPUTime: 1.66s
% Computational Cost: add. (778->153), mult. (1055->234), div. (0->0), fcn. (896->6), ass. (0->100)
t424 = Icges(4,4) - Icges(6,4) - Icges(5,5);
t423 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t422 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t339 = qJ(2) + qJ(3);
t338 = cos(t339);
t421 = t424 * t338;
t337 = sin(t339);
t420 = t424 * t337;
t419 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t418 = Icges(4,6) - Icges(5,6) + Icges(6,6);
t417 = t422 * t337 - t421;
t416 = t423 * t338 - t420;
t415 = rSges(6,1) + pkin(4);
t414 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t341 = sin(qJ(1));
t343 = cos(qJ(1));
t413 = -t417 * t341 - t418 * t343;
t412 = -t418 * t341 + t417 * t343;
t411 = t416 * t341 - t419 * t343;
t410 = t419 * t341 + t416 * t343;
t409 = -t422 * t338 - t420;
t408 = t423 * t337 + t421;
t407 = t418 * t337 - t419 * t338;
t406 = rSges(6,3) + qJ(5);
t405 = rSges(6,2) * t337 + t415 * t338;
t380 = qJD(2) + qJD(3);
t325 = t380 * t341;
t326 = t380 * t343;
t404 = (t413 * t337 - t411 * t338) * t326 + (t412 * t337 + t410 * t338) * t325 + (t409 * t337 + t408 * t338) * qJD(1);
t403 = (t407 * t341 + t414 * t343) * t326 + (t414 * t341 - t407 * t343) * t325 + (t419 * t337 + t418 * t338) * qJD(1);
t342 = cos(qJ(2));
t397 = pkin(2) * t342;
t340 = sin(qJ(2));
t395 = Icges(3,4) * t340;
t394 = Icges(3,4) * t342;
t272 = -pkin(7) * t343 + t397 * t341;
t273 = pkin(7) * t341 + t397 * t343;
t382 = qJD(2) * t343;
t383 = qJD(2) * t341;
t387 = t272 * t383 + t273 * t382;
t335 = pkin(1) * t341 - pkin(6) * t343;
t386 = -t272 - t335;
t385 = t405 * t341 + t406 * t343;
t384 = -t406 * t341 + t405 * t343;
t381 = qJD(4) * t337;
t379 = pkin(2) * qJD(2) * t340;
t371 = pkin(3) * t338 + qJ(4) * t337;
t307 = t371 * t341;
t378 = -t307 + t386;
t377 = -rSges(6,2) * t338 + t415 * t337;
t376 = t343 * t379;
t375 = rSges(3,1) * t342 - rSges(3,2) * t340;
t374 = rSges(4,1) * t338 - rSges(4,2) * t337;
t373 = rSges(5,1) * t338 + rSges(5,3) * t337;
t370 = Icges(3,1) * t342 - t395;
t366 = -Icges(3,2) * t340 + t394;
t362 = Icges(3,5) * t342 - Icges(3,6) * t340;
t301 = -Icges(3,6) * t343 + t366 * t341;
t303 = -Icges(3,5) * t343 + t370 * t341;
t358 = t301 * t340 - t303 * t342;
t302 = Icges(3,6) * t341 + t366 * t343;
t304 = Icges(3,5) * t341 + t370 * t343;
t357 = -t302 * t340 + t304 * t342;
t328 = Icges(3,2) * t342 + t395;
t329 = Icges(3,1) * t340 + t394;
t356 = -t328 * t340 + t329 * t342;
t355 = -qJD(4) * t338 + t325 * t307 + t387;
t324 = qJD(1) * (pkin(1) * t343 + pkin(6) * t341);
t354 = qJD(1) * t273 - t341 * t379 + t324;
t320 = pkin(3) * t337 - qJ(4) * t338;
t353 = -t326 * t320 + t343 * t381 - t376;
t308 = t371 * t343;
t349 = qJD(1) * t308 + t341 * t381 + t354;
t332 = rSges(2,1) * t343 - rSges(2,2) * t341;
t331 = rSges(2,1) * t341 + rSges(2,2) * t343;
t330 = rSges(3,1) * t340 + rSges(3,2) * t342;
t327 = Icges(3,5) * t340 + Icges(3,6) * t342;
t323 = rSges(4,1) * t337 + rSges(4,2) * t338;
t322 = rSges(5,1) * t337 - rSges(5,3) * t338;
t306 = rSges(3,3) * t341 + t375 * t343;
t305 = -rSges(3,3) * t343 + t375 * t341;
t300 = Icges(3,3) * t341 + t362 * t343;
t299 = -Icges(3,3) * t343 + t362 * t341;
t297 = rSges(4,3) * t341 + t374 * t343;
t296 = rSges(5,2) * t341 + t373 * t343;
t294 = -rSges(4,3) * t343 + t374 * t341;
t293 = -rSges(5,2) * t343 + t373 * t341;
t266 = qJD(1) * t306 - t330 * t383 + t324;
t265 = -t330 * t382 + (-t305 - t335) * qJD(1);
t264 = (t305 * t341 + t306 * t343) * qJD(2);
t263 = qJD(1) * t297 - t323 * t325 + t354;
t262 = -t376 - t323 * t326 + (-t294 + t386) * qJD(1);
t261 = t294 * t325 + t297 * t326 + t387;
t260 = qJD(1) * t296 + (-t320 - t322) * t325 + t349;
t259 = -t322 * t326 + (-t293 + t378) * qJD(1) + t353;
t258 = qJD(5) * t343 + t384 * qJD(1) + (-t320 - t377) * t325 + t349;
t257 = -qJD(5) * t341 - t377 * t326 + (t378 - t385) * qJD(1) + t353;
t256 = t293 * t325 - (-t296 - t308) * t326 + t355;
t255 = t385 * t325 - (-t308 - t384) * t326 + t355;
t1 = m(3) * (t264 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + ((t341 * t327 + t356 * t343) * qJD(1) + (t341 ^ 2 * t300 + (t358 * t343 + (-t299 + t357) * t341) * t343) * qJD(2)) * t383 / 0.2e1 - ((-t343 * t327 + t356 * t341) * qJD(1) + (t343 ^ 2 * t299 + (t357 * t341 + (-t300 + t358) * t343) * t341) * qJD(2)) * t382 / 0.2e1 + m(4) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(5) * (t256 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + m(6) * (t255 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + (m(2) * (t331 ^ 2 + t332 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (t403 * t341 + t404 * t343) * t325 / 0.2e1 - (t404 * t341 - t403 * t343) * t326 / 0.2e1 + (((t302 * t342 + t304 * t340) * t341 - (t301 * t342 + t303 * t340) * t343) * qJD(2) - (t411 * t337 + t413 * t338) * t326 + (t410 * t337 - t412 * t338) * t325 + (t342 * t328 + t340 * t329 + t408 * t337 - t409 * t338) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
