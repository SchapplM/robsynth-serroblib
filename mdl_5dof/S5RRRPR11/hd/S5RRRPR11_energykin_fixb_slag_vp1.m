% Calculate kinetic energy for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR11_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR11_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:23
% EndTime: 2019-12-31 21:32:25
% DurationCPUTime: 2.05s
% Computational Cost: add. (897->240), mult. (2161->377), div. (0->0), fcn. (2334->8), ass. (0->119)
t408 = Icges(4,1) + Icges(5,1);
t407 = -Icges(4,4) + Icges(5,5);
t406 = Icges(5,4) + Icges(4,5);
t405 = Icges(4,2) + Icges(5,3);
t404 = -Icges(5,6) + Icges(4,6);
t403 = -Icges(4,3) - Icges(5,2);
t359 = sin(qJ(3));
t363 = cos(qJ(3));
t365 = cos(qJ(1));
t361 = sin(qJ(1));
t364 = cos(qJ(2));
t385 = t361 * t364;
t332 = t359 * t385 + t363 * t365;
t333 = -t359 * t365 + t363 * t385;
t360 = sin(qJ(2));
t387 = t360 * t361;
t402 = t405 * t332 + t407 * t333 - t404 * t387;
t384 = t364 * t365;
t334 = t359 * t384 - t361 * t363;
t335 = t359 * t361 + t363 * t384;
t386 = t360 * t365;
t401 = t405 * t334 + t407 * t335 - t404 * t386;
t400 = -t404 * t332 + t406 * t333 - t403 * t387;
t399 = -t404 * t334 + t406 * t335 - t403 * t386;
t398 = t407 * t332 + t408 * t333 + t406 * t387;
t397 = t407 * t334 + t408 * t335 + t406 * t386;
t396 = t404 * t364 + (t405 * t359 + t407 * t363) * t360;
t395 = t403 * t364 + (-t404 * t359 + t406 * t363) * t360;
t394 = -t406 * t364 + (t407 * t359 + t408 * t363) * t360;
t389 = Icges(3,4) * t360;
t388 = Icges(3,4) * t364;
t378 = pkin(2) * t364 + pkin(7) * t360;
t337 = t378 * t361;
t338 = t378 * t365;
t357 = qJD(2) * t361;
t382 = qJD(2) * t365;
t383 = t337 * t357 + t338 * t382;
t381 = qJD(3) * t360;
t339 = t365 * t381 + t357;
t380 = qJD(5) * t360;
t340 = t361 * t381 - t382;
t302 = pkin(3) * t333 + qJ(4) * t332;
t379 = qJD(4) * t360 * t359 + t339 * t302 + t383;
t377 = rSges(3,1) * t364 - rSges(3,2) * t360;
t376 = Icges(3,1) * t364 - t389;
t375 = -Icges(3,2) * t360 + t388;
t374 = Icges(3,5) * t364 - Icges(3,6) * t360;
t315 = -Icges(3,6) * t365 + t361 * t375;
t319 = -Icges(3,5) * t365 + t361 * t376;
t373 = t315 * t360 - t319 * t364;
t316 = Icges(3,6) * t361 + t365 * t375;
t320 = Icges(3,5) * t361 + t365 * t376;
t372 = -t316 * t360 + t320 * t364;
t345 = Icges(3,2) * t364 + t389;
t346 = Icges(3,1) * t360 + t388;
t371 = -t345 * t360 + t346 * t364;
t343 = qJD(1) * (pkin(1) * t365 + pkin(6) * t361);
t350 = pkin(2) * t360 - pkin(7) * t364;
t370 = qJD(1) * t338 - t350 * t357 + t343;
t303 = pkin(3) * t335 + qJ(4) * t334;
t355 = -qJD(3) * t364 + qJD(1);
t369 = qJD(4) * t332 + t355 * t303 + t370;
t351 = pkin(1) * t361 - pkin(6) * t365;
t368 = (-t337 - t351) * qJD(1) - t350 * t382;
t336 = (pkin(3) * t363 + qJ(4) * t359) * t360;
t367 = qJD(4) * t334 + t340 * t336 + t368;
t362 = cos(qJ(5));
t358 = sin(qJ(5));
t349 = rSges(2,1) * t365 - rSges(2,2) * t361;
t348 = rSges(2,1) * t361 + rSges(2,2) * t365;
t347 = rSges(3,1) * t360 + rSges(3,2) * t364;
t344 = Icges(3,5) * t360 + Icges(3,6) * t364;
t342 = qJD(1) + (-qJD(3) + qJD(5)) * t364;
t341 = pkin(4) * t360 * t363 + pkin(8) * t364;
t328 = (t358 * t359 + t362 * t363) * t360;
t327 = (-t358 * t363 + t359 * t362) * t360;
t324 = rSges(3,3) * t361 + t365 * t377;
t323 = -rSges(3,3) * t365 + t361 * t377;
t322 = -rSges(4,3) * t364 + (rSges(4,1) * t363 - rSges(4,2) * t359) * t360;
t321 = -rSges(5,2) * t364 + (rSges(5,1) * t363 + rSges(5,3) * t359) * t360;
t312 = Icges(3,3) * t361 + t365 * t374;
t311 = -Icges(3,3) * t365 + t361 * t374;
t308 = -t361 * t380 + t340;
t307 = -t365 * t380 + t339;
t306 = pkin(4) * t335 - pkin(8) * t386;
t305 = pkin(4) * t333 - pkin(8) * t387;
t301 = t334 * t358 + t335 * t362;
t300 = t334 * t362 - t335 * t358;
t299 = t332 * t358 + t333 * t362;
t298 = t332 * t362 - t333 * t358;
t297 = rSges(4,1) * t335 - rSges(4,2) * t334 + rSges(4,3) * t386;
t296 = rSges(5,1) * t335 + rSges(5,2) * t386 + rSges(5,3) * t334;
t295 = rSges(4,1) * t333 - rSges(4,2) * t332 + rSges(4,3) * t387;
t294 = rSges(5,1) * t333 + rSges(5,2) * t387 + rSges(5,3) * t332;
t280 = rSges(6,1) * t328 + rSges(6,2) * t327 + rSges(6,3) * t364;
t279 = Icges(6,1) * t328 + Icges(6,4) * t327 + Icges(6,5) * t364;
t278 = Icges(6,4) * t328 + Icges(6,2) * t327 + Icges(6,6) * t364;
t277 = Icges(6,5) * t328 + Icges(6,6) * t327 + Icges(6,3) * t364;
t276 = qJD(1) * t324 - t347 * t357 + t343;
t275 = -t347 * t382 + (-t323 - t351) * qJD(1);
t273 = (t323 * t361 + t324 * t365) * qJD(2);
t272 = rSges(6,1) * t301 + rSges(6,2) * t300 - rSges(6,3) * t386;
t271 = rSges(6,1) * t299 + rSges(6,2) * t298 - rSges(6,3) * t387;
t270 = Icges(6,1) * t301 + Icges(6,4) * t300 - Icges(6,5) * t386;
t269 = Icges(6,1) * t299 + Icges(6,4) * t298 - Icges(6,5) * t387;
t268 = Icges(6,4) * t301 + Icges(6,2) * t300 - Icges(6,6) * t386;
t267 = Icges(6,4) * t299 + Icges(6,2) * t298 - Icges(6,6) * t387;
t266 = Icges(6,5) * t301 + Icges(6,6) * t300 - Icges(6,3) * t386;
t265 = Icges(6,5) * t299 + Icges(6,6) * t298 - Icges(6,3) * t387;
t264 = t297 * t355 - t322 * t339 + t370;
t263 = -t295 * t355 + t322 * t340 + t368;
t262 = t295 * t339 - t297 * t340 + t383;
t261 = t296 * t355 + (-t321 - t336) * t339 + t369;
t260 = t321 * t340 + (-t294 - t302) * t355 + t367;
t259 = t294 * t339 + (-t296 - t303) * t340 + t379;
t258 = t272 * t342 - t280 * t307 + t306 * t355 + (-t336 - t341) * t339 + t369;
t257 = -t271 * t342 + t280 * t308 + t340 * t341 + (-t302 - t305) * t355 + t367;
t256 = t271 * t307 - t272 * t308 + t305 * t339 + (-t303 - t306) * t340 + t379;
t1 = m(3) * (t273 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + ((t361 * t344 + t365 * t371) * qJD(1) + (t361 ^ 2 * t312 + (t373 * t365 + (-t311 + t372) * t361) * t365) * qJD(2)) * t357 / 0.2e1 - ((-t365 * t344 + t361 * t371) * qJD(1) + (t365 ^ 2 * t311 + (t372 * t361 + (-t312 + t373) * t365) * t361) * qJD(2)) * t382 / 0.2e1 + qJD(1) * ((t364 * t345 + t360 * t346) * qJD(1) + ((t316 * t364 + t320 * t360) * t361 - (t315 * t364 + t319 * t360) * t365) * qJD(2)) / 0.2e1 + m(4) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + m(5) * (t259 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + m(6) * (t256 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + t307 * ((-t266 * t386 + t300 * t268 + t301 * t270) * t307 + (-t265 * t386 + t267 * t300 + t269 * t301) * t308 + (-t277 * t386 + t278 * t300 + t279 * t301) * t342) / 0.2e1 + t308 * ((-t266 * t387 + t268 * t298 + t270 * t299) * t307 + (-t265 * t387 + t298 * t267 + t299 * t269) * t308 + (-t277 * t387 + t278 * t298 + t279 * t299) * t342) / 0.2e1 + t342 * ((t266 * t364 + t268 * t327 + t270 * t328) * t307 + (t265 * t364 + t267 * t327 + t269 * t328) * t308 + (t364 * t277 + t327 * t278 + t328 * t279) * t342) / 0.2e1 + ((t396 * t334 + t394 * t335 + t395 * t386) * t355 + (t402 * t334 + t398 * t335 + t400 * t386) * t340 + (t401 * t334 + t397 * t335 + t399 * t386) * t339) * t339 / 0.2e1 + ((t396 * t332 + t394 * t333 + t395 * t387) * t355 + (t402 * t332 + t398 * t333 + t400 * t387) * t340 + (t401 * t332 + t397 * t333 + t399 * t387) * t339) * t340 / 0.2e1 + ((-t399 * t339 - t400 * t340 - t395 * t355) * t364 + ((t396 * t359 + t394 * t363) * t355 + (t402 * t359 + t398 * t363) * t340 + (t401 * t359 + t397 * t363) * t339) * t360) * t355 / 0.2e1 + (m(2) * (t348 ^ 2 + t349 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
