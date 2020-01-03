% Calculate kinetic energy for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:12
% EndTime: 2019-12-31 21:22:15
% DurationCPUTime: 2.34s
% Computational Cost: add. (1279->277), mult. (1863->434), div. (0->0), fcn. (1872->10), ass. (0->134)
t424 = -Icges(5,3) - Icges(4,3);
t374 = qJ(3) + pkin(9);
t368 = sin(t374);
t369 = cos(t374);
t381 = cos(qJ(1));
t378 = sin(qJ(1));
t380 = cos(qJ(2));
t407 = t380 * t378;
t321 = -t368 * t407 - t369 * t381;
t322 = -t368 * t381 + t369 * t407;
t376 = sin(qJ(3));
t379 = cos(qJ(3));
t340 = -t376 * t407 - t379 * t381;
t410 = t376 * t381;
t341 = t379 * t407 - t410;
t377 = sin(qJ(2));
t409 = t377 * t378;
t423 = Icges(4,5) * t341 + Icges(5,5) * t322 + Icges(4,6) * t340 + Icges(5,6) * t321 - t424 * t409;
t406 = t380 * t381;
t323 = -t368 * t406 + t369 * t378;
t324 = t368 * t378 + t369 * t406;
t342 = -t376 * t406 + t378 * t379;
t411 = t376 * t378;
t343 = t379 * t406 + t411;
t408 = t377 * t381;
t422 = Icges(4,5) * t343 + Icges(5,5) * t324 + Icges(4,6) * t342 + Icges(5,6) * t323 - t424 * t408;
t421 = t424 * t380 + (Icges(4,5) * t379 + Icges(5,5) * t369 - Icges(4,6) * t376 - Icges(5,6) * t368) * t377;
t415 = t379 * pkin(3);
t413 = Icges(3,4) * t377;
t412 = Icges(3,4) * t380;
t397 = pkin(2) * t380 + pkin(7) * t377;
t344 = t397 * t378;
t345 = t397 * t381;
t371 = qJD(2) * t378;
t402 = qJD(2) * t381;
t405 = t344 * t371 + t345 * t402;
t404 = pkin(4) * t369;
t401 = qJD(3) * t377;
t346 = t381 * t401 + t371;
t400 = qJD(4) * t377;
t399 = qJD(5) * t377;
t398 = pkin(4) * t368;
t347 = t378 * t401 - t402;
t396 = rSges(3,1) * t380 - rSges(3,2) * t377;
t395 = Icges(3,1) * t380 - t413;
t394 = -Icges(3,2) * t377 + t412;
t393 = Icges(3,5) * t380 - Icges(3,6) * t377;
t329 = -Icges(3,6) * t381 + t378 * t394;
t332 = -Icges(3,5) * t381 + t378 * t395;
t392 = t329 * t377 - t332 * t380;
t330 = Icges(3,6) * t378 + t381 * t394;
t333 = Icges(3,5) * t378 + t381 * t395;
t391 = -t330 * t377 + t333 * t380;
t353 = Icges(3,2) * t380 + t413;
t354 = Icges(3,1) * t377 + t412;
t390 = -t353 * t377 + t354 * t380;
t351 = qJD(1) * (pkin(1) * t381 + pkin(6) * t378);
t358 = pkin(2) * t377 - pkin(7) * t380;
t389 = qJD(1) * t345 - t358 * t371 + t351;
t387 = qJ(4) * t377 + t380 * t415;
t300 = -pkin(3) * t410 + t378 * t387;
t388 = -qJD(4) * t380 + t346 * t300 + t405;
t301 = pkin(3) * t411 + t381 * t387;
t364 = -qJD(3) * t380 + qJD(1);
t386 = t364 * t301 + t378 * t400 + t389;
t359 = pkin(1) * t378 - pkin(6) * t381;
t385 = (-t344 - t359) * qJD(1) - t358 * t402;
t384 = pkin(8) * t377 + t380 * t404;
t310 = -qJ(4) * t380 + t377 * t415;
t383 = t347 * t310 + t381 * t400 + t385;
t370 = qJ(5) + t374;
t366 = cos(t370);
t365 = sin(t370);
t357 = rSges(2,1) * t381 - rSges(2,2) * t378;
t356 = rSges(2,1) * t378 + rSges(2,2) * t381;
t355 = rSges(3,1) * t377 + rSges(3,2) * t380;
t352 = Icges(3,5) * t377 + Icges(3,6) * t380;
t349 = qJD(1) + (-qJD(3) - qJD(5)) * t380;
t336 = rSges(3,3) * t378 + t381 * t396;
t335 = -rSges(3,3) * t381 + t378 * t396;
t334 = -rSges(4,3) * t380 + (rSges(4,1) * t379 - rSges(4,2) * t376) * t377;
t331 = -Icges(4,5) * t380 + (Icges(4,1) * t379 - Icges(4,4) * t376) * t377;
t328 = -Icges(4,6) * t380 + (Icges(4,4) * t379 - Icges(4,2) * t376) * t377;
t327 = Icges(3,3) * t378 + t381 * t393;
t326 = -Icges(3,3) * t381 + t378 * t393;
t320 = t378 * t399 + t347;
t319 = t381 * t399 + t346;
t318 = t365 * t378 + t366 * t406;
t317 = -t365 * t406 + t366 * t378;
t316 = -t365 * t381 + t366 * t407;
t315 = -t365 * t407 - t366 * t381;
t314 = -rSges(5,3) * t380 + (rSges(5,1) * t369 - rSges(5,2) * t368) * t377;
t313 = -Icges(5,5) * t380 + (Icges(5,1) * t369 - Icges(5,4) * t368) * t377;
t312 = -Icges(5,6) * t380 + (Icges(5,4) * t369 - Icges(5,2) * t368) * t377;
t309 = -rSges(6,3) * t380 + (rSges(6,1) * t366 - rSges(6,2) * t365) * t377;
t308 = -Icges(6,5) * t380 + (Icges(6,1) * t366 - Icges(6,4) * t365) * t377;
t307 = -Icges(6,6) * t380 + (Icges(6,4) * t366 - Icges(6,2) * t365) * t377;
t306 = -Icges(6,3) * t380 + (Icges(6,5) * t366 - Icges(6,6) * t365) * t377;
t304 = -pkin(8) * t380 + t377 * t404;
t303 = rSges(4,1) * t343 + rSges(4,2) * t342 + rSges(4,3) * t408;
t302 = rSges(4,1) * t341 + rSges(4,2) * t340 + rSges(4,3) * t409;
t299 = Icges(4,1) * t343 + Icges(4,4) * t342 + Icges(4,5) * t408;
t298 = Icges(4,1) * t341 + Icges(4,4) * t340 + Icges(4,5) * t409;
t297 = Icges(4,4) * t343 + Icges(4,2) * t342 + Icges(4,6) * t408;
t296 = Icges(4,4) * t341 + Icges(4,2) * t340 + Icges(4,6) * t409;
t293 = qJD(1) * t336 - t355 * t371 + t351;
t292 = -t355 * t402 + (-t335 - t359) * qJD(1);
t291 = (t335 * t378 + t336 * t381) * qJD(2);
t289 = rSges(5,1) * t324 + rSges(5,2) * t323 + rSges(5,3) * t408;
t288 = rSges(5,1) * t322 + rSges(5,2) * t321 + rSges(5,3) * t409;
t287 = Icges(5,1) * t324 + Icges(5,4) * t323 + Icges(5,5) * t408;
t286 = Icges(5,1) * t322 + Icges(5,4) * t321 + Icges(5,5) * t409;
t285 = Icges(5,4) * t324 + Icges(5,2) * t323 + Icges(5,6) * t408;
t284 = Icges(5,4) * t322 + Icges(5,2) * t321 + Icges(5,6) * t409;
t280 = rSges(6,1) * t318 + rSges(6,2) * t317 + rSges(6,3) * t408;
t279 = rSges(6,1) * t316 + rSges(6,2) * t315 + rSges(6,3) * t409;
t278 = Icges(6,1) * t318 + Icges(6,4) * t317 + Icges(6,5) * t408;
t277 = Icges(6,1) * t316 + Icges(6,4) * t315 + Icges(6,5) * t409;
t276 = Icges(6,4) * t318 + Icges(6,2) * t317 + Icges(6,6) * t408;
t275 = Icges(6,4) * t316 + Icges(6,2) * t315 + Icges(6,6) * t409;
t274 = Icges(6,5) * t318 + Icges(6,6) * t317 + Icges(6,3) * t408;
t273 = Icges(6,5) * t316 + Icges(6,6) * t315 + Icges(6,3) * t409;
t272 = t378 * t398 + t381 * t384;
t271 = t378 * t384 - t381 * t398;
t270 = t303 * t364 - t334 * t346 + t389;
t269 = -t302 * t364 + t334 * t347 + t385;
t268 = t302 * t346 - t303 * t347 + t405;
t267 = t289 * t364 + (-t310 - t314) * t346 + t386;
t266 = t314 * t347 + (-t288 - t300) * t364 + t383;
t265 = t288 * t346 + (-t289 - t301) * t347 + t388;
t264 = t272 * t364 + t280 * t349 - t309 * t319 + (-t304 - t310) * t346 + t386;
t263 = -t279 * t349 + t304 * t347 + t309 * t320 + (-t271 - t300) * t364 + t383;
t262 = t271 * t346 + t279 * t319 - t280 * t320 + (-t272 - t301) * t347 + t388;
t1 = m(3) * (t291 ^ 2 + t292 ^ 2 + t293 ^ 2) / 0.2e1 + ((t378 * t352 + t381 * t390) * qJD(1) + (t378 ^ 2 * t327 + (t392 * t381 + (-t326 + t391) * t378) * t381) * qJD(2)) * t371 / 0.2e1 - ((-t381 * t352 + t378 * t390) * qJD(1) + (t381 ^ 2 * t326 + (t391 * t378 + (-t327 + t392) * t381) * t378) * qJD(2)) * t402 / 0.2e1 + qJD(1) * ((t380 * t353 + t377 * t354) * qJD(1) + ((t330 * t380 + t333 * t377) * t378 - (t329 * t380 + t332 * t377) * t381) * qJD(2)) / 0.2e1 + m(4) * (t268 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + m(5) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + m(6) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + t319 * ((t274 * t408 + t276 * t317 + t278 * t318) * t319 + (t273 * t408 + t275 * t317 + t277 * t318) * t320 + (t306 * t408 + t307 * t317 + t308 * t318) * t349) / 0.2e1 + t320 * ((t274 * t409 + t276 * t315 + t278 * t316) * t319 + (t273 * t409 + t275 * t315 + t277 * t316) * t320 + (t306 * t409 + t307 * t315 + t308 * t316) * t349) / 0.2e1 + t349 * ((-t273 * t320 - t274 * t319 - t306 * t349) * t380 + ((-t276 * t365 + t278 * t366) * t319 + (-t275 * t365 + t277 * t366) * t320 + (-t307 * t365 + t308 * t366) * t349) * t377) / 0.2e1 + ((t312 * t323 + t313 * t324 + t328 * t342 + t331 * t343 + t421 * t408) * t364 + (t284 * t323 + t286 * t324 + t296 * t342 + t298 * t343 + t423 * t408) * t347 + (t285 * t323 + t287 * t324 + t297 * t342 + t299 * t343 + t422 * t408) * t346) * t346 / 0.2e1 + ((t312 * t321 + t313 * t322 + t328 * t340 + t331 * t341 + t421 * t409) * t364 + (t284 * t321 + t286 * t322 + t296 * t340 + t298 * t341 + t423 * t409) * t347 + (t285 * t321 + t287 * t322 + t297 * t340 + t299 * t341 + t422 * t409) * t346) * t347 / 0.2e1 + ((-t422 * t346 - t423 * t347 - t421 * t364) * t380 + ((-t312 * t368 + t313 * t369 - t328 * t376 + t331 * t379) * t364 + (-t284 * t368 + t286 * t369 - t296 * t376 + t298 * t379) * t347 + (-t285 * t368 + t287 * t369 - t297 * t376 + t299 * t379) * t346) * t377) * t364 / 0.2e1 + (m(2) * (t356 ^ 2 + t357 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
