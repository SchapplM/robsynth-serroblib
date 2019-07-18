% Calculate kinetic energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:07
% EndTime: 2019-07-18 17:17:08
% DurationCPUTime: 1.70s
% Computational Cost: add. (1212->247), mult. (1491->418), div. (0->0), fcn. (1446->10), ass. (0->135)
t343 = sin(qJ(1));
t390 = t343 ^ 2;
t346 = cos(qJ(1));
t389 = t346 ^ 2;
t344 = cos(qJ(4));
t386 = pkin(3) * t344;
t342 = sin(qJ(2));
t385 = Icges(3,4) * t342;
t345 = cos(qJ(2));
t384 = Icges(3,4) * t345;
t340 = qJ(2) + qJ(3);
t336 = sin(t340);
t383 = Icges(4,4) * t336;
t338 = cos(t340);
t382 = Icges(4,4) * t338;
t381 = t336 * t343;
t380 = t336 * t346;
t379 = t338 * t343;
t378 = t338 * t346;
t341 = sin(qJ(4));
t377 = t341 * t343;
t376 = t341 * t346;
t375 = t343 * t344;
t374 = t344 * t346;
t373 = (t389 + t390) * pkin(1) * qJD(2) * t345;
t334 = qJD(2) * t343;
t318 = qJD(3) * t343 + t334;
t372 = qJD(1) * t345;
t371 = qJD(2) * t346;
t370 = qJD(4) * t336;
t369 = qJD(5) * t336;
t301 = t346 * t370 + t318;
t367 = t338 * t386;
t319 = (-qJD(2) - qJD(3)) * t346;
t366 = pkin(2) * t338 + pkin(5) * t336;
t365 = rSges(3,1) * t345 - rSges(3,2) * t342;
t364 = rSges(4,1) * t338 - rSges(4,2) * t336;
t363 = Icges(3,1) * t345 - t385;
t362 = Icges(4,1) * t338 - t383;
t361 = -Icges(3,2) * t342 + t384;
t360 = -Icges(4,2) * t336 + t382;
t359 = Icges(3,5) * t345 - Icges(3,6) * t342;
t358 = Icges(4,5) * t338 - Icges(4,6) * t336;
t297 = -Icges(3,6) * t346 + t343 * t361;
t299 = -Icges(3,5) * t346 + t343 * t363;
t357 = t297 * t342 - t299 * t345;
t298 = Icges(3,6) * t343 + t346 * t361;
t300 = Icges(3,5) * t343 + t346 * t363;
t356 = -t298 * t342 + t300 * t345;
t321 = Icges(3,2) * t345 + t385;
t322 = Icges(3,1) * t342 + t384;
t355 = -t321 * t342 + t322 * t345;
t302 = t343 * t370 + t319;
t354 = (-t334 * t342 + t346 * t372) * pkin(1);
t305 = t366 * t343;
t306 = t366 * t346;
t353 = t318 * t305 - t306 * t319 + t373;
t352 = (Icges(4,5) * t336 + Icges(4,6) * t338) * qJD(1) + (-Icges(4,3) * t346 + t343 * t358) * t319 + (Icges(4,3) * t343 + t346 * t358) * t318;
t351 = (-t342 * t371 - t343 * t372) * pkin(1);
t317 = pkin(2) * t336 - pkin(5) * t338;
t350 = qJD(1) * t306 - t317 * t318 + t354;
t349 = -qJD(1) * t305 + t319 * t317 + t351;
t284 = -Icges(4,6) * t346 + t343 * t360;
t285 = Icges(4,6) * t343 + t346 * t360;
t286 = -Icges(4,5) * t346 + t343 * t362;
t287 = Icges(4,5) * t343 + t346 * t362;
t314 = Icges(4,2) * t338 + t383;
t315 = Icges(4,1) * t336 + t382;
t348 = (-t285 * t336 + t287 * t338) * t318 + (-t284 * t336 + t286 * t338) * t319 + (-t314 * t336 + t315 * t338) * qJD(1);
t339 = qJ(4) + qJ(5);
t337 = cos(t339);
t335 = sin(t339);
t328 = -qJD(4) * t338 + qJD(1);
t325 = rSges(2,1) * t346 - rSges(2,2) * t343;
t324 = rSges(2,1) * t343 + rSges(2,2) * t346;
t323 = rSges(3,1) * t342 + rSges(3,2) * t345;
t320 = Icges(3,5) * t342 + Icges(3,6) * t345;
t316 = rSges(4,1) * t336 + rSges(4,2) * t338;
t312 = qJD(1) + (-qJD(4) - qJD(5)) * t338;
t311 = t338 * t374 + t377;
t310 = -t338 * t376 + t375;
t309 = t338 * t375 - t376;
t308 = -t338 * t377 - t374;
t307 = t386 * t336;
t304 = rSges(3,3) * t343 + t346 * t365;
t303 = -rSges(3,3) * t346 + t343 * t365;
t296 = Icges(3,3) * t343 + t346 * t359;
t295 = -Icges(3,3) * t346 + t343 * t359;
t293 = t335 * t343 + t337 * t378;
t292 = -t335 * t378 + t337 * t343;
t291 = -t335 * t346 + t337 * t379;
t290 = -t335 * t379 - t337 * t346;
t289 = rSges(4,3) * t343 + t346 * t364;
t288 = -rSges(4,3) * t346 + t343 * t364;
t280 = -rSges(5,3) * t338 + (rSges(5,1) * t344 - rSges(5,2) * t341) * t336;
t279 = -Icges(5,5) * t338 + (Icges(5,1) * t344 - Icges(5,4) * t341) * t336;
t278 = -Icges(5,6) * t338 + (Icges(5,4) * t344 - Icges(5,2) * t341) * t336;
t277 = -Icges(5,3) * t338 + (Icges(5,5) * t344 - Icges(5,6) * t341) * t336;
t276 = -rSges(6,3) * t338 + (rSges(6,1) * t337 - rSges(6,2) * t335) * t336;
t275 = t343 * t369 + t302;
t274 = t346 * t369 + t301;
t272 = -Icges(6,5) * t338 + (Icges(6,1) * t337 - Icges(6,4) * t335) * t336;
t271 = -Icges(6,6) * t338 + (Icges(6,4) * t337 - Icges(6,2) * t335) * t336;
t270 = -Icges(6,3) * t338 + (Icges(6,5) * t337 - Icges(6,6) * t335) * t336;
t269 = pkin(3) * t377 + t346 * t367;
t268 = -pkin(3) * t376 + t343 * t367;
t267 = -qJD(1) * t303 - t323 * t371;
t266 = qJD(1) * t304 - t323 * t334;
t265 = rSges(5,1) * t311 + rSges(5,2) * t310 + rSges(5,3) * t380;
t264 = rSges(5,1) * t309 + rSges(5,2) * t308 + rSges(5,3) * t381;
t263 = Icges(5,1) * t311 + Icges(5,4) * t310 + Icges(5,5) * t380;
t262 = Icges(5,1) * t309 + Icges(5,4) * t308 + Icges(5,5) * t381;
t261 = Icges(5,4) * t311 + Icges(5,2) * t310 + Icges(5,6) * t380;
t260 = Icges(5,4) * t309 + Icges(5,2) * t308 + Icges(5,6) * t381;
t259 = Icges(5,5) * t311 + Icges(5,6) * t310 + Icges(5,3) * t380;
t258 = Icges(5,5) * t309 + Icges(5,6) * t308 + Icges(5,3) * t381;
t257 = (t303 * t343 + t304 * t346) * qJD(2);
t256 = rSges(6,1) * t293 + rSges(6,2) * t292 + rSges(6,3) * t380;
t255 = rSges(6,1) * t291 + rSges(6,2) * t290 + rSges(6,3) * t381;
t254 = Icges(6,1) * t293 + Icges(6,4) * t292 + Icges(6,5) * t380;
t253 = Icges(6,1) * t291 + Icges(6,4) * t290 + Icges(6,5) * t381;
t252 = Icges(6,4) * t293 + Icges(6,2) * t292 + Icges(6,6) * t380;
t251 = Icges(6,4) * t291 + Icges(6,2) * t290 + Icges(6,6) * t381;
t250 = Icges(6,5) * t293 + Icges(6,6) * t292 + Icges(6,3) * t380;
t249 = Icges(6,5) * t291 + Icges(6,6) * t290 + Icges(6,3) * t381;
t248 = -qJD(1) * t288 + t316 * t319 + t351;
t247 = qJD(1) * t289 - t316 * t318 + t354;
t246 = t288 * t318 - t289 * t319 + t373;
t245 = -t264 * t328 + t280 * t302 + t349;
t244 = t265 * t328 - t280 * t301 + t350;
t243 = t264 * t301 - t265 * t302 + t353;
t242 = -t255 * t312 - t268 * t328 + t275 * t276 + t302 * t307 + t349;
t241 = t256 * t312 + t269 * t328 - t274 * t276 - t301 * t307 + t350;
t240 = t255 * t274 - t256 * t275 + t268 * t301 - t269 * t302 + t353;
t1 = m(3) * (t257 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + ((t343 * t320 + t346 * t355) * qJD(1) + (t390 * t296 + (t357 * t346 + (-t295 + t356) * t343) * t346) * qJD(2)) * t334 / 0.2e1 - ((-t346 * t320 + t343 * t355) * qJD(1) + (t389 * t295 + (t356 * t343 + (-t296 + t357) * t346) * t343) * qJD(2)) * t371 / 0.2e1 + m(4) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + t318 * (t352 * t343 + t348 * t346) / 0.2e1 + t319 * (t348 * t343 - t352 * t346) / 0.2e1 + m(5) * (t243 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + t301 * ((t259 * t380 + t310 * t261 + t311 * t263) * t301 + (t258 * t380 + t260 * t310 + t262 * t311) * t302 + (t277 * t380 + t278 * t310 + t279 * t311) * t328) / 0.2e1 + t302 * ((t259 * t381 + t261 * t308 + t263 * t309) * t301 + (t258 * t381 + t308 * t260 + t309 * t262) * t302 + (t277 * t381 + t278 * t308 + t279 * t309) * t328) / 0.2e1 + t328 * ((-t258 * t302 - t259 * t301 - t277 * t328) * t338 + ((-t261 * t341 + t263 * t344) * t301 + (-t260 * t341 + t262 * t344) * t302 + (-t278 * t341 + t279 * t344) * t328) * t336) / 0.2e1 + m(6) * (t240 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + t274 * ((t250 * t380 + t292 * t252 + t293 * t254) * t274 + (t249 * t380 + t251 * t292 + t253 * t293) * t275 + (t270 * t380 + t271 * t292 + t272 * t293) * t312) / 0.2e1 + t275 * ((t250 * t381 + t252 * t290 + t254 * t291) * t274 + (t249 * t381 + t290 * t251 + t291 * t253) * t275 + (t270 * t381 + t271 * t290 + t272 * t291) * t312) / 0.2e1 + t312 * ((-t249 * t275 - t250 * t274 - t270 * t312) * t338 + ((-t252 * t335 + t254 * t337) * t274 + (-t251 * t335 + t253 * t337) * t275 + (-t271 * t335 + t272 * t337) * t312) * t336) / 0.2e1 + (m(2) * (t324 ^ 2 + t325 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t298 * t345 + t300 * t342) * t343 - (t297 * t345 + t299 * t342) * t346) * qJD(2) + (t285 * t338 + t287 * t336) * t318 + (t284 * t338 + t286 * t336) * t319 + (t338 * t314 + t336 * t315 + t345 * t321 + t342 * t322) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
