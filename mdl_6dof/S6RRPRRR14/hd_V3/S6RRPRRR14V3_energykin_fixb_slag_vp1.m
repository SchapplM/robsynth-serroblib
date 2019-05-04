% Calculate kinetic energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR14V3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:23
% EndTime: 2019-04-12 15:03:26
% DurationCPUTime: 3.68s
% Computational Cost: add. (1233->277), mult. (3049->448), div. (0->0), fcn. (3504->10), ass. (0->138)
t459 = Icges(3,4) - Icges(4,5);
t458 = Icges(3,1) + Icges(4,1);
t457 = Icges(3,2) + Icges(4,3);
t391 = cos(qJ(2));
t456 = t459 * t391;
t387 = sin(qJ(2));
t455 = t459 * t387;
t454 = Icges(4,4) + Icges(3,5);
t453 = Icges(3,6) - Icges(4,6);
t452 = t457 * t387 - t456;
t451 = t458 * t391 - t455;
t450 = Icges(4,2) + Icges(3,3);
t388 = sin(qJ(1));
t392 = cos(qJ(1));
t449 = t452 * t388 + t453 * t392;
t448 = -t453 * t388 + t452 * t392;
t447 = -t451 * t388 + t454 * t392;
t446 = t454 * t388 + t451 * t392;
t445 = -t457 * t391 - t455;
t444 = t458 * t387 + t456;
t443 = -t453 * t387 + t454 * t391;
t420 = qJ(3) * qJD(2);
t442 = qJD(3) * t387 + t391 * t420;
t441 = t443 * t388 - t392 * t450;
t440 = t388 * t450 + t443 * t392;
t439 = t454 * t387 + t453 * t391;
t438 = t387 * t445 + t391 * t444;
t437 = t387 * t448 + t391 * t446;
t436 = -t387 * t449 + t391 * t447;
t435 = t388 ^ 2;
t434 = t392 ^ 2;
t432 = cos(qJ(5));
t386 = sin(qJ(4));
t426 = t386 * t387;
t425 = t387 * t388;
t424 = t387 * t392;
t423 = t391 * t388;
t422 = t391 * t392;
t421 = t442 * t392;
t383 = qJD(2) * t388;
t417 = qJD(4) * t387;
t357 = t392 * t417 + t383;
t419 = qJD(2) * t392;
t414 = qJ(3) * qJD(1) * t387;
t416 = t442 * t388 + t392 * t414;
t390 = cos(qJ(4));
t354 = t386 * t422 - t388 * t390;
t320 = qJD(5) * t354 + t357;
t415 = t387 * t432;
t358 = t388 * t417 - t419;
t381 = -qJD(4) * t391 + qJD(1);
t352 = t386 * t423 + t390 * t392;
t321 = qJD(5) * t352 + t358;
t409 = -qJD(3) * t391 + (t434 + t435) * t387 * t420;
t408 = rSges(3,1) * t391 - rSges(3,2) * t387;
t407 = rSges(4,1) * t391 + rSges(4,3) * t387;
t356 = qJD(5) * t426 + t381;
t394 = -t388 * t414 + t421;
t389 = cos(qJ(6));
t385 = sin(qJ(5));
t384 = sin(qJ(6));
t368 = rSges(2,1) * t392 - rSges(2,2) * t388;
t367 = rSges(2,1) * t388 + rSges(2,2) * t392;
t366 = rSges(3,1) * t387 + rSges(3,2) * t391;
t365 = rSges(4,1) * t387 - rSges(4,3) * t391;
t355 = t386 * t388 + t390 * t422;
t353 = -t386 * t392 + t390 * t423;
t351 = -t391 * t385 + t390 * t415;
t350 = t387 * t390 * t385 + t391 * t432;
t347 = rSges(3,3) * t388 + t392 * t408;
t346 = rSges(4,2) * t388 + t392 * t407;
t345 = -rSges(3,3) * t392 + t388 * t408;
t344 = -rSges(4,2) * t392 + t388 * t407;
t343 = -rSges(5,3) * t391 + (rSges(5,1) * t390 - rSges(5,2) * t386) * t387;
t338 = -Icges(5,5) * t391 + (Icges(5,1) * t390 - Icges(5,4) * t386) * t387;
t333 = -Icges(5,6) * t391 + (Icges(5,4) * t390 - Icges(5,2) * t386) * t387;
t328 = -Icges(5,3) * t391 + (Icges(5,5) * t390 - Icges(5,6) * t386) * t387;
t327 = t355 * t432 + t385 * t424;
t326 = t355 * t385 - t392 * t415;
t325 = t353 * t432 + t385 * t425;
t324 = t353 * t385 - t388 * t415;
t323 = t351 * t389 + t384 * t426;
t322 = -t351 * t384 + t389 * t426;
t319 = qJD(6) * t350 + t356;
t318 = -qJD(1) * t345 - t366 * t419;
t317 = qJD(1) * t347 - t366 * t383;
t316 = rSges(5,1) * t355 - rSges(5,2) * t354 + rSges(5,3) * t424;
t315 = rSges(5,1) * t353 - rSges(5,2) * t352 + rSges(5,3) * t425;
t314 = rSges(6,1) * t351 - rSges(6,2) * t350 + rSges(6,3) * t426;
t313 = Icges(5,1) * t355 - Icges(5,4) * t354 + Icges(5,5) * t424;
t312 = Icges(5,1) * t353 - Icges(5,4) * t352 + Icges(5,5) * t425;
t311 = Icges(6,1) * t351 - Icges(6,4) * t350 + Icges(6,5) * t426;
t310 = Icges(5,4) * t355 - Icges(5,2) * t354 + Icges(5,6) * t424;
t309 = Icges(5,4) * t353 - Icges(5,2) * t352 + Icges(5,6) * t425;
t308 = Icges(6,4) * t351 - Icges(6,2) * t350 + Icges(6,6) * t426;
t307 = Icges(5,5) * t355 - Icges(5,6) * t354 + Icges(5,3) * t424;
t306 = Icges(5,5) * t353 - Icges(5,6) * t352 + Icges(5,3) * t425;
t305 = Icges(6,5) * t351 - Icges(6,6) * t350 + Icges(6,3) * t426;
t304 = (t345 * t388 + t347 * t392) * qJD(2);
t303 = t327 * t389 + t354 * t384;
t302 = -t327 * t384 + t354 * t389;
t301 = t325 * t389 + t352 * t384;
t300 = -t325 * t384 + t352 * t389;
t299 = qJD(6) * t324 + t321;
t298 = qJD(6) * t326 + t320;
t297 = -t365 * t419 + (-qJ(3) * t425 - t344) * qJD(1) + t421;
t296 = qJD(1) * t346 - t365 * t383 + t416;
t295 = (t344 * t388 + t346 * t392) * qJD(2) + t409;
t294 = rSges(6,1) * t327 - rSges(6,2) * t326 + rSges(6,3) * t354;
t293 = rSges(6,1) * t325 - rSges(6,2) * t324 + rSges(6,3) * t352;
t292 = rSges(7,1) * t323 + rSges(7,2) * t322 + rSges(7,3) * t350;
t291 = Icges(6,1) * t327 - Icges(6,4) * t326 + Icges(6,5) * t354;
t290 = Icges(6,1) * t325 - Icges(6,4) * t324 + Icges(6,5) * t352;
t289 = Icges(7,1) * t323 + Icges(7,4) * t322 + Icges(7,5) * t350;
t288 = Icges(6,4) * t327 - Icges(6,2) * t326 + Icges(6,6) * t354;
t287 = Icges(6,4) * t325 - Icges(6,2) * t324 + Icges(6,6) * t352;
t286 = Icges(7,4) * t323 + Icges(7,2) * t322 + Icges(7,6) * t350;
t285 = Icges(6,5) * t327 - Icges(6,6) * t326 + Icges(6,3) * t354;
t284 = Icges(6,5) * t325 - Icges(6,6) * t324 + Icges(6,3) * t352;
t283 = Icges(7,5) * t323 + Icges(7,6) * t322 + Icges(7,3) * t350;
t282 = -t315 * t381 + t343 * t358 + t394;
t281 = t316 * t381 - t343 * t357 + t416;
t280 = t315 * t357 - t316 * t358 + t409;
t279 = rSges(7,1) * t303 + rSges(7,2) * t302 + rSges(7,3) * t326;
t278 = rSges(7,1) * t301 + rSges(7,2) * t300 + rSges(7,3) * t324;
t277 = Icges(7,1) * t303 + Icges(7,4) * t302 + Icges(7,5) * t326;
t276 = Icges(7,1) * t301 + Icges(7,4) * t300 + Icges(7,5) * t324;
t275 = Icges(7,4) * t303 + Icges(7,2) * t302 + Icges(7,6) * t326;
t274 = Icges(7,4) * t301 + Icges(7,2) * t300 + Icges(7,6) * t324;
t273 = Icges(7,5) * t303 + Icges(7,6) * t302 + Icges(7,3) * t326;
t272 = Icges(7,5) * t301 + Icges(7,6) * t300 + Icges(7,3) * t324;
t271 = -t293 * t356 + t314 * t321 + t394;
t270 = t294 * t356 - t314 * t320 + t416;
t269 = t293 * t320 - t294 * t321 + t409;
t268 = -t278 * t319 + t292 * t299 + t394;
t267 = t279 * t319 - t292 * t298 + t416;
t266 = t278 * t298 - t279 * t299 + t409;
t1 = m(6) * (t269 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + m(7) * (t266 ^ 2 + t267 ^ 2 + t268 ^ 2) / 0.2e1 + t357 * ((t307 * t424 - t354 * t310 + t355 * t313) * t357 + (t306 * t424 - t309 * t354 + t312 * t355) * t358 + (t328 * t424 - t333 * t354 + t338 * t355) * t381) / 0.2e1 + t358 * ((t307 * t425 - t310 * t352 + t313 * t353) * t357 + (t306 * t425 - t352 * t309 + t353 * t312) * t358 + (t328 * t425 - t333 * t352 + t338 * t353) * t381) / 0.2e1 + t381 * ((-t306 * t358 - t307 * t357 - t328 * t381) * t391 + ((-t310 * t386 + t313 * t390) * t357 + (-t309 * t386 + t312 * t390) * t358 + (-t333 * t386 + t338 * t390) * t381) * t387) / 0.2e1 + t320 * ((t354 * t285 - t326 * t288 + t327 * t291) * t320 + (t284 * t354 - t287 * t326 + t290 * t327) * t321 + (t305 * t354 - t308 * t326 + t311 * t327) * t356) / 0.2e1 + t321 * ((t285 * t352 - t288 * t324 + t291 * t325) * t320 + (t352 * t284 - t324 * t287 + t325 * t290) * t321 + (t305 * t352 - t308 * t324 + t311 * t325) * t356) / 0.2e1 + t356 * ((t285 * t426 - t288 * t350 + t291 * t351) * t320 + (t284 * t426 - t287 * t350 + t290 * t351) * t321 + (t305 * t426 - t350 * t308 + t351 * t311) * t356) / 0.2e1 + t298 * ((t326 * t273 + t302 * t275 + t303 * t277) * t298 + (t272 * t326 + t274 * t302 + t276 * t303) * t299 + (t283 * t326 + t286 * t302 + t289 * t303) * t319) / 0.2e1 + t299 * ((t273 * t324 + t275 * t300 + t277 * t301) * t298 + (t324 * t272 + t300 * t274 + t301 * t276) * t299 + (t283 * t324 + t286 * t300 + t289 * t301) * t319) / 0.2e1 + t319 * ((t273 * t350 + t275 * t322 + t277 * t323) * t298 + (t272 * t350 + t274 * t322 + t276 * t323) * t299 + (t350 * t283 + t322 * t286 + t323 * t289) * t319) / 0.2e1 + m(3) * (t304 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(4) * (t295 ^ 2 + t296 ^ 2 + t297 ^ 2) / 0.2e1 + m(5) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t367 ^ 2 + t368 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t387 * t447 + t391 * t449) * t392 + (t387 * t446 - t391 * t448) * t388) * qJD(2) + (t444 * t387 - t445 * t391) * qJD(1)) * qJD(1) / 0.2e1 + ((t440 * t435 + (t436 * t392 + (t437 - t441) * t388) * t392) * qJD(2) + (t388 * t439 + t392 * t438) * qJD(1)) * t383 / 0.2e1 - ((t441 * t434 + (t437 * t388 + (t436 - t440) * t392) * t388) * qJD(2) + (t388 * t438 - t392 * t439) * qJD(1)) * t419 / 0.2e1;
T  = t1;
