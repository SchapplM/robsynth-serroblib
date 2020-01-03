% Calculate kinetic energy for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR11_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR11_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:23
% EndTime: 2019-12-31 19:46:26
% DurationCPUTime: 2.67s
% Computational Cost: add. (755->235), mult. (1454->359), div. (0->0), fcn. (1407->8), ass. (0->122)
t457 = Icges(3,4) + Icges(4,6);
t456 = Icges(3,1) + Icges(4,2);
t455 = Icges(3,2) + Icges(4,3);
t383 = cos(qJ(2));
t454 = t457 * t383;
t381 = sin(qJ(2));
t453 = t457 * t381;
t452 = -Icges(4,4) + Icges(3,5);
t451 = Icges(4,5) - Icges(3,6);
t450 = t455 * t381 - t454;
t449 = -t456 * t383 + t453;
t448 = Icges(4,1) + Icges(3,3);
t382 = sin(qJ(1));
t384 = cos(qJ(1));
t447 = -t450 * t382 + t451 * t384;
t446 = t451 * t382 + t450 * t384;
t445 = t449 * t382 + t452 * t384;
t444 = t452 * t382 - t449 * t384;
t443 = -t455 * t383 - t453;
t442 = t456 * t381 + t454;
t441 = t451 * t381 + t452 * t383;
t440 = t441 * t382 - t448 * t384;
t439 = t448 * t382 + t441 * t384;
t438 = t452 * t381 - t451 * t383;
t437 = t443 * t381 + t442 * t383;
t436 = t446 * t381 + t444 * t383;
t435 = t447 * t381 + t445 * t383;
t378 = sin(pkin(8));
t431 = pkin(4) * t378;
t379 = cos(pkin(8));
t429 = pkin(4) * t379;
t424 = t381 * t382;
t423 = t381 * t384;
t422 = t382 * t383;
t421 = t383 * t384;
t403 = pkin(2) * t383 + qJ(3) * t381;
t348 = t403 * t382;
t367 = pkin(1) * t382 - pkin(6) * t384;
t419 = -t348 - t367;
t415 = qJD(3) * t381;
t371 = t384 * t415;
t414 = qJD(4) * t383;
t418 = t384 * t414 + t371;
t417 = qJD(2) * t382;
t416 = qJD(2) * t384;
t413 = qJD(5) * t383;
t349 = t403 * t384;
t355 = qJD(1) * (pkin(1) * t384 + pkin(6) * t382);
t412 = qJD(1) * t349 + t382 * t415 + t355;
t354 = -pkin(3) * t384 + qJ(4) * t422;
t411 = -t354 + t419;
t362 = pkin(2) * t381 - qJ(3) * t383;
t408 = -qJ(4) * t381 - t362;
t407 = qJD(2) * (rSges(4,2) * t381 + rSges(4,3) * t383 - t362);
t406 = -qJD(3) * t383 + t348 * t417 + t349 * t416;
t405 = rSges(3,1) * t383 - rSges(3,2) * t381;
t404 = -rSges(4,2) * t383 + rSges(4,3) * t381;
t353 = pkin(3) * t382 + qJ(4) * t421;
t402 = qJD(1) * t353 + t382 * t414 + t412;
t389 = qJD(2) * (-rSges(5,3) * t381 - (-rSges(5,1) * t378 - rSges(5,2) * t379) * t383 + t408);
t388 = qJD(2) * (-pkin(7) * t381 + t383 * t431 + t408);
t387 = pkin(7) * t383 + t381 * t431;
t386 = qJD(4) * t381 + t353 * t416 + t354 * t417 + t406;
t377 = pkin(8) + qJ(5);
t375 = cos(t377);
t374 = sin(t377);
t372 = qJD(5) * t381 + qJD(1);
t366 = rSges(2,1) * t384 - rSges(2,2) * t382;
t365 = rSges(2,1) * t382 + rSges(2,2) * t384;
t364 = rSges(3,1) * t381 + rSges(3,2) * t383;
t352 = t382 * t413 - t416;
t351 = t384 * t413 + t417;
t347 = t378 * t424 - t379 * t384;
t346 = t378 * t384 + t379 * t424;
t345 = t378 * t423 + t379 * t382;
t344 = -t378 * t382 + t379 * t423;
t339 = -rSges(4,1) * t384 + t404 * t382;
t338 = rSges(4,1) * t382 + t404 * t384;
t337 = rSges(3,3) * t382 + t405 * t384;
t336 = -rSges(3,3) * t384 + t405 * t382;
t321 = t374 * t424 - t375 * t384;
t320 = t374 * t384 + t375 * t424;
t319 = t374 * t423 + t375 * t382;
t318 = -t374 * t382 + t375 * t423;
t316 = Icges(5,5) * t381 + (-Icges(5,1) * t378 - Icges(5,4) * t379) * t383;
t315 = Icges(5,6) * t381 + (-Icges(5,4) * t378 - Icges(5,2) * t379) * t383;
t314 = Icges(5,3) * t381 + (-Icges(5,5) * t378 - Icges(5,6) * t379) * t383;
t313 = rSges(6,3) * t381 + (-rSges(6,1) * t374 - rSges(6,2) * t375) * t383;
t312 = Icges(6,5) * t381 + (-Icges(6,1) * t374 - Icges(6,4) * t375) * t383;
t311 = Icges(6,6) * t381 + (-Icges(6,4) * t374 - Icges(6,2) * t375) * t383;
t310 = Icges(6,3) * t381 + (-Icges(6,5) * t374 - Icges(6,6) * t375) * t383;
t309 = t382 * t387 - t384 * t429;
t308 = t382 * t429 + t384 * t387;
t307 = rSges(5,1) * t347 + rSges(5,2) * t346 + rSges(5,3) * t422;
t306 = rSges(5,1) * t345 + rSges(5,2) * t344 + rSges(5,3) * t421;
t305 = Icges(5,1) * t347 + Icges(5,4) * t346 + Icges(5,5) * t422;
t304 = Icges(5,1) * t345 + Icges(5,4) * t344 + Icges(5,5) * t421;
t303 = Icges(5,4) * t347 + Icges(5,2) * t346 + Icges(5,6) * t422;
t302 = Icges(5,4) * t345 + Icges(5,2) * t344 + Icges(5,6) * t421;
t301 = Icges(5,5) * t347 + Icges(5,6) * t346 + Icges(5,3) * t422;
t300 = Icges(5,5) * t345 + Icges(5,6) * t344 + Icges(5,3) * t421;
t299 = qJD(1) * t337 - t364 * t417 + t355;
t298 = -t364 * t416 + (-t336 - t367) * qJD(1);
t297 = (t336 * t382 + t337 * t384) * qJD(2);
t296 = rSges(6,1) * t321 + rSges(6,2) * t320 + rSges(6,3) * t422;
t295 = rSges(6,1) * t319 + rSges(6,2) * t318 + rSges(6,3) * t421;
t294 = Icges(6,1) * t321 + Icges(6,4) * t320 + Icges(6,5) * t422;
t293 = Icges(6,1) * t319 + Icges(6,4) * t318 + Icges(6,5) * t421;
t292 = Icges(6,4) * t321 + Icges(6,2) * t320 + Icges(6,6) * t422;
t291 = Icges(6,4) * t319 + Icges(6,2) * t318 + Icges(6,6) * t421;
t290 = Icges(6,5) * t321 + Icges(6,6) * t320 + Icges(6,3) * t422;
t289 = Icges(6,5) * t319 + Icges(6,6) * t318 + Icges(6,3) * t421;
t288 = qJD(1) * t338 + t382 * t407 + t412;
t287 = t371 + t384 * t407 + (-t339 + t419) * qJD(1);
t286 = (t338 * t384 + t339 * t382) * qJD(2) + t406;
t285 = qJD(1) * t306 + t382 * t389 + t402;
t284 = t384 * t389 + (-t307 + t411) * qJD(1) + t418;
t283 = (t306 * t384 + t307 * t382) * qJD(2) + t386;
t282 = qJD(1) * t308 + t295 * t372 - t313 * t351 + t382 * t388 + t402;
t281 = -t296 * t372 + t313 * t352 + t384 * t388 + (-t309 + t411) * qJD(1) + t418;
t280 = -t295 * t352 + t296 * t351 + (t308 * t384 + t309 * t382) * qJD(2) + t386;
t1 = m(3) * (t297 ^ 2 + t298 ^ 2 + t299 ^ 2) / 0.2e1 + m(4) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + m(5) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + m(6) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + t351 * ((t289 * t421 + t318 * t291 + t319 * t293) * t351 + (t290 * t421 + t292 * t318 + t294 * t319) * t352 + (t310 * t421 + t311 * t318 + t312 * t319) * t372) / 0.2e1 + t352 * ((t289 * t422 + t291 * t320 + t293 * t321) * t351 + (t290 * t422 + t320 * t292 + t321 * t294) * t352 + (t310 * t422 + t311 * t320 + t312 * t321) * t372) / 0.2e1 + t372 * ((t289 * t351 + t290 * t352 + t310 * t372) * t381 + ((-t291 * t375 - t293 * t374) * t351 + (-t292 * t375 - t294 * t374) * t352 + (-t311 * t375 - t312 * t374) * t372) * t383) / 0.2e1 + (m(2) * (t365 ^ 2 + t366 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((((t303 * t379 + t305 * t378 - t447) * t383 + (-t301 + t445) * t381) * t384 + ((-t302 * t379 - t304 * t378 - t446) * t383 + (t300 + t444) * t381) * t382) * qJD(2) + ((-t315 * t379 - t316 * t378 - t443) * t383 + (t314 + t442) * t381) * qJD(1)) * qJD(1) / 0.2e1 + (((-t301 * t421 - t303 * t344 - t305 * t345 + t435 * t384) * t384 + (t300 * t421 + t344 * t302 + t345 * t304 + (t436 - t440) * t384 + t439 * t382) * t382) * qJD(2) + (t314 * t421 + t315 * t344 + t316 * t345 + t438 * t382 + t437 * t384) * qJD(1)) * t417 / 0.2e1 - (((t300 * t422 + t302 * t346 + t347 * t304 + t436 * t382) * t382 + (-t301 * t422 - t346 * t303 - t347 * t305 + (t435 - t439) * t382 + t440 * t384) * t384) * qJD(2) + (t314 * t422 + t315 * t346 + t316 * t347 + t437 * t382 - t438 * t384) * qJD(1)) * t416 / 0.2e1;
T = t1;
