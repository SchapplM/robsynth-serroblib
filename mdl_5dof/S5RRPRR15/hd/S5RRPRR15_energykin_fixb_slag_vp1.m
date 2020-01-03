% Calculate kinetic energy for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR15_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR15_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR15_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:09
% EndTime: 2019-12-31 20:41:11
% DurationCPUTime: 2.73s
% Computational Cost: add. (791->234), mult. (1514->379), div. (0->0), fcn. (1467->8), ass. (0->126)
t456 = Icges(3,4) + Icges(4,6);
t455 = Icges(3,1) + Icges(4,2);
t454 = -Icges(3,2) - Icges(4,3);
t383 = cos(qJ(2));
t453 = t456 * t383;
t380 = sin(qJ(2));
t452 = t456 * t380;
t451 = -Icges(4,4) + Icges(3,5);
t450 = Icges(4,5) - Icges(3,6);
t449 = t454 * t380 + t453;
t448 = -t455 * t383 + t452;
t447 = Icges(4,1) + Icges(3,3);
t381 = sin(qJ(1));
t384 = cos(qJ(1));
t446 = t449 * t381 + t450 * t384;
t445 = -t450 * t381 + t449 * t384;
t444 = t448 * t381 + t451 * t384;
t443 = t451 * t381 - t448 * t384;
t442 = t454 * t383 - t452;
t441 = t455 * t380 + t453;
t440 = t450 * t380 + t451 * t383;
t439 = t440 * t381 - t447 * t384;
t438 = t447 * t381 + t440 * t384;
t437 = t451 * t380 - t450 * t383;
t436 = t442 * t380 + t441 * t383;
t435 = -t445 * t380 + t443 * t383;
t434 = t446 * t380 + t444 * t383;
t379 = sin(qJ(4));
t430 = pkin(4) * t379;
t382 = cos(qJ(4));
t428 = pkin(4) * t382;
t422 = t380 * t381;
t421 = t380 * t384;
t420 = t381 * t382;
t419 = t381 * t383;
t418 = t382 * t384;
t417 = t383 * t384;
t404 = pkin(2) * t383 + qJ(3) * t380;
t347 = t404 * t381;
t367 = pkin(1) * t381 - pkin(6) * t384;
t416 = -t347 - t367;
t375 = qJD(2) * t381;
t413 = qJD(4) * t383;
t350 = t384 * t413 + t375;
t415 = qJD(2) * t384;
t414 = qJD(3) * t380;
t412 = qJD(5) * t383;
t372 = qJD(4) * t380 + qJD(1);
t348 = t404 * t384;
t355 = qJD(1) * (pkin(1) * t384 + pkin(6) * t381);
t411 = qJD(1) * t348 + t381 * t414 + t355;
t362 = pkin(2) * t380 - qJ(3) * t383;
t408 = qJD(2) * (rSges(4,2) * t380 + rSges(4,3) * t383 - t362);
t351 = t381 * t413 - t415;
t407 = -qJD(3) * t383 + t347 * t375 + t348 * t415;
t406 = rSges(3,1) * t383 - rSges(3,2) * t380;
t405 = -rSges(4,2) * t383 + rSges(4,3) * t380;
t403 = qJD(2) * (-pkin(7) * t380 - t362);
t352 = pkin(3) * t381 + pkin(7) * t417;
t353 = -pkin(3) * t384 + pkin(7) * t419;
t390 = t352 * t415 + t353 * t375 + t407;
t389 = pkin(8) * t383 + t380 * t430;
t388 = qJD(1) * t352 + t381 * t403 + t411;
t371 = t384 * t414;
t387 = t371 + (-t353 + t416) * qJD(1) + t384 * t403;
t378 = qJ(4) + qJ(5);
t377 = cos(t378);
t376 = sin(t378);
t366 = rSges(2,1) * t384 - rSges(2,2) * t381;
t365 = rSges(2,1) * t381 + rSges(2,2) * t384;
t364 = rSges(3,1) * t380 + rSges(3,2) * t383;
t354 = qJD(5) * t380 + t372;
t346 = t379 * t422 - t418;
t345 = t379 * t384 + t380 * t420;
t344 = t379 * t421 + t420;
t343 = -t379 * t381 + t380 * t418;
t339 = pkin(8) * t380 - t383 * t430;
t338 = t376 * t422 - t377 * t384;
t337 = t376 * t384 + t377 * t422;
t336 = t376 * t421 + t377 * t381;
t335 = -t376 * t381 + t377 * t421;
t334 = -rSges(4,1) * t384 + t381 * t405;
t333 = rSges(4,1) * t381 + t384 * t405;
t332 = rSges(3,3) * t381 + t384 * t406;
t331 = rSges(5,3) * t380 + (-rSges(5,1) * t379 - rSges(5,2) * t382) * t383;
t330 = -rSges(3,3) * t384 + t381 * t406;
t319 = Icges(5,5) * t380 + (-Icges(5,1) * t379 - Icges(5,4) * t382) * t383;
t316 = Icges(5,6) * t380 + (-Icges(5,4) * t379 - Icges(5,2) * t382) * t383;
t313 = Icges(5,3) * t380 + (-Icges(5,5) * t379 - Icges(5,6) * t382) * t383;
t312 = t381 * t412 + t351;
t311 = t384 * t412 + t350;
t310 = rSges(6,3) * t380 + (-rSges(6,1) * t376 - rSges(6,2) * t377) * t383;
t309 = Icges(6,5) * t380 + (-Icges(6,1) * t376 - Icges(6,4) * t377) * t383;
t308 = Icges(6,6) * t380 + (-Icges(6,4) * t376 - Icges(6,2) * t377) * t383;
t307 = Icges(6,3) * t380 + (-Icges(6,5) * t376 - Icges(6,6) * t377) * t383;
t306 = t381 * t389 - t384 * t428;
t305 = t381 * t428 + t384 * t389;
t304 = rSges(5,1) * t346 + rSges(5,2) * t345 + rSges(5,3) * t419;
t303 = rSges(5,1) * t344 + rSges(5,2) * t343 + rSges(5,3) * t417;
t302 = Icges(5,1) * t346 + Icges(5,4) * t345 + Icges(5,5) * t419;
t301 = Icges(5,1) * t344 + Icges(5,4) * t343 + Icges(5,5) * t417;
t300 = Icges(5,4) * t346 + Icges(5,2) * t345 + Icges(5,6) * t419;
t299 = Icges(5,4) * t344 + Icges(5,2) * t343 + Icges(5,6) * t417;
t298 = Icges(5,5) * t346 + Icges(5,6) * t345 + Icges(5,3) * t419;
t297 = Icges(5,5) * t344 + Icges(5,6) * t343 + Icges(5,3) * t417;
t296 = qJD(1) * t332 - t364 * t375 + t355;
t295 = -t364 * t415 + (-t330 - t367) * qJD(1);
t294 = (t330 * t381 + t332 * t384) * qJD(2);
t293 = rSges(6,1) * t338 + rSges(6,2) * t337 + rSges(6,3) * t419;
t292 = rSges(6,1) * t336 + rSges(6,2) * t335 + rSges(6,3) * t417;
t291 = Icges(6,1) * t338 + Icges(6,4) * t337 + Icges(6,5) * t419;
t290 = Icges(6,1) * t336 + Icges(6,4) * t335 + Icges(6,5) * t417;
t289 = Icges(6,4) * t338 + Icges(6,2) * t337 + Icges(6,6) * t419;
t288 = Icges(6,4) * t336 + Icges(6,2) * t335 + Icges(6,6) * t417;
t287 = Icges(6,5) * t338 + Icges(6,6) * t337 + Icges(6,3) * t419;
t286 = Icges(6,5) * t336 + Icges(6,6) * t335 + Icges(6,3) * t417;
t285 = qJD(1) * t333 + t381 * t408 + t411;
t284 = t371 + t384 * t408 + (-t334 + t416) * qJD(1);
t283 = (t333 * t384 + t334 * t381) * qJD(2) + t407;
t282 = t303 * t372 - t331 * t350 + t388;
t281 = -t304 * t372 + t331 * t351 + t387;
t280 = -t303 * t351 + t304 * t350 + t390;
t279 = t292 * t354 + t305 * t372 - t310 * t311 - t339 * t350 + t388;
t278 = -t293 * t354 - t306 * t372 + t310 * t312 + t339 * t351 + t387;
t277 = -t292 * t312 + t293 * t311 - t305 * t351 + t306 * t350 + t390;
t1 = m(3) * (t294 ^ 2 + t295 ^ 2 + t296 ^ 2) / 0.2e1 + m(4) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + m(5) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + t350 * ((t297 * t417 + t343 * t299 + t344 * t301) * t350 + (t298 * t417 + t300 * t343 + t302 * t344) * t351 + (t313 * t417 + t316 * t343 + t319 * t344) * t372) / 0.2e1 + t351 * ((t297 * t419 + t299 * t345 + t301 * t346) * t350 + (t298 * t419 + t345 * t300 + t346 * t302) * t351 + (t313 * t419 + t316 * t345 + t319 * t346) * t372) / 0.2e1 + t372 * ((t297 * t350 + t298 * t351 + t313 * t372) * t380 + ((-t299 * t382 - t301 * t379) * t350 + (-t300 * t382 - t302 * t379) * t351 + (-t316 * t382 - t319 * t379) * t372) * t383) / 0.2e1 + m(6) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + t311 * ((t286 * t417 + t335 * t288 + t336 * t290) * t311 + (t287 * t417 + t289 * t335 + t291 * t336) * t312 + (t307 * t417 + t308 * t335 + t309 * t336) * t354) / 0.2e1 + t312 * ((t286 * t419 + t288 * t337 + t290 * t338) * t311 + (t287 * t419 + t337 * t289 + t338 * t291) * t312 + (t307 * t419 + t308 * t337 + t309 * t338) * t354) / 0.2e1 + t354 * ((t286 * t311 + t287 * t312 + t307 * t354) * t380 + ((-t288 * t377 - t290 * t376) * t311 + (-t289 * t377 - t291 * t376) * t312 + (-t308 * t377 - t309 * t376) * t354) * t383) / 0.2e1 + (m(2) * (t365 ^ 2 + t366 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t444 * t380 - t446 * t383) * t384 + (t443 * t380 + t445 * t383) * t381) * qJD(2) + (t441 * t380 - t442 * t383) * qJD(1)) * qJD(1) / 0.2e1 + ((t438 * t381 ^ 2 + (t434 * t384 + (t435 - t439) * t381) * t384) * qJD(2) + (t381 * t437 + t384 * t436) * qJD(1)) * t375 / 0.2e1 - ((t439 * t384 ^ 2 + (t435 * t381 + (t434 - t438) * t384) * t381) * qJD(2) + (t381 * t436 - t437 * t384) * qJD(1)) * t415 / 0.2e1;
T = t1;
