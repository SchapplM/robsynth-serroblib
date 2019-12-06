% Calculate kinetic energy for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:15
% EndTime: 2019-12-05 18:11:16
% DurationCPUTime: 0.97s
% Computational Cost: add. (996->177), mult. (858->289), div. (0->0), fcn. (728->10), ass. (0->108)
t332 = pkin(9) + qJ(3);
t325 = qJ(4) + t332;
t319 = sin(t325);
t381 = pkin(4) * t319;
t334 = cos(pkin(9));
t379 = t334 * pkin(2);
t323 = sin(t332);
t378 = Icges(4,4) * t323;
t324 = cos(t332);
t377 = Icges(4,4) * t324;
t376 = Icges(5,4) * t319;
t320 = cos(t325);
t375 = Icges(5,4) * t320;
t322 = qJ(5) + t325;
t316 = sin(t322);
t374 = Icges(6,4) * t316;
t317 = cos(t322);
t373 = Icges(6,4) * t317;
t336 = sin(qJ(1));
t337 = cos(qJ(1));
t368 = pkin(3) * t324;
t264 = -pkin(7) * t337 + t368 * t336;
t265 = pkin(7) * t336 + t368 * t337;
t327 = qJD(3) * t336;
t365 = qJD(3) * t337;
t371 = t264 * t327 + t265 * t365;
t313 = pkin(1) * t336 - qJ(2) * t337;
t370 = pkin(6) * t337 - t379 * t336 - t313;
t369 = pkin(4) * t320;
t311 = qJD(4) * t336 + t327;
t364 = -qJD(3) - qJD(4);
t363 = pkin(3) * qJD(3) * t323;
t362 = -t264 + t370;
t310 = qJD(1) * (pkin(1) * t337 + qJ(2) * t336);
t361 = -qJD(2) * t337 + qJD(1) * (pkin(6) * t336 + t379 * t337) + t310;
t333 = sin(pkin(9));
t360 = rSges(3,1) * t334 - rSges(3,2) * t333;
t359 = rSges(4,1) * t324 - rSges(4,2) * t323;
t358 = rSges(5,1) * t320 - rSges(5,2) * t319;
t357 = rSges(6,1) * t317 - rSges(6,2) * t316;
t356 = Icges(4,1) * t324 - t378;
t355 = Icges(5,1) * t320 - t376;
t354 = Icges(6,1) * t317 - t374;
t353 = -Icges(4,2) * t323 + t377;
t352 = -Icges(5,2) * t319 + t375;
t351 = -Icges(6,2) * t316 + t373;
t350 = Icges(4,5) * t324 - Icges(4,6) * t323;
t349 = Icges(5,5) * t320 - Icges(5,6) * t319;
t348 = Icges(6,5) * t317 - Icges(6,6) * t316;
t288 = -Icges(4,6) * t337 + t353 * t336;
t290 = -Icges(4,5) * t337 + t356 * t336;
t347 = t288 * t323 - t290 * t324;
t289 = Icges(4,6) * t336 + t353 * t337;
t291 = Icges(4,5) * t336 + t356 * t337;
t346 = -t289 * t323 + t291 * t324;
t306 = Icges(4,2) * t324 + t378;
t307 = Icges(4,1) * t323 + t377;
t345 = -t306 * t323 + t307 * t324;
t328 = qJD(2) * t336;
t344 = -t337 * t363 + t328;
t303 = qJD(5) * t336 + t311;
t304 = (-qJD(5) + t364) * t337;
t343 = (Icges(6,5) * t316 + Icges(6,6) * t317) * qJD(1) + (-Icges(6,3) * t337 + t348 * t336) * t304 + (Icges(6,3) * t336 + t348 * t337) * t303;
t312 = t364 * t337;
t342 = (Icges(5,5) * t319 + Icges(5,6) * t320) * qJD(1) + (-Icges(5,3) * t337 + t349 * t336) * t312 + (Icges(5,3) * t336 + t349 * t337) * t311;
t341 = qJD(1) * t265 - t336 * t363 + t361;
t270 = -Icges(6,6) * t337 + t351 * t336;
t271 = Icges(6,6) * t336 + t351 * t337;
t272 = -Icges(6,5) * t337 + t354 * t336;
t273 = Icges(6,5) * t336 + t354 * t337;
t296 = Icges(6,2) * t317 + t374;
t297 = Icges(6,1) * t316 + t373;
t340 = (-t271 * t316 + t273 * t317) * t303 + (-t270 * t316 + t272 * t317) * t304 + (-t296 * t316 + t297 * t317) * qJD(1);
t279 = -Icges(5,6) * t337 + t352 * t336;
t280 = Icges(5,6) * t336 + t352 * t337;
t281 = -Icges(5,5) * t337 + t355 * t336;
t282 = Icges(5,5) * t336 + t355 * t337;
t300 = Icges(5,2) * t320 + t376;
t301 = Icges(5,1) * t319 + t375;
t339 = (-t280 * t319 + t282 * t320) * t311 + (-t279 * t319 + t281 * t320) * t312 + (-t300 * t319 + t301 * t320) * qJD(1);
t315 = rSges(2,1) * t337 - rSges(2,2) * t336;
t314 = rSges(2,1) * t336 + rSges(2,2) * t337;
t308 = rSges(4,1) * t323 + rSges(4,2) * t324;
t305 = Icges(4,5) * t323 + Icges(4,6) * t324;
t302 = rSges(5,1) * t319 + rSges(5,2) * t320;
t298 = rSges(6,1) * t316 + rSges(6,2) * t317;
t293 = rSges(4,3) * t336 + t359 * t337;
t292 = -rSges(4,3) * t337 + t359 * t336;
t287 = Icges(4,3) * t336 + t350 * t337;
t286 = -Icges(4,3) * t337 + t350 * t336;
t284 = rSges(5,3) * t336 + t358 * t337;
t283 = -rSges(5,3) * t337 + t358 * t336;
t275 = rSges(6,3) * t336 + t357 * t337;
t274 = -rSges(6,3) * t337 + t357 * t336;
t267 = qJD(1) * t336 * rSges(3,3) + t310 + (qJD(1) * t360 - qJD(2)) * t337;
t266 = t328 + (t337 * rSges(3,3) - t360 * t336 - t313) * qJD(1);
t260 = pkin(8) * t336 + t369 * t337;
t259 = -pkin(8) * t337 + t369 * t336;
t258 = (t292 * t336 + t293 * t337) * qJD(3);
t257 = qJD(1) * t293 - t308 * t327 + t361;
t256 = -t308 * t365 + t328 + (-t292 + t370) * qJD(1);
t255 = qJD(1) * t284 - t302 * t311 + t341;
t254 = t302 * t312 + (-t283 + t362) * qJD(1) + t344;
t253 = t283 * t311 - t284 * t312 + t371;
t252 = -t311 * t381 - t298 * t303 + (t260 + t275) * qJD(1) + t341;
t251 = t312 * t381 + t298 * t304 + (-t259 - t274 + t362) * qJD(1) + t344;
t250 = t259 * t311 - t260 * t312 + t274 * t303 - t275 * t304 + t371;
t1 = m(3) * (t266 ^ 2 + t267 ^ 2) / 0.2e1 + m(4) * (t256 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + ((t336 * t305 + t345 * t337) * qJD(1) + (t336 ^ 2 * t287 + (t347 * t337 + (-t286 + t346) * t336) * t337) * qJD(3)) * t327 / 0.2e1 - ((-t337 * t305 + t345 * t336) * qJD(1) + (t337 ^ 2 * t286 + (t346 * t336 + (-t287 + t347) * t337) * t336) * qJD(3)) * t365 / 0.2e1 + m(5) * (t253 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t311 * (t342 * t336 + t339 * t337) / 0.2e1 + t312 * (t339 * t336 - t342 * t337) / 0.2e1 + m(6) * (t250 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + t303 * (t343 * t336 + t340 * t337) / 0.2e1 + t304 * (t340 * t336 - t343 * t337) / 0.2e1 + (m(2) * (t314 ^ 2 + t315 ^ 2) + Icges(2,3) + Icges(3,2) * t334 ^ 2 + (Icges(3,1) * t333 + 0.2e1 * Icges(3,4) * t334) * t333) * qJD(1) ^ 2 / 0.2e1 + (((t289 * t324 + t291 * t323) * t336 - (t288 * t324 + t323 * t290) * t337) * qJD(3) + (t280 * t320 + t282 * t319) * t311 + (t279 * t320 + t281 * t319) * t312 + (t271 * t317 + t273 * t316) * t303 + (t270 * t317 + t272 * t316) * t304 + (t317 * t296 + t316 * t297 + t320 * t300 + t319 * t301 + t324 * t306 + t323 * t307) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
