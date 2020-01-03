% Calculate kinetic energy for
% S5RPRRR9
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:17
% EndTime: 2019-12-31 19:07:18
% DurationCPUTime: 1.13s
% Computational Cost: add. (1089->199), mult. (1052->329), div. (0->0), fcn. (978->10), ass. (0->113)
t341 = cos(pkin(9));
t385 = t341 * pkin(2);
t339 = pkin(9) + qJ(3);
t331 = sin(t339);
t384 = Icges(4,4) * t331;
t332 = cos(t339);
t383 = Icges(4,4) * t332;
t333 = qJ(4) + t339;
t328 = sin(t333);
t382 = Icges(5,4) * t328;
t329 = cos(t333);
t381 = Icges(5,4) * t329;
t344 = sin(qJ(1));
t380 = t328 * t344;
t346 = cos(qJ(1));
t379 = t328 * t346;
t343 = sin(qJ(5));
t378 = t343 * t344;
t377 = t343 * t346;
t345 = cos(qJ(5));
t376 = t344 * t345;
t375 = t345 * t346;
t371 = pkin(3) * t332;
t277 = -pkin(7) * t346 + t344 * t371;
t278 = pkin(7) * t344 + t346 * t371;
t335 = qJD(3) * t344;
t369 = qJD(3) * t346;
t373 = t277 * t335 + t278 * t369;
t325 = pkin(1) * t344 - qJ(2) * t346;
t372 = pkin(6) * t346 - t344 * t385 - t325;
t323 = qJD(4) * t344 + t335;
t368 = qJD(5) * t328;
t367 = pkin(3) * qJD(3) * t331;
t366 = -t277 + t372;
t324 = (-qJD(3) - qJD(4)) * t346;
t365 = pkin(4) * t329 + pkin(8) * t328;
t321 = qJD(1) * (pkin(1) * t346 + qJ(2) * t344);
t364 = -qJD(2) * t346 + qJD(1) * (pkin(6) * t344 + t346 * t385) + t321;
t340 = sin(pkin(9));
t363 = rSges(3,1) * t341 - rSges(3,2) * t340;
t362 = rSges(4,1) * t332 - rSges(4,2) * t331;
t361 = rSges(5,1) * t329 - rSges(5,2) * t328;
t360 = Icges(4,1) * t332 - t384;
t359 = Icges(5,1) * t329 - t382;
t358 = -Icges(4,2) * t331 + t383;
t357 = -Icges(5,2) * t328 + t381;
t356 = Icges(4,5) * t332 - Icges(4,6) * t331;
t355 = Icges(5,5) * t329 - Icges(5,6) * t328;
t297 = -Icges(4,6) * t346 + t344 * t358;
t299 = -Icges(4,5) * t346 + t344 * t360;
t354 = t297 * t331 - t299 * t332;
t298 = Icges(4,6) * t344 + t346 * t358;
t300 = Icges(4,5) * t344 + t346 * t360;
t353 = -t298 * t331 + t300 * t332;
t317 = Icges(4,2) * t332 + t384;
t318 = Icges(4,1) * t331 + t383;
t352 = -t317 * t331 + t318 * t332;
t336 = qJD(2) * t344;
t351 = -t346 * t367 + t336;
t350 = (Icges(5,5) * t328 + Icges(5,6) * t329) * qJD(1) + (-Icges(5,3) * t346 + t344 * t355) * t324 + (Icges(5,3) * t344 + t346 * t355) * t323;
t349 = qJD(1) * t278 - t344 * t367 + t364;
t288 = -Icges(5,6) * t346 + t344 * t357;
t289 = Icges(5,6) * t344 + t346 * t357;
t290 = -Icges(5,5) * t346 + t344 * t359;
t291 = Icges(5,5) * t344 + t346 * t359;
t312 = Icges(5,2) * t329 + t382;
t313 = Icges(5,1) * t328 + t381;
t348 = (-t289 * t328 + t291 * t329) * t323 + (-t288 * t328 + t290 * t329) * t324 + (-t312 * t328 + t313 * t329) * qJD(1);
t327 = rSges(2,1) * t346 - rSges(2,2) * t344;
t326 = rSges(2,1) * t344 + rSges(2,2) * t346;
t322 = -qJD(5) * t329 + qJD(1);
t319 = rSges(4,1) * t331 + rSges(4,2) * t332;
t316 = Icges(4,5) * t331 + Icges(4,6) * t332;
t315 = pkin(4) * t328 - pkin(8) * t329;
t314 = rSges(5,1) * t328 + rSges(5,2) * t329;
t310 = t329 * t375 + t378;
t309 = -t329 * t377 + t376;
t308 = t329 * t376 - t377;
t307 = -t329 * t378 - t375;
t306 = t344 * t368 + t324;
t305 = t346 * t368 + t323;
t304 = t365 * t346;
t303 = t365 * t344;
t302 = rSges(4,3) * t344 + t346 * t362;
t301 = -rSges(4,3) * t346 + t344 * t362;
t296 = Icges(4,3) * t344 + t346 * t356;
t295 = -Icges(4,3) * t346 + t344 * t356;
t293 = rSges(5,3) * t344 + t346 * t361;
t292 = -rSges(5,3) * t346 + t344 * t361;
t284 = -rSges(6,3) * t329 + (rSges(6,1) * t345 - rSges(6,2) * t343) * t328;
t283 = -Icges(6,5) * t329 + (Icges(6,1) * t345 - Icges(6,4) * t343) * t328;
t282 = -Icges(6,6) * t329 + (Icges(6,4) * t345 - Icges(6,2) * t343) * t328;
t281 = -Icges(6,3) * t329 + (Icges(6,5) * t345 - Icges(6,6) * t343) * t328;
t280 = qJD(1) * t344 * rSges(3,3) + t321 + (qJD(1) * t363 - qJD(2)) * t346;
t279 = t336 + (t346 * rSges(3,3) - t344 * t363 - t325) * qJD(1);
t273 = rSges(6,1) * t310 + rSges(6,2) * t309 + rSges(6,3) * t379;
t272 = rSges(6,1) * t308 + rSges(6,2) * t307 + rSges(6,3) * t380;
t271 = Icges(6,1) * t310 + Icges(6,4) * t309 + Icges(6,5) * t379;
t270 = Icges(6,1) * t308 + Icges(6,4) * t307 + Icges(6,5) * t380;
t269 = Icges(6,4) * t310 + Icges(6,2) * t309 + Icges(6,6) * t379;
t268 = Icges(6,4) * t308 + Icges(6,2) * t307 + Icges(6,6) * t380;
t267 = Icges(6,5) * t310 + Icges(6,6) * t309 + Icges(6,3) * t379;
t266 = Icges(6,5) * t308 + Icges(6,6) * t307 + Icges(6,3) * t380;
t265 = (t301 * t344 + t302 * t346) * qJD(3);
t264 = qJD(1) * t302 - t319 * t335 + t364;
t263 = -t319 * t369 + t336 + (-t301 + t372) * qJD(1);
t262 = qJD(1) * t293 - t314 * t323 + t349;
t261 = t314 * t324 + (-t292 + t366) * qJD(1) + t351;
t260 = t292 * t323 - t293 * t324 + t373;
t259 = qJD(1) * t304 + t273 * t322 - t284 * t305 - t315 * t323 + t349;
t258 = -t272 * t322 + t284 * t306 + t315 * t324 + (-t303 + t366) * qJD(1) + t351;
t257 = t272 * t305 - t273 * t306 + t303 * t323 - t304 * t324 + t373;
t1 = m(3) * (t279 ^ 2 + t280 ^ 2) / 0.2e1 + m(4) * (t263 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + ((t344 * t316 + t346 * t352) * qJD(1) + (t344 ^ 2 * t296 + (t354 * t346 + (-t295 + t353) * t344) * t346) * qJD(3)) * t335 / 0.2e1 - ((-t346 * t316 + t344 * t352) * qJD(1) + (t346 ^ 2 * t295 + (t353 * t344 + (-t296 + t354) * t346) * t344) * qJD(3)) * t369 / 0.2e1 + m(5) * (t260 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + t323 * (t350 * t344 + t348 * t346) / 0.2e1 + t324 * (t348 * t344 - t350 * t346) / 0.2e1 + m(6) * (t257 ^ 2 + t258 ^ 2 + t259 ^ 2) / 0.2e1 + t305 * ((t267 * t379 + t309 * t269 + t310 * t271) * t305 + (t266 * t379 + t268 * t309 + t270 * t310) * t306 + (t281 * t379 + t282 * t309 + t283 * t310) * t322) / 0.2e1 + t306 * ((t267 * t380 + t269 * t307 + t271 * t308) * t305 + (t266 * t380 + t307 * t268 + t308 * t270) * t306 + (t281 * t380 + t282 * t307 + t283 * t308) * t322) / 0.2e1 + t322 * ((-t266 * t306 - t267 * t305 - t281 * t322) * t329 + ((-t269 * t343 + t271 * t345) * t305 + (-t268 * t343 + t270 * t345) * t306 + (-t282 * t343 + t283 * t345) * t322) * t328) / 0.2e1 + (((t298 * t332 + t300 * t331) * t344 - (t297 * t332 + t299 * t331) * t346) * qJD(3) + (t289 * t329 + t291 * t328) * t323 + (t288 * t329 + t290 * t328) * t324 + (t329 * t312 + t328 * t313 + t332 * t317 + t331 * t318) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t326 ^ 2 + t327 ^ 2) + Icges(2,3) + Icges(3,2) * t341 ^ 2 + (Icges(3,1) * t340 + 0.2e1 * Icges(3,4) * t341) * t340) * qJD(1) ^ 2 / 0.2e1;
T = t1;
