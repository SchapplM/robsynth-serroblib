% Calculate kinetic energy for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR14_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR14_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR14_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:29
% EndTime: 2019-12-31 18:34:30
% DurationCPUTime: 1.27s
% Computational Cost: add. (658->193), mult. (968->300), div. (0->0), fcn. (894->8), ass. (0->104)
t380 = Icges(4,3) + Icges(5,3);
t322 = qJ(3) + pkin(8);
t317 = sin(t322);
t318 = cos(t322);
t325 = sin(qJ(3));
t328 = cos(qJ(3));
t379 = Icges(4,5) * t325 + Icges(5,5) * t317 + Icges(4,6) * t328 + Icges(5,6) * t318;
t326 = sin(qJ(1));
t329 = cos(qJ(1));
t378 = t379 * t326 + t380 * t329;
t377 = t380 * t326 - t379 * t329;
t376 = Icges(4,5) * t328 + Icges(5,5) * t318 - Icges(4,6) * t325 - Icges(5,6) * t317;
t362 = Icges(5,4) * t318;
t302 = -Icges(5,2) * t317 + t362;
t363 = Icges(5,4) * t317;
t303 = Icges(5,1) * t318 - t363;
t364 = Icges(4,4) * t328;
t308 = -Icges(4,2) * t325 + t364;
t365 = Icges(4,4) * t325;
t309 = Icges(4,1) * t328 - t365;
t375 = t302 * t318 + t303 * t317 + t308 * t328 + t309 * t325;
t339 = Icges(5,2) * t318 + t363;
t276 = Icges(5,6) * t326 - t339 * t329;
t341 = Icges(5,1) * t317 + t362;
t278 = Icges(5,5) * t326 - t341 * t329;
t340 = Icges(4,2) * t328 + t365;
t286 = Icges(4,6) * t326 - t340 * t329;
t342 = Icges(4,1) * t325 + t364;
t288 = Icges(4,5) * t326 - t342 * t329;
t374 = t276 * t318 + t278 * t317 + t286 * t328 + t288 * t325;
t275 = Icges(5,6) * t329 + t339 * t326;
t277 = Icges(5,5) * t329 + t341 * t326;
t285 = Icges(4,6) * t329 + t340 * t326;
t287 = Icges(4,5) * t329 + t342 * t326;
t373 = -t275 * t318 - t277 * t317 - t285 * t328 - t287 * t325;
t369 = pkin(3) * t325;
t368 = pkin(3) * t328;
t361 = t318 * t326;
t360 = t318 * t329;
t324 = sin(qJ(5));
t359 = t326 * t324;
t327 = cos(qJ(5));
t358 = t326 * t327;
t357 = t329 * t324;
t356 = t329 * t327;
t306 = qJD(1) * (t329 * pkin(1) + t326 * qJ(2));
t355 = qJD(1) * t329 * pkin(6) + t306;
t354 = qJD(3) * t326;
t353 = qJD(3) * t329;
t352 = qJD(5) * t318;
t321 = qJD(2) * t326;
t351 = qJD(4) * t329 + t354 * t368 + t321;
t310 = t326 * pkin(1) - t329 * qJ(2);
t348 = -t326 * pkin(6) - t310;
t294 = qJ(4) * t329 + t326 * t369;
t347 = qJD(1) * t294 + qJD(4) * t326 + t355;
t293 = qJ(4) * t326 - t329 * t369;
t346 = -t293 + t348;
t345 = pkin(4) * t317 - pkin(7) * t318;
t344 = rSges(4,1) * t325 + rSges(4,2) * t328;
t343 = rSges(5,1) * t317 + rSges(5,2) * t318;
t314 = qJD(5) * t317 + qJD(1);
t313 = t329 * rSges(2,1) - t326 * rSges(2,2);
t312 = t328 * rSges(4,1) - t325 * rSges(4,2);
t311 = t326 * rSges(2,1) + t329 * rSges(2,2);
t305 = t318 * pkin(4) + t317 * pkin(7);
t304 = t318 * rSges(5,1) - t317 * rSges(5,2);
t300 = -t326 * t352 + t353;
t299 = t329 * t352 + t354;
t298 = -t317 * t356 + t359;
t297 = t317 * t357 + t358;
t296 = t317 * t358 + t357;
t295 = -t317 * t359 + t356;
t292 = t345 * t329;
t291 = t345 * t326;
t290 = t326 * rSges(4,3) - t344 * t329;
t289 = t329 * rSges(4,3) + t344 * t326;
t281 = t293 * t353;
t280 = t326 * rSges(5,3) - t343 * t329;
t279 = t329 * rSges(5,3) + t343 * t326;
t272 = t317 * rSges(6,3) + (rSges(6,1) * t327 - rSges(6,2) * t324) * t318;
t271 = Icges(6,5) * t317 + (Icges(6,1) * t327 - Icges(6,4) * t324) * t318;
t270 = Icges(6,6) * t317 + (Icges(6,4) * t327 - Icges(6,2) * t324) * t318;
t269 = Icges(6,3) * t317 + (Icges(6,5) * t327 - Icges(6,6) * t324) * t318;
t268 = t306 - qJD(2) * t329 + qJD(1) * (-t329 * rSges(3,2) + t326 * rSges(3,3));
t267 = t321 + (t326 * rSges(3,2) + t329 * rSges(3,3) - t310) * qJD(1);
t266 = (-t289 * t326 + t290 * t329) * qJD(3);
t265 = t298 * rSges(6,1) + t297 * rSges(6,2) + rSges(6,3) * t360;
t264 = t296 * rSges(6,1) + t295 * rSges(6,2) - rSges(6,3) * t361;
t263 = Icges(6,1) * t298 + Icges(6,4) * t297 + Icges(6,5) * t360;
t262 = Icges(6,1) * t296 + Icges(6,4) * t295 - Icges(6,5) * t361;
t261 = Icges(6,4) * t298 + Icges(6,2) * t297 + Icges(6,6) * t360;
t260 = Icges(6,4) * t296 + Icges(6,2) * t295 - Icges(6,6) * t361;
t259 = Icges(6,5) * t298 + Icges(6,6) * t297 + Icges(6,3) * t360;
t258 = Icges(6,5) * t296 + Icges(6,6) * t295 - Icges(6,3) * t361;
t257 = qJD(1) * t289 + (-qJD(3) * t312 - qJD(2)) * t329 + t355;
t256 = t312 * t354 + t321 + (-t290 + t348) * qJD(1);
t255 = qJD(1) * t279 + (-qJD(2) + (-t304 - t368) * qJD(3)) * t329 + t347;
t254 = t304 * t354 + (-t280 + t346) * qJD(1) + t351;
t253 = t281 + (t280 * t329 + (-t279 - t294) * t326) * qJD(3);
t252 = qJD(1) * t291 + t314 * t264 - t300 * t272 + (-qJD(2) + (-t305 - t368) * qJD(3)) * t329 + t347;
t251 = t305 * t354 - t314 * t265 + t299 * t272 + (t292 + t346) * qJD(1) + t351;
t250 = -t299 * t264 + t300 * t265 + t281 + (-t292 * t329 + (-t291 - t294) * t326) * qJD(3);
t1 = m(3) * (t267 ^ 2 + t268 ^ 2) / 0.2e1 + m(4) * (t256 ^ 2 + t257 ^ 2 + t266 ^ 2) / 0.2e1 + m(5) * (t253 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + m(6) * (t250 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + t300 * ((-t258 * t361 + t295 * t260 + t296 * t262) * t300 + (-t259 * t361 + t295 * t261 + t296 * t263) * t299 + (-t269 * t361 + t295 * t270 + t296 * t271) * t314) / 0.2e1 + t299 * ((t258 * t360 + t297 * t260 + t298 * t262) * t300 + (t259 * t360 + t297 * t261 + t298 * t263) * t299 + (t269 * t360 + t297 * t270 + t298 * t271) * t314) / 0.2e1 + t314 * ((t258 * t300 + t259 * t299 + t269 * t314) * t317 + ((-t260 * t324 + t262 * t327) * t300 + (-t261 * t324 + t263 * t327) * t299 + (-t270 * t324 + t271 * t327) * t314) * t318) / 0.2e1 + (((-t317 * t275 + t318 * t277 - t325 * t285 + t328 * t287) * t329 + (-t317 * t276 + t318 * t278 - t325 * t286 + t328 * t288) * t326) * qJD(3) + (-t317 * t302 + t318 * t303 - t325 * t308 + t328 * t309) * qJD(1)) * qJD(1) / 0.2e1 + ((t377 * t326 ^ 2 + (t373 * t329 + (-t374 + t378) * t326) * t329) * qJD(3) + (t326 * t376 - t329 * t375) * qJD(1)) * t354 / 0.2e1 + ((t378 * t329 ^ 2 + (t374 * t326 + (-t373 + t377) * t329) * t326) * qJD(3) + (t326 * t375 + t329 * t376) * qJD(1)) * t353 / 0.2e1 + (m(2) * (t311 ^ 2 + t313 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
