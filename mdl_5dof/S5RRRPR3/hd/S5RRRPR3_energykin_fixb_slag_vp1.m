% Calculate kinetic energy for
% S5RRRPR3
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
% m [6x1]
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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:25
% EndTime: 2022-01-20 11:42:26
% DurationCPUTime: 1.00s
% Computational Cost: add. (1014->164), mult. (764->260), div. (0->0), fcn. (634->10), ass. (0->103)
t369 = Icges(4,3) + Icges(5,3);
t305 = qJ(3) + pkin(9);
t297 = sin(t305);
t298 = cos(t305);
t308 = sin(qJ(3));
t310 = cos(qJ(3));
t368 = Icges(4,5) * t310 + Icges(5,5) * t298 - Icges(4,6) * t308 - Icges(5,6) * t297;
t306 = qJ(1) + qJ(2);
t300 = sin(t306);
t301 = cos(t306);
t367 = t300 * t368 - t301 * t369;
t366 = t300 * t369 + t301 * t368;
t365 = Icges(4,5) * t308 + Icges(5,5) * t297 + Icges(4,6) * t310 + Icges(5,6) * t298;
t351 = Icges(5,4) * t297;
t279 = Icges(5,2) * t298 + t351;
t350 = Icges(5,4) * t298;
t280 = Icges(5,1) * t297 + t350;
t353 = Icges(4,4) * t308;
t287 = Icges(4,2) * t310 + t353;
t352 = Icges(4,4) * t310;
t288 = Icges(4,1) * t308 + t352;
t364 = -t279 * t297 + t280 * t298 - t287 * t308 + t288 * t310;
t327 = -Icges(5,2) * t297 + t350;
t258 = Icges(5,6) * t300 + t301 * t327;
t330 = Icges(5,1) * t298 - t351;
t260 = Icges(5,5) * t300 + t301 * t330;
t328 = -Icges(4,2) * t308 + t352;
t266 = Icges(4,6) * t300 + t301 * t328;
t331 = Icges(4,1) * t310 - t353;
t268 = Icges(4,5) * t300 + t301 * t331;
t363 = -t258 * t297 + t260 * t298 - t266 * t308 + t268 * t310;
t257 = -Icges(5,6) * t301 + t300 * t327;
t259 = -Icges(5,5) * t301 + t300 * t330;
t265 = -Icges(4,6) * t301 + t300 * t328;
t267 = -Icges(4,5) * t301 + t300 * t331;
t362 = t257 * t297 - t259 * t298 + t265 * t308 - t267 * t310;
t357 = pkin(3) * t308;
t356 = t310 * pkin(3);
t354 = pkin(1) * qJD(1);
t299 = qJ(5) + t305;
t293 = sin(t299);
t349 = Icges(6,4) * t293;
t294 = cos(t299);
t348 = Icges(6,4) * t294;
t251 = -qJ(4) * t301 + t300 * t356;
t252 = qJ(4) * t300 + t301 * t356;
t341 = qJD(3) * t301;
t342 = qJD(3) * t300;
t347 = t251 * t342 + t252 * t341;
t284 = t300 * pkin(2) - t301 * pkin(7);
t346 = -t251 - t284;
t311 = cos(qJ(1));
t296 = t311 * t354;
t304 = qJD(1) + qJD(2);
t345 = t304 * (t301 * pkin(2) + t300 * pkin(7)) + t296;
t344 = pkin(4) * t298;
t340 = qJD(3) + qJD(5);
t309 = sin(qJ(1));
t339 = t309 * t354;
t336 = qJD(4) * t300 - t339;
t335 = rSges(4,1) * t310 - rSges(4,2) * t308;
t334 = rSges(5,1) * t298 - rSges(5,2) * t297;
t333 = rSges(6,1) * t294 - rSges(6,2) * t293;
t332 = qJD(3) * (-rSges(5,1) * t297 - rSges(5,2) * t298 - t357);
t329 = Icges(6,1) * t294 - t349;
t326 = -Icges(6,2) * t293 + t348;
t323 = Icges(6,5) * t294 - Icges(6,6) * t293;
t316 = -qJD(4) * t301 + t252 * t304 + t345;
t315 = qJD(3) * (-pkin(4) * t297 - t357);
t282 = t340 * t300;
t283 = t340 * t301;
t314 = -(-Icges(6,3) * t301 + t300 * t323) * t283 + (Icges(6,3) * t300 + t301 * t323) * t282 + (Icges(6,5) * t293 + Icges(6,6) * t294) * t304;
t247 = -Icges(6,6) * t301 + t300 * t326;
t248 = Icges(6,6) * t300 + t301 * t326;
t249 = -Icges(6,5) * t301 + t300 * t329;
t250 = Icges(6,5) * t300 + t301 * t329;
t275 = Icges(6,2) * t294 + t349;
t276 = Icges(6,1) * t293 + t348;
t313 = (-t248 * t293 + t250 * t294) * t282 - (-t247 * t293 + t249 * t294) * t283 + (-t275 * t293 + t276 * t294) * t304;
t291 = rSges(2,1) * t311 - rSges(2,2) * t309;
t290 = rSges(2,1) * t309 + rSges(2,2) * t311;
t289 = rSges(4,1) * t308 + rSges(4,2) * t310;
t277 = rSges(6,1) * t293 + rSges(6,2) * t294;
t272 = t296 + t304 * (rSges(3,1) * t301 - rSges(3,2) * t300);
t271 = -t339 - t304 * (rSges(3,1) * t300 + rSges(3,2) * t301);
t270 = t300 * rSges(4,3) + t301 * t335;
t269 = -t301 * rSges(4,3) + t300 * t335;
t262 = t300 * rSges(5,3) + t301 * t334;
t261 = -t301 * rSges(5,3) + t300 * t334;
t254 = t300 * rSges(6,3) + t301 * t333;
t253 = -t301 * rSges(6,3) + t300 * t333;
t241 = pkin(8) * t300 + t301 * t344;
t240 = -pkin(8) * t301 + t300 * t344;
t239 = (t269 * t300 + t270 * t301) * qJD(3);
t238 = t270 * t304 - t289 * t342 + t345;
t237 = -t339 - t289 * t341 + (-t269 - t284) * t304;
t236 = t304 * t262 + t300 * t332 + t316;
t235 = t301 * t332 + (-t261 + t346) * t304 + t336;
t234 = (t261 * t300 + t262 * t301) * qJD(3) + t347;
t233 = -t282 * t277 + (t241 + t254) * t304 + t300 * t315 + t316;
t232 = -t283 * t277 + t301 * t315 + (-t240 - t253 + t346) * t304 + t336;
t231 = t282 * t253 + t283 * t254 + (t240 * t300 + t241 * t301) * qJD(3) + t347;
t1 = m(3) * (t271 ^ 2 + t272 ^ 2) / 0.2e1 + t304 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t237 ^ 2 + t238 ^ 2 + t239 ^ 2) / 0.2e1 + m(5) * (t234 ^ 2 + t235 ^ 2 + t236 ^ 2) / 0.2e1 + m(6) * (t231 ^ 2 + t232 ^ 2 + t233 ^ 2) / 0.2e1 + t282 * (t314 * t300 + t313 * t301) / 0.2e1 - t283 * (t313 * t300 - t314 * t301) / 0.2e1 + (m(2) * (t290 ^ 2 + t291 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t300 * t365 + t301 * t364) * t304 + (t366 * t300 ^ 2 + (t362 * t301 + (t363 - t367) * t300) * t301) * qJD(3)) * t342 / 0.2e1 - ((t300 * t364 - t301 * t365) * t304 + (t367 * t301 ^ 2 + (t363 * t300 + (t362 - t366) * t301) * t300) * qJD(3)) * t341 / 0.2e1 + ((t248 * t294 + t250 * t293) * t282 - (t247 * t294 + t249 * t293) * t283 + ((-t257 * t298 - t259 * t297 - t265 * t310 - t267 * t308) * t301 + (t258 * t298 + t260 * t297 + t266 * t310 + t268 * t308) * t300) * qJD(3) + (t294 * t275 + t293 * t276 + t298 * t279 + t297 * t280 + t310 * t287 + t308 * t288) * t304) * t304 / 0.2e1;
T = t1;
