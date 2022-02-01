% Calculate kinetic energy for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:37
% EndTime: 2022-01-20 12:07:38
% DurationCPUTime: 0.90s
% Computational Cost: add. (1062->169), mult. (788->278), div. (0->0), fcn. (658->10), ass. (0->107)
t305 = qJ(3) + qJ(4);
t297 = sin(t305);
t352 = pkin(4) * t297;
t309 = cos(qJ(3));
t351 = t309 * pkin(3);
t349 = pkin(1) * qJD(1);
t307 = sin(qJ(3));
t348 = Icges(4,4) * t307;
t347 = Icges(4,4) * t309;
t346 = Icges(5,4) * t297;
t299 = cos(t305);
t345 = Icges(5,4) * t299;
t301 = qJ(5) + t305;
t293 = sin(t301);
t344 = Icges(6,4) * t293;
t294 = cos(t301);
t343 = Icges(6,4) * t294;
t306 = qJ(1) + qJ(2);
t298 = sin(t306);
t300 = cos(t306);
t242 = -pkin(8) * t300 + t351 * t298;
t243 = pkin(8) * t298 + t351 * t300;
t292 = qJD(3) * t298;
t337 = qJD(3) * t300;
t342 = t242 * t292 + t243 * t337;
t283 = pkin(2) * t298 - pkin(7) * t300;
t341 = -t242 - t283;
t310 = cos(qJ(1));
t296 = t310 * t349;
t303 = qJD(1) + qJD(2);
t340 = t303 * (pkin(2) * t300 + pkin(7) * t298) + t296;
t339 = pkin(4) * t299;
t277 = qJD(4) * t298 + t292;
t336 = -qJD(3) - qJD(4);
t335 = pkin(3) * qJD(3) * t307;
t308 = sin(qJ(1));
t334 = t308 * t349;
t333 = rSges(4,1) * t309 - rSges(4,2) * t307;
t332 = rSges(5,1) * t299 - rSges(5,2) * t297;
t331 = rSges(6,1) * t294 - rSges(6,2) * t293;
t330 = Icges(4,1) * t309 - t348;
t329 = Icges(5,1) * t299 - t346;
t328 = Icges(6,1) * t294 - t344;
t327 = -Icges(4,2) * t307 + t347;
t326 = -Icges(5,2) * t297 + t345;
t325 = -Icges(6,2) * t293 + t343;
t324 = Icges(4,5) * t309 - Icges(4,6) * t307;
t323 = Icges(5,5) * t299 - Icges(5,6) * t297;
t322 = Icges(6,5) * t294 - Icges(6,6) * t293;
t262 = -Icges(4,6) * t300 + t298 * t327;
t264 = -Icges(4,5) * t300 + t298 * t330;
t321 = t262 * t307 - t264 * t309;
t263 = Icges(4,6) * t298 + t300 * t327;
t265 = Icges(4,5) * t298 + t300 * t330;
t320 = -t263 * t307 + t265 * t309;
t286 = Icges(4,2) * t309 + t348;
t287 = Icges(4,1) * t307 + t347;
t319 = -t286 * t307 + t287 * t309;
t318 = t303 * t243 - t298 * t335 + t340;
t270 = qJD(5) * t298 + t277;
t271 = (-qJD(5) + t336) * t300;
t317 = (-Icges(6,3) * t300 + t298 * t322) * t271 + (Icges(6,3) * t298 + t300 * t322) * t270 + (Icges(6,5) * t293 + Icges(6,6) * t294) * t303;
t278 = t336 * t300;
t316 = (-Icges(5,3) * t300 + t298 * t323) * t278 + (Icges(5,3) * t298 + t300 * t323) * t277 + (Icges(5,5) * t297 + Icges(5,6) * t299) * t303;
t315 = -t300 * t335 - t334;
t246 = -Icges(6,6) * t300 + t298 * t325;
t247 = Icges(6,6) * t298 + t300 * t325;
t248 = -Icges(6,5) * t300 + t298 * t328;
t249 = Icges(6,5) * t298 + t300 * t328;
t274 = Icges(6,2) * t294 + t344;
t275 = Icges(6,1) * t293 + t343;
t314 = (-t247 * t293 + t249 * t294) * t270 + (-t246 * t293 + t248 * t294) * t271 + (-t274 * t293 + t275 * t294) * t303;
t254 = -Icges(5,6) * t300 + t298 * t326;
t255 = Icges(5,6) * t298 + t300 * t326;
t256 = -Icges(5,5) * t300 + t298 * t329;
t257 = Icges(5,5) * t298 + t300 * t329;
t280 = Icges(5,2) * t299 + t346;
t281 = Icges(5,1) * t297 + t345;
t313 = (-t255 * t297 + t257 * t299) * t277 + (-t254 * t297 + t256 * t299) * t278 + (-t280 * t297 + t281 * t299) * t303;
t290 = rSges(2,1) * t310 - rSges(2,2) * t308;
t289 = rSges(2,1) * t308 + rSges(2,2) * t310;
t288 = rSges(4,1) * t307 + rSges(4,2) * t309;
t285 = Icges(4,5) * t307 + Icges(4,6) * t309;
t282 = rSges(5,1) * t297 + rSges(5,2) * t299;
t276 = rSges(6,1) * t293 + rSges(6,2) * t294;
t269 = t296 + t303 * (rSges(3,1) * t300 - rSges(3,2) * t298);
t268 = -t334 - t303 * (rSges(3,1) * t298 + rSges(3,2) * t300);
t267 = rSges(4,3) * t298 + t300 * t333;
t266 = -rSges(4,3) * t300 + t298 * t333;
t261 = Icges(4,3) * t298 + t300 * t324;
t260 = -Icges(4,3) * t300 + t298 * t324;
t259 = rSges(5,3) * t298 + t300 * t332;
t258 = -rSges(5,3) * t300 + t298 * t332;
t251 = rSges(6,3) * t298 + t300 * t331;
t250 = -rSges(6,3) * t300 + t298 * t331;
t238 = pkin(9) * t298 + t339 * t300;
t237 = -pkin(9) * t300 + t339 * t298;
t236 = (t266 * t298 + t267 * t300) * qJD(3);
t235 = t267 * t303 - t288 * t292 + t340;
t234 = -t334 - t288 * t337 + (-t266 - t283) * t303;
t233 = t259 * t303 - t277 * t282 + t318;
t232 = t278 * t282 + (-t258 + t341) * t303 + t315;
t231 = t258 * t277 - t259 * t278 + t342;
t230 = -t277 * t352 - t270 * t276 + (t238 + t251) * t303 + t318;
t229 = t278 * t352 + t271 * t276 + (-t237 - t250 + t341) * t303 + t315;
t228 = t237 * t277 - t238 * t278 + t250 * t270 - t251 * t271 + t342;
t1 = m(3) * (t268 ^ 2 + t269 ^ 2) / 0.2e1 + t303 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t234 ^ 2 + t235 ^ 2 + t236 ^ 2) / 0.2e1 + ((t285 * t298 + t300 * t319) * t303 + (t298 ^ 2 * t261 + (t321 * t300 + (-t260 + t320) * t298) * t300) * qJD(3)) * t292 / 0.2e1 - ((-t285 * t300 + t298 * t319) * t303 + (t260 * t300 ^ 2 + (t320 * t298 + (-t261 + t321) * t300) * t298) * qJD(3)) * t337 / 0.2e1 + m(5) * (t231 ^ 2 + t232 ^ 2 + t233 ^ 2) / 0.2e1 + t277 * (t298 * t316 + t300 * t313) / 0.2e1 + t278 * (t298 * t313 - t300 * t316) / 0.2e1 + m(6) * (t228 ^ 2 + t229 ^ 2 + t230 ^ 2) / 0.2e1 + t270 * (t298 * t317 + t300 * t314) / 0.2e1 + t271 * (t298 * t314 - t300 * t317) / 0.2e1 + (m(2) * (t289 ^ 2 + t290 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t263 * t309 + t265 * t307) * t298 - (t262 * t309 + t264 * t307) * t300) * qJD(3) + (t255 * t299 + t257 * t297) * t277 + (t254 * t299 + t256 * t297) * t278 + (t247 * t294 + t249 * t293) * t270 + (t246 * t294 + t248 * t293) * t271 + (t274 * t294 + t275 * t293 + t280 * t299 + t281 * t297 + t286 * t309 + t287 * t307) * t303) * t303 / 0.2e1;
T = t1;
