% Calculate kinetic energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:48
% EndTime: 2022-01-20 11:16:48
% DurationCPUTime: 0.84s
% Computational Cost: add. (981->173), mult. (1025->296), div. (0->0), fcn. (1032->10), ass. (0->89)
t314 = cos(qJ(4));
t343 = pkin(4) * t314;
t341 = pkin(1) * qJD(1);
t309 = qJ(1) + qJ(2);
t303 = sin(t309);
t310 = sin(pkin(9));
t340 = t303 * t310;
t311 = cos(pkin(9));
t339 = t303 * t311;
t312 = sin(qJ(4));
t338 = t303 * t312;
t305 = cos(t309);
t337 = t305 * t310;
t336 = t305 * t311;
t335 = t305 * t312;
t334 = t311 * t312;
t333 = t311 * t314;
t315 = cos(qJ(1));
t301 = t315 * t341;
t307 = qJD(1) + qJD(2);
t332 = t307 * (pkin(2) * t305 + qJ(3) * t303) + t301;
t331 = qJD(4) * t310;
t330 = qJD(4) + qJD(5);
t313 = sin(qJ(1));
t329 = t313 * t341;
t323 = pkin(3) * t311 + pkin(7) * t310;
t328 = t307 * t323 * t305 + t332;
t327 = t303 * t331;
t325 = t330 * t310;
t324 = qJD(3) * t303 - t329;
t322 = rSges(4,1) * t311 - rSges(4,2) * t310;
t288 = -t303 * t334 - t305 * t314;
t289 = t303 * t333 - t335;
t290 = t303 * t314 - t305 * t334;
t291 = t305 * t333 + t338;
t321 = (Icges(5,5) * t289 + Icges(5,6) * t288 + Icges(5,3) * t340) * t303 + (Icges(5,5) * t291 + Icges(5,6) * t290 + Icges(5,3) * t337) * t305;
t320 = t321 * t310;
t319 = pkin(8) * t310 + t343 * t311;
t295 = pkin(2) * t303 - qJ(3) * t305;
t318 = (-t323 * t303 - t295) * t307 + t324;
t308 = qJ(4) + qJ(5);
t304 = cos(t308);
t302 = sin(t308);
t298 = -qJD(4) * t311 + t307;
t297 = rSges(2,1) * t315 - rSges(2,2) * t313;
t296 = rSges(2,1) * t313 + rSges(2,2) * t315;
t294 = -t330 * t311 + t307;
t287 = t305 * t325;
t286 = t303 * t325;
t285 = -t311 * rSges(5,3) + (rSges(5,1) * t314 - rSges(5,2) * t312) * t310;
t284 = -Icges(5,5) * t311 + (Icges(5,1) * t314 - Icges(5,4) * t312) * t310;
t283 = -Icges(5,6) * t311 + (Icges(5,4) * t314 - Icges(5,2) * t312) * t310;
t282 = -Icges(5,3) * t311 + (Icges(5,5) * t314 - Icges(5,6) * t312) * t310;
t280 = t302 * t303 + t304 * t336;
t279 = -t302 * t336 + t303 * t304;
t278 = -t302 * t305 + t304 * t339;
t277 = -t302 * t339 - t304 * t305;
t276 = t301 + t307 * (rSges(3,1) * t305 - rSges(3,2) * t303);
t275 = -t329 - t307 * (rSges(3,1) * t303 + rSges(3,2) * t305);
t274 = -t311 * rSges(6,3) + (rSges(6,1) * t304 - rSges(6,2) * t302) * t310;
t273 = -Icges(6,5) * t311 + (Icges(6,1) * t304 - Icges(6,4) * t302) * t310;
t272 = -Icges(6,6) * t311 + (Icges(6,4) * t304 - Icges(6,2) * t302) * t310;
t271 = -Icges(6,3) * t311 + (Icges(6,5) * t304 - Icges(6,6) * t302) * t310;
t270 = -pkin(8) * t311 + t343 * t310;
t269 = rSges(5,1) * t291 + rSges(5,2) * t290 + rSges(5,3) * t337;
t268 = rSges(5,1) * t289 + rSges(5,2) * t288 + rSges(5,3) * t340;
t267 = pkin(4) * t338 + t319 * t305;
t266 = -pkin(4) * t335 + t319 * t303;
t265 = Icges(5,1) * t291 + Icges(5,4) * t290 + Icges(5,5) * t337;
t264 = Icges(5,1) * t289 + Icges(5,4) * t288 + Icges(5,5) * t340;
t263 = Icges(5,4) * t291 + Icges(5,2) * t290 + Icges(5,6) * t337;
t262 = Icges(5,4) * t289 + Icges(5,2) * t288 + Icges(5,6) * t340;
t259 = rSges(6,1) * t280 + rSges(6,2) * t279 + rSges(6,3) * t337;
t258 = rSges(6,1) * t278 + rSges(6,2) * t277 + rSges(6,3) * t340;
t257 = t307 * t303 * rSges(4,3) + (t307 * t322 - qJD(3)) * t305 + t332;
t256 = (t305 * rSges(4,3) - t322 * t303 - t295) * t307 + t324;
t255 = Icges(6,1) * t280 + Icges(6,4) * t279 + Icges(6,5) * t337;
t254 = Icges(6,1) * t278 + Icges(6,4) * t277 + Icges(6,5) * t340;
t253 = Icges(6,4) * t280 + Icges(6,2) * t279 + Icges(6,6) * t337;
t252 = Icges(6,4) * t278 + Icges(6,2) * t277 + Icges(6,6) * t340;
t251 = Icges(6,5) * t280 + Icges(6,6) * t279 + Icges(6,3) * t337;
t250 = Icges(6,5) * t278 + Icges(6,6) * t277 + Icges(6,3) * t340;
t249 = (t268 * t305 - t269 * t303) * t331;
t248 = t269 * t298 + (-t285 * t331 - qJD(3)) * t305 + t328;
t247 = -t268 * t298 + t285 * t327 + t318;
t246 = t259 * t294 + t267 * t298 - t274 * t287 + (-t270 * t331 - qJD(3)) * t305 + t328;
t245 = -t258 * t294 - t266 * t298 + t270 * t327 + t274 * t286 + t318;
t244 = t258 * t287 - t259 * t286 + (t266 * t305 - t267 * t303) * t331;
t1 = m(3) * (t275 ^ 2 + t276 ^ 2) / 0.2e1 + m(4) * (t256 ^ 2 + t257 ^ 2) / 0.2e1 + m(5) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + t298 * ((-t311 * t282 + (-t283 * t312 + t284 * t314) * t310) * t298 + (((-t263 * t312 + t265 * t314) * t305 + (-t262 * t312 + t264 * t314) * t303) * t310 - t321 * t311) * t331) / 0.2e1 + m(6) * (t244 ^ 2 + t245 ^ 2 + t246 ^ 2) / 0.2e1 + t287 * ((t251 * t337 + t279 * t253 + t280 * t255) * t287 + (t250 * t337 + t252 * t279 + t254 * t280) * t286 + (t271 * t337 + t272 * t279 + t273 * t280) * t294) / 0.2e1 + t286 * ((t251 * t340 + t253 * t277 + t255 * t278) * t287 + (t250 * t340 + t277 * t252 + t278 * t254) * t286 + (t271 * t340 + t272 * t277 + t273 * t278) * t294) / 0.2e1 + t294 * ((-t250 * t286 - t251 * t287 - t271 * t294) * t311 + ((-t253 * t302 + t255 * t304) * t287 + (-t252 * t302 + t254 * t304) * t286 + (-t272 * t302 + t273 * t304) * t294) * t310) / 0.2e1 + (Icges(3,3) + Icges(4,2) * t311 ^ 2 + (Icges(4,1) * t310 + 0.2e1 * Icges(4,4) * t311) * t310) * t307 ^ 2 / 0.2e1 + (m(2) * (t296 ^ 2 + t297 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (t305 * ((t282 * t337 + t283 * t290 + t284 * t291) * t298 + ((t262 * t290 + t264 * t291) * t303 + (t290 * t263 + t291 * t265 + t320) * t305) * t331) + t303 * ((t282 * t340 + t283 * t288 + t284 * t289) * t298 + ((t263 * t288 + t265 * t289) * t305 + (t288 * t262 + t289 * t264 + t320) * t303) * t331)) * t331 / 0.2e1;
T = t1;
