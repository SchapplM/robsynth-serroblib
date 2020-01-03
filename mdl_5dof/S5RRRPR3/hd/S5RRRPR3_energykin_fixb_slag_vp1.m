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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:08:49
% EndTime: 2020-01-03 12:08:50
% DurationCPUTime: 1.07s
% Computational Cost: add. (1014->166), mult. (764->258), div. (0->0), fcn. (634->10), ass. (0->103)
t376 = -Icges(4,3) - Icges(5,3);
t312 = qJ(3) + pkin(9);
t304 = sin(t312);
t305 = cos(t312);
t315 = sin(qJ(3));
t317 = cos(qJ(3));
t375 = Icges(4,5) * t317 + Icges(5,5) * t305 - Icges(4,6) * t315 - Icges(5,6) * t304;
t313 = qJ(1) + qJ(2);
t307 = sin(t313);
t308 = cos(t313);
t374 = t375 * t307 + t376 * t308;
t373 = t376 * t307 - t375 * t308;
t372 = -Icges(4,5) * t315 - Icges(5,5) * t304 - Icges(4,6) * t317 - Icges(5,6) * t305;
t356 = Icges(5,4) * t304;
t285 = Icges(5,2) * t305 + t356;
t355 = Icges(5,4) * t305;
t286 = Icges(5,1) * t304 + t355;
t358 = Icges(4,4) * t315;
t294 = Icges(4,2) * t317 + t358;
t357 = Icges(4,4) * t317;
t295 = Icges(4,1) * t315 + t357;
t371 = t285 * t304 - t286 * t305 + t294 * t315 - t295 * t317;
t335 = -Icges(5,2) * t304 + t355;
t264 = -Icges(5,6) * t307 - t308 * t335;
t338 = Icges(5,1) * t305 - t356;
t266 = -Icges(5,5) * t307 - t308 * t338;
t336 = -Icges(4,2) * t315 + t357;
t272 = -Icges(4,6) * t307 - t308 * t336;
t339 = Icges(4,1) * t317 - t358;
t274 = -Icges(4,5) * t307 - t308 * t339;
t370 = t264 * t304 - t266 * t305 + t272 * t315 - t274 * t317;
t263 = -Icges(5,6) * t308 + t307 * t335;
t265 = -Icges(5,5) * t308 + t307 * t338;
t271 = -Icges(4,6) * t308 + t307 * t336;
t273 = -Icges(4,5) * t308 + t307 * t339;
t369 = -t263 * t304 + t265 * t305 - t271 * t315 + t273 * t317;
t306 = qJ(5) + t312;
t299 = sin(t306);
t300 = cos(t306);
t353 = Icges(6,4) * t300;
t334 = -Icges(6,2) * t299 + t353;
t253 = -Icges(6,6) * t308 + t307 * t334;
t254 = -Icges(6,6) * t307 - t308 * t334;
t354 = Icges(6,4) * t299;
t337 = Icges(6,1) * t300 - t354;
t255 = -Icges(6,5) * t308 + t307 * t337;
t256 = -Icges(6,5) * t307 - t308 * t337;
t281 = Icges(6,2) * t300 + t354;
t282 = Icges(6,1) * t299 + t353;
t346 = -qJD(3) - qJD(5);
t288 = t346 * t307;
t289 = t346 * t308;
t311 = qJD(1) + qJD(2);
t368 = (t281 * t299 - t282 * t300) * t311 + (t253 * t299 - t255 * t300) * t289 + (t254 * t299 - t256 * t300) * t288;
t363 = pkin(3) * t315;
t362 = pkin(4) * t304;
t361 = t317 * pkin(3);
t359 = pkin(1) * qJD(1);
t258 = -qJ(4) * t307 - t308 * t361;
t290 = -t308 * pkin(2) - t307 * pkin(7);
t352 = -t258 - t290;
t316 = sin(qJ(1));
t302 = t316 * t359;
t351 = t311 * (t307 * pkin(2) - t308 * pkin(7)) + t302;
t350 = pkin(4) * t305;
t348 = qJD(3) * t307;
t347 = qJD(3) * t308;
t318 = cos(qJ(1));
t303 = t318 * t359;
t343 = -qJD(4) * t308 + t303;
t342 = rSges(4,1) * t317 - rSges(4,2) * t315;
t341 = rSges(5,1) * t305 - rSges(5,2) * t304;
t340 = rSges(6,1) * t300 - rSges(6,2) * t299;
t331 = Icges(6,5) * t300 - Icges(6,6) * t299;
t257 = -qJ(4) * t308 + t307 * t361;
t321 = -qJD(4) * t307 + t311 * t257 + t347 * t363 + t351;
t320 = -(-Icges(6,3) * t308 + t307 * t331) * t289 - (-Icges(6,3) * t307 - t308 * t331) * t288 - (Icges(6,5) * t299 + Icges(6,6) * t300) * t311;
t298 = -t318 * rSges(2,1) + t316 * rSges(2,2);
t297 = t316 * rSges(2,1) + t318 * rSges(2,2);
t296 = t315 * rSges(4,1) + t317 * rSges(4,2);
t287 = t304 * rSges(5,1) + t305 * rSges(5,2);
t283 = t299 * rSges(6,1) + t300 * rSges(6,2);
t278 = t303 - t311 * (-t308 * rSges(3,1) + t307 * rSges(3,2));
t277 = t302 + t311 * (t307 * rSges(3,1) + t308 * rSges(3,2));
t276 = -t307 * rSges(4,3) - t308 * t342;
t275 = -t308 * rSges(4,3) + t307 * t342;
t268 = -t307 * rSges(5,3) - t308 * t341;
t267 = -t308 * rSges(5,3) + t307 * t341;
t260 = -t307 * rSges(6,3) - t308 * t340;
t259 = -t308 * rSges(6,3) + t307 * t340;
t249 = t257 * t348;
t248 = -pkin(8) * t307 - t308 * t350;
t247 = -pkin(8) * t308 + t307 * t350;
t246 = (t275 * t307 - t276 * t308) * qJD(3);
t245 = -t296 * t348 + t303 + (-t276 - t290) * t311;
t244 = t311 * t275 + t296 * t347 + t351;
t243 = (-t287 - t363) * t348 + (-t268 + t352) * t311 + t343;
t242 = t311 * t267 + t287 * t347 + t321;
t241 = t249 + (t267 * t307 + (-t258 - t268) * t308) * qJD(3);
t240 = t288 * t283 + (-t362 - t363) * t348 + (-t248 - t260 + t352) * t311 + t343;
t239 = t347 * t362 - t289 * t283 + (t247 + t259) * t311 + t321;
t238 = -t288 * t259 + t289 * t260 + t249 + (t247 * t307 + (-t248 - t258) * t308) * qJD(3);
t1 = m(3) * (t277 ^ 2 + t278 ^ 2) / 0.2e1 + t311 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t244 ^ 2 + t245 ^ 2 + t246 ^ 2) / 0.2e1 + m(5) * (t241 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + m(6) * (t238 ^ 2 + t239 ^ 2 + t240 ^ 2) / 0.2e1 + t289 * (-t368 * t307 + t320 * t308) / 0.2e1 + t288 * (t320 * t307 + t368 * t308) / 0.2e1 + (m(2) * (t297 ^ 2 + t298 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 - ((t372 * t307 + t371 * t308) * t311 + (t373 * t307 ^ 2 + (t369 * t308 + (-t370 + t374) * t307) * t308) * qJD(3)) * t348 / 0.2e1 - ((-t371 * t307 + t372 * t308) * t311 + (t374 * t308 ^ 2 + (t370 * t307 + (-t369 + t373) * t308) * t307) * qJD(3)) * t347 / 0.2e1 + ((t300 * t253 + t299 * t255) * t289 + (t300 * t254 + t299 * t256) * t288 + ((-t305 * t263 - t304 * t265 - t317 * t271 - t315 * t273) * t308 + (-t305 * t264 - t304 * t266 - t317 * t272 - t315 * t274) * t307) * qJD(3) + (t300 * t281 + t299 * t282 + t305 * t285 + t304 * t286 + t317 * t294 + t315 * t295) * t311) * t311 / 0.2e1;
T = t1;
