% Calculate kinetic energy for
% S5RRPRR5
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:29
% EndTime: 2020-01-03 12:03:30
% DurationCPUTime: 0.97s
% Computational Cost: add. (809->133), mult. (557->224), div. (0->0), fcn. (452->10), ass. (0->84)
t294 = qJD(1) + qJD(2);
t296 = sin(pkin(9));
t297 = cos(pkin(9));
t339 = -qJD(3) + t294 * (rSges(4,1) * t297 - rSges(4,2) * t296);
t295 = qJ(1) + qJ(2);
t288 = sin(t295);
t289 = cos(t295);
t293 = pkin(9) + qJ(4);
t287 = qJ(5) + t293;
t280 = sin(t287);
t281 = cos(t287);
t329 = Icges(6,4) * t281;
t312 = -Icges(6,2) * t280 + t329;
t249 = -Icges(6,6) * t289 + t312 * t288;
t250 = -Icges(6,6) * t288 - t312 * t289;
t330 = Icges(6,4) * t280;
t314 = Icges(6,1) * t281 - t330;
t251 = -Icges(6,5) * t289 + t314 * t288;
t252 = -Icges(6,5) * t288 - t314 * t289;
t267 = Icges(6,2) * t281 + t330;
t268 = Icges(6,1) * t280 + t329;
t321 = -qJD(4) - qJD(5);
t274 = t321 * t288;
t275 = t321 * t289;
t338 = (t267 * t280 - t268 * t281) * t294 + (t249 * t280 - t251 * t281) * t275 + (t250 * t280 - t252 * t281) * t274;
t334 = t297 * pkin(3);
t333 = pkin(1) * qJD(1);
t285 = sin(t293);
t332 = Icges(5,4) * t285;
t286 = cos(t293);
t331 = Icges(5,4) * t286;
t276 = -t289 * pkin(2) - t288 * qJ(3);
t327 = pkin(7) * t288 + t334 * t289 - t276;
t299 = sin(qJ(1));
t283 = t299 * t333;
t326 = t294 * (t288 * pkin(2) - t289 * qJ(3)) + t283;
t325 = pkin(4) * t286;
t323 = qJD(4) * t288;
t322 = qJD(4) * t289;
t320 = pkin(4) * qJD(4) * t285;
t300 = cos(qJ(1));
t284 = t300 * t333;
t319 = -qJD(3) * t289 + t284;
t317 = rSges(5,1) * t286 - rSges(5,2) * t285;
t316 = rSges(6,1) * t281 - rSges(6,2) * t280;
t315 = Icges(5,1) * t286 - t332;
t313 = -Icges(5,2) * t285 + t331;
t311 = Icges(5,5) * t286 - Icges(5,6) * t285;
t310 = Icges(6,5) * t281 - Icges(6,6) * t280;
t257 = -Icges(5,6) * t289 + t313 * t288;
t259 = -Icges(5,5) * t289 + t315 * t288;
t307 = -t257 * t285 + t259 * t286;
t258 = -Icges(5,6) * t288 - t313 * t289;
t260 = -Icges(5,5) * t288 - t315 * t289;
t306 = t258 * t285 - t260 * t286;
t271 = Icges(5,2) * t286 + t332;
t272 = Icges(5,1) * t285 + t331;
t304 = t271 * t285 - t272 * t286;
t303 = -qJD(3) * t288 + t294 * (-pkin(7) * t289 + t334 * t288) + t326;
t302 = -(-Icges(6,3) * t289 + t310 * t288) * t275 - (-Icges(6,3) * t288 - t310 * t289) * t274 - (Icges(6,5) * t280 + Icges(6,6) * t281) * t294;
t279 = -t300 * rSges(2,1) + t299 * rSges(2,2);
t278 = t299 * rSges(2,1) + t300 * rSges(2,2);
t273 = t285 * rSges(5,1) + t286 * rSges(5,2);
t270 = Icges(5,5) * t285 + Icges(5,6) * t286;
t269 = t280 * rSges(6,1) + t281 * rSges(6,2);
t264 = t284 - t294 * (-t289 * rSges(3,1) + t288 * rSges(3,2));
t263 = t283 + t294 * (t288 * rSges(3,1) + t289 * rSges(3,2));
t262 = -t288 * rSges(5,3) - t317 * t289;
t261 = -t289 * rSges(5,3) + t317 * t288;
t256 = -Icges(5,3) * t288 - t311 * t289;
t255 = -Icges(5,3) * t289 + t311 * t288;
t254 = -t288 * rSges(6,3) - t316 * t289;
t253 = -t289 * rSges(6,3) + t316 * t288;
t244 = -pkin(8) * t288 - t325 * t289;
t243 = -pkin(8) * t289 + t325 * t288;
t242 = t284 + (t288 * rSges(4,3) - t276) * t294 + t339 * t289;
t241 = -t294 * t289 * rSges(4,3) + t339 * t288 + t326;
t240 = (t261 * t288 - t262 * t289) * qJD(4);
t239 = -t273 * t323 + (-t262 + t327) * t294 + t319;
t238 = t294 * t261 + t273 * t322 + t303;
t237 = -t288 * t320 + t274 * t269 + (-t244 - t254 + t327) * t294 + t319;
t236 = t289 * t320 - t275 * t269 + (t243 + t253) * t294 + t303;
t235 = -t274 * t253 + t275 * t254 + (t243 * t288 - t244 * t289) * qJD(4);
t1 = m(3) * (t263 ^ 2 + t264 ^ 2) / 0.2e1 + m(4) * (t241 ^ 2 + t242 ^ 2) / 0.2e1 + m(5) * (t238 ^ 2 + t239 ^ 2 + t240 ^ 2) / 0.2e1 - ((-t289 * t270 - t304 * t288) * t294 + (t289 ^ 2 * t255 + (t306 * t288 + (t256 - t307) * t289) * t288) * qJD(4)) * t322 / 0.2e1 - ((-t288 * t270 + t304 * t289) * t294 + (t288 ^ 2 * t256 + (t307 * t289 + (t255 - t306) * t288) * t289) * qJD(4)) * t323 / 0.2e1 + m(6) * (t235 ^ 2 + t236 ^ 2 + t237 ^ 2) / 0.2e1 + t275 * (-t338 * t288 + t302 * t289) / 0.2e1 + t274 * (t302 * t288 + t338 * t289) / 0.2e1 + (Icges(3,3) + Icges(4,2) * t297 ^ 2 + (Icges(4,1) * t296 + 0.2e1 * Icges(4,4) * t297) * t296) * t294 ^ 2 / 0.2e1 + ((-(t286 * t257 + t285 * t259) * t289 - (t286 * t258 + t285 * t260) * t288) * qJD(4) + (t281 * t249 + t280 * t251) * t275 + (t281 * t250 + t280 * t252) * t274 + (t281 * t267 + t280 * t268 + t286 * t271 + t285 * t272) * t294) * t294 / 0.2e1 + (m(2) * (t278 ^ 2 + t279 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
