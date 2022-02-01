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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:02:10
% EndTime: 2022-01-20 11:02:10
% DurationCPUTime: 0.69s
% Computational Cost: add. (809->133), mult. (557->224), div. (0->0), fcn. (452->10), ass. (0->84)
t289 = cos(pkin(9));
t325 = t289 * pkin(3);
t324 = pkin(1) * qJD(1);
t285 = pkin(9) + qJ(4);
t277 = sin(t285);
t323 = Icges(5,4) * t277;
t278 = cos(t285);
t322 = Icges(5,4) * t278;
t279 = qJ(5) + t285;
t273 = sin(t279);
t321 = Icges(6,4) * t273;
t274 = cos(t279);
t320 = Icges(6,4) * t274;
t287 = qJ(1) + qJ(2);
t280 = sin(t287);
t281 = cos(t287);
t268 = pkin(2) * t280 - qJ(3) * t281;
t318 = pkin(7) * t281 - t325 * t280 - t268;
t292 = cos(qJ(1));
t276 = t292 * t324;
t286 = qJD(1) + qJD(2);
t317 = t286 * (pkin(2) * t281 + qJ(3) * t280) + t276;
t316 = pkin(4) * t278;
t314 = qJD(4) * t280;
t313 = qJD(4) * t281;
t312 = qJD(4) + qJD(5);
t311 = pkin(4) * qJD(4) * t277;
t291 = sin(qJ(1));
t310 = t291 * t324;
t309 = qJD(3) * t280 - t310;
t288 = sin(pkin(9));
t308 = rSges(4,1) * t289 - rSges(4,2) * t288;
t307 = rSges(5,1) * t278 - rSges(5,2) * t277;
t306 = rSges(6,1) * t274 - rSges(6,2) * t273;
t305 = Icges(5,1) * t278 - t323;
t304 = Icges(6,1) * t274 - t321;
t303 = -Icges(5,2) * t277 + t322;
t302 = -Icges(6,2) * t273 + t320;
t301 = Icges(5,5) * t278 - Icges(5,6) * t277;
t300 = Icges(6,5) * t274 - Icges(6,6) * t273;
t249 = -Icges(5,6) * t281 + t303 * t280;
t251 = -Icges(5,5) * t281 + t305 * t280;
t299 = t249 * t277 - t251 * t278;
t250 = Icges(5,6) * t280 + t303 * t281;
t252 = Icges(5,5) * t280 + t305 * t281;
t298 = -t250 * t277 + t252 * t278;
t263 = Icges(5,2) * t278 + t323;
t264 = Icges(5,1) * t277 + t322;
t297 = -t263 * t277 + t264 * t278;
t296 = -qJD(3) * t281 + t286 * (pkin(7) * t280 + t325 * t281) + t317;
t266 = t312 * t280;
t267 = t312 * t281;
t295 = -(-Icges(6,3) * t281 + t300 * t280) * t267 + (Icges(6,3) * t280 + t300 * t281) * t266 + (Icges(6,5) * t273 + Icges(6,6) * t274) * t286;
t241 = -Icges(6,6) * t281 + t302 * t280;
t242 = Icges(6,6) * t280 + t302 * t281;
t243 = -Icges(6,5) * t281 + t304 * t280;
t244 = Icges(6,5) * t280 + t304 * t281;
t259 = Icges(6,2) * t274 + t321;
t260 = Icges(6,1) * t273 + t320;
t294 = (-t242 * t273 + t244 * t274) * t266 - (-t241 * t273 + t243 * t274) * t267 + (-t259 * t273 + t260 * t274) * t286;
t271 = rSges(2,1) * t292 - rSges(2,2) * t291;
t270 = rSges(2,1) * t291 + rSges(2,2) * t292;
t265 = rSges(5,1) * t277 + rSges(5,2) * t278;
t262 = Icges(5,5) * t277 + Icges(5,6) * t278;
t261 = rSges(6,1) * t273 + rSges(6,2) * t274;
t256 = t276 + t286 * (rSges(3,1) * t281 - rSges(3,2) * t280);
t255 = -t310 - t286 * (rSges(3,1) * t280 + rSges(3,2) * t281);
t254 = rSges(5,3) * t280 + t307 * t281;
t253 = -rSges(5,3) * t281 + t307 * t280;
t248 = Icges(5,3) * t280 + t301 * t281;
t247 = -Icges(5,3) * t281 + t301 * t280;
t246 = rSges(6,3) * t280 + t306 * t281;
t245 = -rSges(6,3) * t281 + t306 * t280;
t236 = pkin(8) * t280 + t316 * t281;
t235 = -pkin(8) * t281 + t316 * t280;
t234 = t286 * t280 * rSges(4,3) + (t286 * t308 - qJD(3)) * t281 + t317;
t233 = (t281 * rSges(4,3) - t308 * t280 - t268) * t286 + t309;
t232 = (t253 * t280 + t254 * t281) * qJD(4);
t231 = t254 * t286 - t265 * t314 + t296;
t230 = -t265 * t313 + (-t253 + t318) * t286 + t309;
t229 = -t280 * t311 - t261 * t266 + (t236 + t246) * t286 + t296;
t228 = -t281 * t311 - t261 * t267 + (-t235 - t245 + t318) * t286 + t309;
t227 = t245 * t266 + t246 * t267 + (t235 * t280 + t236 * t281) * qJD(4);
t1 = m(3) * (t255 ^ 2 + t256 ^ 2) / 0.2e1 + m(4) * (t233 ^ 2 + t234 ^ 2) / 0.2e1 + m(5) * (t230 ^ 2 + t231 ^ 2 + t232 ^ 2) / 0.2e1 + ((t280 * t262 + t297 * t281) * t286 + (t280 ^ 2 * t248 + (t299 * t281 + (-t247 + t298) * t280) * t281) * qJD(4)) * t314 / 0.2e1 - ((-t281 * t262 + t297 * t280) * t286 + (t281 ^ 2 * t247 + (t298 * t280 + (-t248 + t299) * t281) * t280) * qJD(4)) * t313 / 0.2e1 + m(6) * (t227 ^ 2 + t228 ^ 2 + t229 ^ 2) / 0.2e1 + t266 * (t295 * t280 + t294 * t281) / 0.2e1 - t267 * (t294 * t280 - t295 * t281) / 0.2e1 + (Icges(3,3) + Icges(4,2) * t289 ^ 2 + (Icges(4,1) * t288 + 0.2e1 * Icges(4,4) * t289) * t288) * t286 ^ 2 / 0.2e1 + (((t250 * t278 + t252 * t277) * t280 - (t249 * t278 + t251 * t277) * t281) * qJD(4) + (t242 * t274 + t244 * t273) * t266 - (t241 * t274 + t243 * t273) * t267 + (t274 * t259 + t273 * t260 + t278 * t263 + t277 * t264) * t286) * t286 / 0.2e1 + (m(2) * (t270 ^ 2 + t271 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
