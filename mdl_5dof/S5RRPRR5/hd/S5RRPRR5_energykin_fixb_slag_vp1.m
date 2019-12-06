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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:33:36
% EndTime: 2019-12-05 18:33:36
% DurationCPUTime: 0.61s
% Computational Cost: add. (809->132), mult. (557->223), div. (0->0), fcn. (452->10), ass. (0->84)
t297 = qJ(1) + qJ(2);
t290 = sin(t297);
t291 = cos(t297);
t295 = pkin(9) + qJ(4);
t289 = qJ(5) + t295;
t284 = sin(t289);
t285 = cos(t289);
t333 = Icges(6,4) * t285;
t315 = -Icges(6,2) * t284 + t333;
t251 = Icges(6,6) * t291 - t315 * t290;
t252 = Icges(6,6) * t290 + t315 * t291;
t334 = Icges(6,4) * t284;
t317 = Icges(6,1) * t285 - t334;
t253 = Icges(6,5) * t291 - t317 * t290;
t254 = Icges(6,5) * t290 + t317 * t291;
t269 = Icges(6,2) * t285 + t334;
t270 = Icges(6,1) * t284 + t333;
t326 = qJD(4) + qJD(5);
t276 = t326 * t290;
t277 = t326 * t291;
t296 = qJD(1) + qJD(2);
t342 = (t269 * t284 - t270 * t285) * t296 + (t251 * t284 - t253 * t285) * t277 + (t252 * t284 - t254 * t285) * t276;
t299 = cos(pkin(9));
t338 = t299 * pkin(3);
t337 = pkin(1) * qJD(1);
t287 = sin(t295);
t336 = Icges(5,4) * t287;
t288 = cos(t295);
t335 = Icges(5,4) * t288;
t278 = pkin(2) * t291 + qJ(3) * t290;
t331 = -pkin(7) * t290 - t338 * t291 - t278;
t330 = pkin(4) * t288;
t328 = qJD(4) * t290;
t327 = qJD(4) * t291;
t325 = pkin(4) * qJD(4) * t287;
t301 = sin(qJ(1));
t324 = t301 * t337;
t302 = cos(qJ(1));
t323 = t302 * t337;
t322 = qJD(3) * t291 - t323;
t298 = sin(pkin(9));
t321 = -rSges(4,1) * t299 + rSges(4,2) * t298;
t320 = rSges(5,1) * t288 - rSges(5,2) * t287;
t319 = rSges(6,1) * t285 - rSges(6,2) * t284;
t318 = Icges(5,1) * t288 - t336;
t316 = -Icges(5,2) * t287 + t335;
t314 = Icges(5,5) * t288 - Icges(5,6) * t287;
t313 = Icges(6,5) * t285 - Icges(6,6) * t284;
t259 = Icges(5,6) * t291 - t316 * t290;
t261 = Icges(5,5) * t291 - t318 * t290;
t310 = -t259 * t287 + t261 * t288;
t260 = Icges(5,6) * t290 + t316 * t291;
t262 = Icges(5,5) * t290 + t318 * t291;
t309 = t260 * t287 - t262 * t288;
t273 = Icges(5,2) * t288 + t336;
t274 = Icges(5,1) * t287 + t335;
t307 = t273 * t287 - t274 * t288;
t306 = t296 * (-pkin(2) * t290 + qJ(3) * t291) + qJD(3) * t290 - t324;
t305 = t296 * (pkin(7) * t291 - t338 * t290) + t306;
t304 = (Icges(6,3) * t291 - t313 * t290) * t277 + (Icges(6,3) * t290 + t313 * t291) * t276 + (Icges(6,5) * t284 + Icges(6,6) * t285) * t296;
t281 = rSges(2,1) * t302 - rSges(2,2) * t301;
t280 = -rSges(2,1) * t301 - rSges(2,2) * t302;
t275 = rSges(5,1) * t287 + rSges(5,2) * t288;
t272 = Icges(5,5) * t287 + Icges(5,6) * t288;
t271 = rSges(6,1) * t284 + rSges(6,2) * t285;
t266 = -t323 - t296 * (rSges(3,1) * t291 - rSges(3,2) * t290);
t265 = -t324 + t296 * (-rSges(3,1) * t290 - rSges(3,2) * t291);
t264 = rSges(5,3) * t290 + t320 * t291;
t263 = rSges(5,3) * t291 - t320 * t290;
t258 = Icges(5,3) * t290 + t314 * t291;
t257 = Icges(5,3) * t291 - t314 * t290;
t256 = rSges(6,3) * t290 + t319 * t291;
t255 = rSges(6,3) * t291 - t319 * t290;
t246 = pkin(8) * t290 + t330 * t291;
t245 = pkin(8) * t291 - t330 * t290;
t244 = (-t290 * rSges(4,3) + t321 * t291 - t278) * t296 + t322;
t243 = t296 * (rSges(4,3) * t291 + t321 * t290) + t306;
t242 = (-t263 * t290 + t264 * t291) * qJD(4);
t241 = t275 * t328 + (-t264 + t331) * t296 + t322;
t240 = t263 * t296 - t275 * t327 + t305;
t239 = t290 * t325 + t271 * t276 + (-t246 - t256 + t331) * t296 + t322;
t238 = -t291 * t325 - t271 * t277 + (t245 + t255) * t296 + t305;
t237 = -t255 * t276 + t256 * t277 + (-t245 * t290 + t246 * t291) * qJD(4);
t1 = m(3) * (t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(4) * (t243 ^ 2 + t244 ^ 2) / 0.2e1 + m(5) * (t240 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + ((t291 * t272 + t307 * t290) * t296 + (t291 ^ 2 * t257 + (t309 * t290 + (t258 - t310) * t291) * t290) * qJD(4)) * t327 / 0.2e1 + ((t290 * t272 - t307 * t291) * t296 + (t290 ^ 2 * t258 + (t310 * t291 + (t257 - t309) * t290) * t291) * qJD(4)) * t328 / 0.2e1 + m(6) * (t237 ^ 2 + t238 ^ 2 + t239 ^ 2) / 0.2e1 + t277 * (t342 * t290 + t304 * t291) / 0.2e1 + t276 * (t304 * t290 - t342 * t291) / 0.2e1 + (Icges(3,3) + Icges(4,2) * t299 ^ 2 + (Icges(4,1) * t298 + 0.2e1 * Icges(4,4) * t299) * t298) * t296 ^ 2 / 0.2e1 + (((t259 * t288 + t261 * t287) * t291 + (t260 * t288 + t262 * t287) * t290) * qJD(4) + (t251 * t285 + t253 * t284) * t277 + (t252 * t285 + t254 * t284) * t276 + (t285 * t269 + t284 * t270 + t288 * t273 + t287 * t274) * t296) * t296 / 0.2e1 + (m(2) * (t280 ^ 2 + t281 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
