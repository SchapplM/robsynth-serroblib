% Calculate kinetic energy for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:25
% EndTime: 2020-01-03 11:27:26
% DurationCPUTime: 0.70s
% Computational Cost: add. (783->136), mult. (559->222), div. (0->0), fcn. (452->10), ass. (0->84)
t295 = qJ(1) + pkin(8);
t288 = sin(t295);
t290 = cos(t295);
t294 = pkin(9) + qJ(4);
t291 = qJ(5) + t294;
t282 = sin(t291);
t283 = cos(t291);
t330 = Icges(6,4) * t283;
t313 = -Icges(6,2) * t282 + t330;
t251 = -Icges(6,6) * t290 + t313 * t288;
t252 = -Icges(6,6) * t288 - t313 * t290;
t331 = Icges(6,4) * t282;
t315 = Icges(6,1) * t283 - t331;
t253 = -Icges(6,5) * t290 + t315 * t288;
t254 = -Icges(6,5) * t288 - t315 * t290;
t268 = Icges(6,2) * t283 + t331;
t269 = Icges(6,1) * t282 + t330;
t322 = -qJD(4) - qJD(5);
t272 = t322 * t288;
t273 = t322 * t290;
t338 = (t268 * t282 - t269 * t283) * qJD(1) + (t251 * t282 - t253 * t283) * t273 + (t252 * t282 - t254 * t283) * t272;
t297 = cos(pkin(9));
t335 = t297 * pkin(3);
t334 = pkin(1) * qJD(1);
t287 = sin(t294);
t333 = Icges(5,4) * t287;
t289 = cos(t294);
t332 = Icges(5,4) * t289;
t278 = -t290 * pkin(2) - t288 * qJ(3);
t328 = pkin(6) * t288 + t335 * t290 - t278;
t299 = sin(qJ(1));
t285 = t299 * t334;
t327 = qJD(1) * (t288 * pkin(2) - t290 * qJ(3)) + t285;
t326 = pkin(4) * t289;
t324 = qJD(4) * t288;
t323 = qJD(4) * t290;
t321 = pkin(4) * qJD(4) * t287;
t300 = cos(qJ(1));
t286 = t300 * t334;
t320 = -qJD(3) * t290 + t286;
t296 = sin(pkin(9));
t319 = rSges(4,1) * t297 - rSges(4,2) * t296;
t318 = rSges(5,1) * t289 - rSges(5,2) * t287;
t317 = rSges(6,1) * t283 - rSges(6,2) * t282;
t316 = Icges(5,1) * t289 - t333;
t314 = -Icges(5,2) * t287 + t332;
t312 = Icges(5,5) * t289 - Icges(5,6) * t287;
t311 = Icges(6,5) * t283 - Icges(6,6) * t282;
t259 = -Icges(5,6) * t290 + t314 * t288;
t261 = -Icges(5,5) * t290 + t316 * t288;
t308 = -t259 * t287 + t261 * t289;
t260 = -Icges(5,6) * t288 - t314 * t290;
t262 = -Icges(5,5) * t288 - t316 * t290;
t307 = t260 * t287 - t262 * t289;
t275 = Icges(5,2) * t289 + t333;
t276 = Icges(5,1) * t287 + t332;
t305 = t275 * t287 - t276 * t289;
t304 = -qJD(3) * t288 + qJD(1) * (-pkin(6) * t290 + t335 * t288) + t327;
t303 = -(Icges(6,5) * t282 + Icges(6,6) * t283) * qJD(1) - (-Icges(6,3) * t290 + t311 * t288) * t273 - (-Icges(6,3) * t288 - t311 * t290) * t272;
t301 = qJD(2) ^ 2;
t281 = -t300 * rSges(2,1) + t299 * rSges(2,2);
t280 = t299 * rSges(2,1) + t300 * rSges(2,2);
t277 = t287 * rSges(5,1) + t289 * rSges(5,2);
t274 = Icges(5,5) * t287 + Icges(5,6) * t289;
t270 = t282 * rSges(6,1) + t283 * rSges(6,2);
t266 = t286 - qJD(1) * (-t290 * rSges(3,1) + t288 * rSges(3,2));
t265 = t285 + qJD(1) * (t288 * rSges(3,1) + t290 * rSges(3,2));
t264 = -t288 * rSges(5,3) - t318 * t290;
t263 = -t290 * rSges(5,3) + t318 * t288;
t258 = -Icges(5,3) * t288 - t312 * t290;
t257 = -Icges(5,3) * t290 + t312 * t288;
t256 = -t288 * rSges(6,3) - t317 * t290;
t255 = -t290 * rSges(6,3) + t317 * t288;
t246 = -pkin(7) * t288 - t326 * t290;
t245 = -pkin(7) * t290 + t326 * t288;
t244 = (t288 * rSges(4,3) + t319 * t290 - t278) * qJD(1) + t320;
t243 = -qJD(1) * t290 * rSges(4,3) + (qJD(1) * t319 - qJD(3)) * t288 + t327;
t242 = qJD(2) + (t263 * t288 - t264 * t290) * qJD(4);
t241 = -t277 * t324 + (-t264 + t328) * qJD(1) + t320;
t240 = qJD(1) * t263 + t277 * t323 + t304;
t239 = -t288 * t321 + t272 * t270 + (-t246 - t256 + t328) * qJD(1) + t320;
t238 = t290 * t321 - t273 * t270 + (t245 + t255) * qJD(1) + t304;
t237 = -t272 * t255 + t273 * t256 + qJD(2) + (t245 * t288 - t246 * t290) * qJD(4);
t1 = m(3) * (t265 ^ 2 + t266 ^ 2 + t301) / 0.2e1 + m(4) * (t243 ^ 2 + t244 ^ 2 + t301) / 0.2e1 + m(5) * (t240 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 - ((-t290 * t274 - t305 * t288) * qJD(1) + (t290 ^ 2 * t257 + (t307 * t288 + (t258 - t308) * t290) * t288) * qJD(4)) * t323 / 0.2e1 - ((-t288 * t274 + t305 * t290) * qJD(1) + (t288 ^ 2 * t258 + (t308 * t290 + (t257 - t307) * t288) * t290) * qJD(4)) * t324 / 0.2e1 + m(6) * (t237 ^ 2 + t238 ^ 2 + t239 ^ 2) / 0.2e1 + t273 * (-t338 * t288 + t303 * t290) / 0.2e1 + t272 * (t303 * t288 + t338 * t290) / 0.2e1 + ((-(t289 * t259 + t287 * t261) * t290 - (t289 * t260 + t287 * t262) * t288) * qJD(4) + (t283 * t251 + t282 * t253) * t273 + (t283 * t252 + t282 * t254) * t272 + (t283 * t268 + t282 * t269 + t289 * t275 + t287 * t276) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t280 ^ 2 + t281 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t297 ^ 2 + (Icges(4,1) * t296 + 0.2e1 * Icges(4,4) * t297) * t296) * qJD(1) ^ 2 / 0.2e1;
T = t1;
