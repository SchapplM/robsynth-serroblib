% Calculate kinetic energy for
% S5RRPRR4
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:47
% EndTime: 2020-01-03 12:01:48
% DurationCPUTime: 0.60s
% Computational Cost: add. (803->127), mult. (514->211), div. (0->0), fcn. (408->10), ass. (0->82)
t267 = qJ(1) + qJ(2);
t259 = pkin(9) + t267;
t254 = sin(t259);
t255 = cos(t259);
t266 = qJ(4) + qJ(5);
t260 = sin(t266);
t262 = cos(t266);
t296 = Icges(6,4) * t262;
t283 = -Icges(6,2) * t260 + t296;
t222 = -Icges(6,6) * t255 + t254 * t283;
t223 = -Icges(6,6) * t254 - t255 * t283;
t297 = Icges(6,4) * t260;
t285 = Icges(6,1) * t262 - t297;
t224 = -Icges(6,5) * t255 + t254 * t285;
t225 = -Icges(6,5) * t254 - t255 * t285;
t291 = -qJD(4) - qJD(5);
t239 = t291 * t254;
t240 = t291 * t255;
t243 = Icges(6,2) * t262 + t297;
t244 = Icges(6,1) * t260 + t296;
t265 = qJD(1) + qJD(2);
t307 = (t243 * t260 - t244 * t262) * t265 + (t222 * t260 - t224 * t262) * t240 + (t223 * t260 - t225 * t262) * t239;
t303 = pkin(2) * t265;
t270 = cos(qJ(4));
t302 = pkin(4) * t270;
t300 = pkin(1) * qJD(1);
t268 = sin(qJ(4));
t299 = Icges(5,4) * t268;
t298 = Icges(5,4) * t270;
t269 = sin(qJ(1));
t257 = t269 * t300;
t261 = sin(t267);
t295 = t261 * t303 + t257;
t271 = cos(qJ(1));
t258 = t271 * t300;
t263 = cos(t267);
t294 = t263 * t303 + t258;
t293 = qJD(4) * t254;
t292 = qJD(4) * t255;
t290 = pkin(4) * qJD(4) * t268;
t289 = t265 * (pkin(3) * t254 - pkin(7) * t255) + t295;
t288 = rSges(5,1) * t270 - rSges(5,2) * t268;
t287 = rSges(6,1) * t262 - rSges(6,2) * t260;
t286 = Icges(5,1) * t270 - t299;
t284 = -Icges(5,2) * t268 + t298;
t282 = Icges(5,5) * t270 - Icges(5,6) * t268;
t281 = Icges(6,5) * t262 - Icges(6,6) * t260;
t230 = -Icges(5,6) * t255 + t254 * t284;
t232 = -Icges(5,5) * t255 + t254 * t286;
t278 = -t230 * t268 + t232 * t270;
t231 = -Icges(5,6) * t254 - t255 * t284;
t233 = -Icges(5,5) * t254 - t255 * t286;
t277 = t231 * t268 - t233 * t270;
t247 = Icges(5,2) * t270 + t299;
t248 = Icges(5,1) * t268 + t298;
t275 = t247 * t268 - t248 * t270;
t274 = -(-Icges(6,3) * t255 + t254 * t281) * t240 - (-Icges(6,3) * t254 - t255 * t281) * t239 - (Icges(6,5) * t260 + Icges(6,6) * t262) * t265;
t253 = -rSges(2,1) * t271 + rSges(2,2) * t269;
t252 = rSges(2,1) * t269 + rSges(2,2) * t271;
t251 = rSges(5,1) * t268 + rSges(5,2) * t270;
t246 = Icges(5,5) * t268 + Icges(5,6) * t270;
t245 = rSges(6,1) * t260 + rSges(6,2) * t262;
t241 = -pkin(3) * t255 - pkin(7) * t254;
t237 = t258 - t265 * (-rSges(3,1) * t263 + rSges(3,2) * t261);
t236 = t257 + t265 * (rSges(3,1) * t261 + rSges(3,2) * t263);
t235 = -rSges(5,3) * t254 - t255 * t288;
t234 = -rSges(5,3) * t255 + t254 * t288;
t229 = -Icges(5,3) * t254 - t255 * t282;
t228 = -Icges(5,3) * t255 + t254 * t282;
t227 = -rSges(6,3) * t254 - t255 * t287;
t226 = -rSges(6,3) * t255 + t254 * t287;
t219 = -t265 * (-rSges(4,1) * t255 + rSges(4,2) * t254) + t294;
t218 = t265 * (rSges(4,1) * t254 + rSges(4,2) * t255) + t295;
t217 = -pkin(8) * t254 - t255 * t302;
t216 = -pkin(8) * t255 + t254 * t302;
t215 = qJD(3) + (t234 * t254 - t235 * t255) * qJD(4);
t214 = -t251 * t293 + (-t235 - t241) * t265 + t294;
t213 = t234 * t265 + t251 * t292 + t289;
t212 = -t254 * t290 + t239 * t245 + (-t217 - t227 - t241) * t265 + t294;
t211 = t255 * t290 - t240 * t245 + (t216 + t226) * t265 + t289;
t210 = -t226 * t239 + t227 * t240 + qJD(3) + (t216 * t254 - t217 * t255) * qJD(4);
t1 = m(3) * (t236 ^ 2 + t237 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t218 ^ 2 + t219 ^ 2) / 0.2e1 + m(5) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 - ((-t255 * t246 - t254 * t275) * t265 + (t255 ^ 2 * t228 + (t277 * t254 + (t229 - t278) * t255) * t254) * qJD(4)) * t292 / 0.2e1 - ((-t254 * t246 + t255 * t275) * t265 + (t254 ^ 2 * t229 + (t278 * t255 + (t228 - t277) * t254) * t255) * qJD(4)) * t293 / 0.2e1 + m(6) * (t210 ^ 2 + t211 ^ 2 + t212 ^ 2) / 0.2e1 + t240 * (-t307 * t254 + t274 * t255) / 0.2e1 + t239 * (t274 * t254 + t307 * t255) / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t265 ^ 2 / 0.2e1 + ((-(t230 * t270 + t232 * t268) * t255 - (t231 * t270 + t233 * t268) * t254) * qJD(4) + (t222 * t262 + t224 * t260) * t240 + (t223 * t262 + t225 * t260) * t239 + (t262 * t243 + t260 * t244 + t270 * t247 + t268 * t248) * t265) * t265 / 0.2e1 + (m(2) * (t252 ^ 2 + t253 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
