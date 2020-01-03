% Calculate kinetic energy for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR4_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR4_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:40
% EndTime: 2019-12-31 17:25:41
% DurationCPUTime: 0.96s
% Computational Cost: add. (675->179), mult. (961->302), div. (0->0), fcn. (904->8), ass. (0->103)
t269 = cos(qJ(2));
t305 = pkin(2) * t269;
t266 = sin(qJ(2));
t303 = Icges(3,4) * t266;
t302 = Icges(3,4) * t269;
t264 = qJ(2) + qJ(3);
t262 = sin(t264);
t301 = Icges(4,4) * t262;
t263 = cos(t264);
t300 = Icges(4,4) * t263;
t267 = sin(qJ(1));
t299 = t262 * t267;
t270 = cos(qJ(1));
t298 = t262 * t270;
t265 = sin(qJ(4));
t297 = t265 * t267;
t296 = t265 * t270;
t268 = cos(qJ(4));
t295 = t267 * t268;
t294 = t268 * t270;
t217 = -pkin(6) * t270 + t267 * t305;
t218 = pkin(6) * t267 + t270 * t305;
t261 = qJD(2) * t267;
t291 = qJD(2) * t270;
t293 = t217 * t261 + t218 * t291;
t257 = pkin(1) * t267 - pkin(5) * t270;
t292 = -t217 - t257;
t249 = qJD(3) * t267 + t261;
t290 = qJD(4) * t262;
t289 = pkin(2) * qJD(2) * t266;
t250 = (-qJD(2) - qJD(3)) * t270;
t288 = t270 * t289;
t287 = pkin(3) * t263 + pkin(7) * t262;
t286 = rSges(3,1) * t269 - rSges(3,2) * t266;
t285 = rSges(4,1) * t263 - rSges(4,2) * t262;
t284 = Icges(3,1) * t269 - t303;
t283 = Icges(4,1) * t263 - t301;
t282 = -Icges(3,2) * t266 + t302;
t281 = -Icges(4,2) * t262 + t300;
t280 = Icges(3,5) * t269 - Icges(3,6) * t266;
t279 = Icges(4,5) * t263 - Icges(4,6) * t262;
t229 = -Icges(3,6) * t270 + t267 * t282;
t231 = -Icges(3,5) * t270 + t267 * t284;
t278 = t229 * t266 - t231 * t269;
t230 = Icges(3,6) * t267 + t270 * t282;
t232 = Icges(3,5) * t267 + t270 * t284;
t277 = -t230 * t266 + t232 * t269;
t252 = Icges(3,2) * t269 + t303;
t253 = Icges(3,1) * t266 + t302;
t276 = -t252 * t266 + t253 * t269;
t248 = qJD(1) * (pkin(1) * t270 + pkin(5) * t267);
t275 = qJD(1) * t218 - t267 * t289 + t248;
t274 = (Icges(4,5) * t262 + Icges(4,6) * t263) * qJD(1) + (-Icges(4,3) * t270 + t267 * t279) * t250 + (Icges(4,3) * t267 + t270 * t279) * t249;
t221 = -Icges(4,6) * t270 + t267 * t281;
t222 = Icges(4,6) * t267 + t270 * t281;
t223 = -Icges(4,5) * t270 + t267 * t283;
t224 = Icges(4,5) * t267 + t270 * t283;
t244 = Icges(4,2) * t263 + t301;
t245 = Icges(4,1) * t262 + t300;
t273 = (-t222 * t262 + t224 * t263) * t249 + (-t221 * t262 + t223 * t263) * t250 + (-t244 * t262 + t245 * t263) * qJD(1);
t258 = -qJD(4) * t263 + qJD(1);
t256 = rSges(2,1) * t270 - rSges(2,2) * t267;
t255 = rSges(2,1) * t267 + rSges(2,2) * t270;
t254 = rSges(3,1) * t266 + rSges(3,2) * t269;
t251 = Icges(3,5) * t266 + Icges(3,6) * t269;
t247 = pkin(3) * t262 - pkin(7) * t263;
t246 = rSges(4,1) * t262 + rSges(4,2) * t263;
t242 = t263 * t294 + t297;
t241 = -t263 * t296 + t295;
t240 = t263 * t295 - t296;
t239 = -t263 * t297 - t294;
t238 = t287 * t270;
t237 = t287 * t267;
t236 = rSges(3,3) * t267 + t270 * t286;
t235 = -rSges(3,3) * t270 + t267 * t286;
t234 = t267 * t290 + t250;
t233 = t270 * t290 + t249;
t228 = Icges(3,3) * t267 + t270 * t280;
t227 = -Icges(3,3) * t270 + t267 * t280;
t226 = rSges(4,3) * t267 + t270 * t285;
t225 = -rSges(4,3) * t270 + t267 * t285;
t216 = -rSges(5,3) * t263 + (rSges(5,1) * t268 - rSges(5,2) * t265) * t262;
t215 = -Icges(5,5) * t263 + (Icges(5,1) * t268 - Icges(5,4) * t265) * t262;
t214 = -Icges(5,6) * t263 + (Icges(5,4) * t268 - Icges(5,2) * t265) * t262;
t213 = -Icges(5,3) * t263 + (Icges(5,5) * t268 - Icges(5,6) * t265) * t262;
t209 = qJD(1) * t236 - t254 * t261 + t248;
t208 = -t254 * t291 + (-t235 - t257) * qJD(1);
t207 = rSges(5,1) * t242 + rSges(5,2) * t241 + rSges(5,3) * t298;
t206 = rSges(5,1) * t240 + rSges(5,2) * t239 + rSges(5,3) * t299;
t205 = Icges(5,1) * t242 + Icges(5,4) * t241 + Icges(5,5) * t298;
t204 = Icges(5,1) * t240 + Icges(5,4) * t239 + Icges(5,5) * t299;
t203 = Icges(5,4) * t242 + Icges(5,2) * t241 + Icges(5,6) * t298;
t202 = Icges(5,4) * t240 + Icges(5,2) * t239 + Icges(5,6) * t299;
t201 = Icges(5,5) * t242 + Icges(5,6) * t241 + Icges(5,3) * t298;
t200 = Icges(5,5) * t240 + Icges(5,6) * t239 + Icges(5,3) * t299;
t199 = (t235 * t267 + t236 * t270) * qJD(2);
t198 = qJD(1) * t226 - t246 * t249 + t275;
t197 = -t288 + t246 * t250 + (-t225 + t292) * qJD(1);
t196 = t225 * t249 - t226 * t250 + t293;
t195 = qJD(1) * t238 + t207 * t258 - t216 * t233 - t247 * t249 + t275;
t194 = -t288 - t206 * t258 + t216 * t234 + t247 * t250 + (-t237 + t292) * qJD(1);
t193 = t206 * t233 - t207 * t234 + t237 * t249 - t238 * t250 + t293;
t1 = m(3) * (t199 ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + ((t267 * t251 + t276 * t270) * qJD(1) + (t267 ^ 2 * t228 + (t278 * t270 + (-t227 + t277) * t267) * t270) * qJD(2)) * t261 / 0.2e1 - ((-t270 * t251 + t267 * t276) * qJD(1) + (t270 ^ 2 * t227 + (t277 * t267 + (-t228 + t278) * t270) * t267) * qJD(2)) * t291 / 0.2e1 + m(4) * (t196 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + t249 * (t274 * t267 + t273 * t270) / 0.2e1 + t250 * (t273 * t267 - t274 * t270) / 0.2e1 + m(5) * (t193 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + t233 * ((t201 * t298 + t203 * t241 + t205 * t242) * t233 + (t200 * t298 + t202 * t241 + t242 * t204) * t234 + (t213 * t298 + t214 * t241 + t215 * t242) * t258) / 0.2e1 + t234 * ((t201 * t299 + t203 * t239 + t205 * t240) * t233 + (t200 * t299 + t202 * t239 + t204 * t240) * t234 + (t213 * t299 + t214 * t239 + t215 * t240) * t258) / 0.2e1 + t258 * ((-t200 * t234 - t201 * t233 - t213 * t258) * t263 + ((-t203 * t265 + t205 * t268) * t233 + (-t202 * t265 + t204 * t268) * t234 + (-t214 * t265 + t215 * t268) * t258) * t262) / 0.2e1 + (m(2) * (t255 ^ 2 + t256 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t230 * t269 + t232 * t266) * t267 - (t229 * t269 + t266 * t231) * t270) * qJD(2) + (t222 * t263 + t224 * t262) * t249 + (t221 * t263 + t223 * t262) * t250 + (t263 * t244 + t262 * t245 + t269 * t252 + t266 * t253) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
