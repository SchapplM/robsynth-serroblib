% Calculate kinetic energy for
% S4RRRR5
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
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:21
% EndTime: 2019-12-31 17:27:22
% DurationCPUTime: 1.21s
% Computational Cost: add. (670->204), mult. (1227->344), div. (0->0), fcn. (1226->8), ass. (0->104)
t276 = cos(qJ(3));
t305 = pkin(3) * t276;
t274 = sin(qJ(2));
t303 = Icges(3,4) * t274;
t277 = cos(qJ(2));
t302 = Icges(3,4) * t277;
t273 = sin(qJ(3));
t275 = sin(qJ(1));
t301 = t273 * t275;
t278 = cos(qJ(1));
t300 = t273 * t278;
t299 = t274 * t275;
t298 = t274 * t278;
t297 = t275 * t277;
t296 = t277 * t278;
t291 = pkin(2) * t277 + pkin(6) * t274;
t251 = t291 * t275;
t252 = t291 * t278;
t269 = qJD(2) * t275;
t294 = qJD(2) * t278;
t295 = t251 * t269 + t252 * t294;
t293 = qJD(3) * t274;
t253 = t278 * t293 + t269;
t292 = qJD(4) * t274;
t254 = t275 * t293 - t294;
t290 = rSges(3,1) * t277 - rSges(3,2) * t274;
t289 = Icges(3,1) * t277 - t303;
t288 = -Icges(3,2) * t274 + t302;
t287 = Icges(3,5) * t277 - Icges(3,6) * t274;
t232 = -Icges(3,6) * t278 + t275 * t288;
t235 = -Icges(3,5) * t278 + t275 * t289;
t286 = t232 * t274 - t235 * t277;
t233 = Icges(3,6) * t275 + t278 * t288;
t236 = Icges(3,5) * t275 + t278 * t289;
t285 = -t233 * t274 + t236 * t277;
t258 = Icges(3,2) * t277 + t303;
t259 = Icges(3,1) * t274 + t302;
t284 = -t258 * t274 + t259 * t277;
t256 = qJD(1) * (pkin(1) * t278 + pkin(5) * t275);
t263 = pkin(2) * t274 - pkin(6) * t277;
t283 = qJD(1) * t252 - t263 * t269 + t256;
t282 = pkin(7) * t274 + t277 * t305;
t264 = pkin(1) * t275 - pkin(5) * t278;
t281 = (-t251 - t264) * qJD(1) - t263 * t294;
t272 = qJ(3) + qJ(4);
t271 = cos(t272);
t270 = sin(t272);
t267 = -qJD(3) * t277 + qJD(1);
t262 = rSges(2,1) * t278 - rSges(2,2) * t275;
t261 = rSges(2,1) * t275 + rSges(2,2) * t278;
t260 = rSges(3,1) * t274 + rSges(3,2) * t277;
t257 = Icges(3,5) * t274 + Icges(3,6) * t277;
t255 = qJD(1) + (-qJD(3) - qJD(4)) * t277;
t250 = t276 * t296 + t301;
t249 = -t273 * t296 + t275 * t276;
t248 = t276 * t297 - t300;
t247 = -t273 * t297 - t276 * t278;
t243 = t270 * t275 + t271 * t296;
t242 = -t270 * t296 + t271 * t275;
t241 = -t270 * t278 + t271 * t297;
t240 = -t270 * t297 - t271 * t278;
t239 = rSges(3,3) * t275 + t278 * t290;
t238 = -rSges(3,3) * t278 + t275 * t290;
t237 = -rSges(4,3) * t277 + (rSges(4,1) * t276 - rSges(4,2) * t273) * t274;
t234 = -Icges(4,5) * t277 + (Icges(4,1) * t276 - Icges(4,4) * t273) * t274;
t231 = -Icges(4,6) * t277 + (Icges(4,4) * t276 - Icges(4,2) * t273) * t274;
t230 = Icges(3,3) * t275 + t278 * t287;
t229 = -Icges(3,3) * t278 + t275 * t287;
t228 = -Icges(4,3) * t277 + (Icges(4,5) * t276 - Icges(4,6) * t273) * t274;
t227 = t275 * t292 + t254;
t226 = t278 * t292 + t253;
t225 = -rSges(5,3) * t277 + (rSges(5,1) * t271 - rSges(5,2) * t270) * t274;
t224 = -Icges(5,5) * t277 + (Icges(5,1) * t271 - Icges(5,4) * t270) * t274;
t223 = -Icges(5,6) * t277 + (Icges(5,4) * t271 - Icges(5,2) * t270) * t274;
t222 = -Icges(5,3) * t277 + (Icges(5,5) * t271 - Icges(5,6) * t270) * t274;
t221 = -pkin(7) * t277 + t274 * t305;
t220 = pkin(3) * t301 + t278 * t282;
t219 = -pkin(3) * t300 + t275 * t282;
t218 = rSges(4,1) * t250 + rSges(4,2) * t249 + rSges(4,3) * t298;
t217 = rSges(4,1) * t248 + rSges(4,2) * t247 + rSges(4,3) * t299;
t216 = Icges(4,1) * t250 + Icges(4,4) * t249 + Icges(4,5) * t298;
t215 = Icges(4,1) * t248 + Icges(4,4) * t247 + Icges(4,5) * t299;
t214 = Icges(4,4) * t250 + Icges(4,2) * t249 + Icges(4,6) * t298;
t213 = Icges(4,4) * t248 + Icges(4,2) * t247 + Icges(4,6) * t299;
t212 = Icges(4,5) * t250 + Icges(4,6) * t249 + Icges(4,3) * t298;
t211 = Icges(4,5) * t248 + Icges(4,6) * t247 + Icges(4,3) * t299;
t210 = qJD(1) * t239 - t260 * t269 + t256;
t209 = -t260 * t294 + (-t238 - t264) * qJD(1);
t208 = (t238 * t275 + t239 * t278) * qJD(2);
t207 = rSges(5,1) * t243 + rSges(5,2) * t242 + rSges(5,3) * t298;
t206 = rSges(5,1) * t241 + rSges(5,2) * t240 + rSges(5,3) * t299;
t205 = Icges(5,1) * t243 + Icges(5,4) * t242 + Icges(5,5) * t298;
t204 = Icges(5,1) * t241 + Icges(5,4) * t240 + Icges(5,5) * t299;
t203 = Icges(5,4) * t243 + Icges(5,2) * t242 + Icges(5,6) * t298;
t202 = Icges(5,4) * t241 + Icges(5,2) * t240 + Icges(5,6) * t299;
t201 = Icges(5,5) * t243 + Icges(5,6) * t242 + Icges(5,3) * t298;
t200 = Icges(5,5) * t241 + Icges(5,6) * t240 + Icges(5,3) * t299;
t199 = t218 * t267 - t237 * t253 + t283;
t198 = -t217 * t267 + t237 * t254 + t281;
t197 = t217 * t253 - t218 * t254 + t295;
t196 = t207 * t255 + t220 * t267 - t221 * t253 - t225 * t226 + t283;
t195 = -t206 * t255 - t219 * t267 + t221 * t254 + t225 * t227 + t281;
t194 = t206 * t226 - t207 * t227 + t219 * t253 - t220 * t254 + t295;
t1 = m(3) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + ((t275 * t257 + t278 * t284) * qJD(1) + (t275 ^ 2 * t230 + (t286 * t278 + (-t229 + t285) * t275) * t278) * qJD(2)) * t269 / 0.2e1 - ((-t278 * t257 + t284 * t275) * qJD(1) + (t278 ^ 2 * t229 + (t285 * t275 + (-t230 + t286) * t278) * t275) * qJD(2)) * t294 / 0.2e1 + qJD(1) * ((t277 * t258 + t274 * t259) * qJD(1) + ((t233 * t277 + t236 * t274) * t275 - (t232 * t277 + t274 * t235) * t278) * qJD(2)) / 0.2e1 + m(4) * (t197 ^ 2 + t198 ^ 2 + t199 ^ 2) / 0.2e1 + t253 * ((t212 * t298 + t249 * t214 + t250 * t216) * t253 + (t211 * t298 + t213 * t249 + t215 * t250) * t254 + (t228 * t298 + t231 * t249 + t234 * t250) * t267) / 0.2e1 + t254 * ((t212 * t299 + t214 * t247 + t216 * t248) * t253 + (t211 * t299 + t247 * t213 + t248 * t215) * t254 + (t228 * t299 + t231 * t247 + t234 * t248) * t267) / 0.2e1 + t267 * ((-t211 * t254 - t212 * t253 - t228 * t267) * t277 + ((-t214 * t273 + t216 * t276) * t253 + (-t213 * t273 + t215 * t276) * t254 + (-t231 * t273 + t234 * t276) * t267) * t274) / 0.2e1 + m(5) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + t226 * ((t201 * t298 + t242 * t203 + t243 * t205) * t226 + (t200 * t298 + t202 * t242 + t204 * t243) * t227 + (t222 * t298 + t223 * t242 + t224 * t243) * t255) / 0.2e1 + t227 * ((t201 * t299 + t203 * t240 + t205 * t241) * t226 + (t200 * t299 + t240 * t202 + t241 * t204) * t227 + (t222 * t299 + t223 * t240 + t224 * t241) * t255) / 0.2e1 + t255 * ((-t200 * t227 - t201 * t226 - t222 * t255) * t277 + ((-t203 * t270 + t205 * t271) * t226 + (-t202 * t270 + t204 * t271) * t227 + (-t223 * t270 + t224 * t271) * t255) * t274) / 0.2e1 + (m(2) * (t261 ^ 2 + t262 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
