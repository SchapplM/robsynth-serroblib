% Calculate kinetic energy for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR9_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR9_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:20
% EndTime: 2019-12-31 17:09:21
% DurationCPUTime: 1.19s
% Computational Cost: add. (634->204), mult. (1167->331), div. (0->0), fcn. (1166->8), ass. (0->102)
t270 = cos(pkin(7));
t306 = pkin(3) * t270;
t272 = sin(qJ(2));
t305 = Icges(3,4) * t272;
t274 = cos(qJ(2));
t304 = Icges(3,4) * t274;
t269 = sin(pkin(7));
t273 = sin(qJ(1));
t303 = t269 * t273;
t275 = cos(qJ(1));
t302 = t269 * t275;
t301 = t272 * t273;
t300 = t272 * t275;
t299 = t273 * t274;
t298 = t274 * t275;
t284 = pkin(2) * t274 + qJ(3) * t272;
t249 = t284 * t273;
t261 = pkin(1) * t273 - pkin(5) * t275;
t296 = -t249 - t261;
t295 = qJD(2) * t273;
t294 = qJD(2) * t275;
t293 = qJD(3) * t272;
t292 = qJD(4) * t272;
t250 = t284 * t275;
t253 = qJD(1) * (pkin(1) * t275 + pkin(5) * t273);
t291 = qJD(1) * t250 + t273 * t293 + t253;
t257 = pkin(2) * t272 - qJ(3) * t274;
t288 = qJD(2) * (pkin(6) * t274 - t272 * t306 - t257);
t287 = qJD(2) * (rSges(4,3) * t274 - (rSges(4,1) * t270 - rSges(4,2) * t269) * t272 - t257);
t286 = -qJD(3) * t274 + t249 * t295 + t250 * t294;
t285 = rSges(3,1) * t274 - rSges(3,2) * t272;
t283 = Icges(3,1) * t274 - t305;
t282 = -Icges(3,2) * t272 + t304;
t281 = Icges(3,5) * t274 - Icges(3,6) * t272;
t236 = -Icges(3,6) * t275 + t273 * t282;
t238 = -Icges(3,5) * t275 + t273 * t283;
t280 = t236 * t272 - t238 * t274;
t237 = Icges(3,6) * t273 + t275 * t282;
t239 = Icges(3,5) * t273 + t275 * t283;
t279 = -t237 * t272 + t239 * t274;
t255 = Icges(3,2) * t274 + t305;
t256 = Icges(3,1) * t272 + t304;
t278 = -t255 * t272 + t256 * t274;
t277 = pkin(6) * t272 + t274 * t306;
t268 = pkin(7) + qJ(4);
t267 = cos(t268);
t266 = sin(t268);
t264 = -qJD(4) * t274 + qJD(1);
t263 = t275 * t293;
t260 = rSges(2,1) * t275 - rSges(2,2) * t273;
t259 = rSges(2,1) * t273 + rSges(2,2) * t275;
t258 = rSges(3,1) * t272 + rSges(3,2) * t274;
t254 = Icges(3,5) * t272 + Icges(3,6) * t274;
t252 = t273 * t292 - t294;
t251 = t275 * t292 + t295;
t248 = t270 * t298 + t303;
t247 = -t269 * t298 + t270 * t273;
t246 = t270 * t299 - t302;
t245 = -t269 * t299 - t270 * t275;
t243 = rSges(3,3) * t273 + t275 * t285;
t242 = -rSges(3,3) * t275 + t273 * t285;
t235 = Icges(3,3) * t273 + t275 * t281;
t234 = -Icges(3,3) * t275 + t273 * t281;
t233 = t266 * t273 + t267 * t298;
t232 = -t266 * t298 + t267 * t273;
t231 = -t266 * t275 + t267 * t299;
t230 = -t266 * t299 - t267 * t275;
t228 = -Icges(4,5) * t274 + (Icges(4,1) * t270 - Icges(4,4) * t269) * t272;
t227 = -Icges(4,6) * t274 + (Icges(4,4) * t270 - Icges(4,2) * t269) * t272;
t226 = -Icges(4,3) * t274 + (Icges(4,5) * t270 - Icges(4,6) * t269) * t272;
t225 = -rSges(5,3) * t274 + (rSges(5,1) * t267 - rSges(5,2) * t266) * t272;
t224 = -Icges(5,5) * t274 + (Icges(5,1) * t267 - Icges(5,4) * t266) * t272;
t223 = -Icges(5,6) * t274 + (Icges(5,4) * t267 - Icges(5,2) * t266) * t272;
t222 = -Icges(5,3) * t274 + (Icges(5,5) * t267 - Icges(5,6) * t266) * t272;
t220 = pkin(3) * t303 + t275 * t277;
t219 = -pkin(3) * t302 + t273 * t277;
t218 = rSges(4,1) * t248 + rSges(4,2) * t247 + rSges(4,3) * t300;
t217 = rSges(4,1) * t246 + rSges(4,2) * t245 + rSges(4,3) * t301;
t216 = Icges(4,1) * t248 + Icges(4,4) * t247 + Icges(4,5) * t300;
t215 = Icges(4,1) * t246 + Icges(4,4) * t245 + Icges(4,5) * t301;
t214 = Icges(4,4) * t248 + Icges(4,2) * t247 + Icges(4,6) * t300;
t213 = Icges(4,4) * t246 + Icges(4,2) * t245 + Icges(4,6) * t301;
t212 = Icges(4,5) * t248 + Icges(4,6) * t247 + Icges(4,3) * t300;
t211 = Icges(4,5) * t246 + Icges(4,6) * t245 + Icges(4,3) * t301;
t210 = qJD(1) * t243 - t258 * t295 + t253;
t209 = -t258 * t294 + (-t242 - t261) * qJD(1);
t208 = (t242 * t273 + t243 * t275) * qJD(2);
t207 = rSges(5,1) * t233 + rSges(5,2) * t232 + rSges(5,3) * t300;
t206 = rSges(5,1) * t231 + rSges(5,2) * t230 + rSges(5,3) * t301;
t205 = Icges(5,1) * t233 + Icges(5,4) * t232 + Icges(5,5) * t300;
t204 = Icges(5,1) * t231 + Icges(5,4) * t230 + Icges(5,5) * t301;
t203 = Icges(5,4) * t233 + Icges(5,2) * t232 + Icges(5,6) * t300;
t202 = Icges(5,4) * t231 + Icges(5,2) * t230 + Icges(5,6) * t301;
t201 = Icges(5,5) * t233 + Icges(5,6) * t232 + Icges(5,3) * t300;
t200 = Icges(5,5) * t231 + Icges(5,6) * t230 + Icges(5,3) * t301;
t199 = qJD(1) * t218 + t273 * t287 + t291;
t198 = t263 + t275 * t287 + (-t217 + t296) * qJD(1);
t197 = (t217 * t273 + t218 * t275) * qJD(2) + t286;
t196 = qJD(1) * t220 + t207 * t264 - t225 * t251 + t273 * t288 + t291;
t195 = -t206 * t264 + t225 * t252 + t263 + t275 * t288 + (-t219 + t296) * qJD(1);
t194 = t206 * t251 - t207 * t252 + (t219 * t273 + t220 * t275) * qJD(2) + t286;
t1 = m(3) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + m(4) * (t197 ^ 2 + t198 ^ 2 + t199 ^ 2) / 0.2e1 + m(5) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + t251 * ((t201 * t300 + t232 * t203 + t233 * t205) * t251 + (t200 * t300 + t202 * t232 + t204 * t233) * t252 + (t222 * t300 + t223 * t232 + t224 * t233) * t264) / 0.2e1 + t252 * ((t201 * t301 + t203 * t230 + t205 * t231) * t251 + (t200 * t301 + t230 * t202 + t231 * t204) * t252 + (t222 * t301 + t223 * t230 + t224 * t231) * t264) / 0.2e1 + t264 * ((-t200 * t252 - t201 * t251 - t222 * t264) * t274 + ((-t203 * t266 + t205 * t267) * t251 + (-t202 * t266 + t204 * t267) * t252 + (-t223 * t266 + t224 * t267) * t264) * t272) / 0.2e1 + (m(2) * (t259 ^ 2 + t260 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t237 * t274 + t239 * t272) * t273 - (t236 * t274 + t238 * t272) * t275 + (t211 * t275 - t212 * t273) * t274 + ((-t214 * t269 + t216 * t270) * t273 - (-t213 * t269 + t215 * t270) * t275) * t272) * qJD(2) + ((t255 - t226) * t274 + (-t227 * t269 + t228 * t270 + t256) * t272) * qJD(1)) * qJD(1) / 0.2e1 + (((-t211 * t300 - t213 * t247 - t215 * t248 + t280 * t275) * t275 + ((-t234 + t279) * t275 + t212 * t300 + t214 * t247 + t216 * t248 + t235 * t273) * t273) * qJD(2) + (t226 * t300 + t227 * t247 + t228 * t248 + t273 * t254 + t275 * t278) * qJD(1)) * t295 / 0.2e1 - (((-t211 * t301 - t213 * t245 - t215 * t246 + t234 * t275) * t275 + ((-t235 + t280) * t275 + t212 * t301 + t214 * t245 + t216 * t246 + t279 * t273) * t273) * qJD(2) + (t226 * t301 + t227 * t245 + t228 * t246 - t275 * t254 + t273 * t278) * qJD(1)) * t294 / 0.2e1;
T = t1;
