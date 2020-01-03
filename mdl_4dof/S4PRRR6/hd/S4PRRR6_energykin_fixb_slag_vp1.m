% Calculate kinetic energy for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR6_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR6_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:38
% EndTime: 2019-12-31 16:34:39
% DurationCPUTime: 1.18s
% Computational Cost: add. (601->160), mult. (1108->279), div. (0->0), fcn. (1130->8), ass. (0->86)
t271 = cos(pkin(7));
t304 = t271 ^ 2;
t270 = sin(pkin(7));
t305 = t270 ^ 2;
t306 = t304 + t305;
t308 = qJD(2) * t306;
t307 = t270 * t271;
t303 = qJD(2) ^ 2;
t274 = cos(qJ(3));
t302 = pkin(3) * t274;
t272 = sin(qJ(3));
t300 = t270 * t272;
t273 = sin(qJ(2));
t299 = t270 * t273;
t275 = cos(qJ(2));
t298 = t270 * t275;
t297 = t271 * t272;
t296 = t271 * t273;
t295 = t271 * t275;
t294 = t272 * t275;
t293 = t274 * t275;
t266 = qJD(2) * t270;
t291 = qJD(3) * t273;
t257 = t271 * t291 + t266;
t292 = qJD(2) * t271;
t290 = qJD(3) * t275;
t289 = qJD(4) * t273;
t262 = pkin(2) * t273 - pkin(5) * t275;
t288 = t262 * t266;
t287 = t262 * t292;
t286 = qJD(1) + (pkin(2) * t275 + pkin(5) * t273) * t308;
t258 = t270 * t291 - t292;
t283 = Icges(3,5) * t275 - Icges(3,6) * t273;
t278 = pkin(6) * t273 + t275 * t302;
t269 = qJ(3) + qJ(4);
t268 = cos(t269);
t267 = sin(t269);
t261 = rSges(3,1) * t273 + rSges(3,2) * t275;
t259 = (-qJD(3) - qJD(4)) * t275;
t256 = t271 * t293 + t300;
t255 = t270 * t274 - t271 * t294;
t254 = t270 * t293 - t297;
t253 = -t270 * t294 - t271 * t274;
t252 = -rSges(4,3) * t275 + (rSges(4,1) * t274 - rSges(4,2) * t272) * t273;
t251 = -Icges(4,5) * t275 + (Icges(4,1) * t274 - Icges(4,4) * t272) * t273;
t250 = -Icges(4,6) * t275 + (Icges(4,4) * t274 - Icges(4,2) * t272) * t273;
t249 = -Icges(4,3) * t275 + (Icges(4,5) * t274 - Icges(4,6) * t272) * t273;
t248 = t267 * t270 + t268 * t295;
t247 = -t267 * t295 + t268 * t270;
t246 = -t267 * t271 + t268 * t298;
t245 = -t267 * t298 - t268 * t271;
t238 = Icges(3,3) * t270 + t271 * t283;
t237 = -Icges(3,3) * t271 + t270 * t283;
t236 = t270 * t289 + t258;
t235 = t271 * t289 + t257;
t234 = -rSges(5,3) * t275 + (rSges(5,1) * t268 - rSges(5,2) * t267) * t273;
t233 = -Icges(5,5) * t275 + (Icges(5,1) * t268 - Icges(5,4) * t267) * t273;
t232 = -Icges(5,6) * t275 + (Icges(5,4) * t268 - Icges(5,2) * t267) * t273;
t231 = -Icges(5,3) * t275 + (Icges(5,5) * t268 - Icges(5,6) * t267) * t273;
t230 = -pkin(6) * t275 + t273 * t302;
t229 = pkin(3) * t300 + t271 * t278;
t228 = -pkin(3) * t297 + t270 * t278;
t227 = rSges(4,1) * t256 + rSges(4,2) * t255 + rSges(4,3) * t296;
t226 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t299;
t225 = Icges(4,1) * t256 + Icges(4,4) * t255 + Icges(4,5) * t296;
t224 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t299;
t223 = Icges(4,4) * t256 + Icges(4,2) * t255 + Icges(4,6) * t296;
t222 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t299;
t221 = Icges(4,5) * t256 + Icges(4,6) * t255 + Icges(4,3) * t296;
t220 = Icges(4,5) * t254 + Icges(4,6) * t253 + Icges(4,3) * t299;
t219 = rSges(5,1) * t248 + rSges(5,2) * t247 + rSges(5,3) * t296;
t218 = rSges(5,1) * t246 + rSges(5,2) * t245 + rSges(5,3) * t299;
t217 = Icges(5,1) * t248 + Icges(5,4) * t247 + Icges(5,5) * t296;
t216 = Icges(5,1) * t246 + Icges(5,4) * t245 + Icges(5,5) * t299;
t215 = Icges(5,4) * t248 + Icges(5,2) * t247 + Icges(5,6) * t296;
t214 = Icges(5,4) * t246 + Icges(5,2) * t245 + Icges(5,6) * t299;
t213 = Icges(5,5) * t248 + Icges(5,6) * t247 + Icges(5,3) * t296;
t212 = Icges(5,5) * t246 + Icges(5,6) * t245 + Icges(5,3) * t299;
t211 = qJD(1) + (rSges(3,1) * t275 - rSges(3,2) * t273) * t308;
t210 = t226 * t290 + t252 * t258 - t287;
t209 = -t227 * t290 - t252 * t257 - t288;
t208 = t226 * t257 - t227 * t258 + t286;
t207 = -t218 * t259 + t228 * t290 + t230 * t258 + t234 * t236 - t287;
t206 = t219 * t259 - t229 * t290 - t230 * t257 - t234 * t235 - t288;
t205 = t218 * t235 - t219 * t236 + t228 * t257 - t229 * t258 + t286;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t306 * t303 * t261 ^ 2 + t211 ^ 2) / 0.2e1 + t303 * t270 * (-t237 * t307 + t305 * t238) / 0.2e1 - t303 * t271 * (t304 * t237 - t238 * t307) / 0.2e1 + m(4) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + t257 * ((t221 * t296 + t255 * t223 + t256 * t225) * t257 + (t220 * t296 + t222 * t255 + t224 * t256) * t258 - (t249 * t296 + t250 * t255 + t251 * t256) * t290) / 0.2e1 + t258 * ((t221 * t299 + t223 * t253 + t225 * t254) * t257 + (t220 * t299 + t253 * t222 + t254 * t224) * t258 - (t249 * t299 + t250 * t253 + t251 * t254) * t290) / 0.2e1 - ((-t220 * t258 - t221 * t257 + t249 * t290) * t275 + ((-t223 * t272 + t225 * t274) * t257 + (-t222 * t272 + t224 * t274) * t258 - (-t250 * t272 + t251 * t274) * t290) * t273) * t290 / 0.2e1 + m(5) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + t235 * ((t213 * t296 + t247 * t215 + t248 * t217) * t235 + (t212 * t296 + t214 * t247 + t216 * t248) * t236 + (t231 * t296 + t232 * t247 + t233 * t248) * t259) / 0.2e1 + t236 * ((t213 * t299 + t215 * t245 + t217 * t246) * t235 + (t212 * t299 + t245 * t214 + t246 * t216) * t236 + (t231 * t299 + t232 * t245 + t233 * t246) * t259) / 0.2e1 + t259 * ((-t212 * t236 - t213 * t235 - t231 * t259) * t275 + ((-t215 * t267 + t217 * t268) * t235 + (-t214 * t267 + t216 * t268) * t236 + (-t232 * t267 + t233 * t268) * t259) * t273) / 0.2e1;
T = t1;
