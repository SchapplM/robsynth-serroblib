% Calculate kinetic energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:24:52
% EndTime: 2019-07-18 13:24:54
% DurationCPUTime: 1.05s
% Computational Cost: add. (593->192), mult. (1436->326), div. (0->0), fcn. (1566->8), ass. (0->95)
t263 = sin(qJ(3));
t289 = Icges(4,4) * t263;
t267 = cos(qJ(3));
t288 = Icges(4,4) * t267;
t262 = sin(qJ(4));
t287 = t262 * t263;
t264 = sin(qJ(1));
t286 = t263 * t264;
t266 = cos(qJ(4));
t285 = t263 * t266;
t268 = cos(qJ(1));
t284 = t263 * t268;
t283 = t267 * t264;
t282 = t267 * t268;
t280 = qJ(2) * qJD(1);
t281 = qJD(2) * t264 + t268 * t280;
t259 = qJD(3) * t264;
t278 = qJD(4) * t263;
t245 = t268 * t278 + t259;
t279 = qJD(3) * t268;
t277 = -qJD(2) * t268 + t264 * t280;
t246 = t264 * t278 - t279;
t255 = -qJD(4) * t267 + qJD(1);
t276 = rSges(4,1) * t267 - rSges(4,2) * t263;
t275 = Icges(4,1) * t267 - t289;
t274 = -Icges(4,2) * t263 + t288;
t273 = Icges(4,5) * t267 - Icges(4,6) * t263;
t230 = -Icges(4,6) * t268 + t264 * t274;
t233 = -Icges(4,5) * t268 + t264 * t275;
t272 = t230 * t263 - t233 * t267;
t231 = Icges(4,6) * t264 + t268 * t274;
t234 = Icges(4,5) * t264 + t268 * t275;
t271 = -t231 * t263 + t234 * t267;
t248 = Icges(4,2) * t267 + t289;
t249 = Icges(4,1) * t263 + t288;
t270 = -t248 * t263 + t249 * t267;
t265 = cos(qJ(5));
t261 = sin(qJ(5));
t252 = rSges(2,1) * t268 - rSges(2,2) * t264;
t251 = rSges(2,1) * t264 + rSges(2,2) * t268;
t250 = rSges(4,1) * t263 + rSges(4,2) * t267;
t247 = Icges(4,5) * t263 + Icges(4,6) * t267;
t244 = qJD(5) * t287 + t255;
t243 = t262 * t264 + t266 * t282;
t242 = t262 * t282 - t264 * t266;
t241 = -t262 * t268 + t266 * t283;
t240 = t262 * t283 + t266 * t268;
t239 = -t261 * t267 + t265 * t285;
t238 = -t261 * t285 - t265 * t267;
t237 = rSges(4,3) * t264 + t268 * t276;
t236 = -rSges(4,3) * t268 + t264 * t276;
t235 = -rSges(5,3) * t267 + (rSges(5,1) * t266 - rSges(5,2) * t262) * t263;
t232 = -Icges(5,5) * t267 + (Icges(5,1) * t266 - Icges(5,4) * t262) * t263;
t229 = -Icges(5,6) * t267 + (Icges(5,4) * t266 - Icges(5,2) * t262) * t263;
t228 = Icges(4,3) * t264 + t268 * t273;
t227 = -Icges(4,3) * t268 + t264 * t273;
t226 = -Icges(5,3) * t267 + (Icges(5,5) * t266 - Icges(5,6) * t262) * t263;
t225 = -qJD(1) * (rSges(3,1) * t264 - rSges(3,3) * t268) + t281;
t224 = qJD(1) * (rSges(3,1) * t268 + rSges(3,3) * t264) + t277;
t223 = t243 * t265 + t261 * t284;
t222 = -t243 * t261 + t265 * t284;
t221 = t241 * t265 + t261 * t286;
t220 = -t241 * t261 + t265 * t286;
t219 = qJD(5) * t240 + t246;
t218 = qJD(5) * t242 + t245;
t217 = rSges(5,1) * t243 - rSges(5,2) * t242 + rSges(5,3) * t284;
t216 = rSges(5,1) * t241 - rSges(5,2) * t240 + rSges(5,3) * t286;
t215 = rSges(6,1) * t239 + rSges(6,2) * t238 + rSges(6,3) * t287;
t214 = Icges(5,1) * t243 - Icges(5,4) * t242 + Icges(5,5) * t284;
t213 = Icges(5,1) * t241 - Icges(5,4) * t240 + Icges(5,5) * t286;
t212 = Icges(6,1) * t239 + Icges(6,4) * t238 + Icges(6,5) * t287;
t211 = Icges(5,4) * t243 - Icges(5,2) * t242 + Icges(5,6) * t284;
t210 = Icges(5,4) * t241 - Icges(5,2) * t240 + Icges(5,6) * t286;
t209 = Icges(6,4) * t239 + Icges(6,2) * t238 + Icges(6,6) * t287;
t208 = Icges(5,5) * t243 - Icges(5,6) * t242 + Icges(5,3) * t284;
t207 = Icges(5,5) * t241 - Icges(5,6) * t240 + Icges(5,3) * t286;
t206 = Icges(6,5) * t239 + Icges(6,6) * t238 + Icges(6,3) * t287;
t205 = -qJD(1) * t236 - t250 * t279 + t281;
t204 = qJD(1) * t237 - t250 * t259 + t277;
t203 = (t236 * t264 + t237 * t268) * qJD(3);
t202 = rSges(6,1) * t223 + rSges(6,2) * t222 + rSges(6,3) * t242;
t201 = rSges(6,1) * t221 + rSges(6,2) * t220 + rSges(6,3) * t240;
t200 = Icges(6,1) * t223 + Icges(6,4) * t222 + Icges(6,5) * t242;
t199 = Icges(6,1) * t221 + Icges(6,4) * t220 + Icges(6,5) * t240;
t198 = Icges(6,4) * t223 + Icges(6,2) * t222 + Icges(6,6) * t242;
t197 = Icges(6,4) * t221 + Icges(6,2) * t220 + Icges(6,6) * t240;
t196 = Icges(6,5) * t223 + Icges(6,6) * t222 + Icges(6,3) * t242;
t195 = Icges(6,5) * t221 + Icges(6,6) * t220 + Icges(6,3) * t240;
t194 = -t216 * t255 + t235 * t246 + t281;
t193 = t217 * t255 - t235 * t245 + t277;
t192 = t216 * t245 - t217 * t246;
t191 = -t201 * t244 + t215 * t219 + t281;
t190 = t202 * t244 - t215 * t218 + t277;
t189 = t201 * t218 - t202 * t219;
t1 = m(3) * (t224 ^ 2 + t225 ^ 2) / 0.2e1 + m(4) * (t203 ^ 2 + t204 ^ 2 + t205 ^ 2) / 0.2e1 + ((t264 * t247 + t268 * t270) * qJD(1) + (t264 ^ 2 * t228 + (t272 * t268 + (-t227 + t271) * t264) * t268) * qJD(3)) * t259 / 0.2e1 - ((-t268 * t247 + t264 * t270) * qJD(1) + (t268 ^ 2 * t227 + (t271 * t264 + (-t228 + t272) * t268) * t264) * qJD(3)) * t279 / 0.2e1 + qJD(1) * ((t267 * t248 + t263 * t249) * qJD(1) + ((t231 * t267 + t234 * t263) * t264 - (t230 * t267 + t233 * t263) * t268) * qJD(3)) / 0.2e1 + m(5) * (t192 ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + t245 * ((t208 * t284 - t242 * t211 + t243 * t214) * t245 + (t207 * t284 - t210 * t242 + t213 * t243) * t246 + (t226 * t284 - t229 * t242 + t232 * t243) * t255) / 0.2e1 + t246 * ((t208 * t286 - t211 * t240 + t214 * t241) * t245 + (t207 * t286 - t240 * t210 + t241 * t213) * t246 + (t226 * t286 - t229 * t240 + t232 * t241) * t255) / 0.2e1 + t255 * ((-t207 * t246 - t208 * t245 - t226 * t255) * t267 + ((-t211 * t262 + t214 * t266) * t245 + (-t210 * t262 + t213 * t266) * t246 + (-t229 * t262 + t232 * t266) * t255) * t263) / 0.2e1 + m(6) * (t189 ^ 2 + t190 ^ 2 + t191 ^ 2) / 0.2e1 + t218 * ((t242 * t196 + t222 * t198 + t223 * t200) * t218 + (t195 * t242 + t197 * t222 + t199 * t223) * t219 + (t206 * t242 + t209 * t222 + t212 * t223) * t244) / 0.2e1 + t219 * ((t196 * t240 + t198 * t220 + t200 * t221) * t218 + (t240 * t195 + t220 * t197 + t221 * t199) * t219 + (t206 * t240 + t209 * t220 + t212 * t221) * t244) / 0.2e1 + t244 * ((t196 * t287 + t198 * t238 + t200 * t239) * t218 + (t195 * t287 + t197 * t238 + t199 * t239) * t219 + (t206 * t287 + t238 * t209 + t239 * t212) * t244) / 0.2e1 + (m(2) * (t251 ^ 2 + t252 ^ 2) + Icges(2,3) + Icges(3,2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
