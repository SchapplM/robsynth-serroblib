% Calculate kinetic energy for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:42:20
% EndTime: 2019-03-09 00:42:24
% DurationCPUTime: 4.20s
% Computational Cost: add. (3795->443), mult. (6268->660), div. (0->0), fcn. (7460->14), ass. (0->196)
t316 = cos(pkin(6));
t365 = pkin(7) * t316;
t321 = cos(qJ(3));
t364 = pkin(3) * t321;
t320 = cos(qJ(5));
t363 = pkin(5) * t320;
t313 = sin(pkin(12));
t360 = Icges(2,4) * t313;
t315 = cos(pkin(12));
t319 = sin(qJ(2));
t322 = cos(qJ(2));
t349 = t316 * t322;
t273 = t313 * t319 - t315 * t349;
t317 = sin(qJ(5));
t359 = t273 * t317;
t275 = t313 * t349 + t315 * t319;
t358 = t275 * t317;
t314 = sin(pkin(6));
t357 = t313 * t314;
t356 = t314 * t315;
t318 = sin(qJ(3));
t355 = t314 * t318;
t354 = t314 * t319;
t353 = t314 * t321;
t352 = t314 * t322;
t351 = t316 * t318;
t350 = t316 * t319;
t348 = qJ(3) + qJ(4);
t347 = qJD(2) * t314;
t346 = V_base(5) * qJ(1) + V_base(1);
t342 = qJD(1) + V_base(3);
t341 = t313 * t355;
t340 = t315 * t355;
t339 = t317 * t352;
t289 = t313 * t347 + V_base(4);
t300 = qJD(2) * t316 + V_base(6);
t338 = cos(t348);
t253 = qJD(3) * t275 + t289;
t337 = t314 * t338;
t222 = qJD(4) * t275 + t253;
t288 = -t315 * t347 + V_base(5);
t276 = -t313 * t350 + t315 * t322;
t309 = sin(t348);
t247 = t276 * t309 - t313 * t337;
t195 = qJD(5) * t247 + t222;
t252 = qJD(3) * t273 + t288;
t283 = pkin(1) * t313 - pkin(7) * t356;
t336 = -t283 * V_base(6) + V_base(5) * t365 + t346;
t284 = pkin(1) * t315 + pkin(7) * t357;
t335 = V_base(4) * t283 - t284 * V_base(5) + t342;
t221 = qJD(4) * t273 + t252;
t261 = (-qJD(3) - qJD(4)) * t352 + t300;
t274 = t313 * t322 + t315 * t350;
t245 = t274 * t309 + t315 * t337;
t194 = qJD(5) * t245 + t221;
t334 = V_base(6) * t284 + V_base(2) + (-qJ(1) - t365) * V_base(4);
t267 = t309 * t354 - t316 * t338;
t223 = qJD(5) * t267 + t261;
t241 = t274 * pkin(2) + t273 * pkin(8);
t282 = (pkin(2) * t319 - pkin(8) * t322) * t314;
t333 = -t241 * t300 + t288 * t282 + t336;
t242 = t276 * pkin(2) + t275 * pkin(8);
t332 = t289 * t241 - t242 * t288 + t335;
t331 = t300 * t242 - t282 * t289 + t334;
t192 = -pkin(3) * t340 + pkin(9) * t273 + t274 * t364;
t238 = pkin(3) * t351 + (-pkin(9) * t322 + t319 * t364) * t314;
t277 = -qJD(3) * t352 + t300;
t330 = -t192 * t277 + t252 * t238 + t333;
t193 = pkin(3) * t341 + pkin(9) * t275 + t276 * t364;
t329 = t253 * t192 - t193 * t252 + t332;
t328 = t277 * t193 - t238 * t253 + t331;
t246 = t274 * t338 - t309 * t356;
t205 = pkin(4) * t246 + pkin(10) * t245;
t268 = t316 * t309 + t319 * t337;
t237 = pkin(4) * t268 + pkin(10) * t267;
t327 = -t205 * t261 + t221 * t237 + t330;
t248 = t276 * t338 + t309 * t357;
t206 = pkin(4) * t248 + pkin(10) * t247;
t326 = t222 * t205 - t206 * t221 + t329;
t325 = t261 * t206 - t222 * t237 + t328;
t312 = qJ(5) + qJ(6);
t310 = cos(t312);
t308 = sin(t312);
t307 = Icges(2,4) * t315;
t298 = rSges(2,1) * t315 - rSges(2,2) * t313;
t297 = rSges(2,1) * t313 + rSges(2,2) * t315;
t296 = Icges(2,1) * t315 - t360;
t295 = Icges(2,1) * t313 + t307;
t294 = -Icges(2,2) * t313 + t307;
t293 = Icges(2,2) * t315 + t360;
t287 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t286 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t285 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t281 = t319 * t353 + t351;
t280 = t316 * t321 - t318 * t354;
t266 = rSges(3,3) * t316 + (rSges(3,1) * t319 + rSges(3,2) * t322) * t314;
t265 = Icges(3,5) * t316 + (Icges(3,1) * t319 + Icges(3,4) * t322) * t314;
t264 = Icges(3,6) * t316 + (Icges(3,4) * t319 + Icges(3,2) * t322) * t314;
t263 = Icges(3,3) * t316 + (Icges(3,5) * t319 + Icges(3,6) * t322) * t314;
t260 = V_base(5) * rSges(2,3) - t297 * V_base(6) + t346;
t259 = t298 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t257 = t276 * t321 + t341;
t256 = -t276 * t318 + t313 * t353;
t255 = t274 * t321 - t340;
t254 = -t274 * t318 - t315 * t353;
t251 = t297 * V_base(4) - t298 * V_base(5) + t342;
t250 = t268 * t320 - t339;
t249 = -t268 * t317 - t320 * t352;
t240 = t268 * t310 - t308 * t352;
t239 = -t268 * t308 - t310 * t352;
t236 = rSges(4,1) * t281 + rSges(4,2) * t280 - rSges(4,3) * t352;
t235 = Icges(4,1) * t281 + Icges(4,4) * t280 - Icges(4,5) * t352;
t234 = Icges(4,4) * t281 + Icges(4,2) * t280 - Icges(4,6) * t352;
t233 = Icges(4,5) * t281 + Icges(4,6) * t280 - Icges(4,3) * t352;
t232 = rSges(3,1) * t276 - rSges(3,2) * t275 + rSges(3,3) * t357;
t231 = rSges(3,1) * t274 - rSges(3,2) * t273 - rSges(3,3) * t356;
t230 = Icges(3,1) * t276 - Icges(3,4) * t275 + Icges(3,5) * t357;
t229 = Icges(3,1) * t274 - Icges(3,4) * t273 - Icges(3,5) * t356;
t228 = Icges(3,4) * t276 - Icges(3,2) * t275 + Icges(3,6) * t357;
t227 = Icges(3,4) * t274 - Icges(3,2) * t273 - Icges(3,6) * t356;
t226 = Icges(3,5) * t276 - Icges(3,6) * t275 + Icges(3,3) * t357;
t225 = Icges(3,5) * t274 - Icges(3,6) * t273 - Icges(3,3) * t356;
t219 = rSges(5,1) * t268 - rSges(5,2) * t267 - rSges(5,3) * t352;
t218 = Icges(5,1) * t268 - Icges(5,4) * t267 - Icges(5,5) * t352;
t217 = Icges(5,4) * t268 - Icges(5,2) * t267 - Icges(5,6) * t352;
t216 = Icges(5,5) * t268 - Icges(5,6) * t267 - Icges(5,3) * t352;
t215 = t248 * t320 + t358;
t214 = -t248 * t317 + t275 * t320;
t213 = t246 * t320 + t359;
t212 = -t246 * t317 + t273 * t320;
t211 = t248 * t310 + t275 * t308;
t210 = -t248 * t308 + t275 * t310;
t209 = t246 * t310 + t273 * t308;
t208 = -t246 * t308 + t273 * t310;
t207 = qJD(6) * t267 + t223;
t203 = rSges(4,1) * t257 + rSges(4,2) * t256 + rSges(4,3) * t275;
t202 = rSges(4,1) * t255 + rSges(4,2) * t254 + rSges(4,3) * t273;
t201 = Icges(4,1) * t257 + Icges(4,4) * t256 + Icges(4,5) * t275;
t200 = Icges(4,1) * t255 + Icges(4,4) * t254 + Icges(4,5) * t273;
t199 = Icges(4,4) * t257 + Icges(4,2) * t256 + Icges(4,6) * t275;
t198 = Icges(4,4) * t255 + Icges(4,2) * t254 + Icges(4,6) * t273;
t197 = Icges(4,5) * t257 + Icges(4,6) * t256 + Icges(4,3) * t275;
t196 = Icges(4,5) * t255 + Icges(4,6) * t254 + Icges(4,3) * t273;
t191 = rSges(5,1) * t248 - rSges(5,2) * t247 + rSges(5,3) * t275;
t190 = rSges(5,1) * t246 - rSges(5,2) * t245 + rSges(5,3) * t273;
t189 = Icges(5,1) * t248 - Icges(5,4) * t247 + Icges(5,5) * t275;
t188 = Icges(5,1) * t246 - Icges(5,4) * t245 + Icges(5,5) * t273;
t187 = Icges(5,4) * t248 - Icges(5,2) * t247 + Icges(5,6) * t275;
t186 = Icges(5,4) * t246 - Icges(5,2) * t245 + Icges(5,6) * t273;
t185 = Icges(5,5) * t248 - Icges(5,6) * t247 + Icges(5,3) * t275;
t184 = Icges(5,5) * t246 - Icges(5,6) * t245 + Icges(5,3) * t273;
t183 = rSges(6,1) * t250 + rSges(6,2) * t249 + rSges(6,3) * t267;
t182 = Icges(6,1) * t250 + Icges(6,4) * t249 + Icges(6,5) * t267;
t181 = Icges(6,4) * t250 + Icges(6,2) * t249 + Icges(6,6) * t267;
t180 = Icges(6,5) * t250 + Icges(6,6) * t249 + Icges(6,3) * t267;
t177 = rSges(7,1) * t240 + rSges(7,2) * t239 + rSges(7,3) * t267;
t176 = Icges(7,1) * t240 + Icges(7,4) * t239 + Icges(7,5) * t267;
t175 = Icges(7,4) * t240 + Icges(7,2) * t239 + Icges(7,6) * t267;
t174 = Icges(7,5) * t240 + Icges(7,6) * t239 + Icges(7,3) * t267;
t173 = -pkin(5) * t339 + pkin(11) * t267 + t268 * t363;
t171 = -t231 * t300 + t266 * t288 + t336;
t170 = t232 * t300 - t266 * t289 + t334;
t168 = qJD(6) * t247 + t195;
t167 = qJD(6) * t245 + t194;
t165 = t231 * t289 - t232 * t288 + t335;
t164 = rSges(6,1) * t215 + rSges(6,2) * t214 + rSges(6,3) * t247;
t163 = rSges(6,1) * t213 + rSges(6,2) * t212 + rSges(6,3) * t245;
t162 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t247;
t161 = Icges(6,1) * t213 + Icges(6,4) * t212 + Icges(6,5) * t245;
t160 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t247;
t159 = Icges(6,4) * t213 + Icges(6,2) * t212 + Icges(6,6) * t245;
t158 = Icges(6,5) * t215 + Icges(6,6) * t214 + Icges(6,3) * t247;
t157 = Icges(6,5) * t213 + Icges(6,6) * t212 + Icges(6,3) * t245;
t156 = rSges(7,1) * t211 + rSges(7,2) * t210 + rSges(7,3) * t247;
t155 = rSges(7,1) * t209 + rSges(7,2) * t208 + rSges(7,3) * t245;
t154 = Icges(7,1) * t211 + Icges(7,4) * t210 + Icges(7,5) * t247;
t153 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t245;
t152 = Icges(7,4) * t211 + Icges(7,2) * t210 + Icges(7,6) * t247;
t151 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t245;
t150 = Icges(7,5) * t211 + Icges(7,6) * t210 + Icges(7,3) * t247;
t149 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t245;
t148 = pkin(5) * t358 + pkin(11) * t247 + t248 * t363;
t147 = pkin(5) * t359 + pkin(11) * t245 + t246 * t363;
t146 = -t202 * t277 + t236 * t252 + t333;
t145 = t203 * t277 - t236 * t253 + t331;
t144 = t202 * t253 - t203 * t252 + t332;
t143 = -t190 * t261 + t219 * t221 + t330;
t142 = t191 * t261 - t219 * t222 + t328;
t141 = t190 * t222 - t191 * t221 + t329;
t140 = -t163 * t223 + t183 * t194 + t327;
t139 = t164 * t223 - t183 * t195 + t325;
t138 = t163 * t195 - t164 * t194 + t326;
t137 = -t147 * t223 - t155 * t207 + t167 * t177 + t173 * t194 + t327;
t136 = t148 * t223 + t156 * t207 - t168 * t177 - t173 * t195 + t325;
t135 = t147 * t195 - t148 * t194 + t155 * t168 - t156 * t167 + t326;
t1 = ((t293 * t315 + t295 * t313 + Icges(1,2)) * V_base(5) + (t294 * t315 + t296 * t313 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + m(3) * (t165 ^ 2 + t170 ^ 2 + t171 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + t300 * ((t225 * t288 + t226 * t289 + t263 * t300) * t316 + ((t228 * t322 + t230 * t319) * t289 + (t227 * t322 + t229 * t319) * t288 + (t264 * t322 + t265 * t319) * t300) * t314) / 0.2e1 + ((Icges(2,5) * t315 - Icges(2,6) * t313 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t313 + Icges(2,6) * t315 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + t167 * ((t150 * t245 + t152 * t208 + t154 * t209) * t168 + (t245 * t149 + t208 * t151 + t209 * t153) * t167 + (t174 * t245 + t175 * t208 + t176 * t209) * t207) / 0.2e1 + t194 * ((t158 * t245 + t160 * t212 + t162 * t213) * t195 + (t245 * t157 + t212 * t159 + t213 * t161) * t194 + (t180 * t245 + t181 * t212 + t182 * t213) * t223) / 0.2e1 + t168 * ((t247 * t150 + t210 * t152 + t211 * t154) * t168 + (t149 * t247 + t151 * t210 + t153 * t211) * t167 + (t174 * t247 + t175 * t210 + t176 * t211) * t207) / 0.2e1 + t195 * ((t247 * t158 + t214 * t160 + t215 * t162) * t195 + (t157 * t247 + t159 * t214 + t161 * t215) * t194 + (t180 * t247 + t181 * t214 + t182 * t215) * t223) / 0.2e1 + m(2) * (t251 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + t277 * ((-t197 * t352 + t199 * t280 + t201 * t281) * t253 + (-t196 * t352 + t198 * t280 + t200 * t281) * t252 + (-t233 * t352 + t280 * t234 + t281 * t235) * t277) / 0.2e1 + t207 * ((t150 * t267 + t152 * t239 + t154 * t240) * t168 + (t149 * t267 + t151 * t239 + t153 * t240) * t167 + (t267 * t174 + t239 * t175 + t240 * t176) * t207) / 0.2e1 + t261 * ((-t185 * t352 - t187 * t267 + t189 * t268) * t222 + (-t184 * t352 - t186 * t267 + t188 * t268) * t221 + (-t216 * t352 - t267 * t217 + t268 * t218) * t261) / 0.2e1 + t223 * ((t158 * t267 + t160 * t249 + t162 * t250) * t195 + (t157 * t267 + t159 * t249 + t161 * t250) * t194 + (t267 * t180 + t249 * t181 + t250 * t182) * t223) / 0.2e1 + t221 * ((t185 * t273 - t187 * t245 + t189 * t246) * t222 + (t184 * t273 - t245 * t186 + t246 * t188) * t221 + (t216 * t273 - t217 * t245 + t218 * t246) * t261) / 0.2e1 + t288 * ((-t226 * t356 - t228 * t273 + t230 * t274) * t289 + (-t225 * t356 - t273 * t227 + t274 * t229) * t288 + (-t263 * t356 - t264 * t273 + t265 * t274) * t300) / 0.2e1 + t289 * ((t226 * t357 - t275 * t228 + t276 * t230) * t289 + (t225 * t357 - t227 * t275 + t229 * t276) * t288 + (t263 * t357 - t264 * t275 + t265 * t276) * t300) / 0.2e1 + t222 * ((t275 * t185 - t247 * t187 + t248 * t189) * t222 + (t184 * t275 - t186 * t247 + t188 * t248) * t221 + (t216 * t275 - t217 * t247 + t218 * t248) * t261) / 0.2e1 + t252 * ((t197 * t273 + t199 * t254 + t201 * t255) * t253 + (t273 * t196 + t254 * t198 + t255 * t200) * t252 + (t233 * t273 + t234 * t254 + t235 * t255) * t277) / 0.2e1 + t253 * ((t275 * t197 + t256 * t199 + t257 * t201) * t253 + (t196 * t275 + t198 * t256 + t200 * t257) * t252 + (t233 * t275 + t234 * t256 + t235 * t257) * t277) / 0.2e1 + m(1) * (t285 ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1 + ((-t293 * t313 + t295 * t315 + Icges(1,4)) * V_base(5) + (-t313 * t294 + t315 * t296 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1;
T  = t1;
