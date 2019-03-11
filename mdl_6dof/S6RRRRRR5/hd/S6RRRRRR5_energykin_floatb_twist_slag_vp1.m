% Calculate kinetic energy for
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:56:33
% EndTime: 2019-03-10 03:56:37
% DurationCPUTime: 4.07s
% Computational Cost: add. (3797->444), mult. (5704->663), div. (0->0), fcn. (6640->14), ass. (0->196)
t321 = cos(pkin(6));
t369 = pkin(8) * t321;
t327 = cos(qJ(3));
t368 = t327 * pkin(3);
t325 = sin(qJ(1));
t366 = Icges(2,4) * t325;
t320 = sin(pkin(6));
t324 = sin(qJ(2));
t365 = t320 * t324;
t364 = t320 * t325;
t363 = t320 * t327;
t328 = cos(qJ(2));
t362 = t320 * t328;
t329 = cos(qJ(1));
t361 = t320 * t329;
t323 = sin(qJ(3));
t360 = t321 * t323;
t359 = t324 * t325;
t358 = t324 * t329;
t357 = t325 * t328;
t356 = t328 * t329;
t319 = qJ(3) + qJ(4);
t314 = cos(t319);
t355 = pkin(4) * t314;
t353 = qJD(2) * t320;
t352 = -qJD(3) - qJD(4);
t351 = V_base(5) * pkin(7) + V_base(1);
t348 = t323 * t364;
t347 = t323 * t361;
t346 = qJ(5) + t319;
t295 = t325 * t353 + V_base(4);
t312 = V_base(6) + qJD(1);
t313 = sin(t319);
t345 = pkin(4) * t313;
t282 = t321 * t357 + t358;
t252 = qJD(3) * t282 + t295;
t344 = cos(t346);
t296 = qJD(2) * t321 + t312;
t223 = qJD(4) * t282 + t252;
t343 = t320 * t344;
t294 = -t329 * t353 + V_base(5);
t206 = qJD(5) * t282 + t223;
t286 = pkin(1) * t325 - pkin(8) * t361;
t342 = -t286 * t312 + V_base(5) * t369 + t351;
t287 = pkin(1) * t329 + pkin(8) * t364;
t341 = V_base(4) * t286 - t287 * V_base(5) + V_base(3);
t280 = -t321 * t356 + t359;
t251 = qJD(3) * t280 + t294;
t222 = qJD(4) * t280 + t251;
t340 = t312 * t287 + V_base(2) + (-pkin(7) - t369) * V_base(4);
t205 = qJD(5) * t280 + t222;
t281 = t321 * t358 + t357;
t240 = t281 * pkin(2) + t280 * pkin(9);
t285 = (pkin(2) * t324 - pkin(9) * t328) * t320;
t339 = -t240 * t296 + t294 * t285 + t342;
t283 = -t321 * t359 + t356;
t241 = t283 * pkin(2) + t282 * pkin(9);
t338 = t295 * t240 - t241 * t294 + t341;
t246 = (-qJD(5) + t352) * t362 + t296;
t337 = t296 * t241 - t285 * t295 + t340;
t191 = -pkin(3) * t347 + pkin(10) * t280 + t281 * t368;
t235 = pkin(3) * t360 + (-pkin(10) * t328 + t324 * t368) * t320;
t276 = -qJD(3) * t362 + t296;
t336 = -t191 * t276 + t251 * t235 + t339;
t192 = pkin(3) * t348 + pkin(10) * t282 + t283 * t368;
t335 = t252 * t191 - t192 * t251 + t338;
t334 = t276 * t192 - t235 * t252 + t337;
t165 = pkin(11) * t280 + t281 * t355 - t345 * t361;
t204 = t345 * t321 + (-pkin(11) * t328 + t324 * t355) * t320;
t261 = t352 * t362 + t296;
t333 = -t165 * t261 + t222 * t204 + t336;
t166 = pkin(11) * t282 + t283 * t355 + t345 * t364;
t332 = t223 * t165 - t166 * t222 + t335;
t331 = t261 * t166 - t204 * t223 + t334;
t326 = cos(qJ(6));
t322 = sin(qJ(6));
t315 = Icges(2,4) * t329;
t309 = sin(t346);
t304 = rSges(2,1) * t329 - rSges(2,2) * t325;
t303 = rSges(2,1) * t325 + rSges(2,2) * t329;
t302 = Icges(2,1) * t329 - t366;
t301 = Icges(2,1) * t325 + t315;
t300 = -Icges(2,2) * t325 + t315;
t299 = Icges(2,2) * t329 + t366;
t292 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t291 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t290 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t279 = t324 * t363 + t360;
t278 = t321 * t327 - t323 * t365;
t269 = t313 * t321 + t314 * t365;
t268 = -t313 * t365 + t314 * t321;
t267 = rSges(3,3) * t321 + (rSges(3,1) * t324 + rSges(3,2) * t328) * t320;
t266 = Icges(3,5) * t321 + (Icges(3,1) * t324 + Icges(3,4) * t328) * t320;
t265 = Icges(3,6) * t321 + (Icges(3,4) * t324 + Icges(3,2) * t328) * t320;
t264 = Icges(3,3) * t321 + (Icges(3,5) * t324 + Icges(3,6) * t328) * t320;
t263 = t321 * t309 + t324 * t343;
t262 = t309 * t365 - t321 * t344;
t260 = V_base(5) * rSges(2,3) - t303 * t312 + t351;
t259 = t304 * t312 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t257 = t303 * V_base(4) - t304 * V_base(5) + V_base(3);
t256 = t283 * t327 + t348;
t255 = -t283 * t323 + t325 * t363;
t254 = t281 * t327 - t347;
t253 = -t281 * t323 - t327 * t361;
t250 = t283 * t314 + t313 * t364;
t249 = -t283 * t313 + t314 * t364;
t248 = t281 * t314 - t313 * t361;
t247 = -t281 * t313 - t314 * t361;
t245 = t283 * t344 + t309 * t364;
t244 = t283 * t309 - t325 * t343;
t243 = t281 * t344 - t309 * t361;
t242 = t281 * t309 + t329 * t343;
t239 = t263 * t326 - t322 * t362;
t238 = -t263 * t322 - t326 * t362;
t237 = rSges(3,1) * t283 - rSges(3,2) * t282 + rSges(3,3) * t364;
t236 = rSges(3,1) * t281 - rSges(3,2) * t280 - rSges(3,3) * t361;
t234 = Icges(3,1) * t283 - Icges(3,4) * t282 + Icges(3,5) * t364;
t233 = Icges(3,1) * t281 - Icges(3,4) * t280 - Icges(3,5) * t361;
t232 = Icges(3,4) * t283 - Icges(3,2) * t282 + Icges(3,6) * t364;
t231 = Icges(3,4) * t281 - Icges(3,2) * t280 - Icges(3,6) * t361;
t230 = Icges(3,5) * t283 - Icges(3,6) * t282 + Icges(3,3) * t364;
t229 = Icges(3,5) * t281 - Icges(3,6) * t280 - Icges(3,3) * t361;
t228 = rSges(4,1) * t279 + rSges(4,2) * t278 - rSges(4,3) * t362;
t227 = Icges(4,1) * t279 + Icges(4,4) * t278 - Icges(4,5) * t362;
t226 = Icges(4,4) * t279 + Icges(4,2) * t278 - Icges(4,6) * t362;
t225 = Icges(4,5) * t279 + Icges(4,6) * t278 - Icges(4,3) * t362;
t220 = pkin(5) * t263 + pkin(12) * t262;
t219 = rSges(5,1) * t269 + rSges(5,2) * t268 - rSges(5,3) * t362;
t218 = Icges(5,1) * t269 + Icges(5,4) * t268 - Icges(5,5) * t362;
t217 = Icges(5,4) * t269 + Icges(5,2) * t268 - Icges(5,6) * t362;
t216 = Icges(5,5) * t269 + Icges(5,6) * t268 - Icges(5,3) * t362;
t215 = rSges(6,1) * t263 - rSges(6,2) * t262 - rSges(6,3) * t362;
t214 = Icges(6,1) * t263 - Icges(6,4) * t262 - Icges(6,5) * t362;
t213 = Icges(6,4) * t263 - Icges(6,2) * t262 - Icges(6,6) * t362;
t212 = Icges(6,5) * t263 - Icges(6,6) * t262 - Icges(6,3) * t362;
t211 = qJD(6) * t262 + t246;
t210 = t245 * t326 + t282 * t322;
t209 = -t245 * t322 + t282 * t326;
t208 = t243 * t326 + t280 * t322;
t207 = -t243 * t322 + t280 * t326;
t202 = pkin(5) * t245 + pkin(12) * t244;
t201 = pkin(5) * t243 + pkin(12) * t242;
t200 = rSges(4,1) * t256 + rSges(4,2) * t255 + rSges(4,3) * t282;
t199 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t280;
t198 = Icges(4,1) * t256 + Icges(4,4) * t255 + Icges(4,5) * t282;
t197 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t280;
t196 = Icges(4,4) * t256 + Icges(4,2) * t255 + Icges(4,6) * t282;
t195 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t280;
t194 = Icges(4,5) * t256 + Icges(4,6) * t255 + Icges(4,3) * t282;
t193 = Icges(4,5) * t254 + Icges(4,6) * t253 + Icges(4,3) * t280;
t190 = rSges(5,1) * t250 + rSges(5,2) * t249 + rSges(5,3) * t282;
t189 = rSges(5,1) * t248 + rSges(5,2) * t247 + rSges(5,3) * t280;
t188 = Icges(5,1) * t250 + Icges(5,4) * t249 + Icges(5,5) * t282;
t187 = Icges(5,1) * t248 + Icges(5,4) * t247 + Icges(5,5) * t280;
t186 = Icges(5,4) * t250 + Icges(5,2) * t249 + Icges(5,6) * t282;
t185 = Icges(5,4) * t248 + Icges(5,2) * t247 + Icges(5,6) * t280;
t184 = Icges(5,5) * t250 + Icges(5,6) * t249 + Icges(5,3) * t282;
t183 = Icges(5,5) * t248 + Icges(5,6) * t247 + Icges(5,3) * t280;
t182 = rSges(6,1) * t245 - rSges(6,2) * t244 + rSges(6,3) * t282;
t181 = rSges(6,1) * t243 - rSges(6,2) * t242 + rSges(6,3) * t280;
t180 = Icges(6,1) * t245 - Icges(6,4) * t244 + Icges(6,5) * t282;
t179 = Icges(6,1) * t243 - Icges(6,4) * t242 + Icges(6,5) * t280;
t178 = Icges(6,4) * t245 - Icges(6,2) * t244 + Icges(6,6) * t282;
t177 = Icges(6,4) * t243 - Icges(6,2) * t242 + Icges(6,6) * t280;
t176 = Icges(6,5) * t245 - Icges(6,6) * t244 + Icges(6,3) * t282;
t175 = Icges(6,5) * t243 - Icges(6,6) * t242 + Icges(6,3) * t280;
t174 = rSges(7,1) * t239 + rSges(7,2) * t238 + rSges(7,3) * t262;
t173 = Icges(7,1) * t239 + Icges(7,4) * t238 + Icges(7,5) * t262;
t172 = Icges(7,4) * t239 + Icges(7,2) * t238 + Icges(7,6) * t262;
t171 = Icges(7,5) * t239 + Icges(7,6) * t238 + Icges(7,3) * t262;
t170 = qJD(6) * t244 + t206;
t169 = qJD(6) * t242 + t205;
t164 = -t236 * t296 + t267 * t294 + t342;
t163 = t237 * t296 - t267 * t295 + t340;
t160 = t236 * t295 - t237 * t294 + t341;
t159 = rSges(7,1) * t210 + rSges(7,2) * t209 + rSges(7,3) * t244;
t158 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t242;
t157 = Icges(7,1) * t210 + Icges(7,4) * t209 + Icges(7,5) * t244;
t156 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t242;
t155 = Icges(7,4) * t210 + Icges(7,2) * t209 + Icges(7,6) * t244;
t154 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t242;
t153 = Icges(7,5) * t210 + Icges(7,6) * t209 + Icges(7,3) * t244;
t152 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t242;
t150 = -t199 * t276 + t228 * t251 + t339;
t149 = t200 * t276 - t228 * t252 + t337;
t148 = t199 * t252 - t200 * t251 + t338;
t147 = -t189 * t261 + t219 * t222 + t336;
t146 = t190 * t261 - t219 * t223 + t334;
t145 = t189 * t223 - t190 * t222 + t335;
t144 = -t181 * t246 + t205 * t215 + t333;
t143 = t182 * t246 - t206 * t215 + t331;
t142 = t181 * t206 - t182 * t205 + t332;
t141 = -t158 * t211 + t169 * t174 - t201 * t246 + t205 * t220 + t333;
t140 = t159 * t211 - t170 * t174 + t202 * t246 - t206 * t220 + t331;
t139 = t158 * t170 - t159 * t169 + t201 * t206 - t202 * t205 + t332;
t1 = ((t329 * t299 + t325 * t301 + Icges(1,2)) * V_base(5) + (t300 * t329 + t302 * t325 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + m(4) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(5) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(6) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(7) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + ((Icges(2,5) * t325 + Icges(2,6) * t329) * V_base(5) + (Icges(2,5) * t329 - Icges(2,6) * t325) * V_base(4) + Icges(2,3) * t312 / 0.2e1) * t312 + t296 * ((t229 * t294 + t230 * t295 + t264 * t296) * t321 + ((t232 * t328 + t234 * t324) * t295 + (t231 * t328 + t233 * t324) * t294 + (t265 * t328 + t266 * t324) * t296) * t320) / 0.2e1 + ((-t299 * t325 + t301 * t329 + Icges(1,4)) * V_base(5) + (-t325 * t300 + t329 * t302 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + m(3) * (t160 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + t294 * ((-t230 * t361 - t232 * t280 + t234 * t281) * t295 + (-t229 * t361 - t280 * t231 + t281 * t233) * t294 + (-t264 * t361 - t265 * t280 + t266 * t281) * t296) / 0.2e1 + t276 * ((-t194 * t362 + t196 * t278 + t198 * t279) * t252 + (-t193 * t362 + t195 * t278 + t197 * t279) * t251 + (-t225 * t362 + t278 * t226 + t279 * t227) * t276) / 0.2e1 + t261 * ((-t184 * t362 + t186 * t268 + t188 * t269) * t223 + (-t183 * t362 + t185 * t268 + t187 * t269) * t222 + (-t216 * t362 + t268 * t217 + t269 * t218) * t261) / 0.2e1 + t246 * ((-t176 * t362 - t178 * t262 + t180 * t263) * t206 + (-t175 * t362 - t177 * t262 + t179 * t263) * t205 + (-t212 * t362 - t262 * t213 + t263 * t214) * t246) / 0.2e1 + t295 * ((t230 * t364 - t282 * t232 + t283 * t234) * t295 + (t229 * t364 - t231 * t282 + t233 * t283) * t294 + (t264 * t364 - t265 * t282 + t266 * t283) * t296) / 0.2e1 + t169 * ((t153 * t242 + t155 * t207 + t157 * t208) * t170 + (t242 * t152 + t207 * t154 + t208 * t156) * t169 + (t171 * t242 + t172 * t207 + t173 * t208) * t211) / 0.2e1 + t170 * ((t244 * t153 + t209 * t155 + t210 * t157) * t170 + (t152 * t244 + t154 * t209 + t156 * t210) * t169 + (t171 * t244 + t172 * t209 + t173 * t210) * t211) / 0.2e1 + m(2) * (t257 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + t211 * ((t153 * t262 + t155 * t238 + t157 * t239) * t170 + (t152 * t262 + t154 * t238 + t156 * t239) * t169 + (t262 * t171 + t238 * t172 + t239 * t173) * t211) / 0.2e1 + t205 * ((t176 * t280 - t178 * t242 + t180 * t243) * t206 + (t280 * t175 - t242 * t177 + t243 * t179) * t205 + (t212 * t280 - t213 * t242 + t214 * t243) * t246) / 0.2e1 + t222 * ((t184 * t280 + t186 * t247 + t188 * t248) * t223 + (t280 * t183 + t247 * t185 + t248 * t187) * t222 + (t216 * t280 + t217 * t247 + t218 * t248) * t261) / 0.2e1 + t251 * ((t194 * t280 + t196 * t253 + t198 * t254) * t252 + (t280 * t193 + t253 * t195 + t254 * t197) * t251 + (t225 * t280 + t226 * t253 + t227 * t254) * t276) / 0.2e1 + t206 * ((t282 * t176 - t244 * t178 + t245 * t180) * t206 + (t175 * t282 - t177 * t244 + t179 * t245) * t205 + (t212 * t282 - t213 * t244 + t214 * t245) * t246) / 0.2e1 + t223 * ((t282 * t184 + t249 * t186 + t250 * t188) * t223 + (t183 * t282 + t185 * t249 + t187 * t250) * t222 + (t216 * t282 + t217 * t249 + t218 * t250) * t261) / 0.2e1 + t252 * ((t282 * t194 + t255 * t196 + t256 * t198) * t252 + (t193 * t282 + t195 * t255 + t197 * t256) * t251 + (t225 * t282 + t226 * t255 + t227 * t256) * t276) / 0.2e1 + m(1) * (t290 ^ 2 + t291 ^ 2 + t292 ^ 2) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
