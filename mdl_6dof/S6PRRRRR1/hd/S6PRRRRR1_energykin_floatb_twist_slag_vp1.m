% Calculate kinetic energy for
% S6PRRRRR1
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
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:37:27
% EndTime: 2019-03-09 00:37:31
% DurationCPUTime: 4.22s
% Computational Cost: add. (3737->444), mult. (5704->662), div. (0->0), fcn. (6640->14), ass. (0->195)
t322 = cos(pkin(6));
t369 = pkin(7) * t322;
t327 = cos(qJ(3));
t368 = t327 * pkin(3);
t319 = sin(pkin(12));
t366 = Icges(2,4) * t319;
t320 = sin(pkin(6));
t365 = t319 * t320;
t321 = cos(pkin(12));
t364 = t320 * t321;
t324 = sin(qJ(3));
t363 = t320 * t324;
t325 = sin(qJ(2));
t362 = t320 * t325;
t361 = t320 * t327;
t328 = cos(qJ(2));
t360 = t320 * t328;
t359 = t322 * t324;
t358 = t322 * t325;
t357 = t322 * t328;
t318 = qJ(3) + qJ(4);
t314 = cos(t318);
t356 = pkin(4) * t314;
t354 = qJD(2) * t320;
t353 = -qJD(3) - qJD(4);
t352 = V_base(5) * qJ(1) + V_base(1);
t348 = qJD(1) + V_base(3);
t347 = t319 * t363;
t346 = t321 * t363;
t345 = qJ(5) + t318;
t294 = t319 * t354 + V_base(4);
t306 = qJD(2) * t322 + V_base(6);
t313 = sin(t318);
t344 = pkin(4) * t313;
t278 = t319 * t357 + t321 * t325;
t253 = qJD(3) * t278 + t294;
t343 = cos(t345);
t223 = qJD(4) * t278 + t253;
t342 = t320 * t343;
t293 = -t321 * t354 + V_base(5);
t206 = qJD(5) * t278 + t223;
t276 = t319 * t325 - t321 * t357;
t252 = qJD(3) * t276 + t293;
t286 = pkin(1) * t319 - pkin(7) * t364;
t341 = -t286 * V_base(6) + V_base(5) * t369 + t352;
t287 = pkin(1) * t321 + pkin(7) * t365;
t340 = V_base(4) * t286 - t287 * V_base(5) + t348;
t222 = qJD(4) * t276 + t252;
t205 = qJD(5) * t276 + t222;
t339 = V_base(6) * t287 + V_base(2) + (-qJ(1) - t369) * V_base(4);
t250 = (-qJD(5) + t353) * t360 + t306;
t277 = t319 * t328 + t321 * t358;
t238 = t277 * pkin(2) + t276 * pkin(8);
t285 = (pkin(2) * t325 - pkin(8) * t328) * t320;
t338 = -t238 * t306 + t293 * t285 + t341;
t279 = -t319 * t358 + t321 * t328;
t239 = t279 * pkin(2) + t278 * pkin(8);
t337 = t294 * t238 - t239 * t293 + t340;
t336 = t306 * t239 - t285 * t294 + t339;
t191 = -pkin(3) * t346 + pkin(9) * t276 + t277 * t368;
t237 = pkin(3) * t359 + (-pkin(9) * t328 + t325 * t368) * t320;
t280 = -qJD(3) * t360 + t306;
t335 = -t191 * t280 + t252 * t237 + t338;
t192 = pkin(3) * t347 + pkin(9) * t278 + t279 * t368;
t334 = t253 * t191 - t192 * t252 + t337;
t333 = t280 * t192 - t237 * t253 + t336;
t163 = pkin(10) * t276 + t277 * t356 - t344 * t364;
t204 = t344 * t322 + (-pkin(10) * t328 + t325 * t356) * t320;
t261 = t353 * t360 + t306;
t332 = -t163 * t261 + t222 * t204 + t335;
t164 = pkin(10) * t278 + t279 * t356 + t344 * t365;
t331 = t223 * t163 - t164 * t222 + t334;
t330 = t261 * t164 - t204 * t223 + t333;
t326 = cos(qJ(6));
t323 = sin(qJ(6));
t312 = Icges(2,4) * t321;
t309 = sin(t345);
t303 = rSges(2,1) * t321 - rSges(2,2) * t319;
t302 = rSges(2,1) * t319 + rSges(2,2) * t321;
t301 = Icges(2,1) * t321 - t366;
t300 = Icges(2,1) * t319 + t312;
t299 = -Icges(2,2) * t319 + t312;
t298 = Icges(2,2) * t321 + t366;
t292 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t291 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t290 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t284 = t325 * t361 + t359;
t283 = t322 * t327 - t324 * t362;
t269 = t313 * t322 + t314 * t362;
t268 = -t313 * t362 + t314 * t322;
t267 = rSges(3,3) * t322 + (rSges(3,1) * t325 + rSges(3,2) * t328) * t320;
t266 = Icges(3,5) * t322 + (Icges(3,1) * t325 + Icges(3,4) * t328) * t320;
t265 = Icges(3,6) * t322 + (Icges(3,4) * t325 + Icges(3,2) * t328) * t320;
t264 = Icges(3,3) * t322 + (Icges(3,5) * t325 + Icges(3,6) * t328) * t320;
t263 = t322 * t309 + t325 * t342;
t262 = t309 * t362 - t322 * t343;
t260 = V_base(5) * rSges(2,3) - t302 * V_base(6) + t352;
t259 = t303 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t257 = t279 * t327 + t347;
t256 = -t279 * t324 + t319 * t361;
t255 = t277 * t327 - t346;
t254 = -t277 * t324 - t321 * t361;
t251 = t302 * V_base(4) - t303 * V_base(5) + t348;
t249 = t279 * t314 + t313 * t365;
t248 = -t279 * t313 + t314 * t365;
t247 = t277 * t314 - t313 * t364;
t246 = -t277 * t313 - t314 * t364;
t245 = t263 * t326 - t323 * t360;
t244 = -t263 * t323 - t326 * t360;
t243 = t279 * t343 + t309 * t365;
t242 = t279 * t309 - t319 * t342;
t241 = t277 * t343 - t309 * t364;
t240 = t277 * t309 + t321 * t342;
t236 = rSges(4,1) * t284 + rSges(4,2) * t283 - rSges(4,3) * t360;
t235 = Icges(4,1) * t284 + Icges(4,4) * t283 - Icges(4,5) * t360;
t234 = Icges(4,4) * t284 + Icges(4,2) * t283 - Icges(4,6) * t360;
t233 = Icges(4,5) * t284 + Icges(4,6) * t283 - Icges(4,3) * t360;
t232 = rSges(3,1) * t279 - rSges(3,2) * t278 + rSges(3,3) * t365;
t231 = rSges(3,1) * t277 - rSges(3,2) * t276 - rSges(3,3) * t364;
t230 = Icges(3,1) * t279 - Icges(3,4) * t278 + Icges(3,5) * t365;
t229 = Icges(3,1) * t277 - Icges(3,4) * t276 - Icges(3,5) * t364;
t228 = Icges(3,4) * t279 - Icges(3,2) * t278 + Icges(3,6) * t365;
t227 = Icges(3,4) * t277 - Icges(3,2) * t276 - Icges(3,6) * t364;
t226 = Icges(3,5) * t279 - Icges(3,6) * t278 + Icges(3,3) * t365;
t225 = Icges(3,5) * t277 - Icges(3,6) * t276 - Icges(3,3) * t364;
t221 = pkin(5) * t263 + pkin(11) * t262;
t219 = rSges(5,1) * t269 + rSges(5,2) * t268 - rSges(5,3) * t360;
t218 = Icges(5,1) * t269 + Icges(5,4) * t268 - Icges(5,5) * t360;
t217 = Icges(5,4) * t269 + Icges(5,2) * t268 - Icges(5,6) * t360;
t216 = Icges(5,5) * t269 + Icges(5,6) * t268 - Icges(5,3) * t360;
t215 = qJD(6) * t262 + t250;
t214 = rSges(6,1) * t263 - rSges(6,2) * t262 - rSges(6,3) * t360;
t213 = Icges(6,1) * t263 - Icges(6,4) * t262 - Icges(6,5) * t360;
t212 = Icges(6,4) * t263 - Icges(6,2) * t262 - Icges(6,6) * t360;
t211 = Icges(6,5) * t263 - Icges(6,6) * t262 - Icges(6,3) * t360;
t210 = t243 * t326 + t278 * t323;
t209 = -t243 * t323 + t278 * t326;
t208 = t241 * t326 + t276 * t323;
t207 = -t241 * t323 + t276 * t326;
t202 = pkin(5) * t243 + pkin(11) * t242;
t201 = pkin(5) * t241 + pkin(11) * t240;
t200 = rSges(4,1) * t257 + rSges(4,2) * t256 + rSges(4,3) * t278;
t199 = rSges(4,1) * t255 + rSges(4,2) * t254 + rSges(4,3) * t276;
t198 = Icges(4,1) * t257 + Icges(4,4) * t256 + Icges(4,5) * t278;
t197 = Icges(4,1) * t255 + Icges(4,4) * t254 + Icges(4,5) * t276;
t196 = Icges(4,4) * t257 + Icges(4,2) * t256 + Icges(4,6) * t278;
t195 = Icges(4,4) * t255 + Icges(4,2) * t254 + Icges(4,6) * t276;
t194 = Icges(4,5) * t257 + Icges(4,6) * t256 + Icges(4,3) * t278;
t193 = Icges(4,5) * t255 + Icges(4,6) * t254 + Icges(4,3) * t276;
t190 = rSges(5,1) * t249 + rSges(5,2) * t248 + rSges(5,3) * t278;
t189 = rSges(5,1) * t247 + rSges(5,2) * t246 + rSges(5,3) * t276;
t188 = Icges(5,1) * t249 + Icges(5,4) * t248 + Icges(5,5) * t278;
t187 = Icges(5,1) * t247 + Icges(5,4) * t246 + Icges(5,5) * t276;
t186 = Icges(5,4) * t249 + Icges(5,2) * t248 + Icges(5,6) * t278;
t185 = Icges(5,4) * t247 + Icges(5,2) * t246 + Icges(5,6) * t276;
t184 = Icges(5,5) * t249 + Icges(5,6) * t248 + Icges(5,3) * t278;
t183 = Icges(5,5) * t247 + Icges(5,6) * t246 + Icges(5,3) * t276;
t182 = rSges(6,1) * t243 - rSges(6,2) * t242 + rSges(6,3) * t278;
t181 = rSges(6,1) * t241 - rSges(6,2) * t240 + rSges(6,3) * t276;
t180 = Icges(6,1) * t243 - Icges(6,4) * t242 + Icges(6,5) * t278;
t179 = Icges(6,1) * t241 - Icges(6,4) * t240 + Icges(6,5) * t276;
t178 = Icges(6,4) * t243 - Icges(6,2) * t242 + Icges(6,6) * t278;
t177 = Icges(6,4) * t241 - Icges(6,2) * t240 + Icges(6,6) * t276;
t176 = Icges(6,5) * t243 - Icges(6,6) * t242 + Icges(6,3) * t278;
t175 = Icges(6,5) * t241 - Icges(6,6) * t240 + Icges(6,3) * t276;
t174 = rSges(7,1) * t245 + rSges(7,2) * t244 + rSges(7,3) * t262;
t173 = Icges(7,1) * t245 + Icges(7,4) * t244 + Icges(7,5) * t262;
t172 = Icges(7,4) * t245 + Icges(7,2) * t244 + Icges(7,6) * t262;
t171 = Icges(7,5) * t245 + Icges(7,6) * t244 + Icges(7,3) * t262;
t170 = qJD(6) * t242 + t206;
t169 = qJD(6) * t240 + t205;
t167 = -t231 * t306 + t267 * t293 + t341;
t166 = t232 * t306 - t267 * t294 + t339;
t160 = t231 * t294 - t232 * t293 + t340;
t159 = rSges(7,1) * t210 + rSges(7,2) * t209 + rSges(7,3) * t242;
t158 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t240;
t157 = Icges(7,1) * t210 + Icges(7,4) * t209 + Icges(7,5) * t242;
t156 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t240;
t155 = Icges(7,4) * t210 + Icges(7,2) * t209 + Icges(7,6) * t242;
t154 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t240;
t153 = Icges(7,5) * t210 + Icges(7,6) * t209 + Icges(7,3) * t242;
t152 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t240;
t150 = -t199 * t280 + t236 * t252 + t338;
t149 = t200 * t280 - t236 * t253 + t336;
t148 = t199 * t253 - t200 * t252 + t337;
t147 = -t189 * t261 + t219 * t222 + t335;
t146 = t190 * t261 - t219 * t223 + t333;
t145 = t189 * t223 - t190 * t222 + t334;
t144 = -t181 * t250 + t205 * t214 + t332;
t143 = t182 * t250 - t206 * t214 + t330;
t142 = t181 * t206 - t182 * t205 + t331;
t141 = -t158 * t215 + t169 * t174 - t201 * t250 + t205 * t221 + t332;
t140 = t159 * t215 - t170 * t174 + t202 * t250 - t206 * t221 + t330;
t139 = t158 * t170 - t159 * t169 + t201 * t206 - t202 * t205 + t331;
t1 = ((-t298 * t319 + t300 * t321 + Icges(1,4)) * V_base(5) + (-t299 * t319 + t301 * t321 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t298 * t321 + t300 * t319 + Icges(1,2)) * V_base(5) + (t299 * t321 + t301 * t319 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t306 * ((t225 * t293 + t226 * t294 + t264 * t306) * t322 + ((t228 * t328 + t230 * t325) * t294 + (t227 * t328 + t229 * t325) * t293 + (t265 * t328 + t266 * t325) * t306) * t320) / 0.2e1 + m(5) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(4) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(6) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(7) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + ((Icges(2,5) * t321 - Icges(2,6) * t319 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t319 + Icges(2,6) * t321 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + m(3) * (t160 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + t169 * ((t153 * t240 + t155 * t207 + t157 * t208) * t170 + (t240 * t152 + t207 * t154 + t208 * t156) * t169 + (t171 * t240 + t172 * t207 + t173 * t208) * t215) / 0.2e1 + t170 * ((t242 * t153 + t209 * t155 + t210 * t157) * t170 + (t152 * t242 + t154 * t209 + t156 * t210) * t169 + (t171 * t242 + t172 * t209 + t173 * t210) * t215) / 0.2e1 + m(2) * (t251 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + t215 * ((t153 * t262 + t155 * t244 + t157 * t245) * t170 + (t152 * t262 + t154 * t244 + t156 * t245) * t169 + (t262 * t171 + t244 * t172 + t245 * t173) * t215) / 0.2e1 + t205 * ((t176 * t276 - t178 * t240 + t180 * t241) * t206 + (t276 * t175 - t240 * t177 + t241 * t179) * t205 + (t211 * t276 - t212 * t240 + t213 * t241) * t250) / 0.2e1 + t222 * ((t184 * t276 + t186 * t246 + t188 * t247) * t223 + (t276 * t183 + t246 * t185 + t247 * t187) * t222 + (t216 * t276 + t217 * t246 + t218 * t247) * t261) / 0.2e1 + t206 * ((t278 * t176 - t242 * t178 + t243 * t180) * t206 + (t175 * t278 - t177 * t242 + t179 * t243) * t205 + (t211 * t278 - t212 * t242 + t213 * t243) * t250) / 0.2e1 + t223 * ((t278 * t184 + t248 * t186 + t249 * t188) * t223 + (t183 * t278 + t185 * t248 + t187 * t249) * t222 + (t216 * t278 + t217 * t248 + t218 * t249) * t261) / 0.2e1 + t280 * ((-t194 * t360 + t196 * t283 + t198 * t284) * t253 + (-t193 * t360 + t195 * t283 + t197 * t284) * t252 + (-t233 * t360 + t234 * t283 + t235 * t284) * t280) / 0.2e1 + t252 * ((t194 * t276 + t196 * t254 + t198 * t255) * t253 + (t193 * t276 + t195 * t254 + t197 * t255) * t252 + (t233 * t276 + t234 * t254 + t235 * t255) * t280) / 0.2e1 + t253 * ((t194 * t278 + t196 * t256 + t198 * t257) * t253 + (t193 * t278 + t195 * t256 + t197 * t257) * t252 + (t233 * t278 + t234 * t256 + t235 * t257) * t280) / 0.2e1 + t261 * ((-t184 * t360 + t186 * t268 + t188 * t269) * t223 + (-t183 * t360 + t185 * t268 + t187 * t269) * t222 + (-t216 * t360 + t217 * t268 + t218 * t269) * t261) / 0.2e1 + t250 * ((-t176 * t360 - t178 * t262 + t180 * t263) * t206 + (-t175 * t360 - t177 * t262 + t179 * t263) * t205 + (-t211 * t360 - t212 * t262 + t213 * t263) * t250) / 0.2e1 + m(1) * (t290 ^ 2 + t291 ^ 2 + t292 ^ 2) / 0.2e1 + t293 * ((-t226 * t364 - t228 * t276 + t230 * t277) * t294 + (-t225 * t364 - t227 * t276 + t229 * t277) * t293 + (-t264 * t364 - t265 * t276 + t266 * t277) * t306) / 0.2e1 + t294 * ((t226 * t365 - t228 * t278 + t230 * t279) * t294 + (t225 * t365 - t227 * t278 + t229 * t279) * t293 + (t264 * t365 - t265 * t278 + t266 * t279) * t306) / 0.2e1;
T  = t1;
