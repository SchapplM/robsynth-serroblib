% Calculate kinetic energy for
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:33
% EndTime: 2019-03-08 23:15:37
% DurationCPUTime: 4.71s
% Computational Cost: add. (3483->437), mult. (6882->628), div. (0->0), fcn. (8332->14), ass. (0->187)
t367 = Icges(5,3) + Icges(6,3);
t312 = sin(pkin(11));
t314 = cos(pkin(11));
t321 = cos(qJ(2));
t315 = cos(pkin(6));
t319 = sin(qJ(2));
t346 = t315 * t319;
t268 = t312 * t321 + t314 * t346;
t313 = sin(pkin(6));
t318 = sin(qJ(3));
t348 = t313 * t318;
t358 = cos(qJ(3));
t251 = t268 * t358 - t314 * t348;
t345 = t315 * t321;
t267 = t312 * t319 - t314 * t345;
t311 = qJ(4) + pkin(12);
t304 = sin(t311);
t305 = cos(t311);
t212 = -t251 * t304 + t267 * t305;
t213 = t251 * t305 + t267 * t304;
t317 = sin(qJ(4));
t320 = cos(qJ(4));
t216 = -t251 * t317 + t267 * t320;
t352 = t267 * t317;
t217 = t251 * t320 + t352;
t335 = t313 * t358;
t250 = t268 * t318 + t314 * t335;
t365 = Icges(5,5) * t217 + Icges(6,5) * t213 + Icges(5,6) * t216 + Icges(6,6) * t212 + t367 * t250;
t270 = -t312 * t346 + t314 * t321;
t253 = t270 * t358 + t312 * t348;
t269 = t312 * t345 + t314 * t319;
t214 = -t253 * t304 + t269 * t305;
t215 = t253 * t305 + t269 * t304;
t218 = -t253 * t317 + t269 * t320;
t351 = t269 * t317;
t219 = t253 * t320 + t351;
t252 = t270 * t318 - t312 * t335;
t364 = Icges(5,5) * t219 + Icges(6,5) * t215 + Icges(5,6) * t218 + Icges(6,6) * t214 + t367 * t252;
t275 = t315 * t318 + t319 * t335;
t347 = t313 * t321;
t241 = -t275 * t304 - t305 * t347;
t242 = t275 * t305 - t304 * t347;
t254 = -t275 * t317 - t320 * t347;
t336 = t317 * t347;
t255 = t275 * t320 - t336;
t274 = -t315 * t358 + t319 * t348;
t363 = Icges(5,5) * t255 + Icges(6,5) * t242 + Icges(5,6) * t254 + Icges(6,6) * t241 + t367 * t274;
t356 = pkin(7) * t315;
t355 = t320 * pkin(4);
t353 = Icges(2,4) * t312;
t350 = t312 * t313;
t349 = t313 * t314;
t344 = pkin(5) * t305;
t342 = qJD(2) * t313;
t341 = V_base(5) * qJ(1) + V_base(1);
t337 = qJD(1) + V_base(3);
t285 = t312 * t342 + V_base(4);
t296 = qJD(2) * t315 + V_base(6);
t334 = pkin(5) * t304;
t249 = qJD(3) * t269 + t285;
t207 = qJD(4) * t252 + t249;
t284 = -t314 * t342 + V_base(5);
t248 = qJD(3) * t267 + t284;
t271 = -qJD(3) * t347 + t296;
t277 = pkin(1) * t312 - pkin(7) * t349;
t333 = -t277 * V_base(6) + V_base(5) * t356 + t341;
t278 = pkin(1) * t314 + pkin(7) * t350;
t332 = V_base(4) * t277 - t278 * V_base(5) + t337;
t206 = qJD(4) * t250 + t248;
t240 = qJD(4) * t274 + t271;
t331 = V_base(6) * t278 + V_base(2) + (-qJ(1) - t356) * V_base(4);
t235 = pkin(2) * t268 + pkin(8) * t267;
t276 = (pkin(2) * t319 - pkin(8) * t321) * t313;
t330 = -t235 * t296 + t284 * t276 + t333;
t236 = pkin(2) * t270 + pkin(8) * t269;
t329 = t285 * t235 - t236 * t284 + t332;
t328 = t296 * t236 - t276 * t285 + t331;
t204 = pkin(3) * t251 + pkin(9) * t250;
t239 = t275 * pkin(3) + t274 * pkin(9);
t327 = -t204 * t271 + t248 * t239 + t330;
t205 = pkin(3) * t253 + pkin(9) * t252;
t326 = t249 * t204 - t205 * t248 + t329;
t325 = t271 * t205 - t239 * t249 + t328;
t189 = -pkin(4) * t336 + qJ(5) * t274 + t275 * t355;
t324 = qJD(5) * t252 + t206 * t189 + t327;
t155 = pkin(4) * t352 + qJ(5) * t250 + t251 * t355;
t323 = qJD(5) * t274 + t207 * t155 + t326;
t156 = pkin(4) * t351 + qJ(5) * t252 + t253 * t355;
t322 = qJD(5) * t250 + t240 * t156 + t325;
t307 = qJ(6) + t311;
t306 = Icges(2,4) * t314;
t301 = cos(t307);
t300 = sin(t307);
t293 = rSges(2,1) * t314 - rSges(2,2) * t312;
t292 = rSges(2,1) * t312 + rSges(2,2) * t314;
t291 = Icges(2,1) * t314 - t353;
t290 = Icges(2,1) * t312 + t306;
t289 = -Icges(2,2) * t312 + t306;
t288 = Icges(2,2) * t314 + t353;
t282 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t281 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t280 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t262 = t315 * rSges(3,3) + (rSges(3,1) * t319 + rSges(3,2) * t321) * t313;
t261 = Icges(3,5) * t315 + (Icges(3,1) * t319 + Icges(3,4) * t321) * t313;
t260 = Icges(3,6) * t315 + (Icges(3,4) * t319 + Icges(3,2) * t321) * t313;
t259 = Icges(3,3) * t315 + (Icges(3,5) * t319 + Icges(3,6) * t321) * t313;
t258 = V_base(5) * rSges(2,3) - t292 * V_base(6) + t341;
t257 = t293 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t247 = t292 * V_base(4) - t293 * V_base(5) + t337;
t238 = t275 * t301 - t300 * t347;
t237 = -t275 * t300 - t301 * t347;
t234 = t275 * rSges(4,1) - t274 * rSges(4,2) - rSges(4,3) * t347;
t233 = Icges(4,1) * t275 - Icges(4,4) * t274 - Icges(4,5) * t347;
t232 = Icges(4,4) * t275 - Icges(4,2) * t274 - Icges(4,6) * t347;
t231 = Icges(4,5) * t275 - Icges(4,6) * t274 - Icges(4,3) * t347;
t230 = rSges(3,1) * t270 - rSges(3,2) * t269 + rSges(3,3) * t350;
t229 = rSges(3,1) * t268 - rSges(3,2) * t267 - rSges(3,3) * t349;
t228 = Icges(3,1) * t270 - Icges(3,4) * t269 + Icges(3,5) * t350;
t227 = Icges(3,1) * t268 - Icges(3,4) * t267 - Icges(3,5) * t349;
t226 = Icges(3,4) * t270 - Icges(3,2) * t269 + Icges(3,6) * t350;
t225 = Icges(3,4) * t268 - Icges(3,2) * t267 - Icges(3,6) * t349;
t224 = Icges(3,5) * t270 - Icges(3,6) * t269 + Icges(3,3) * t350;
t223 = Icges(3,5) * t268 - Icges(3,6) * t267 - Icges(3,3) * t349;
t220 = qJD(6) * t274 + t240;
t211 = t253 * t301 + t269 * t300;
t210 = -t253 * t300 + t269 * t301;
t209 = t251 * t301 + t267 * t300;
t208 = -t251 * t300 + t267 * t301;
t202 = rSges(5,1) * t255 + rSges(5,2) * t254 + rSges(5,3) * t274;
t200 = Icges(5,1) * t255 + Icges(5,4) * t254 + Icges(5,5) * t274;
t199 = Icges(5,4) * t255 + Icges(5,2) * t254 + Icges(5,6) * t274;
t197 = rSges(4,1) * t253 - rSges(4,2) * t252 + rSges(4,3) * t269;
t196 = rSges(4,1) * t251 - rSges(4,2) * t250 + rSges(4,3) * t267;
t195 = Icges(4,1) * t253 - Icges(4,4) * t252 + Icges(4,5) * t269;
t194 = Icges(4,1) * t251 - Icges(4,4) * t250 + Icges(4,5) * t267;
t193 = Icges(4,4) * t253 - Icges(4,2) * t252 + Icges(4,6) * t269;
t192 = Icges(4,4) * t251 - Icges(4,2) * t250 + Icges(4,6) * t267;
t191 = Icges(4,5) * t253 - Icges(4,6) * t252 + Icges(4,3) * t269;
t190 = Icges(4,5) * t251 - Icges(4,6) * t250 + Icges(4,3) * t267;
t188 = rSges(6,1) * t242 + rSges(6,2) * t241 + rSges(6,3) * t274;
t187 = Icges(6,1) * t242 + Icges(6,4) * t241 + Icges(6,5) * t274;
t186 = Icges(6,4) * t242 + Icges(6,2) * t241 + Icges(6,6) * t274;
t184 = qJD(6) * t252 + t207;
t183 = qJD(6) * t250 + t206;
t181 = rSges(7,1) * t238 + rSges(7,2) * t237 + rSges(7,3) * t274;
t180 = Icges(7,1) * t238 + Icges(7,4) * t237 + Icges(7,5) * t274;
t179 = Icges(7,4) * t238 + Icges(7,2) * t237 + Icges(7,6) * t274;
t178 = Icges(7,5) * t238 + Icges(7,6) * t237 + Icges(7,3) * t274;
t177 = -t229 * t296 + t262 * t284 + t333;
t176 = t230 * t296 - t262 * t285 + t331;
t175 = pkin(10) * t274 + t275 * t344 - t334 * t347;
t174 = rSges(5,1) * t219 + rSges(5,2) * t218 + rSges(5,3) * t252;
t173 = rSges(5,1) * t217 + rSges(5,2) * t216 + rSges(5,3) * t250;
t172 = Icges(5,1) * t219 + Icges(5,4) * t218 + Icges(5,5) * t252;
t171 = Icges(5,1) * t217 + Icges(5,4) * t216 + Icges(5,5) * t250;
t170 = Icges(5,4) * t219 + Icges(5,2) * t218 + Icges(5,6) * t252;
t169 = Icges(5,4) * t217 + Icges(5,2) * t216 + Icges(5,6) * t250;
t166 = t229 * t285 - t230 * t284 + t332;
t164 = rSges(6,1) * t215 + rSges(6,2) * t214 + rSges(6,3) * t252;
t163 = rSges(6,1) * t213 + rSges(6,2) * t212 + rSges(6,3) * t250;
t162 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t252;
t161 = Icges(6,1) * t213 + Icges(6,4) * t212 + Icges(6,5) * t250;
t160 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t252;
t159 = Icges(6,4) * t213 + Icges(6,2) * t212 + Icges(6,6) * t250;
t154 = rSges(7,1) * t211 + rSges(7,2) * t210 + rSges(7,3) * t252;
t153 = rSges(7,1) * t209 + rSges(7,2) * t208 + rSges(7,3) * t250;
t152 = Icges(7,1) * t211 + Icges(7,4) * t210 + Icges(7,5) * t252;
t151 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t250;
t150 = Icges(7,4) * t211 + Icges(7,2) * t210 + Icges(7,6) * t252;
t149 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t250;
t148 = Icges(7,5) * t211 + Icges(7,6) * t210 + Icges(7,3) * t252;
t147 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t250;
t145 = pkin(10) * t252 + t253 * t344 + t269 * t334;
t144 = pkin(10) * t250 + t251 * t344 + t267 * t334;
t142 = -t196 * t271 + t234 * t248 + t330;
t141 = t197 * t271 - t234 * t249 + t328;
t140 = t196 * t249 - t197 * t248 + t329;
t139 = -t173 * t240 + t202 * t206 + t327;
t138 = t174 * t240 - t202 * t207 + t325;
t137 = t173 * t207 - t174 * t206 + t326;
t136 = t188 * t206 + (-t155 - t163) * t240 + t324;
t135 = t164 * t240 + (-t188 - t189) * t207 + t322;
t134 = t163 * t207 + (-t156 - t164) * t206 + t323;
t133 = t324 - t153 * t220 + t175 * t206 + t181 * t183 + (-t144 - t155) * t240;
t132 = t145 * t240 + t154 * t220 - t181 * t184 + (-t175 - t189) * t207 + t322;
t131 = t323 + (-t145 - t156) * t206 + t144 * t207 + t153 * t184 - t154 * t183;
t1 = t296 * ((t223 * t284 + t224 * t285 + t259 * t296) * t315 + ((t226 * t321 + t228 * t319) * t285 + (t225 * t321 + t227 * t319) * t284 + (t260 * t321 + t261 * t319) * t296) * t313) / 0.2e1 + t271 * ((-t191 * t347 - t274 * t193 + t275 * t195) * t249 + (-t190 * t347 - t274 * t192 + t275 * t194) * t248 + (-t231 * t347 - t274 * t232 + t275 * t233) * t271) / 0.2e1 + t284 * ((-t224 * t349 - t226 * t267 + t228 * t268) * t285 + (-t223 * t349 - t267 * t225 + t268 * t227) * t284 + (-t259 * t349 - t260 * t267 + t261 * t268) * t296) / 0.2e1 + t285 * ((t224 * t350 - t269 * t226 + t270 * t228) * t285 + (t223 * t350 - t225 * t269 + t227 * t270) * t284 + (t259 * t350 - t260 * t269 + t261 * t270) * t296) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(3) * (t166 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + t183 * ((t148 * t250 + t150 * t208 + t152 * t209) * t184 + (t250 * t147 + t208 * t149 + t209 * t151) * t183 + (t178 * t250 + t179 * t208 + t180 * t209) * t220) / 0.2e1 + t184 * ((t148 * t252 + t210 * t150 + t211 * t152) * t184 + (t147 * t252 + t149 * t210 + t151 * t211) * t183 + (t178 * t252 + t179 * t210 + t180 * t211) * t220) / 0.2e1 + m(2) * (t247 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + t248 * ((t191 * t267 - t193 * t250 + t195 * t251) * t249 + (t267 * t190 - t250 * t192 + t251 * t194) * t248 + (t231 * t267 - t232 * t250 + t233 * t251) * t271) / 0.2e1 + t249 * ((t269 * t191 - t252 * t193 + t253 * t195) * t249 + (t190 * t269 - t192 * t252 + t194 * t253) * t248 + (t231 * t269 - t232 * t252 + t233 * t253) * t271) / 0.2e1 + t220 * ((t148 * t274 + t150 * t237 + t152 * t238) * t184 + (t147 * t274 + t149 * t237 + t151 * t238) * t183 + (t274 * t178 + t237 * t179 + t238 * t180) * t220) / 0.2e1 + m(1) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + ((t186 * t212 + t187 * t213 + t199 * t216 + t200 * t217 + t250 * t363) * t240 + (t160 * t212 + t162 * t213 + t170 * t216 + t172 * t217 + t250 * t364) * t207 + (t212 * t159 + t213 * t161 + t216 * t169 + t217 * t171 + t365 * t250) * t206) * t206 / 0.2e1 + ((t186 * t214 + t187 * t215 + t199 * t218 + t200 * t219 + t252 * t363) * t240 + (t214 * t160 + t162 * t215 + t218 * t170 + t219 * t172 + t364 * t252) * t207 + (t159 * t214 + t161 * t215 + t169 * t218 + t171 * t219 + t252 * t365) * t206) * t207 / 0.2e1 + ((t241 * t186 + t242 * t187 + t254 * t199 + t255 * t200 + t363 * t274) * t240 + (t160 * t241 + t162 * t242 + t170 * t254 + t172 * t255 + t274 * t364) * t207 + (t159 * t241 + t161 * t242 + t169 * t254 + t171 * t255 + t274 * t365) * t206) * t240 / 0.2e1 + ((-t288 * t312 + t290 * t314 + Icges(1,4)) * V_base(5) + (-t312 * t289 + t314 * t291 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t314 * t288 + t312 * t290 + Icges(1,2)) * V_base(5) + (t289 * t314 + t291 * t312 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t314 - Icges(2,6) * t312 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t312 + Icges(2,6) * t314 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
