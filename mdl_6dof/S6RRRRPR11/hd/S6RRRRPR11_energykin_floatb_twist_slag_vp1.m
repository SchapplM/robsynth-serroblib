% Calculate kinetic energy for
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:12:15
% EndTime: 2019-03-09 23:12:19
% DurationCPUTime: 4.69s
% Computational Cost: add. (3543->437), mult. (6882->630), div. (0->0), fcn. (8332->14), ass. (0->188)
t365 = Icges(5,3) + Icges(6,3);
t314 = cos(pkin(6));
t319 = sin(qJ(1));
t321 = cos(qJ(2));
t345 = t319 * t321;
t318 = sin(qJ(2));
t322 = cos(qJ(1));
t346 = t318 * t322;
t272 = t314 * t346 + t345;
t317 = sin(qJ(3));
t313 = sin(pkin(6));
t348 = t313 * t322;
t358 = cos(qJ(3));
t252 = t272 * t358 - t317 * t348;
t344 = t321 * t322;
t347 = t318 * t319;
t271 = -t314 * t344 + t347;
t312 = qJ(4) + pkin(12);
t304 = sin(t312);
t305 = cos(t312);
t212 = -t252 * t304 + t271 * t305;
t213 = t252 * t305 + t271 * t304;
t316 = sin(qJ(4));
t320 = cos(qJ(4));
t217 = -t252 * t316 + t271 * t320;
t352 = t271 * t316;
t218 = t252 * t320 + t352;
t336 = t313 * t358;
t251 = t272 * t317 + t322 * t336;
t364 = Icges(5,5) * t218 + Icges(6,5) * t213 + Icges(5,6) * t217 + Icges(6,6) * t212 + t365 * t251;
t274 = -t314 * t347 + t344;
t350 = t313 * t319;
t254 = t274 * t358 + t317 * t350;
t273 = t314 * t345 + t346;
t214 = -t254 * t304 + t273 * t305;
t215 = t254 * t305 + t273 * t304;
t219 = -t254 * t316 + t273 * t320;
t351 = t273 * t316;
t220 = t254 * t320 + t351;
t253 = t274 * t317 - t319 * t336;
t363 = Icges(5,5) * t220 + Icges(6,5) * t215 + Icges(5,6) * t219 + Icges(6,6) * t214 + t365 * t253;
t270 = t314 * t317 + t318 * t336;
t349 = t313 * t321;
t241 = -t270 * t304 - t305 * t349;
t242 = t270 * t305 - t304 * t349;
t247 = -t270 * t316 - t320 * t349;
t337 = t316 * t349;
t248 = t270 * t320 - t337;
t269 = t313 * t317 * t318 - t314 * t358;
t362 = Icges(5,5) * t248 + Icges(6,5) * t242 + Icges(5,6) * t247 + Icges(6,6) * t241 + t365 * t269;
t356 = pkin(8) * t314;
t355 = t320 * pkin(4);
t353 = Icges(2,4) * t319;
t343 = pkin(5) * t305;
t341 = qJD(2) * t313;
t340 = V_base(5) * pkin(7) + V_base(1);
t285 = t319 * t341 + V_base(4);
t307 = V_base(6) + qJD(1);
t335 = pkin(5) * t304;
t250 = qJD(3) * t273 + t285;
t286 = qJD(2) * t314 + t307;
t207 = qJD(4) * t253 + t250;
t284 = -t322 * t341 + V_base(5);
t278 = t319 * pkin(1) - pkin(8) * t348;
t334 = -t278 * t307 + V_base(5) * t356 + t340;
t279 = pkin(1) * t322 + pkin(8) * t350;
t333 = V_base(4) * t278 - t279 * V_base(5) + V_base(3);
t249 = qJD(3) * t271 + t284;
t206 = qJD(4) * t251 + t249;
t267 = -qJD(3) * t349 + t286;
t332 = t307 * t279 + V_base(2) + (-pkin(7) - t356) * V_base(4);
t238 = qJD(4) * t269 + t267;
t239 = pkin(2) * t272 + pkin(9) * t271;
t276 = (pkin(2) * t318 - pkin(9) * t321) * t313;
t331 = -t239 * t286 + t284 * t276 + t334;
t240 = pkin(2) * t274 + pkin(9) * t273;
t330 = t285 * t239 - t240 * t284 + t333;
t329 = t286 * t240 - t276 * t285 + t332;
t204 = pkin(3) * t252 + pkin(10) * t251;
t237 = pkin(3) * t270 + pkin(10) * t269;
t328 = -t204 * t267 + t249 * t237 + t331;
t205 = pkin(3) * t254 + pkin(10) * t253;
t327 = t250 * t204 - t205 * t249 + t330;
t189 = -pkin(4) * t337 + qJ(5) * t269 + t270 * t355;
t326 = qJD(5) * t253 + t206 * t189 + t328;
t155 = pkin(4) * t352 + qJ(5) * t251 + t252 * t355;
t325 = qJD(5) * t269 + t207 * t155 + t327;
t324 = t267 * t205 - t237 * t250 + t329;
t156 = pkin(4) * t351 + qJ(5) * t253 + t254 * t355;
t323 = qJD(5) * t251 + t238 * t156 + t324;
t308 = Icges(2,4) * t322;
t306 = qJ(6) + t312;
t301 = cos(t306);
t300 = sin(t306);
t294 = rSges(2,1) * t322 - t319 * rSges(2,2);
t293 = t319 * rSges(2,1) + rSges(2,2) * t322;
t292 = Icges(2,1) * t322 - t353;
t291 = Icges(2,1) * t319 + t308;
t290 = -Icges(2,2) * t319 + t308;
t289 = Icges(2,2) * t322 + t353;
t282 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t281 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t280 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t262 = rSges(3,3) * t314 + (rSges(3,1) * t318 + rSges(3,2) * t321) * t313;
t261 = Icges(3,5) * t314 + (Icges(3,1) * t318 + Icges(3,4) * t321) * t313;
t260 = Icges(3,6) * t314 + (Icges(3,4) * t318 + Icges(3,2) * t321) * t313;
t259 = Icges(3,3) * t314 + (Icges(3,5) * t318 + Icges(3,6) * t321) * t313;
t258 = V_base(5) * rSges(2,3) - t293 * t307 + t340;
t257 = t294 * t307 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t255 = t293 * V_base(4) - t294 * V_base(5) + V_base(3);
t236 = t270 * t301 - t300 * t349;
t235 = -t270 * t300 - t301 * t349;
t234 = rSges(3,1) * t274 - rSges(3,2) * t273 + rSges(3,3) * t350;
t233 = t272 * rSges(3,1) - t271 * rSges(3,2) - rSges(3,3) * t348;
t232 = Icges(3,1) * t274 - Icges(3,4) * t273 + Icges(3,5) * t350;
t231 = Icges(3,1) * t272 - Icges(3,4) * t271 - Icges(3,5) * t348;
t230 = Icges(3,4) * t274 - Icges(3,2) * t273 + Icges(3,6) * t350;
t229 = Icges(3,4) * t272 - Icges(3,2) * t271 - Icges(3,6) * t348;
t228 = Icges(3,5) * t274 - Icges(3,6) * t273 + Icges(3,3) * t350;
t227 = Icges(3,5) * t272 - Icges(3,6) * t271 - Icges(3,3) * t348;
t226 = rSges(4,1) * t270 - rSges(4,2) * t269 - rSges(4,3) * t349;
t225 = Icges(4,1) * t270 - Icges(4,4) * t269 - Icges(4,5) * t349;
t224 = Icges(4,4) * t270 - Icges(4,2) * t269 - Icges(4,6) * t349;
t223 = Icges(4,5) * t270 - Icges(4,6) * t269 - Icges(4,3) * t349;
t216 = qJD(6) * t269 + t238;
t211 = t254 * t301 + t273 * t300;
t210 = -t254 * t300 + t273 * t301;
t209 = t252 * t301 + t271 * t300;
t208 = -t252 * t300 + t271 * t301;
t202 = rSges(4,1) * t254 - rSges(4,2) * t253 + rSges(4,3) * t273;
t201 = rSges(4,1) * t252 - rSges(4,2) * t251 + rSges(4,3) * t271;
t199 = Icges(4,1) * t254 - Icges(4,4) * t253 + Icges(4,5) * t273;
t198 = Icges(4,1) * t252 - Icges(4,4) * t251 + Icges(4,5) * t271;
t197 = Icges(4,4) * t254 - Icges(4,2) * t253 + Icges(4,6) * t273;
t196 = Icges(4,4) * t252 - Icges(4,2) * t251 + Icges(4,6) * t271;
t195 = Icges(4,5) * t254 - Icges(4,6) * t253 + Icges(4,3) * t273;
t194 = Icges(4,5) * t252 - Icges(4,6) * t251 + Icges(4,3) * t271;
t193 = rSges(5,1) * t248 + rSges(5,2) * t247 + rSges(5,3) * t269;
t192 = Icges(5,1) * t248 + Icges(5,4) * t247 + Icges(5,5) * t269;
t191 = Icges(5,4) * t248 + Icges(5,2) * t247 + Icges(5,6) * t269;
t188 = rSges(6,1) * t242 + rSges(6,2) * t241 + rSges(6,3) * t269;
t187 = Icges(6,1) * t242 + Icges(6,4) * t241 + Icges(6,5) * t269;
t186 = Icges(6,4) * t242 + Icges(6,2) * t241 + Icges(6,6) * t269;
t184 = qJD(6) * t253 + t207;
t183 = qJD(6) * t251 + t206;
t181 = rSges(7,1) * t236 + rSges(7,2) * t235 + rSges(7,3) * t269;
t180 = Icges(7,1) * t236 + Icges(7,4) * t235 + Icges(7,5) * t269;
t179 = Icges(7,4) * t236 + Icges(7,2) * t235 + Icges(7,6) * t269;
t178 = Icges(7,5) * t236 + Icges(7,6) * t235 + Icges(7,3) * t269;
t177 = -t233 * t286 + t262 * t284 + t334;
t176 = t234 * t286 - t262 * t285 + t332;
t175 = pkin(11) * t269 + t270 * t343 - t335 * t349;
t174 = rSges(5,1) * t220 + rSges(5,2) * t219 + rSges(5,3) * t253;
t173 = rSges(5,1) * t218 + rSges(5,2) * t217 + rSges(5,3) * t251;
t172 = Icges(5,1) * t220 + Icges(5,4) * t219 + Icges(5,5) * t253;
t171 = Icges(5,1) * t218 + Icges(5,4) * t217 + Icges(5,5) * t251;
t170 = Icges(5,4) * t220 + Icges(5,2) * t219 + Icges(5,6) * t253;
t169 = Icges(5,4) * t218 + Icges(5,2) * t217 + Icges(5,6) * t251;
t166 = t233 * t285 - t234 * t284 + t333;
t165 = rSges(6,1) * t215 + rSges(6,2) * t214 + rSges(6,3) * t253;
t164 = rSges(6,1) * t213 + rSges(6,2) * t212 + rSges(6,3) * t251;
t163 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t253;
t162 = Icges(6,1) * t213 + Icges(6,4) * t212 + Icges(6,5) * t251;
t161 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t253;
t160 = Icges(6,4) * t213 + Icges(6,2) * t212 + Icges(6,6) * t251;
t154 = rSges(7,1) * t211 + rSges(7,2) * t210 + rSges(7,3) * t253;
t153 = rSges(7,1) * t209 + rSges(7,2) * t208 + rSges(7,3) * t251;
t152 = Icges(7,1) * t211 + Icges(7,4) * t210 + Icges(7,5) * t253;
t151 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t251;
t150 = Icges(7,4) * t211 + Icges(7,2) * t210 + Icges(7,6) * t253;
t149 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t251;
t148 = Icges(7,5) * t211 + Icges(7,6) * t210 + Icges(7,3) * t253;
t147 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t251;
t145 = pkin(11) * t253 + t254 * t343 + t273 * t335;
t144 = pkin(11) * t251 + t252 * t343 + t271 * t335;
t142 = -t201 * t267 + t226 * t249 + t331;
t141 = t202 * t267 - t226 * t250 + t329;
t140 = t201 * t250 - t202 * t249 + t330;
t139 = -t173 * t238 + t193 * t206 + t328;
t138 = t174 * t238 - t193 * t207 + t324;
t137 = t173 * t207 - t174 * t206 + t327;
t136 = t188 * t206 + (-t155 - t164) * t238 + t326;
t135 = t165 * t238 + (-t188 - t189) * t207 + t323;
t134 = t164 * t207 + (-t156 - t165) * t206 + t325;
t133 = (-t144 - t155) * t238 - t153 * t216 + t175 * t206 + t181 * t183 + t326;
t132 = t145 * t238 + t154 * t216 - t181 * t184 + (-t175 - t189) * t207 + t323;
t131 = t144 * t207 + t153 * t184 - t154 * t183 + (-t145 - t156) * t206 + t325;
t1 = m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + t267 * ((-t195 * t349 - t197 * t269 + t199 * t270) * t250 + (-t194 * t349 - t196 * t269 + t198 * t270) * t249 + (-t223 * t349 - t224 * t269 + t225 * t270) * t267) / 0.2e1 + t285 * ((t228 * t350 - t230 * t273 + t232 * t274) * t285 + (t227 * t350 - t229 * t273 + t231 * t274) * t284 + (t259 * t350 - t260 * t273 + t261 * t274) * t286) / 0.2e1 + t284 * ((-t228 * t348 - t271 * t230 + t272 * t232) * t285 + (-t227 * t348 - t271 * t229 + t272 * t231) * t284 + (-t259 * t348 - t271 * t260 + t272 * t261) * t286) / 0.2e1 + t286 * ((t227 * t284 + t228 * t285 + t259 * t286) * t314 + ((t230 * t321 + t232 * t318) * t285 + (t229 * t321 + t231 * t318) * t284 + (t260 * t321 + t261 * t318) * t286) * t313) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(3) * (t166 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + t183 * ((t148 * t251 + t150 * t208 + t152 * t209) * t184 + (t147 * t251 + t149 * t208 + t151 * t209) * t183 + (t178 * t251 + t179 * t208 + t180 * t209) * t216) / 0.2e1 + t184 * ((t148 * t253 + t150 * t210 + t152 * t211) * t184 + (t147 * t253 + t149 * t210 + t151 * t211) * t183 + (t178 * t253 + t179 * t210 + t180 * t211) * t216) / 0.2e1 + m(2) * (t255 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + t216 * ((t148 * t269 + t150 * t235 + t152 * t236) * t184 + (t147 * t269 + t149 * t235 + t151 * t236) * t183 + (t178 * t269 + t179 * t235 + t180 * t236) * t216) / 0.2e1 + t249 * ((t195 * t271 - t197 * t251 + t199 * t252) * t250 + (t194 * t271 - t196 * t251 + t198 * t252) * t249 + (t223 * t271 - t224 * t251 + t225 * t252) * t267) / 0.2e1 + t250 * ((t195 * t273 - t197 * t253 + t199 * t254) * t250 + (t194 * t273 - t196 * t253 + t198 * t254) * t249 + (t223 * t273 - t224 * t253 + t225 * t254) * t267) / 0.2e1 + m(1) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + ((t186 * t212 + t187 * t213 + t191 * t217 + t192 * t218 + t251 * t362) * t238 + (t161 * t212 + t163 * t213 + t170 * t217 + t172 * t218 + t251 * t363) * t207 + (t160 * t212 + t162 * t213 + t169 * t217 + t171 * t218 + t364 * t251) * t206) * t206 / 0.2e1 + ((t186 * t214 + t187 * t215 + t191 * t219 + t192 * t220 + t253 * t362) * t238 + (t161 * t214 + t163 * t215 + t170 * t219 + t172 * t220 + t363 * t253) * t207 + (t160 * t214 + t162 * t215 + t169 * t219 + t171 * t220 + t253 * t364) * t206) * t207 / 0.2e1 + ((t186 * t241 + t187 * t242 + t191 * t247 + t192 * t248 + t362 * t269) * t238 + (t161 * t241 + t163 * t242 + t170 * t247 + t172 * t248 + t269 * t363) * t207 + (t160 * t241 + t162 * t242 + t169 * t247 + t171 * t248 + t269 * t364) * t206) * t238 / 0.2e1 + ((-t319 * t289 + t291 * t322 + Icges(1,4)) * V_base(5) + (-t319 * t290 + t292 * t322 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t289 * t322 + t319 * t291 + Icges(1,2)) * V_base(5) + (t290 * t322 + t319 * t292 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t319 + Icges(2,6) * t322) * V_base(5) + (Icges(2,5) * t322 - Icges(2,6) * t319) * V_base(4) + Icges(2,3) * t307 / 0.2e1) * t307;
T  = t1;
