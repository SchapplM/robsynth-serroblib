% Calculate kinetic energy for
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:12
% EndTime: 2019-03-08 22:16:17
% DurationCPUTime: 4.77s
% Computational Cost: add. (3411->437), mult. (6720->628), div. (0->0), fcn. (8125->14), ass. (0->187)
t367 = Icges(4,2) + Icges(5,3);
t313 = sin(pkin(11));
t316 = cos(pkin(11));
t321 = cos(qJ(2));
t317 = cos(pkin(6));
t320 = sin(qJ(2));
t347 = t317 * t320;
t268 = t313 * t321 + t316 * t347;
t314 = sin(pkin(6));
t319 = sin(qJ(3));
t349 = t314 * t319;
t358 = cos(qJ(3));
t253 = t268 * t358 - t316 * t349;
t346 = t317 * t321;
t267 = t313 * t320 - t316 * t346;
t312 = sin(pkin(12));
t315 = cos(pkin(12));
t216 = -t253 * t312 + t267 * t315;
t353 = t267 * t312;
t217 = t253 * t315 + t353;
t335 = t314 * t358;
t252 = t268 * t319 + t316 * t335;
t364 = -Icges(4,4) * t253 + Icges(5,5) * t217 - Icges(4,6) * t267 + Icges(5,6) * t216 + t367 * t252;
t270 = -t313 * t347 + t316 * t321;
t255 = t270 * t358 + t313 * t349;
t269 = t313 * t346 + t316 * t320;
t218 = -t255 * t312 + t269 * t315;
t352 = t269 * t312;
t219 = t255 * t315 + t352;
t254 = t270 * t319 - t313 * t335;
t363 = -Icges(4,4) * t255 + Icges(5,5) * t219 - Icges(4,6) * t269 + Icges(5,6) * t218 + t367 * t254;
t275 = t317 * t319 + t320 * t335;
t348 = t314 * t321;
t250 = -t275 * t312 - t315 * t348;
t336 = t312 * t348;
t251 = t275 * t315 - t336;
t274 = -t317 * t358 + t320 * t349;
t362 = -Icges(4,4) * t275 + Icges(5,5) * t251 + Icges(4,6) * t348 + Icges(5,6) * t250 + t367 * t274;
t356 = pkin(7) * t317;
t355 = t315 * pkin(4);
t354 = Icges(2,4) * t313;
t351 = t313 * t314;
t350 = t314 * t316;
t311 = pkin(12) + qJ(5);
t305 = cos(t311);
t344 = pkin(5) * t305;
t342 = qJD(2) * t314;
t341 = V_base(5) * qJ(1) + V_base(1);
t337 = qJD(1) + V_base(3);
t285 = t313 * t342 + V_base(4);
t296 = qJD(2) * t317 + V_base(6);
t304 = sin(t311);
t334 = pkin(5) * t304;
t249 = qJD(3) * t269 + t285;
t207 = qJD(5) * t254 + t249;
t284 = -t316 * t342 + V_base(5);
t248 = qJD(3) * t267 + t284;
t271 = -qJD(3) * t348 + t296;
t277 = pkin(1) * t313 - pkin(7) * t350;
t333 = -t277 * V_base(6) + V_base(5) * t356 + t341;
t278 = pkin(1) * t316 + pkin(7) * t351;
t332 = V_base(4) * t277 - t278 * V_base(5) + t337;
t206 = qJD(5) * t252 + t248;
t240 = qJD(5) * t274 + t271;
t331 = V_base(6) * t278 + V_base(2) + (-qJ(1) - t356) * V_base(4);
t235 = pkin(2) * t268 + pkin(8) * t267;
t276 = (pkin(2) * t320 - pkin(8) * t321) * t314;
t330 = -t235 * t296 + t284 * t276 + t333;
t236 = pkin(2) * t270 + pkin(8) * t269;
t329 = t285 * t235 - t236 * t284 + t332;
t328 = t296 * t236 - t276 * t285 + t331;
t237 = t275 * pkin(3) + t274 * qJ(4);
t327 = qJD(4) * t254 + t248 * t237 + t330;
t204 = pkin(3) * t253 + qJ(4) * t252;
t326 = qJD(4) * t274 + t249 * t204 + t329;
t205 = pkin(3) * t255 + qJ(4) * t254;
t325 = qJD(4) * t252 + t271 * t205 + t328;
t155 = pkin(4) * t353 + pkin(9) * t252 + t253 * t355;
t189 = -pkin(4) * t336 + pkin(9) * t274 + t275 * t355;
t324 = t248 * t189 + (-t155 - t204) * t271 + t327;
t156 = pkin(4) * t352 + pkin(9) * t254 + t255 * t355;
t323 = t249 * t155 + (-t156 - t205) * t248 + t326;
t322 = t271 * t156 + (-t189 - t237) * t249 + t325;
t307 = qJ(6) + t311;
t306 = Icges(2,4) * t316;
t301 = cos(t307);
t300 = sin(t307);
t293 = rSges(2,1) * t316 - rSges(2,2) * t313;
t292 = rSges(2,1) * t313 + rSges(2,2) * t316;
t291 = Icges(2,1) * t316 - t354;
t290 = Icges(2,1) * t313 + t306;
t289 = -Icges(2,2) * t313 + t306;
t288 = Icges(2,2) * t316 + t354;
t283 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t282 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t281 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t262 = t317 * rSges(3,3) + (rSges(3,1) * t320 + rSges(3,2) * t321) * t314;
t261 = Icges(3,5) * t317 + (Icges(3,1) * t320 + Icges(3,4) * t321) * t314;
t260 = Icges(3,6) * t317 + (Icges(3,4) * t320 + Icges(3,2) * t321) * t314;
t259 = Icges(3,3) * t317 + (Icges(3,5) * t320 + Icges(3,6) * t321) * t314;
t258 = V_base(5) * rSges(2,3) - t292 * V_base(6) + t341;
t257 = t293 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t247 = t292 * V_base(4) - t293 * V_base(5) + t337;
t242 = t275 * t305 - t304 * t348;
t241 = -t275 * t304 - t305 * t348;
t239 = t275 * t301 - t300 * t348;
t238 = -t275 * t300 - t301 * t348;
t234 = t275 * rSges(4,1) - t274 * rSges(4,2) - rSges(4,3) * t348;
t233 = Icges(4,1) * t275 - Icges(4,4) * t274 - Icges(4,5) * t348;
t231 = Icges(4,5) * t275 - Icges(4,6) * t274 - Icges(4,3) * t348;
t230 = rSges(3,1) * t270 - rSges(3,2) * t269 + rSges(3,3) * t351;
t229 = rSges(3,1) * t268 - rSges(3,2) * t267 - rSges(3,3) * t350;
t228 = Icges(3,1) * t270 - Icges(3,4) * t269 + Icges(3,5) * t351;
t227 = Icges(3,1) * t268 - Icges(3,4) * t267 - Icges(3,5) * t350;
t226 = Icges(3,4) * t270 - Icges(3,2) * t269 + Icges(3,6) * t351;
t225 = Icges(3,4) * t268 - Icges(3,2) * t267 - Icges(3,6) * t350;
t224 = Icges(3,5) * t270 - Icges(3,6) * t269 + Icges(3,3) * t351;
t223 = Icges(3,5) * t268 - Icges(3,6) * t267 - Icges(3,3) * t350;
t220 = qJD(6) * t274 + t240;
t215 = t255 * t305 + t269 * t304;
t214 = -t255 * t304 + t269 * t305;
t213 = t253 * t305 + t267 * t304;
t212 = -t253 * t304 + t267 * t305;
t211 = t255 * t301 + t269 * t300;
t210 = -t255 * t300 + t269 * t301;
t209 = t253 * t301 + t267 * t300;
t208 = -t253 * t300 + t267 * t301;
t201 = rSges(5,1) * t251 + rSges(5,2) * t250 + rSges(5,3) * t274;
t200 = rSges(4,1) * t255 - rSges(4,2) * t254 + rSges(4,3) * t269;
t199 = rSges(4,1) * t253 - rSges(4,2) * t252 + rSges(4,3) * t267;
t198 = Icges(5,1) * t251 + Icges(5,4) * t250 + Icges(5,5) * t274;
t197 = Icges(5,4) * t251 + Icges(5,2) * t250 + Icges(5,6) * t274;
t195 = Icges(4,1) * t255 - Icges(4,4) * t254 + Icges(4,5) * t269;
t194 = Icges(4,1) * t253 - Icges(4,4) * t252 + Icges(4,5) * t267;
t191 = Icges(4,5) * t255 - Icges(4,6) * t254 + Icges(4,3) * t269;
t190 = Icges(4,5) * t253 - Icges(4,6) * t252 + Icges(4,3) * t267;
t188 = rSges(6,1) * t242 + rSges(6,2) * t241 + rSges(6,3) * t274;
t187 = Icges(6,1) * t242 + Icges(6,4) * t241 + Icges(6,5) * t274;
t186 = Icges(6,4) * t242 + Icges(6,2) * t241 + Icges(6,6) * t274;
t185 = Icges(6,5) * t242 + Icges(6,6) * t241 + Icges(6,3) * t274;
t184 = qJD(6) * t254 + t207;
t183 = qJD(6) * t252 + t206;
t182 = rSges(7,1) * t239 + rSges(7,2) * t238 + rSges(7,3) * t274;
t180 = Icges(7,1) * t239 + Icges(7,4) * t238 + Icges(7,5) * t274;
t179 = Icges(7,4) * t239 + Icges(7,2) * t238 + Icges(7,6) * t274;
t178 = Icges(7,5) * t239 + Icges(7,6) * t238 + Icges(7,3) * t274;
t177 = -t229 * t296 + t262 * t284 + t333;
t176 = t230 * t296 - t262 * t285 + t331;
t174 = pkin(10) * t274 + t275 * t344 - t334 * t348;
t173 = rSges(5,1) * t219 + rSges(5,2) * t218 + rSges(5,3) * t254;
t172 = rSges(5,1) * t217 + rSges(5,2) * t216 + rSges(5,3) * t252;
t171 = Icges(5,1) * t219 + Icges(5,4) * t218 + Icges(5,5) * t254;
t170 = Icges(5,1) * t217 + Icges(5,4) * t216 + Icges(5,5) * t252;
t169 = Icges(5,4) * t219 + Icges(5,2) * t218 + Icges(5,6) * t254;
t168 = Icges(5,4) * t217 + Icges(5,2) * t216 + Icges(5,6) * t252;
t165 = t229 * t285 - t230 * t284 + t332;
t164 = rSges(6,1) * t215 + rSges(6,2) * t214 + rSges(6,3) * t254;
t163 = rSges(6,1) * t213 + rSges(6,2) * t212 + rSges(6,3) * t252;
t162 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t254;
t161 = Icges(6,1) * t213 + Icges(6,4) * t212 + Icges(6,5) * t252;
t160 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t254;
t159 = Icges(6,4) * t213 + Icges(6,2) * t212 + Icges(6,6) * t252;
t158 = Icges(6,5) * t215 + Icges(6,6) * t214 + Icges(6,3) * t254;
t157 = Icges(6,5) * t213 + Icges(6,6) * t212 + Icges(6,3) * t252;
t154 = rSges(7,1) * t211 + rSges(7,2) * t210 + rSges(7,3) * t254;
t153 = rSges(7,1) * t209 + rSges(7,2) * t208 + rSges(7,3) * t252;
t152 = Icges(7,1) * t211 + Icges(7,4) * t210 + Icges(7,5) * t254;
t151 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t252;
t150 = Icges(7,4) * t211 + Icges(7,2) * t210 + Icges(7,6) * t254;
t149 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t252;
t148 = Icges(7,5) * t211 + Icges(7,6) * t210 + Icges(7,3) * t254;
t147 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t252;
t144 = pkin(10) * t254 + t255 * t344 + t269 * t334;
t143 = pkin(10) * t252 + t253 * t344 + t267 * t334;
t142 = -t199 * t271 + t234 * t248 + t330;
t141 = t200 * t271 - t234 * t249 + t328;
t140 = t199 * t249 - t200 * t248 + t329;
t139 = t201 * t248 + (-t172 - t204) * t271 + t327;
t138 = t173 * t271 + (-t201 - t237) * t249 + t325;
t137 = t172 * t249 + (-t173 - t205) * t248 + t326;
t136 = -t163 * t240 + t188 * t206 + t324;
t135 = t164 * t240 - t188 * t207 + t322;
t134 = t163 * t207 - t164 * t206 + t323;
t133 = -t143 * t240 - t153 * t220 + t174 * t206 + t182 * t183 + t324;
t132 = t144 * t240 + t154 * t220 - t174 * t207 - t182 * t184 + t322;
t131 = t143 * t207 - t144 * t206 + t153 * t184 - t154 * t183 + t323;
t1 = t296 * ((t223 * t284 + t224 * t285 + t259 * t296) * t317 + ((t226 * t321 + t228 * t320) * t285 + (t225 * t321 + t227 * t320) * t284 + (t260 * t321 + t261 * t320) * t296) * t314) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + t284 * ((-t224 * t350 - t226 * t267 + t228 * t268) * t285 + (-t223 * t350 - t267 * t225 + t268 * t227) * t284 + (-t259 * t350 - t260 * t267 + t261 * t268) * t296) / 0.2e1 + t285 * ((t224 * t351 - t269 * t226 + t270 * t228) * t285 + (t223 * t351 - t225 * t269 + t227 * t270) * t284 + (t259 * t351 - t260 * t269 + t261 * t270) * t296) / 0.2e1 + m(3) * (t165 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + t183 * ((t148 * t252 + t150 * t208 + t152 * t209) * t184 + (t252 * t147 + t149 * t208 + t209 * t151) * t183 + (t178 * t252 + t179 * t208 + t180 * t209) * t220) / 0.2e1 + t206 * ((t158 * t252 + t160 * t212 + t162 * t213) * t207 + (t252 * t157 + t212 * t159 + t161 * t213) * t206 + (t185 * t252 + t186 * t212 + t187 * t213) * t240) / 0.2e1 + t184 * ((t254 * t148 + t210 * t150 + t211 * t152) * t184 + (t147 * t254 + t149 * t210 + t151 * t211) * t183 + (t178 * t254 + t179 * t210 + t180 * t211) * t220) / 0.2e1 + t207 * ((t254 * t158 + t160 * t214 + t215 * t162) * t207 + (t157 * t254 + t159 * t214 + t161 * t215) * t206 + (t185 * t254 + t186 * t214 + t187 * t215) * t240) / 0.2e1 + m(2) * (t247 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + t220 * ((t148 * t274 + t150 * t238 + t152 * t239) * t184 + (t147 * t274 + t149 * t238 + t151 * t239) * t183 + (t274 * t178 + t238 * t179 + t239 * t180) * t220) / 0.2e1 + t240 * ((t158 * t274 + t160 * t241 + t162 * t242) * t207 + (t157 * t274 + t159 * t241 + t161 * t242) * t206 + (t185 * t274 + t241 * t186 + t242 * t187) * t240) / 0.2e1 + m(1) * (t281 ^ 2 + t282 ^ 2 + t283 ^ 2) / 0.2e1 + ((t197 * t216 + t198 * t217 + t231 * t267 + t233 * t253 + t252 * t362) * t271 + (t169 * t216 + t171 * t217 + t191 * t267 + t195 * t253 + t252 * t363) * t249 + (t216 * t168 + t217 * t170 + t267 * t190 + t253 * t194 + t364 * t252) * t248) * t248 / 0.2e1 + ((t197 * t218 + t198 * t219 + t231 * t269 + t233 * t255 + t254 * t362) * t271 + (t218 * t169 + t171 * t219 + t269 * t191 + t255 * t195 + t363 * t254) * t249 + (t168 * t218 + t170 * t219 + t190 * t269 + t194 * t255 + t254 * t364) * t248) * t249 / 0.2e1 + ((t250 * t197 + t251 * t198 - t231 * t348 + t275 * t233 + t362 * t274) * t271 + (t169 * t250 + t171 * t251 - t191 * t348 + t275 * t195 + t274 * t363) * t249 + (t168 * t250 + t170 * t251 - t190 * t348 + t275 * t194 + t274 * t364) * t248) * t271 / 0.2e1 + ((-t288 * t313 + t290 * t316 + Icges(1,4)) * V_base(5) + (-t313 * t289 + t316 * t291 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t316 * t288 + t313 * t290 + Icges(1,2)) * V_base(5) + (t289 * t316 + t291 * t313 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t313 + Icges(2,6) * t316 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t316 - Icges(2,6) * t313 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
