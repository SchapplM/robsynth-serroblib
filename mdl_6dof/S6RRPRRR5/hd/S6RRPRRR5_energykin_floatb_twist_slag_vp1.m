% Calculate kinetic energy for
% S6RRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:37:56
% EndTime: 2019-03-09 13:38:01
% DurationCPUTime: 4.31s
% Computational Cost: add. (4041->439), mult. (9164->640), div. (0->0), fcn. (11536->14), ass. (0->197)
t327 = sin(qJ(1));
t329 = cos(qJ(1));
t326 = sin(qJ(2));
t364 = sin(pkin(12));
t365 = cos(pkin(12));
t372 = cos(qJ(2));
t288 = -t326 * t364 + t372 * t365;
t323 = cos(pkin(6));
t338 = t323 * t288;
t339 = t326 * t365 + t364 * t372;
t259 = -t327 * t339 + t329 * t338;
t279 = t339 * t323;
t260 = t279 * t329 + t288 * t327;
t322 = sin(pkin(6));
t358 = t322 * t329;
t211 = Icges(4,5) * t260 + Icges(4,6) * t259 - Icges(4,3) * t358;
t348 = t329 * t372;
t357 = t327 * t326;
t282 = t323 * t348 - t357;
t349 = t327 * t372;
t356 = t329 * t326;
t283 = t323 * t356 + t349;
t246 = Icges(3,5) * t283 + Icges(3,6) * t282 - Icges(3,3) * t358;
t378 = t211 + t246;
t261 = -t327 * t338 - t329 * t339;
t262 = -t279 * t327 + t288 * t329;
t359 = t322 * t327;
t212 = Icges(4,5) * t262 + Icges(4,6) * t261 + Icges(4,3) * t359;
t284 = -t323 * t349 - t356;
t285 = -t323 * t357 + t348;
t247 = Icges(3,5) * t285 + Icges(3,6) * t284 + Icges(3,3) * t359;
t377 = t212 + t247;
t277 = t288 * t322;
t278 = t339 * t322;
t242 = Icges(4,5) * t278 + Icges(4,6) * t277 + Icges(4,3) * t323;
t273 = Icges(3,3) * t323 + (Icges(3,5) * t326 + Icges(3,6) * t372) * t322;
t376 = t242 + t273;
t371 = cos(qJ(4));
t370 = pkin(2) * t326;
t369 = pkin(8) * t323;
t368 = pkin(2) * t372;
t328 = cos(qJ(5));
t367 = pkin(5) * t328;
t363 = Icges(2,4) * t327;
t324 = sin(qJ(5));
t362 = t259 * t324;
t361 = t261 * t324;
t360 = t277 * t324;
t355 = qJD(2) * t322;
t354 = qJD(3) * t322;
t353 = V_base(5) * pkin(7) + V_base(1);
t350 = t322 * t371;
t296 = t327 * t355 + V_base(4);
t316 = V_base(6) + qJD(1);
t347 = -qJ(3) * t322 + t323 * t370;
t236 = -qJD(4) * t261 + t296;
t297 = qJD(2) * t323 + t316;
t325 = sin(qJ(4));
t239 = t262 * t325 - t327 * t350;
t199 = qJD(5) * t239 + t236;
t264 = -qJD(4) * t277 + t297;
t295 = -t329 * t355 + V_base(5);
t290 = pkin(1) * t327 - pkin(8) * t358;
t344 = -t290 * t316 + V_base(5) * t369 + t353;
t265 = t278 * t325 - t323 * t371;
t226 = qJD(5) * t265 + t264;
t291 = pkin(1) * t329 + pkin(8) * t359;
t343 = V_base(4) * t290 - t291 * V_base(5) + V_base(3);
t235 = -qJD(4) * t259 + t295;
t237 = t260 * t325 + t329 * t350;
t198 = qJD(5) * t237 + t235;
t289 = qJ(3) * t323 + t322 * t370;
t342 = t295 * t289 + t327 * t354 + t344;
t257 = t327 * t368 + t329 * t347;
t341 = qJD(3) * t323 + t296 * t257 + t343;
t340 = t316 * t291 + V_base(2) + (-pkin(7) - t369) * V_base(4);
t223 = pkin(3) * t260 - pkin(9) * t259;
t254 = pkin(3) * t278 - pkin(9) * t277;
t337 = t295 * t254 + (-t223 - t257) * t297 + t342;
t224 = pkin(3) * t262 - pkin(9) * t261;
t258 = -t327 * t347 + t329 * t368;
t336 = t296 * t223 + (-t224 - t258) * t295 + t341;
t335 = t297 * t258 - t329 * t354 + t340;
t238 = t260 * t371 - t325 * t358;
t196 = t238 * pkin(4) + t237 * pkin(10);
t266 = t278 * t371 + t323 * t325;
t225 = t266 * pkin(4) + t265 * pkin(10);
t334 = -t196 * t264 + t235 * t225 + t337;
t240 = t262 * t371 + t325 * t359;
t197 = t240 * pkin(4) + t239 * pkin(10);
t333 = t236 * t196 - t197 * t235 + t336;
t332 = t297 * t224 + (-t254 - t289) * t296 + t335;
t331 = t264 * t197 - t225 * t236 + t332;
t321 = qJ(5) + qJ(6);
t319 = Icges(2,4) * t329;
t318 = cos(t321);
t317 = sin(t321);
t305 = rSges(2,1) * t329 - rSges(2,2) * t327;
t304 = rSges(2,1) * t327 + rSges(2,2) * t329;
t303 = Icges(2,1) * t329 - t363;
t302 = Icges(2,1) * t327 + t319;
t301 = -Icges(2,2) * t327 + t319;
t300 = Icges(2,2) * t329 + t363;
t294 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t293 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t292 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t276 = t323 * rSges(3,3) + (rSges(3,1) * t326 + rSges(3,2) * t372) * t322;
t275 = Icges(3,5) * t323 + (Icges(3,1) * t326 + Icges(3,4) * t372) * t322;
t274 = Icges(3,6) * t323 + (Icges(3,4) * t326 + Icges(3,2) * t372) * t322;
t269 = V_base(5) * rSges(2,3) - t304 * t316 + t353;
t268 = t305 * t316 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t267 = t304 * V_base(4) - t305 * V_base(5) + V_base(3);
t253 = rSges(3,1) * t285 + rSges(3,2) * t284 + rSges(3,3) * t359;
t252 = rSges(3,1) * t283 + rSges(3,2) * t282 - rSges(3,3) * t358;
t251 = Icges(3,1) * t285 + Icges(3,4) * t284 + Icges(3,5) * t359;
t250 = Icges(3,1) * t283 + Icges(3,4) * t282 - Icges(3,5) * t358;
t249 = Icges(3,4) * t285 + Icges(3,2) * t284 + Icges(3,6) * t359;
t248 = Icges(3,4) * t283 + Icges(3,2) * t282 - Icges(3,6) * t358;
t245 = rSges(4,1) * t278 + rSges(4,2) * t277 + rSges(4,3) * t323;
t244 = Icges(4,1) * t278 + Icges(4,4) * t277 + Icges(4,5) * t323;
t243 = Icges(4,4) * t278 + Icges(4,2) * t277 + Icges(4,6) * t323;
t230 = t266 * t328 - t360;
t229 = -t266 * t324 - t277 * t328;
t228 = t266 * t318 - t277 * t317;
t227 = -t266 * t317 - t277 * t318;
t222 = rSges(5,1) * t266 - rSges(5,2) * t265 - rSges(5,3) * t277;
t221 = Icges(5,1) * t266 - Icges(5,4) * t265 - Icges(5,5) * t277;
t220 = Icges(5,4) * t266 - Icges(5,2) * t265 - Icges(5,6) * t277;
t219 = Icges(5,5) * t266 - Icges(5,6) * t265 - Icges(5,3) * t277;
t218 = rSges(4,1) * t262 + rSges(4,2) * t261 + rSges(4,3) * t359;
t217 = rSges(4,1) * t260 + rSges(4,2) * t259 - rSges(4,3) * t358;
t216 = Icges(4,1) * t262 + Icges(4,4) * t261 + Icges(4,5) * t359;
t215 = Icges(4,1) * t260 + Icges(4,4) * t259 - Icges(4,5) * t358;
t214 = Icges(4,4) * t262 + Icges(4,2) * t261 + Icges(4,6) * t359;
t213 = Icges(4,4) * t260 + Icges(4,2) * t259 - Icges(4,6) * t358;
t209 = qJD(6) * t265 + t226;
t207 = t240 * t328 - t361;
t206 = -t240 * t324 - t261 * t328;
t205 = t238 * t328 - t362;
t204 = -t238 * t324 - t259 * t328;
t203 = t240 * t318 - t261 * t317;
t202 = -t240 * t317 - t261 * t318;
t201 = t238 * t318 - t259 * t317;
t200 = -t238 * t317 - t259 * t318;
t195 = -t252 * t297 + t276 * t295 + t344;
t194 = t253 * t297 - t276 * t296 + t340;
t191 = t252 * t296 - t253 * t295 + t343;
t190 = rSges(6,1) * t230 + rSges(6,2) * t229 + rSges(6,3) * t265;
t189 = Icges(6,1) * t230 + Icges(6,4) * t229 + Icges(6,5) * t265;
t188 = Icges(6,4) * t230 + Icges(6,2) * t229 + Icges(6,6) * t265;
t187 = Icges(6,5) * t230 + Icges(6,6) * t229 + Icges(6,3) * t265;
t186 = rSges(5,1) * t240 - rSges(5,2) * t239 - rSges(5,3) * t261;
t185 = rSges(5,1) * t238 - rSges(5,2) * t237 - rSges(5,3) * t259;
t184 = Icges(5,1) * t240 - Icges(5,4) * t239 - Icges(5,5) * t261;
t183 = Icges(5,1) * t238 - Icges(5,4) * t237 - Icges(5,5) * t259;
t182 = Icges(5,4) * t240 - Icges(5,2) * t239 - Icges(5,6) * t261;
t181 = Icges(5,4) * t238 - Icges(5,2) * t237 - Icges(5,6) * t259;
t180 = Icges(5,5) * t240 - Icges(5,6) * t239 - Icges(5,3) * t261;
t179 = Icges(5,5) * t238 - Icges(5,6) * t237 - Icges(5,3) * t259;
t178 = -pkin(5) * t360 + pkin(11) * t265 + t266 * t367;
t177 = rSges(7,1) * t228 + rSges(7,2) * t227 + rSges(7,3) * t265;
t176 = Icges(7,1) * t228 + Icges(7,4) * t227 + Icges(7,5) * t265;
t175 = Icges(7,4) * t228 + Icges(7,2) * t227 + Icges(7,6) * t265;
t174 = Icges(7,5) * t228 + Icges(7,6) * t227 + Icges(7,3) * t265;
t173 = qJD(6) * t239 + t199;
t172 = qJD(6) * t237 + t198;
t170 = rSges(6,1) * t207 + rSges(6,2) * t206 + rSges(6,3) * t239;
t169 = rSges(6,1) * t205 + rSges(6,2) * t204 + rSges(6,3) * t237;
t168 = Icges(6,1) * t207 + Icges(6,4) * t206 + Icges(6,5) * t239;
t167 = Icges(6,1) * t205 + Icges(6,4) * t204 + Icges(6,5) * t237;
t166 = Icges(6,4) * t207 + Icges(6,2) * t206 + Icges(6,6) * t239;
t165 = Icges(6,4) * t205 + Icges(6,2) * t204 + Icges(6,6) * t237;
t164 = Icges(6,5) * t207 + Icges(6,6) * t206 + Icges(6,3) * t239;
t163 = Icges(6,5) * t205 + Icges(6,6) * t204 + Icges(6,3) * t237;
t162 = rSges(7,1) * t203 + rSges(7,2) * t202 + rSges(7,3) * t239;
t161 = rSges(7,1) * t201 + rSges(7,2) * t200 + rSges(7,3) * t237;
t160 = Icges(7,1) * t203 + Icges(7,4) * t202 + Icges(7,5) * t239;
t159 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t237;
t158 = Icges(7,4) * t203 + Icges(7,2) * t202 + Icges(7,6) * t239;
t157 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t237;
t156 = Icges(7,5) * t203 + Icges(7,6) * t202 + Icges(7,3) * t239;
t155 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t237;
t154 = -pkin(5) * t361 + pkin(11) * t239 + t240 * t367;
t153 = -pkin(5) * t362 + pkin(11) * t237 + t238 * t367;
t152 = t245 * t295 + (-t217 - t257) * t297 + t342;
t151 = t218 * t297 + (-t245 - t289) * t296 + t335;
t150 = t217 * t296 + (-t218 - t258) * t295 + t341;
t149 = -t185 * t264 + t222 * t235 + t337;
t148 = t186 * t264 - t222 * t236 + t332;
t147 = t185 * t236 - t186 * t235 + t336;
t146 = -t169 * t226 + t190 * t198 + t334;
t145 = t170 * t226 - t190 * t199 + t331;
t144 = t169 * t199 - t170 * t198 + t333;
t143 = -t153 * t226 - t161 * t209 + t172 * t177 + t178 * t198 + t334;
t142 = t154 * t226 + t162 * t209 - t173 * t177 - t178 * t199 + t331;
t141 = t153 * t199 - t154 * t198 + t161 * t173 - t162 * t172 + t333;
t1 = m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(3) * (t191 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + t172 * ((t156 * t237 + t158 * t200 + t160 * t201) * t173 + (t155 * t237 + t157 * t200 + t159 * t201) * t172 + (t174 * t237 + t175 * t200 + t176 * t201) * t209) / 0.2e1 + t198 * ((t164 * t237 + t166 * t204 + t168 * t205) * t199 + (t163 * t237 + t165 * t204 + t167 * t205) * t198 + (t187 * t237 + t188 * t204 + t189 * t205) * t226) / 0.2e1 + t173 * ((t156 * t239 + t158 * t202 + t160 * t203) * t173 + (t155 * t239 + t157 * t202 + t159 * t203) * t172 + (t174 * t239 + t175 * t202 + t176 * t203) * t209) / 0.2e1 + t199 * ((t164 * t239 + t166 * t206 + t168 * t207) * t199 + (t163 * t239 + t165 * t206 + t167 * t207) * t198 + (t187 * t239 + t188 * t206 + t189 * t207) * t226) / 0.2e1 + t235 * ((-t180 * t259 - t182 * t237 + t184 * t238) * t236 + (-t179 * t259 - t181 * t237 + t183 * t238) * t235 + (-t219 * t259 - t220 * t237 + t221 * t238) * t264) / 0.2e1 + t236 * ((-t180 * t261 - t182 * t239 + t184 * t240) * t236 + (-t179 * t261 - t181 * t239 + t183 * t240) * t235 + (-t219 * t261 - t220 * t239 + t221 * t240) * t264) / 0.2e1 + t209 * ((t156 * t265 + t158 * t227 + t160 * t228) * t173 + (t155 * t265 + t157 * t227 + t159 * t228) * t172 + (t174 * t265 + t175 * t227 + t176 * t228) * t209) / 0.2e1 + t226 * ((t164 * t265 + t166 * t229 + t168 * t230) * t199 + (t163 * t265 + t165 * t229 + t167 * t230) * t198 + (t187 * t265 + t188 * t229 + t189 * t230) * t226) / 0.2e1 + m(2) * (t267 ^ 2 + t268 ^ 2 + t269 ^ 2) / 0.2e1 + t264 * ((-t180 * t277 - t182 * t265 + t184 * t266) * t236 + (-t179 * t277 - t181 * t265 + t183 * t266) * t235 + (-t219 * t277 - t220 * t265 + t221 * t266) * t264) / 0.2e1 + m(1) * (t292 ^ 2 + t293 ^ 2 + t294 ^ 2) / 0.2e1 + ((t243 * t259 + t244 * t260 + t274 * t282 + t275 * t283 - t358 * t376) * t297 + (t214 * t259 + t216 * t260 + t249 * t282 + t251 * t283 - t358 * t377) * t296 + (t213 * t259 + t215 * t260 + t248 * t282 + t250 * t283 - t378 * t358) * t295) * t295 / 0.2e1 + ((t243 * t261 + t244 * t262 + t274 * t284 + t275 * t285 + t359 * t376) * t297 + (t214 * t261 + t216 * t262 + t249 * t284 + t251 * t285 + t377 * t359) * t296 + (t213 * t261 + t215 * t262 + t248 * t284 + t250 * t285 + t359 * t378) * t295) * t296 / 0.2e1 + ((t246 * t295 + t247 * t296 + t273 * t297) * t323 + ((t249 * t372 + t251 * t326) * t296 + (t248 * t372 + t250 * t326) * t295 + (t274 * t372 + t275 * t326) * t297) * t322 + (t212 * t323 + t214 * t277 + t216 * t278) * t296 + (t211 * t323 + t213 * t277 + t215 * t278) * t295 + (t242 * t323 + t243 * t277 + t244 * t278) * t297) * t297 / 0.2e1 + ((-t300 * t327 + t302 * t329 + Icges(1,4)) * V_base(5) + (-t301 * t327 + t303 * t329 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t300 * t329 + t302 * t327 + Icges(1,2)) * V_base(5) + (t301 * t329 + t303 * t327 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t327 + Icges(2,6) * t329) * V_base(5) + (Icges(2,5) * t329 - Icges(2,6) * t327) * V_base(4) + Icges(2,3) * t316 / 0.2e1) * t316;
T  = t1;
