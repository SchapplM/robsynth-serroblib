% Calculate kinetic energy for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:32
% EndTime: 2019-03-09 16:11:36
% DurationCPUTime: 4.97s
% Computational Cost: add. (3081->386), mult. (7256->539), div. (0->0), fcn. (8949->12), ass. (0->166)
t361 = Icges(5,1) + Icges(6,1);
t360 = -Icges(5,4) + Icges(6,5);
t359 = Icges(6,4) + Icges(5,5);
t358 = Icges(5,2) + Icges(6,3);
t357 = -Icges(5,6) + Icges(6,6);
t356 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t309 = sin(qJ(2));
t310 = sin(qJ(1));
t312 = cos(qJ(2));
t313 = cos(qJ(1));
t342 = cos(pkin(6));
t326 = t313 * t342;
t275 = t309 * t326 + t310 * t312;
t308 = sin(qJ(3));
t306 = sin(pkin(6));
t337 = t306 * t313;
t343 = cos(qJ(3));
t254 = t275 * t343 - t308 * t337;
t274 = t309 * t310 - t312 * t326;
t305 = sin(pkin(11));
t341 = cos(pkin(11));
t223 = t254 * t305 - t274 * t341;
t224 = t254 * t341 + t274 * t305;
t329 = t306 * t343;
t253 = t275 * t308 + t313 * t329;
t355 = t358 * t223 + t360 * t224 + t357 * t253;
t327 = t310 * t342;
t277 = -t309 * t327 + t313 * t312;
t339 = t306 * t310;
t256 = t277 * t343 + t308 * t339;
t276 = t313 * t309 + t312 * t327;
t225 = t256 * t305 - t276 * t341;
t226 = t256 * t341 + t276 * t305;
t255 = t277 * t308 - t310 * t329;
t354 = t358 * t225 + t360 * t226 + t357 * t255;
t353 = t360 * t223 + t361 * t224 + t359 * t253;
t352 = t360 * t225 + t361 * t226 + t359 * t255;
t273 = t308 * t342 + t309 * t329;
t338 = t306 * t312;
t249 = t273 * t305 + t338 * t341;
t250 = t273 * t341 - t305 * t338;
t272 = t306 * t308 * t309 - t342 * t343;
t351 = t358 * t249 + t360 * t250 + t357 * t272;
t350 = t360 * t249 + t361 * t250 + t359 * t272;
t349 = -Icges(4,4) * t254 - Icges(4,6) * t274 + t357 * t223 + t359 * t224 + t356 * t253;
t348 = -Icges(4,4) * t256 - Icges(4,6) * t276 + t357 * t225 + t359 * t226 + t356 * t255;
t347 = -Icges(4,4) * t273 + Icges(4,6) * t338 + t357 * t249 + t359 * t250 + t356 * t272;
t340 = Icges(2,4) * t310;
t190 = pkin(4) * t224 + qJ(5) * t223;
t217 = pkin(3) * t254 + qJ(4) * t253;
t336 = -t190 - t217;
t191 = pkin(4) * t226 + qJ(5) * t225;
t218 = pkin(3) * t256 + qJ(4) * t255;
t335 = -t191 - t218;
t216 = pkin(4) * t250 + qJ(5) * t249;
t242 = pkin(3) * t273 + qJ(4) * t272;
t334 = -t216 - t242;
t333 = qJD(2) * t306;
t332 = V_base(5) * pkin(7) + V_base(1);
t328 = t342 * pkin(8);
t286 = t310 * t333 + V_base(4);
t302 = V_base(6) + qJD(1);
t252 = qJD(3) * t276 + t286;
t287 = qJD(2) * t342 + t302;
t285 = -t313 * t333 + V_base(5);
t280 = t310 * pkin(1) - pkin(8) * t337;
t325 = -t280 * t302 + V_base(5) * t328 + t332;
t281 = pkin(1) * t313 + pkin(8) * t339;
t324 = V_base(4) * t280 - t281 * V_base(5) + V_base(3);
t251 = qJD(3) * t274 + t285;
t270 = -qJD(3) * t338 + t287;
t244 = pkin(2) * t275 + pkin(9) * t274;
t279 = (pkin(2) * t309 - pkin(9) * t312) * t306;
t323 = -t244 * t287 + t285 * t279 + t325;
t245 = pkin(2) * t277 + pkin(9) * t276;
t322 = t286 * t244 - t245 * t285 + t324;
t321 = t302 * t281 + V_base(2) + (-t328 - pkin(7)) * V_base(4);
t320 = qJD(4) * t255 + t251 * t242 + t323;
t319 = qJD(4) * t272 + t252 * t217 + t322;
t318 = qJD(5) * t225 + t251 * t216 + t320;
t317 = t287 * t245 - t286 * t279 + t321;
t316 = qJD(5) * t249 + t252 * t190 + t319;
t315 = qJD(4) * t253 + t270 * t218 + t317;
t314 = qJD(5) * t223 + t270 * t191 + t315;
t311 = cos(qJ(6));
t307 = sin(qJ(6));
t303 = Icges(2,4) * t313;
t295 = rSges(2,1) * t313 - t310 * rSges(2,2);
t294 = t310 * rSges(2,1) + rSges(2,2) * t313;
t293 = Icges(2,1) * t313 - t340;
t292 = Icges(2,1) * t310 + t303;
t291 = -Icges(2,2) * t310 + t303;
t290 = Icges(2,2) * t313 + t340;
t284 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t283 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t282 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t264 = t342 * rSges(3,3) + (rSges(3,1) * t309 + rSges(3,2) * t312) * t306;
t263 = Icges(3,5) * t342 + (Icges(3,1) * t309 + Icges(3,4) * t312) * t306;
t262 = Icges(3,6) * t342 + (Icges(3,4) * t309 + Icges(3,2) * t312) * t306;
t261 = Icges(3,3) * t342 + (Icges(3,5) * t309 + Icges(3,6) * t312) * t306;
t260 = V_base(5) * rSges(2,3) - t294 * t302 + t332;
t259 = t295 * t302 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t257 = t294 * V_base(4) - t295 * V_base(5) + V_base(3);
t243 = -qJD(6) * t272 + t270;
t241 = rSges(3,1) * t277 - rSges(3,2) * t276 + rSges(3,3) * t339;
t240 = t275 * rSges(3,1) - t274 * rSges(3,2) - rSges(3,3) * t337;
t239 = Icges(3,1) * t277 - Icges(3,4) * t276 + Icges(3,5) * t339;
t238 = Icges(3,1) * t275 - Icges(3,4) * t274 - Icges(3,5) * t337;
t237 = Icges(3,4) * t277 - Icges(3,2) * t276 + Icges(3,6) * t339;
t236 = Icges(3,4) * t275 - Icges(3,2) * t274 - Icges(3,6) * t337;
t235 = Icges(3,5) * t277 - Icges(3,6) * t276 + Icges(3,3) * t339;
t234 = Icges(3,5) * t275 - Icges(3,6) * t274 - Icges(3,3) * t337;
t233 = rSges(4,1) * t273 - rSges(4,2) * t272 - rSges(4,3) * t338;
t232 = Icges(4,1) * t273 - Icges(4,4) * t272 - Icges(4,5) * t338;
t230 = Icges(4,5) * t273 - Icges(4,6) * t272 - Icges(4,3) * t338;
t227 = pkin(5) * t250 - pkin(10) * t272;
t220 = -qJD(6) * t255 + t252;
t219 = -qJD(6) * t253 + t251;
t214 = t249 * t307 + t250 * t311;
t213 = t249 * t311 - t250 * t307;
t212 = rSges(4,1) * t256 - rSges(4,2) * t255 + rSges(4,3) * t276;
t211 = rSges(4,1) * t254 - rSges(4,2) * t253 + rSges(4,3) * t274;
t210 = Icges(4,1) * t256 - Icges(4,4) * t255 + Icges(4,5) * t276;
t209 = Icges(4,1) * t254 - Icges(4,4) * t253 + Icges(4,5) * t274;
t206 = Icges(4,5) * t256 - Icges(4,6) * t255 + Icges(4,3) * t276;
t205 = Icges(4,5) * t254 - Icges(4,6) * t253 + Icges(4,3) * t274;
t203 = rSges(5,1) * t250 - rSges(5,2) * t249 + rSges(5,3) * t272;
t202 = rSges(6,1) * t250 + rSges(6,2) * t272 + rSges(6,3) * t249;
t195 = pkin(5) * t226 - pkin(10) * t255;
t194 = pkin(5) * t224 - pkin(10) * t253;
t189 = t225 * t307 + t226 * t311;
t188 = t225 * t311 - t226 * t307;
t187 = t223 * t307 + t224 * t311;
t186 = t223 * t311 - t224 * t307;
t185 = -t240 * t287 + t264 * t285 + t325;
t184 = t287 * t241 - t286 * t264 + t321;
t181 = rSges(5,1) * t226 - rSges(5,2) * t225 + rSges(5,3) * t255;
t180 = rSges(6,1) * t226 + rSges(6,2) * t255 + rSges(6,3) * t225;
t179 = rSges(5,1) * t224 - rSges(5,2) * t223 + rSges(5,3) * t253;
t178 = rSges(6,1) * t224 + rSges(6,2) * t253 + rSges(6,3) * t223;
t165 = t240 * t286 - t241 * t285 + t324;
t164 = rSges(7,1) * t214 + rSges(7,2) * t213 - rSges(7,3) * t272;
t163 = Icges(7,1) * t214 + Icges(7,4) * t213 - Icges(7,5) * t272;
t162 = Icges(7,4) * t214 + Icges(7,2) * t213 - Icges(7,6) * t272;
t161 = Icges(7,5) * t214 + Icges(7,6) * t213 - Icges(7,3) * t272;
t160 = rSges(7,1) * t189 + rSges(7,2) * t188 - rSges(7,3) * t255;
t159 = rSges(7,1) * t187 + rSges(7,2) * t186 - rSges(7,3) * t253;
t158 = Icges(7,1) * t189 + Icges(7,4) * t188 - Icges(7,5) * t255;
t157 = Icges(7,1) * t187 + Icges(7,4) * t186 - Icges(7,5) * t253;
t156 = Icges(7,4) * t189 + Icges(7,2) * t188 - Icges(7,6) * t255;
t155 = Icges(7,4) * t187 + Icges(7,2) * t186 - Icges(7,6) * t253;
t154 = Icges(7,5) * t189 + Icges(7,6) * t188 - Icges(7,3) * t255;
t153 = Icges(7,5) * t187 + Icges(7,6) * t186 - Icges(7,3) * t253;
t152 = -t211 * t270 + t233 * t251 + t323;
t151 = t270 * t212 - t252 * t233 + t317;
t150 = t211 * t252 - t212 * t251 + t322;
t149 = t203 * t251 + (-t179 - t217) * t270 + t320;
t148 = t270 * t181 + (-t203 - t242) * t252 + t315;
t147 = t179 * t252 + (-t181 - t218) * t251 + t319;
t146 = t202 * t251 + (-t178 + t336) * t270 + t318;
t145 = t270 * t180 + (-t202 + t334) * t252 + t314;
t144 = t178 * t252 + (-t180 + t335) * t251 + t316;
t143 = t318 - t159 * t243 + t164 * t219 + t227 * t251 + (-t194 + t336) * t270;
t142 = t243 * t160 - t220 * t164 + t270 * t195 + (-t227 + t334) * t252 + t314;
t141 = t159 * t220 - t160 * t219 + t194 * t252 + (-t195 + t335) * t251 + t316;
t1 = t286 * ((t235 * t339 - t237 * t276 + t239 * t277) * t286 + (t234 * t339 - t236 * t276 + t238 * t277) * t285 + (t261 * t339 - t262 * t276 + t263 * t277) * t287) / 0.2e1 + t285 * ((-t235 * t337 - t274 * t237 + t275 * t239) * t286 + (-t234 * t337 - t274 * t236 + t275 * t238) * t285 + (-t261 * t337 - t274 * t262 + t275 * t263) * t287) / 0.2e1 + t287 * (((t237 * t312 + t239 * t309) * t286 + (t236 * t312 + t238 * t309) * t285 + (t262 * t312 + t263 * t309) * t287) * t306 + (t234 * t285 + t235 * t286 + t261 * t287) * t342) / 0.2e1 + m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(3) * (t165 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + t219 * ((-t154 * t253 + t156 * t186 + t158 * t187) * t220 + (-t153 * t253 + t155 * t186 + t157 * t187) * t219 + (-t161 * t253 + t162 * t186 + t163 * t187) * t243) / 0.2e1 + t220 * ((-t154 * t255 + t156 * t188 + t158 * t189) * t220 + (-t153 * t255 + t155 * t188 + t157 * t189) * t219 + (-t161 * t255 + t162 * t188 + t163 * t189) * t243) / 0.2e1 + m(2) * (t257 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + t243 * ((-t154 * t272 + t156 * t213 + t158 * t214) * t220 + (-t153 * t272 + t155 * t213 + t157 * t214) * t219 + (-t161 * t272 + t162 * t213 + t163 * t214) * t243) / 0.2e1 + m(1) * (t282 ^ 2 + t283 ^ 2 + t284 ^ 2) / 0.2e1 + ((-t310 * t290 + t292 * t313 + Icges(1,4)) * V_base(5) + (-t310 * t291 + t293 * t313 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t290 * t313 + t310 * t292 + Icges(1,2)) * V_base(5) + (t291 * t313 + t310 * t293 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t223 * t351 + t224 * t350 + t230 * t274 + t232 * t254 + t253 * t347) * t270 + (t206 * t274 + t210 * t254 + t223 * t354 + t224 * t352 + t253 * t348) * t252 + (t205 * t274 + t209 * t254 + t355 * t223 + t353 * t224 + t349 * t253) * t251) * t251 / 0.2e1 + ((t225 * t351 + t226 * t350 + t230 * t276 + t232 * t256 + t255 * t347) * t270 + (t206 * t276 + t210 * t256 + t354 * t225 + t352 * t226 + t348 * t255) * t252 + (t205 * t276 + t209 * t256 + t225 * t355 + t353 * t226 + t349 * t255) * t251) * t252 / 0.2e1 + ((-t230 * t338 + t232 * t273 + t351 * t249 + t350 * t250 + t347 * t272) * t270 + (-t206 * t338 + t210 * t273 + t249 * t354 + t250 * t352 + t272 * t348) * t252 + (-t205 * t338 + t209 * t273 + t249 * t355 + t353 * t250 + t349 * t272) * t251) * t270 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t310 + Icges(2,6) * t313) * V_base(5) + (Icges(2,5) * t313 - Icges(2,6) * t310) * V_base(4) + Icges(2,3) * t302 / 0.2e1) * t302;
T  = t1;
