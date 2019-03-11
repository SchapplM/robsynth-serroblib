% Calculate kinetic energy for
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:50:58
% EndTime: 2019-03-09 23:51:04
% DurationCPUTime: 6.11s
% Computational Cost: add. (3249->392), mult. (7634->569), div. (0->0), fcn. (9432->12), ass. (0->169)
t361 = Icges(5,1) + Icges(6,1);
t360 = -Icges(5,4) + Icges(6,5);
t359 = Icges(6,4) + Icges(5,5);
t358 = Icges(5,2) + Icges(6,3);
t357 = Icges(6,2) + Icges(5,3);
t356 = -Icges(5,6) + Icges(6,6);
t312 = sin(qJ(2));
t313 = sin(qJ(1));
t315 = cos(qJ(2));
t316 = cos(qJ(1));
t341 = cos(pkin(6));
t329 = t316 * t341;
t278 = t312 * t329 + t313 * t315;
t311 = sin(qJ(3));
t308 = sin(pkin(6));
t337 = t308 * t316;
t343 = cos(qJ(3));
t257 = t278 * t343 - t311 * t337;
t277 = t312 * t313 - t315 * t329;
t310 = sin(qJ(4));
t342 = cos(qJ(4));
t226 = t257 * t310 - t277 * t342;
t227 = t257 * t342 + t277 * t310;
t332 = t308 * t343;
t256 = t278 * t311 + t316 * t332;
t355 = t358 * t226 + t360 * t227 + t356 * t256;
t330 = t313 * t341;
t280 = -t312 * t330 + t316 * t315;
t339 = t308 * t313;
t259 = t280 * t343 + t311 * t339;
t279 = t316 * t312 + t315 * t330;
t228 = t259 * t310 - t279 * t342;
t229 = t259 * t342 + t279 * t310;
t258 = t280 * t311 - t313 * t332;
t354 = t358 * t228 + t360 * t229 + t356 * t258;
t353 = t356 * t226 + t359 * t227 + t357 * t256;
t352 = t356 * t228 + t359 * t229 + t357 * t258;
t351 = t360 * t226 + t361 * t227 + t359 * t256;
t350 = t360 * t228 + t361 * t229 + t359 * t258;
t276 = t311 * t341 + t312 * t332;
t338 = t308 * t315;
t252 = t276 * t310 + t338 * t342;
t253 = t276 * t342 - t310 * t338;
t275 = t308 * t311 * t312 - t341 * t343;
t349 = t358 * t252 + t360 * t253 + t356 * t275;
t348 = t356 * t252 + t359 * t253 + t357 * t275;
t347 = t360 * t252 + t361 * t253 + t359 * t275;
t340 = Icges(2,4) * t313;
t336 = qJD(2) * t308;
t335 = V_base(5) * pkin(7) + V_base(1);
t331 = t341 * pkin(8);
t289 = t313 * t336 + V_base(4);
t305 = V_base(6) + qJD(1);
t255 = qJD(3) * t279 + t289;
t290 = qJD(2) * t341 + t305;
t222 = qJD(4) * t258 + t255;
t288 = -t316 * t336 + V_base(5);
t283 = t313 * pkin(1) - pkin(8) * t337;
t328 = -t283 * t305 + V_base(5) * t331 + t335;
t284 = pkin(1) * t316 + pkin(8) * t339;
t327 = V_base(4) * t283 - t284 * V_base(5) + V_base(3);
t254 = qJD(3) * t277 + t288;
t221 = qJD(4) * t256 + t254;
t273 = -qJD(3) * t338 + t290;
t246 = qJD(4) * t275 + t273;
t247 = pkin(2) * t278 + pkin(9) * t277;
t282 = (pkin(2) * t312 - pkin(9) * t315) * t308;
t326 = -t247 * t290 + t288 * t282 + t328;
t248 = pkin(2) * t280 + pkin(9) * t279;
t325 = t289 * t247 - t248 * t288 + t327;
t324 = t305 * t284 + V_base(2) + (-t331 - pkin(7)) * V_base(4);
t219 = pkin(3) * t257 + pkin(10) * t256;
t245 = pkin(3) * t276 + pkin(10) * t275;
t323 = -t219 * t273 + t254 * t245 + t326;
t220 = pkin(3) * t259 + pkin(10) * t258;
t322 = t255 * t219 - t220 * t254 + t325;
t321 = t290 * t248 - t289 * t282 + t324;
t218 = pkin(4) * t253 + qJ(5) * t252;
t320 = qJD(5) * t228 + t221 * t218 + t323;
t191 = pkin(4) * t227 + qJ(5) * t226;
t319 = qJD(5) * t252 + t222 * t191 + t322;
t318 = t273 * t220 - t255 * t245 + t321;
t192 = pkin(4) * t229 + qJ(5) * t228;
t317 = qJD(5) * t226 + t246 * t192 + t318;
t314 = cos(qJ(6));
t309 = sin(qJ(6));
t306 = Icges(2,4) * t316;
t298 = rSges(2,1) * t316 - t313 * rSges(2,2);
t297 = t313 * rSges(2,1) + rSges(2,2) * t316;
t296 = Icges(2,1) * t316 - t340;
t295 = Icges(2,1) * t313 + t306;
t294 = -Icges(2,2) * t313 + t306;
t293 = Icges(2,2) * t316 + t340;
t287 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t286 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t285 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t267 = t341 * rSges(3,3) + (rSges(3,1) * t312 + rSges(3,2) * t315) * t308;
t266 = Icges(3,5) * t341 + (Icges(3,1) * t312 + Icges(3,4) * t315) * t308;
t265 = Icges(3,6) * t341 + (Icges(3,4) * t312 + Icges(3,2) * t315) * t308;
t264 = Icges(3,3) * t341 + (Icges(3,5) * t312 + Icges(3,6) * t315) * t308;
t263 = V_base(5) * rSges(2,3) - t297 * t305 + t335;
t262 = t298 * t305 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t260 = t297 * V_base(4) - t298 * V_base(5) + V_base(3);
t244 = rSges(3,1) * t280 - rSges(3,2) * t279 + rSges(3,3) * t339;
t243 = t278 * rSges(3,1) - t277 * rSges(3,2) - rSges(3,3) * t337;
t242 = Icges(3,1) * t280 - Icges(3,4) * t279 + Icges(3,5) * t339;
t241 = Icges(3,1) * t278 - Icges(3,4) * t277 - Icges(3,5) * t337;
t240 = Icges(3,4) * t280 - Icges(3,2) * t279 + Icges(3,6) * t339;
t239 = Icges(3,4) * t278 - Icges(3,2) * t277 - Icges(3,6) * t337;
t238 = Icges(3,5) * t280 - Icges(3,6) * t279 + Icges(3,3) * t339;
t237 = Icges(3,5) * t278 - Icges(3,6) * t277 - Icges(3,3) * t337;
t236 = rSges(4,1) * t276 - rSges(4,2) * t275 - rSges(4,3) * t338;
t235 = Icges(4,1) * t276 - Icges(4,4) * t275 - Icges(4,5) * t338;
t234 = Icges(4,4) * t276 - Icges(4,2) * t275 - Icges(4,6) * t338;
t233 = Icges(4,5) * t276 - Icges(4,6) * t275 - Icges(4,3) * t338;
t230 = pkin(5) * t253 - pkin(11) * t275;
t223 = -qJD(6) * t275 + t246;
t216 = t252 * t309 + t253 * t314;
t215 = t252 * t314 - t253 * t309;
t214 = rSges(4,1) * t259 - rSges(4,2) * t258 + rSges(4,3) * t279;
t213 = rSges(4,1) * t257 - rSges(4,2) * t256 + rSges(4,3) * t277;
t211 = Icges(4,1) * t259 - Icges(4,4) * t258 + Icges(4,5) * t279;
t210 = Icges(4,1) * t257 - Icges(4,4) * t256 + Icges(4,5) * t277;
t209 = Icges(4,4) * t259 - Icges(4,2) * t258 + Icges(4,6) * t279;
t208 = Icges(4,4) * t257 - Icges(4,2) * t256 + Icges(4,6) * t277;
t207 = Icges(4,5) * t259 - Icges(4,6) * t258 + Icges(4,3) * t279;
t206 = Icges(4,5) * t257 - Icges(4,6) * t256 + Icges(4,3) * t277;
t205 = rSges(5,1) * t253 - rSges(5,2) * t252 + rSges(5,3) * t275;
t204 = rSges(6,1) * t253 + rSges(6,2) * t275 + rSges(6,3) * t252;
t197 = pkin(5) * t229 - pkin(11) * t258;
t196 = pkin(5) * t227 - pkin(11) * t256;
t195 = -qJD(6) * t258 + t222;
t194 = -qJD(6) * t256 + t221;
t190 = t228 * t309 + t229 * t314;
t189 = t228 * t314 - t229 * t309;
t188 = t226 * t309 + t227 * t314;
t187 = t226 * t314 - t227 * t309;
t186 = -t243 * t290 + t267 * t288 + t328;
t185 = t290 * t244 - t289 * t267 + t324;
t183 = rSges(5,1) * t229 - rSges(5,2) * t228 + rSges(5,3) * t258;
t182 = rSges(6,1) * t229 + rSges(6,2) * t258 + rSges(6,3) * t228;
t181 = rSges(5,1) * t227 - rSges(5,2) * t226 + rSges(5,3) * t256;
t180 = rSges(6,1) * t227 + rSges(6,2) * t256 + rSges(6,3) * t226;
t167 = t243 * t289 - t244 * t288 + t327;
t165 = rSges(7,1) * t216 + rSges(7,2) * t215 - rSges(7,3) * t275;
t164 = Icges(7,1) * t216 + Icges(7,4) * t215 - Icges(7,5) * t275;
t163 = Icges(7,4) * t216 + Icges(7,2) * t215 - Icges(7,6) * t275;
t162 = Icges(7,5) * t216 + Icges(7,6) * t215 - Icges(7,3) * t275;
t160 = rSges(7,1) * t190 + rSges(7,2) * t189 - rSges(7,3) * t258;
t159 = rSges(7,1) * t188 + rSges(7,2) * t187 - rSges(7,3) * t256;
t158 = Icges(7,1) * t190 + Icges(7,4) * t189 - Icges(7,5) * t258;
t157 = Icges(7,1) * t188 + Icges(7,4) * t187 - Icges(7,5) * t256;
t156 = Icges(7,4) * t190 + Icges(7,2) * t189 - Icges(7,6) * t258;
t155 = Icges(7,4) * t188 + Icges(7,2) * t187 - Icges(7,6) * t256;
t154 = Icges(7,5) * t190 + Icges(7,6) * t189 - Icges(7,3) * t258;
t153 = Icges(7,5) * t188 + Icges(7,6) * t187 - Icges(7,3) * t256;
t152 = -t213 * t273 + t236 * t254 + t326;
t151 = t273 * t214 - t255 * t236 + t321;
t150 = t213 * t255 - t214 * t254 + t325;
t149 = -t181 * t246 + t205 * t221 + t323;
t148 = t246 * t183 - t222 * t205 + t318;
t147 = t181 * t222 - t183 * t221 + t322;
t146 = t204 * t221 + (-t180 - t191) * t246 + t320;
t145 = t246 * t182 + (-t204 - t218) * t222 + t317;
t144 = t180 * t222 + (-t182 - t192) * t221 + t319;
t143 = -t159 * t223 + t165 * t194 + t221 * t230 + t320 + (-t191 - t196) * t246;
t142 = t223 * t160 - t195 * t165 + t246 * t197 + (-t218 - t230) * t222 + t317;
t141 = t159 * t195 - t160 * t194 + t196 * t222 + (-t192 - t197) * t221 + t319;
t1 = t288 * ((-t238 * t337 - t277 * t240 + t278 * t242) * t289 + (-t237 * t337 - t277 * t239 + t278 * t241) * t288 + (-t264 * t337 - t277 * t265 + t278 * t266) * t290) / 0.2e1 + t273 * ((-t207 * t338 - t209 * t275 + t211 * t276) * t255 + (-t206 * t338 - t208 * t275 + t210 * t276) * t254 + (-t233 * t338 - t234 * t275 + t235 * t276) * t273) / 0.2e1 + t289 * ((t238 * t339 - t240 * t279 + t242 * t280) * t289 + (t237 * t339 - t239 * t279 + t241 * t280) * t288 + (t264 * t339 - t265 * t279 + t266 * t280) * t290) / 0.2e1 + t255 * ((t207 * t279 - t209 * t258 + t211 * t259) * t255 + (t206 * t279 - t208 * t258 + t210 * t259) * t254 + (t233 * t279 - t234 * t258 + t235 * t259) * t273) / 0.2e1 + m(1) * (t285 ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1 + t223 * ((-t154 * t275 + t156 * t215 + t158 * t216) * t195 + (-t153 * t275 + t155 * t215 + t157 * t216) * t194 + (-t275 * t162 + t215 * t163 + t216 * t164) * t223) / 0.2e1 + t254 * ((t207 * t277 - t209 * t256 + t211 * t257) * t255 + (t277 * t206 - t208 * t256 + t257 * t210) * t254 + (t233 * t277 - t234 * t256 + t235 * t257) * t273) / 0.2e1 + m(2) * (t260 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + t195 * ((-t258 * t154 + t189 * t156 + t190 * t158) * t195 + (-t153 * t258 + t155 * t189 + t157 * t190) * t194 + (-t162 * t258 + t163 * t189 + t164 * t190) * t223) / 0.2e1 + t194 * ((-t154 * t256 + t156 * t187 + t158 * t188) * t195 + (-t256 * t153 + t187 * t155 + t188 * t157) * t194 + (-t162 * t256 + t163 * t187 + t164 * t188) * t223) / 0.2e1 + m(3) * (t167 ^ 2 + t185 ^ 2 + t186 ^ 2) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + t290 * (((t240 * t315 + t242 * t312) * t289 + (t239 * t315 + t241 * t312) * t288 + (t265 * t315 + t266 * t312) * t290) * t308 + (t237 * t288 + t238 * t289 + t264 * t290) * t341) / 0.2e1 + ((t226 * t349 + t227 * t347 + t256 * t348) * t246 + (t226 * t354 + t227 * t350 + t256 * t352) * t222 + (t355 * t226 + t351 * t227 + t353 * t256) * t221) * t221 / 0.2e1 + ((t228 * t349 + t229 * t347 + t258 * t348) * t246 + (t354 * t228 + t350 * t229 + t352 * t258) * t222 + (t228 * t355 + t351 * t229 + t353 * t258) * t221) * t222 / 0.2e1 + ((t349 * t252 + t347 * t253 + t348 * t275) * t246 + (t252 * t354 + t253 * t350 + t275 * t352) * t222 + (t252 * t355 + t351 * t253 + t353 * t275) * t221) * t246 / 0.2e1 + ((-t313 * t293 + t295 * t316 + Icges(1,4)) * V_base(5) + (-t313 * t294 + t296 * t316 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t293 * t316 + t313 * t295 + Icges(1,2)) * V_base(5) + (t294 * t316 + t313 * t296 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t313 + Icges(2,6) * t316) * V_base(5) + (Icges(2,5) * t316 - Icges(2,6) * t313) * V_base(4) + Icges(2,3) * t305 / 0.2e1) * t305;
T  = t1;
