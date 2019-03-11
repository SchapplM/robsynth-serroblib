% Calculate kinetic energy for
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:21:09
% EndTime: 2019-03-09 10:21:16
% DurationCPUTime: 6.54s
% Computational Cost: add. (3877->436), mult. (7918->616), div. (0->0), fcn. (9820->14), ass. (0->191)
t394 = -Icges(5,3) - Icges(6,3);
t335 = sin(qJ(2));
t375 = sin(pkin(11));
t377 = cos(pkin(6));
t354 = t377 * t375;
t376 = cos(pkin(11));
t355 = t377 * t376;
t381 = cos(qJ(2));
t290 = t335 * t355 + t354 * t381;
t299 = -t335 * t375 + t381 * t376;
t336 = sin(qJ(1));
t339 = cos(qJ(1));
t269 = t290 * t339 + t336 * t299;
t369 = qJ(4) + pkin(12);
t327 = sin(t369);
t331 = sin(pkin(6));
t358 = cos(t369);
t353 = t331 * t358;
t236 = t269 * t327 + t339 * t353;
t372 = t331 * t339;
t237 = t269 * t358 - t327 * t372;
t334 = sin(qJ(4));
t338 = cos(qJ(4));
t244 = -t269 * t334 - t338 * t372;
t364 = t334 * t372;
t245 = t269 * t338 - t364;
t298 = -t335 * t376 - t375 * t381;
t345 = -t335 * t354 + t355 * t381;
t268 = t336 * t298 + t339 * t345;
t393 = Icges(5,5) * t245 + Icges(6,5) * t237 + Icges(5,6) * t244 - Icges(6,6) * t236 + t394 * t268;
t271 = -t336 * t290 + t299 * t339;
t238 = t271 * t327 - t336 * t353;
t373 = t331 * t336;
t239 = t271 * t358 + t327 * t373;
t246 = -t271 * t334 + t338 * t373;
t365 = t334 * t373;
t247 = t271 * t338 + t365;
t270 = t298 * t339 - t336 * t345;
t392 = Icges(5,5) * t247 + Icges(6,5) * t239 + Icges(5,6) * t246 - Icges(6,6) * t238 + t394 * t270;
t289 = t298 * t331;
t272 = -t289 * t327 - t358 * t377;
t273 = -t289 * t358 + t327 * t377;
t275 = t289 * t334 + t338 * t377;
t359 = t377 * t334;
t276 = -t289 * t338 + t359;
t288 = t299 * t331;
t391 = Icges(5,5) * t276 + Icges(6,5) * t273 + Icges(5,6) * t275 - Icges(6,6) * t272 + t394 * t288;
t218 = Icges(4,5) * t269 + Icges(4,6) * t268 - Icges(4,3) * t372;
t357 = t377 * t381;
t293 = -t336 * t335 + t339 * t357;
t361 = t335 * t377;
t294 = t336 * t381 + t339 * t361;
t253 = Icges(3,5) * t294 + Icges(3,6) * t293 - Icges(3,3) * t372;
t390 = t218 + t253;
t219 = Icges(4,5) * t271 + Icges(4,6) * t270 + Icges(4,3) * t373;
t295 = -t339 * t335 - t336 * t357;
t296 = -t336 * t361 + t339 * t381;
t254 = Icges(3,5) * t296 + Icges(3,6) * t295 + Icges(3,3) * t373;
t389 = t219 + t254;
t249 = -Icges(4,5) * t289 + Icges(4,6) * t288 + Icges(4,3) * t377;
t284 = Icges(3,3) * t377 + (Icges(3,5) * t335 + Icges(3,6) * t381) * t331;
t388 = t249 + t284;
t380 = pkin(2) * t381;
t379 = pkin(4) * t338;
t374 = Icges(2,4) * t336;
t371 = qJD(2) * t331;
t370 = qJD(3) * t331;
t368 = V_base(5) * pkin(7) + V_base(1);
t363 = t377 * pkin(8);
t308 = t336 * t371 + V_base(4);
t328 = V_base(6) + qJD(1);
t362 = pkin(2) * t361 - qJ(3) * t331;
t243 = -qJD(4) * t270 + t308;
t309 = qJD(2) * t377 + t328;
t274 = -qJD(4) * t288 + t309;
t307 = -t339 * t371 + V_base(5);
t301 = t336 * pkin(1) - pkin(8) * t372;
t352 = -t301 * t328 + V_base(5) * t363 + t368;
t302 = pkin(1) * t339 + pkin(8) * t373;
t351 = V_base(4) * t301 - t302 * V_base(5) + V_base(3);
t242 = -qJD(4) * t268 + t307;
t300 = t331 * t335 * pkin(2) + qJ(3) * t377;
t350 = t307 * t300 + t336 * t370 + t352;
t266 = t336 * t380 + t339 * t362;
t349 = qJD(3) * t377 + t308 * t266 + t351;
t348 = t328 * t302 + V_base(2) + (-t363 - pkin(7)) * V_base(4);
t230 = t269 * pkin(3) - t268 * pkin(9);
t261 = -t289 * pkin(3) - t288 * pkin(9);
t347 = t307 * t261 + (-t230 - t266) * t309 + t350;
t231 = pkin(3) * t271 - pkin(9) * t270;
t267 = -t336 * t362 + t339 * t380;
t346 = t308 * t230 + (-t231 - t267) * t307 + t349;
t344 = t309 * t267 - t339 * t370 + t348;
t211 = pkin(4) * t359 - qJ(5) * t288 - t289 * t379;
t343 = -qJD(5) * t270 + t242 * t211 + t347;
t177 = -pkin(4) * t364 - qJ(5) * t268 + t269 * t379;
t342 = -qJD(5) * t288 + t243 * t177 + t346;
t341 = t309 * t231 + (-t261 - t300) * t308 + t344;
t178 = pkin(4) * t365 - qJ(5) * t270 + t271 * t379;
t340 = -qJD(5) * t268 + t274 * t178 + t341;
t337 = cos(qJ(6));
t333 = sin(qJ(6));
t329 = Icges(2,4) * t339;
t317 = rSges(2,1) * t339 - t336 * rSges(2,2);
t316 = t336 * rSges(2,1) + rSges(2,2) * t339;
t315 = Icges(2,1) * t339 - t374;
t314 = Icges(2,1) * t336 + t329;
t313 = -Icges(2,2) * t336 + t329;
t312 = Icges(2,2) * t339 + t374;
t305 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t304 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t303 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t287 = t377 * rSges(3,3) + (rSges(3,1) * t335 + rSges(3,2) * t381) * t331;
t286 = Icges(3,5) * t377 + (Icges(3,1) * t335 + Icges(3,4) * t381) * t331;
t285 = Icges(3,6) * t377 + (Icges(3,4) * t335 + Icges(3,2) * t381) * t331;
t279 = V_base(5) * rSges(2,3) - t316 * t328 + t368;
t278 = t317 * t328 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t277 = t316 * V_base(4) - t317 * V_base(5) + V_base(3);
t260 = rSges(3,1) * t296 + rSges(3,2) * t295 + rSges(3,3) * t373;
t259 = t294 * rSges(3,1) + t293 * rSges(3,2) - rSges(3,3) * t372;
t258 = Icges(3,1) * t296 + Icges(3,4) * t295 + Icges(3,5) * t373;
t257 = Icges(3,1) * t294 + Icges(3,4) * t293 - Icges(3,5) * t372;
t256 = Icges(3,4) * t296 + Icges(3,2) * t295 + Icges(3,6) * t373;
t255 = Icges(3,4) * t294 + Icges(3,2) * t293 - Icges(3,6) * t372;
t252 = -t289 * rSges(4,1) + t288 * rSges(4,2) + rSges(4,3) * t377;
t251 = -Icges(4,1) * t289 + Icges(4,4) * t288 + Icges(4,5) * t377;
t250 = -Icges(4,4) * t289 + Icges(4,2) * t288 + Icges(4,6) * t377;
t235 = t273 * t337 - t288 * t333;
t234 = -t273 * t333 - t288 * t337;
t233 = qJD(6) * t272 + t274;
t232 = pkin(5) * t273 + pkin(10) * t272;
t229 = rSges(5,1) * t276 + rSges(5,2) * t275 - rSges(5,3) * t288;
t228 = Icges(5,1) * t276 + Icges(5,4) * t275 - Icges(5,5) * t288;
t227 = Icges(5,4) * t276 + Icges(5,2) * t275 - Icges(5,6) * t288;
t225 = rSges(4,1) * t271 + rSges(4,2) * t270 + rSges(4,3) * t373;
t224 = t269 * rSges(4,1) + t268 * rSges(4,2) - rSges(4,3) * t372;
t223 = Icges(4,1) * t271 + Icges(4,4) * t270 + Icges(4,5) * t373;
t222 = Icges(4,1) * t269 + Icges(4,4) * t268 - Icges(4,5) * t372;
t221 = Icges(4,4) * t271 + Icges(4,2) * t270 + Icges(4,6) * t373;
t220 = Icges(4,4) * t269 + Icges(4,2) * t268 - Icges(4,6) * t372;
t217 = rSges(6,1) * t273 - rSges(6,2) * t272 - rSges(6,3) * t288;
t216 = Icges(6,1) * t273 - Icges(6,4) * t272 - Icges(6,5) * t288;
t215 = Icges(6,4) * t273 - Icges(6,2) * t272 - Icges(6,6) * t288;
t210 = t239 * t337 - t270 * t333;
t209 = -t239 * t333 - t270 * t337;
t208 = t237 * t337 - t268 * t333;
t207 = -t237 * t333 - t268 * t337;
t206 = qJD(6) * t238 + t243;
t205 = qJD(6) * t236 + t242;
t204 = -t259 * t309 + t287 * t307 + t352;
t203 = t309 * t260 - t308 * t287 + t348;
t202 = pkin(5) * t239 + pkin(10) * t238;
t201 = pkin(5) * t237 + pkin(10) * t236;
t200 = t259 * t308 - t260 * t307 + t351;
t199 = rSges(5,1) * t247 + rSges(5,2) * t246 - rSges(5,3) * t270;
t198 = rSges(5,1) * t245 + rSges(5,2) * t244 - rSges(5,3) * t268;
t196 = Icges(5,1) * t247 + Icges(5,4) * t246 - Icges(5,5) * t270;
t195 = Icges(5,1) * t245 + Icges(5,4) * t244 - Icges(5,5) * t268;
t194 = Icges(5,4) * t247 + Icges(5,2) * t246 - Icges(5,6) * t270;
t193 = Icges(5,4) * t245 + Icges(5,2) * t244 - Icges(5,6) * t268;
t190 = rSges(6,1) * t239 - rSges(6,2) * t238 - rSges(6,3) * t270;
t189 = rSges(6,1) * t237 - rSges(6,2) * t236 - rSges(6,3) * t268;
t188 = Icges(6,1) * t239 - Icges(6,4) * t238 - Icges(6,5) * t270;
t187 = Icges(6,1) * t237 - Icges(6,4) * t236 - Icges(6,5) * t268;
t186 = Icges(6,4) * t239 - Icges(6,2) * t238 - Icges(6,6) * t270;
t185 = Icges(6,4) * t237 - Icges(6,2) * t236 - Icges(6,6) * t268;
t182 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t272;
t181 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t272;
t180 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t272;
t179 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t272;
t174 = rSges(7,1) * t210 + rSges(7,2) * t209 + rSges(7,3) * t238;
t173 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t236;
t172 = Icges(7,1) * t210 + Icges(7,4) * t209 + Icges(7,5) * t238;
t171 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t236;
t170 = Icges(7,4) * t210 + Icges(7,2) * t209 + Icges(7,6) * t238;
t169 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t236;
t168 = Icges(7,5) * t210 + Icges(7,6) * t209 + Icges(7,3) * t238;
t167 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t236;
t166 = t252 * t307 + (-t224 - t266) * t309 + t350;
t165 = t309 * t225 + (-t252 - t300) * t308 + t344;
t164 = t224 * t308 + (-t225 - t267) * t307 + t349;
t163 = -t198 * t274 + t229 * t242 + t347;
t162 = t274 * t199 - t243 * t229 + t341;
t161 = t198 * t243 - t199 * t242 + t346;
t160 = t217 * t242 + (-t177 - t189) * t274 + t343;
t159 = t274 * t190 + (-t211 - t217) * t243 + t340;
t158 = t189 * t243 + (-t178 - t190) * t242 + t342;
t157 = (-t177 - t201) * t274 - t173 * t233 + t182 * t205 + t232 * t242 + t343;
t156 = t233 * t174 - t206 * t182 + t274 * t202 + (-t211 - t232) * t243 + t340;
t155 = t173 * t206 - t174 * t205 + t201 * t243 + (-t178 - t202) * t242 + t342;
t1 = m(3) * (t200 ^ 2 + t203 ^ 2 + t204 ^ 2) / 0.2e1 + m(4) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(6) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(5) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(7) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + t205 * ((t168 * t236 + t170 * t207 + t172 * t208) * t206 + (t236 * t167 + t207 * t169 + t208 * t171) * t205 + (t179 * t236 + t180 * t207 + t181 * t208) * t233) / 0.2e1 + t206 * ((t238 * t168 + t209 * t170 + t210 * t172) * t206 + (t167 * t238 + t169 * t209 + t171 * t210) * t205 + (t179 * t238 + t180 * t209 + t181 * t210) * t233) / 0.2e1 + t233 * ((t168 * t272 + t170 * t234 + t172 * t235) * t206 + (t167 * t272 + t169 * t234 + t171 * t235) * t205 + (t272 * t179 + t234 * t180 + t235 * t181) * t233) / 0.2e1 + m(2) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + m(1) * (t303 ^ 2 + t304 ^ 2 + t305 ^ 2) / 0.2e1 + ((-t215 * t236 + t216 * t237 + t227 * t244 + t228 * t245 - t268 * t391) * t274 + (-t186 * t236 + t188 * t237 + t194 * t244 + t196 * t245 - t268 * t392) * t243 + (-t236 * t185 + t237 * t187 + t244 * t193 + t245 * t195 - t268 * t393) * t242) * t242 / 0.2e1 + ((-t215 * t238 + t216 * t239 + t227 * t246 + t228 * t247 - t270 * t391) * t274 + (-t238 * t186 + t239 * t188 + t246 * t194 + t247 * t196 - t270 * t392) * t243 + (-t185 * t238 + t187 * t239 + t193 * t246 + t195 * t247 - t270 * t393) * t242) * t243 / 0.2e1 + ((-t272 * t215 + t273 * t216 + t275 * t227 + t276 * t228 - t288 * t391) * t274 + (-t186 * t272 + t188 * t273 + t194 * t275 + t196 * t276 - t288 * t392) * t243 + (-t185 * t272 + t187 * t273 + t193 * t275 + t195 * t276 - t288 * t393) * t242) * t274 / 0.2e1 + ((t268 * t250 + t269 * t251 + t293 * t285 + t294 * t286 - t372 * t388) * t309 + (t268 * t221 + t269 * t223 + t293 * t256 + t294 * t258 - t372 * t389) * t308 + (t268 * t220 + t269 * t222 + t293 * t255 + t294 * t257 - t372 * t390) * t307) * t307 / 0.2e1 + ((t250 * t270 + t251 * t271 + t285 * t295 + t286 * t296 + t373 * t388) * t309 + (t270 * t221 + t271 * t223 + t295 * t256 + t296 * t258 + t373 * t389) * t308 + (t220 * t270 + t222 * t271 + t255 * t295 + t257 * t296 + t373 * t390) * t307) * t308 / 0.2e1 + ((t219 * t377 + t288 * t221 - t289 * t223) * t308 + (t218 * t377 + t288 * t220 - t289 * t222) * t307 + (t249 * t377 + t288 * t250 - t289 * t251) * t309 + ((t256 * t381 + t258 * t335) * t308 + (t255 * t381 + t257 * t335) * t307 + (t285 * t381 + t286 * t335) * t309) * t331 + (t253 * t307 + t254 * t308 + t284 * t309) * t377) * t309 / 0.2e1 + ((-t336 * t312 + t314 * t339 + Icges(1,4)) * V_base(5) + (-t336 * t313 + t339 * t315 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t339 * t312 + t336 * t314 + Icges(1,2)) * V_base(5) + (t313 * t339 + t336 * t315 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t336 + Icges(2,6) * t339) * V_base(5) + (Icges(2,5) * t339 - Icges(2,6) * t336) * V_base(4) + Icges(2,3) * t328 / 0.2e1) * t328;
T  = t1;
