% Calculate kinetic energy for
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:25
% EndTime: 2019-03-09 08:54:30
% DurationCPUTime: 5.72s
% Computational Cost: add. (3787->431), mult. (7711->608), div. (0->0), fcn. (9550->14), ass. (0->189)
t395 = -Icges(4,2) - Icges(5,3);
t394 = Icges(4,3) + Icges(3,3);
t336 = sin(qJ(2));
t379 = sin(pkin(11));
t380 = cos(pkin(11));
t384 = cos(qJ(2));
t298 = -t336 * t380 - t379 * t384;
t337 = sin(qJ(1));
t339 = cos(qJ(1));
t381 = cos(pkin(6));
t354 = t381 * t379;
t355 = t381 * t380;
t345 = -t336 * t354 + t355 * t384;
t268 = t337 * t298 + t339 * t345;
t290 = t336 * t355 + t354 * t384;
t299 = -t336 * t379 + t384 * t380;
t269 = t290 * t339 + t337 * t299;
t357 = t381 * t384;
t293 = -t337 * t336 + t339 * t357;
t361 = t336 * t381;
t294 = t337 * t384 + t339 * t361;
t332 = sin(pkin(6));
t376 = t332 * t339;
t393 = Icges(3,5) * t294 + Icges(4,5) * t269 + Icges(3,6) * t293 + Icges(4,6) * t268 - t394 * t376;
t270 = t298 * t339 - t337 * t345;
t271 = -t337 * t290 + t299 * t339;
t295 = -t339 * t336 - t337 * t357;
t296 = -t337 * t361 + t339 * t384;
t377 = t332 * t337;
t392 = Icges(3,5) * t296 + Icges(4,5) * t271 + Icges(3,6) * t295 + Icges(4,6) * t270 + t394 * t377;
t331 = sin(pkin(12));
t333 = cos(pkin(12));
t241 = -t269 * t331 - t333 * t376;
t364 = t331 * t376;
t242 = t269 * t333 - t364;
t391 = -Icges(4,4) * t269 + Icges(5,5) * t242 + Icges(4,6) * t376 + Icges(5,6) * t241 + t395 * t268;
t243 = -t271 * t331 + t333 * t377;
t365 = t331 * t377;
t244 = t271 * t333 + t365;
t390 = -Icges(4,4) * t271 + Icges(5,5) * t244 - Icges(4,6) * t377 + Icges(5,6) * t243 + t395 * t270;
t288 = t299 * t332;
t289 = t298 * t332;
t389 = -Icges(4,5) * t289 + Icges(4,6) * t288 + (Icges(3,5) * t336 + Icges(3,6) * t384) * t332 + t394 * t381;
t275 = t289 * t331 + t333 * t381;
t359 = t381 * t331;
t276 = -t289 * t333 + t359;
t388 = Icges(4,4) * t289 + Icges(5,5) * t276 - Icges(4,6) * t381 + Icges(5,6) * t275 + t395 * t288;
t383 = pkin(2) * t384;
t382 = pkin(4) * t333;
t378 = Icges(2,4) * t337;
t230 = t269 * pkin(3) - t268 * qJ(4);
t362 = pkin(2) * t361 - qJ(3) * t332;
t266 = t337 * t383 + t339 * t362;
t374 = -t230 - t266;
t231 = pkin(3) * t271 - qJ(4) * t270;
t267 = -t337 * t362 + t339 * t383;
t373 = -t231 - t267;
t261 = -t289 * pkin(3) - t288 * qJ(4);
t300 = t332 * t336 * pkin(2) + qJ(3) * t381;
t372 = -t261 - t300;
t371 = qJD(2) * t332;
t370 = qJD(3) * t332;
t369 = pkin(12) + qJ(5);
t368 = V_base(5) * pkin(7) + V_base(1);
t363 = t381 * pkin(8);
t308 = t337 * t371 + V_base(4);
t328 = V_base(6) + qJD(1);
t358 = cos(t369);
t247 = -qJD(5) * t270 + t308;
t309 = qJD(2) * t381 + t328;
t353 = t332 * t358;
t274 = -qJD(5) * t288 + t309;
t307 = -t339 * t371 + V_base(5);
t301 = t337 * pkin(1) - pkin(8) * t376;
t352 = -t301 * t328 + V_base(5) * t363 + t368;
t302 = pkin(1) * t339 + pkin(8) * t377;
t351 = V_base(4) * t301 - t302 * V_base(5) + V_base(3);
t246 = -qJD(5) * t268 + t307;
t350 = t307 * t300 + t337 * t370 + t352;
t349 = qJD(3) * t381 + t308 * t266 + t351;
t348 = -qJD(4) * t270 + t307 * t261 + t350;
t347 = t328 * t302 + V_base(2) + (-t363 - pkin(7)) * V_base(4);
t346 = -qJD(4) * t288 + t308 * t230 + t349;
t344 = t309 * t267 - t339 * t370 + t347;
t177 = -pkin(4) * t364 - pkin(9) * t268 + t269 * t382;
t211 = pkin(4) * t359 - pkin(9) * t288 - t289 * t382;
t343 = t307 * t211 + (-t177 + t374) * t309 + t348;
t178 = pkin(4) * t365 - pkin(9) * t270 + t271 * t382;
t342 = t308 * t177 + (-t178 + t373) * t307 + t346;
t341 = -qJD(4) * t268 + t309 * t231 + t344;
t340 = t309 * t178 + (-t211 + t372) * t308 + t341;
t338 = cos(qJ(6));
t335 = sin(qJ(6));
t329 = Icges(2,4) * t339;
t327 = sin(t369);
t317 = rSges(2,1) * t339 - t337 * rSges(2,2);
t316 = t337 * rSges(2,1) + rSges(2,2) * t339;
t315 = Icges(2,1) * t339 - t378;
t314 = Icges(2,1) * t337 + t329;
t313 = -Icges(2,2) * t337 + t329;
t312 = Icges(2,2) * t339 + t378;
t305 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t304 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t303 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t287 = t381 * rSges(3,3) + (rSges(3,1) * t336 + rSges(3,2) * t384) * t332;
t286 = Icges(3,5) * t381 + (Icges(3,1) * t336 + Icges(3,4) * t384) * t332;
t285 = Icges(3,6) * t381 + (Icges(3,4) * t336 + Icges(3,2) * t384) * t332;
t279 = V_base(5) * rSges(2,3) - t316 * t328 + t368;
t278 = t317 * t328 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t277 = t316 * V_base(4) - t317 * V_base(5) + V_base(3);
t273 = -t289 * t358 + t327 * t381;
t272 = -t289 * t327 - t358 * t381;
t260 = rSges(3,1) * t296 + rSges(3,2) * t295 + rSges(3,3) * t377;
t259 = t294 * rSges(3,1) + t293 * rSges(3,2) - rSges(3,3) * t376;
t258 = Icges(3,1) * t296 + Icges(3,4) * t295 + Icges(3,5) * t377;
t257 = Icges(3,1) * t294 + Icges(3,4) * t293 - Icges(3,5) * t376;
t256 = Icges(3,4) * t296 + Icges(3,2) * t295 + Icges(3,6) * t377;
t255 = Icges(3,4) * t294 + Icges(3,2) * t293 - Icges(3,6) * t376;
t252 = -t289 * rSges(4,1) + t288 * rSges(4,2) + rSges(4,3) * t381;
t251 = -Icges(4,1) * t289 + Icges(4,4) * t288 + Icges(4,5) * t381;
t239 = t271 * t358 + t327 * t377;
t238 = t271 * t327 - t337 * t353;
t237 = t269 * t358 - t327 * t376;
t236 = t269 * t327 + t339 * t353;
t235 = t273 * t338 - t288 * t335;
t234 = -t273 * t335 - t288 * t338;
t233 = qJD(6) * t272 + t274;
t232 = pkin(5) * t273 + pkin(10) * t272;
t229 = rSges(5,1) * t276 + rSges(5,2) * t275 - rSges(5,3) * t288;
t228 = Icges(5,1) * t276 + Icges(5,4) * t275 - Icges(5,5) * t288;
t227 = Icges(5,4) * t276 + Icges(5,2) * t275 - Icges(5,6) * t288;
t225 = rSges(4,1) * t271 + rSges(4,2) * t270 + rSges(4,3) * t377;
t224 = t269 * rSges(4,1) + t268 * rSges(4,2) - rSges(4,3) * t376;
t223 = Icges(4,1) * t271 + Icges(4,4) * t270 + Icges(4,5) * t377;
t222 = Icges(4,1) * t269 + Icges(4,4) * t268 - Icges(4,5) * t376;
t217 = rSges(6,1) * t273 - rSges(6,2) * t272 - rSges(6,3) * t288;
t216 = Icges(6,1) * t273 - Icges(6,4) * t272 - Icges(6,5) * t288;
t215 = Icges(6,4) * t273 - Icges(6,2) * t272 - Icges(6,6) * t288;
t214 = Icges(6,5) * t273 - Icges(6,6) * t272 - Icges(6,3) * t288;
t210 = t239 * t338 - t270 * t335;
t209 = -t239 * t335 - t270 * t338;
t208 = t237 * t338 - t268 * t335;
t207 = -t237 * t335 - t268 * t338;
t205 = qJD(6) * t238 + t247;
t204 = qJD(6) * t236 + t246;
t203 = -t259 * t309 + t287 * t307 + t352;
t202 = t309 * t260 - t308 * t287 + t347;
t201 = pkin(5) * t239 + pkin(10) * t238;
t200 = pkin(5) * t237 + pkin(10) * t236;
t199 = t259 * t308 - t260 * t307 + t351;
t198 = rSges(5,1) * t244 + rSges(5,2) * t243 - rSges(5,3) * t270;
t197 = rSges(5,1) * t242 + rSges(5,2) * t241 - rSges(5,3) * t268;
t196 = Icges(5,1) * t244 + Icges(5,4) * t243 - Icges(5,5) * t270;
t195 = Icges(5,1) * t242 + Icges(5,4) * t241 - Icges(5,5) * t268;
t194 = Icges(5,4) * t244 + Icges(5,2) * t243 - Icges(5,6) * t270;
t193 = Icges(5,4) * t242 + Icges(5,2) * t241 - Icges(5,6) * t268;
t190 = rSges(6,1) * t239 - rSges(6,2) * t238 - rSges(6,3) * t270;
t189 = rSges(6,1) * t237 - rSges(6,2) * t236 - rSges(6,3) * t268;
t188 = Icges(6,1) * t239 - Icges(6,4) * t238 - Icges(6,5) * t270;
t187 = Icges(6,1) * t237 - Icges(6,4) * t236 - Icges(6,5) * t268;
t186 = Icges(6,4) * t239 - Icges(6,2) * t238 - Icges(6,6) * t270;
t185 = Icges(6,4) * t237 - Icges(6,2) * t236 - Icges(6,6) * t268;
t184 = Icges(6,5) * t239 - Icges(6,6) * t238 - Icges(6,3) * t270;
t183 = Icges(6,5) * t237 - Icges(6,6) * t236 - Icges(6,3) * t268;
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
t163 = t229 * t307 + (-t197 + t374) * t309 + t348;
t162 = t309 * t198 + (-t229 + t372) * t308 + t341;
t161 = t197 * t308 + (-t198 + t373) * t307 + t346;
t160 = -t189 * t274 + t217 * t246 + t343;
t159 = t274 * t190 - t247 * t217 + t340;
t158 = t189 * t247 - t190 * t246 + t342;
t157 = -t173 * t233 + t182 * t204 - t200 * t274 + t232 * t246 + t343;
t156 = t233 * t174 - t205 * t182 + t274 * t201 - t247 * t232 + t340;
t155 = t173 * t205 - t174 * t204 + t200 * t247 - t201 * t246 + t342;
t1 = m(7) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(6) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(5) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(4) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(3) * (t199 ^ 2 + t202 ^ 2 + t203 ^ 2) / 0.2e1 + t204 * ((t168 * t236 + t170 * t207 + t172 * t208) * t205 + (t236 * t167 + t207 * t169 + t208 * t171) * t204 + (t179 * t236 + t180 * t207 + t181 * t208) * t233) / 0.2e1 + t205 * ((t238 * t168 + t209 * t170 + t210 * t172) * t205 + (t167 * t238 + t169 * t209 + t171 * t210) * t204 + (t179 * t238 + t180 * t209 + t181 * t210) * t233) / 0.2e1 + t233 * ((t168 * t272 + t170 * t234 + t172 * t235) * t205 + (t167 * t272 + t169 * t234 + t171 * t235) * t204 + (t272 * t179 + t234 * t180 + t235 * t181) * t233) / 0.2e1 + t246 * ((-t184 * t268 - t186 * t236 + t188 * t237) * t247 + (-t268 * t183 - t236 * t185 + t237 * t187) * t246 + (-t214 * t268 - t215 * t236 + t216 * t237) * t274) / 0.2e1 + t247 * ((-t270 * t184 - t238 * t186 + t239 * t188) * t247 + (-t183 * t270 - t185 * t238 + t187 * t239) * t246 + (-t214 * t270 - t215 * t238 + t216 * t239) * t274) / 0.2e1 + m(2) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + t274 * ((-t184 * t288 - t186 * t272 + t188 * t273) * t247 + (-t183 * t288 - t185 * t272 + t187 * t273) * t246 + (-t288 * t214 - t272 * t215 + t273 * t216) * t274) / 0.2e1 + m(1) * (t303 ^ 2 + t304 ^ 2 + t305 ^ 2) / 0.2e1 + ((-t337 * t312 + t314 * t339 + Icges(1,4)) * V_base(5) + (-t337 * t313 + t339 * t315 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t339 * t312 + t337 * t314 + Icges(1,2)) * V_base(5) + (t313 * t339 + t337 * t315 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t227 * t241 + t228 * t242 + t269 * t251 - t268 * t388 + t293 * t285 + t294 * t286 - t376 * t389) * t309 + (t194 * t241 + t196 * t242 + t269 * t223 + t293 * t256 + t294 * t258 - t268 * t390 - t376 * t392) * t308 + (t241 * t193 + t242 * t195 + t269 * t222 + t293 * t255 + t294 * t257 - t391 * t268 - t393 * t376) * t307) * t307 / 0.2e1 + ((t227 * t243 + t228 * t244 + t251 * t271 - t270 * t388 + t285 * t295 + t286 * t296 + t377 * t389) * t309 + (t243 * t194 + t244 * t196 + t271 * t223 + t295 * t256 + t296 * t258 - t390 * t270 + t392 * t377) * t308 + (t193 * t243 + t195 * t244 + t222 * t271 + t255 * t295 + t257 * t296 - t391 * t270 + t377 * t393) * t307) * t308 / 0.2e1 + ((-t289 * t251 + (t285 * t384 + t286 * t336) * t332 + t275 * t227 + t276 * t228 + t389 * t381 - t388 * t288) * t309 + (-t289 * t223 + (t256 * t384 + t258 * t336) * t332 + t194 * t275 + t196 * t276 + t392 * t381 - t390 * t288) * t308 + (-t289 * t222 + (t255 * t384 + t257 * t336) * t332 + t193 * t275 + t195 * t276 + t393 * t381 - t391 * t288) * t307) * t309 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t337 + Icges(2,6) * t339) * V_base(5) + (Icges(2,5) * t339 - Icges(2,6) * t337) * V_base(4) + Icges(2,3) * t328 / 0.2e1) * t328;
T  = t1;
