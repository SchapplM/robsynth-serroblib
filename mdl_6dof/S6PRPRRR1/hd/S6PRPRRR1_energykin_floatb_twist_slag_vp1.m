% Calculate kinetic energy for
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:15
% EndTime: 2019-03-08 20:22:21
% DurationCPUTime: 6.22s
% Computational Cost: add. (3937->441), mult. (8194->644), div. (0->0), fcn. (10180->14), ass. (0->198)
t334 = sin(qJ(2));
t377 = sin(pkin(12));
t378 = cos(pkin(12));
t383 = cos(qJ(2));
t297 = -t334 * t378 - t377 * t383;
t329 = sin(pkin(11));
t331 = cos(pkin(11));
t379 = cos(pkin(6));
t351 = t379 * t377;
t352 = t379 * t378;
t345 = -t334 * t351 + t352 * t383;
t267 = t329 * t297 + t331 * t345;
t289 = t334 * t352 + t351 * t383;
t298 = -t334 * t377 + t383 * t378;
t268 = t289 * t331 + t298 * t329;
t330 = sin(pkin(6));
t374 = t330 * t331;
t212 = Icges(4,5) * t268 + Icges(4,6) * t267 - Icges(4,3) * t374;
t355 = t379 * t383;
t291 = -t329 * t334 + t331 * t355;
t358 = t334 * t379;
t292 = t329 * t383 + t331 * t358;
t252 = Icges(3,5) * t292 + Icges(3,6) * t291 - Icges(3,3) * t374;
t391 = t212 + t252;
t269 = t297 * t331 - t329 * t345;
t270 = -t289 * t329 + t298 * t331;
t375 = t329 * t330;
t213 = Icges(4,5) * t270 + Icges(4,6) * t269 + Icges(4,3) * t375;
t293 = -t329 * t355 - t331 * t334;
t294 = -t329 * t358 + t331 * t383;
t253 = Icges(3,5) * t294 + Icges(3,6) * t293 + Icges(3,3) * t375;
t390 = t213 + t253;
t287 = t298 * t330;
t288 = t297 * t330;
t248 = -Icges(4,5) * t288 + Icges(4,6) * t287 + Icges(4,3) * t379;
t283 = Icges(3,3) * t379 + (Icges(3,5) * t334 + Icges(3,6) * t383) * t330;
t389 = t248 + t283;
t382 = pkin(2) * t383;
t336 = cos(qJ(4));
t381 = pkin(4) * t336;
t376 = Icges(2,4) * t329;
t333 = sin(qJ(4));
t373 = t330 * t333;
t372 = t330 * t336;
t371 = qJ(4) + qJ(5);
t370 = qJD(2) * t330;
t369 = qJD(3) * t330;
t368 = V_base(5) * qJ(1) + V_base(1);
t364 = qJD(1) + V_base(3);
t363 = t329 * t373;
t362 = t331 * t373;
t361 = t379 * pkin(7);
t306 = t329 * t370 + V_base(4);
t317 = qJD(2) * t379 + V_base(6);
t360 = cos(t371);
t359 = pkin(2) * t358 - qJ(3) * t330;
t356 = t379 * t333;
t241 = -qJD(4) * t269 + t306;
t274 = -qJD(4) * t287 + t317;
t353 = t330 * t360;
t210 = -qJD(5) * t269 + t241;
t247 = -qJD(5) * t287 + t274;
t305 = -t331 * t370 + V_base(5);
t240 = -qJD(4) * t267 + t305;
t300 = pkin(1) * t329 - pkin(7) * t374;
t350 = -t300 * V_base(6) + V_base(5) * t361 + t368;
t301 = pkin(1) * t331 + pkin(7) * t375;
t349 = V_base(4) * t300 - t301 * V_base(5) + t364;
t209 = -qJD(5) * t267 + t240;
t299 = t330 * t334 * pkin(2) + qJ(3) * t379;
t348 = t305 * t299 + t329 * t369 + t350;
t260 = t329 * t382 + t331 * t359;
t347 = qJD(3) * t379 + t306 * t260 + t349;
t346 = V_base(6) * t301 + V_base(2) + (-t361 - qJ(1)) * V_base(4);
t225 = t268 * pkin(3) - t267 * pkin(8);
t266 = -t288 * pkin(3) - t287 * pkin(8);
t344 = t305 * t266 + (-t225 - t260) * t317 + t348;
t226 = t270 * pkin(3) - t269 * pkin(8);
t261 = -t329 * t359 + t331 * t382;
t343 = t306 * t225 + (-t226 - t261) * t305 + t347;
t342 = t317 * t261 - t331 * t369 + t346;
t173 = -pkin(4) * t362 - pkin(9) * t267 + t268 * t381;
t207 = pkin(4) * t356 - pkin(9) * t287 - t288 * t381;
t341 = -t173 * t274 + t240 * t207 + t344;
t174 = pkin(4) * t363 - pkin(9) * t269 + t270 * t381;
t340 = t241 * t173 - t174 * t240 + t343;
t339 = t317 * t226 + (-t266 - t299) * t306 + t342;
t338 = t274 * t174 - t241 * t207 + t339;
t335 = cos(qJ(6));
t332 = sin(qJ(6));
t327 = sin(t371);
t326 = Icges(2,4) * t331;
t315 = rSges(2,1) * t331 - rSges(2,2) * t329;
t314 = rSges(2,1) * t329 + rSges(2,2) * t331;
t313 = Icges(2,1) * t331 - t376;
t312 = Icges(2,1) * t329 + t326;
t311 = -Icges(2,2) * t329 + t326;
t310 = Icges(2,2) * t331 + t376;
t304 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t303 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t302 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t286 = t379 * rSges(3,3) + (rSges(3,1) * t334 + rSges(3,2) * t383) * t330;
t285 = Icges(3,5) * t379 + (Icges(3,1) * t334 + Icges(3,4) * t383) * t330;
t284 = Icges(3,6) * t379 + (Icges(3,4) * t334 + Icges(3,2) * t383) * t330;
t278 = V_base(5) * rSges(2,3) - t314 * V_base(6) + t368;
t277 = t315 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t276 = -t288 * t336 + t356;
t275 = t288 * t333 + t336 * t379;
t273 = t314 * V_base(4) - t315 * V_base(5) + t364;
t272 = -t288 * t360 + t327 * t379;
t271 = -t288 * t327 - t360 * t379;
t259 = rSges(3,1) * t294 + rSges(3,2) * t293 + rSges(3,3) * t375;
t258 = rSges(3,1) * t292 + rSges(3,2) * t291 - rSges(3,3) * t374;
t257 = Icges(3,1) * t294 + Icges(3,4) * t293 + Icges(3,5) * t375;
t256 = Icges(3,1) * t292 + Icges(3,4) * t291 - Icges(3,5) * t374;
t255 = Icges(3,4) * t294 + Icges(3,2) * t293 + Icges(3,6) * t375;
t254 = Icges(3,4) * t292 + Icges(3,2) * t291 - Icges(3,6) * t374;
t251 = -t288 * rSges(4,1) + t287 * rSges(4,2) + rSges(4,3) * t379;
t250 = -Icges(4,1) * t288 + Icges(4,4) * t287 + Icges(4,5) * t379;
t249 = -Icges(4,4) * t288 + Icges(4,2) * t287 + Icges(4,6) * t379;
t245 = t270 * t336 + t363;
t244 = -t270 * t333 + t329 * t372;
t243 = t268 * t336 - t362;
t242 = -t268 * t333 - t331 * t372;
t237 = t270 * t360 + t327 * t375;
t236 = t270 * t327 - t329 * t353;
t235 = t268 * t360 - t327 * t374;
t234 = t268 * t327 + t331 * t353;
t233 = t272 * t335 - t287 * t332;
t232 = -t272 * t332 - t287 * t335;
t231 = pkin(5) * t272 + pkin(10) * t271;
t230 = rSges(5,1) * t276 + rSges(5,2) * t275 - rSges(5,3) * t287;
t229 = Icges(5,1) * t276 + Icges(5,4) * t275 - Icges(5,5) * t287;
t228 = Icges(5,4) * t276 + Icges(5,2) * t275 - Icges(5,6) * t287;
t227 = Icges(5,5) * t276 + Icges(5,6) * t275 - Icges(5,3) * t287;
t224 = qJD(6) * t271 + t247;
t223 = rSges(6,1) * t272 - rSges(6,2) * t271 - rSges(6,3) * t287;
t222 = Icges(6,1) * t272 - Icges(6,4) * t271 - Icges(6,5) * t287;
t221 = Icges(6,4) * t272 - Icges(6,2) * t271 - Icges(6,6) * t287;
t220 = Icges(6,5) * t272 - Icges(6,6) * t271 - Icges(6,3) * t287;
t219 = rSges(4,1) * t270 + rSges(4,2) * t269 + rSges(4,3) * t375;
t218 = rSges(4,1) * t268 + rSges(4,2) * t267 - rSges(4,3) * t374;
t217 = Icges(4,1) * t270 + Icges(4,4) * t269 + Icges(4,5) * t375;
t216 = Icges(4,1) * t268 + Icges(4,4) * t267 - Icges(4,5) * t374;
t215 = Icges(4,4) * t270 + Icges(4,2) * t269 + Icges(4,6) * t375;
t214 = Icges(4,4) * t268 + Icges(4,2) * t267 - Icges(4,6) * t374;
t206 = t237 * t335 - t269 * t332;
t205 = -t237 * t332 - t269 * t335;
t204 = t235 * t335 - t267 * t332;
t203 = -t235 * t332 - t267 * t335;
t202 = -t258 * t317 + t286 * t305 + t350;
t201 = t317 * t259 - t306 * t286 + t346;
t200 = pkin(5) * t237 + pkin(10) * t236;
t199 = pkin(5) * t235 + pkin(10) * t234;
t198 = t258 * t306 - t259 * t305 + t349;
t196 = rSges(5,1) * t245 + rSges(5,2) * t244 - rSges(5,3) * t269;
t195 = rSges(5,1) * t243 + rSges(5,2) * t242 - rSges(5,3) * t267;
t194 = Icges(5,1) * t245 + Icges(5,4) * t244 - Icges(5,5) * t269;
t193 = Icges(5,1) * t243 + Icges(5,4) * t242 - Icges(5,5) * t267;
t192 = Icges(5,4) * t245 + Icges(5,2) * t244 - Icges(5,6) * t269;
t191 = Icges(5,4) * t243 + Icges(5,2) * t242 - Icges(5,6) * t267;
t190 = Icges(5,5) * t245 + Icges(5,6) * t244 - Icges(5,3) * t269;
t189 = Icges(5,5) * t243 + Icges(5,6) * t242 - Icges(5,3) * t267;
t188 = qJD(6) * t236 + t210;
t187 = qJD(6) * t234 + t209;
t186 = rSges(7,1) * t233 + rSges(7,2) * t232 + rSges(7,3) * t271;
t185 = Icges(7,1) * t233 + Icges(7,4) * t232 + Icges(7,5) * t271;
t184 = Icges(7,4) * t233 + Icges(7,2) * t232 + Icges(7,6) * t271;
t183 = Icges(7,5) * t233 + Icges(7,6) * t232 + Icges(7,3) * t271;
t182 = rSges(6,1) * t237 - rSges(6,2) * t236 - rSges(6,3) * t269;
t181 = rSges(6,1) * t235 - rSges(6,2) * t234 - rSges(6,3) * t267;
t180 = Icges(6,1) * t237 - Icges(6,4) * t236 - Icges(6,5) * t269;
t179 = Icges(6,1) * t235 - Icges(6,4) * t234 - Icges(6,5) * t267;
t178 = Icges(6,4) * t237 - Icges(6,2) * t236 - Icges(6,6) * t269;
t177 = Icges(6,4) * t235 - Icges(6,2) * t234 - Icges(6,6) * t267;
t176 = Icges(6,5) * t237 - Icges(6,6) * t236 - Icges(6,3) * t269;
t175 = Icges(6,5) * t235 - Icges(6,6) * t234 - Icges(6,3) * t267;
t170 = rSges(7,1) * t206 + rSges(7,2) * t205 + rSges(7,3) * t236;
t169 = rSges(7,1) * t204 + rSges(7,2) * t203 + rSges(7,3) * t234;
t168 = Icges(7,1) * t206 + Icges(7,4) * t205 + Icges(7,5) * t236;
t167 = Icges(7,1) * t204 + Icges(7,4) * t203 + Icges(7,5) * t234;
t166 = Icges(7,4) * t206 + Icges(7,2) * t205 + Icges(7,6) * t236;
t165 = Icges(7,4) * t204 + Icges(7,2) * t203 + Icges(7,6) * t234;
t164 = Icges(7,5) * t206 + Icges(7,6) * t205 + Icges(7,3) * t236;
t163 = Icges(7,5) * t204 + Icges(7,6) * t203 + Icges(7,3) * t234;
t162 = t251 * t305 + (-t218 - t260) * t317 + t348;
t161 = t317 * t219 + (-t251 - t299) * t306 + t342;
t160 = t218 * t306 + (-t219 - t261) * t305 + t347;
t159 = -t195 * t274 + t230 * t240 + t344;
t158 = t274 * t196 - t241 * t230 + t339;
t157 = t195 * t241 - t196 * t240 + t343;
t156 = -t181 * t247 + t209 * t223 + t341;
t155 = t247 * t182 - t210 * t223 + t338;
t154 = t181 * t210 - t182 * t209 + t340;
t153 = -t169 * t224 + t186 * t187 - t199 * t247 + t209 * t231 + t341;
t152 = t224 * t170 - t188 * t186 + t247 * t200 - t210 * t231 + t338;
t151 = t169 * t188 - t170 * t187 + t199 * t210 - t200 * t209 + t340;
t1 = m(7) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(6) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(3) * (t198 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + t187 * ((t164 * t234 + t166 * t203 + t168 * t204) * t188 + (t163 * t234 + t165 * t203 + t167 * t204) * t187 + (t183 * t234 + t184 * t203 + t185 * t204) * t224) / 0.2e1 + t188 * ((t164 * t236 + t166 * t205 + t168 * t206) * t188 + (t163 * t236 + t165 * t205 + t167 * t206) * t187 + (t183 * t236 + t184 * t205 + t185 * t206) * t224) / 0.2e1 + t209 * ((-t176 * t267 - t178 * t234 + t180 * t235) * t210 + (-t175 * t267 - t177 * t234 + t179 * t235) * t209 + (-t220 * t267 - t221 * t234 + t222 * t235) * t247) / 0.2e1 + t210 * ((-t176 * t269 - t178 * t236 + t237 * t180) * t210 + (-t175 * t269 - t177 * t236 + t179 * t237) * t209 + (-t220 * t269 - t221 * t236 + t222 * t237) * t247) / 0.2e1 + t224 * ((t164 * t271 + t166 * t232 + t168 * t233) * t188 + (t163 * t271 + t165 * t232 + t167 * t233) * t187 + (t183 * t271 + t184 * t232 + t185 * t233) * t224) / 0.2e1 + t240 * ((-t190 * t267 + t192 * t242 + t194 * t243) * t241 + (-t189 * t267 + t191 * t242 + t193 * t243) * t240 + (-t227 * t267 + t228 * t242 + t229 * t243) * t274) / 0.2e1 + t241 * ((-t190 * t269 + t192 * t244 + t194 * t245) * t241 + (-t189 * t269 + t191 * t244 + t193 * t245) * t240 + (-t227 * t269 + t228 * t244 + t229 * t245) * t274) / 0.2e1 + m(2) * (t273 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + t247 * ((-t176 * t287 - t178 * t271 + t180 * t272) * t210 + (-t175 * t287 - t177 * t271 + t179 * t272) * t209 + (-t220 * t287 - t221 * t271 + t222 * t272) * t247) / 0.2e1 + t274 * ((-t190 * t287 + t192 * t275 + t194 * t276) * t241 + (-t189 * t287 + t191 * t275 + t193 * t276) * t240 + (-t227 * t287 + t228 * t275 + t229 * t276) * t274) / 0.2e1 + m(1) * (t302 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + ((t249 * t267 + t250 * t268 + t284 * t291 + t285 * t292 - t389 * t374) * t317 + (t215 * t267 + t217 * t268 + t255 * t291 + t257 * t292 - t390 * t374) * t306 + (t214 * t267 + t216 * t268 + t254 * t291 + t256 * t292 - t391 * t374) * t305) * t305 / 0.2e1 + ((t249 * t269 + t250 * t270 + t284 * t293 + t285 * t294 + t389 * t375) * t317 + (t215 * t269 + t217 * t270 + t255 * t293 + t257 * t294 + t390 * t375) * t306 + (t214 * t269 + t216 * t270 + t254 * t293 + t256 * t294 + t391 * t375) * t305) * t306 / 0.2e1 + ((t213 * t379 + t287 * t215 - t288 * t217) * t306 + (t212 * t379 + t287 * t214 - t288 * t216) * t305 + (t248 * t379 + t287 * t249 - t288 * t250) * t317 + ((t255 * t383 + t257 * t334) * t306 + (t254 * t383 + t256 * t334) * t305 + (t284 * t383 + t285 * t334) * t317) * t330 + (t252 * t305 + t253 * t306 + t283 * t317) * t379) * t317 / 0.2e1 + ((-t310 * t329 + t312 * t331 + Icges(1,4)) * V_base(5) + (-t311 * t329 + t313 * t331 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t310 * t331 + t312 * t329 + Icges(1,2)) * V_base(5) + (t311 * t331 + t313 * t329 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t331 - Icges(2,6) * t329 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t329 + Icges(2,6) * t331 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
