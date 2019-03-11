% Calculate kinetic energy for
% S6RRRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:16:12
% EndTime: 2019-03-10 02:16:16
% DurationCPUTime: 4.48s
% Computational Cost: add. (3380->386), mult. (6925->562), div. (0->0), fcn. (8407->12), ass. (0->175)
t387 = Icges(6,1) + Icges(7,1);
t386 = -Icges(6,4) + Icges(7,5);
t385 = Icges(7,4) + Icges(6,5);
t384 = Icges(6,2) + Icges(7,3);
t383 = Icges(7,2) + Icges(6,3);
t382 = -Icges(6,6) + Icges(7,6);
t381 = rSges(7,1) + pkin(5);
t380 = rSges(7,3) + qJ(6);
t322 = cos(pkin(6));
t326 = sin(qJ(1));
t328 = cos(qJ(2));
t355 = t326 * t328;
t325 = sin(qJ(2));
t329 = cos(qJ(1));
t356 = t325 * t329;
t289 = t322 * t356 + t355;
t324 = sin(qJ(3));
t321 = sin(pkin(6));
t358 = t321 * t329;
t367 = cos(qJ(3));
t268 = t289 * t367 - t324 * t358;
t354 = t328 * t329;
t357 = t325 * t326;
t288 = -t322 * t354 + t357;
t353 = qJ(4) + qJ(5);
t318 = sin(t353);
t343 = cos(t353);
t232 = t268 * t318 - t288 * t343;
t233 = t268 * t343 + t288 * t318;
t344 = t321 * t367;
t267 = t289 * t324 + t329 * t344;
t379 = t384 * t232 + t386 * t233 + t382 * t267;
t291 = -t322 * t357 + t354;
t360 = t321 * t326;
t270 = t291 * t367 + t324 * t360;
t290 = t322 * t355 + t356;
t234 = t270 * t318 - t290 * t343;
t235 = t270 * t343 + t290 * t318;
t269 = t291 * t324 - t326 * t344;
t378 = t384 * t234 + t386 * t235 + t382 * t269;
t377 = t382 * t232 + t385 * t233 + t383 * t267;
t376 = t382 * t234 + t385 * t235 + t383 * t269;
t375 = t386 * t232 + t387 * t233 + t385 * t267;
t374 = t386 * t234 + t387 * t235 + t385 * t269;
t287 = t322 * t324 + t325 * t344;
t359 = t321 * t328;
t259 = t287 * t318 + t343 * t359;
t260 = t287 * t343 - t318 * t359;
t286 = t321 * t324 * t325 - t322 * t367;
t373 = t384 * t259 + t386 * t260 + t382 * t286;
t372 = t382 * t259 + t385 * t260 + t383 * t286;
t371 = t386 * t259 + t387 * t260 + t385 * t286;
t366 = pkin(8) * t322;
t327 = cos(qJ(4));
t365 = pkin(4) * t327;
t363 = Icges(2,4) * t326;
t323 = sin(qJ(4));
t362 = t288 * t323;
t361 = t290 * t323;
t352 = rSges(7,2) * t267 + t380 * t232 + t233 * t381;
t351 = rSges(7,2) * t269 + t380 * t234 + t235 * t381;
t350 = rSges(7,2) * t286 + t380 * t259 + t260 * t381;
t349 = qJD(2) * t321;
t348 = V_base(5) * pkin(7) + V_base(1);
t345 = t323 * t359;
t300 = t326 * t349 + V_base(4);
t317 = V_base(6) + qJD(1);
t266 = qJD(3) * t290 + t300;
t301 = qJD(2) * t322 + t317;
t231 = qJD(4) * t269 + t266;
t299 = -t329 * t349 + V_base(5);
t294 = pkin(1) * t326 - pkin(8) * t358;
t342 = -t294 * t317 + V_base(5) * t366 + t348;
t295 = pkin(1) * t329 + pkin(8) * t360;
t341 = V_base(4) * t294 - t295 * V_base(5) + V_base(3);
t265 = qJD(3) * t288 + t299;
t230 = qJD(4) * t267 + t265;
t284 = -qJD(3) * t359 + t301;
t340 = t317 * t295 + V_base(2) + (-pkin(7) - t366) * V_base(4);
t256 = qJD(4) * t286 + t284;
t257 = pkin(2) * t289 + pkin(9) * t288;
t293 = (pkin(2) * t325 - pkin(9) * t328) * t321;
t339 = -t257 * t301 + t299 * t293 + t342;
t258 = pkin(2) * t291 + pkin(9) * t290;
t338 = t300 * t257 - t258 * t299 + t341;
t337 = t301 * t258 - t293 * t300 + t340;
t228 = t268 * pkin(3) + t267 * pkin(10);
t255 = t287 * pkin(3) + t286 * pkin(10);
t336 = -t228 * t284 + t265 * t255 + t339;
t229 = t270 * pkin(3) + t269 * pkin(10);
t335 = t266 * t228 - t229 * t265 + t338;
t334 = t284 * t229 - t255 * t266 + t337;
t169 = pkin(4) * t362 + pkin(11) * t267 + t268 * t365;
t212 = -pkin(4) * t345 + pkin(11) * t286 + t287 * t365;
t333 = -t169 * t256 + t230 * t212 + t336;
t170 = pkin(4) * t361 + pkin(11) * t269 + t270 * t365;
t332 = t231 * t169 - t170 * t230 + t335;
t331 = t256 * t170 - t212 * t231 + t334;
t319 = Icges(2,4) * t329;
t309 = rSges(2,1) * t329 - rSges(2,2) * t326;
t308 = rSges(2,1) * t326 + rSges(2,2) * t329;
t307 = Icges(2,1) * t329 - t363;
t306 = Icges(2,1) * t326 + t319;
t305 = -Icges(2,2) * t326 + t319;
t304 = Icges(2,2) * t329 + t363;
t298 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t297 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t296 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t280 = rSges(3,3) * t322 + (rSges(3,1) * t325 + rSges(3,2) * t328) * t321;
t279 = Icges(3,5) * t322 + (Icges(3,1) * t325 + Icges(3,4) * t328) * t321;
t278 = Icges(3,6) * t322 + (Icges(3,4) * t325 + Icges(3,2) * t328) * t321;
t277 = Icges(3,3) * t322 + (Icges(3,5) * t325 + Icges(3,6) * t328) * t321;
t274 = V_base(5) * rSges(2,3) - t308 * t317 + t348;
t273 = t309 * t317 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t271 = t308 * V_base(4) - t309 * V_base(5) + V_base(3);
t264 = t287 * t327 - t345;
t263 = -t287 * t323 - t327 * t359;
t254 = rSges(3,1) * t291 - rSges(3,2) * t290 + rSges(3,3) * t360;
t253 = rSges(3,1) * t289 - rSges(3,2) * t288 - rSges(3,3) * t358;
t252 = Icges(3,1) * t291 - Icges(3,4) * t290 + Icges(3,5) * t360;
t251 = Icges(3,1) * t289 - Icges(3,4) * t288 - Icges(3,5) * t358;
t250 = Icges(3,4) * t291 - Icges(3,2) * t290 + Icges(3,6) * t360;
t249 = Icges(3,4) * t289 - Icges(3,2) * t288 - Icges(3,6) * t358;
t248 = Icges(3,5) * t291 - Icges(3,6) * t290 + Icges(3,3) * t360;
t247 = Icges(3,5) * t289 - Icges(3,6) * t288 - Icges(3,3) * t358;
t246 = rSges(4,1) * t287 - rSges(4,2) * t286 - rSges(4,3) * t359;
t245 = Icges(4,1) * t287 - Icges(4,4) * t286 - Icges(4,5) * t359;
t244 = Icges(4,4) * t287 - Icges(4,2) * t286 - Icges(4,6) * t359;
t243 = Icges(4,5) * t287 - Icges(4,6) * t286 - Icges(4,3) * t359;
t240 = t270 * t327 + t361;
t239 = -t270 * t323 + t290 * t327;
t238 = t268 * t327 + t362;
t237 = -t268 * t323 + t288 * t327;
t236 = qJD(5) * t286 + t256;
t225 = rSges(4,1) * t270 - rSges(4,2) * t269 + rSges(4,3) * t290;
t224 = rSges(4,1) * t268 - rSges(4,2) * t267 + rSges(4,3) * t288;
t222 = Icges(4,1) * t270 - Icges(4,4) * t269 + Icges(4,5) * t290;
t221 = Icges(4,1) * t268 - Icges(4,4) * t267 + Icges(4,5) * t288;
t220 = Icges(4,4) * t270 - Icges(4,2) * t269 + Icges(4,6) * t290;
t219 = Icges(4,4) * t268 - Icges(4,2) * t267 + Icges(4,6) * t288;
t218 = Icges(4,5) * t270 - Icges(4,6) * t269 + Icges(4,3) * t290;
t217 = Icges(4,5) * t268 - Icges(4,6) * t267 + Icges(4,3) * t288;
t216 = rSges(5,1) * t264 + rSges(5,2) * t263 + rSges(5,3) * t286;
t215 = Icges(5,1) * t264 + Icges(5,4) * t263 + Icges(5,5) * t286;
t214 = Icges(5,4) * t264 + Icges(5,2) * t263 + Icges(5,6) * t286;
t213 = Icges(5,5) * t264 + Icges(5,6) * t263 + Icges(5,3) * t286;
t211 = rSges(6,1) * t260 - rSges(6,2) * t259 + rSges(6,3) * t286;
t203 = qJD(5) * t269 + t231;
t202 = qJD(5) * t267 + t230;
t198 = -t253 * t301 + t280 * t299 + t342;
t197 = t254 * t301 - t280 * t300 + t340;
t196 = rSges(5,1) * t240 + rSges(5,2) * t239 + rSges(5,3) * t269;
t195 = rSges(5,1) * t238 + rSges(5,2) * t237 + rSges(5,3) * t267;
t194 = Icges(5,1) * t240 + Icges(5,4) * t239 + Icges(5,5) * t269;
t193 = Icges(5,1) * t238 + Icges(5,4) * t237 + Icges(5,5) * t267;
t192 = Icges(5,4) * t240 + Icges(5,2) * t239 + Icges(5,6) * t269;
t191 = Icges(5,4) * t238 + Icges(5,2) * t237 + Icges(5,6) * t267;
t190 = Icges(5,5) * t240 + Icges(5,6) * t239 + Icges(5,3) * t269;
t189 = Icges(5,5) * t238 + Icges(5,6) * t237 + Icges(5,3) * t267;
t188 = t253 * t300 - t254 * t299 + t341;
t187 = rSges(6,1) * t235 - rSges(6,2) * t234 + rSges(6,3) * t269;
t185 = rSges(6,1) * t233 - rSges(6,2) * t232 + rSges(6,3) * t267;
t166 = -t224 * t284 + t246 * t265 + t339;
t165 = t225 * t284 - t246 * t266 + t337;
t164 = t224 * t266 - t225 * t265 + t338;
t163 = -t195 * t256 + t216 * t230 + t336;
t162 = t196 * t256 - t216 * t231 + t334;
t161 = t195 * t231 - t196 * t230 + t335;
t160 = -t185 * t236 + t202 * t211 + t333;
t159 = t187 * t236 - t203 * t211 + t331;
t158 = t185 * t203 - t187 * t202 + t332;
t157 = qJD(6) * t234 + t202 * t350 - t236 * t352 + t333;
t156 = qJD(6) * t232 - t203 * t350 + t236 * t351 + t331;
t155 = qJD(6) * t259 - t202 * t351 + t203 * t352 + t332;
t1 = t301 * ((t247 * t299 + t248 * t300 + t277 * t301) * t322 + ((t250 * t328 + t252 * t325) * t300 + (t249 * t328 + t251 * t325) * t299 + (t278 * t328 + t279 * t325) * t301) * t321) / 0.2e1 + m(1) * (t296 ^ 2 + t297 ^ 2 + t298 ^ 2) / 0.2e1 + t265 * ((t218 * t288 - t220 * t267 + t222 * t268) * t266 + (t217 * t288 - t219 * t267 + t221 * t268) * t265 + (t243 * t288 - t244 * t267 + t245 * t268) * t284) / 0.2e1 + t266 * ((t218 * t290 - t220 * t269 + t222 * t270) * t266 + (t217 * t290 - t219 * t269 + t221 * t270) * t265 + (t243 * t290 - t244 * t269 + t245 * t270) * t284) / 0.2e1 + t256 * ((t190 * t286 + t192 * t263 + t194 * t264) * t231 + (t189 * t286 + t191 * t263 + t193 * t264) * t230 + (t213 * t286 + t214 * t263 + t215 * t264) * t256) / 0.2e1 + m(2) * (t271 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + m(3) * (t188 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + t230 * ((t190 * t267 + t192 * t237 + t194 * t238) * t231 + (t267 * t189 + t237 * t191 + t238 * t193) * t230 + (t213 * t267 + t214 * t237 + t215 * t238) * t256) / 0.2e1 + t231 * ((t269 * t190 + t239 * t192 + t240 * t194) * t231 + (t189 * t269 + t191 * t239 + t193 * t240) * t230 + (t213 * t269 + t214 * t239 + t215 * t240) * t256) / 0.2e1 + t299 * ((-t248 * t358 - t250 * t288 + t252 * t289) * t300 + (-t247 * t358 - t249 * t288 + t251 * t289) * t299 + (-t277 * t358 - t278 * t288 + t279 * t289) * t301) / 0.2e1 + t284 * ((-t218 * t359 - t220 * t286 + t222 * t287) * t266 + (-t217 * t359 - t219 * t286 + t221 * t287) * t265 + (-t243 * t359 - t244 * t286 + t245 * t287) * t284) / 0.2e1 + t300 * ((t248 * t360 - t250 * t290 + t252 * t291) * t300 + (t247 * t360 - t249 * t290 + t251 * t291) * t299 + (t277 * t360 - t278 * t290 + t279 * t291) * t301) / 0.2e1 + m(4) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(5) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(6) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(7) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + ((t232 * t373 + t233 * t371 + t267 * t372) * t236 + (t232 * t378 + t233 * t374 + t267 * t376) * t203 + (t379 * t232 + t375 * t233 + t377 * t267) * t202) * t202 / 0.2e1 + ((t234 * t373 + t235 * t371 + t269 * t372) * t236 + (t378 * t234 + t374 * t235 + t376 * t269) * t203 + (t234 * t379 + t375 * t235 + t377 * t269) * t202) * t203 / 0.2e1 + ((t259 * t373 + t260 * t371 + t286 * t372) * t236 + (t259 * t378 + t260 * t374 + t286 * t376) * t203 + (t259 * t379 + t375 * t260 + t377 * t286) * t202) * t236 / 0.2e1 + ((-t304 * t326 + t306 * t329 + Icges(1,4)) * V_base(5) + (-t305 * t326 + t307 * t329 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t304 * t329 + t306 * t326 + Icges(1,2)) * V_base(5) + (t305 * t329 + t307 * t326 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t326 + Icges(2,6) * t329) * V_base(5) + (Icges(2,5) * t329 - Icges(2,6) * t326) * V_base(4) + Icges(2,3) * t317 / 0.2e1) * t317;
T  = t1;
