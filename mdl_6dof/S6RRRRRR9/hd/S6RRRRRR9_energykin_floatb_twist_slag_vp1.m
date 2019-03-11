% Calculate kinetic energy for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:37
% EndTime: 2019-03-10 05:16:42
% DurationCPUTime: 4.56s
% Computational Cost: add. (5955->447), mult. (15086->675), div. (0->0), fcn. (19332->16), ass. (0->199)
t411 = cos(qJ(3));
t410 = cos(qJ(4));
t360 = cos(pkin(6));
t409 = pkin(9) * t360;
t366 = cos(qJ(5));
t408 = pkin(5) * t366;
t406 = cos(pkin(7));
t405 = sin(pkin(7));
t365 = sin(qJ(1));
t404 = Icges(2,4) * t365;
t367 = cos(qJ(2));
t368 = cos(qJ(1));
t394 = t367 * t368;
t364 = sin(qJ(2));
t397 = t364 * t365;
t325 = t360 * t394 - t397;
t395 = t365 * t367;
t396 = t364 * t368;
t326 = t360 * t396 + t395;
t363 = sin(qJ(3));
t359 = sin(pkin(6));
t386 = t411 * t405;
t385 = t359 * t386;
t387 = t406 * t411;
t287 = -t325 * t387 + t326 * t363 + t368 * t385;
t361 = sin(qJ(5));
t403 = t287 * t361;
t327 = -t360 * t395 - t396;
t328 = -t360 * t397 + t394;
t289 = -t327 * t387 + t328 * t363 - t365 * t385;
t402 = t289 * t361;
t400 = t359 * t364;
t309 = -t359 * t367 * t387 - t360 * t386 + t363 * t400;
t401 = t309 * t361;
t399 = t359 * t365;
t398 = t359 * t368;
t393 = qJD(2) * t359;
t392 = V_base(5) * pkin(8) + V_base(1);
t338 = t365 * t393 + V_base(4);
t353 = V_base(6) + qJD(1);
t389 = t359 * t406;
t388 = t359 * t405;
t380 = -t327 * t405 + t365 * t389;
t302 = qJD(3) * t380 + t338;
t339 = qJD(2) * t360 + t353;
t263 = qJD(4) * t289 + t302;
t379 = t360 * t406 - t367 * t388;
t311 = qJD(3) * t379 + t339;
t337 = -t368 * t393 + V_base(5);
t290 = t328 * t411 + (t327 * t406 + t365 * t388) * t363;
t362 = sin(qJ(4));
t272 = t290 * t362 - t380 * t410;
t228 = qJD(5) * t272 + t263;
t330 = pkin(1) * t365 - pkin(9) * t398;
t384 = -t330 * t353 + V_base(5) * t409 + t392;
t279 = qJD(4) * t309 + t311;
t331 = pkin(1) * t368 + pkin(9) * t399;
t383 = V_base(4) * t330 - t331 * V_base(5) + V_base(3);
t381 = -t325 * t405 - t368 * t389;
t301 = qJD(3) * t381 + t337;
t310 = t360 * t405 * t363 + (t406 * t363 * t367 + t364 * t411) * t359;
t285 = t310 * t362 - t379 * t410;
t255 = qJD(5) * t285 + t279;
t262 = qJD(4) * t287 + t301;
t382 = t353 * t331 + V_base(2) + (-pkin(8) - t409) * V_base(4);
t288 = t326 * t411 + (t325 * t406 - t368 * t388) * t363;
t270 = t288 * t362 - t381 * t410;
t227 = qJD(5) * t270 + t262;
t291 = t326 * pkin(2) + pkin(10) * t381;
t315 = pkin(2) * t400 + pkin(10) * t379;
t378 = -t291 * t339 + t337 * t315 + t384;
t292 = t328 * pkin(2) + pkin(10) * t380;
t377 = t338 * t291 - t292 * t337 + t383;
t376 = t339 * t292 - t315 * t338 + t382;
t259 = pkin(3) * t288 + pkin(11) * t287;
t278 = pkin(3) * t310 + pkin(11) * t309;
t375 = -t259 * t311 + t301 * t278 + t378;
t260 = pkin(3) * t290 + pkin(11) * t289;
t374 = t302 * t259 - t260 * t301 + t377;
t373 = t311 * t260 - t278 * t302 + t376;
t271 = t288 * t410 + t362 * t381;
t229 = t271 * pkin(4) + t270 * pkin(12);
t286 = t310 * t410 + t362 * t379;
t258 = t286 * pkin(4) + t285 * pkin(12);
t372 = -t229 * t279 + t262 * t258 + t375;
t273 = t290 * t410 + t362 * t380;
t230 = t273 * pkin(4) + t272 * pkin(12);
t371 = t263 * t229 - t230 * t262 + t374;
t370 = t279 * t230 - t258 * t263 + t373;
t358 = qJ(5) + qJ(6);
t356 = Icges(2,4) * t368;
t355 = cos(t358);
t354 = sin(t358);
t347 = rSges(2,1) * t368 - rSges(2,2) * t365;
t346 = rSges(2,1) * t365 + rSges(2,2) * t368;
t345 = Icges(2,1) * t368 - t404;
t344 = Icges(2,1) * t365 + t356;
t343 = -Icges(2,2) * t365 + t356;
t342 = Icges(2,2) * t368 + t404;
t336 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t335 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t334 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t321 = rSges(3,3) * t360 + (rSges(3,1) * t364 + rSges(3,2) * t367) * t359;
t320 = Icges(3,5) * t360 + (Icges(3,1) * t364 + Icges(3,4) * t367) * t359;
t319 = Icges(3,6) * t360 + (Icges(3,4) * t364 + Icges(3,2) * t367) * t359;
t318 = Icges(3,3) * t360 + (Icges(3,5) * t364 + Icges(3,6) * t367) * t359;
t314 = V_base(5) * rSges(2,3) - t346 * t353 + t392;
t313 = t347 * t353 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t312 = t346 * V_base(4) - t347 * V_base(5) + V_base(3);
t300 = rSges(3,1) * t328 + rSges(3,2) * t327 + rSges(3,3) * t399;
t299 = rSges(3,1) * t326 + rSges(3,2) * t325 - rSges(3,3) * t398;
t298 = Icges(3,1) * t328 + Icges(3,4) * t327 + Icges(3,5) * t399;
t297 = Icges(3,1) * t326 + Icges(3,4) * t325 - Icges(3,5) * t398;
t296 = Icges(3,4) * t328 + Icges(3,2) * t327 + Icges(3,6) * t399;
t295 = Icges(3,4) * t326 + Icges(3,2) * t325 - Icges(3,6) * t398;
t294 = Icges(3,5) * t328 + Icges(3,6) * t327 + Icges(3,3) * t399;
t293 = Icges(3,5) * t326 + Icges(3,6) * t325 - Icges(3,3) * t398;
t277 = t310 * rSges(4,1) - t309 * rSges(4,2) + rSges(4,3) * t379;
t276 = Icges(4,1) * t310 - Icges(4,4) * t309 + Icges(4,5) * t379;
t275 = Icges(4,4) * t310 - Icges(4,2) * t309 + Icges(4,6) * t379;
t274 = Icges(4,5) * t310 - Icges(4,6) * t309 + Icges(4,3) * t379;
t267 = t286 * t366 + t401;
t266 = -t286 * t361 + t309 * t366;
t265 = t286 * t355 + t309 * t354;
t264 = -t286 * t354 + t309 * t355;
t257 = -t299 * t339 + t321 * t337 + t384;
t256 = t300 * t339 - t321 * t338 + t382;
t253 = t290 * rSges(4,1) - t289 * rSges(4,2) + rSges(4,3) * t380;
t252 = t288 * rSges(4,1) - t287 * rSges(4,2) + rSges(4,3) * t381;
t251 = t299 * t338 - t300 * t337 + t383;
t250 = Icges(4,1) * t290 - Icges(4,4) * t289 + Icges(4,5) * t380;
t249 = Icges(4,1) * t288 - Icges(4,4) * t287 + Icges(4,5) * t381;
t248 = Icges(4,4) * t290 - Icges(4,2) * t289 + Icges(4,6) * t380;
t247 = Icges(4,4) * t288 - Icges(4,2) * t287 + Icges(4,6) * t381;
t246 = Icges(4,5) * t290 - Icges(4,6) * t289 + Icges(4,3) * t380;
t245 = Icges(4,5) * t288 - Icges(4,6) * t287 + Icges(4,3) * t381;
t244 = rSges(5,1) * t286 - rSges(5,2) * t285 + rSges(5,3) * t309;
t243 = Icges(5,1) * t286 - Icges(5,4) * t285 + Icges(5,5) * t309;
t242 = Icges(5,4) * t286 - Icges(5,2) * t285 + Icges(5,6) * t309;
t241 = Icges(5,5) * t286 - Icges(5,6) * t285 + Icges(5,3) * t309;
t240 = t273 * t366 + t402;
t239 = -t273 * t361 + t289 * t366;
t238 = t271 * t366 + t403;
t237 = -t271 * t361 + t287 * t366;
t235 = t273 * t355 + t289 * t354;
t234 = -t273 * t354 + t289 * t355;
t233 = t271 * t355 + t287 * t354;
t232 = -t271 * t354 + t287 * t355;
t231 = qJD(6) * t285 + t255;
t225 = rSges(5,1) * t273 - rSges(5,2) * t272 + rSges(5,3) * t289;
t224 = rSges(5,1) * t271 - rSges(5,2) * t270 + rSges(5,3) * t287;
t223 = Icges(5,1) * t273 - Icges(5,4) * t272 + Icges(5,5) * t289;
t222 = Icges(5,1) * t271 - Icges(5,4) * t270 + Icges(5,5) * t287;
t221 = Icges(5,4) * t273 - Icges(5,2) * t272 + Icges(5,6) * t289;
t220 = Icges(5,4) * t271 - Icges(5,2) * t270 + Icges(5,6) * t287;
t219 = Icges(5,5) * t273 - Icges(5,6) * t272 + Icges(5,3) * t289;
t218 = Icges(5,5) * t271 - Icges(5,6) * t270 + Icges(5,3) * t287;
t216 = rSges(6,1) * t267 + rSges(6,2) * t266 + rSges(6,3) * t285;
t215 = Icges(6,1) * t267 + Icges(6,4) * t266 + Icges(6,5) * t285;
t214 = Icges(6,4) * t267 + Icges(6,2) * t266 + Icges(6,6) * t285;
t213 = Icges(6,5) * t267 + Icges(6,6) * t266 + Icges(6,3) * t285;
t212 = rSges(7,1) * t265 + rSges(7,2) * t264 + rSges(7,3) * t285;
t211 = Icges(7,1) * t265 + Icges(7,4) * t264 + Icges(7,5) * t285;
t210 = Icges(7,4) * t265 + Icges(7,2) * t264 + Icges(7,6) * t285;
t209 = Icges(7,5) * t265 + Icges(7,6) * t264 + Icges(7,3) * t285;
t208 = pkin(5) * t401 + pkin(13) * t285 + t286 * t408;
t207 = qJD(6) * t272 + t228;
t206 = qJD(6) * t270 + t227;
t204 = rSges(6,1) * t240 + rSges(6,2) * t239 + rSges(6,3) * t272;
t203 = rSges(6,1) * t238 + rSges(6,2) * t237 + rSges(6,3) * t270;
t202 = Icges(6,1) * t240 + Icges(6,4) * t239 + Icges(6,5) * t272;
t201 = Icges(6,1) * t238 + Icges(6,4) * t237 + Icges(6,5) * t270;
t200 = Icges(6,4) * t240 + Icges(6,2) * t239 + Icges(6,6) * t272;
t199 = Icges(6,4) * t238 + Icges(6,2) * t237 + Icges(6,6) * t270;
t198 = Icges(6,5) * t240 + Icges(6,6) * t239 + Icges(6,3) * t272;
t197 = Icges(6,5) * t238 + Icges(6,6) * t237 + Icges(6,3) * t270;
t196 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t272;
t195 = rSges(7,1) * t233 + rSges(7,2) * t232 + rSges(7,3) * t270;
t194 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t272;
t193 = Icges(7,1) * t233 + Icges(7,4) * t232 + Icges(7,5) * t270;
t192 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t272;
t191 = Icges(7,4) * t233 + Icges(7,2) * t232 + Icges(7,6) * t270;
t190 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t272;
t189 = Icges(7,5) * t233 + Icges(7,6) * t232 + Icges(7,3) * t270;
t188 = -t252 * t311 + t277 * t301 + t378;
t187 = t253 * t311 - t277 * t302 + t376;
t186 = pkin(5) * t402 + pkin(13) * t272 + t273 * t408;
t185 = pkin(5) * t403 + pkin(13) * t270 + t271 * t408;
t184 = t252 * t302 - t253 * t301 + t377;
t183 = -t224 * t279 + t244 * t262 + t375;
t182 = t225 * t279 - t244 * t263 + t373;
t181 = t224 * t263 - t225 * t262 + t374;
t180 = -t203 * t255 + t216 * t227 + t372;
t179 = t204 * t255 - t216 * t228 + t370;
t178 = t203 * t228 - t204 * t227 + t371;
t177 = -t185 * t255 - t195 * t231 + t206 * t212 + t208 * t227 + t372;
t176 = t186 * t255 + t196 * t231 - t207 * t212 - t208 * t228 + t370;
t175 = t185 * t228 - t186 * t227 + t195 * t207 - t196 * t206 + t371;
t1 = t311 * ((t246 * t379 - t309 * t248 + t310 * t250) * t302 + (t245 * t379 - t309 * t247 + t310 * t249) * t301 + (t274 * t379 - t309 * t275 + t310 * t276) * t311) / 0.2e1 + t302 * ((t246 * t380 - t289 * t248 + t290 * t250) * t302 + (t245 * t380 - t289 * t247 + t290 * t249) * t301 + (t274 * t380 - t289 * t275 + t290 * t276) * t311) / 0.2e1 + t301 * ((t246 * t381 - t287 * t248 + t288 * t250) * t302 + (t245 * t381 - t287 * t247 + t288 * t249) * t301 + (t274 * t381 - t287 * t275 + t288 * t276) * t311) / 0.2e1 + ((Icges(2,5) * t365 + Icges(2,6) * t368) * V_base(5) + (Icges(2,5) * t368 - Icges(2,6) * t365) * V_base(4) + Icges(2,3) * t353 / 0.2e1) * t353 + ((t342 * t368 + t344 * t365 + Icges(1,2)) * V_base(5) + (t343 * t368 + t345 * t365 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(1) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1 + m(2) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + t279 * ((t219 * t309 - t221 * t285 + t223 * t286) * t263 + (t218 * t309 - t220 * t285 + t222 * t286) * t262 + (t241 * t309 - t242 * t285 + t243 * t286) * t279) / 0.2e1 + t262 * ((t219 * t287 - t221 * t270 + t223 * t271) * t263 + (t218 * t287 - t220 * t270 + t222 * t271) * t262 + (t241 * t287 - t242 * t270 + t243 * t271) * t279) / 0.2e1 + t263 * ((t219 * t289 - t221 * t272 + t223 * t273) * t263 + (t218 * t289 - t220 * t272 + t222 * t273) * t262 + (t241 * t289 - t242 * t272 + t243 * t273) * t279) / 0.2e1 + t255 * ((t198 * t285 + t200 * t266 + t202 * t267) * t228 + (t197 * t285 + t199 * t266 + t201 * t267) * t227 + (t213 * t285 + t214 * t266 + t215 * t267) * t255) / 0.2e1 + t231 * ((t190 * t285 + t192 * t264 + t194 * t265) * t207 + (t189 * t285 + t191 * t264 + t193 * t265) * t206 + (t285 * t209 + t264 * t210 + t265 * t211) * t231) / 0.2e1 + t206 * ((t190 * t270 + t192 * t232 + t194 * t233) * t207 + (t270 * t189 + t232 * t191 + t233 * t193) * t206 + (t209 * t270 + t210 * t232 + t211 * t233) * t231) / 0.2e1 + t228 * ((t272 * t198 + t239 * t200 + t240 * t202) * t228 + (t197 * t272 + t199 * t239 + t201 * t240) * t227 + (t213 * t272 + t214 * t239 + t215 * t240) * t255) / 0.2e1 + t207 * ((t272 * t190 + t234 * t192 + t235 * t194) * t207 + (t189 * t272 + t191 * t234 + t193 * t235) * t206 + (t209 * t272 + t210 * t234 + t211 * t235) * t231) / 0.2e1 + t227 * ((t198 * t270 + t200 * t237 + t202 * t238) * t228 + (t270 * t197 + t237 * t199 + t238 * t201) * t227 + (t213 * t270 + t214 * t237 + t215 * t238) * t255) / 0.2e1 + m(3) * (t251 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + ((-t342 * t365 + t344 * t368 + Icges(1,4)) * V_base(5) + (-t343 * t365 + t345 * t368 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + t337 * ((-t294 * t398 + t296 * t325 + t298 * t326) * t338 + (-t293 * t398 + t295 * t325 + t297 * t326) * t337 + (-t318 * t398 + t319 * t325 + t320 * t326) * t339) / 0.2e1 + t338 * ((t294 * t399 + t296 * t327 + t298 * t328) * t338 + (t293 * t399 + t295 * t327 + t297 * t328) * t337 + (t318 * t399 + t319 * t327 + t320 * t328) * t339) / 0.2e1 + m(7) * (t175 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + m(6) * (t178 ^ 2 + t179 ^ 2 + t180 ^ 2) / 0.2e1 + m(5) * (t181 ^ 2 + t182 ^ 2 + t183 ^ 2) / 0.2e1 + m(4) * (t184 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + t339 * ((t293 * t337 + t294 * t338 + t318 * t339) * t360 + ((t296 * t367 + t298 * t364) * t338 + (t295 * t367 + t297 * t364) * t337 + (t319 * t367 + t320 * t364) * t339) * t359) / 0.2e1;
T  = t1;
