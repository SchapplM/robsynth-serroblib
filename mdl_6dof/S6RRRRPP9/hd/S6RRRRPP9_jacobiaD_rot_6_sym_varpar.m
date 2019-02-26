% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:16
% EndTime: 2019-02-26 22:30:18
% DurationCPUTime: 2.37s
% Computational Cost: add. (9759->192), mult. (28191->370), div. (951->12), fcn. (35841->13), ass. (0->146)
t302 = sin(qJ(1));
t298 = cos(pkin(6));
t325 = qJD(2) * t298 + qJD(1);
t301 = sin(qJ(2));
t355 = t301 * t302;
t333 = t298 * t355;
t348 = qJD(2) * t301;
t305 = cos(qJ(2));
t306 = cos(qJ(1));
t350 = t305 * t306;
t262 = -qJD(1) * t333 - t302 * t348 + t325 * t350;
t300 = sin(qJ(3));
t304 = cos(qJ(3));
t353 = t302 * t305;
t354 = t301 * t306;
t285 = t298 * t354 + t353;
t297 = sin(pkin(6));
t357 = t297 * t306;
t317 = t285 * t300 + t304 * t357;
t349 = qJD(1) * t297;
t331 = t302 * t349;
t242 = t317 * qJD(3) - t262 * t304 - t300 * t331;
t334 = t300 * t357;
t277 = -t285 * t304 + t334;
t332 = t298 * t350;
t284 = -t332 + t355;
t299 = sin(qJ(4));
t303 = cos(qJ(4));
t253 = t277 * t303 - t284 * t299;
t286 = t298 * t353 + t354;
t261 = t286 * qJD(1) + t285 * qJD(2);
t226 = -t253 * qJD(4) - t242 * t299 - t261 * t303;
t364 = t277 * t299;
t386 = -t284 * t303 - t364;
t287 = -t333 + t350;
t359 = t297 * t300;
t279 = t287 * t304 + t302 * t359;
t254 = -t279 * t299 + t286 * t303;
t260 = -t285 * qJD(1) - t286 * qJD(2);
t330 = t306 * t349;
t238 = t279 * qJD(3) + t260 * t300 - t304 * t330;
t358 = t297 * t304;
t316 = -t287 * t300 + t302 * t358;
t268 = 0.1e1 / t316;
t269 = 0.1e1 / t316 ^ 2;
t270 = t268 * t269;
t372 = t238 * t270;
t384 = 0.2e1 * t254 * t372;
t282 = t298 * t304 - t301 * t359;
t347 = qJD(2) * t305;
t328 = t297 * t347;
t272 = t282 * qJD(3) + t304 * t328;
t283 = t298 * t300 + t301 * t358;
t351 = t303 * t305;
t273 = -t283 * t299 - t297 * t351;
t329 = t297 * t348;
t244 = t273 * qJD(4) + t272 * t303 + t299 * t329;
t356 = t299 * t305;
t318 = -t283 * t303 + t297 * t356;
t266 = 0.1e1 / t318 ^ 2;
t383 = t244 * t266;
t265 = 0.1e1 / t318;
t369 = t253 * t266;
t319 = -t265 * t317 - t282 * t369;
t382 = t303 * t319;
t236 = atan2(t253, -t318);
t229 = sin(t236);
t230 = cos(t236);
t223 = t229 * t253 - t230 * t318;
t220 = 0.1e1 / t223;
t221 = 0.1e1 / t223 ^ 2;
t381 = 0.2e1 * t253;
t361 = t286 * t299;
t255 = t279 * t303 + t361;
t380 = 0.2e1 * t255;
t248 = t255 ^ 2;
t219 = t221 * t248 + 0.1e1;
t239 = t316 * qJD(3) + t260 * t304 + t300 * t330;
t259 = -qJD(1) * t332 - t306 * t347 + t325 * t355;
t344 = qJD(4) * t299;
t225 = -t259 * t299 - t279 * t344 + (qJD(4) * t286 + t239) * t303;
t375 = t221 * t255;
t246 = t253 ^ 2;
t235 = t246 * t266 + 0.1e1;
t231 = 0.1e1 / t235;
t366 = t261 * t299;
t227 = t366 + qJD(4) * t364 + (qJD(4) * t284 - t242) * t303;
t322 = t227 * t265 - t244 * t369;
t211 = t322 * t231;
t324 = t229 * t318 + t230 * t253;
t206 = t324 * t211 - t229 * t227 + t230 * t244;
t222 = t220 * t221;
t378 = t206 * t222;
t379 = (t225 * t375 - t248 * t378) / t219 ^ 2;
t371 = t265 * t383;
t377 = (-t227 * t369 + t246 * t371) / t235 ^ 2;
t376 = t221 * t225;
t374 = t229 * t255;
t373 = t230 * t255;
t370 = t253 * t265;
t368 = t254 * t269;
t367 = t254 * t279;
t363 = t316 * t303;
t360 = t286 * t300;
t352 = t303 * t304;
t346 = qJD(3) * t300;
t345 = qJD(3) * t304;
t343 = qJD(4) * t304;
t342 = 0.2e1 * t379;
t224 = -t255 * qJD(4) - t239 * t299 - t259 * t303;
t247 = t254 ^ 2;
t237 = t247 * t269 + 0.1e1;
t341 = 0.2e1 * (t224 * t368 + t247 * t372) / t237 ^ 2;
t340 = -0.2e1 * t377;
t339 = t222 * t380;
t338 = t265 * t377;
t337 = t221 * t374;
t336 = t221 * t373;
t327 = t206 * t339;
t326 = t371 * t381;
t323 = -qJD(4) * t287 - t259 * t304;
t321 = -t265 * t386 - t273 * t369;
t256 = -t284 * t352 + t285 * t299;
t280 = (t299 * t301 + t304 * t351) * t297;
t320 = t256 * t265 - t280 * t369;
t313 = -t229 + (t230 * t370 + t229) * t231;
t240 = qJD(3) * t334 - t262 * t300 - t285 * t345 + t304 * t331;
t271 = -t283 * qJD(3) - t300 * t328;
t258 = -t286 * t352 + t287 * t299;
t257 = t287 * t303 + t304 * t361;
t245 = ((qJD(2) - t343) * t356 + (-t305 * t346 + (-qJD(2) * t304 + qJD(4)) * t301) * t303) * t297;
t243 = t318 * qJD(4) - t272 * t299 + t303 * t329;
t233 = 0.1e1 / t237;
t228 = (t284 * t343 + t262) * t299 + (qJD(4) * t285 - t261 * t304 + t284 * t346) * t303;
t217 = 0.1e1 / t219;
t216 = t231 * t382;
t215 = t320 * t231;
t214 = t321 * t231;
t210 = t313 * t255;
t209 = (t229 * t317 + t230 * t282) * t303 + t324 * t216;
t208 = t324 * t215 - t229 * t256 + t230 * t280;
t207 = t324 * t214 + t229 * t386 + t230 * t273;
t205 = t320 * t340 + (-t280 * t326 + t228 * t265 + (t227 * t280 + t244 * t256 - t245 * t253) * t266) * t231;
t203 = t321 * t340 + (-t273 * t326 - t226 * t265 + (t227 * t273 - t243 * t253 - t244 * t386) * t266) * t231;
t202 = t340 * t382 + (-t319 * t344 + (-t282 * t326 + t240 * t265 + (t227 * t282 - t244 * t317 - t253 * t271) * t266) * t303) * t231;
t1 = [-t338 * t380 + (t225 * t265 + t255 * t383) * t231, t205, t202, t203, 0, 0; -0.2e1 * t253 * t220 * t379 + ((qJD(4) * t386 + t242 * t303 - t366) * t220 + (-t206 * t253 - t210 * t225) * t221) * t217 + (t210 * t221 * t342 + (0.2e1 * t210 * t378 - (-t211 * t231 * t370 + t340) * t337 - (-t338 * t381 - t211 + (t211 - t322) * t231) * t336 - t313 * t376) * t217) * t255 (t208 * t375 - t220 * t258) * t342 + (t208 * t327 + (-t258 * t206 - t208 * t225 - (t205 * t253 - t215 * t227 + t245 + (t215 * t318 - t256) * t211) * t373 - (t205 * t318 - t215 * t244 - t228 + (-t215 * t253 - t280) * t211) * t374) * t221 + ((t286 * t343 + t260) * t299 + (t286 * t346 - t323) * t303) * t220) * t217 (t209 * t375 - t220 * t363) * t342 + (-t209 * t376 + (-t238 * t303 - t316 * t344) * t220 + (t209 * t339 - t221 * t363) * t206 - (-t282 * t344 + t202 * t253 - t216 * t227 + t271 * t303 + (t216 * t318 + t303 * t317) * t211) * t336 - (-t317 * t344 + t202 * t318 - t216 * t244 - t240 * t303 + (-t216 * t253 - t282 * t303) * t211) * t337) * t217 (t207 * t375 - t220 * t254) * t342 + (t207 * t327 + t224 * t220 + (-t254 * t206 - t207 * t225 - (t203 * t253 - t214 * t227 + t243 + (t214 * t318 + t386) * t211) * t373 - (t203 * t318 - t214 * t244 + t226 + (-t214 * t253 - t273) * t211) * t374) * t221) * t217, 0, 0; (t268 * t386 - t317 * t368) * t341 + (-t226 * t268 + t317 * t384 + (t224 * t317 - t238 * t386 - t240 * t254) * t269) * t233 (t257 * t268 - t360 * t368) * t341 + (-(t260 * t303 + t323 * t299) * t268 + (-(-t299 * t346 + t303 * t343) * t268 + t300 * t384) * t286 + (t224 * t360 - t257 * t238 + (-t259 * t300 + t286 * t345) * t254) * t269) * t233 (-t268 * t299 * t316 + t269 * t367) * t341 + (qJD(4) * t268 * t363 + (-t224 * t279 - t239 * t254) * t269 + (-0.2e1 * t270 * t367 + (t269 * t316 - t268) * t299) * t238) * t233, -t255 * t268 * t341 + (t238 * t255 * t269 + t225 * t268) * t233, 0, 0;];
JaD_rot  = t1;
