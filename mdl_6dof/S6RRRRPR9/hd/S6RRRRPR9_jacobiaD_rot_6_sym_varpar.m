% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR9
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:13
% EndTime: 2019-02-26 22:35:15
% DurationCPUTime: 1.75s
% Computational Cost: add. (12483->154), mult. (18084->301), div. (989->12), fcn. (22873->13), ass. (0->135)
t310 = qJ(3) + qJ(4);
t306 = sin(t310);
t312 = cos(pkin(6));
t314 = cos(qJ(2));
t384 = sin(qJ(1));
t343 = t384 * t314;
t313 = sin(qJ(2));
t315 = cos(qJ(1));
t362 = t315 * t313;
t326 = -t312 * t362 - t343;
t307 = cos(t310);
t311 = sin(pkin(6));
t363 = t311 * t315;
t348 = t307 * t363;
t277 = -t306 * t326 + t348;
t365 = t311 * t313;
t288 = t306 * t365 - t307 * t312;
t266 = atan2(-t277, t288);
t261 = sin(t266);
t262 = cos(t266);
t244 = -t261 * t277 + t262 * t288;
t242 = 0.1e1 / t244 ^ 2;
t344 = t384 * t313;
t336 = t312 * t344;
t361 = t315 * t314;
t293 = -t336 + t361;
t345 = t311 * t384;
t282 = t293 * t306 - t307 * t345;
t276 = t282 ^ 2;
t240 = t242 * t276 + 0.1e1;
t325 = -t312 * t343 - t362;
t272 = t326 * qJD(1) + t325 * qJD(2);
t309 = qJD(3) + qJD(4);
t332 = t309 * t345 + t272;
t342 = qJD(1) * t363;
t366 = t307 * t309;
t248 = t293 * t366 + t332 * t306 - t307 * t342;
t377 = t242 * t282;
t275 = t277 ^ 2;
t286 = 0.1e1 / t288 ^ 2;
t265 = t275 * t286 + 0.1e1;
t263 = 0.1e1 / t265;
t341 = qJD(2) * t384;
t274 = -qJD(1) * t336 - t313 * t341 + (qJD(2) * t312 + qJD(1)) * t361;
t299 = t306 * t363;
t340 = t384 * qJD(1);
t335 = t311 * t340;
t250 = t274 * t306 - t309 * t299 - t307 * t335 - t326 * t366;
t359 = qJD(2) * t314;
t328 = t309 * t312 + t311 * t359;
t347 = t309 * t365;
t269 = t328 * t306 + t307 * t347;
t285 = 0.1e1 / t288;
t369 = t277 * t286;
t331 = -t250 * t285 + t269 * t369;
t232 = t331 * t263;
t333 = -t261 * t288 - t262 * t277;
t227 = t333 * t232 - t261 * t250 + t262 * t269;
t241 = 0.1e1 / t244;
t243 = t241 * t242;
t382 = t227 * t243;
t357 = 0.2e1 * (t248 * t377 - t276 * t382) / t240 ^ 2;
t388 = t269 * t286;
t346 = t312 * t361;
t290 = -t344 + t346;
t364 = t311 * t314;
t327 = -t285 * t290 + t364 * t369;
t387 = t306 * t327;
t251 = (t309 * t326 + t335) * t306 + t274 * t307 - t309 * t348;
t283 = t293 * t307 + t306 * t345;
t308 = pkin(12) + qJ(6);
t304 = sin(t308);
t305 = cos(t308);
t260 = t283 * t305 - t304 * t325;
t254 = 0.1e1 / t260;
t255 = 0.1e1 / t260 ^ 2;
t386 = -0.2e1 * t277;
t385 = 0.2e1 * t282;
t249 = t332 * t307 + (-t293 * t309 + t342) * t306;
t271 = -qJD(1) * t346 - t315 * t359 + (t312 * t341 + t340) * t313;
t236 = t260 * qJD(6) + t249 * t304 + t271 * t305;
t259 = t283 * t304 + t305 * t325;
t253 = t259 ^ 2;
t247 = t253 * t255 + 0.1e1;
t375 = t255 * t259;
t358 = qJD(6) * t259;
t237 = t249 * t305 - t271 * t304 - t358;
t379 = t237 * t254 * t255;
t381 = (t236 * t375 - t253 * t379) / t247 ^ 2;
t371 = t285 * t388;
t380 = (t250 * t369 - t275 * t371) / t265 ^ 2;
t378 = t242 * t248;
t376 = t254 * t304;
t374 = t259 * t305;
t373 = t261 * t282;
t372 = t262 * t282;
t370 = t277 * t285;
t368 = t325 * t306;
t367 = t325 * t307;
t360 = qJD(2) * t313;
t356 = -0.2e1 * t381;
t355 = 0.2e1 * t381;
t354 = -0.2e1 * t380;
t353 = t243 * t385;
t352 = t285 * t380;
t351 = t259 * t379;
t350 = t242 * t373;
t349 = t242 * t372;
t339 = 0.2e1 * t351;
t338 = t371 * t386;
t279 = -t307 * t326 - t299;
t334 = -qJD(6) * t367 + t272;
t258 = -t279 * t305 + t290 * t304;
t257 = -t279 * t304 - t290 * t305;
t330 = t255 * t374 - t376;
t289 = t306 * t312 + t307 * t365;
t329 = -t279 * t285 + t289 * t369;
t323 = -t261 + (t262 * t370 + t261) * t263;
t322 = qJD(6) * t293 + t271 * t307 - t309 * t368;
t273 = t325 * qJD(1) + t326 * qJD(2);
t270 = -t306 * t347 + t328 * t307;
t268 = t293 * t304 + t305 * t367;
t267 = -t293 * t305 + t304 * t367;
t245 = 0.1e1 / t247;
t238 = 0.1e1 / t240;
t235 = t263 * t387;
t234 = t329 * t263;
t231 = t323 * t282;
t229 = (-t261 * t290 + t262 * t364) * t306 + t333 * t235;
t228 = t333 * t234 - t261 * t279 + t262 * t289;
t225 = t329 * t354 + (t289 * t338 - t251 * t285 + (t250 * t289 + t269 * t279 + t270 * t277) * t286) * t263;
t224 = t354 * t387 + (t327 * t366 + (t338 * t364 - t273 * t285 + (t269 * t290 + (t250 * t314 - t277 * t360) * t311) * t286) * t306) * t263;
t223 = t330 * t282 * t356 + (t330 * t248 + ((-qJD(6) * t254 - 0.2e1 * t351) * t305 + (t236 * t305 + (t237 - t358) * t304) * t255) * t282) * t245;
t222 = (t228 * t377 - t241 * t283) * t357 + (t228 * t227 * t353 + t249 * t241 + (-t283 * t227 - t228 * t248 - (-t225 * t277 - t234 * t250 + t270 + (-t234 * t288 - t279) * t232) * t372 - (-t225 * t288 - t234 * t269 - t251 + (t234 * t277 - t289) * t232) * t373) * t242) * t238;
t1 = [t352 * t385 + (-t248 * t285 + t282 * t388) * t263, t224, t225, t225, 0, 0; t277 * t241 * t357 + (-t250 * t241 + (t227 * t277 - t231 * t248) * t242) * t238 + (t231 * t242 * t357 + (0.2e1 * t231 * t382 - (-t232 * t263 * t370 + t354) * t350 - (t352 * t386 - t232 + (t232 - t331) * t263) * t349 - t323 * t378) * t238) * t282 (t229 * t377 - t241 * t368) * t357 + (-t229 * t378 + (t271 * t306 + t325 * t366) * t241 + (t229 * t353 - t242 * t368) * t227 - (-t224 * t277 - t235 * t250 + (-t306 * t360 + t314 * t366) * t311 + (-t235 * t288 - t290 * t306) * t232) * t349 - (-t290 * t366 - t224 * t288 - t235 * t269 - t273 * t306 + (t235 * t277 - t306 * t364) * t232) * t350) * t238, t222, t222, 0, 0; (-t254 * t257 + t258 * t375) * t355 + ((t258 * qJD(6) - t251 * t304 - t273 * t305) * t254 + t258 * t339 + (-t257 * t237 - (-t257 * qJD(6) - t251 * t305 + t273 * t304) * t259 - t258 * t236) * t255) * t245 (-t254 * t267 + t268 * t375) * t355 + (t268 * t339 - t334 * t254 * t305 + t322 * t376 + (-t334 * t259 * t304 - t268 * t236 - t267 * t237 - t322 * t374) * t255) * t245, t223, t223, 0, t356 + 0.2e1 * (t236 * t245 * t255 + (-t245 * t379 - t255 * t381) * t259) * t259;];
JaD_rot  = t1;
