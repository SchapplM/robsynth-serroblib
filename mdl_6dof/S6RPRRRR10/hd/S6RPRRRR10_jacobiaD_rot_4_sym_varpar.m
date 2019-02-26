% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR10
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR10_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiaD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:50
% EndTime: 2019-02-26 21:19:52
% DurationCPUTime: 1.34s
% Computational Cost: add. (4837->111), mult. (15324->224), div. (448->12), fcn. (19446->15), ass. (0->108)
t347 = sin(pkin(13));
t352 = cos(pkin(6));
t320 = t352 * t347;
t350 = cos(pkin(13));
t353 = sin(qJ(1));
t355 = cos(qJ(1));
t276 = t355 * t320 + t353 * t350;
t287 = sin(qJ(3));
t348 = sin(pkin(7));
t349 = sin(pkin(6));
t318 = t349 * t348;
t311 = t355 * t318;
t283 = t287 * t311;
t322 = t352 * t350;
t275 = -t355 * t322 + t353 * t347;
t351 = cos(pkin(7));
t326 = t275 * t351;
t354 = cos(qJ(3));
t362 = -t276 * t354 + t287 * t326 + t283;
t304 = t354 * t311;
t324 = t351 * t354;
t361 = -t275 * t324 - t304;
t249 = t276 * t287 - t361;
t317 = t349 * t347;
t319 = t351 * t349;
t359 = t350 * t319 + t352 * t348;
t295 = -t287 * t317 + t359 * t354;
t246 = atan2(-t249, -t295);
t241 = sin(t246);
t242 = cos(t246);
t222 = -t241 * t249 - t242 * t295;
t220 = 0.1e1 / t222 ^ 2;
t277 = -t353 * t320 + t355 * t350;
t300 = t353 * t322 + t355 * t347;
t298 = t300 * t351;
t308 = t353 * t318;
t303 = t354 * t308;
t358 = -t277 * t287 - t354 * t298 + t303;
t248 = t358 ^ 2;
t218 = t248 * t220 + 0.1e1;
t255 = t277 * t354 + (-t298 + t308) * t287;
t272 = t276 * qJD(1);
t297 = qJD(1) * t326;
t228 = -qJD(1) * t304 + t255 * qJD(3) - t272 * t287 - t354 * t297;
t339 = t228 * t220;
t219 = 0.1e1 / t222;
t247 = t249 ^ 2;
t262 = 0.1e1 / t295 ^ 2;
t245 = t247 * t262 + 0.1e1;
t243 = 0.1e1 / t245;
t261 = 0.1e1 / t295;
t265 = t359 * t287 + t354 * t317;
t257 = t265 * qJD(3);
t335 = t257 * t262;
t273 = t300 * qJD(1);
t274 = t277 * qJD(1);
t357 = qJD(1) * t303 + t362 * qJD(3) - t273 * t324 - t274 * t287;
t314 = t249 * t335 - t261 * t357;
t213 = t314 * t243;
t316 = t241 * t295 - t242 * t249;
t209 = t316 * t213 + t241 * t357 + t242 * t257;
t360 = t209 * t220;
t345 = t219 * t360;
t331 = 0.2e1 * (-t248 * t345 - t339 * t358) / t218 ^ 2;
t231 = (qJD(1) * t308 - qJD(3) * t276 - t351 * t273) * t287 + t274 * t354 + t361 * qJD(3);
t309 = t353 * t319;
t267 = t300 * t348 + t309;
t286 = sin(qJ(4));
t288 = cos(qJ(4));
t240 = t255 * t288 + t267 * t286;
t234 = 0.1e1 / t240;
t235 = 0.1e1 / t240 ^ 2;
t356 = -0.2e1 * t358;
t229 = qJD(1) * t283 + t358 * qJD(3) - t272 * t354 + t287 * t297;
t266 = -t275 * t348 + t355 * t319;
t258 = t266 * qJD(1);
t223 = t240 * qJD(4) + t229 * t286 - t258 * t288;
t239 = t255 * t286 - t267 * t288;
t233 = t239 ^ 2;
t227 = t233 * t235 + 0.1e1;
t338 = t235 * t239;
t332 = qJD(4) * t239;
t224 = t229 * t288 + t258 * t286 - t332;
t340 = t224 * t234 * t235;
t344 = (t223 * t338 - t233 * t340) / t227 ^ 2;
t334 = t261 * t335;
t342 = (-t249 * t262 * t357 + t247 * t334) / t245 ^ 2;
t341 = t220 * t358;
t337 = t249 * t261;
t336 = t249 * t265;
t330 = -0.2e1 * t344;
t329 = -0.2e1 * t342;
t328 = t261 * t342;
t327 = t239 * t340;
t325 = t345 * t356;
t238 = t266 * t286 + t288 * t362;
t237 = -t266 * t288 + t286 * t362;
t313 = -t286 * t234 + t288 * t338;
t312 = -t261 * t362 + t262 * t336;
t306 = -t241 + (-t242 * t337 + t241) * t243;
t259 = -qJD(1) * t309 - t273 * t348;
t256 = t295 * qJD(3);
t225 = 0.1e1 / t227;
t216 = 0.1e1 / t218;
t214 = t312 * t243;
t210 = t316 * t214 + t241 * t362 + t242 * t265;
t208 = t312 * t329 + (0.2e1 * t334 * t336 + t231 * t261 + (t249 * t256 - t257 * t362 - t265 * t357) * t262) * t243;
t1 = [-t328 * t356 + (t228 * t261 - t335 * t358) * t243, 0, t208, 0, 0, 0; t249 * t219 * t331 + (t357 * t219 + t249 * t360 + (t306 * t228 - ((t213 * t243 * t337 + t329) * t241 + (0.2e1 * t249 * t328 - t213 + (t213 - t314) * t243) * t242) * t358) * t341) * t216 - (-t341 * t331 + (-t339 + t325) * t216) * t306 * t358, 0 (-t210 * t341 - t219 * t255) * t331 + (t210 * t325 + t229 * t219 + (-t255 * t209 - t210 * t228 - (-(-t208 * t249 + t214 * t357 + t256 + (t214 * t295 + t362) * t213) * t242 - (t208 * t295 - t214 * t257 - t231 + (t214 * t249 - t265) * t213) * t241) * t358) * t220) * t216, 0, 0, 0; 0.2e1 * (-t234 * t237 + t238 * t338) * t344 + ((t238 * qJD(4) - t231 * t286 - t259 * t288) * t234 + 0.2e1 * t238 * t327 + (-t237 * t224 - (-t237 * qJD(4) - t231 * t288 + t259 * t286) * t239 - t238 * t223) * t235) * t225, 0, -t313 * t358 * t330 + (t313 * t228 - ((-qJD(4) * t234 - 0.2e1 * t327) * t288 + (t223 * t288 + (t224 - t332) * t286) * t235) * t358) * t225, t330 + 0.2e1 * (t223 * t235 * t225 + (-t225 * t340 - t235 * t344) * t239) * t239, 0, 0;];
JaD_rot  = t1;
