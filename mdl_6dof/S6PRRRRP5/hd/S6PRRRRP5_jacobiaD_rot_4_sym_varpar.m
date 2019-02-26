% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP5_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobiaD_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:29
% EndTime: 2019-02-26 20:17:30
% DurationCPUTime: 1.16s
% Computational Cost: add. (5341->130), mult. (16984->252), div. (538->12), fcn. (21480->15), ass. (0->122)
t287 = sin(qJ(2));
t290 = cos(qJ(2));
t354 = cos(pkin(12));
t355 = cos(pkin(6));
t320 = t355 * t354;
t353 = sin(pkin(12));
t306 = -t353 * t287 + t290 * t320;
t271 = t306 * qJD(2);
t282 = sin(pkin(7));
t283 = sin(pkin(6));
t323 = t354 * t283 * t282;
t358 = -qJD(3) * t323 + t271;
t286 = sin(qJ(3));
t305 = -t287 * t320 - t353 * t290;
t300 = t305 * qJD(2);
t289 = cos(qJ(3));
t301 = t306 * t289;
t357 = qJD(3) * t301 + t286 * t300;
t284 = cos(pkin(7));
t302 = t306 * t286;
t297 = t284 * t302 - t289 * t305;
t338 = t284 * t289;
t225 = t297 * qJD(3) + t358 * t286 - t300 * t338;
t343 = t305 * t286;
t249 = -t284 * t301 + t289 * t323 - t343;
t247 = t249 ^ 2;
t334 = t289 * t290;
t337 = t286 * t287;
t311 = t284 * t334 - t337;
t327 = t355 * t282;
t262 = -t311 * t283 - t289 * t327;
t260 = 0.1e1 / t262 ^ 2;
t241 = t247 * t260 + 0.1e1;
t344 = t249 * t260;
t335 = t287 * t289;
t336 = t286 * t290;
t309 = t284 * t336 + t335;
t310 = t284 * t335 + t336;
t321 = qJD(3) * t327;
t245 = t286 * t321 + (t310 * qJD(2) + t309 * qJD(3)) * t283;
t259 = 0.1e1 / t262;
t345 = t245 * t259 * t260;
t356 = -0.2e1 * (t225 * t344 - t247 * t345) / t241 ^ 2;
t242 = atan2(-t249, t262);
t237 = sin(t242);
t238 = cos(t242);
t219 = -t237 * t249 + t238 * t262;
t216 = 0.1e1 / t219;
t319 = t355 * t353;
t304 = t287 * t319 - t354 * t290;
t303 = t354 * t287 + t290 * t319;
t328 = t283 * t353;
t322 = t282 * t328;
t307 = -t284 * t303 + t322;
t253 = t307 * t286 - t289 * t304;
t264 = t282 * t303 + t284 * t328;
t285 = sin(qJ(4));
t288 = cos(qJ(4));
t236 = t253 * t288 + t264 * t285;
t232 = 0.1e1 / t236;
t217 = 0.1e1 / t219 ^ 2;
t233 = 0.1e1 / t236 ^ 2;
t239 = 0.1e1 / t241;
t209 = (-t225 * t259 + t245 * t344) * t239;
t318 = -t237 * t262 - t238 * t249;
t205 = t318 * t209 - t237 * t225 + t238 * t245;
t352 = t205 * t216 * t217;
t272 = t303 * qJD(2);
t273 = t304 * qJD(2);
t339 = t284 * t286;
t342 = t304 * t286;
t228 = t273 * t339 - t272 * t289 + (t307 * t289 + t342) * qJD(3);
t340 = t282 * t288;
t220 = t236 * qJD(4) + t228 * t285 + t273 * t340;
t235 = t253 * t285 - t264 * t288;
t231 = t235 ^ 2;
t224 = t231 * t233 + 0.1e1;
t348 = t233 * t235;
t333 = qJD(4) * t235;
t341 = t282 * t285;
t221 = t228 * t288 - t273 * t341 - t333;
t349 = t221 * t232 * t233;
t351 = (t220 * t348 - t231 * t349) / t224 ^ 2;
t252 = -t289 * t322 + t303 * t338 - t342;
t350 = t217 * t252;
t347 = t237 * t252;
t346 = t238 * t252;
t248 = t252 ^ 2;
t215 = t248 * t217 + 0.1e1;
t227 = t253 * qJD(3) - t272 * t286 - t273 * t338;
t332 = 0.2e1 * (t227 * t350 - t248 * t352) / t215 ^ 2;
t331 = -0.2e1 * t351;
t330 = t235 * t349;
t329 = qJD(3) * t343;
t325 = -0.2e1 * t249 * t345;
t324 = 0.2e1 * t252 * t352;
t315 = -t285 * t232 + t288 * t348;
t251 = -t286 * t323 + t297;
t263 = t309 * t283 + t286 * t327;
t314 = -t251 * t259 + t263 * t344;
t255 = -t305 * t338 + t302;
t270 = t310 * t283;
t313 = -t255 * t259 + t270 * t344;
t257 = -t289 * t303 + t304 * t339;
t312 = -t257 * t285 - t304 * t340;
t244 = t257 * t288 - t304 * t341;
t256 = -t286 * t303 - t304 * t338;
t308 = -t284 * t337 + t334;
t254 = (t311 * qJD(2) + t308 * qJD(3)) * t283;
t246 = t289 * t321 + (t308 * qJD(2) + t311 * qJD(3)) * t283;
t230 = -t256 * qJD(3) + t272 * t339 + t273 * t289;
t229 = t271 * t338 + t284 * t329 + t357;
t226 = t357 * t284 + t358 * t289 + t329;
t222 = 0.1e1 / t224;
t213 = 0.1e1 / t215;
t211 = t313 * t239;
t210 = t314 * t239;
t207 = t318 * t211 - t237 * t255 + t238 * t270;
t206 = t318 * t210 - t237 * t251 + t238 * t263;
t204 = t313 * t356 + (t270 * t325 - t229 * t259 + (t225 * t270 + t245 * t255 + t249 * t254) * t260) * t239;
t203 = t314 * t356 + (t263 * t325 - t226 * t259 + (t225 * t263 + t245 * t251 + t246 * t249) * t260) * t239;
t1 = [0, t204, t203, 0, 0, 0; 0 (t207 * t350 - t216 * t256) * t332 + ((t257 * qJD(3) - t272 * t338 + t273 * t286) * t216 + t207 * t324 + (-t256 * t205 - t207 * t227 - (-t204 * t249 - t211 * t225 + t254 + (-t211 * t262 - t255) * t209) * t346 - (-t204 * t262 - t211 * t245 - t229 + (t211 * t249 - t270) * t209) * t347) * t217) * t213 (t206 * t350 - t216 * t253) * t332 + (t206 * t324 + t228 * t216 + (-t253 * t205 - t206 * t227 - (-t203 * t249 - t210 * t225 + t246 + (-t210 * t262 - t251) * t209) * t346 - (-t203 * t262 - t210 * t245 - t226 + (t210 * t249 - t263) * t209) * t347) * t217) * t213, 0, 0, 0; 0, 0.2e1 * (t232 * t312 + t244 * t348) * t351 + ((t244 * qJD(4) + t230 * t285 + t272 * t340) * t232 + 0.2e1 * t244 * t330 + (t312 * t221 - (t312 * qJD(4) + t230 * t288 - t272 * t341) * t235 - t244 * t220) * t233) * t222, t315 * t252 * t331 + (t315 * t227 + ((-qJD(4) * t232 - 0.2e1 * t330) * t288 + (t220 * t288 + (t221 - t333) * t285) * t233) * t252) * t222, t331 + 0.2e1 * (t220 * t233 * t222 + (-t222 * t349 - t233 * t351) * t235) * t235, 0, 0;];
JaD_rot  = t1;
