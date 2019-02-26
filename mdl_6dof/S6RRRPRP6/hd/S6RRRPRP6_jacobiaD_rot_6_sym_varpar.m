% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:04
% EndTime: 2019-02-26 22:12:06
% DurationCPUTime: 1.48s
% Computational Cost: add. (8428->149), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->128)
t273 = cos(pkin(6));
t275 = sin(qJ(2));
t347 = sin(qJ(1));
t308 = t347 * t275;
t298 = t273 * t308;
t303 = qJD(2) * t347;
t277 = cos(qJ(2));
t278 = cos(qJ(1));
t324 = t278 * t277;
t272 = sin(pkin(6));
t327 = t272 * t278;
t352 = -qJD(1) * t298 - t275 * t303 + (qJD(2) * t273 + qJD(1)) * t324 - qJD(3) * t327;
t271 = qJ(3) + pkin(11);
t269 = sin(t271);
t270 = cos(t271);
t307 = t347 * t277;
t325 = t278 * t275;
t290 = -t273 * t325 - t307;
t242 = -t269 * t290 + t270 * t327;
t329 = t272 * t275;
t252 = t269 * t329 - t273 * t270;
t231 = atan2(-t242, t252);
t218 = sin(t231);
t219 = cos(t231);
t209 = -t218 * t242 + t219 * t252;
t207 = 0.1e1 / t209 ^ 2;
t258 = -t298 + t324;
t309 = t272 * t347;
t289 = -t258 * t269 + t270 * t309;
t241 = t289 ^ 2;
t203 = t241 * t207 + 0.1e1;
t288 = -t273 * t307 - t325;
t235 = t290 * qJD(1) + t288 * qJD(2);
t248 = t258 * t270 + t269 * t309;
t306 = qJD(1) * t327;
t213 = t248 * qJD(3) + t235 * t269 - t270 * t306;
t340 = t207 * t289;
t240 = t242 ^ 2;
t250 = 0.1e1 / t252 ^ 2;
t230 = t240 * t250 + 0.1e1;
t224 = 0.1e1 / t230;
t302 = t347 * qJD(1);
t297 = t272 * t302;
t321 = qJD(3) * t270;
t215 = t269 * t352 - t270 * t297 - t290 * t321;
t253 = t273 * t269 + t270 * t329;
t322 = qJD(2) * t277;
t305 = t272 * t322;
t238 = t253 * qJD(3) + t269 * t305;
t249 = 0.1e1 / t252;
t333 = t242 * t250;
t294 = -t215 * t249 + t238 * t333;
t197 = t294 * t224;
t295 = -t218 * t252 - t219 * t242;
t192 = t295 * t197 - t218 * t215 + t219 * t238;
t206 = 0.1e1 / t209;
t208 = t206 * t207;
t345 = t192 * t208;
t319 = 0.2e1 * (-t213 * t340 - t241 * t345) / t203 ^ 2;
t351 = t238 * t250;
t310 = t273 * t324;
t255 = -t308 + t310;
t328 = t272 * t277;
t291 = -t249 * t255 + t328 * t333;
t350 = t269 * t291;
t216 = (qJD(3) * t290 + t297) * t269 + t352 * t270;
t276 = cos(qJ(5));
t274 = sin(qJ(5));
t331 = t288 * t274;
t229 = t248 * t276 - t331;
t221 = 0.1e1 / t229;
t222 = 0.1e1 / t229 ^ 2;
t349 = -0.2e1 * t242;
t348 = -0.2e1 * t289;
t214 = t289 * qJD(3) + t235 * t270 + t269 * t306;
t234 = -qJD(1) * t310 - t278 * t322 + (t273 * t303 + t302) * t275;
t204 = t229 * qJD(5) + t214 * t274 + t234 * t276;
t330 = t288 * t276;
t228 = t248 * t274 + t330;
t220 = t228 ^ 2;
t212 = t220 * t222 + 0.1e1;
t337 = t222 * t228;
t320 = qJD(5) * t228;
t205 = t214 * t276 - t234 * t274 - t320;
t342 = t205 * t221 * t222;
t344 = (t204 * t337 - t220 * t342) / t212 ^ 2;
t335 = t249 * t351;
t343 = (t215 * t333 - t240 * t335) / t230 ^ 2;
t341 = t207 * t213;
t339 = t218 * t289;
t338 = t219 * t289;
t336 = t228 * t276;
t334 = t242 * t249;
t332 = t288 * t269;
t326 = t274 * t221;
t323 = qJD(2) * t275;
t318 = -0.2e1 * t344;
t317 = 0.2e1 * t344;
t316 = -0.2e1 * t343;
t315 = t208 * t348;
t314 = t249 * t343;
t313 = t207 * t339;
t312 = t207 * t338;
t311 = t228 * t342;
t301 = 0.2e1 * t311;
t300 = t335 * t349;
t244 = -t269 * t327 - t270 * t290;
t296 = -qJD(5) * t270 * t288 + t235;
t227 = -t244 * t276 + t255 * t274;
t226 = -t244 * t274 - t255 * t276;
t293 = t222 * t336 - t326;
t292 = -t244 * t249 + t253 * t333;
t286 = -t218 + (t219 * t334 + t218) * t224;
t285 = -qJD(3) * t332 + qJD(5) * t258 + t234 * t270;
t239 = -t252 * qJD(3) + t270 * t305;
t236 = t288 * qJD(1) + t290 * qJD(2);
t233 = t258 * t274 + t270 * t330;
t232 = -t258 * t276 + t270 * t331;
t210 = 0.1e1 / t212;
t201 = 0.1e1 / t203;
t200 = t224 * t350;
t198 = t292 * t224;
t196 = t286 * t289;
t194 = (-t218 * t255 + t219 * t328) * t269 + t295 * t200;
t193 = t295 * t198 - t218 * t244 + t219 * t253;
t191 = t292 * t316 + (t253 * t300 - t216 * t249 + (t215 * t253 + t238 * t244 + t239 * t242) * t250) * t224;
t189 = t316 * t350 + (t291 * t321 + (t300 * t328 - t236 * t249 + (t238 * t255 + (t215 * t277 - t242 * t323) * t272) * t250) * t269) * t224;
t1 = [t314 * t348 + (-t213 * t249 - t289 * t351) * t224, t189, t191, 0, 0, 0; t242 * t206 * t319 + (-t215 * t206 + (t192 * t242 + t196 * t213) * t207) * t201 - (-t196 * t207 * t319 + (-0.2e1 * t196 * t345 + (-t197 * t224 * t334 + t316) * t313 + (t314 * t349 - t197 + (t197 - t294) * t224) * t312 - t286 * t341) * t201) * t289 (-t194 * t340 - t206 * t332) * t319 + (-t194 * t341 + (t234 * t269 + t288 * t321) * t206 + (t194 * t315 - t207 * t332) * t192 + (-t189 * t242 - t200 * t215 + (-t269 * t323 + t277 * t321) * t272 + (-t200 * t252 - t255 * t269) * t197) * t312 + (-t255 * t321 - t189 * t252 - t200 * t238 - t236 * t269 + (t200 * t242 - t269 * t328) * t197) * t313) * t201 (-t193 * t340 - t206 * t248) * t319 + (t193 * t192 * t315 + t214 * t206 + (-t248 * t192 - t193 * t213 + (-t191 * t242 - t198 * t215 + t239 + (-t198 * t252 - t244) * t197) * t338 + (-t191 * t252 - t198 * t238 - t216 + (t198 * t242 - t253) * t197) * t339) * t207) * t201, 0, 0, 0; (-t221 * t226 + t227 * t337) * t317 + ((t227 * qJD(5) - t216 * t274 - t236 * t276) * t221 + t227 * t301 + (-t226 * t205 - (-t226 * qJD(5) - t216 * t276 + t236 * t274) * t228 - t227 * t204) * t222) * t210 (-t221 * t232 + t233 * t337) * t317 + (t233 * t301 - t296 * t221 * t276 + t285 * t326 + (-t296 * t228 * t274 - t233 * t204 - t232 * t205 - t285 * t336) * t222) * t210, -t293 * t289 * t318 + (t293 * t213 - ((-qJD(5) * t221 - 0.2e1 * t311) * t276 + (t204 * t276 + (t205 - t320) * t274) * t222) * t289) * t210, 0, t318 + 0.2e1 * (t204 * t222 * t210 + (-t210 * t342 - t222 * t344) * t228) * t228, 0;];
JaD_rot  = t1;
