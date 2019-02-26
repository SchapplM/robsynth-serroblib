% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:51
% EndTime: 2019-02-26 22:05:52
% DurationCPUTime: 1.42s
% Computational Cost: add. (8849->150), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->128)
t276 = cos(pkin(6));
t277 = sin(qJ(2));
t347 = sin(qJ(1));
t309 = t347 * t277;
t299 = t276 * t309;
t303 = t347 * qJD(2);
t278 = cos(qJ(2));
t279 = cos(qJ(1));
t325 = t279 * t278;
t275 = sin(pkin(6));
t327 = t275 * t279;
t352 = -qJD(1) * t299 - t277 * t303 + (qJD(2) * t276 + qJD(1)) * t325 - qJD(3) * t327;
t274 = qJ(3) + pkin(11);
t270 = sin(t274);
t272 = cos(t274);
t308 = t347 * t278;
t326 = t279 * t277;
t291 = -t276 * t326 - t308;
t242 = -t270 * t291 + t272 * t327;
t329 = t275 * t277;
t253 = t270 * t329 - t276 * t272;
t231 = atan2(-t242, t253);
t226 = sin(t231);
t227 = cos(t231);
t209 = -t226 * t242 + t227 * t253;
t207 = 0.1e1 / t209 ^ 2;
t258 = -t299 + t325;
t310 = t275 * t347;
t290 = -t258 * t270 + t272 * t310;
t241 = t290 ^ 2;
t203 = t241 * t207 + 0.1e1;
t289 = -t276 * t308 - t326;
t235 = t291 * qJD(1) + t289 * qJD(2);
t248 = t258 * t272 + t270 * t310;
t307 = qJD(1) * t327;
t213 = t248 * qJD(3) + t235 * t270 - t272 * t307;
t340 = t213 * t207;
t240 = t242 ^ 2;
t251 = 0.1e1 / t253 ^ 2;
t230 = t240 * t251 + 0.1e1;
t228 = 0.1e1 / t230;
t304 = t347 * qJD(1);
t298 = t275 * t304;
t322 = qJD(3) * t272;
t215 = t352 * t270 - t272 * t298 - t291 * t322;
t254 = t276 * t270 + t272 * t329;
t323 = qJD(2) * t278;
t306 = t275 * t323;
t238 = t254 * qJD(3) + t270 * t306;
t250 = 0.1e1 / t253;
t334 = t242 * t251;
t295 = -t215 * t250 + t238 * t334;
t197 = t295 * t228;
t296 = -t226 * t253 - t227 * t242;
t192 = t296 * t197 - t226 * t215 + t227 * t238;
t206 = 0.1e1 / t209;
t208 = t206 * t207;
t345 = t192 * t208;
t320 = 0.2e1 * (-t241 * t345 - t290 * t340) / t203 ^ 2;
t351 = t238 * t251;
t311 = t276 * t325;
t255 = -t309 + t311;
t328 = t275 * t278;
t292 = -t250 * t255 + t328 * t334;
t350 = t270 * t292;
t216 = (qJD(3) * t291 + t298) * t270 + t352 * t272;
t273 = pkin(12) + qJ(6);
t269 = sin(t273);
t271 = cos(t273);
t225 = t248 * t271 - t269 * t289;
t219 = 0.1e1 / t225;
t220 = 0.1e1 / t225 ^ 2;
t349 = -0.2e1 * t242;
t348 = -0.2e1 * t290;
t214 = t290 * qJD(3) + t235 * t272 + t270 * t307;
t234 = -qJD(1) * t311 - t279 * t323 + (t276 * t303 + t304) * t277;
t204 = t225 * qJD(6) + t214 * t269 + t234 * t271;
t224 = t248 * t269 + t271 * t289;
t218 = t224 ^ 2;
t212 = t218 * t220 + 0.1e1;
t339 = t220 * t224;
t321 = qJD(6) * t224;
t205 = t214 * t271 - t234 * t269 - t321;
t342 = t205 * t219 * t220;
t344 = (t204 * t339 - t218 * t342) / t212 ^ 2;
t336 = t250 * t351;
t343 = (t215 * t334 - t240 * t336) / t230 ^ 2;
t341 = t207 * t290;
t338 = t226 * t290;
t337 = t227 * t290;
t335 = t242 * t250;
t333 = t289 * t270;
t332 = t289 * t272;
t331 = t269 * t219;
t330 = t271 * t224;
t324 = qJD(2) * t277;
t319 = -0.2e1 * t344;
t318 = 0.2e1 * t344;
t317 = -0.2e1 * t343;
t316 = t208 * t348;
t315 = t250 * t343;
t314 = t207 * t338;
t313 = t207 * t337;
t312 = t224 * t342;
t302 = 0.2e1 * t312;
t301 = t336 * t349;
t244 = -t270 * t327 - t272 * t291;
t297 = -qJD(6) * t332 + t235;
t223 = -t244 * t271 + t255 * t269;
t222 = -t244 * t269 - t255 * t271;
t294 = t220 * t330 - t331;
t293 = -t244 * t250 + t254 * t334;
t287 = -t226 + (t227 * t335 + t226) * t228;
t286 = -qJD(3) * t333 + qJD(6) * t258 + t234 * t272;
t239 = -t253 * qJD(3) + t272 * t306;
t236 = t289 * qJD(1) + t291 * qJD(2);
t233 = t258 * t269 + t271 * t332;
t232 = -t258 * t271 + t269 * t332;
t210 = 0.1e1 / t212;
t201 = 0.1e1 / t203;
t200 = t228 * t350;
t198 = t293 * t228;
t196 = t287 * t290;
t194 = (-t226 * t255 + t227 * t328) * t270 + t296 * t200;
t193 = t296 * t198 - t226 * t244 + t227 * t254;
t191 = t293 * t317 + (t254 * t301 - t216 * t250 + (t215 * t254 + t238 * t244 + t239 * t242) * t251) * t228;
t189 = t317 * t350 + (t292 * t322 + (t301 * t328 - t236 * t250 + (t238 * t255 + (t215 * t278 - t242 * t324) * t275) * t251) * t270) * t228;
t1 = [t315 * t348 + (-t213 * t250 - t290 * t351) * t228, t189, t191, 0, 0, 0; t242 * t206 * t320 + (-t215 * t206 + (t192 * t242 + t196 * t213) * t207) * t201 - (-t196 * t207 * t320 + (-0.2e1 * t196 * t345 + (-t197 * t228 * t335 + t317) * t314 + (t315 * t349 - t197 + (t197 - t295) * t228) * t313 - t287 * t340) * t201) * t290 (-t194 * t341 - t206 * t333) * t320 + (-t194 * t340 + (t234 * t270 + t289 * t322) * t206 + (t194 * t316 - t207 * t333) * t192 + (-t189 * t242 - t200 * t215 + (-t270 * t324 + t278 * t322) * t275 + (-t200 * t253 - t255 * t270) * t197) * t313 + (-t255 * t322 - t189 * t253 - t200 * t238 - t236 * t270 + (t200 * t242 - t270 * t328) * t197) * t314) * t201 (-t193 * t341 - t206 * t248) * t320 + (t193 * t192 * t316 + t214 * t206 + (-t248 * t192 - t193 * t213 + (-t191 * t242 - t198 * t215 + t239 + (-t198 * t253 - t244) * t197) * t337 + (-t191 * t253 - t198 * t238 - t216 + (t198 * t242 - t254) * t197) * t338) * t207) * t201, 0, 0, 0; (-t219 * t222 + t223 * t339) * t318 + ((t223 * qJD(6) - t216 * t269 - t236 * t271) * t219 + t223 * t302 + (-t222 * t205 - (-t222 * qJD(6) - t216 * t271 + t236 * t269) * t224 - t223 * t204) * t220) * t210 (-t219 * t232 + t233 * t339) * t318 + (t233 * t302 - t297 * t219 * t271 + t286 * t331 + (-t297 * t224 * t269 - t233 * t204 - t232 * t205 - t286 * t330) * t220) * t210, -t294 * t290 * t319 + (t294 * t213 - ((-qJD(6) * t219 - 0.2e1 * t312) * t271 + (t204 * t271 + (t205 - t321) * t269) * t220) * t290) * t210, 0, 0, t319 + 0.2e1 * (t204 * t220 * t210 + (-t210 * t342 - t220 * t344) * t224) * t224;];
JaD_rot  = t1;
