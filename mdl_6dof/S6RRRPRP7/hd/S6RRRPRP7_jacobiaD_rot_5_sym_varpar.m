% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP7
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

function JaD_rot = S6RRRPRP7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:39
% EndTime: 2019-02-26 22:12:41
% DurationCPUTime: 1.53s
% Computational Cost: add. (8428->149), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->128)
t271 = cos(pkin(6));
t273 = sin(qJ(2));
t345 = sin(qJ(1));
t306 = t345 * t273;
t296 = t271 * t306;
t300 = t345 * qJD(2);
t275 = cos(qJ(2));
t276 = cos(qJ(1));
t322 = t276 * t275;
t270 = sin(pkin(6));
t326 = t270 * t276;
t350 = -qJD(1) * t296 - t273 * t300 + (qJD(2) * t271 + qJD(1)) * t322 - qJD(3) * t326;
t269 = qJ(3) + pkin(11);
t267 = sin(t269);
t268 = cos(t269);
t305 = t345 * t275;
t323 = t276 * t273;
t288 = -t271 * t323 - t305;
t240 = -t267 * t288 + t268 * t326;
t328 = t270 * t273;
t250 = t267 * t328 - t271 * t268;
t229 = atan2(-t240, t250);
t216 = sin(t229);
t217 = cos(t229);
t207 = -t216 * t240 + t217 * t250;
t205 = 0.1e1 / t207 ^ 2;
t256 = -t296 + t322;
t307 = t270 * t345;
t287 = -t256 * t267 + t268 * t307;
t239 = t287 ^ 2;
t201 = t239 * t205 + 0.1e1;
t286 = -t271 * t305 - t323;
t233 = t288 * qJD(1) + t286 * qJD(2);
t246 = t256 * t268 + t267 * t307;
t304 = qJD(1) * t326;
t211 = t246 * qJD(3) + t233 * t267 - t268 * t304;
t338 = t211 * t205;
t238 = t240 ^ 2;
t248 = 0.1e1 / t250 ^ 2;
t228 = t238 * t248 + 0.1e1;
t222 = 0.1e1 / t228;
t301 = t345 * qJD(1);
t295 = t270 * t301;
t319 = qJD(3) * t268;
t213 = t350 * t267 - t268 * t295 - t288 * t319;
t251 = t271 * t267 + t268 * t328;
t320 = qJD(2) * t275;
t303 = t270 * t320;
t236 = t251 * qJD(3) + t267 * t303;
t247 = 0.1e1 / t250;
t332 = t240 * t248;
t292 = -t213 * t247 + t236 * t332;
t195 = t292 * t222;
t293 = -t216 * t250 - t217 * t240;
t190 = t293 * t195 - t216 * t213 + t217 * t236;
t204 = 0.1e1 / t207;
t206 = t204 * t205;
t343 = t190 * t206;
t317 = 0.2e1 * (-t239 * t343 - t287 * t338) / t201 ^ 2;
t349 = t236 * t248;
t308 = t271 * t322;
t253 = -t306 + t308;
t327 = t270 * t275;
t289 = -t247 * t253 + t327 * t332;
t348 = t267 * t289;
t214 = (qJD(3) * t288 + t295) * t267 + t350 * t268;
t274 = cos(qJ(5));
t272 = sin(qJ(5));
t330 = t286 * t272;
t227 = t246 * t274 - t330;
t219 = 0.1e1 / t227;
t220 = 0.1e1 / t227 ^ 2;
t347 = -0.2e1 * t240;
t346 = -0.2e1 * t287;
t212 = t287 * qJD(3) + t233 * t268 + t267 * t304;
t232 = -qJD(1) * t308 - t276 * t320 + (t271 * t300 + t301) * t273;
t202 = t227 * qJD(5) + t212 * t272 + t232 * t274;
t329 = t286 * t274;
t226 = t246 * t272 + t329;
t218 = t226 ^ 2;
t210 = t218 * t220 + 0.1e1;
t335 = t220 * t226;
t318 = qJD(5) * t226;
t203 = t212 * t274 - t232 * t272 - t318;
t340 = t203 * t219 * t220;
t342 = (t202 * t335 - t218 * t340) / t210 ^ 2;
t334 = t247 * t349;
t341 = (t213 * t332 - t238 * t334) / t228 ^ 2;
t339 = t205 * t287;
t337 = t216 * t287;
t336 = t217 * t287;
t333 = t240 * t247;
t331 = t286 * t267;
t325 = t272 * t219;
t324 = t274 * t226;
t321 = qJD(2) * t273;
t316 = -0.2e1 * t342;
t315 = 0.2e1 * t342;
t314 = -0.2e1 * t341;
t313 = t206 * t346;
t312 = t247 * t341;
t311 = t205 * t337;
t310 = t205 * t336;
t309 = t226 * t340;
t299 = 0.2e1 * t309;
t298 = t334 * t347;
t242 = -t267 * t326 - t268 * t288;
t294 = -qJD(5) * t268 * t286 + t233;
t225 = -t242 * t274 + t253 * t272;
t224 = -t242 * t272 - t253 * t274;
t291 = t220 * t324 - t325;
t290 = -t242 * t247 + t251 * t332;
t284 = -t216 + (t217 * t333 + t216) * t222;
t283 = -qJD(3) * t331 + qJD(5) * t256 + t232 * t268;
t237 = -t250 * qJD(3) + t268 * t303;
t234 = t286 * qJD(1) + t288 * qJD(2);
t231 = t256 * t272 + t268 * t329;
t230 = -t256 * t274 + t268 * t330;
t208 = 0.1e1 / t210;
t199 = 0.1e1 / t201;
t198 = t222 * t348;
t196 = t290 * t222;
t194 = t284 * t287;
t192 = (-t216 * t253 + t217 * t327) * t267 + t293 * t198;
t191 = t293 * t196 - t216 * t242 + t217 * t251;
t189 = t290 * t314 + (t251 * t298 - t214 * t247 + (t213 * t251 + t236 * t242 + t237 * t240) * t248) * t222;
t187 = t314 * t348 + (t289 * t319 + (t298 * t327 - t234 * t247 + (t236 * t253 + (t213 * t275 - t240 * t321) * t270) * t248) * t267) * t222;
t1 = [t312 * t346 + (-t211 * t247 - t287 * t349) * t222, t187, t189, 0, 0, 0; t240 * t204 * t317 + (-t213 * t204 + (t190 * t240 + t194 * t211) * t205) * t199 - (-t194 * t205 * t317 + (-0.2e1 * t194 * t343 + (-t195 * t222 * t333 + t314) * t311 + (t312 * t347 - t195 + (t195 - t292) * t222) * t310 - t284 * t338) * t199) * t287 (-t192 * t339 - t204 * t331) * t317 + (-t192 * t338 + (t232 * t267 + t286 * t319) * t204 + (t192 * t313 - t205 * t331) * t190 + (-t187 * t240 - t198 * t213 + (-t267 * t321 + t275 * t319) * t270 + (-t198 * t250 - t253 * t267) * t195) * t310 + (-t253 * t319 - t187 * t250 - t198 * t236 - t234 * t267 + (t198 * t240 - t267 * t327) * t195) * t311) * t199 (-t191 * t339 - t204 * t246) * t317 + (t191 * t190 * t313 + t212 * t204 + (-t246 * t190 - t191 * t211 + (-t189 * t240 - t196 * t213 + t237 + (-t196 * t250 - t242) * t195) * t336 + (-t189 * t250 - t196 * t236 - t214 + (t196 * t240 - t251) * t195) * t337) * t205) * t199, 0, 0, 0; (-t219 * t224 + t225 * t335) * t315 + ((t225 * qJD(5) - t214 * t272 - t234 * t274) * t219 + t225 * t299 + (-t224 * t203 - (-t224 * qJD(5) - t214 * t274 + t234 * t272) * t226 - t225 * t202) * t220) * t208 (-t219 * t230 + t231 * t335) * t315 + (t231 * t299 - t294 * t219 * t274 + t283 * t325 + (-t294 * t226 * t272 - t231 * t202 - t230 * t203 - t283 * t324) * t220) * t208, -t291 * t287 * t316 + (t291 * t211 - ((-qJD(5) * t219 - 0.2e1 * t309) * t274 + (t202 * t274 + (t203 - t318) * t272) * t220) * t287) * t208, 0, t316 + 0.2e1 * (t202 * t220 * t208 + (-t208 * t340 - t220 * t342) * t226) * t226, 0;];
JaD_rot  = t1;
