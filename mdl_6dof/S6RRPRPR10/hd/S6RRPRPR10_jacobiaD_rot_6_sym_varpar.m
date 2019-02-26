% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:01
% EndTime: 2019-02-26 21:43:02
% DurationCPUTime: 1.42s
% Computational Cost: add. (8428->151), mult. (13478->298), div. (726->12), fcn. (17045->13), ass. (0->129)
t270 = sin(qJ(1));
t273 = cos(qJ(1));
t267 = cos(pkin(6));
t290 = qJD(2) * t267 + qJD(1);
t269 = sin(qJ(2));
t317 = t270 * t269;
t299 = t267 * t317;
t266 = sin(pkin(6));
t307 = qJD(4) * t266;
t311 = qJD(2) * t269;
t272 = cos(qJ(2));
t313 = t273 * t272;
t344 = -qJD(1) * t299 - t270 * t311 - t273 * t307 + t290 * t313;
t253 = -t299 + t313;
t265 = pkin(11) + qJ(4);
t263 = sin(t265);
t264 = cos(t265);
t321 = t266 * t270;
t241 = t253 * t263 - t264 * t321;
t271 = cos(qJ(6));
t314 = t273 * t269;
t316 = t270 * t272;
t252 = t267 * t316 + t314;
t268 = sin(qJ(6));
t324 = t252 * t268;
t286 = t241 * t271 - t324;
t343 = t286 * qJD(6);
t251 = t267 * t314 + t316;
t319 = t266 * t273;
t237 = t251 * t264 - t263 * t319;
t322 = t266 * t269;
t249 = t267 * t263 + t264 * t322;
t224 = atan2(-t237, t249);
t211 = sin(t224);
t212 = cos(t224);
t202 = -t211 * t237 + t212 * t249;
t200 = 0.1e1 / t202 ^ 2;
t242 = t253 * t264 + t263 * t321;
t235 = t242 ^ 2;
t196 = t235 * t200 + 0.1e1;
t229 = -t251 * qJD(1) - t252 * qJD(2);
t312 = qJD(1) * t266;
t296 = t273 * t312;
t309 = qJD(4) * t263;
t207 = t263 * t296 - t253 * t309 + (t270 * t307 + t229) * t264;
t332 = t207 * t200;
t234 = t237 ^ 2;
t246 = 0.1e1 / t249 ^ 2;
t223 = t234 * t246 + 0.1e1;
t217 = 0.1e1 / t223;
t291 = t344 * t264;
t297 = t270 * t312;
t209 = -t251 * t309 + t263 * t297 + t291;
t248 = -t263 * t322 + t267 * t264;
t310 = qJD(2) * t272;
t295 = t266 * t310;
t233 = t248 * qJD(4) + t264 * t295;
t245 = 0.1e1 / t249;
t326 = t237 * t246;
t285 = -t209 * t245 + t233 * t326;
t190 = t285 * t217;
t288 = -t211 * t249 - t212 * t237;
t185 = t288 * t190 - t211 * t209 + t212 * t233;
t199 = 0.1e1 / t202;
t201 = t199 * t200;
t337 = t185 * t201;
t306 = 0.2e1 * (-t235 * t337 + t242 * t332) / t196 ^ 2;
t342 = t233 * t246;
t298 = t267 * t313;
t250 = t298 - t317;
t320 = t266 * t272;
t282 = -t245 * t250 + t320 * t326;
t341 = t264 * t282;
t323 = t252 * t271;
t222 = t241 * t268 + t323;
t214 = 0.1e1 / t222;
t215 = 0.1e1 / t222 ^ 2;
t340 = -0.2e1 * t237;
t339 = 0.2e1 * t242;
t206 = t242 * qJD(4) + t229 * t263 - t264 * t296;
t228 = -qJD(1) * t298 - t273 * t310 + t290 * t317;
t197 = t222 * qJD(6) - t206 * t271 - t228 * t268;
t213 = t286 ^ 2;
t205 = t213 * t215 + 0.1e1;
t329 = t215 * t286;
t198 = t206 * t268 - t228 * t271 + t343;
t334 = t198 * t214 * t215;
t336 = (-t197 * t329 - t213 * t334) / t205 ^ 2;
t328 = t245 * t342;
t335 = (t209 * t326 - t234 * t328) / t223 ^ 2;
t333 = t200 * t242;
t331 = t211 * t242;
t330 = t212 * t242;
t327 = t237 * t245;
t325 = t252 * t264;
t318 = t268 * t286;
t315 = t271 * t214;
t308 = qJD(4) * t264;
t305 = 0.2e1 * t336;
t304 = -0.2e1 * t335;
t303 = t201 * t339;
t302 = t245 * t335;
t301 = t200 * t331;
t300 = t200 * t330;
t293 = -0.2e1 * t286 * t334;
t292 = t328 * t340;
t289 = -qJD(6) * t252 * t263 + t229;
t236 = t251 * t263 + t264 * t319;
t287 = -t236 * t271 - t250 * t268;
t220 = -t236 * t268 + t250 * t271;
t284 = -t215 * t318 + t315;
t283 = t236 * t245 + t248 * t326;
t281 = -t211 + (t212 * t327 + t211) * t217;
t208 = t251 * t308 + t263 * t344 - t264 * t297;
t280 = qJD(6) * t253 - t228 * t263 + t252 * t308;
t232 = -t249 * qJD(4) - t263 * t295;
t230 = -t252 * qJD(1) - t251 * qJD(2);
t226 = t253 * t271 - t263 * t324;
t225 = t253 * t268 + t263 * t323;
t203 = 0.1e1 / t205;
t194 = 0.1e1 / t196;
t193 = t217 * t341;
t191 = t283 * t217;
t189 = t281 * t242;
t187 = (-t211 * t250 + t212 * t320) * t264 + t288 * t193;
t186 = t288 * t191 + t211 * t236 + t212 * t248;
t184 = t283 * t304 + (t248 * t292 + t208 * t245 + (t209 * t248 + t232 * t237 - t233 * t236) * t246) * t217;
t182 = t304 * t341 + (-t282 * t309 + (t292 * t320 - t230 * t245 + (t233 * t250 + (t209 * t272 - t237 * t311) * t266) * t246) * t264) * t217;
t1 = [t302 * t339 + (-t207 * t245 + t242 * t342) * t217, t182, 0, t184, 0, 0; t237 * t199 * t306 + (((qJD(4) * t251 - t297) * t263 - t291) * t199 + (t185 * t237 - t189 * t207) * t200) * t194 + (t189 * t200 * t306 + (0.2e1 * t189 * t337 - (-t190 * t217 * t327 + t304) * t301 - (t302 * t340 - t190 + (t190 - t285) * t217) * t300 - t281 * t332) * t194) * t242 (t187 * t333 + t199 * t325) * t306 + (-t187 * t332 + (t228 * t264 + t252 * t309) * t199 + (t187 * t303 + t200 * t325) * t185 - (-t182 * t237 - t193 * t209 + (-t264 * t311 - t272 * t309) * t266 + (-t193 * t249 - t250 * t264) * t190) * t300 - (t250 * t309 - t182 * t249 - t193 * t233 - t230 * t264 + (t193 * t237 - t264 * t320) * t190) * t301) * t194, 0 (t186 * t333 + t199 * t241) * t306 + (t186 * t185 * t303 - t206 * t199 + (t241 * t185 - t186 * t207 - (-t184 * t237 - t191 * t209 + t232 + (-t191 * t249 + t236) * t190) * t330 - (-t184 * t249 - t191 * t233 + t208 + (t191 * t237 - t248) * t190) * t331) * t200) * t194, 0, 0; (t214 * t287 - t220 * t329) * t305 + ((t220 * qJD(6) + t208 * t271 + t230 * t268) * t214 + t220 * t293 + (t287 * t198 + (t287 * qJD(6) - t208 * t268 + t230 * t271) * t286 - t220 * t197) * t215) * t203 (-t214 * t225 - t226 * t329) * t305 + (t226 * t293 + t289 * t214 * t268 + t280 * t315 + (t271 * t286 * t289 - t226 * t197 - t225 * t198 - t280 * t318) * t215) * t203, 0, t284 * t242 * t305 + (-t284 * t207 + ((qJD(6) * t214 + t293) * t268 + (-t197 * t268 + (t198 + t343) * t271) * t215) * t242) * t203, 0, -0.2e1 * t336 - 0.2e1 * (t197 * t215 * t203 - (-t203 * t334 - t215 * t336) * t286) * t286;];
JaD_rot  = t1;
