% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPPRRR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiaD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:44
% EndTime: 2019-02-26 19:38:45
% DurationCPUTime: 1.13s
% Computational Cost: add. (7555->76), mult. (22578->165), div. (281->12), fcn. (29955->19), ass. (0->93)
t289 = sin(pkin(14));
t339 = sin(pkin(13));
t340 = sin(pkin(12));
t318 = t340 * t339;
t345 = cos(pkin(13));
t346 = cos(pkin(12));
t324 = t346 * t345;
t349 = cos(pkin(6));
t308 = t349 * t324 - t318;
t320 = t340 * t345;
t323 = t346 * t339;
t309 = t349 * t323 + t320;
t343 = sin(pkin(6));
t348 = cos(pkin(7));
t325 = t348 * t343;
t341 = sin(pkin(8));
t342 = sin(pkin(7));
t344 = cos(pkin(14));
t347 = cos(pkin(8));
t322 = t343 * t342;
t354 = t308 * t348 - t346 * t322;
t356 = (-t309 * t289 + t354 * t344) * t347 + (-t308 * t342 - t346 * t325) * t341;
t311 = -t349 * t318 + t324;
t310 = -t349 * t320 - t323;
t319 = t340 * t343;
t353 = t310 * t348 + t342 * t319;
t301 = t311 * t289 - t353 * t344;
t304 = -t310 * t342 + t348 * t319;
t355 = t301 * t347 - t304 * t341;
t352 = t345 * t325 + t349 * t342;
t321 = t343 * t339;
t351 = t347 * (-t289 * t321 + t352 * t344) + t341 * (-t345 * t322 + t349 * t348);
t282 = t354 * t289 + t309 * t344;
t291 = sin(qJ(4));
t350 = cos(qJ(4));
t266 = t282 * t291 - t356 * t350;
t287 = t352 * t289 + t344 * t321;
t307 = -t287 * t291 + t351 * t350;
t259 = atan2(-t266, -t307);
t254 = sin(t259);
t255 = cos(t259);
t242 = -t254 * t266 - t255 * t307;
t239 = 0.1e1 / t242;
t283 = t353 * t289 + t311 * t344;
t270 = t283 * t350 - t355 * t291;
t279 = t301 * t341 + t304 * t347;
t290 = sin(qJ(5));
t292 = cos(qJ(5));
t253 = t270 * t292 + t279 * t290;
t249 = 0.1e1 / t253;
t274 = 0.1e1 / t307;
t240 = 0.1e1 / t242 ^ 2;
t250 = 0.1e1 / t253 ^ 2;
t275 = 0.1e1 / t307 ^ 2;
t264 = t266 ^ 2;
t258 = t264 * t275 + 0.1e1;
t256 = 0.1e1 / t258;
t268 = t282 * t350 + t356 * t291;
t261 = t268 * qJD(4);
t278 = t287 * t350 + t351 * t291;
t272 = t278 * qJD(4);
t333 = t266 * t275;
t233 = (t261 * t274 + t272 * t333) * t256;
t317 = t254 * t307 - t255 * t266;
t230 = t317 * t233 - t254 * t261 + t255 * t272;
t338 = t230 * t239 * t240;
t269 = t283 * t291 + t355 * t350;
t337 = t240 * t269;
t252 = t270 * t290 - t279 * t292;
t248 = t252 ^ 2;
t245 = t248 * t250 + 0.1e1;
t262 = t269 * qJD(4);
t246 = t253 * qJD(5) - t262 * t290;
t334 = t250 * t252;
t330 = qJD(5) * t252;
t247 = -t262 * t292 - t330;
t335 = t247 * t249 * t250;
t336 = 0.1e1 / t245 ^ 2 * (t246 * t334 - t248 * t335);
t332 = t266 * t278;
t331 = t272 * t274 * t275;
t329 = -0.2e1 * t336;
t316 = -t249 * t290 + t292 * t334;
t315 = t268 * t274 + t275 * t332;
t271 = t307 * qJD(4);
t265 = t269 ^ 2;
t263 = t270 * qJD(4);
t260 = t266 * qJD(4);
t243 = 0.1e1 / t245;
t237 = t240 * t265 + 0.1e1;
t234 = t315 * t256;
t231 = t317 * t234 - t254 * t268 + t255 * t278;
t229 = -0.2e1 * t315 / t258 ^ 2 * (t261 * t333 + t264 * t331) + (0.2e1 * t331 * t332 - t260 * t274 + (t261 * t278 + t266 * t271 + t268 * t272) * t275) * t256;
t1 = [0, 0, 0, t229, 0, 0; 0, 0, 0, 0.2e1 * (t231 * t337 - t239 * t270) / t237 ^ 2 * (t263 * t337 - t265 * t338) + (-t262 * t239 + (-t270 * t230 - t231 * t263) * t240 + (0.2e1 * t231 * t338 + (-(-t229 * t266 - t234 * t261 + t271 + (t234 * t307 - t268) * t233) * t255 - (t229 * t307 - t234 * t272 + t260 + (t234 * t266 - t278) * t233) * t254) * t240) * t269) / t237, 0, 0; 0, 0, 0, t316 * t269 * t329 + (t316 * t263 + ((-qJD(5) * t249 - 0.2e1 * t252 * t335) * t292 + (t246 * t292 + (t247 - t330) * t290) * t250) * t269) * t243, t329 + 0.2e1 * (t243 * t246 * t250 + (-t243 * t335 - t250 * t336) * t252) * t252, 0;];
JaD_rot  = t1;
