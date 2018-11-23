% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR10_jacobiaD_rot_3_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_rot_3_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_rot_3_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiaD_rot_3_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:17
% EndTime: 2018-11-23 11:27:19
% DurationCPUTime: 1.27s
% Computational Cost: add. (7877->134), mult. (10089->265), div. (449->12), fcn. (10309->21), ass. (0->127)
t350 = pkin(6) + qJ(2);
t339 = sin(t350);
t374 = -qJD(2) / 0.2e1;
t291 = t339 * t374;
t351 = pkin(6) - qJ(2);
t340 = sin(t351);
t331 = qJD(2) * t340;
t275 = t331 / 0.2e1 + t291;
t313 = cos(qJ(2));
t314 = cos(qJ(1));
t301 = cos(t350) / 0.2e1;
t304 = cos(t351);
t284 = t304 / 0.2e1 + t301;
t310 = sin(qJ(2));
t311 = sin(qJ(1));
t328 = t311 * t284 + t310 * t314;
t355 = qJD(2) * t311;
t248 = t328 * qJD(1) - t275 * t314 + t313 * t355;
t307 = sin(pkin(7));
t308 = sin(pkin(6));
t373 = cos(pkin(7));
t342 = t308 * t373;
t335 = qJD(1) * t342;
t235 = -t248 * t307 - t311 * t335;
t263 = -t284 * t314 + t311 * t310;
t323 = -t263 * t307 + t314 * t342;
t251 = t323 ^ 2;
t333 = t339 / 0.2e1;
t262 = -(t333 + t340 / 0.2e1) * t307 + cos(pkin(6)) * t373;
t260 = 0.1e1 / t262 ^ 2;
t240 = t251 * t260 + 0.1e1;
t259 = 0.1e1 / t262;
t261 = t259 * t260;
t277 = qJD(2) * t301 + t304 * t374;
t360 = t277 * t307;
t369 = (t235 * t260 * t323 + t251 * t261 * t360) / t240 ^ 2;
t379 = -0.2e1 * t369;
t280 = t333 - t340 / 0.2e1;
t268 = t311 * t280 - t313 * t314;
t349 = pkin(7) - qJ(3);
t338 = sin(t349);
t332 = t338 / 0.2e1;
t306 = pkin(7) + qJ(3);
t302 = sin(t306);
t376 = t302 / 0.2e1;
t278 = t376 + t332;
t341 = cos(t349);
t345 = cos(t306) / 0.2e1;
t282 = t345 + t341 / 0.2e1;
t309 = sin(qJ(3));
t359 = t308 * t311;
t231 = -t268 * t309 - t278 * t359 + t282 * t328;
t225 = t231 ^ 2;
t279 = t376 - t338 / 0.2e1;
t281 = t345 - t341 / 0.2e1;
t312 = cos(qJ(3));
t324 = -t268 * t312 - t279 * t328 - t281 * t359;
t227 = 0.1e1 / t324 ^ 2;
t378 = t225 * t227;
t329 = t280 * t314 + t311 * t313;
t283 = t301 - t304 / 0.2e1;
t362 = t323 * t283;
t377 = t307 * (-t259 * t329 + t260 * t362);
t276 = t284 * qJD(2);
t354 = qJD(2) * t314;
t246 = t329 * qJD(1) + t311 * t276 + t310 * t354;
t243 = atan2(t323, t262);
t236 = sin(t243);
t237 = cos(t243);
t221 = t236 * t323 + t237 * t262;
t218 = 0.1e1 / t221;
t226 = 0.1e1 / t324;
t219 = 0.1e1 / t221 ^ 2;
t255 = -t307 * t328 - t311 * t342;
t252 = t255 ^ 2;
t216 = t219 * t252 + 0.1e1;
t245 = t263 * qJD(1) - t311 * t275 - t313 * t354;
t234 = t245 * t307 - t314 * t335;
t366 = t234 * t219;
t238 = 0.1e1 / t240;
t343 = t260 * t360;
t325 = t235 * t259 + t323 * t343;
t210 = t325 * t238;
t330 = -t236 * t262 + t237 * t323;
t206 = t330 * t210 + t235 * t236 - t237 * t360;
t371 = t206 * t218 * t219;
t372 = (-t252 * t371 + t255 * t366) / t216 ^ 2;
t270 = t278 * qJD(3);
t272 = t282 * qJD(3);
t353 = qJD(3) * t309;
t356 = qJD(1) * t314;
t213 = t268 * t353 + t245 * t279 - t328 * t272 - t246 * t312 + (t270 * t311 - t281 * t356) * t308;
t228 = t226 * t227;
t370 = t213 * t228;
t368 = t219 * t255;
t367 = t227 * t231;
t365 = t236 * t255;
t364 = t237 * t255;
t363 = t323 * t259;
t361 = t277 * t307 ^ 2;
t358 = t308 * t314;
t357 = qJD(1) * t311;
t352 = qJD(3) * t312;
t348 = -0.2e1 * t372;
t347 = -0.2e1 * t371;
t224 = 0.1e1 + t378;
t271 = (t332 - t302 / 0.2e1) * qJD(3);
t273 = t281 * qJD(3);
t212 = -t268 * t352 - t245 * t282 - t246 * t309 + t328 * t271 + (-t273 * t311 - t278 * t356) * t308;
t344 = t212 * t367;
t346 = 0.2e1 * (-t225 * t370 + t344) / t224 ^ 2;
t337 = 0.2e1 * t231 * t370;
t336 = t259 * t379;
t322 = t236 + (t237 * t363 - t236) * t238;
t321 = t268 * qJD(1) - t276 * t314 + t310 * t355;
t274 = t291 - t331 / 0.2e1;
t242 = t268 * t279 - t312 * t328;
t241 = -t268 * t282 - t309 * t328;
t230 = t263 * t279 - t281 * t358 - t312 * t329;
t229 = -t263 * t282 - t278 * t358 - t309 * t329;
t222 = 0.1e1 / t224;
t214 = 0.1e1 / t216;
t211 = t238 * t377;
t209 = t322 * t255;
t207 = (-t236 * t329 - t237 * t283) * t307 + t330 * t211;
t205 = t377 * t379 + (0.2e1 * t261 * t361 * t362 + t321 * t259 * t307 + (-t329 * t361 + (t235 * t283 + t274 * t323) * t307) * t260) * t238;
t1 = [t255 * t336 + (t234 * t259 + t255 * t343) * t238, t205, 0, 0, 0, 0; t323 * t218 * t348 + (t235 * t218 + (-t206 * t323 + t209 * t234) * t219) * t214 + ((t209 * t347 + t322 * t366) * t214 + (t209 * t348 + ((-t210 * t238 * t363 + 0.2e1 * t369) * t365 + (t323 * t336 + t210 + (-t210 + t325) * t238) * t364) * t214) * t219) * t255, 0.2e1 * (t218 * t268 * t307 - t207 * t368) * t372 + ((t330 * t205 + (-t221 * t210 + t235 * t237) * t211) * t368 + (t255 * t347 + t366) * t207 + (-t246 * t218 + (t268 * t206 + (-t210 * t329 - t274) * t364 + (t210 * t283 + t211 * t277 + t321) * t365) * t219) * t307) * t214, 0, 0, 0, 0; (-t226 * t229 + t230 * t367) * t346 + ((-t329 * t352 - t248 * t282 + t321 * t309 - t263 * t271 + (-t273 * t314 + t278 * t357) * t308) * t226 + t230 * t337 + (-t229 * t213 - (t329 * t353 + t248 * t279 + t321 * t312 + t263 * t272 + (t270 * t314 + t281 * t357) * t308) * t231 - t230 * t212) * t227) * t222 (-t226 * t241 + t242 * t367) * t346 + ((t245 * t309 - t246 * t282 - t268 * t271 - t328 * t352) * t226 + t242 * t337 + (-t241 * t213 - (t245 * t312 + t246 * t279 + t268 * t272 + t328 * t353) * t231 - t242 * t212) * t227) * t222 (-t226 * t324 - t378) * t346 + (0.2e1 * t344 + (-0.2e1 * t225 * t228 - t227 * t324 + t226) * t213) * t222, 0, 0, 0;];
JaD_rot  = t1;
