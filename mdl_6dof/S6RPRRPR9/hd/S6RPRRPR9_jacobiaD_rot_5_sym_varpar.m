% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:29
% EndTime: 2019-02-26 21:05:30
% DurationCPUTime: 1.44s
% Computational Cost: add. (5146->112), mult. (15324->225), div. (448->12), fcn. (19446->15), ass. (0->109)
t358 = sin(pkin(12));
t363 = cos(pkin(6));
t332 = t363 * t358;
t361 = cos(pkin(12));
t364 = sin(qJ(1));
t366 = cos(qJ(1));
t287 = t366 * t332 + t364 * t361;
t300 = sin(qJ(3));
t359 = sin(pkin(7));
t360 = sin(pkin(6));
t330 = t360 * t359;
t323 = t366 * t330;
t294 = t300 * t323;
t334 = t363 * t361;
t286 = -t366 * t334 + t364 * t358;
t362 = cos(pkin(7));
t338 = t286 * t362;
t365 = cos(qJ(3));
t372 = -t287 * t365 + t300 * t338 + t294;
t316 = t365 * t323;
t336 = t362 * t365;
t371 = -t286 * t336 - t316;
t329 = t360 * t358;
t331 = t362 * t360;
t370 = t361 * t331 + t363 * t359;
t276 = t370 * t300 + t365 * t329;
t268 = t276 * qJD(3);
t307 = -t300 * t329 + t370 * t365;
t273 = 0.1e1 / t307 ^ 2;
t347 = t268 * t273;
t288 = -t364 * t332 + t366 * t361;
t312 = t364 * t334 + t366 * t358;
t310 = t312 * t362;
t320 = t364 * t330;
t315 = t365 * t320;
t369 = -t288 * t300 - t365 * t310 + t315;
t284 = t312 * qJD(1);
t285 = t288 * qJD(1);
t242 = (qJD(1) * t320 - qJD(3) * t287 - t362 * t284) * t300 + t285 * t365 + t371 * qJD(3);
t368 = qJD(1) * t315 + t372 * qJD(3) - t284 * t336 - t285 * t300;
t260 = t287 * t300 - t371;
t257 = atan2(-t260, -t307);
t252 = sin(t257);
t253 = cos(t257);
t233 = -t252 * t260 - t253 * t307;
t230 = 0.1e1 / t233;
t266 = t288 * t365 + (-t310 + t320) * t300;
t321 = t364 * t331;
t278 = t312 * t359 + t321;
t299 = qJ(4) + pkin(13);
t297 = sin(t299);
t298 = cos(t299);
t251 = t266 * t298 + t278 * t297;
t245 = 0.1e1 / t251;
t272 = 0.1e1 / t307;
t231 = 0.1e1 / t233 ^ 2;
t246 = 0.1e1 / t251 ^ 2;
t367 = -0.2e1 * t369;
t259 = t369 ^ 2;
t229 = t231 * t259 + 0.1e1;
t283 = t287 * qJD(1);
t309 = qJD(1) * t338;
t239 = -qJD(1) * t316 + t266 * qJD(3) - t283 * t300 - t365 * t309;
t353 = t231 * t369;
t258 = t260 ^ 2;
t256 = t258 * t273 + 0.1e1;
t254 = 0.1e1 / t256;
t326 = t260 * t347 - t272 * t368;
t224 = t326 * t254;
t328 = t252 * t307 - t253 * t260;
t220 = t328 * t224 + t252 * t368 + t253 * t268;
t356 = t220 * t230 * t231;
t357 = (-t239 * t353 - t259 * t356) / t229 ^ 2;
t346 = t272 * t347;
t355 = (-t260 * t273 * t368 + t258 * t346) / t256 ^ 2;
t227 = 0.1e1 / t229;
t354 = t227 * t231;
t240 = qJD(1) * t294 + t369 * qJD(3) - t283 * t365 + t300 * t309;
t277 = -t286 * t359 + t366 * t331;
t270 = t277 * qJD(1);
t250 = t266 * t297 - t278 * t298;
t344 = qJD(4) * t250;
t235 = t240 * t298 + t270 * t297 - t344;
t352 = t235 * t245 * t246;
t234 = t251 * qJD(4) + t240 * t297 - t270 * t298;
t244 = t250 ^ 2;
t238 = t244 * t246 + 0.1e1;
t350 = t246 * t250;
t351 = 0.1e1 / t238 ^ 2 * (t234 * t350 - t244 * t352);
t349 = t260 * t272;
t348 = t260 * t276;
t343 = 0.2e1 * t357;
t342 = -0.2e1 * t355;
t341 = -0.2e1 * t351;
t340 = t272 * t355;
t339 = t250 * t352;
t337 = t356 * t367;
t249 = t277 * t297 + t298 * t372;
t248 = -t277 * t298 + t297 * t372;
t325 = -t245 * t297 + t298 * t350;
t324 = -t272 * t372 + t273 * t348;
t318 = -t252 + (-t253 * t349 + t252) * t254;
t271 = -qJD(1) * t321 - t284 * t359;
t267 = t307 * qJD(3);
t236 = 0.1e1 / t238;
t225 = t324 * t254;
t221 = t328 * t225 + t252 * t372 + t253 * t276;
t219 = t324 * t342 + (0.2e1 * t346 * t348 + t242 * t272 + (t260 * t267 - t268 * t372 - t276 * t368) * t273) * t254;
t1 = [-t340 * t367 + (t239 * t272 - t347 * t369) * t254, 0, t219, 0, 0, 0; -(-t220 * t354 - 0.2e1 * t230 * t357) * t260 + (t368 * t230 + (t318 * t239 - ((t224 * t254 * t349 + t342) * t252 + (0.2e1 * t260 * t340 - t224 + (t224 - t326) * t254) * t253) * t369) * t353) * t227 - (t227 * t337 - t239 * t354 - t353 * t343) * t318 * t369, 0 (-t221 * t353 - t230 * t266) * t343 + (t221 * t337 + t240 * t230 + (-t266 * t220 - t221 * t239 - (-(-t219 * t260 + t225 * t368 + t267 + (t225 * t307 + t372) * t224) * t253 - (t219 * t307 - t225 * t268 - t242 + (t225 * t260 - t276) * t224) * t252) * t369) * t231) * t227, 0, 0, 0; 0.2e1 * (-t245 * t248 + t249 * t350) * t351 + ((t249 * qJD(4) - t242 * t297 - t271 * t298) * t245 + 0.2e1 * t249 * t339 + (-t248 * t235 - (-t248 * qJD(4) - t242 * t298 + t271 * t297) * t250 - t249 * t234) * t246) * t236, 0, -t325 * t369 * t341 + (t325 * t239 - ((-qJD(4) * t245 - 0.2e1 * t339) * t298 + (t234 * t298 + (t235 - t344) * t297) * t246) * t369) * t236, t341 + 0.2e1 * (t234 * t246 * t236 + (-t236 * t352 - t246 * t351) * t250) * t250, 0, 0;];
JaD_rot  = t1;
