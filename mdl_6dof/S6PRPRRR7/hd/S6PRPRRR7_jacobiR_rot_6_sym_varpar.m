% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:41
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiR_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:41:41
% EndTime: 2019-02-22 09:41:42
% DurationCPUTime: 0.40s
% Computational Cost: add. (653->83), mult. (1935->179), div. (0->0), fcn. (2618->18), ass. (0->88)
t307 = sin(pkin(8));
t312 = cos(pkin(8));
t305 = sin(pkin(14));
t309 = sin(pkin(6));
t310 = cos(pkin(14));
t318 = sin(qJ(2));
t313 = cos(pkin(7));
t322 = cos(qJ(2));
t345 = t313 * t322;
t308 = sin(pkin(7));
t314 = cos(pkin(6));
t350 = t308 * t314;
t326 = t310 * t350 + (-t305 * t318 + t310 * t345) * t309;
t352 = t308 * t309;
t333 = t314 * t313 - t322 * t352;
t357 = t333 * t307 + t326 * t312;
t306 = sin(pkin(13));
t311 = cos(pkin(13));
t344 = t314 * t318;
t303 = -t306 * t344 + t311 * t322;
t343 = t314 * t322;
t302 = -t306 * t343 - t311 * t318;
t334 = t302 * t313 + t306 * t352;
t327 = t303 * t305 - t334 * t310;
t349 = t309 * t313;
t335 = -t302 * t308 + t306 * t349;
t356 = -t335 * t307 + t327 * t312;
t301 = t306 * t322 + t311 * t344;
t300 = -t306 * t318 + t311 * t343;
t336 = t300 * t313 - t311 * t352;
t328 = t301 * t305 - t336 * t310;
t337 = -t300 * t308 - t311 * t349;
t355 = -t337 * t307 + t328 * t312;
t354 = t305 * t313;
t353 = t308 * t307;
t351 = t308 * t312;
t348 = t309 * t318;
t347 = t310 * t313;
t346 = t313 * t318;
t315 = sin(qJ(6));
t320 = cos(qJ(5));
t342 = t315 * t320;
t319 = cos(qJ(6));
t341 = t319 * t320;
t340 = t308 * t348;
t288 = -t300 * t305 - t301 * t347;
t339 = t288 * t312 + t301 * t353;
t290 = -t302 * t305 - t303 * t347;
t338 = t290 * t312 + t303 * t353;
t298 = (-t305 * t322 - t310 * t346) * t309;
t332 = t298 * t312 + t307 * t340;
t321 = cos(qJ(4));
t317 = sin(qJ(4));
t316 = sin(qJ(5));
t299 = (-t305 * t346 + t310 * t322) * t309;
t296 = t310 * t348 + (t309 * t345 + t350) * t305;
t292 = -t298 * t307 + t312 * t340;
t291 = t302 * t310 - t303 * t354;
t289 = t300 * t310 - t301 * t354;
t287 = t303 * t310 + t334 * t305;
t286 = t301 * t310 + t336 * t305;
t285 = -t326 * t307 + t333 * t312;
t282 = t299 * t321 + t332 * t317;
t281 = t299 * t317 - t332 * t321;
t280 = -t290 * t307 + t303 * t351;
t279 = -t288 * t307 + t301 * t351;
t278 = t327 * t307 + t335 * t312;
t277 = t328 * t307 + t337 * t312;
t276 = t296 * t321 + t357 * t317;
t275 = t296 * t317 - t357 * t321;
t274 = t282 * t320 + t292 * t316;
t273 = t291 * t321 + t338 * t317;
t272 = t291 * t317 - t338 * t321;
t271 = t289 * t321 + t339 * t317;
t270 = t289 * t317 - t339 * t321;
t269 = t287 * t321 - t356 * t317;
t268 = t287 * t317 + t356 * t321;
t267 = t286 * t321 - t355 * t317;
t266 = t286 * t317 + t355 * t321;
t265 = t276 * t320 + t285 * t316;
t264 = -t276 * t316 + t285 * t320;
t263 = t273 * t320 + t280 * t316;
t262 = t271 * t320 + t279 * t316;
t261 = t269 * t320 + t278 * t316;
t260 = -t269 * t316 + t278 * t320;
t259 = t267 * t320 + t277 * t316;
t258 = -t267 * t316 + t277 * t320;
t1 = [0, t263 * t319 + t272 * t315, 0, -t268 * t341 + t269 * t315, t260 * t319, -t261 * t315 + t268 * t319; 0, t262 * t319 + t270 * t315, 0, -t266 * t341 + t267 * t315, t258 * t319, -t259 * t315 + t266 * t319; 0, t274 * t319 + t281 * t315, 0, -t275 * t341 + t276 * t315, t264 * t319, -t265 * t315 + t275 * t319; 0, -t263 * t315 + t272 * t319, 0, t268 * t342 + t269 * t319, -t260 * t315, -t261 * t319 - t268 * t315; 0, -t262 * t315 + t270 * t319, 0, t266 * t342 + t267 * t319, -t258 * t315, -t259 * t319 - t266 * t315; 0, -t274 * t315 + t281 * t319, 0, t275 * t342 + t276 * t319, -t264 * t315, -t265 * t319 - t275 * t315; 0, t273 * t316 - t280 * t320, 0, -t268 * t316, t261, 0; 0, t271 * t316 - t279 * t320, 0, -t266 * t316, t259, 0; 0, t282 * t316 - t292 * t320, 0, -t275 * t316, t265, 0;];
JR_rot  = t1;
