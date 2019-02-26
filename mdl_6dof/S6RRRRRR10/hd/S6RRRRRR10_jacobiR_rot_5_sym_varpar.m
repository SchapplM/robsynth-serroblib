% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR10_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiR_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:48
% EndTime: 2019-02-26 22:52:48
% DurationCPUTime: 0.63s
% Computational Cost: add. (495->79), mult. (1476->167), div. (0->0), fcn. (2007->16), ass. (0->85)
t321 = sin(qJ(2));
t322 = sin(qJ(1));
t326 = cos(qJ(2));
t327 = cos(qJ(1));
t358 = cos(pkin(6));
t337 = t327 * t358;
t306 = t322 * t321 - t326 * t337;
t307 = t321 * t337 + t322 * t326;
t320 = sin(qJ(3));
t325 = cos(qJ(3));
t314 = sin(pkin(7));
t315 = sin(pkin(6));
t349 = t315 * t327;
t339 = t314 * t349;
t317 = cos(pkin(7));
t346 = t317 * t320;
t289 = t306 * t346 - t307 * t325 + t320 * t339;
t319 = sin(qJ(4));
t324 = cos(qJ(4));
t288 = t325 * (t306 * t317 + t339) + t307 * t320;
t300 = -t306 * t314 + t317 * t349;
t313 = sin(pkin(8));
t316 = cos(pkin(8));
t334 = t288 * t316 + t300 * t313;
t268 = t289 * t324 + t319 * t334;
t277 = t288 * t313 - t300 * t316;
t318 = sin(qJ(5));
t323 = cos(qJ(5));
t366 = t268 * t318 + t277 * t323;
t365 = t268 * t323 - t277 * t318;
t266 = t289 * t319 - t324 * t334;
t338 = t322 * t358;
t308 = -t327 * t321 - t326 * t338;
t350 = t315 * t322;
t302 = -t308 * t314 + t317 * t350;
t357 = t302 * t313;
t355 = t313 * t314;
t354 = t313 * t318;
t353 = t313 * t323;
t352 = t314 * t315;
t351 = t314 * t316;
t348 = t316 * t319;
t347 = t316 * t324;
t345 = t317 * t325;
t344 = t320 * t321;
t343 = t320 * t326;
t342 = t321 * t325;
t341 = t325 * t326;
t340 = t321 * t352;
t336 = t358 * t314;
t298 = t325 * t336 + (t317 * t341 - t344) * t315;
t305 = t317 * t358 - t326 * t352;
t333 = t298 * t316 + t305 * t313;
t292 = t306 * t320 - t307 * t345;
t332 = t292 * t316 + t307 * t355;
t309 = -t321 * t338 + t327 * t326;
t294 = -t308 * t320 - t309 * t345;
t331 = t294 * t316 + t309 * t355;
t329 = t308 * t317 + t314 * t350;
t303 = (-t317 * t342 - t343) * t315;
t328 = t303 * t316 + t313 * t340;
t304 = (-t317 * t344 + t341) * t315;
t299 = t320 * t336 + (t317 * t343 + t342) * t315;
t296 = -t303 * t313 + t316 * t340;
t295 = t308 * t325 - t309 * t346;
t293 = -t306 * t325 - t307 * t346;
t291 = t309 * t325 + t320 * t329;
t290 = -t309 * t320 + t325 * t329;
t285 = -t298 * t313 + t305 * t316;
t283 = -t294 * t313 + t309 * t351;
t282 = -t292 * t313 + t307 * t351;
t281 = t304 * t324 + t319 * t328;
t280 = t298 * t324 - t299 * t348;
t279 = -t290 * t313 + t302 * t316;
t276 = t299 * t324 + t319 * t333;
t275 = -t299 * t319 + t324 * t333;
t274 = t290 * t324 - t291 * t348;
t273 = -t288 * t324 + t289 * t348;
t272 = t295 * t324 + t319 * t331;
t271 = t293 * t324 + t319 * t332;
t270 = t291 * t324 + (t290 * t316 + t357) * t319;
t269 = -t290 * t347 + t291 * t319 - t324 * t357;
t265 = t270 * t323 + t279 * t318;
t264 = -t270 * t318 + t279 * t323;
t1 = [t365, t272 * t323 + t283 * t318, t274 * t323 + t291 * t354, -t269 * t323, t264, 0; t265, t271 * t323 + t282 * t318, t273 * t323 - t289 * t354, t266 * t323, t366, 0; 0, t281 * t323 + t296 * t318, t280 * t323 + t299 * t354, t275 * t323, -t276 * t318 + t285 * t323, 0; -t366, -t272 * t318 + t283 * t323, -t274 * t318 + t291 * t353, t269 * t318, -t265, 0; t264, -t271 * t318 + t282 * t323, -t273 * t318 - t289 * t353, -t266 * t318, t365, 0; 0, -t281 * t318 + t296 * t323, -t280 * t318 + t299 * t353, -t275 * t318, -t276 * t323 - t285 * t318, 0; t266, t295 * t319 - t324 * t331, t290 * t319 + t291 * t347, t270, 0, 0; t269, t293 * t319 - t324 * t332, -t288 * t319 - t289 * t347, -t268, 0, 0; 0, t304 * t319 - t324 * t328, t298 * t319 + t299 * t347, t276, 0, 0;];
JR_rot  = t1;
