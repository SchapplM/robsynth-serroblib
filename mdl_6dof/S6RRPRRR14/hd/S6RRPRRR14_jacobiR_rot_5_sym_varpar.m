% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiR_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:46
% EndTime: 2019-02-26 22:55:46
% DurationCPUTime: 0.53s
% Computational Cost: add. (400->65), mult. (1183->139), div. (0->0), fcn. (1610->16), ass. (0->76)
t298 = sin(qJ(2));
t299 = sin(qJ(1));
t302 = cos(qJ(2));
t303 = cos(qJ(1));
t329 = cos(pkin(6));
t313 = t303 * t329;
t284 = t298 * t313 + t299 * t302;
t289 = sin(pkin(14));
t293 = cos(pkin(14));
t283 = t299 * t298 - t302 * t313;
t291 = sin(pkin(7));
t295 = cos(pkin(7));
t292 = sin(pkin(6));
t319 = t292 * t303;
t306 = t283 * t295 + t291 * t319;
t266 = -t284 * t293 + t306 * t289;
t297 = sin(qJ(4));
t301 = cos(qJ(4));
t265 = t284 * t289 + t306 * t293;
t277 = -t283 * t291 + t295 * t319;
t290 = sin(pkin(8));
t294 = cos(pkin(8));
t310 = t265 * t294 + t277 * t290;
t248 = t266 * t301 + t310 * t297;
t255 = t265 * t290 - t277 * t294;
t296 = sin(qJ(5));
t300 = cos(qJ(5));
t338 = t248 * t296 + t255 * t300;
t337 = t248 * t300 - t255 * t296;
t246 = t266 * t297 - t310 * t301;
t314 = t299 * t329;
t286 = -t298 * t314 + t303 * t302;
t285 = -t303 * t298 - t302 * t314;
t320 = t292 * t299;
t305 = t285 * t295 + t291 * t320;
t267 = -t286 * t289 + t305 * t293;
t279 = -t285 * t291 + t295 * t320;
t330 = t267 * t294 + t279 * t290;
t324 = t289 * t295;
t323 = t290 * t291;
t322 = t291 * t294;
t321 = t292 * t298;
t318 = t293 * t295;
t317 = t295 * t298;
t316 = t295 * t302;
t315 = t291 * t321;
t312 = t329 * t291;
t275 = t293 * t312 + (-t289 * t298 + t293 * t316) * t292;
t282 = -t292 * t302 * t291 + t329 * t295;
t309 = t275 * t294 + t282 * t290;
t269 = t283 * t289 - t284 * t318;
t308 = t269 * t294 + t284 * t323;
t271 = -t285 * t289 - t286 * t318;
t307 = t271 * t294 + t286 * t323;
t280 = (-t289 * t302 - t293 * t317) * t292;
t304 = t280 * t294 + t290 * t315;
t281 = (-t289 * t317 + t293 * t302) * t292;
t276 = t293 * t321 + (t292 * t316 + t312) * t289;
t273 = -t280 * t290 + t294 * t315;
t272 = t285 * t293 - t286 * t324;
t270 = -t283 * t293 - t284 * t324;
t268 = t286 * t293 + t305 * t289;
t262 = -t275 * t290 + t282 * t294;
t260 = -t271 * t290 + t286 * t322;
t259 = -t269 * t290 + t284 * t322;
t258 = t281 * t301 + t304 * t297;
t257 = -t267 * t290 + t279 * t294;
t254 = t276 * t301 + t309 * t297;
t253 = -t276 * t297 + t309 * t301;
t252 = t272 * t301 + t307 * t297;
t251 = t270 * t301 + t308 * t297;
t250 = t268 * t301 + t330 * t297;
t249 = t268 * t297 - t330 * t301;
t245 = t250 * t300 + t257 * t296;
t244 = -t250 * t296 + t257 * t300;
t1 = [t337, t252 * t300 + t260 * t296, 0, -t249 * t300, t244, 0; t245, t251 * t300 + t259 * t296, 0, t246 * t300, t338, 0; 0, t258 * t300 + t273 * t296, 0, t253 * t300, -t254 * t296 + t262 * t300, 0; -t338, -t252 * t296 + t260 * t300, 0, t249 * t296, -t245, 0; t244, -t251 * t296 + t259 * t300, 0, -t246 * t296, t337, 0; 0, -t258 * t296 + t273 * t300, 0, -t253 * t296, -t254 * t300 - t262 * t296, 0; t246, t272 * t297 - t307 * t301, 0, t250, 0, 0; t249, t270 * t297 - t308 * t301, 0, -t248, 0, 0; 0, t281 * t297 - t304 * t301, 0, t254, 0, 0;];
JR_rot  = t1;
