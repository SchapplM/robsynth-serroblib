% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:32
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP12_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:32:43
% EndTime: 2019-02-22 12:32:43
% DurationCPUTime: 0.41s
% Computational Cost: add. (290->65), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->69)
t295 = sin(qJ(2));
t296 = sin(qJ(1));
t300 = cos(qJ(2));
t301 = cos(qJ(1));
t325 = cos(pkin(6));
t304 = t301 * t325;
t281 = t295 * t304 + t296 * t300;
t294 = sin(qJ(3));
t299 = cos(qJ(3));
t280 = t296 * t295 - t300 * t304;
t291 = cos(pkin(7));
t289 = sin(pkin(7));
t290 = sin(pkin(6));
t318 = t290 * t301;
t307 = t289 * t318;
t303 = t280 * t291 + t307;
t262 = -t281 * t299 + t303 * t294;
t273 = -t280 * t289 + t291 * t318;
t293 = sin(qJ(4));
t298 = cos(qJ(4));
t252 = t262 * t298 + t273 * t293;
t292 = sin(qJ(5));
t329 = t252 * t292;
t297 = cos(qJ(5));
t328 = t252 * t297;
t250 = t262 * t293 - t273 * t298;
t324 = t281 * t294;
t322 = t289 * t290;
t321 = t289 * t293;
t320 = t289 * t298;
t319 = t290 * t296;
t317 = t291 * t294;
t316 = t291 * t299;
t315 = t292 * t298;
t314 = t294 * t295;
t313 = t294 * t300;
t312 = t295 * t299;
t311 = t297 * t298;
t310 = t299 * t300;
t309 = t295 * t322;
t308 = t289 * t319;
t306 = t289 * t325;
t305 = t296 * t325;
t282 = -t301 * t295 - t300 * t305;
t302 = -t282 * t289 + t291 * t319;
t283 = -t295 * t305 + t301 * t300;
t279 = t325 * t291 - t300 * t322;
t278 = (-t291 * t314 + t310) * t290;
t277 = (t291 * t312 + t313) * t290;
t272 = t294 * t306 + (t291 * t313 + t312) * t290;
t271 = -t299 * t306 + (-t291 * t310 + t314) * t290;
t269 = t278 * t298 + t293 * t309;
t268 = t282 * t299 - t283 * t317;
t267 = t282 * t294 + t283 * t316;
t266 = -t280 * t299 - t281 * t317;
t265 = -t280 * t294 + t281 * t316;
t264 = t283 * t299 + (t282 * t291 + t308) * t294;
t263 = -t282 * t316 + t283 * t294 - t299 * t308;
t261 = -t303 * t299 - t324;
t259 = t280 * t316 + t299 * t307 + t324;
t258 = t272 * t298 + t279 * t293;
t257 = -t272 * t293 + t279 * t298;
t256 = t268 * t298 + t283 * t321;
t255 = t266 * t298 + t281 * t321;
t254 = t264 * t298 + t302 * t293;
t253 = t264 * t293 - t302 * t298;
t249 = t254 * t297 + t263 * t292;
t248 = t254 * t292 - t263 * t297;
t1 = [t261 * t292 + t328, t256 * t297 + t267 * t292, -t263 * t311 + t264 * t292, -t253 * t297, -t248, 0; t249, t255 * t297 + t265 * t292, -t259 * t311 - t262 * t292, t250 * t297, t259 * t297 + t329, 0; 0, t269 * t297 + t277 * t292, -t271 * t311 + t272 * t292, t257 * t297, -t258 * t292 + t271 * t297, 0; t250, t268 * t293 - t283 * t320, -t263 * t293, t254, 0, 0; t253, t266 * t293 - t281 * t320, -t259 * t293, -t252, 0, 0; 0, t278 * t293 - t298 * t309, -t271 * t293, t258, 0, 0; -t261 * t297 + t329, t256 * t292 - t267 * t297, -t263 * t315 - t264 * t297, -t253 * t292, t249, 0; t248, t255 * t292 - t265 * t297, -t259 * t315 + t262 * t297, t250 * t292, t259 * t292 - t328, 0; 0, t269 * t292 - t277 * t297, -t271 * t315 - t272 * t297, t257 * t292, t258 * t297 + t271 * t292, 0;];
JR_rot  = t1;
