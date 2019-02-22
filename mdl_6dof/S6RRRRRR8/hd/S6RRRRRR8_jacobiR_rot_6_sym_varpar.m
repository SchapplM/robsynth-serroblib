% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:37
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:37:43
% EndTime: 2019-02-22 12:37:44
% DurationCPUTime: 0.44s
% Computational Cost: add. (437->71), mult. (1029->135), div. (0->0), fcn. (1421->14), ass. (0->76)
t292 = sin(qJ(2));
t295 = cos(qJ(2));
t296 = cos(qJ(1));
t323 = cos(pkin(6));
t302 = t296 * t323;
t324 = sin(qJ(1));
t276 = t292 * t302 + t324 * t295;
t291 = sin(qJ(3));
t294 = cos(qJ(3));
t275 = t324 * t292 - t295 * t302;
t289 = cos(pkin(7));
t287 = sin(pkin(7));
t288 = sin(pkin(6));
t312 = t288 * t296;
t304 = t287 * t312;
t298 = t275 * t289 + t304;
t257 = -t276 * t294 + t298 * t291;
t268 = -t275 * t287 + t289 * t312;
t286 = qJ(4) + qJ(5);
t284 = sin(t286);
t285 = cos(t286);
t246 = t257 * t285 + t268 * t284;
t290 = sin(qJ(6));
t328 = t246 * t290;
t293 = cos(qJ(6));
t327 = t246 * t293;
t244 = t257 * t284 - t268 * t285;
t322 = t244 * t290;
t299 = t323 * t324;
t277 = -t296 * t292 - t295 * t299;
t278 = -t292 * t299 + t296 * t295;
t303 = t288 * t324;
t300 = t287 * t303;
t259 = t278 * t294 + (t277 * t289 + t300) * t291;
t297 = -t277 * t287 + t289 * t303;
t247 = t259 * t284 - t297 * t285;
t321 = t247 * t290;
t301 = t323 * t287;
t307 = t292 * t294;
t308 = t291 * t295;
t267 = t291 * t301 + (t289 * t308 + t307) * t288;
t313 = t287 * t288;
t274 = t323 * t289 - t295 * t313;
t252 = -t267 * t284 + t274 * t285;
t320 = t252 * t290;
t319 = t276 * t291;
t317 = t284 * t287;
t316 = t285 * t287;
t315 = t285 * t290;
t314 = t285 * t293;
t311 = t289 * t291;
t310 = t289 * t294;
t309 = t291 * t292;
t306 = t294 * t295;
t305 = t292 * t313;
t273 = (-t289 * t309 + t306) * t288;
t272 = (t289 * t307 + t308) * t288;
t266 = -t294 * t301 + (-t289 * t306 + t309) * t288;
t264 = t273 * t285 + t284 * t305;
t263 = t277 * t294 - t278 * t311;
t262 = t277 * t291 + t278 * t310;
t261 = -t275 * t294 - t276 * t311;
t260 = -t275 * t291 + t276 * t310;
t258 = -t277 * t310 + t278 * t291 - t294 * t300;
t256 = -t298 * t294 - t319;
t254 = t275 * t310 + t294 * t304 + t319;
t253 = t267 * t285 + t274 * t284;
t251 = t252 * t293;
t250 = t263 * t285 + t278 * t317;
t249 = t261 * t285 + t276 * t317;
t248 = t259 * t285 + t297 * t284;
t243 = t247 * t293;
t242 = t244 * t293;
t241 = t248 * t293 + t258 * t290;
t240 = -t248 * t290 + t258 * t293;
t1 = [t256 * t290 + t327, t250 * t293 + t262 * t290, -t258 * t314 + t259 * t290, -t243, -t243, t240; t241, t249 * t293 + t260 * t290, -t254 * t314 - t257 * t290, t242, t242, t254 * t293 + t328; 0, t264 * t293 + t272 * t290, -t266 * t314 + t267 * t290, t251, t251, -t253 * t290 + t266 * t293; t256 * t293 - t328, -t250 * t290 + t262 * t293, t258 * t315 + t259 * t293, t321, t321, -t241; t240, -t249 * t290 + t260 * t293, t254 * t315 - t257 * t293, -t322, -t322, -t254 * t290 + t327; 0, -t264 * t290 + t272 * t293, t266 * t315 + t267 * t293, -t320, -t320, -t253 * t293 - t266 * t290; t244, t263 * t284 - t278 * t316, -t258 * t284, t248, t248, 0; t247, t261 * t284 - t276 * t316, -t254 * t284, -t246, -t246, 0; 0, t273 * t284 - t285 * t305, -t266 * t284, t253, t253, 0;];
JR_rot  = t1;
