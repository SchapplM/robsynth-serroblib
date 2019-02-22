% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:05
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobiR_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:05:18
% EndTime: 2019-02-22 10:05:19
% DurationCPUTime: 0.27s
% Computational Cost: add. (391->75), mult. (1178->165), div. (0->0), fcn. (1599->16), ass. (0->80)
t264 = sin(pkin(8));
t265 = sin(pkin(7));
t304 = t264 * t265;
t271 = sin(qJ(5));
t303 = t264 * t271;
t275 = cos(qJ(5));
t302 = t264 * t275;
t266 = sin(pkin(6));
t301 = t265 * t266;
t268 = cos(pkin(8));
t300 = t265 * t268;
t270 = cos(pkin(6));
t299 = t265 * t270;
t269 = cos(pkin(7));
t298 = t266 * t269;
t272 = sin(qJ(4));
t297 = t268 * t272;
t276 = cos(qJ(4));
t296 = t268 * t276;
t273 = sin(qJ(3));
t295 = t269 * t273;
t277 = cos(qJ(3));
t294 = t269 * t277;
t274 = sin(qJ(2));
t293 = t270 * t274;
t278 = cos(qJ(2));
t292 = t270 * t278;
t291 = t273 * t274;
t290 = t273 * t278;
t289 = t274 * t277;
t288 = t277 * t278;
t267 = cos(pkin(14));
t287 = t267 * t301;
t286 = t274 * t301;
t263 = sin(pkin(14));
t257 = -t263 * t274 + t267 * t292;
t258 = t263 * t278 + t267 * t293;
t241 = -t258 * t273 + (t257 * t269 - t287) * t277;
t252 = -t257 * t265 - t267 * t298;
t285 = t241 * t268 + t252 * t264;
t260 = -t263 * t293 + t267 * t278;
t259 = -t263 * t292 - t267 * t274;
t280 = t259 * t269 + t263 * t301;
t243 = -t260 * t273 + t280 * t277;
t253 = -t259 * t265 + t263 * t298;
t284 = t243 * t268 + t253 * t264;
t250 = t277 * t299 + (t269 * t288 - t291) * t266;
t256 = t270 * t269 - t278 * t301;
t283 = t250 * t268 + t256 * t264;
t245 = -t257 * t273 - t258 * t294;
t282 = t245 * t268 + t258 * t304;
t247 = -t259 * t273 - t260 * t294;
t281 = t247 * t268 + t260 * t304;
t254 = (-t269 * t289 - t290) * t266;
t279 = t254 * t268 + t264 * t286;
t255 = (-t269 * t291 + t288) * t266;
t251 = t273 * t299 + (t269 * t290 + t289) * t266;
t249 = -t254 * t264 + t268 * t286;
t248 = t259 * t277 - t260 * t295;
t246 = t257 * t277 - t258 * t295;
t244 = t260 * t277 + t280 * t273;
t242 = t257 * t295 + t258 * t277 - t273 * t287;
t240 = -t250 * t264 + t256 * t268;
t239 = t255 * t276 + t279 * t272;
t238 = -t247 * t264 + t260 * t300;
t237 = -t245 * t264 + t258 * t300;
t236 = t250 * t276 - t251 * t297;
t235 = -t243 * t264 + t253 * t268;
t234 = -t241 * t264 + t252 * t268;
t233 = t251 * t276 + t283 * t272;
t232 = -t251 * t272 + t283 * t276;
t231 = t243 * t276 - t244 * t297;
t230 = t241 * t276 - t242 * t297;
t229 = t248 * t276 + t281 * t272;
t228 = t246 * t276 + t282 * t272;
t227 = t244 * t276 + t284 * t272;
t226 = -t244 * t272 + t284 * t276;
t225 = t242 * t276 + t285 * t272;
t224 = -t242 * t272 + t285 * t276;
t1 = [0, t229 * t275 + t238 * t271, t231 * t275 + t244 * t303, t226 * t275, -t227 * t271 + t235 * t275, 0; 0, t228 * t275 + t237 * t271, t230 * t275 + t242 * t303, t224 * t275, -t225 * t271 + t234 * t275, 0; 0, t239 * t275 + t249 * t271, t236 * t275 + t251 * t303, t232 * t275, -t233 * t271 + t240 * t275, 0; 0, -t229 * t271 + t238 * t275, -t231 * t271 + t244 * t302, -t226 * t271, -t227 * t275 - t235 * t271, 0; 0, -t228 * t271 + t237 * t275, -t230 * t271 + t242 * t302, -t224 * t271, -t225 * t275 - t234 * t271, 0; 0, -t239 * t271 + t249 * t275, -t236 * t271 + t251 * t302, -t232 * t271, -t233 * t275 - t240 * t271, 0; 0, t248 * t272 - t281 * t276, t243 * t272 + t244 * t296, t227, 0, 0; 0, t246 * t272 - t282 * t276, t241 * t272 + t242 * t296, t225, 0, 0; 0, t255 * t272 - t279 * t276, t250 * t272 + t251 * t296, t233, 0, 0;];
JR_rot  = t1;
