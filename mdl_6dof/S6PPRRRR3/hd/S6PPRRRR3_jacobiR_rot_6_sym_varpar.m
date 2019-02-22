% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:27
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiR_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:27:12
% EndTime: 2019-02-22 09:27:12
% DurationCPUTime: 0.27s
% Computational Cost: add. (652->75), mult. (1920->156), div. (0->0), fcn. (2602->18), ass. (0->73)
t274 = sin(pkin(13));
t282 = cos(pkin(6));
t310 = t274 * t282;
t275 = sin(pkin(8));
t284 = sin(qJ(5));
t309 = t275 * t284;
t288 = cos(qJ(5));
t308 = t275 * t288;
t276 = sin(pkin(7));
t277 = sin(pkin(6));
t307 = t276 * t277;
t306 = t276 * t282;
t281 = cos(pkin(7));
t305 = t277 * t281;
t279 = cos(pkin(13));
t304 = t279 * t282;
t280 = cos(pkin(8));
t285 = sin(qJ(4));
t303 = t280 * t285;
t289 = cos(qJ(4));
t302 = t280 * t289;
t286 = sin(qJ(3));
t301 = t281 * t286;
t283 = sin(qJ(6));
t300 = t283 * t288;
t287 = cos(qJ(6));
t299 = t287 * t288;
t298 = t279 * t307;
t273 = sin(pkin(14));
t278 = cos(pkin(14));
t268 = -t274 * t273 + t278 * t304;
t297 = -t268 * t276 - t279 * t305;
t270 = -t279 * t273 - t278 * t310;
t296 = -t270 * t276 + t274 * t305;
t295 = t270 * t281 + t274 * t307;
t294 = -t278 * t307 + t282 * t281;
t293 = t297 * t275;
t292 = t296 * t275;
t291 = t294 * t275;
t290 = cos(qJ(3));
t271 = -t273 * t310 + t279 * t278;
t269 = t273 * t304 + t274 * t278;
t266 = t286 * t306 + (t273 * t290 + t278 * t301) * t277;
t265 = t290 * t306 + (t278 * t281 * t290 - t273 * t286) * t277;
t261 = -t265 * t275 + t294 * t280;
t260 = t271 * t290 + t295 * t286;
t259 = -t271 * t286 + t295 * t290;
t258 = t268 * t301 + t269 * t290 - t286 * t298;
t257 = -t269 * t286 + (t268 * t281 - t298) * t290;
t254 = t265 * t289 - t266 * t303;
t253 = t265 * t285 + t266 * t302;
t252 = -t259 * t275 + t296 * t280;
t251 = -t257 * t275 + t297 * t280;
t250 = t266 * t289 + (t265 * t280 + t291) * t285;
t249 = -t265 * t302 + t266 * t285 - t289 * t291;
t248 = t254 * t288 + t266 * t309;
t247 = t259 * t289 - t260 * t303;
t246 = t259 * t285 + t260 * t302;
t245 = t257 * t289 - t258 * t303;
t244 = t257 * t285 + t258 * t302;
t243 = t250 * t288 + t261 * t284;
t242 = -t250 * t284 + t261 * t288;
t241 = t260 * t289 + (t259 * t280 + t292) * t285;
t240 = -t259 * t302 + t260 * t285 - t289 * t292;
t239 = t258 * t289 + (t257 * t280 + t293) * t285;
t238 = -t257 * t302 + t258 * t285 - t289 * t293;
t237 = t247 * t288 + t260 * t309;
t236 = t245 * t288 + t258 * t309;
t235 = t241 * t288 + t252 * t284;
t234 = -t241 * t284 + t252 * t288;
t233 = t239 * t288 + t251 * t284;
t232 = -t239 * t284 + t251 * t288;
t1 = [0, 0, t237 * t287 + t246 * t283, -t240 * t299 + t241 * t283, t234 * t287, -t235 * t283 + t240 * t287; 0, 0, t236 * t287 + t244 * t283, -t238 * t299 + t239 * t283, t232 * t287, -t233 * t283 + t238 * t287; 0, 0, t248 * t287 + t253 * t283, -t249 * t299 + t250 * t283, t242 * t287, -t243 * t283 + t249 * t287; 0, 0, -t237 * t283 + t246 * t287, t240 * t300 + t241 * t287, -t234 * t283, -t235 * t287 - t240 * t283; 0, 0, -t236 * t283 + t244 * t287, t238 * t300 + t239 * t287, -t232 * t283, -t233 * t287 - t238 * t283; 0, 0, -t248 * t283 + t253 * t287, t249 * t300 + t250 * t287, -t242 * t283, -t243 * t287 - t249 * t283; 0, 0, t247 * t284 - t260 * t308, -t240 * t284, t235, 0; 0, 0, t245 * t284 - t258 * t308, -t238 * t284, t233, 0; 0, 0, t254 * t284 - t266 * t308, -t249 * t284, t243, 0;];
JR_rot  = t1;
