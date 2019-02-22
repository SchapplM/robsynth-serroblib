% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR9
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
% Datum: 2019-02-22 12:38
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:38:22
% EndTime: 2019-02-22 12:38:23
% DurationCPUTime: 0.44s
% Computational Cost: add. (417->68), mult. (1045->135), div. (0->0), fcn. (1442->14), ass. (0->74)
t280 = sin(qJ(2));
t281 = sin(qJ(1));
t284 = cos(qJ(2));
t285 = cos(qJ(1));
t309 = cos(pkin(6));
t289 = t285 * t309;
t264 = t280 * t289 + t281 * t284;
t279 = sin(qJ(3));
t283 = cos(qJ(3));
t263 = t281 * t280 - t284 * t289;
t277 = cos(pkin(7));
t275 = sin(pkin(7));
t276 = sin(pkin(6));
t300 = t276 * t285;
t291 = t275 * t300;
t287 = t263 * t277 + t291;
t245 = -t264 * t283 + t279 * t287;
t256 = -t263 * t275 + t277 * t300;
t278 = sin(qJ(4));
t282 = cos(qJ(4));
t235 = t245 * t282 + t256 * t278;
t274 = qJ(5) + qJ(6);
t272 = sin(t274);
t313 = t235 * t272;
t273 = cos(t274);
t312 = t235 * t273;
t233 = t245 * t278 - t256 * t282;
t308 = t264 * t279;
t306 = t272 * t282;
t305 = t273 * t282;
t304 = t275 * t276;
t303 = t275 * t278;
t302 = t275 * t282;
t301 = t276 * t281;
t299 = t277 * t279;
t298 = t277 * t283;
t297 = t279 * t280;
t296 = t279 * t284;
t295 = t280 * t283;
t294 = t283 * t284;
t293 = t280 * t304;
t292 = t275 * t301;
t290 = t281 * t309;
t288 = t309 * t275;
t265 = -t285 * t280 - t284 * t290;
t286 = -t265 * t275 + t277 * t301;
t266 = -t280 * t290 + t285 * t284;
t262 = t277 * t309 - t284 * t304;
t261 = (-t277 * t297 + t294) * t276;
t260 = (t277 * t295 + t296) * t276;
t255 = t279 * t288 + (t277 * t296 + t295) * t276;
t254 = -t283 * t288 + (-t277 * t294 + t297) * t276;
t252 = t261 * t282 + t278 * t293;
t251 = t265 * t283 - t266 * t299;
t250 = t265 * t279 + t266 * t298;
t249 = -t263 * t283 - t264 * t299;
t248 = -t263 * t279 + t264 * t298;
t247 = t266 * t283 + (t265 * t277 + t292) * t279;
t246 = -t265 * t298 + t266 * t279 - t283 * t292;
t244 = -t283 * t287 - t308;
t242 = t263 * t298 + t283 * t291 + t308;
t241 = t255 * t282 + t262 * t278;
t240 = -t255 * t278 + t262 * t282;
t239 = t251 * t282 + t266 * t303;
t238 = t249 * t282 + t264 * t303;
t237 = t247 * t282 + t278 * t286;
t236 = t247 * t278 - t282 * t286;
t232 = -t241 * t273 - t254 * t272;
t231 = -t241 * t272 + t254 * t273;
t230 = t237 * t273 + t246 * t272;
t229 = -t237 * t272 + t246 * t273;
t228 = -t242 * t272 + t312;
t227 = t242 * t273 + t313;
t1 = [t244 * t272 + t312, t239 * t273 + t250 * t272, -t246 * t305 + t247 * t272, -t236 * t273, t229, t229; t230, t238 * t273 + t248 * t272, -t242 * t305 - t245 * t272, t233 * t273, t227, t227; 0, t252 * t273 + t260 * t272, -t254 * t305 + t255 * t272, t240 * t273, t231, t231; t244 * t273 - t313, -t239 * t272 + t250 * t273, t246 * t306 + t247 * t273, t236 * t272, -t230, -t230; t229, -t238 * t272 + t248 * t273, t242 * t306 - t245 * t273, -t233 * t272, t228, t228; 0, -t252 * t272 + t260 * t273, t254 * t306 + t255 * t273, -t240 * t272, t232, t232; t233, t251 * t278 - t266 * t302, -t246 * t278, t237, 0, 0; t236, t249 * t278 - t264 * t302, -t242 * t278, -t235, 0, 0; 0, t261 * t278 - t282 * t293, -t254 * t278, t241, 0, 0;];
JR_rot  = t1;
