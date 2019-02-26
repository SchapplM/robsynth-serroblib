% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:38
% EndTime: 2019-02-26 22:20:38
% DurationCPUTime: 0.35s
% Computational Cost: add. (461->63), mult. (1310->135), div. (0->0), fcn. (1807->16), ass. (0->65)
t270 = sin(pkin(13));
t276 = sin(qJ(3));
t298 = cos(pkin(13));
t300 = cos(qJ(3));
t266 = -t270 * t300 - t276 * t298;
t271 = sin(pkin(7));
t257 = t266 * t271;
t273 = cos(pkin(7));
t259 = t266 * t273;
t277 = sin(qJ(2));
t278 = sin(qJ(1));
t281 = cos(qJ(2));
t282 = cos(qJ(1));
t299 = cos(pkin(6));
t288 = t282 * t299;
t261 = t278 * t277 - t281 * t288;
t262 = t277 * t288 + t278 * t281;
t285 = -t276 * t270 + t298 * t300;
t272 = sin(pkin(6));
t293 = t272 * t282;
t238 = -t257 * t293 - t261 * t259 - t262 * t285;
t253 = -t261 * t271 + t273 * t293;
t275 = sin(qJ(5));
t280 = cos(qJ(5));
t227 = t238 * t280 + t253 * t275;
t274 = sin(qJ(6));
t279 = cos(qJ(6));
t256 = t285 * t271;
t258 = t285 * t273;
t286 = -t256 * t293 - t261 * t258 + t262 * t266;
t304 = t227 * t274 - t279 * t286;
t303 = t227 * t279 + t274 * t286;
t225 = t238 * t275 - t253 * t280;
t297 = t271 * t272;
t296 = t271 * t275;
t295 = t271 * t280;
t294 = t272 * t278;
t292 = t274 * t280;
t291 = t279 * t280;
t290 = t277 * t297;
t289 = t278 * t299;
t263 = -t282 * t277 - t281 * t289;
t287 = -t263 * t271 + t273 * t294;
t264 = -t277 * t289 + t282 * t281;
t284 = -t257 * t294 - t263 * t259 + t264 * t285;
t283 = -t299 * t257 + (-t259 * t281 + t277 * t285) * t272;
t260 = t273 * t299 - t281 * t297;
t251 = (t259 * t277 + t281 * t285) * t272;
t250 = (t258 * t277 - t266 * t281) * t272;
t249 = t251 * t280 + t275 * t290;
t248 = t264 * t259 + t263 * t285;
t247 = t264 * t258 - t263 * t266;
t246 = t262 * t259 - t261 * t285;
t245 = t262 * t258 + t261 * t266;
t243 = t299 * t256 + (t258 * t281 + t266 * t277) * t272;
t240 = t256 * t294 + t263 * t258 + t264 * t266;
t233 = t248 * t280 + t264 * t296;
t232 = t246 * t280 + t262 * t296;
t231 = t260 * t275 + t280 * t283;
t230 = t260 * t280 - t275 * t283;
t229 = t275 * t287 + t280 * t284;
t228 = t275 * t284 - t280 * t287;
t224 = t229 * t279 - t240 * t274;
t223 = -t229 * t274 - t240 * t279;
t1 = [t303, t233 * t279 + t247 * t274, t240 * t291 + t274 * t284, 0, -t228 * t279, t223; t224, t232 * t279 + t245 * t274, -t238 * t274 + t286 * t291, 0, t225 * t279, t304; 0, t249 * t279 + t250 * t274, t243 * t291 + t274 * t283, 0, t230 * t279, -t231 * t274 - t243 * t279; -t304, -t233 * t274 + t247 * t279, -t240 * t292 + t279 * t284, 0, t228 * t274, -t224; t223, -t232 * t274 + t245 * t279, -t238 * t279 - t286 * t292, 0, -t225 * t274, t303; 0, -t249 * t274 + t250 * t279, -t243 * t292 + t279 * t283, 0, -t230 * t274, -t231 * t279 + t243 * t274; t225, t248 * t275 - t264 * t295, t240 * t275, 0, t229, 0; t228, t246 * t275 - t262 * t295, t286 * t275, 0, -t227, 0; 0, t251 * t275 - t280 * t290, t243 * t275, 0, t231, 0;];
JR_rot  = t1;
