% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
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
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:25
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR14_jacobig_rot_6_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobig_rot_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobig_rot_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:25:35
% EndTime: 2019-01-03 10:25:35
% DurationCPUTime: 0.09s
% Computational Cost: add. (102->33), mult. (296->69), div. (0->0), fcn. (409->16), ass. (0->46)
t274 = sin(pkin(7));
t279 = cos(pkin(6));
t300 = t274 * t279;
t278 = cos(pkin(7));
t286 = cos(qJ(2));
t299 = t278 * t286;
t275 = sin(pkin(6));
t283 = sin(qJ(1));
t298 = t283 * t275;
t282 = sin(qJ(2));
t297 = t283 * t282;
t296 = t283 * t286;
t287 = cos(qJ(1));
t295 = t287 * t275;
t294 = t287 * t282;
t293 = t287 * t286;
t269 = t279 * t294 + t296;
t272 = sin(pkin(14));
t276 = cos(pkin(14));
t268 = t279 * t293 - t297;
t289 = t268 * t278 - t274 * t295;
t259 = -t269 * t272 + t289 * t276;
t265 = -t268 * t274 - t278 * t295;
t273 = sin(pkin(8));
t277 = cos(pkin(8));
t292 = t259 * t277 + t265 * t273;
t271 = -t279 * t297 + t293;
t270 = -t279 * t296 - t294;
t288 = t270 * t278 + t274 * t298;
t261 = -t271 * t272 + t288 * t276;
t266 = -t270 * t274 + t278 * t298;
t291 = t261 * t277 + t266 * t273;
t263 = t276 * t300 + (-t272 * t282 + t276 * t299) * t275;
t267 = -t275 * t286 * t274 + t279 * t278;
t290 = t263 * t277 + t267 * t273;
t285 = cos(qJ(4));
t284 = cos(qJ(5));
t281 = sin(qJ(4));
t280 = sin(qJ(5));
t264 = t275 * t282 * t276 + (t275 * t299 + t300) * t272;
t262 = t271 * t276 + t288 * t272;
t260 = t269 * t276 + t289 * t272;
t258 = -t263 * t273 + t267 * t277;
t257 = -t261 * t273 + t266 * t277;
t256 = -t259 * t273 + t265 * t277;
t1 = [0, t298, 0, t257, t262 * t281 - t291 * t285 (t262 * t285 + t291 * t281) * t280 - t257 * t284; 0, -t295, 0, t256, t260 * t281 - t292 * t285 (t260 * t285 + t292 * t281) * t280 - t256 * t284; 1, t279, 0, t258, t264 * t281 - t290 * t285 (t264 * t285 + t290 * t281) * t280 - t258 * t284;];
Jg_rot  = t1;
