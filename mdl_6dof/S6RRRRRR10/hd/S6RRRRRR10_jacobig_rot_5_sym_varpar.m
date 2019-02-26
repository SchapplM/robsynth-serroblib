% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR10_jacobig_rot_5_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobig_rot_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobig_rot_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:19
% EndTime: 2018-11-23 11:27:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (208->45), mult. (212->65), div. (0->0), fcn. (227->25), ass. (0->50)
t304 = sin(pkin(6));
t311 = sin(qJ(1));
t316 = t311 * t304;
t314 = cos(qJ(1));
t315 = t314 * t304;
t313 = cos(qJ(2));
t312 = cos(qJ(3));
t310 = sin(qJ(2));
t309 = sin(qJ(3));
t308 = sin(qJ(4));
t307 = cos(pkin(6));
t306 = cos(pkin(7));
t305 = cos(pkin(8));
t303 = sin(pkin(7));
t302 = sin(pkin(8));
t301 = pkin(6) - qJ(2);
t300 = pkin(6) + qJ(2);
t299 = pkin(7) - qJ(3);
t298 = pkin(7) + qJ(3);
t297 = pkin(8) - qJ(4);
t296 = pkin(8) + qJ(4);
t295 = cos(t300);
t294 = cos(t298);
t293 = sin(t301);
t292 = sin(t299);
t291 = cos(t301) / 0.2e1;
t290 = cos(t299) / 0.2e1;
t289 = sin(t300) / 0.2e1;
t288 = sin(t298) / 0.2e1;
t287 = t291 - t295 / 0.2e1;
t286 = t291 + t295 / 0.2e1;
t285 = t290 - t294 / 0.2e1;
t284 = t290 + t294 / 0.2e1;
t283 = cos(t297) / 0.2e1 + cos(t296) / 0.2e1;
t282 = t289 - t293 / 0.2e1;
t281 = t289 + t293 / 0.2e1;
t280 = t288 - t292 / 0.2e1;
t279 = t288 + t292 / 0.2e1;
t278 = sin(t296) / 0.2e1 + sin(t297) / 0.2e1;
t277 = -t311 * t282 + t314 * t313;
t276 = -t311 * t286 - t314 * t310;
t275 = t314 * t282 + t311 * t313;
t274 = t314 * t286 - t311 * t310;
t273 = -t281 * t303 + t307 * t306;
t272 = -t276 * t303 + t306 * t316;
t271 = -t274 * t303 - t306 * t315;
t270 = t307 * t279 + t281 * t284 - t287 * t309;
t269 = t276 * t284 - t277 * t309 + t279 * t316;
t268 = t274 * t284 - t275 * t309 - t279 * t315;
t1 = [0, t316, t272, -t269 * t302 + t272 * t305 (t276 * t280 + t277 * t312 + t285 * t316) * t308 - t269 * t283 - t272 * t278, 0; 0, -t315, t271, -t268 * t302 + t271 * t305 (t274 * t280 + t275 * t312 - t285 * t315) * t308 - t268 * t283 - t271 * t278, 0; 1, t307, t273, -t270 * t302 + t273 * t305 (t281 * t280 + t307 * t285 + t287 * t312) * t308 - t270 * t283 - t273 * t278, 0;];
Jg_rot  = t1;
