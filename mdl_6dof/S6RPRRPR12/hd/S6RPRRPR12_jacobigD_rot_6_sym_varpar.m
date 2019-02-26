% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRPR12_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:17
% EndTime: 2019-02-26 21:07:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (71->39), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->42)
t270 = sin(pkin(7));
t274 = cos(pkin(6));
t299 = t270 * t274;
t271 = sin(pkin(6));
t277 = sin(qJ(1));
t298 = t271 * t277;
t280 = cos(qJ(1));
t297 = t271 * t280;
t273 = cos(pkin(7));
t276 = sin(qJ(3));
t296 = t273 * t276;
t269 = sin(pkin(12));
t295 = t277 * t269;
t272 = cos(pkin(12));
t294 = t277 * t272;
t293 = t280 * t269;
t292 = t280 * t272;
t291 = qJD(1) * t271;
t278 = cos(qJ(4));
t290 = qJD(3) * t278;
t289 = t277 * t291;
t288 = t280 * t291;
t287 = t270 * t289;
t286 = t270 * t288;
t265 = t274 * t292 - t295;
t285 = t265 * t273 - t270 * t297;
t267 = -t274 * t294 - t293;
t284 = t267 * t273 + t270 * t298;
t266 = t274 * t293 + t294;
t268 = -t274 * t295 + t292;
t279 = cos(qJ(3));
t283 = t266 * t279 + t285 * t276;
t282 = t268 * t279 + t284 * t276;
t281 = t276 * t299 + (t269 * t279 + t272 * t296) * t271;
t275 = sin(qJ(4));
t264 = t268 * qJD(1);
t263 = t267 * qJD(1);
t262 = t266 * qJD(1);
t261 = t265 * qJD(1);
t260 = -t263 * t270 + t273 * t289;
t259 = t261 * t270 + t273 * t288;
t1 = [0, 0, t259, -t262 * t276 + (t261 * t273 - t286) * t279 + t282 * qJD(3), 0 (-t261 * t296 - t262 * t279 + t276 * t286) * t278 + t259 * t275 + (-t282 * t275 + (-t267 * t270 + t273 * t298) * t278) * qJD(4) + (-t268 * t276 + t284 * t279) * t290; 0, 0, t260, t264 * t276 + (-t263 * t273 - t287) * t279 + t283 * qJD(3), 0 (t263 * t296 + t264 * t279 + t276 * t287) * t278 + t260 * t275 + (-t283 * t275 + (-t265 * t270 - t273 * t297) * t278) * qJD(4) + (-t266 * t276 + t285 * t279) * t290; 0, 0, 0, t281 * qJD(3), 0 (-t281 * t275 + (-t271 * t272 * t270 + t274 * t273) * t278) * qJD(4) + (t279 * t299 + (t272 * t273 * t279 - t269 * t276) * t271) * t290;];
JgD_rot  = t1;
