% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRPR9_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:30
% EndTime: 2019-02-26 21:05:30
% DurationCPUTime: 0.17s
% Computational Cost: add. (82->40), mult. (247->86), div. (0->0), fcn. (270->12), ass. (0->42)
t274 = sin(pkin(7));
t278 = cos(pkin(6));
t300 = t274 * t278;
t275 = sin(pkin(6));
t280 = sin(qJ(1));
t299 = t275 * t280;
t282 = cos(qJ(1));
t298 = t275 * t282;
t277 = cos(pkin(7));
t279 = sin(qJ(3));
t297 = t277 * t279;
t296 = t278 * t282;
t273 = sin(pkin(12));
t295 = t280 * t273;
t276 = cos(pkin(12));
t294 = t280 * t276;
t293 = qJD(1) * t275;
t272 = qJ(4) + pkin(13);
t270 = sin(t272);
t292 = qJD(3) * t270;
t291 = t280 * t293;
t290 = t282 * t293;
t289 = t274 * t291;
t288 = t274 * t290;
t266 = t276 * t296 - t295;
t287 = t266 * t277 - t274 * t298;
t268 = -t273 * t282 - t278 * t294;
t286 = t268 * t277 + t274 * t299;
t267 = t273 * t296 + t294;
t269 = t276 * t282 - t278 * t295;
t281 = cos(qJ(3));
t285 = t267 * t281 + t287 * t279;
t284 = t269 * t281 + t286 * t279;
t283 = t279 * t300 + (t273 * t281 + t276 * t297) * t275;
t271 = cos(t272);
t265 = t269 * qJD(1);
t264 = t268 * qJD(1);
t263 = t267 * qJD(1);
t262 = t266 * qJD(1);
t261 = -t264 * t274 + t277 * t291;
t260 = t262 * t274 + t277 * t290;
t1 = [0, 0, t260, -t263 * t279 + (t262 * t277 - t288) * t281 + t284 * qJD(3), 0 (-t262 * t297 - t263 * t281 + t279 * t288) * t270 - t260 * t271 + (t284 * t271 + (-t268 * t274 + t277 * t299) * t270) * qJD(4) + (-t269 * t279 + t286 * t281) * t292; 0, 0, t261, t265 * t279 + (-t264 * t277 - t289) * t281 + t285 * qJD(3), 0 (t264 * t297 + t265 * t281 + t279 * t289) * t270 - t261 * t271 + (t285 * t271 + (-t266 * t274 - t277 * t298) * t270) * qJD(4) + (-t267 * t279 + t287 * t281) * t292; 0, 0, 0, t283 * qJD(3), 0 (t283 * t271 + (-t275 * t276 * t274 + t277 * t278) * t270) * qJD(4) + (t281 * t300 + (t276 * t277 * t281 - t273 * t279) * t275) * t292;];
JgD_rot  = t1;
