% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRPRR11_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:32
% EndTime: 2019-02-26 20:54:32
% DurationCPUTime: 0.17s
% Computational Cost: add. (82->40), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->43)
t275 = sin(pkin(7));
t279 = cos(pkin(6));
t302 = t275 * t279;
t276 = sin(pkin(6));
t281 = sin(qJ(1));
t301 = t276 * t281;
t283 = cos(qJ(1));
t300 = t276 * t283;
t278 = cos(pkin(7));
t280 = sin(qJ(3));
t299 = t278 * t280;
t274 = sin(pkin(12));
t298 = t281 * t274;
t277 = cos(pkin(12));
t297 = t281 * t277;
t296 = t283 * t274;
t295 = t283 * t277;
t294 = qJD(1) * t276;
t273 = pkin(13) + qJ(5);
t271 = sin(t273);
t293 = qJD(3) * t271;
t292 = t281 * t294;
t291 = t283 * t294;
t290 = t275 * t292;
t289 = t275 * t291;
t267 = t279 * t295 - t298;
t288 = t267 * t278 - t275 * t300;
t269 = -t279 * t297 - t296;
t287 = t269 * t278 + t275 * t301;
t268 = t279 * t296 + t297;
t270 = -t279 * t298 + t295;
t282 = cos(qJ(3));
t286 = t268 * t282 + t280 * t288;
t285 = t270 * t282 + t280 * t287;
t284 = t280 * t302 + (t274 * t282 + t277 * t299) * t276;
t272 = cos(t273);
t266 = t270 * qJD(1);
t265 = t269 * qJD(1);
t264 = t268 * qJD(1);
t263 = t267 * qJD(1);
t262 = -t265 * t275 + t278 * t292;
t261 = t263 * t275 + t278 * t291;
t1 = [0, 0, t261, 0, -t264 * t280 + (t263 * t278 - t289) * t282 + t285 * qJD(3) (-t263 * t299 - t264 * t282 + t280 * t289) * t271 - t261 * t272 + (t285 * t272 + (-t269 * t275 + t278 * t301) * t271) * qJD(5) + (-t270 * t280 + t282 * t287) * t293; 0, 0, t262, 0, t266 * t280 + (-t265 * t278 - t290) * t282 + t286 * qJD(3) (t265 * t299 + t266 * t282 + t280 * t290) * t271 - t262 * t272 + (t286 * t272 + (-t267 * t275 - t278 * t300) * t271) * qJD(5) + (-t268 * t280 + t282 * t288) * t293; 0, 0, 0, 0, t284 * qJD(3) (t284 * t272 + (-t276 * t277 * t275 + t279 * t278) * t271) * qJD(5) + (t282 * t302 + (t277 * t278 * t282 - t274 * t280) * t276) * t293;];
JgD_rot  = t1;
