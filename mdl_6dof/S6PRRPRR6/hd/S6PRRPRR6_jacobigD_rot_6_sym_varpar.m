% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR6_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:07
% EndTime: 2019-02-26 20:07:07
% DurationCPUTime: 0.18s
% Computational Cost: add. (82->41), mult. (240->92), div. (0->0), fcn. (263->12), ass. (0->41)
t266 = pkin(13) + qJ(5);
t265 = cos(t266);
t268 = sin(pkin(7));
t296 = t265 * t268;
t269 = sin(pkin(6));
t295 = t268 * t269;
t274 = sin(qJ(2));
t294 = t268 * t274;
t271 = cos(pkin(7));
t293 = t269 * t271;
t273 = sin(qJ(3));
t292 = t271 * t273;
t275 = cos(qJ(3));
t291 = t271 * t275;
t272 = cos(pkin(6));
t290 = t272 * t274;
t276 = cos(qJ(2));
t289 = t272 * t276;
t288 = t273 * t274;
t287 = t273 * t276;
t286 = t274 * t275;
t285 = t275 * t276;
t264 = sin(t266);
t284 = qJD(3) * t264;
t283 = qJD(3) * t268;
t267 = sin(pkin(12));
t270 = cos(pkin(12));
t260 = -t267 * t274 + t270 * t289;
t282 = t260 * t271 - t270 * t295;
t262 = -t267 * t289 - t270 * t274;
t281 = t262 * t271 + t267 * t295;
t261 = t267 * t276 + t270 * t290;
t280 = t267 * t290 - t270 * t276;
t279 = t271 * t287 + t286;
t278 = t261 * t275 + t282 * t273;
t277 = t281 * t273 - t275 * t280;
t259 = t280 * qJD(2);
t258 = t262 * qJD(2);
t257 = t261 * qJD(2);
t256 = t260 * qJD(2);
t1 = [0, 0, -t259 * t268, 0, t277 * qJD(3) + t258 * t273 - t259 * t291 (t258 * t275 + t259 * t292) * t264 + t259 * t296 + (t277 * t265 + (-t262 * t268 + t267 * t293) * t264) * qJD(5) + (t273 * t280 + t281 * t275) * t284; 0, 0, t257 * t268, 0, t278 * qJD(3) + t256 * t273 + t257 * t291 (t256 * t275 - t257 * t292) * t264 - t257 * t296 + (t278 * t265 + (-t260 * t268 - t270 * t293) * t264) * qJD(5) + (-t261 * t273 + t282 * t275) * t284; 0, 0, t269 * qJD(2) * t294, 0, t272 * t273 * t283 + (t279 * qJD(3) + (t271 * t286 + t287) * qJD(2)) * t269 (t275 * t264 * t283 + (t264 * t271 + t273 * t296) * qJD(5)) * t272 + ((-t268 * t276 * t264 + t279 * t265) * qJD(5) + (t271 * t285 - t288) * t284 + ((-t271 * t288 + t285) * t264 - t265 * t294) * qJD(2)) * t269;];
JgD_rot  = t1;
