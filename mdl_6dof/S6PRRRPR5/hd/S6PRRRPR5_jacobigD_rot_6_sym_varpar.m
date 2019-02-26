% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR5_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:55
% EndTime: 2019-02-26 20:12:55
% DurationCPUTime: 0.17s
% Computational Cost: add. (82->41), mult. (240->92), div. (0->0), fcn. (263->12), ass. (0->41)
t267 = qJ(4) + pkin(13);
t266 = cos(t267);
t269 = sin(pkin(7));
t297 = t269 * t266;
t270 = sin(pkin(6));
t296 = t269 * t270;
t275 = sin(qJ(2));
t295 = t269 * t275;
t272 = cos(pkin(7));
t294 = t270 * t272;
t274 = sin(qJ(3));
t293 = t272 * t274;
t276 = cos(qJ(3));
t292 = t272 * t276;
t273 = cos(pkin(6));
t291 = t273 * t275;
t277 = cos(qJ(2));
t290 = t273 * t277;
t289 = t274 * t275;
t288 = t274 * t277;
t287 = t275 * t276;
t286 = t276 * t277;
t265 = sin(t267);
t285 = qJD(3) * t265;
t284 = qJD(3) * t269;
t268 = sin(pkin(12));
t271 = cos(pkin(12));
t261 = -t268 * t275 + t271 * t290;
t283 = t261 * t272 - t271 * t296;
t263 = -t268 * t290 - t271 * t275;
t282 = t263 * t272 + t268 * t296;
t262 = t268 * t277 + t271 * t291;
t281 = t268 * t291 - t271 * t277;
t280 = t272 * t288 + t287;
t279 = t262 * t276 + t283 * t274;
t278 = t282 * t274 - t276 * t281;
t260 = t281 * qJD(2);
t259 = t263 * qJD(2);
t258 = t262 * qJD(2);
t257 = t261 * qJD(2);
t1 = [0, 0, -t260 * t269, t278 * qJD(3) + t259 * t274 - t260 * t292, 0 (t259 * t276 + t260 * t293) * t265 + t260 * t297 + (t278 * t266 + (-t263 * t269 + t268 * t294) * t265) * qJD(4) + (t274 * t281 + t282 * t276) * t285; 0, 0, t258 * t269, t279 * qJD(3) + t257 * t274 + t258 * t292, 0 (t257 * t276 - t258 * t293) * t265 - t258 * t297 + (t279 * t266 + (-t261 * t269 - t271 * t294) * t265) * qJD(4) + (-t262 * t274 + t283 * t276) * t285; 0, 0, t270 * qJD(2) * t295, t273 * t274 * t284 + (t280 * qJD(3) + (t272 * t287 + t288) * qJD(2)) * t270, 0 (t276 * t265 * t284 + (t265 * t272 + t274 * t297) * qJD(4)) * t273 + ((-t269 * t277 * t265 + t280 * t266) * qJD(4) + (t272 * t286 - t289) * t285 + ((-t272 * t289 + t286) * t265 - t266 * t295) * qJD(2)) * t270;];
JgD_rot  = t1;
