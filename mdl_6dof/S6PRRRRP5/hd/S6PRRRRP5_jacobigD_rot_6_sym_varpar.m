% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRP5_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:29
% EndTime: 2019-02-26 20:17:29
% DurationCPUTime: 0.16s
% Computational Cost: add. (70->40), mult. (240->91), div. (0->0), fcn. (263->12), ass. (0->38)
t266 = sin(pkin(7));
t267 = sin(pkin(6));
t294 = t266 * t267;
t274 = cos(qJ(4));
t293 = t266 * t274;
t269 = cos(pkin(7));
t292 = t267 * t269;
t272 = sin(qJ(3));
t291 = t269 * t272;
t275 = cos(qJ(3));
t290 = t269 * t275;
t270 = cos(pkin(6));
t273 = sin(qJ(2));
t289 = t270 * t273;
t276 = cos(qJ(2));
t288 = t270 * t276;
t287 = t272 * t273;
t286 = t272 * t276;
t285 = t273 * t275;
t284 = t275 * t276;
t271 = sin(qJ(4));
t283 = qJD(3) * t271;
t265 = sin(pkin(12));
t268 = cos(pkin(12));
t261 = -t265 * t273 + t268 * t288;
t282 = t261 * t269 - t268 * t294;
t263 = -t265 * t288 - t268 * t273;
t281 = t263 * t269 + t265 * t294;
t262 = t265 * t276 + t268 * t289;
t280 = t265 * t289 - t268 * t276;
t279 = t269 * t286 + t285;
t278 = t262 * t275 + t272 * t282;
t277 = t272 * t281 - t275 * t280;
t260 = t280 * qJD(2);
t259 = t263 * qJD(2);
t258 = t262 * qJD(2);
t257 = t261 * qJD(2);
t1 = [0, 0, -t260 * t266, qJD(3) * t277 + t259 * t272 - t260 * t290 (t259 * t275 + t260 * t291) * t271 + t260 * t293 + (t277 * t274 + (-t263 * t266 + t265 * t292) * t271) * qJD(4) + (t272 * t280 + t275 * t281) * t283, 0; 0, 0, t258 * t266, qJD(3) * t278 + t257 * t272 + t258 * t290 (t257 * t275 - t258 * t291) * t271 - t258 * t293 + (t278 * t274 + (-t261 * t266 - t268 * t292) * t271) * qJD(4) + (-t262 * t272 + t275 * t282) * t283, 0; 0, 0, qJD(2) * t273 * t294, t270 * t266 * qJD(3) * t272 + (t279 * qJD(3) + (t269 * t285 + t286) * qJD(2)) * t267 (t266 * t275 * t283 + (t269 * t271 + t272 * t293) * qJD(4)) * t270 + ((-t266 * t276 * t271 + t274 * t279) * qJD(4) + (t269 * t284 - t287) * t283 + ((-t269 * t287 + t284) * t271 - t273 * t293) * qJD(2)) * t267, 0;];
JgD_rot  = t1;
