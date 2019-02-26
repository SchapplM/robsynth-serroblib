% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRP11_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:32
% EndTime: 2019-02-26 21:13:33
% DurationCPUTime: 0.17s
% Computational Cost: add. (71->39), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->42)
t273 = sin(pkin(7));
t277 = cos(pkin(6));
t302 = t273 * t277;
t274 = sin(pkin(6));
t280 = sin(qJ(1));
t301 = t274 * t280;
t283 = cos(qJ(1));
t300 = t274 * t283;
t276 = cos(pkin(7));
t279 = sin(qJ(3));
t299 = t276 * t279;
t272 = sin(pkin(12));
t298 = t280 * t272;
t275 = cos(pkin(12));
t297 = t280 * t275;
t296 = t283 * t272;
t295 = t283 * t275;
t294 = qJD(1) * t274;
t278 = sin(qJ(4));
t293 = qJD(3) * t278;
t292 = t280 * t294;
t291 = t283 * t294;
t290 = t273 * t292;
t289 = t273 * t291;
t268 = t277 * t295 - t298;
t288 = t268 * t276 - t273 * t300;
t270 = -t277 * t297 - t296;
t287 = t270 * t276 + t273 * t301;
t269 = t277 * t296 + t297;
t271 = -t277 * t298 + t295;
t282 = cos(qJ(3));
t286 = t269 * t282 + t288 * t279;
t285 = t271 * t282 + t287 * t279;
t284 = t279 * t302 + (t272 * t282 + t275 * t299) * t274;
t281 = cos(qJ(4));
t267 = t271 * qJD(1);
t266 = t270 * qJD(1);
t265 = t269 * qJD(1);
t264 = t268 * qJD(1);
t263 = -t266 * t273 + t276 * t292;
t262 = t264 * t273 + t276 * t291;
t1 = [0, 0, t262, -t265 * t279 + (t264 * t276 - t289) * t282 + t285 * qJD(3) (-t264 * t299 - t265 * t282 + t279 * t289) * t278 - t262 * t281 + (t285 * t281 + (-t270 * t273 + t276 * t301) * t278) * qJD(4) + (-t271 * t279 + t287 * t282) * t293, 0; 0, 0, t263, t267 * t279 + (-t266 * t276 - t290) * t282 + t286 * qJD(3) (t266 * t299 + t267 * t282 + t279 * t290) * t278 - t263 * t281 + (t286 * t281 + (-t268 * t273 - t276 * t300) * t278) * qJD(4) + (-t269 * t279 + t288 * t282) * t293, 0; 0, 0, 0, t284 * qJD(3) (t284 * t281 + (-t274 * t275 * t273 + t277 * t276) * t278) * qJD(4) + (t282 * t302 + (t275 * t276 * t282 - t272 * t279) * t274) * t293, 0;];
JgD_rot  = t1;
