% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR11
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
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRPR11_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:45
% EndTime: 2019-02-26 21:06:45
% DurationCPUTime: 0.14s
% Computational Cost: add. (71->39), mult. (247->86), div. (0->0), fcn. (270->12), ass. (0->41)
t275 = sin(pkin(7));
t279 = cos(pkin(6));
t303 = t275 * t279;
t276 = sin(pkin(6));
t282 = sin(qJ(1));
t302 = t276 * t282;
t285 = cos(qJ(1));
t301 = t276 * t285;
t278 = cos(pkin(7));
t281 = sin(qJ(3));
t300 = t278 * t281;
t299 = t279 * t285;
t274 = sin(pkin(12));
t298 = t282 * t274;
t277 = cos(pkin(12));
t297 = t282 * t277;
t296 = qJD(1) * t276;
t280 = sin(qJ(4));
t295 = qJD(3) * t280;
t294 = t282 * t296;
t293 = t285 * t296;
t292 = t275 * t294;
t291 = t275 * t293;
t270 = t277 * t299 - t298;
t290 = t270 * t278 - t275 * t301;
t272 = -t274 * t285 - t279 * t297;
t289 = t272 * t278 + t275 * t302;
t271 = t274 * t299 + t297;
t273 = t277 * t285 - t279 * t298;
t284 = cos(qJ(3));
t288 = t271 * t284 + t290 * t281;
t287 = t273 * t284 + t289 * t281;
t286 = t281 * t303 + (t274 * t284 + t277 * t300) * t276;
t283 = cos(qJ(4));
t269 = t273 * qJD(1);
t268 = t272 * qJD(1);
t267 = t271 * qJD(1);
t266 = t270 * qJD(1);
t265 = -t268 * t275 + t278 * t294;
t264 = t266 * t275 + t278 * t293;
t1 = [0, 0, t264, -t267 * t281 + (t266 * t278 - t291) * t284 + t287 * qJD(3), 0 (-t266 * t300 - t267 * t284 + t281 * t291) * t280 - t264 * t283 + (t287 * t283 + (-t272 * t275 + t278 * t302) * t280) * qJD(4) + (-t273 * t281 + t289 * t284) * t295; 0, 0, t265, t269 * t281 + (-t268 * t278 - t292) * t284 + t288 * qJD(3), 0 (t268 * t300 + t269 * t284 + t281 * t292) * t280 - t265 * t283 + (t288 * t283 + (-t270 * t275 - t278 * t301) * t280) * qJD(4) + (-t271 * t281 + t290 * t284) * t295; 0, 0, 0, t286 * qJD(3), 0 (t286 * t283 + (-t275 * t276 * t277 + t278 * t279) * t280) * qJD(4) + (t284 * t303 + (t277 * t278 * t284 - t274 * t281) * t276) * t295;];
JgD_rot  = t1;
