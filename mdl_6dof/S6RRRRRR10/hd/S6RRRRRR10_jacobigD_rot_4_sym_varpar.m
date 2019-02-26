% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR10_jacobigD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobigD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:57
% EndTime: 2019-02-26 22:52:57
% DurationCPUTime: 0.12s
% Computational Cost: add. (49->27), mult. (173->62), div. (0->0), fcn. (177->12), ass. (0->34)
t270 = sin(pkin(7));
t271 = sin(pkin(6));
t297 = t271 * t270;
t273 = cos(pkin(7));
t278 = cos(qJ(3));
t296 = t273 * t278;
t275 = sin(qJ(3));
t279 = cos(qJ(2));
t295 = t275 * t279;
t276 = sin(qJ(2));
t294 = t276 * t278;
t277 = sin(qJ(1));
t293 = t277 * t276;
t292 = t277 * t279;
t280 = cos(qJ(1));
t291 = t280 * t276;
t290 = t280 * t279;
t289 = qJD(1) * t271;
t269 = sin(pkin(8));
t288 = qJD(3) * t269;
t287 = t277 * t289;
t286 = t280 * t289;
t285 = t270 * t278 * t289;
t274 = cos(pkin(6));
t284 = t274 * t290 - t293;
t283 = -t274 * t292 - t291;
t282 = t274 * t291 + t292;
t281 = t274 * t293 - t290;
t272 = cos(pkin(8));
t268 = t283 * qJD(1) - t282 * qJD(2);
t267 = -t284 * qJD(1) + t281 * qJD(2);
t266 = -t268 * t270 + t273 * t287;
t265 = -t267 * t270 + t273 * t286;
t1 = [0, t286, t265 -(-(-t282 * qJD(1) + t283 * qJD(2)) * t275 + t267 * t296 + t280 * t285) * t269 + t265 * t272 - (t281 * t278 + (-t283 * t273 - t277 * t297) * t275) * t288, 0, 0; 0, t287, t266 -(-(-t281 * qJD(1) + t284 * qJD(2)) * t275 + t268 * t296 + t277 * t285) * t269 + t266 * t272 - (-t282 * t278 + (-t284 * t273 + t280 * t297) * t275) * t288, 0, 0; 0, 0, qJD(2) * t276 * t297, t274 * t270 * t275 * t288 + (-(-t273 * t295 - t294) * t288 + (-(-t273 * t294 - t295) * t269 + t276 * t270 * t272) * qJD(2)) * t271, 0, 0;];
JgD_rot  = t1;
