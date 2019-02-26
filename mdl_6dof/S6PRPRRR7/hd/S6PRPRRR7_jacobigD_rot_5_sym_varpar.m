% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR7_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobigD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:31
% EndTime: 2019-02-26 19:57:31
% DurationCPUTime: 0.21s
% Computational Cost: add. (68->41), mult. (251->100), div. (0->0), fcn. (276->14), ass. (0->41)
t284 = sin(pkin(14));
t292 = cos(pkin(7));
t313 = t284 * t292;
t286 = sin(pkin(8));
t287 = sin(pkin(7));
t312 = t287 * t286;
t288 = sin(pkin(6));
t311 = t287 * t288;
t291 = cos(pkin(8));
t310 = t287 * t291;
t293 = cos(pkin(6));
t309 = t287 * t293;
t295 = sin(qJ(2));
t308 = t287 * t295;
t307 = t288 * t292;
t289 = cos(pkin(14));
t306 = t289 * t292;
t305 = t292 * t295;
t297 = cos(qJ(2));
t304 = t292 * t297;
t303 = t293 * t295;
t302 = t293 * t297;
t301 = qJD(2) * t288;
t285 = sin(pkin(13));
t290 = cos(pkin(13));
t280 = -t285 * t295 + t290 * t302;
t300 = t280 * t292 - t290 * t311;
t282 = -t285 * t302 - t290 * t295;
t299 = t282 * t292 + t285 * t311;
t281 = t285 * t297 + t290 * t303;
t298 = t285 * t303 - t290 * t297;
t296 = cos(qJ(4));
t294 = sin(qJ(4));
t279 = t298 * qJD(2);
t278 = t282 * qJD(2);
t277 = t281 * qJD(2);
t276 = t280 * qJD(2);
t275 = (-t284 * t297 - t289 * t305) * t301;
t274 = -t278 * t284 + t279 * t306;
t273 = -t276 * t284 - t277 * t306;
t1 = [0, 0, 0, -t274 * t286 - t279 * t310 (t278 * t289 + t279 * t313) * t294 + (-t274 * t291 + t279 * t312) * t296 + ((t284 * t299 - t289 * t298) * t296 + ((t284 * t298 + t289 * t299) * t291 + (-t282 * t287 + t285 * t307) * t286) * t294) * qJD(4), 0; 0, 0, 0, -t273 * t286 + t277 * t310 (t276 * t289 - t277 * t313) * t294 + (-t273 * t291 - t277 * t312) * t296 + ((t281 * t289 + t284 * t300) * t296 + ((-t281 * t284 + t289 * t300) * t291 + (-t280 * t287 - t290 * t307) * t286) * t294) * qJD(4), 0; 0, 0, 0, t291 * t301 * t308 - t275 * t286, -t275 * t291 * t296 + ((t288 * t295 * t289 + (t288 * t304 + t309) * t284) * t296 + ((t289 * t309 + (-t284 * t295 + t289 * t304) * t288) * t291 + (t293 * t292 - t297 * t311) * t286) * t294) * qJD(4) + ((-t284 * t305 + t289 * t297) * t294 - t286 * t296 * t308) * t301, 0;];
JgD_rot  = t1;
