% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRR11_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:20:26
% EndTime: 2019-02-26 21:20:27
% DurationCPUTime: 0.15s
% Computational Cost: add. (118->39), mult. (402->86), div. (0->0), fcn. (444->12), ass. (0->44)
t286 = sin(pkin(7));
t290 = cos(pkin(6));
t314 = t286 * t290;
t287 = sin(pkin(6));
t293 = sin(qJ(1));
t313 = t287 * t293;
t296 = cos(qJ(1));
t312 = t287 * t296;
t289 = cos(pkin(7));
t292 = sin(qJ(3));
t311 = t289 * t292;
t310 = t290 * t296;
t285 = sin(pkin(13));
t309 = t293 * t285;
t288 = cos(pkin(13));
t308 = t293 * t288;
t307 = qJD(1) * t287;
t291 = sin(qJ(4));
t306 = qJD(3) * t291;
t305 = t293 * t307;
t304 = t296 * t307;
t303 = t286 * t305;
t302 = t286 * t304;
t281 = t288 * t310 - t309;
t301 = t281 * t289 - t286 * t312;
t283 = -t285 * t296 - t290 * t308;
t300 = t283 * t289 + t286 * t313;
t282 = t285 * t310 + t308;
t284 = t288 * t296 - t290 * t309;
t295 = cos(qJ(3));
t299 = t282 * t295 + t301 * t292;
t298 = t284 * t295 + t300 * t292;
t297 = t292 * t314 + (t285 * t295 + t288 * t311) * t287;
t294 = cos(qJ(4));
t280 = t284 * qJD(1);
t279 = t283 * qJD(1);
t278 = t282 * qJD(1);
t277 = t281 * qJD(1);
t276 = -t279 * t286 + t289 * t305;
t275 = t277 * t286 + t289 * t304;
t274 = (t297 * t294 + (-t286 * t287 * t288 + t289 * t290) * t291) * qJD(4) + (t295 * t314 + (t288 * t289 * t295 - t285 * t292) * t287) * t306;
t273 = (t279 * t311 + t280 * t295 + t292 * t303) * t291 - t276 * t294 + (t299 * t294 + (-t281 * t286 - t289 * t312) * t291) * qJD(4) + (-t282 * t292 + t301 * t295) * t306;
t272 = (-t277 * t311 - t278 * t295 + t292 * t302) * t291 - t275 * t294 + (t298 * t294 + (-t283 * t286 + t289 * t313) * t291) * qJD(4) + (-t284 * t292 + t300 * t295) * t306;
t1 = [0, 0, t275, -t278 * t292 + (t277 * t289 - t302) * t295 + t298 * qJD(3), t272, t272; 0, 0, t276, t280 * t292 + (-t279 * t289 - t303) * t295 + t299 * qJD(3), t273, t273; 0, 0, 0, t297 * qJD(3), t274, t274;];
JgD_rot  = t1;
