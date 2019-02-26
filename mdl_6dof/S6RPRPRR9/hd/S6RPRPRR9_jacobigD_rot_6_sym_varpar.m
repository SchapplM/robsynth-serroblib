% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR9
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
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRPRR9_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:25
% EndTime: 2019-02-26 20:53:25
% DurationCPUTime: 0.22s
% Computational Cost: add. (110->49), mult. (369->102), div. (0->0), fcn. (413->14), ass. (0->46)
t292 = sin(pkin(13));
t296 = cos(pkin(13));
t301 = sin(qJ(3));
t304 = cos(qJ(3));
t307 = t301 * t292 - t304 * t296;
t318 = t307 * qJD(3);
t295 = sin(pkin(6));
t302 = sin(qJ(1));
t317 = t295 * t302;
t305 = cos(qJ(1));
t316 = t295 * t305;
t293 = sin(pkin(12));
t315 = t302 * t293;
t297 = cos(pkin(12));
t314 = t302 * t297;
t313 = t305 * t293;
t312 = t305 * t297;
t311 = qJD(1) * t302;
t310 = qJD(1) * t305;
t298 = cos(pkin(7));
t309 = qJD(1) * t295 * t298;
t308 = t304 * t292 + t301 * t296;
t299 = cos(pkin(6));
t284 = t299 * t312 - t315;
t286 = -t299 * t314 - t313;
t285 = t299 * t313 + t314;
t287 = -t299 * t315 + t312;
t289 = t308 * qJD(3);
t303 = cos(qJ(5));
t300 = sin(qJ(5));
t294 = sin(pkin(7));
t283 = t287 * qJD(1);
t282 = t286 * qJD(1);
t281 = t285 * qJD(1);
t280 = t284 * qJD(1);
t279 = t308 * t298;
t278 = t307 * t298;
t277 = t308 * t294;
t276 = t307 * t294;
t275 = t298 * t318;
t274 = t298 * t289;
t273 = t294 * t318;
t272 = t294 * t289;
t271 = -t282 * t294 + t302 * t309;
t270 = t280 * t294 + t305 * t309;
t1 = [0, 0, t270, 0, t286 * t274 - t280 * t278 - t281 * t308 - t287 * t318 + (t272 * t302 + t276 * t310) * t295 (-t286 * t275 - t280 * t279 + t281 * t307 - t287 * t289 + (-t273 * t302 + t277 * t310) * t295) * t300 - t270 * t303 + ((t277 * t317 + t286 * t279 - t287 * t307) * t303 + (-t286 * t294 + t298 * t317) * t300) * qJD(5); 0, 0, t271, 0, t284 * t274 + t282 * t278 + t283 * t308 - t285 * t318 + (-t272 * t305 + t276 * t311) * t295 (-t284 * t275 + t282 * t279 - t283 * t307 - t285 * t289 + (t273 * t305 + t277 * t311) * t295) * t300 - t271 * t303 + ((-t277 * t316 + t284 * t279 - t285 * t307) * t303 + (-t284 * t294 - t298 * t316) * t300) * qJD(5); 0, 0, 0, 0, t299 * t272 + (t274 * t297 - t293 * t318) * t295 (-t299 * t273 + (-t275 * t297 - t289 * t293) * t295) * t300 + ((t299 * t277 + (t279 * t297 - t293 * t307) * t295) * t303 + (-t295 * t297 * t294 + t299 * t298) * t300) * qJD(5);];
JgD_rot  = t1;
