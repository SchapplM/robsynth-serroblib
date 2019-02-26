% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR5_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:11
% EndTime: 2019-02-26 20:21:11
% DurationCPUTime: 0.18s
% Computational Cost: add. (116->40), mult. (392->91), div. (0->0), fcn. (434->12), ass. (0->41)
t285 = sin(pkin(7));
t286 = sin(pkin(6));
t313 = t285 * t286;
t293 = cos(qJ(4));
t312 = t285 * t293;
t288 = cos(pkin(7));
t311 = t286 * t288;
t291 = sin(qJ(3));
t310 = t288 * t291;
t294 = cos(qJ(3));
t309 = t288 * t294;
t289 = cos(pkin(6));
t292 = sin(qJ(2));
t308 = t289 * t292;
t295 = cos(qJ(2));
t307 = t289 * t295;
t306 = t291 * t292;
t305 = t291 * t295;
t304 = t292 * t294;
t303 = t294 * t295;
t290 = sin(qJ(4));
t302 = qJD(3) * t290;
t284 = sin(pkin(13));
t287 = cos(pkin(13));
t280 = -t284 * t292 + t287 * t307;
t301 = t280 * t288 - t287 * t313;
t282 = -t284 * t307 - t287 * t292;
t300 = t282 * t288 + t284 * t313;
t281 = t284 * t295 + t287 * t308;
t299 = t284 * t308 - t287 * t295;
t298 = t288 * t305 + t304;
t297 = t281 * t294 + t301 * t291;
t296 = t300 * t291 - t294 * t299;
t279 = t299 * qJD(2);
t278 = t282 * qJD(2);
t277 = t281 * qJD(2);
t276 = t280 * qJD(2);
t275 = (t285 * t294 * t302 + (t288 * t290 + t291 * t312) * qJD(4)) * t289 + ((-t285 * t295 * t290 + t298 * t293) * qJD(4) + (t288 * t303 - t306) * t302 + ((-t288 * t306 + t303) * t290 - t292 * t312) * qJD(2)) * t286;
t274 = (t278 * t294 + t279 * t310) * t290 + t279 * t312 + (t296 * t293 + (-t282 * t285 + t284 * t311) * t290) * qJD(4) + (t291 * t299 + t300 * t294) * t302;
t273 = (t276 * t294 - t277 * t310) * t290 - t277 * t312 + (t297 * t293 + (-t280 * t285 - t287 * t311) * t290) * qJD(4) + (-t281 * t291 + t301 * t294) * t302;
t1 = [0, 0, -t279 * t285, t296 * qJD(3) + t278 * t291 - t279 * t309, t274, t274; 0, 0, t277 * t285, t297 * qJD(3) + t276 * t291 + t277 * t309, t273, t273; 0, 0, qJD(2) * t292 * t313, t289 * t285 * qJD(3) * t291 + (t298 * qJD(3) + (t288 * t304 + t305) * qJD(2)) * t286, t275, t275;];
JgD_rot  = t1;
