% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR4
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
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:37
% EndTime: 2019-02-26 20:20:37
% DurationCPUTime: 0.17s
% Computational Cost: add. (108->42), mult. (313->90), div. (0->0), fcn. (340->12), ass. (0->46)
t296 = sin(pkin(13));
t299 = cos(pkin(13));
t305 = cos(qJ(2));
t301 = cos(pkin(6));
t303 = sin(qJ(2));
t319 = t301 * t303;
t289 = t296 * t305 + t299 * t319;
t285 = t289 * qJD(2);
t297 = sin(pkin(7));
t326 = t285 * t297;
t309 = t296 * t319 - t299 * t305;
t287 = t309 * qJD(2);
t325 = t287 * t297;
t298 = sin(pkin(6));
t324 = t297 * t298;
t323 = t297 * t303;
t300 = cos(pkin(7));
t322 = t298 * t300;
t302 = sin(qJ(3));
t321 = t300 * t302;
t304 = cos(qJ(3));
t320 = t300 * t304;
t318 = t301 * t305;
t317 = t302 * t303;
t316 = t302 * t305;
t315 = t303 * t304;
t314 = t304 * t305;
t295 = qJ(4) + qJ(5);
t292 = sin(t295);
t313 = qJD(3) * t292;
t312 = qJD(3) * t297;
t288 = -t296 * t303 + t299 * t318;
t311 = t288 * t300 - t299 * t324;
t290 = -t296 * t318 - t299 * t303;
t310 = t290 * t300 + t296 * t324;
t308 = t300 * t316 + t315;
t307 = t289 * t304 + t302 * t311;
t306 = t302 * t310 - t304 * t309;
t294 = qJD(4) + qJD(5);
t293 = cos(t295);
t286 = t290 * qJD(2);
t284 = t288 * qJD(2);
t283 = t301 * t302 * t312 + (t308 * qJD(3) + (t300 * t315 + t316) * qJD(2)) * t298;
t282 = qJD(3) * t306 + t286 * t302 - t287 * t320;
t281 = qJD(3) * t307 + t284 * t302 + t285 * t320;
t1 = [0, 0, -t325, t282, t282 (t294 * t306 + t325) * t293 + (t287 * t321 + t286 * t304 + (-t290 * t297 + t296 * t322) * t294) * t292 + (t302 * t309 + t304 * t310) * t313; 0, 0, t326, t281, t281 (t294 * t307 - t326) * t293 + (-t285 * t321 + t284 * t304 + (-t288 * t297 - t299 * t322) * t294) * t292 + (-t289 * t302 + t304 * t311) * t313; 0, 0, t298 * qJD(2) * t323, t283, t283 (t297 * t302 * t294 * t293 + (t294 * t300 + t304 * t312) * t292) * t301 + ((-t297 * t305 * t292 + t293 * t308) * t294 + (t300 * t314 - t317) * t313 + ((-t300 * t317 + t314) * t292 - t293 * t323) * qJD(2)) * t298;];
JgD_rot  = t1;
