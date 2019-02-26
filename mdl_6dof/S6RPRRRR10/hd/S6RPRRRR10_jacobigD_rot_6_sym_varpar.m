% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR10
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

function JgD_rot = S6RPRRRR10_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:51
% EndTime: 2019-02-26 21:19:51
% DurationCPUTime: 0.17s
% Computational Cost: add. (108->39), mult. (321->80), div. (0->0), fcn. (348->12), ass. (0->47)
t303 = sin(pkin(7));
t307 = cos(pkin(6));
t330 = t303 * t307;
t304 = sin(pkin(6));
t309 = sin(qJ(1));
t329 = t304 * t309;
t311 = cos(qJ(1));
t328 = t304 * t311;
t305 = cos(pkin(13));
t306 = cos(pkin(7));
t327 = t305 * t306;
t302 = sin(pkin(13));
t326 = t309 * t302;
t325 = t309 * t305;
t324 = t311 * t302;
t323 = t311 * t305;
t322 = qJD(1) * t304;
t301 = qJ(4) + qJ(5);
t298 = sin(t301);
t321 = qJD(3) * t298;
t320 = t309 * t322;
t319 = t311 * t322;
t294 = t307 * t323 - t326;
t318 = t294 * t306 - t303 * t328;
t296 = -t307 * t325 - t324;
t317 = t296 * t306 + t303 * t329;
t295 = t307 * t324 + t325;
t297 = -t307 * t326 + t323;
t290 = t294 * qJD(1);
t316 = -t290 * t306 + t303 * t319;
t292 = t296 * qJD(1);
t315 = t292 * t306 + t303 * t320;
t308 = sin(qJ(3));
t310 = cos(qJ(3));
t314 = t295 * t310 + t318 * t308;
t313 = t297 * t310 + t317 * t308;
t312 = t308 * t330 + (t302 * t310 + t308 * t327) * t304;
t300 = qJD(4) + qJD(5);
t299 = cos(t301);
t293 = t297 * qJD(1);
t291 = t295 * qJD(1);
t289 = -t292 * t303 + t306 * t320;
t288 = t290 * t303 + t306 * t319;
t287 = t312 * qJD(3);
t286 = t314 * qJD(3) + t293 * t308 - t315 * t310;
t285 = t313 * qJD(3) - t291 * t308 - t316 * t310;
t1 = [0, 0, t288, t285, t285 (t313 * t300 - t288) * t299 + (-t291 * t310 + (-t296 * t303 + t306 * t329) * t300 + t316 * t308) * t298 + (-t297 * t308 + t317 * t310) * t321; 0, 0, t289, t286, t286 (t314 * t300 - t289) * t299 + (t293 * t310 + (-t294 * t303 - t306 * t328) * t300 + t315 * t308) * t298 + (-t295 * t308 + t318 * t310) * t321; 0, 0, 0, t287, t287 (t312 * t299 + (-t304 * t305 * t303 + t307 * t306) * t298) * t300 + (t310 * t330 + (-t302 * t308 + t310 * t327) * t304) * t321;];
JgD_rot  = t1;
