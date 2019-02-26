% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR15_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:16
% EndTime: 2019-02-26 22:24:16
% DurationCPUTime: 0.23s
% Computational Cost: add. (100->49), mult. (332->103), div. (0->0), fcn. (355->12), ass. (0->47)
t307 = sin(qJ(3));
t311 = cos(qJ(3));
t305 = cos(pkin(6));
t312 = cos(qJ(2));
t313 = cos(qJ(1));
t324 = t312 * t313;
t308 = sin(qJ(2));
t309 = sin(qJ(1));
t327 = t309 * t308;
t314 = t305 * t327 - t324;
t326 = t309 * t312;
t328 = t308 * t313;
t300 = -t305 * t326 - t328;
t302 = sin(pkin(7));
t304 = cos(pkin(7));
t303 = sin(pkin(6));
t334 = t303 * t309;
t316 = t300 * t304 + t302 * t334;
t340 = t307 * t314 + t316 * t311;
t299 = t305 * t328 + t326;
t298 = t305 * t324 - t327;
t333 = t303 * t313;
t317 = -t298 * t304 + t302 * t333;
t339 = t299 * t307 + t317 * t311;
t336 = t302 * t308;
t335 = t302 * t311;
t332 = t304 * t311;
t331 = t307 * t308;
t330 = t307 * t312;
t329 = t308 * t311;
t325 = t311 * t312;
t323 = qJD(1) * t303;
t310 = cos(qJ(5));
t322 = qJD(3) * t310;
t321 = t309 * t323;
t320 = t313 * t323;
t319 = t302 * t321;
t318 = t302 * t320;
t315 = t304 * t325 - t331;
t306 = sin(qJ(5));
t297 = -t314 * qJD(1) + t298 * qJD(2);
t296 = t300 * qJD(1) - t299 * qJD(2);
t295 = -t299 * qJD(1) + t300 * qJD(2);
t294 = -t298 * qJD(1) + t314 * qJD(2);
t293 = -t296 * t302 + t304 * t321;
t292 = -t294 * t302 + t304 * t320;
t1 = [0, t320, t292, 0, t295 * t311 + (t294 * t304 + t318) * t307 + t340 * qJD(3), t292 * t306 - (-t294 * t332 + t295 * t307 - t311 * t318) * t310 + ((-t300 * t302 + t304 * t334) * t310 - t340 * t306) * qJD(5) - (t316 * t307 - t311 * t314) * t322; 0, t321, t293, 0, t297 * t311 + (t296 * t304 + t319) * t307 - t339 * qJD(3), t293 * t306 - (-t296 * t332 + t297 * t307 - t311 * t319) * t310 + ((-t298 * t302 - t304 * t333) * t310 + t339 * t306) * qJD(5) - (t299 * t311 - t317 * t307) * t322; 0, 0, t303 * qJD(2) * t336, 0, t305 * qJD(3) * t335 + (t315 * qJD(3) + (-t304 * t331 + t325) * qJD(2)) * t303 (-t302 * t307 * t322 + (t304 * t310 - t306 * t335) * qJD(5)) * t305 + ((-t302 * t312 * t310 - t315 * t306) * qJD(5) - (t304 * t330 + t329) * t322 + (t306 * t336 - (t304 * t329 + t330) * t310) * qJD(2)) * t303;];
JgD_rot  = t1;
