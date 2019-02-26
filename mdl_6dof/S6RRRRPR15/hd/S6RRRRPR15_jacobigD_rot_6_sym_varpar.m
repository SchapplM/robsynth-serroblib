% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR15_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:48
% EndTime: 2019-02-26 22:38:48
% DurationCPUTime: 0.24s
% Computational Cost: add. (100->49), mult. (332->103), div. (0->0), fcn. (355->12), ass. (0->47)
t297 = sin(pkin(7));
t302 = sin(qJ(3));
t333 = t297 * t302;
t303 = sin(qJ(2));
t332 = t297 * t303;
t298 = sin(pkin(6));
t304 = sin(qJ(1));
t331 = t298 * t304;
t308 = cos(qJ(1));
t330 = t298 * t308;
t299 = cos(pkin(7));
t329 = t299 * t302;
t328 = t302 * t303;
t307 = cos(qJ(2));
t327 = t302 * t307;
t306 = cos(qJ(3));
t326 = t303 * t306;
t325 = t304 * t303;
t324 = t304 * t307;
t323 = t306 * t307;
t322 = t308 * t303;
t321 = t308 * t307;
t320 = qJD(1) * t298;
t305 = cos(qJ(4));
t319 = qJD(3) * t305;
t318 = t304 * t320;
t317 = t308 * t320;
t316 = t297 * t318;
t315 = t297 * t317;
t300 = cos(pkin(6));
t293 = t300 * t321 - t325;
t314 = t293 * t299 - t297 * t330;
t295 = -t300 * t324 - t322;
t313 = t295 * t299 + t297 * t331;
t312 = t299 * t327 + t326;
t294 = t300 * t322 + t324;
t311 = t300 * t325 - t321;
t310 = t294 * t306 + t314 * t302;
t309 = t313 * t302 - t306 * t311;
t301 = sin(qJ(4));
t292 = -t311 * qJD(1) + t293 * qJD(2);
t291 = t295 * qJD(1) - t294 * qJD(2);
t290 = -t294 * qJD(1) + t295 * qJD(2);
t289 = -t293 * qJD(1) + t311 * qJD(2);
t288 = -t291 * t297 + t299 * t318;
t287 = -t289 * t297 + t299 * t317;
t1 = [0, t317, t287, t290 * t302 + (-t289 * t299 - t315) * t306 + t309 * qJD(3), 0 (t289 * t329 + t290 * t306 + t302 * t315) * t305 + t287 * t301 + (-t309 * t301 + (-t295 * t297 + t299 * t331) * t305) * qJD(4) + (t302 * t311 + t313 * t306) * t319; 0, t318, t288, t292 * t302 + (-t291 * t299 - t316) * t306 + t310 * qJD(3), 0 (t291 * t329 + t292 * t306 + t302 * t316) * t305 + t288 * t301 + (-t310 * t301 + (-t293 * t297 - t299 * t330) * t305) * qJD(4) + (-t294 * t302 + t314 * t306) * t319; 0, 0, t298 * qJD(2) * t332, t300 * qJD(3) * t333 + (t312 * qJD(3) + (t299 * t326 + t327) * qJD(2)) * t298, 0 (t297 * t306 * t319 + (t299 * t305 - t301 * t333) * qJD(4)) * t300 + ((-t297 * t307 * t305 - t312 * t301) * qJD(4) + (t299 * t323 - t328) * t319 + ((-t299 * t328 + t323) * t305 + t301 * t332) * qJD(2)) * t298;];
JgD_rot  = t1;
