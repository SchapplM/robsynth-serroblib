% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRP12_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:19
% EndTime: 2019-02-26 21:14:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (71->39), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->42)
t292 = sin(pkin(7));
t296 = cos(pkin(6));
t321 = t292 * t296;
t293 = sin(pkin(6));
t299 = sin(qJ(1));
t320 = t293 * t299;
t302 = cos(qJ(1));
t319 = t293 * t302;
t295 = cos(pkin(7));
t298 = sin(qJ(3));
t318 = t295 * t298;
t291 = sin(pkin(12));
t317 = t299 * t291;
t294 = cos(pkin(12));
t316 = t299 * t294;
t315 = t302 * t291;
t314 = t302 * t294;
t313 = qJD(1) * t293;
t297 = sin(qJ(4));
t312 = qJD(3) * t297;
t311 = t299 * t313;
t310 = t302 * t313;
t309 = t292 * t311;
t308 = t292 * t310;
t287 = t296 * t314 - t317;
t307 = t287 * t295 - t292 * t319;
t289 = -t296 * t316 - t315;
t306 = t289 * t295 + t292 * t320;
t288 = t296 * t315 + t316;
t290 = -t296 * t317 + t314;
t301 = cos(qJ(3));
t305 = t288 * t301 + t307 * t298;
t304 = t290 * t301 + t306 * t298;
t303 = t298 * t321 + (t291 * t301 + t294 * t318) * t293;
t300 = cos(qJ(4));
t286 = t290 * qJD(1);
t285 = t289 * qJD(1);
t284 = t288 * qJD(1);
t283 = t287 * qJD(1);
t282 = -t285 * t292 + t295 * t311;
t281 = t283 * t292 + t295 * t310;
t1 = [0, 0, t281, -t284 * t298 + (t283 * t295 - t308) * t301 + t304 * qJD(3) (-t283 * t318 - t284 * t301 + t298 * t308) * t297 - t281 * t300 + (t304 * t300 + (-t289 * t292 + t295 * t320) * t297) * qJD(4) + (-t290 * t298 + t306 * t301) * t312, 0; 0, 0, t282, t286 * t298 + (-t285 * t295 - t309) * t301 + t305 * qJD(3) (t285 * t318 + t286 * t301 + t298 * t309) * t297 - t282 * t300 + (t305 * t300 + (-t287 * t292 - t295 * t319) * t297) * qJD(4) + (-t288 * t298 + t307 * t301) * t312, 0; 0, 0, 0, t303 * qJD(3) (t303 * t300 + (-t293 * t294 * t292 + t296 * t295) * t297) * qJD(4) + (t301 * t321 + (t294 * t295 * t301 - t291 * t298) * t293) * t312, 0;];
JgD_rot  = t1;
