% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRP6_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:03
% EndTime: 2019-02-26 20:18:03
% DurationCPUTime: 0.14s
% Computational Cost: add. (70->40), mult. (240->91), div. (0->0), fcn. (263->12), ass. (0->38)
t292 = sin(pkin(7));
t293 = sin(pkin(6));
t320 = t292 * t293;
t300 = cos(qJ(4));
t319 = t292 * t300;
t295 = cos(pkin(7));
t318 = t293 * t295;
t298 = sin(qJ(3));
t317 = t295 * t298;
t301 = cos(qJ(3));
t316 = t295 * t301;
t296 = cos(pkin(6));
t299 = sin(qJ(2));
t315 = t296 * t299;
t302 = cos(qJ(2));
t314 = t296 * t302;
t313 = t298 * t299;
t312 = t298 * t302;
t311 = t299 * t301;
t310 = t301 * t302;
t297 = sin(qJ(4));
t309 = qJD(3) * t297;
t291 = sin(pkin(12));
t294 = cos(pkin(12));
t287 = -t291 * t299 + t294 * t314;
t308 = t287 * t295 - t294 * t320;
t289 = -t291 * t314 - t294 * t299;
t307 = t289 * t295 + t291 * t320;
t288 = t291 * t302 + t294 * t315;
t306 = t291 * t315 - t294 * t302;
t305 = t295 * t312 + t311;
t304 = t288 * t301 + t308 * t298;
t303 = t307 * t298 - t301 * t306;
t286 = t306 * qJD(2);
t285 = t289 * qJD(2);
t284 = t288 * qJD(2);
t283 = t287 * qJD(2);
t1 = [0, 0, -t286 * t292, t303 * qJD(3) + t285 * t298 - t286 * t316 (t285 * t301 + t286 * t317) * t297 + t286 * t319 + (t303 * t300 + (-t289 * t292 + t291 * t318) * t297) * qJD(4) + (t298 * t306 + t307 * t301) * t309, 0; 0, 0, t284 * t292, t304 * qJD(3) + t283 * t298 + t284 * t316 (t283 * t301 - t284 * t317) * t297 - t284 * t319 + (t304 * t300 + (-t287 * t292 - t294 * t318) * t297) * qJD(4) + (-t288 * t298 + t308 * t301) * t309, 0; 0, 0, qJD(2) * t299 * t320, t296 * t292 * qJD(3) * t298 + (t305 * qJD(3) + (t295 * t311 + t312) * qJD(2)) * t293 (t292 * t301 * t309 + (t295 * t297 + t298 * t319) * qJD(4)) * t296 + ((-t292 * t302 * t297 + t305 * t300) * qJD(4) + (t295 * t310 - t313) * t309 + ((-t295 * t313 + t310) * t297 - t299 * t319) * qJD(2)) * t293, 0;];
JgD_rot  = t1;
