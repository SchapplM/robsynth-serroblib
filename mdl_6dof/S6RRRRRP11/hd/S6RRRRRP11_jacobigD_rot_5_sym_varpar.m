% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRP11_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:45
% EndTime: 2019-02-26 22:45:45
% DurationCPUTime: 0.20s
% Computational Cost: add. (100->49), mult. (332->104), div. (0->0), fcn. (355->12), ass. (0->46)
t301 = sin(pkin(7));
t309 = cos(qJ(4));
t336 = t301 * t309;
t302 = sin(pkin(6));
t308 = sin(qJ(1));
t335 = t302 * t308;
t312 = cos(qJ(1));
t334 = t302 * t312;
t303 = cos(pkin(7));
t306 = sin(qJ(3));
t333 = t303 * t306;
t307 = sin(qJ(2));
t332 = t306 * t307;
t311 = cos(qJ(2));
t331 = t306 * t311;
t310 = cos(qJ(3));
t330 = t307 * t310;
t329 = t308 * t307;
t328 = t308 * t311;
t327 = t310 * t311;
t326 = t312 * t307;
t325 = t312 * t311;
t324 = qJD(1) * t302;
t305 = sin(qJ(4));
t323 = qJD(3) * t305;
t322 = t308 * t324;
t321 = t312 * t324;
t320 = t301 * t322;
t319 = t301 * t321;
t304 = cos(pkin(6));
t297 = t304 * t325 - t329;
t318 = t297 * t303 - t301 * t334;
t299 = -t304 * t328 - t326;
t317 = t299 * t303 + t301 * t335;
t316 = t303 * t331 + t330;
t298 = t304 * t326 + t328;
t315 = t304 * t329 - t325;
t314 = t298 * t310 + t318 * t306;
t313 = t317 * t306 - t310 * t315;
t296 = -t315 * qJD(1) + t297 * qJD(2);
t295 = t299 * qJD(1) - t298 * qJD(2);
t294 = -t298 * qJD(1) + t299 * qJD(2);
t293 = -t297 * qJD(1) + t315 * qJD(2);
t292 = -t295 * t301 + t303 * t322;
t291 = -t293 * t301 + t303 * t321;
t1 = [0, t321, t291, t294 * t306 + (-t293 * t303 - t319) * t310 + t313 * qJD(3) (t293 * t333 + t294 * t310 + t306 * t319) * t305 - t291 * t309 + (t313 * t309 + (-t299 * t301 + t303 * t335) * t305) * qJD(4) + (t306 * t315 + t317 * t310) * t323, 0; 0, t322, t292, t296 * t306 + (-t295 * t303 - t320) * t310 + t314 * qJD(3) (t295 * t333 + t296 * t310 + t306 * t320) * t305 - t292 * t309 + (t314 * t309 + (-t297 * t301 - t303 * t334) * t305) * qJD(4) + (-t298 * t306 + t318 * t310) * t323, 0; 0, 0, t302 * qJD(2) * t307 * t301, t304 * t301 * qJD(3) * t306 + (t316 * qJD(3) + (t303 * t330 + t331) * qJD(2)) * t302 (t301 * t310 * t323 + (t303 * t305 + t306 * t336) * qJD(4)) * t304 + ((-t301 * t311 * t305 + t316 * t309) * qJD(4) + (t303 * t327 - t332) * t323 + ((-t303 * t332 + t327) * t305 - t307 * t336) * qJD(2)) * t302, 0;];
JgD_rot  = t1;
