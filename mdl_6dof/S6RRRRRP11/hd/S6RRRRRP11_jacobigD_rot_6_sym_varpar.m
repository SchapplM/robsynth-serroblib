% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JgD_rot = S6RRRRRP11_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:40
% EndTime: 2019-02-26 22:45:40
% DurationCPUTime: 0.18s
% Computational Cost: add. (100->49), mult. (332->104), div. (0->0), fcn. (355->12), ass. (0->46)
t302 = sin(pkin(7));
t310 = cos(qJ(4));
t337 = t302 * t310;
t303 = sin(pkin(6));
t309 = sin(qJ(1));
t336 = t303 * t309;
t313 = cos(qJ(1));
t335 = t303 * t313;
t304 = cos(pkin(7));
t307 = sin(qJ(3));
t334 = t304 * t307;
t308 = sin(qJ(2));
t333 = t307 * t308;
t312 = cos(qJ(2));
t332 = t307 * t312;
t311 = cos(qJ(3));
t331 = t308 * t311;
t330 = t309 * t308;
t329 = t309 * t312;
t328 = t311 * t312;
t327 = t313 * t308;
t326 = t313 * t312;
t325 = qJD(1) * t303;
t306 = sin(qJ(4));
t324 = qJD(3) * t306;
t323 = t309 * t325;
t322 = t313 * t325;
t321 = t302 * t323;
t320 = t302 * t322;
t305 = cos(pkin(6));
t298 = t305 * t326 - t330;
t319 = t298 * t304 - t302 * t335;
t300 = -t305 * t329 - t327;
t318 = t300 * t304 + t302 * t336;
t317 = t304 * t332 + t331;
t299 = t305 * t327 + t329;
t316 = t305 * t330 - t326;
t315 = t299 * t311 + t319 * t307;
t314 = t318 * t307 - t311 * t316;
t297 = -t316 * qJD(1) + t298 * qJD(2);
t296 = t300 * qJD(1) - t299 * qJD(2);
t295 = -t299 * qJD(1) + t300 * qJD(2);
t294 = -t298 * qJD(1) + t316 * qJD(2);
t293 = -t296 * t302 + t304 * t323;
t292 = -t294 * t302 + t304 * t322;
t1 = [0, t322, t292, t295 * t307 + (-t294 * t304 - t320) * t311 + t314 * qJD(3) (t294 * t334 + t295 * t311 + t307 * t320) * t306 - t292 * t310 + (t314 * t310 + (-t300 * t302 + t304 * t336) * t306) * qJD(4) + (t307 * t316 + t318 * t311) * t324, 0; 0, t323, t293, t297 * t307 + (-t296 * t304 - t321) * t311 + t315 * qJD(3) (t296 * t334 + t297 * t311 + t307 * t321) * t306 - t293 * t310 + (t315 * t310 + (-t298 * t302 - t304 * t335) * t306) * qJD(4) + (-t299 * t307 + t319 * t311) * t324, 0; 0, 0, t303 * qJD(2) * t308 * t302, t305 * t302 * qJD(3) * t307 + (t317 * qJD(3) + (t304 * t331 + t332) * qJD(2)) * t303 (t302 * t311 * t324 + (t304 * t306 + t307 * t337) * qJD(4)) * t305 + ((-t302 * t312 * t306 + t317 * t310) * qJD(4) + (t304 * t328 - t333) * t324 + ((-t304 * t333 + t328) * t306 - t308 * t337) * qJD(2)) * t303, 0;];
JgD_rot  = t1;
