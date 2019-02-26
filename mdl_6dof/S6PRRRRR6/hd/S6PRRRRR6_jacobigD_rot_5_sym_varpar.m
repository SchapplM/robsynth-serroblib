% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR6
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR6_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobigD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:22:02
% EndTime: 2019-02-26 20:22:02
% DurationCPUTime: 0.24s
% Computational Cost: add. (111->50), mult. (386->112), div. (0->0), fcn. (423->14), ass. (0->51)
t325 = sin(qJ(3));
t328 = cos(qJ(3));
t316 = sin(pkin(14));
t320 = cos(pkin(14));
t329 = cos(qJ(2));
t323 = cos(pkin(6));
t326 = sin(qJ(2));
t345 = t323 * t326;
t334 = t316 * t345 - t320 * t329;
t344 = t323 * t329;
t314 = -t316 * t344 - t320 * t326;
t322 = cos(pkin(7));
t318 = sin(pkin(7));
t319 = sin(pkin(6));
t352 = t318 * t319;
t335 = t314 * t322 + t316 * t352;
t357 = t335 * t325 - t328 * t334;
t313 = t316 * t329 + t320 * t345;
t312 = -t316 * t326 + t320 * t344;
t336 = -t312 * t322 + t320 * t352;
t356 = -t313 * t328 + t336 * t325;
t317 = sin(pkin(8));
t353 = t318 * t317;
t321 = cos(pkin(8));
t351 = t318 * t321;
t350 = t318 * t323;
t349 = t318 * t326;
t348 = t319 * t322;
t347 = t322 * t325;
t346 = t322 * t328;
t343 = t325 * t326;
t342 = t325 * t329;
t341 = t326 * t328;
t340 = t328 * t329;
t324 = sin(qJ(4));
t339 = qJD(3) * t324;
t338 = qJD(3) * t350;
t337 = t319 * qJD(2) * t349;
t333 = t322 * t340 - t343;
t332 = t322 * t342 + t341;
t331 = -t313 * t325 - t336 * t328;
t330 = t325 * t334 + t335 * t328;
t327 = cos(qJ(4));
t311 = t334 * qJD(2);
t310 = t314 * qJD(2);
t309 = t313 * qJD(2);
t308 = t312 * qJD(2);
t307 = -t325 * t338 + (-t332 * qJD(3) + (-t322 * t341 - t342) * qJD(2)) * t319;
t306 = -qJD(3) * t357 - t310 * t325 + t311 * t346;
t305 = qJD(3) * t356 - t308 * t325 - t309 * t346;
t1 = [0, 0, -t311 * t318, -t306 * t317 - t311 * t351 (t310 * t328 + t311 * t347) * t324 + (-t306 * t321 + t311 * t353) * t327 + t330 * t339 + (t357 * t327 + (t330 * t321 + (-t314 * t318 + t316 * t348) * t317) * t324) * qJD(4), 0; 0, 0, t309 * t318, -t305 * t317 + t309 * t351 (t308 * t328 - t309 * t347) * t324 + (-t305 * t321 - t309 * t353) * t327 + t331 * t339 + (-t356 * t327 + (t331 * t321 + (-t312 * t318 - t320 * t348) * t317) * t324) * qJD(4), 0; 0, 0, t337, -t307 * t317 + t321 * t337, t328 * t324 * t338 - t307 * t321 * t327 + (t333 * t339 + ((-t322 * t343 + t340) * t324 - t317 * t327 * t349) * qJD(2)) * t319 + ((t332 * t319 + t325 * t350) * t327 + ((t333 * t319 + t328 * t350) * t321 + (t323 * t322 - t329 * t352) * t317) * t324) * qJD(4), 0;];
JgD_rot  = t1;
