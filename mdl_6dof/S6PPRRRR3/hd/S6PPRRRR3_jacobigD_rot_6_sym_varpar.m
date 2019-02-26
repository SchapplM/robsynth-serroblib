% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRRR3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobigD_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:57
% EndTime: 2019-02-26 19:43:57
% DurationCPUTime: 0.24s
% Computational Cost: add. (195->51), mult. (628->110), div. (0->0), fcn. (751->16), ass. (0->54)
t330 = sin(pkin(14));
t334 = sin(pkin(6));
t342 = sin(qJ(3));
t345 = cos(qJ(3));
t335 = cos(pkin(14));
t338 = cos(pkin(7));
t360 = t335 * t338;
t333 = sin(pkin(7));
t339 = cos(pkin(6));
t362 = t333 * t339;
t322 = (t330 * t345 + t342 * t360) * t334 + t342 * t362;
t336 = cos(pkin(13));
t331 = sin(pkin(13));
t365 = t331 * t339;
t329 = -t330 * t365 + t336 * t335;
t328 = -t336 * t330 - t335 * t365;
t363 = t333 * t334;
t350 = t328 * t338 + t331 * t363;
t318 = t329 * t345 + t350 * t342;
t359 = t336 * t339;
t327 = t330 * t359 + t331 * t335;
t326 = -t331 * t330 + t335 * t359;
t351 = -t326 * t338 + t336 * t363;
t368 = -t327 * t345 + t351 * t342;
t332 = sin(pkin(8));
t343 = cos(qJ(5));
t364 = t332 * t343;
t361 = t334 * t338;
t337 = cos(pkin(8));
t341 = sin(qJ(4));
t358 = t337 * t341;
t344 = cos(qJ(4));
t357 = t337 * t344;
t340 = sin(qJ(5));
t356 = qJD(4) * t340;
t315 = -t327 * t342 - t351 * t345;
t323 = -t326 * t333 - t336 * t361;
t354 = t315 * t337 + t323 * t332;
t317 = -t329 * t342 + t350 * t345;
t324 = -t328 * t333 + t331 * t361;
t353 = t317 * t337 + t324 * t332;
t321 = t345 * t362 + (-t330 * t342 + t345 * t360) * t334;
t325 = -t335 * t363 + t339 * t338;
t352 = t321 * t337 + t325 * t332;
t348 = t354 * t341 - t344 * t368;
t347 = t318 * t344 + t353 * t341;
t346 = t322 * t344 + t352 * t341;
t320 = t322 * qJD(3);
t319 = t321 * qJD(3);
t314 = t318 * qJD(3);
t313 = t317 * qJD(3);
t312 = t368 * qJD(3);
t311 = t315 * qJD(3);
t1 = [0, 0, 0, t314 * t332, t347 * qJD(4) + t313 * t341 + t314 * t357 (t313 * t344 - t314 * t358) * t340 - t314 * t364 + (t347 * t343 + (-t317 * t332 + t324 * t337) * t340) * qJD(5) + (-t318 * t341 + t353 * t344) * t356; 0, 0, 0, -t312 * t332, t348 * qJD(4) + t311 * t341 - t312 * t357 (t311 * t344 + t312 * t358) * t340 + t312 * t364 + (t348 * t343 + (-t315 * t332 + t323 * t337) * t340) * qJD(5) + (t341 * t368 + t354 * t344) * t356; 0, 0, 0, t320 * t332, t346 * qJD(4) + t319 * t341 + t320 * t357 (t319 * t344 - t320 * t358) * t340 - t320 * t364 + (t346 * t343 + (-t321 * t332 + t325 * t337) * t340) * qJD(5) + (-t322 * t341 + t352 * t344) * t356;];
JgD_rot  = t1;
