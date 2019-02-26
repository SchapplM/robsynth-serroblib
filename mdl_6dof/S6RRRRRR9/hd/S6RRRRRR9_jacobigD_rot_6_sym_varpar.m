% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR9_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:08
% EndTime: 2019-02-26 22:52:09
% DurationCPUTime: 0.22s
% Computational Cost: add. (162->49), mult. (528->104), div. (0->0), fcn. (570->12), ass. (0->49)
t327 = sin(pkin(7));
t335 = cos(qJ(4));
t362 = t327 * t335;
t328 = sin(pkin(6));
t334 = sin(qJ(1));
t361 = t328 * t334;
t338 = cos(qJ(1));
t360 = t328 * t338;
t329 = cos(pkin(7));
t332 = sin(qJ(3));
t359 = t329 * t332;
t333 = sin(qJ(2));
t358 = t332 * t333;
t337 = cos(qJ(2));
t357 = t332 * t337;
t336 = cos(qJ(3));
t356 = t333 * t336;
t355 = t334 * t333;
t354 = t334 * t337;
t353 = t336 * t337;
t352 = t338 * t333;
t351 = t338 * t337;
t350 = qJD(1) * t328;
t331 = sin(qJ(4));
t349 = qJD(3) * t331;
t348 = t334 * t350;
t347 = t338 * t350;
t346 = t327 * t348;
t345 = t327 * t347;
t330 = cos(pkin(6));
t323 = t330 * t351 - t355;
t344 = t323 * t329 - t327 * t360;
t325 = -t330 * t354 - t352;
t343 = t325 * t329 + t327 * t361;
t342 = t329 * t357 + t356;
t324 = t330 * t352 + t354;
t341 = t330 * t355 - t351;
t340 = t324 * t336 + t344 * t332;
t339 = t343 * t332 - t336 * t341;
t322 = -t341 * qJD(1) + t323 * qJD(2);
t321 = t325 * qJD(1) - t324 * qJD(2);
t320 = -t324 * qJD(1) + t325 * qJD(2);
t319 = -t323 * qJD(1) + t341 * qJD(2);
t318 = -t321 * t327 + t329 * t348;
t317 = -t319 * t327 + t329 * t347;
t316 = (t327 * t336 * t349 + (t329 * t331 + t332 * t362) * qJD(4)) * t330 + ((-t327 * t337 * t331 + t342 * t335) * qJD(4) + (t329 * t353 - t358) * t349 + ((-t329 * t358 + t353) * t331 - t333 * t362) * qJD(2)) * t328;
t315 = (t321 * t359 + t322 * t336 + t332 * t346) * t331 - t318 * t335 + (t340 * t335 + (-t323 * t327 - t329 * t360) * t331) * qJD(4) + (-t324 * t332 + t344 * t336) * t349;
t314 = (t319 * t359 + t320 * t336 + t332 * t345) * t331 - t317 * t335 + (t339 * t335 + (-t325 * t327 + t329 * t361) * t331) * qJD(4) + (t332 * t341 + t343 * t336) * t349;
t1 = [0, t347, t317, t320 * t332 + (-t319 * t329 - t345) * t336 + t339 * qJD(3), t314, t314; 0, t348, t318, t322 * t332 + (-t321 * t329 - t346) * t336 + t340 * qJD(3), t315, t315; 0, 0, t328 * qJD(2) * t333 * t327, t330 * t327 * qJD(3) * t332 + (t342 * qJD(3) + (t329 * t356 + t357) * qJD(2)) * t328, t316, t316;];
JgD_rot  = t1;
