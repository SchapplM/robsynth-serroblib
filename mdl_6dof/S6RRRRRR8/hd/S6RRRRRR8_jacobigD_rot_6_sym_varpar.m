% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR8
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
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR8_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:24
% EndTime: 2019-02-26 22:51:24
% DurationCPUTime: 0.22s
% Computational Cost: add. (148->49), mult. (433->99), div. (0->0), fcn. (460->12), ass. (0->51)
t346 = sin(pkin(7));
t351 = sin(qJ(2));
t379 = t346 * t351;
t347 = sin(pkin(6));
t352 = sin(qJ(1));
t378 = t347 * t352;
t355 = cos(qJ(1));
t377 = t347 * t355;
t350 = sin(qJ(3));
t376 = t350 * t351;
t354 = cos(qJ(2));
t375 = t350 * t354;
t353 = cos(qJ(3));
t374 = t351 * t353;
t373 = t351 * t355;
t372 = t352 * t351;
t371 = t352 * t354;
t370 = t353 * t354;
t369 = t354 * t355;
t368 = qJD(1) * t347;
t345 = qJ(4) + qJ(5);
t342 = sin(t345);
t367 = qJD(3) * t342;
t366 = qJD(3) * t346;
t365 = t352 * t368;
t364 = t355 * t368;
t349 = cos(pkin(6));
t338 = t349 * t369 - t372;
t348 = cos(pkin(7));
t363 = t338 * t348 - t346 * t377;
t340 = -t349 * t371 - t373;
t362 = t340 * t348 + t346 * t378;
t361 = t348 * t375 + t374;
t339 = t349 * t373 + t371;
t360 = t349 * t372 - t369;
t334 = -qJD(1) * t338 + qJD(2) * t360;
t359 = t334 * t348 + t346 * t364;
t336 = qJD(1) * t340 - qJD(2) * t339;
t358 = t336 * t348 + t346 * t365;
t357 = t339 * t353 + t350 * t363;
t356 = t350 * t362 - t353 * t360;
t344 = qJD(4) + qJD(5);
t343 = cos(t345);
t337 = -qJD(1) * t360 + qJD(2) * t338;
t335 = -qJD(1) * t339 + qJD(2) * t340;
t333 = -t336 * t346 + t348 * t365;
t332 = -t334 * t346 + t348 * t364;
t331 = t349 * t350 * t366 + (t361 * qJD(3) + (t348 * t374 + t375) * qJD(2)) * t347;
t330 = qJD(3) * t357 + t337 * t350 - t353 * t358;
t329 = qJD(3) * t356 + t335 * t350 - t353 * t359;
t1 = [0, t364, t332, t329, t329 (t344 * t356 - t332) * t343 + (t335 * t353 + (-t340 * t346 + t348 * t378) * t344 + t359 * t350) * t342 + (t350 * t360 + t353 * t362) * t367; 0, t365, t333, t330, t330 (t344 * t357 - t333) * t343 + (t337 * t353 + (-t338 * t346 - t348 * t377) * t344 + t358 * t350) * t342 + (-t339 * t350 + t353 * t363) * t367; 0, 0, t347 * qJD(2) * t379, t331, t331 (t346 * t350 * t344 * t343 + (t344 * t348 + t353 * t366) * t342) * t349 + ((-t346 * t354 * t342 + t343 * t361) * t344 + (t348 * t370 - t376) * t367 + ((-t348 * t376 + t370) * t342 - t343 * t379) * qJD(2)) * t347;];
JgD_rot  = t1;
