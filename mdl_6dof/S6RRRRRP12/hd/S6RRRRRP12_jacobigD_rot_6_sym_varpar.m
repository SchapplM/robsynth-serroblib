% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP12
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
% Datum: 2019-02-26 22:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRP12_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:46:24
% EndTime: 2019-02-26 22:46:24
% DurationCPUTime: 0.21s
% Computational Cost: add. (100->49), mult. (332->104), div. (0->0), fcn. (355->12), ass. (0->46)
t340 = sin(pkin(7));
t348 = cos(qJ(4));
t375 = t340 * t348;
t341 = sin(pkin(6));
t347 = sin(qJ(1));
t374 = t341 * t347;
t351 = cos(qJ(1));
t373 = t341 * t351;
t342 = cos(pkin(7));
t345 = sin(qJ(3));
t372 = t342 * t345;
t346 = sin(qJ(2));
t371 = t345 * t346;
t350 = cos(qJ(2));
t370 = t345 * t350;
t349 = cos(qJ(3));
t369 = t346 * t349;
t368 = t347 * t346;
t367 = t347 * t350;
t366 = t349 * t350;
t365 = t351 * t346;
t364 = t351 * t350;
t363 = qJD(1) * t341;
t344 = sin(qJ(4));
t362 = qJD(3) * t344;
t361 = t347 * t363;
t360 = t351 * t363;
t359 = t340 * t361;
t358 = t340 * t360;
t343 = cos(pkin(6));
t336 = t343 * t364 - t368;
t357 = t336 * t342 - t340 * t373;
t338 = -t343 * t367 - t365;
t356 = t338 * t342 + t340 * t374;
t355 = t342 * t370 + t369;
t337 = t343 * t365 + t367;
t354 = t343 * t368 - t364;
t353 = t337 * t349 + t357 * t345;
t352 = t356 * t345 - t349 * t354;
t335 = -t354 * qJD(1) + t336 * qJD(2);
t334 = t338 * qJD(1) - t337 * qJD(2);
t333 = -t337 * qJD(1) + t338 * qJD(2);
t332 = -t336 * qJD(1) + t354 * qJD(2);
t331 = -t334 * t340 + t342 * t361;
t330 = -t332 * t340 + t342 * t360;
t1 = [0, t360, t330, t333 * t345 + (-t332 * t342 - t358) * t349 + t352 * qJD(3) (t332 * t372 + t333 * t349 + t345 * t358) * t344 - t330 * t348 + (t352 * t348 + (-t338 * t340 + t342 * t374) * t344) * qJD(4) + (t345 * t354 + t356 * t349) * t362, 0; 0, t361, t331, t335 * t345 + (-t334 * t342 - t359) * t349 + t353 * qJD(3) (t334 * t372 + t335 * t349 + t345 * t359) * t344 - t331 * t348 + (t353 * t348 + (-t336 * t340 - t342 * t373) * t344) * qJD(4) + (-t337 * t345 + t357 * t349) * t362, 0; 0, 0, t341 * qJD(2) * t346 * t340, t343 * t340 * qJD(3) * t345 + (t355 * qJD(3) + (t342 * t369 + t370) * qJD(2)) * t341 (t340 * t349 * t362 + (t342 * t344 + t345 * t375) * qJD(4)) * t343 + ((-t340 * t350 * t344 + t355 * t348) * qJD(4) + (t342 * t366 - t371) * t362 + ((-t342 * t371 + t366) * t344 - t346 * t375) * qJD(2)) * t341, 0;];
JgD_rot  = t1;
