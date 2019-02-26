% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:03
% EndTime: 2019-02-26 20:10:04
% DurationCPUTime: 0.48s
% Computational Cost: add. (596->104), mult. (1844->182), div. (0->0), fcn. (1930->10), ass. (0->72)
t353 = sin(qJ(2));
t356 = cos(qJ(2));
t349 = sin(pkin(10));
t400 = cos(pkin(6));
t379 = t349 * t400;
t399 = cos(pkin(10));
t339 = -t353 * t379 + t356 * t399;
t352 = sin(qJ(3));
t355 = cos(qJ(3));
t385 = pkin(5) + pkin(9) + r_i_i_C(1);
t403 = -pkin(3) * t352 + t355 * t385;
t401 = r_i_i_C(2) + qJ(5);
t350 = sin(pkin(6));
t398 = t350 * t353;
t397 = t350 * t355;
t351 = sin(qJ(4));
t396 = t351 * t355;
t395 = t351 * t356;
t354 = cos(qJ(4));
t394 = t354 * t355;
t393 = t354 * t356;
t392 = t355 * t356;
t391 = qJD(2) * t353;
t390 = qJD(2) * t355;
t389 = qJD(3) * t352;
t388 = qJD(3) * t356;
t387 = qJD(4) * t354;
t386 = qJD(4) * t355;
t384 = -r_i_i_C(3) - qJ(6) - pkin(4);
t383 = t350 * t391;
t382 = qJD(2) * t350 * t356;
t381 = qJD(4) * t395;
t380 = t352 * t388;
t378 = t350 * t399;
t373 = t400 * t399;
t368 = t356 * t373;
t332 = -qJD(2) * t368 + t349 * t391;
t336 = t349 * t353 - t368;
t375 = t336 * t386 - t332;
t338 = t353 * t399 + t356 * t379;
t334 = t338 * qJD(2);
t374 = t338 * t386 - t334;
t337 = t349 * t356 + t353 * t373;
t362 = -t337 * t355 + t352 * t378;
t372 = t336 * t351 - t354 * t362;
t371 = -t336 * t354 - t351 * t362;
t327 = t349 * t350 * t352 + t339 * t355;
t370 = t327 * t354 + t338 * t351;
t369 = t327 * t351 - t338 * t354;
t367 = -t339 * t352 + t349 * t397;
t341 = t352 * t400 + t353 * t397;
t366 = t341 * t351 + t350 * t393;
t365 = -t355 * pkin(3) - t352 * t385 - pkin(2);
t364 = -t337 * t352 - t355 * t378;
t363 = -t352 * t398 + t355 * t400;
t361 = qJD(3) * t403;
t333 = t337 * qJD(2);
t360 = qJD(4) * t337 - t333 * t355 + t336 * t389;
t335 = t339 * qJD(2);
t359 = qJD(4) * t339 - t335 * t355 + t338 * t389;
t358 = t351 * t401 - t354 * t384 + pkin(3);
t357 = t351 * qJD(5) + t354 * qJD(6) + (t351 * t384 + t354 * t401) * qJD(4);
t329 = qJD(3) * t363 + t355 * t382;
t323 = qJD(3) * t367 - t334 * t355;
t321 = qJD(3) * t364 - t332 * t355;
t317 = -qJD(4) * t366 + t329 * t354 + t351 * t383;
t316 = t329 * t351 + t341 * t387 - t350 * t381 - t354 * t383;
t311 = -qJD(4) * t369 + t323 * t354 + t335 * t351;
t310 = qJD(4) * t370 + t323 * t351 - t335 * t354;
t309 = -qJD(4) * t371 + t321 * t354 + t333 * t351;
t308 = qJD(4) * t372 + t321 * t351 - t333 * t354;
t1 = [0 -(t338 * t394 - t339 * t351) * qJD(6) - (t338 * t396 + t339 * t354) * qJD(5) - t334 * pkin(8) + t401 * (t351 * t359 - t354 * t374) - t384 * (t351 * t374 + t354 * t359) - t338 * t361 + t365 * t335, t385 * t323 + t358 * (-qJD(3) * t327 + t334 * t352) + t357 * t367, qJD(5) * t370 - qJD(6) * t369 + t310 * t384 + t311 * t401, t310, t311; 0 -(t336 * t394 - t337 * t351) * qJD(6) - (t336 * t396 + t337 * t354) * qJD(5) - t332 * pkin(8) + t401 * (t351 * t360 - t354 * t375) - t384 * (t351 * t375 + t354 * t360) - t336 * t361 + t365 * t333, t385 * t321 + t358 * (qJD(3) * t362 + t332 * t352) + t357 * t364, qJD(5) * t372 - qJD(6) * t371 + t308 * t384 + t309 * t401, t308, t309; 0, t384 * (-t351 * t382 - t387 * t398) + (-t401 * ((qJD(2) - t386) * t393 + (t380 + (-qJD(4) + t390) * t353) * t351) + t384 * (t355 * t381 + (t353 * t390 + t380) * t354) - (-t351 * t353 - t354 * t392) * qJD(6) - (-t351 * t392 + t353 * t354) * qJD(5) + t403 * t388 + (t356 * pkin(8) + t353 * t365) * qJD(2)) * t350, t385 * t329 + t358 * (-qJD(3) * t341 - t352 * t382) + t357 * t363, -t366 * qJD(6) - (-t341 * t354 + t350 * t395) * qJD(5) + t401 * t317 + t384 * t316, t316, t317;];
JaD_transl  = t1;
