% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_transl = S6PRRRPP3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:04
% EndTime: 2019-02-26 20:10:04
% DurationCPUTime: 0.36s
% Computational Cost: add. (437->91), mult. (1367->160), div. (0->0), fcn. (1417->10), ass. (0->64)
t353 = sin(qJ(3));
t356 = cos(qJ(3));
t395 = pkin(9) + r_i_i_C(1);
t396 = -pkin(3) * t353 + t395 * t356;
t394 = r_i_i_C(2) - pkin(4);
t392 = r_i_i_C(3) + qJ(5);
t391 = cos(pkin(6));
t350 = sin(pkin(6));
t390 = t350 * t353;
t389 = t350 * t356;
t352 = sin(qJ(4));
t388 = t352 * t356;
t357 = cos(qJ(2));
t387 = t352 * t357;
t355 = cos(qJ(4));
t386 = t355 * t357;
t354 = sin(qJ(2));
t385 = qJD(2) * t354;
t384 = qJD(2) * t357;
t383 = qJD(3) * t353;
t382 = qJD(3) * t357;
t381 = qJD(4) * t356;
t380 = t350 * t385;
t379 = t350 * t384;
t378 = t354 * t391;
t377 = t357 * t391;
t376 = -qJD(2) + t381;
t349 = sin(pkin(10));
t375 = t349 * t378;
t351 = cos(pkin(10));
t374 = t351 * t377;
t336 = -qJD(2) * t374 + t349 * t385;
t340 = t349 * t354 - t374;
t373 = t340 * t381 - t336;
t342 = t349 * t377 + t351 * t354;
t338 = t342 * qJD(2);
t372 = t342 * t381 - t338;
t341 = t349 * t357 + t351 * t378;
t368 = -t341 * t356 + t351 * t390;
t371 = t340 * t352 - t355 * t368;
t343 = t351 * t357 - t375;
t333 = t343 * t356 + t349 * t390;
t370 = t333 * t355 + t342 * t352;
t369 = -t341 * t353 - t351 * t389;
t367 = -t343 * t353 + t349 * t389;
t345 = t391 * t353 + t354 * t389;
t366 = -t345 * t355 + t350 * t387;
t365 = -t356 * pkin(3) - t395 * t353 - pkin(2);
t364 = qJD(3) * t396;
t363 = -t354 * t390 + t391 * t356;
t362 = t392 * t352 - t394 * t355 + pkin(3);
t337 = t341 * qJD(2);
t361 = qJD(4) * t341 - t337 * t356 + t340 * t383;
t339 = -qJD(2) * t375 + t351 * t384;
t360 = qJD(4) * t343 - t339 * t356 + t342 * t383;
t359 = t353 * t382 + (qJD(2) * t356 - qJD(4)) * t354;
t358 = qJD(5) * t352 + (t394 * t352 + t392 * t355) * qJD(4);
t335 = t363 * qJD(3) + t356 * t379;
t329 = t367 * qJD(3) - t338 * t356;
t327 = t369 * qJD(3) - t336 * t356;
t322 = -t366 * qJD(4) + t335 * t352 - t355 * t380;
t316 = t370 * qJD(4) + t329 * t352 - t339 * t355;
t314 = t371 * qJD(4) + t327 * t352 - t337 * t355;
t1 = [0 -(t342 * t388 + t343 * t355) * qJD(5) - t338 * pkin(8) - t394 * (t372 * t352 + t360 * t355) + t392 * (t360 * t352 - t372 * t355) - t342 * t364 + t365 * t339, t395 * t329 + t358 * t367 + t362 * (-t333 * qJD(3) + t338 * t353) t370 * qJD(5) + t392 * (t329 * t355 + t339 * t352 + (-t333 * t352 + t342 * t355) * qJD(4)) + t394 * t316, t316, 0; 0 -(t340 * t388 + t341 * t355) * qJD(5) - t336 * pkin(8) - t394 * (t373 * t352 + t361 * t355) + t392 * (t361 * t352 - t373 * t355) - t340 * t364 + t365 * t337, t395 * t327 + t358 * t369 + t362 * (t368 * qJD(3) + t336 * t353) t371 * qJD(5) + t392 * (t327 * t355 + t337 * t352 + (t340 * t355 + t352 * t368) * qJD(4)) + t394 * t314, t314, 0; 0 (t394 * (t359 * t355 + t376 * t387) - t392 * (t359 * t352 - t376 * t386) - (t354 * t355 - t356 * t387) * qJD(5) + t396 * t382 + (t357 * pkin(8) + t365 * t354) * qJD(2)) * t350, t395 * t335 + t358 * t363 + t362 * (-t345 * qJD(3) - t353 * t379) -t366 * qJD(5) + t392 * (t352 * t380 + t335 * t355 + (-t345 * t352 - t350 * t386) * qJD(4)) + t394 * t322, t322, 0;];
JaD_transl  = t1;
