% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:55:24
% EndTime: 2019-02-26 19:55:25
% DurationCPUTime: 0.44s
% Computational Cost: add. (677->92), mult. (1118->152), div. (0->0), fcn. (1125->13), ass. (0->65)
t346 = pkin(12) + qJ(4);
t342 = sin(t346);
t343 = cos(t346);
t354 = cos(qJ(5));
t348 = qJ(5) + qJ(6);
t344 = sin(t348);
t345 = cos(t348);
t370 = t345 * r_i_i_C(1) - t344 * r_i_i_C(2);
t366 = t354 * pkin(5) + pkin(4) + t370;
t395 = r_i_i_C(3) + pkin(10) + pkin(9);
t401 = (t366 * t342 - t395 * t343) * qJD(4);
t353 = sin(qJ(2));
t355 = cos(qJ(2));
t349 = sin(pkin(11));
t394 = cos(pkin(6));
t379 = t349 * t394;
t393 = cos(pkin(11));
t336 = -t353 * t379 + t393 * t355;
t369 = t344 * r_i_i_C(1) + t345 * r_i_i_C(2);
t347 = qJD(5) + qJD(6);
t352 = sin(qJ(5));
t396 = t352 * pkin(5);
t400 = qJD(5) * t396 + t347 * t369;
t392 = t344 * t347;
t391 = t345 * t347;
t350 = sin(pkin(6));
t390 = t349 * t350;
t389 = t350 * t353;
t388 = t350 * t355;
t367 = t394 * t393;
t334 = t349 * t355 + t353 * t367;
t330 = t334 * qJD(2);
t378 = t350 * t393;
t360 = -t334 * t343 + t342 * t378;
t374 = -t347 * t360 - t330;
t365 = t355 * t367;
t384 = qJD(2) * t353;
t329 = -qJD(2) * t365 + t349 * t384;
t362 = -t334 * t342 - t343 * t378;
t318 = qJD(4) * t362 - t329 * t343;
t333 = t349 * t353 - t365;
t376 = -t333 * t347 - t318;
t387 = (t376 * t344 - t374 * t345) * r_i_i_C(1) + (t374 * t344 + t376 * t345) * r_i_i_C(2);
t326 = t336 * t343 + t342 * t390;
t332 = t336 * qJD(2);
t373 = t326 * t347 - t332;
t335 = t393 * t353 + t355 * t379;
t331 = t335 * qJD(2);
t364 = -t336 * t342 + t343 * t390;
t320 = qJD(4) * t364 - t331 * t343;
t375 = -t335 * t347 - t320;
t386 = (t375 * t344 - t373 * t345) * r_i_i_C(1) + (t373 * t344 + t375 * t345) * r_i_i_C(2);
t328 = t394 * t342 + t343 * t389;
t380 = t350 * t384;
t363 = -t328 * t347 + t380;
t361 = -t342 * t389 + t394 * t343;
t381 = qJD(2) * t388;
t322 = qJD(4) * t361 + t343 * t381;
t368 = t347 * t388 - t322;
t385 = (t344 * t368 + t345 * t363) * r_i_i_C(1) + (-t344 * t363 + t345 * t368) * r_i_i_C(2);
t383 = qJD(5) * t354;
t358 = -t395 * t342 - t366 * t343 - cos(pkin(12)) * pkin(3) - pkin(2);
t357 = t400 * t343 + t401;
t351 = -pkin(8) - qJ(3);
t1 = [0 (-t331 * t344 + t336 * t391) * r_i_i_C(1) + (-t331 * t345 - t336 * t392) * r_i_i_C(2) + t331 * t351 + t336 * qJD(3) + (-t331 * t352 + t336 * t383) * pkin(5) + t358 * t332 + t357 * t335, t332, t395 * t320 - t400 * t364 + t366 * (-qJD(4) * t326 + t331 * t342) (-t320 * t352 + t332 * t354 + (-t326 * t354 - t335 * t352) * qJD(5)) * pkin(5) + t386, t386; 0 (-t329 * t344 + t334 * t391) * r_i_i_C(1) + (-t329 * t345 - t334 * t392) * r_i_i_C(2) + t329 * t351 + t334 * qJD(3) + (-t329 * t352 + t334 * t383) * pkin(5) + t358 * t330 + t357 * t333, t330, t395 * t318 - t400 * t362 + t366 * (qJD(4) * t360 + t329 * t342) (-t318 * t352 + t330 * t354 + (-t333 * t352 + t354 * t360) * qJD(5)) * pkin(5) + t387, t387; 0 ((pkin(5) * t383 + t358 * qJD(2) + t370 * t347 + qJD(3)) * t353 + (-qJD(2) * t351 + (-qJD(5) * t343 + qJD(2)) * t396 - t401 + t369 * (-t343 * t347 + qJD(2))) * t355) * t350, t380, t395 * t322 - t400 * t361 + t366 * (-qJD(4) * t328 - t342 * t381) (t354 * t380 - t322 * t352 + (-t328 * t354 + t352 * t388) * qJD(5)) * pkin(5) + t385, t385;];
JaD_transl  = t1;
