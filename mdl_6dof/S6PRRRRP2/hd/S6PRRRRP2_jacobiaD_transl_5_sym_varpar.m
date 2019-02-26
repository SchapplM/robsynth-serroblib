% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:44
% EndTime: 2019-02-26 20:15:44
% DurationCPUTime: 0.45s
% Computational Cost: add. (593->88), mult. (1038->158), div. (0->0), fcn. (1028->12), ass. (0->61)
t406 = pkin(10) + r_i_i_C(3);
t362 = sin(qJ(5));
t365 = cos(qJ(5));
t379 = r_i_i_C(1) * t365 - r_i_i_C(2) * t362;
t412 = pkin(4) + t379;
t392 = qJD(5) * t365;
t393 = qJD(5) * t362;
t411 = -r_i_i_C(1) * t393 - t392 * r_i_i_C(2);
t359 = qJ(3) + qJ(4);
t356 = sin(t359);
t357 = cos(t359);
t358 = qJD(3) + qJD(4);
t363 = sin(qJ(3));
t410 = -qJD(3) * t363 * pkin(3) - (t356 * t412 - t406 * t357) * t358;
t364 = sin(qJ(2));
t367 = cos(qJ(2));
t360 = sin(pkin(11));
t402 = cos(pkin(6));
t385 = t360 * t402;
t401 = cos(pkin(11));
t349 = -t364 * t385 + t401 * t367;
t408 = t362 * r_i_i_C(1) + t365 * r_i_i_C(2);
t400 = t356 * t358;
t399 = t357 * t358;
t361 = sin(pkin(6));
t398 = t360 * t361;
t397 = t361 * t364;
t396 = t361 * t367;
t395 = qJD(2) * t364;
t394 = qJD(5) * t357;
t389 = t356 * t397;
t388 = t357 * t397;
t387 = t361 * t395;
t386 = qJD(2) * t396;
t384 = t361 * t401;
t380 = t357 * t384;
t348 = t401 * t364 + t367 * t385;
t344 = t348 * qJD(2);
t378 = t358 * t398 - t344;
t377 = t402 * t401;
t375 = t367 * t377;
t374 = t402 * t358 + t386;
t347 = t360 * t367 + t364 * t377;
t366 = cos(qJ(3));
t373 = -pkin(3) * t366 - t406 * t356 - t357 * t412 - pkin(2);
t328 = -t349 * t400 + t378 * t357;
t372 = t411 * (-t349 * t356 + t357 * t398) + t406 * t328 + t412 * (-t349 * t399 - t378 * t356);
t342 = -qJD(2) * t375 + t360 * t395;
t326 = -t342 * t357 - t347 * t400 - t358 * t380;
t371 = t411 * (-t347 * t356 - t380) + t406 * t326 + t412 * (-t347 * t399 + (t358 * t384 + t342) * t356);
t333 = t374 * t357 - t358 * t389;
t370 = t411 * (t402 * t357 - t389) + t406 * t333 + t412 * (-t374 * t356 - t358 * t388);
t369 = t408 * t394 - t410;
t368 = -pkin(9) - pkin(8);
t346 = t360 * t364 - t375;
t345 = t349 * qJD(2);
t343 = t347 * qJD(2);
t341 = t402 * t356 + t388;
t337 = t349 * t357 + t356 * t398;
t335 = t347 * t357 - t356 * t384;
t1 = [0 (-t344 * t362 + t349 * t392) * r_i_i_C(1) + (-t344 * t365 - t349 * t393) * r_i_i_C(2) + t344 * t368 + t373 * t345 + t369 * t348 (t344 * t363 + (-t349 * t366 - t363 * t398) * qJD(3)) * pkin(3) + t372, t372 (-t328 * t362 + t345 * t365) * r_i_i_C(1) + (-t328 * t365 - t345 * t362) * r_i_i_C(2) + ((-t337 * t365 - t348 * t362) * r_i_i_C(1) + (t337 * t362 - t348 * t365) * r_i_i_C(2)) * qJD(5), 0; 0 (-t342 * t362 + t347 * t392) * r_i_i_C(1) + (-t342 * t365 - t347 * t393) * r_i_i_C(2) + t342 * t368 + t373 * t343 + t369 * t346 (t342 * t363 + (-t347 * t366 + t363 * t384) * qJD(3)) * pkin(3) + t371, t371 (-t326 * t362 + t343 * t365) * r_i_i_C(1) + (-t326 * t365 - t343 * t362) * r_i_i_C(2) + ((-t335 * t365 - t346 * t362) * r_i_i_C(1) + (t335 * t362 - t346 * t365) * r_i_i_C(2)) * qJD(5), 0; 0 ((t373 * qJD(2) + t379 * qJD(5)) * t364 + (-qJD(2) * t368 + t408 * (qJD(2) - t394) + t410) * t367) * t361 (-t363 * t386 + (-t402 * t363 - t366 * t397) * qJD(3)) * pkin(3) + t370, t370 (-t333 * t362 + t365 * t387) * r_i_i_C(1) + (-t333 * t365 - t362 * t387) * r_i_i_C(2) + ((-t341 * t365 + t362 * t396) * r_i_i_C(1) + (t341 * t362 + t365 * t396) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
