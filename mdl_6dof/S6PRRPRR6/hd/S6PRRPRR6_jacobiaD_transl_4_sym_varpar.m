% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR6
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR6_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR6_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiaD_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:08
% EndTime: 2019-02-26 20:07:08
% DurationCPUTime: 0.26s
% Computational Cost: add. (272->74), mult. (919->136), div. (0->0), fcn. (944->12), ass. (0->52)
t354 = sin(pkin(12));
t358 = cos(pkin(12));
t362 = sin(qJ(2));
t360 = cos(pkin(6));
t364 = cos(qJ(2));
t383 = t360 * t364;
t346 = -t354 * t362 + t358 * t383;
t361 = sin(qJ(3));
t363 = cos(qJ(3));
t384 = t360 * t362;
t369 = t354 * t384 - t358 * t364;
t359 = cos(pkin(7));
t370 = t354 * t383 + t358 * t362;
t355 = sin(pkin(7));
t356 = sin(pkin(6));
t388 = t355 * t356;
t371 = t354 * t388 - t359 * t370;
t395 = t371 * t361 - t363 * t369;
t394 = t355 * pkin(9);
t393 = r_i_i_C(3) + qJ(4);
t347 = t354 * t364 + t358 * t384;
t392 = t347 * t363;
t353 = sin(pkin(13));
t390 = t353 * t355;
t357 = cos(pkin(13));
t387 = t355 * t357;
t386 = t359 * t361;
t385 = t359 * t363;
t382 = t361 * t362;
t381 = t361 * t364;
t380 = t362 * t363;
t379 = t363 * t364;
t378 = qJD(3) * t355;
t376 = t361 * t378;
t375 = t357 * r_i_i_C(1) - t353 * r_i_i_C(2) + pkin(3);
t374 = -t346 * t361 - t347 * t385;
t373 = -t346 * t359 + t358 * t388;
t372 = t361 * t370 + t369 * t385;
t368 = t359 * t379 - t382;
t367 = t359 * t380 + t381;
t366 = t359 * t381 + t380;
t365 = t359 * t382 - t379;
t345 = t369 * qJD(2);
t344 = t370 * qJD(2);
t343 = t347 * qJD(2);
t342 = t346 * qJD(2);
t336 = t360 * t376 + (t367 * qJD(2) + t366 * qJD(3)) * t356;
t335 = t372 * qJD(3) + t344 * t386 + t345 * t363;
t333 = t374 * qJD(3) - t342 * t386 - t343 * t363;
t330 = t395 * qJD(3) - t344 * t361 - t345 * t385;
t328 = t342 * t361 + t343 * t385 - t358 * t356 * t376 + (t346 * t386 + t392) * qJD(3);
t1 = [0 (t335 * t357 - t344 * t390) * r_i_i_C(1) + (-t335 * t353 - t344 * t387) * r_i_i_C(2) + t335 * pkin(3) - t372 * qJD(4) + t345 * pkin(2) - t344 * t394 + t393 * (-t344 * t385 + t345 * t361 + (-t363 * t370 + t369 * t386) * qJD(3)) t395 * qJD(4) + t393 * (t345 * t386 - t344 * t363 + (t361 * t369 + t371 * t363) * qJD(3)) - t375 * t330, t330, 0, 0; 0 (t333 * t357 + t342 * t390) * r_i_i_C(1) + (-t333 * t353 + t342 * t387) * r_i_i_C(2) + t333 * pkin(3) - t374 * qJD(4) - t343 * pkin(2) + t342 * t394 + t393 * (t342 * t385 - t343 * t361 + (t346 * t363 - t347 * t386) * qJD(3)) -(t373 * t361 - t392) * qJD(4) + t393 * (-t343 * t386 + t342 * t363 + (-t347 * t361 - t373 * t363) * qJD(3)) - t375 * t328, t328, 0, 0; 0 (-t393 * (-t368 * qJD(2) + t365 * qJD(3)) + t375 * (-t366 * qJD(2) - t367 * qJD(3)) + t367 * qJD(4) + (-t362 * pkin(2) + (r_i_i_C(1) * t353 + r_i_i_C(2) * t357 + pkin(9)) * t364 * t355) * qJD(2)) * t356 -(-t360 * t355 * t361 - t366 * t356) * qJD(4) + t393 * (t360 * t363 * t378 + (-t365 * qJD(2) + t368 * qJD(3)) * t356) - t375 * t336, t336, 0, 0;];
JaD_transl  = t1;
