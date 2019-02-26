% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:08
% EndTime: 2019-02-26 21:28:08
% DurationCPUTime: 0.43s
% Computational Cost: add. (687->80), mult. (1076->125), div. (0->0), fcn. (998->10), ass. (0->58)
t357 = qJ(2) + pkin(10);
t356 = cos(t357);
t360 = sin(qJ(5));
t355 = sin(t357);
t407 = cos(qJ(5));
t392 = t355 * t407;
t336 = -t356 * t360 + t392;
t409 = pkin(3) + pkin(4);
t371 = t409 * t355 + sin(qJ(2)) * pkin(2) - qJ(4) * t356;
t366 = -t371 * qJD(2) + t355 * qJD(4);
t362 = sin(qJ(1));
t365 = cos(qJ(1));
t387 = qJD(2) * t407;
t383 = t362 * t387;
t414 = t336 * qJD(1);
t335 = t355 * t360 + t356 * t407;
t397 = qJD(2) * t360;
t416 = qJD(5) * t335 - t355 * t397;
t327 = -t356 * t383 + t362 * t416 - t365 * t414;
t386 = qJD(5) * t407;
t382 = t355 * t386;
t369 = t356 * t397 + t382;
t374 = qJD(1) * t335;
t396 = qJD(5) * t360;
t388 = t356 * t396;
t328 = -t355 * t383 + t365 * t374 + (t369 - t388) * t362;
t359 = sin(qJ(6));
t363 = cos(qJ(6));
t381 = -r_i_i_C(1) * t363 + r_i_i_C(2) * t359;
t378 = pkin(5) - t381;
t408 = -r_i_i_C(3) - pkin(9);
t405 = r_i_i_C(2) * t363;
t413 = qJD(6) * (r_i_i_C(1) * t359 + t405);
t419 = -t336 * t362 * t413 - t327 * t378 - t328 * t408;
t370 = t355 * t387 + t388;
t398 = qJD(2) * t356;
t389 = t365 * t398;
t325 = -t360 * t389 + t362 * t374 + (t370 - t382) * t365;
t334 = t335 * t365;
t401 = t356 * t365;
t326 = -t365 * t355 * t396 + qJD(2) * t334 - t362 * t414 - t386 * t401;
t418 = (t360 * t401 - t365 * t392) * t413 + t408 * t325 + t378 * t326;
t329 = -t356 * t387 + t416;
t417 = t335 * t413 + t408 * t329 - t378 * (t369 - t370);
t415 = -qJ(4) * t355 - t409 * t356 - cos(qJ(2)) * pkin(2);
t404 = -pkin(8) + qJ(3) + pkin(7);
t403 = t328 * t359;
t400 = qJD(1) * t362;
t399 = qJD(1) * t365;
t395 = qJD(6) * t359;
t394 = qJD(6) * t363;
t385 = -t328 * t363 + t359 * t400;
t368 = qJD(3) + (-pkin(1) + t415) * qJD(1);
t367 = qJD(2) * t415 + qJD(4) * t356;
t332 = t335 * t362;
t324 = -t359 * t399 - t325 * t363 + (-t334 * t359 - t362 * t363) * qJD(6);
t323 = -t363 * t399 + t325 * t359 + (-t334 * t363 + t359 * t362) * qJD(6);
t1 = [(t332 * t395 + t385) * r_i_i_C(1) + (t332 * t394 + t403) * r_i_i_C(2) - t328 * pkin(5) + t408 * t327 + (qJD(6) * t381 + t368) * t365 + ((-t404 + t405) * qJD(1) - t366) * t362, t365 * t367 + t371 * t400 - t418, t399, -t355 * t400 + t389, t418, r_i_i_C(1) * t323 - r_i_i_C(2) * t324; -t325 * pkin(5) + t324 * r_i_i_C(1) + t323 * r_i_i_C(2) + t408 * t326 + t368 * t362 + (qJD(1) * t404 + t366) * t365, t362 * t367 - t371 * t399 - t419, t400, t355 * t399 + t362 * t398, t419 (-t363 * t400 - t403) * r_i_i_C(1) + t385 * r_i_i_C(2) + ((-t332 * t363 - t359 * t365) * r_i_i_C(1) + (t332 * t359 - t363 * t365) * r_i_i_C(2)) * qJD(6); 0, t366 - t417, 0, qJD(2) * t355, t417 (t329 * t363 + t336 * t395) * r_i_i_C(2) + (t329 * t359 - t336 * t394) * r_i_i_C(1);];
JaD_transl  = t1;
