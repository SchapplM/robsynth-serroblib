% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:14
% EndTime: 2019-02-26 21:41:14
% DurationCPUTime: 0.48s
% Computational Cost: add. (687->90), mult. (1196->138), div. (0->0), fcn. (1090->10), ass. (0->63)
t398 = qJ(4) + pkin(10);
t355 = sin(t398);
t363 = cos(qJ(2));
t359 = sin(qJ(2));
t392 = cos(t398);
t388 = t359 * t392;
t336 = -t363 * t355 + t388;
t358 = sin(qJ(4));
t393 = pkin(4) * t358 + qJ(3);
t362 = cos(qJ(4));
t412 = t362 * pkin(4) + pkin(2) + pkin(3);
t374 = t412 * t359 - t393 * t363;
t424 = t374 * qJD(2) - t359 * qJD(3);
t386 = qJD(4) * t392;
t387 = qJD(2) * t392;
t402 = qJD(4) * t355;
t423 = t363 * t402 + (-t386 + t387) * t359;
t335 = t359 * t355 + t363 * t392;
t404 = qJD(2) * t359;
t329 = t335 * qJD(4) - t355 * t404 - t363 * t387;
t360 = sin(qJ(1));
t364 = cos(qJ(1));
t419 = t336 * qJD(1);
t327 = t329 * t360 - t419 * t364;
t403 = qJD(2) * t363;
t330 = t355 * t403 - t423;
t372 = qJD(1) * t335;
t328 = t330 * t360 + t364 * t372;
t357 = sin(qJ(6));
t361 = cos(qJ(6));
t390 = -t361 * r_i_i_C(1) + t357 * r_i_i_C(2);
t383 = pkin(5) - t390;
t414 = -r_i_i_C(3) - pkin(9);
t413 = t361 * r_i_i_C(2);
t418 = qJD(6) * (t357 * r_i_i_C(1) + t413);
t422 = -t336 * t360 * t418 - t383 * t327 - t414 * t328;
t394 = t364 * t403;
t325 = -t355 * t394 + t360 * t372 + t364 * t423;
t334 = t335 * t364;
t407 = t363 * t364;
t326 = -t364 * t359 * t402 + qJD(2) * t334 - t419 * t360 - t386 * t407;
t421 = (t355 * t407 - t364 * t388) * t418 + t414 * t325 + t383 * t326;
t420 = t414 * t329 - t383 * t330 + t335 * t418;
t411 = pkin(7) - qJ(5) - pkin(8);
t410 = pkin(4) * qJD(4);
t409 = t328 * t357;
t406 = qJD(1) * t360;
t405 = qJD(1) * t364;
t401 = qJD(6) * t357;
t400 = qJD(6) * t361;
t391 = -t328 * t361 + t357 * t406;
t385 = t358 * t363 - t359 * t362;
t384 = t358 * t359 + t362 * t363;
t376 = t385 * qJD(4);
t375 = -t393 * t359 - t412 * t363;
t368 = -qJD(5) + (-pkin(1) + t375) * qJD(1);
t367 = (qJD(2) - qJD(4)) * t384;
t366 = t375 * qJD(2) + qJD(3) * t363 + t384 * t410;
t365 = -t385 * t410 - t424;
t332 = t335 * t360;
t324 = -t357 * t405 - t325 * t361 + (-t334 * t357 - t360 * t361) * qJD(6);
t323 = -t361 * t405 + t325 * t357 + (-t334 * t361 + t357 * t360) * qJD(6);
t1 = [(t332 * t401 + t391) * r_i_i_C(1) + (t332 * t400 + t409) * r_i_i_C(2) - t328 * pkin(5) + t414 * t327 + (t390 * qJD(6) + t368) * t364 + (pkin(4) * t376 + (-t411 + t413) * qJD(1) + t424) * t360, t366 * t364 + t374 * t406 - t421, -t359 * t406 + t394 (t367 * t364 + t385 * t406) * pkin(4) + t421, -t405, t323 * r_i_i_C(1) - t324 * r_i_i_C(2); -t325 * pkin(5) + t324 * r_i_i_C(1) + t323 * r_i_i_C(2) + t414 * t326 + t368 * t360 + (t411 * qJD(1) + t365) * t364, t366 * t360 - t374 * t405 - t422, t359 * t405 + t360 * t403 (t367 * t360 - t385 * t405) * pkin(4) + t422, -t406 (-t361 * t406 - t409) * r_i_i_C(1) + t391 * r_i_i_C(2) + ((-t332 * t361 - t357 * t364) * r_i_i_C(1) + (t332 * t357 - t361 * t364) * r_i_i_C(2)) * qJD(6); 0, t365 - t420, t404 (-t385 * qJD(2) + t376) * pkin(4) + t420, 0 (t329 * t361 + t336 * t401) * r_i_i_C(2) + (t329 * t357 - t336 * t400) * r_i_i_C(1);];
JaD_transl  = t1;
