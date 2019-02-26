% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR15_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR15_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:22
% EndTime: 2019-02-26 22:24:22
% DurationCPUTime: 0.42s
% Computational Cost: add. (470->85), mult. (1462->141), div. (0->0), fcn. (1486->10), ass. (0->59)
t398 = cos(pkin(6));
t400 = sin(qJ(2));
t401 = sin(qJ(1));
t428 = t401 * t400;
t421 = t398 * t428;
t403 = cos(qJ(2));
t404 = cos(qJ(1));
t424 = t404 * t403;
t384 = -qJD(1) * t421 - qJD(2) * t428 + (qJD(2) * t398 + qJD(1)) * t424;
t399 = sin(qJ(3));
t402 = cos(qJ(3));
t425 = t404 * t400;
t427 = t401 * t403;
t387 = t398 * t425 + t427;
t386 = -t398 * t424 + t428;
t395 = sin(pkin(7));
t397 = cos(pkin(7));
t396 = sin(pkin(6));
t434 = t396 * t404;
t416 = t386 * t397 + t395 * t434;
t405 = -t387 * t402 + t416 * t399;
t409 = t398 * t427 + t425;
t383 = t409 * qJD(1) + t387 * qJD(2);
t420 = qJD(1) * t395 * t396;
t407 = -t383 * t397 + t401 * t420;
t368 = -t405 * qJD(3) + t384 * t399 - t407 * t402;
t406 = t387 * t399 + t416 * t402;
t453 = -t406 * qJD(3) + t384 * t402 + t407 * t399;
t445 = pkin(10) + r_i_i_C(1);
t452 = t445 * t395;
t449 = t445 * t397 + pkin(9);
t408 = t421 - t424;
t435 = t396 * t401;
t414 = t395 * t435 - t397 * t409;
t447 = t414 * t399 - t402 * t408;
t444 = -r_i_i_C(2) + pkin(3);
t443 = r_i_i_C(3) + qJ(4);
t438 = t408 * t399;
t436 = t395 * t398;
t433 = t397 * t399;
t432 = t397 * t402;
t431 = t399 * t400;
t430 = t399 * t403;
t429 = t400 * t402;
t426 = t402 * t403;
t419 = qJD(3) * t436;
t418 = t404 * t420;
t417 = t386 * t399 - t387 * t432;
t415 = t399 * t409 + t408 * t432;
t413 = t397 * t426 - t431;
t412 = t397 * t429 + t430;
t411 = t397 * t430 + t429;
t410 = t397 * t431 - t426;
t382 = t387 * qJD(1) + t409 * qJD(2);
t381 = t386 * qJD(1) + t408 * qJD(2);
t376 = t399 * t419 + (t412 * qJD(2) + t411 * qJD(3)) * t396;
t367 = qJD(3) * t438 + (t381 * t397 + t418) * t399 + (t414 * qJD(3) - t382) * t402;
t366 = t447 * qJD(3) - t381 * t432 - t382 * t399 - t402 * t418;
t1 = [-t406 * qJD(4) - t384 * pkin(2) - t383 * t452 - t444 * t453 - t443 * t368 + (-t404 * pkin(1) - t449 * t435) * qJD(1), -t415 * qJD(4) + t381 * pkin(2) - t382 * t452 + t444 * (t415 * qJD(3) + t381 * t402 + t382 * t433) + t443 * (-t382 * t432 + t381 * t399 + (-t402 * t409 + t408 * t433) * qJD(3)) t447 * qJD(4) - t444 * t366 + t443 * t367, t366, 0, 0; -(t414 * t402 + t438) * qJD(4) - t382 * pkin(2) - t381 * t452 + t444 * t367 + t443 * t366 + (-t401 * pkin(1) + t449 * t434) * qJD(1), -t417 * qJD(4) - t383 * pkin(2) + t384 * t452 + t444 * (t417 * qJD(3) - t383 * t402 - t384 * t433) + t443 * (t384 * t432 - t383 * t399 + (-t386 * t402 - t387 * t433) * qJD(3)) -t405 * qJD(4) - t444 * t368 + t443 * t453, t368, 0, 0; 0 (-t444 * (t411 * qJD(2) + t412 * qJD(3)) - t443 * (-t413 * qJD(2) + t410 * qJD(3)) + t412 * qJD(4) + (-t400 * pkin(2) + t403 * t452) * qJD(2)) * t396 -(-t411 * t396 - t399 * t436) * qJD(4) + t443 * (t402 * t419 + (-t410 * qJD(2) + t413 * qJD(3)) * t396) - t444 * t376, t376, 0, 0;];
JaD_transl  = t1;
