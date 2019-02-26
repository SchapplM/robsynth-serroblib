% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:57:00
% EndTime: 2019-02-26 21:57:00
% DurationCPUTime: 0.56s
% Computational Cost: add. (986->90), mult. (1471->141), div. (0->0), fcn. (1350->10), ass. (0->69)
t389 = sin(qJ(2));
t390 = sin(qJ(1));
t394 = cos(qJ(1));
t393 = cos(qJ(2));
t443 = qJ(4) + qJ(5);
t433 = sin(t443);
t434 = cos(t443);
t367 = t389 * t433 + t393 * t434;
t405 = qJD(1) * t367;
t436 = qJD(4) + qJD(5);
t415 = t436 * t434;
t410 = t394 * t415;
t414 = t436 * t433;
t428 = qJD(2) * t433;
t429 = qJD(2) * t434;
t453 = t389 * t429 + (t414 - t428) * t393;
t354 = -t389 * t410 + t390 * t405 + t453 * t394;
t366 = t367 * t394;
t412 = t389 * t414;
t368 = t389 * t434 - t393 * t433;
t450 = t368 * qJD(1);
t355 = qJD(2) * t366 - t450 * t390 - t393 * t410 - t394 * t412;
t449 = -r_i_i_C(3) - pkin(10);
t387 = sin(qJ(6));
t391 = cos(qJ(6));
t451 = -t391 * r_i_i_C(1) + t387 * r_i_i_C(2) - pkin(5);
t460 = t449 * t354 - t451 * t355;
t361 = -t389 * t428 + t412 + (t415 - t429) * t393;
t356 = t361 * t390 - t450 * t394;
t362 = t389 * t415 - t453;
t357 = t362 * t390 + t394 * t405;
t459 = t451 * t356 - t449 * t357;
t458 = t449 * t361 + t451 * t362;
t388 = sin(qJ(4));
t432 = pkin(4) * t388 + qJ(3);
t392 = cos(qJ(4));
t446 = t392 * pkin(4) + pkin(2) + pkin(3);
t408 = t446 * t389 - t432 * t393;
t454 = t408 * qJD(2) - t389 * qJD(3);
t438 = qJD(6) * t391;
t439 = qJD(6) * t387;
t452 = r_i_i_C(1) * t439 + t438 * r_i_i_C(2);
t447 = t391 * r_i_i_C(2);
t445 = pkin(7) - pkin(9) - pkin(8);
t444 = pkin(4) * qJD(4);
t442 = qJD(1) * t390;
t441 = qJD(1) * t394;
t440 = qJD(2) * t393;
t364 = t367 * t390;
t427 = t364 * t391 + t387 * t394;
t426 = t364 * t387 - t391 * t394;
t425 = t388 * t393 - t389 * t392;
t424 = t388 * t389 + t392 * t393;
t416 = qJD(6) * (-t387 * r_i_i_C(1) - t447);
t413 = t425 * qJD(4);
t409 = -t432 * t389 - t446 * t393;
t406 = qJD(1) * (-pkin(1) + t409);
t402 = (qJD(2) - qJD(4)) * t424;
t363 = t368 * t390;
t400 = -t452 * t363 + t459;
t365 = t368 * t394;
t399 = -t452 * t365 + t460;
t398 = t452 * t367 + t458;
t397 = t409 * qJD(2) + qJD(3) * t393 + t424 * t444;
t396 = -t425 * t444 - t454;
t385 = t387 * t442;
t347 = -t387 * t441 - t354 * t391 + (-t366 * t387 - t390 * t391) * qJD(6);
t346 = -t391 * t441 + t354 * t387 + (-t366 * t391 + t387 * t390) * qJD(6);
t1 = [t385 * r_i_i_C(1) + t451 * t357 + t449 * t356 + (t426 * r_i_i_C(1) + t427 * r_i_i_C(2)) * qJD(6) + t394 * t406 + (pkin(4) * t413 + (-t445 + t447) * qJD(1) + t454) * t390, -t365 * t416 + t397 * t394 + t408 * t442 - t460, -t389 * t442 + t394 * t440 (t402 * t394 + t425 * t442) * pkin(4) + t399, t399, t346 * r_i_i_C(1) - t347 * r_i_i_C(2); -t354 * pkin(5) + t347 * r_i_i_C(1) + t346 * r_i_i_C(2) + t449 * t355 + t390 * t406 + (t445 * qJD(1) + t396) * t394, -t363 * t416 + t397 * t390 - t408 * t441 - t459, t389 * t441 + t390 * t440 (t402 * t390 - t425 * t441) * pkin(4) + t400, t400 (-t357 * t387 - t391 * t442) * r_i_i_C(1) + (-t357 * t391 + t385) * r_i_i_C(2) + (-t427 * r_i_i_C(1) + t426 * r_i_i_C(2)) * qJD(6); 0, t367 * t416 + t396 - t458, qJD(2) * t389 (-t425 * qJD(2) + t413) * pkin(4) + t398, t398 (t361 * t391 + t368 * t439) * r_i_i_C(2) + (t361 * t387 - t368 * t438) * r_i_i_C(1);];
JaD_transl  = t1;
