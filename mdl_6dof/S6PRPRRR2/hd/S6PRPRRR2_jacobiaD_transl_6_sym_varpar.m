% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR2
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
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:13
% EndTime: 2019-02-26 19:54:14
% DurationCPUTime: 0.62s
% Computational Cost: add. (787->97), mult. (2002->173), div. (0->0), fcn. (2172->14), ass. (0->70)
t422 = sin(qJ(4));
t425 = cos(qJ(4));
t415 = qJ(5) + qJ(6);
t412 = sin(t415);
t413 = cos(t415);
t424 = cos(qJ(5));
t438 = t424 * pkin(5) + r_i_i_C(1) * t413 - r_i_i_C(2) * t412 + pkin(4);
t467 = r_i_i_C(3) + pkin(10) + pkin(9);
t414 = qJD(5) + qJD(6);
t421 = sin(qJ(5));
t470 = qJD(5) * t421 * pkin(5) + (t412 * r_i_i_C(1) + t413 * r_i_i_C(2)) * t414;
t473 = t470 * t425 + (t438 * t422 - t467 * t425) * qJD(4);
t420 = cos(pkin(6));
t416 = sin(pkin(12));
t423 = sin(qJ(2));
t465 = cos(pkin(12));
t468 = cos(qJ(2));
t433 = t468 * t416 + t423 * t465;
t400 = t433 * t420;
t443 = t468 * t465;
t455 = qJD(2) * t423;
t471 = -qJD(2) * t443 + t416 * t455;
t432 = -t423 * t416 + t443;
t428 = t467 * t422 + t438 * t425 + pkin(3);
t466 = pkin(2) * qJD(2);
t464 = t412 * t414;
t463 = t413 * t414;
t418 = sin(pkin(6));
t462 = t418 * t422;
t461 = t418 * t425;
t460 = t420 * t423;
t401 = t432 * qJD(2);
t417 = sin(pkin(11));
t419 = cos(pkin(11));
t429 = qJD(2) * t400;
t378 = -t417 * t401 - t419 * t429;
t440 = t419 * t400 + t417 * t432;
t435 = t419 * t462 - t425 * t440;
t446 = -t414 * t435 + t378;
t397 = t471 * t420;
t402 = t433 * qJD(2);
t380 = t419 * t397 + t417 * t402;
t436 = -t419 * t461 - t422 * t440;
t369 = t436 * qJD(4) - t380 * t425;
t431 = t432 * t420;
t385 = -t417 * t433 + t419 * t431;
t449 = t385 * t414 - t369;
t458 = (t449 * t412 - t446 * t413) * r_i_i_C(1) + (t446 * t412 + t449 * t413) * r_i_i_C(2);
t439 = -t417 * t400 + t419 * t432;
t377 = t417 * t462 + t425 * t439;
t381 = -t419 * t401 + t417 * t429;
t445 = t377 * t414 + t381;
t382 = t417 * t397 - t419 * t402;
t434 = t417 * t461 - t422 * t439;
t371 = t434 * qJD(4) + t382 * t425;
t388 = -t417 * t431 - t419 * t433;
t448 = t388 * t414 - t371;
t457 = (t448 * t412 - t445 * t413) * r_i_i_C(1) + (t445 * t412 + t448 * t413) * r_i_i_C(2);
t399 = t433 * t418;
t391 = t399 * t425 + t420 * t422;
t396 = qJD(2) * t399;
t444 = t391 * t414 - t396;
t395 = t471 * t418;
t441 = -t399 * t422 + t420 * t425;
t373 = t441 * qJD(4) - t395 * t425;
t398 = t432 * t418;
t447 = t398 * t414 - t373;
t456 = (t447 * t412 - t444 * t413) * r_i_i_C(1) + (t444 * t412 + t447 * t413) * r_i_i_C(2);
t454 = qJD(5) * t424;
t1 = [0 (t382 * t412 + t439 * t463) * r_i_i_C(1) + (t382 * t413 - t439 * t464) * r_i_i_C(2) + t382 * pkin(8) + (t382 * t421 + t439 * t454) * pkin(5) + (t417 * t460 - t468 * t419) * t466 + t428 * t381 - t473 * t388, 0, t467 * t371 - t470 * t434 + t438 * (-t377 * qJD(4) - t382 * t422) (-t371 * t421 - t381 * t424 + (-t377 * t424 + t388 * t421) * qJD(5)) * pkin(5) + t457, t457; 0 (-t380 * t412 + t440 * t463) * r_i_i_C(1) + (-t380 * t413 - t440 * t464) * r_i_i_C(2) - t380 * pkin(8) + (-t380 * t421 + t440 * t454) * pkin(5) + (-t468 * t417 - t419 * t460) * t466 + t428 * t378 - t473 * t385, 0, t467 * t369 - t470 * t436 + t438 * (t435 * qJD(4) + t380 * t422) (-t369 * t421 - t378 * t424 + (t385 * t421 + t424 * t435) * qJD(5)) * pkin(5) + t458, t458; 0 (-t395 * t412 + t399 * t463) * r_i_i_C(1) + (-t395 * t413 - t399 * t464) * r_i_i_C(2) - t395 * pkin(8) - t418 * pkin(2) * t455 + (-t395 * t421 + t399 * t454) * pkin(5) - t428 * t396 - t473 * t398, 0, t467 * t373 - t470 * t441 + t438 * (-t391 * qJD(4) + t395 * t422) (-t373 * t421 + t396 * t424 + (-t391 * t424 + t398 * t421) * qJD(5)) * pkin(5) + t456, t456;];
JaD_transl  = t1;
