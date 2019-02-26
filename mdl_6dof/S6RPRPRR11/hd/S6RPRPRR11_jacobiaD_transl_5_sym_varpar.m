% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR11_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR11_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:33
% EndTime: 2019-02-26 20:54:33
% DurationCPUTime: 0.44s
% Computational Cost: add. (502->90), mult. (1462->155), div. (0->0), fcn. (1557->14), ass. (0->73)
t420 = cos(pkin(12));
t422 = cos(pkin(6));
t417 = sin(pkin(12));
t425 = sin(qJ(1));
t456 = t425 * t417;
t446 = t422 * t456;
t427 = cos(qJ(1));
t451 = qJD(1) * t427;
t397 = -qJD(1) * t446 + t420 * t451;
t418 = sin(pkin(7));
t419 = sin(pkin(6));
t459 = t419 * t427;
t445 = t418 * t459;
t482 = -qJD(3) * t445 + t397;
t421 = cos(pkin(7));
t453 = t427 * t420;
t434 = t422 * t453 - t456;
t476 = t434 * t421;
t481 = -t476 + t445;
t454 = t427 * t417;
t455 = t425 * t420;
t433 = t422 * t455 + t454;
t396 = t433 * qJD(1);
t400 = t422 * t454 + t455;
t424 = sin(qJ(3));
t426 = cos(qJ(3));
t460 = t419 * t425;
t443 = qJD(1) * t460;
t375 = (-qJD(3) * t400 - t396 * t421 + t418 * t443) * t424 + (qJD(3) * t476 + t482) * t426;
t466 = t396 * t418;
t387 = t421 * t443 + t466;
t415 = pkin(13) + qJ(5);
t413 = sin(t415);
t414 = cos(t415);
t480 = t375 * t413 - t387 * t414;
t479 = -t375 * t414 - t387 * t413;
t475 = t434 * t418 + t421 * t459;
t457 = t421 * t426;
t462 = t418 * t422;
t474 = (-t417 * t424 + t420 * t457) * t419 + t426 * t462;
t458 = t421 * t424;
t437 = -t400 * t426 - t434 * t458;
t380 = t424 * t445 + t437;
t473 = t380 * t414 + t413 * t475;
t472 = -t380 * t413 + t414 * t475;
t470 = t400 * t424 + t481 * t426;
t452 = qJD(1) * t426;
t461 = t419 * t418;
t442 = t452 * t461;
t469 = t437 * qJD(3) - t396 * t457 - t482 * t424 + t425 * t442;
t468 = sin(pkin(13)) * pkin(4);
t467 = r_i_i_C(3) + pkin(10) + qJ(4);
t449 = t419 * qJD(2);
t440 = t413 * r_i_i_C(1) + t414 * r_i_i_C(2);
t412 = cos(pkin(13)) * pkin(4) + pkin(3);
t438 = -t414 * r_i_i_C(1) + t413 * r_i_i_C(2) - t412;
t435 = t418 * t460 - t421 * t433;
t432 = qJD(5) * t440;
t402 = -t446 + t453;
t381 = -t402 * t424 + t435 * t426;
t382 = t402 * t426 + t435 * t424;
t389 = t424 * t462 + (t417 * t426 + t420 * t458) * t419;
t398 = -t420 * t461 + t422 * t421;
t395 = t400 * qJD(1);
t392 = t418 * t433 + t421 * t460;
t385 = t475 * qJD(1);
t384 = t389 * qJD(3);
t383 = t474 * qJD(3);
t373 = t481 * t424 * qJD(1) + t381 * qJD(3) - t395 * t426;
t372 = t382 * qJD(3) - t395 * t424 - t427 * t442 + t452 * t476;
t371 = t373 * t414 + t385 * t413 + (-t382 * t413 + t392 * t414) * qJD(5);
t370 = -t373 * t413 + t385 * t414 + (-t382 * t414 - t392 * t413) * qJD(5);
t1 = [t479 * r_i_i_C(1) + t480 * r_i_i_C(2) - t375 * t412 - t470 * qJD(4) - t387 * t468 - t397 * pkin(2) - pkin(9) * t466 + t427 * t449 + t467 * t469 + (t472 * r_i_i_C(1) - t473 * r_i_i_C(2)) * qJD(5) + (-t427 * pkin(1) + (-pkin(9) * t421 - qJ(2)) * t460) * qJD(1), t419 * t451, t382 * qJD(4) + t438 * t372 + t467 * t373 - t381 * t432, t372, t370 * r_i_i_C(1) - t371 * r_i_i_C(2), 0; t385 * t468 + t425 * t449 - t395 * pkin(2) + t371 * r_i_i_C(1) + t370 * r_i_i_C(2) - t381 * qJD(4) + t373 * t412 + t467 * t372 + (-t425 * pkin(1) + pkin(9) * t475 + qJ(2) * t459) * qJD(1), t443, -qJD(4) * t380 + t467 * t375 + t470 * t432 - t438 * t469, -t469, -t480 * r_i_i_C(1) + t479 * r_i_i_C(2) + (t473 * r_i_i_C(1) + t472 * r_i_i_C(2)) * qJD(5), 0; 0, 0, t389 * qJD(4) + t467 * t383 + t438 * t384 - t474 * t432, t384, -t440 * t383 + ((-t389 * t414 - t398 * t413) * r_i_i_C(1) + (t389 * t413 - t398 * t414) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
