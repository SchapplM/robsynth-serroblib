% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:10
% EndTime: 2019-02-26 21:40:11
% DurationCPUTime: 0.43s
% Computational Cost: add. (661->92), mult. (1960->157), div. (0->0), fcn. (2112->12), ass. (0->65)
t397 = sin(pkin(6));
t405 = cos(qJ(1));
t437 = t397 * t405;
t399 = cos(pkin(6));
t396 = sin(pkin(11));
t404 = cos(qJ(2));
t439 = cos(pkin(11));
t422 = qJD(2) * t439;
t401 = sin(qJ(2));
t430 = qJD(2) * t401;
t444 = t396 * t430 - t404 * t422;
t373 = t444 * t399;
t384 = -t404 * t396 - t401 * t439;
t379 = t384 * t399;
t429 = qJD(2) * t404;
t381 = -t396 * t429 - t401 * t422;
t402 = sin(qJ(1));
t410 = -t401 * t396 + t404 * t439;
t432 = qJD(1) * t402;
t442 = -t379 * t432 - t402 * t381 + (-qJD(1) * t410 + t373) * t405;
t446 = qJD(4) * t437 + t442;
t366 = -t379 * t405 + t402 * t410;
t426 = t397 * t432;
t445 = -qJD(4) * t366 + t426;
t400 = sin(qJ(4));
t403 = cos(qJ(4));
t443 = t445 * t400 - t446 * t403;
t395 = sin(pkin(12));
t398 = cos(pkin(12));
t416 = r_i_i_C(1) * t398 - r_i_i_C(2) * t395 + pkin(4);
t440 = r_i_i_C(3) + qJ(5);
t407 = t440 * t400 + t416 * t403 + pkin(3);
t441 = pkin(2) * t399;
t438 = t397 * t402;
t436 = t401 * t402;
t435 = t401 * t405;
t434 = t402 * t404;
t433 = t404 * t405;
t431 = qJD(1) * t405;
t427 = pkin(2) * t430;
t425 = t397 * t431;
t377 = t384 * t397;
t418 = -t377 * t403 + t399 * t400;
t378 = t410 * t399;
t417 = t378 * t405 + t402 * t384;
t415 = r_i_i_C(1) * t395 + r_i_i_C(2) * t398 + pkin(9);
t368 = t402 * t379 + t405 * t410;
t414 = -t368 * t400 + t403 * t438;
t413 = t368 * t403 + t400 * t438;
t409 = qJD(2) * t384;
t408 = qJD(2) * t410;
t353 = -t446 * t400 - t445 * t403;
t359 = -t366 * qJD(1) + t402 * t373 + t405 * t381;
t406 = t400 * qJD(5) + (-t416 * t400 + t440 * t403) * qJD(4);
t394 = pkin(2) * t404 + pkin(1);
t382 = -qJD(3) * t397 + t429 * t441;
t380 = t401 * t441 + (-pkin(8) - qJ(3)) * t397;
t374 = t399 * t409;
t371 = t444 * t397;
t363 = t418 * qJD(4) - t371 * t400;
t361 = t405 * t374 - t378 * t432 + t384 * t431 - t402 * t408;
t358 = t417 * qJD(1) + t402 * t374 + t405 * t408;
t352 = t414 * qJD(4) + t359 * t403 + t400 * t425;
t351 = t413 * qJD(4) + t359 * t400 - t403 * t425;
t1 = [(t361 * t395 - t398 * t443) * r_i_i_C(1) + (t361 * t398 + t395 * t443) * r_i_i_C(2) - t443 * pkin(4) - (t366 * t400 + t403 * t437) * qJD(5) + t442 * pkin(3) + t361 * pkin(9) + t402 * t427 - t405 * t382 - t440 * t353 + (t380 * t402 - t394 * t405) * qJD(1), t415 * t359 + t406 * (-t402 * t378 + t384 * t405) - t407 * t358 + ((t399 * t436 - t433) * qJD(2) + (-t399 * t433 + t436) * qJD(1)) * pkin(2), t425, t413 * qJD(5) - t416 * t351 + t440 * t352, t351, 0; (t352 * t398 + t358 * t395) * r_i_i_C(1) + (-t352 * t395 + t358 * t398) * r_i_i_C(2) + t352 * pkin(4) - t414 * qJD(5) + t359 * pkin(3) + t358 * pkin(9) - t405 * t427 - t402 * t382 + t440 * t351 + (-t380 * t405 - t394 * t402) * qJD(1), -t415 * t442 + t406 * t417 + t407 * t361 + ((-t399 * t435 - t434) * qJD(2) + (-t399 * t434 - t435) * qJD(1)) * pkin(2), t426 -(-t366 * t403 + t400 * t437) * qJD(5) + t440 * t443 - t416 * t353, t353, 0; 0, -t415 * t371 + (t406 * t410 + t407 * t409 - t427) * t397, 0, t418 * qJD(5) + t440 * (-t371 * t403 + (t377 * t400 + t399 * t403) * qJD(4)) - t416 * t363, t363, 0;];
JaD_transl  = t1;
