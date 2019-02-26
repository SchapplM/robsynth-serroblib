% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRPR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:21
% EndTime: 2019-02-26 19:40:21
% DurationCPUTime: 0.50s
% Computational Cost: add. (761->83), mult. (2148->148), div. (0->0), fcn. (2498->16), ass. (0->72)
t478 = sin(pkin(12));
t479 = sin(pkin(11));
t465 = t479 * t478;
t482 = cos(pkin(12));
t483 = cos(pkin(11));
t471 = t483 * t482;
t485 = cos(pkin(6));
t447 = -t485 * t471 + t465;
t484 = cos(pkin(7));
t445 = t447 * t484;
t480 = sin(pkin(7));
t481 = sin(pkin(6));
t469 = t481 * t480;
t456 = t483 * t469;
t491 = t445 + t456;
t467 = t479 * t482;
t470 = t483 * t478;
t448 = t485 * t467 + t470;
t466 = t479 * t481;
t490 = t448 * t484 - t480 * t466;
t472 = t484 * t481;
t489 = t482 * t472 + t485 * t480;
t424 = t485 * t470 + t467;
t440 = sin(qJ(3));
t487 = cos(qJ(3));
t407 = t424 * t440 + t491 * t487;
t488 = t490 * t487;
t436 = pkin(13) + qJ(6);
t434 = sin(t436);
t435 = cos(t436);
t474 = t434 * r_i_i_C(1) + t435 * r_i_i_C(2);
t458 = qJD(6) * t474;
t468 = t481 * t478;
t415 = t440 * t468 - t489 * t487;
t486 = r_i_i_C(3) + pkin(10) + qJ(5);
t477 = qJD(3) * t440;
t476 = t424 * t487;
t475 = t435 * r_i_i_C(1) - t434 * r_i_i_C(2);
t408 = -t491 * t440 + t476;
t417 = t447 * t480 - t483 * t472;
t439 = sin(qJ(4));
t441 = cos(qJ(4));
t400 = t408 * t441 + t417 * t439;
t464 = -t408 * t439 + t417 * t441;
t425 = -t485 * t465 + t471;
t410 = t425 * t487 - t490 * t440;
t418 = t448 * t480 + t484 * t466;
t402 = t410 * t441 + t418 * t439;
t463 = -t410 * t439 + t418 * t441;
t416 = t489 * t440 + t487 * t468;
t423 = -t482 * t469 + t485 * t484;
t412 = t416 * t441 + t423 * t439;
t462 = -t416 * t439 + t423 * t441;
t461 = cos(pkin(13)) * pkin(5) + pkin(4) + t475;
t459 = qJD(6) * t475;
t453 = -sin(pkin(13)) * pkin(5) - pkin(9) - t474;
t450 = -t486 * t439 - t461 * t441 - pkin(3);
t442 = -t439 * qJD(5) + t441 * t458 + (t461 * t439 - t486 * t441) * qJD(4);
t414 = t416 * qJD(3);
t413 = t415 * qJD(3);
t409 = t425 * t440 + t488;
t406 = t410 * qJD(3);
t405 = t488 * qJD(3) + t425 * t477;
t404 = -t456 * t477 + (-t440 * t445 + t476) * qJD(3);
t403 = t407 * qJD(3);
t398 = t462 * qJD(4) - t413 * t441;
t397 = t412 * qJD(4) - t413 * t439;
t396 = t463 * qJD(4) - t405 * t441;
t395 = t402 * qJD(4) - t405 * t439;
t394 = t464 * qJD(4) - t403 * t441;
t393 = t400 * qJD(4) - t403 * t439;
t1 = [0, 0, t453 * t405 + t450 * t406 + t442 * t409 + t410 * t459, qJD(5) * t402 - t461 * t395 + t486 * t396 - t463 * t458, t395 (-t396 * t434 + t406 * t435) * r_i_i_C(1) + (-t396 * t435 - t406 * t434) * r_i_i_C(2) + ((-t402 * t435 - t409 * t434) * r_i_i_C(1) + (t402 * t434 - t409 * t435) * r_i_i_C(2)) * qJD(6); 0, 0, t453 * t403 + t450 * t404 + t442 * t407 + t408 * t459, qJD(5) * t400 - t461 * t393 + t486 * t394 - t464 * t458, t393 (-t394 * t434 + t404 * t435) * r_i_i_C(1) + (-t394 * t435 - t404 * t434) * r_i_i_C(2) + ((-t400 * t435 - t407 * t434) * r_i_i_C(1) + (t400 * t434 - t407 * t435) * r_i_i_C(2)) * qJD(6); 0, 0, t453 * t413 + t450 * t414 + t442 * t415 + t416 * t459, qJD(5) * t412 - t461 * t397 + t486 * t398 - t462 * t458, t397 (-t398 * t434 + t414 * t435) * r_i_i_C(1) + (-t398 * t435 - t414 * t434) * r_i_i_C(2) + ((-t412 * t435 - t415 * t434) * r_i_i_C(1) + (t412 * t434 - t415 * t435) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
