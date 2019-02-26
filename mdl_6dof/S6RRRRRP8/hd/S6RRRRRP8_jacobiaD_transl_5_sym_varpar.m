% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:50
% EndTime: 2019-02-26 22:43:51
% DurationCPUTime: 0.72s
% Computational Cost: add. (953->124), mult. (1772->205), div. (0->0), fcn. (1760->12), ass. (0->81)
t491 = pkin(11) + r_i_i_C(3);
t435 = sin(qJ(5));
t439 = cos(qJ(5));
t453 = r_i_i_C(1) * t439 - r_i_i_C(2) * t435;
t451 = pkin(4) + t453;
t438 = sin(qJ(1));
t434 = cos(pkin(6));
t455 = qJD(2) * t434 + qJD(1);
t437 = sin(qJ(2));
t476 = t437 * t438;
t462 = t434 * t476;
t471 = qJD(2) * t437;
t441 = cos(qJ(2));
t442 = cos(qJ(1));
t473 = t441 * t442;
t406 = -qJD(1) * t462 - t438 * t471 + t455 * t473;
t431 = qJD(3) + qJD(4);
t433 = sin(pkin(6));
t477 = t433 * t442;
t500 = t431 * t477 - t406;
t467 = qJD(5) * t439;
t468 = qJD(5) * t435;
t499 = -r_i_i_C(1) * t468 - t467 * r_i_i_C(2);
t474 = t438 * t441;
t475 = t437 * t442;
t417 = t434 * t475 + t474;
t432 = qJ(3) + qJ(4);
t429 = sin(t432);
t430 = cos(t432);
t409 = -t417 * t430 + t429 * t477;
t461 = t434 * t473;
t416 = -t461 + t476;
t498 = -t409 * t435 - t416 * t439;
t497 = t409 * t439 - t416 * t435;
t436 = sin(qJ(3));
t496 = -qJD(3) * t436 * pkin(3) - (t451 * t429 - t491 * t430) * t431;
t494 = r_i_i_C(1) * t435 + r_i_i_C(2) * t439;
t440 = cos(qJ(3));
t428 = pkin(3) * t440 + pkin(2);
t492 = t491 * t429 + t451 * t430 + t428;
t418 = t434 * t474 + t475;
t405 = t418 * qJD(1) + t417 * qJD(2);
t487 = t405 * t435;
t486 = t405 * t439;
t419 = -t462 + t473;
t483 = t419 * t430;
t482 = t429 * t431;
t481 = t433 * t437;
t480 = t433 * t438;
t479 = t433 * t440;
t478 = t433 * t441;
t472 = qJD(1) * t433;
t470 = qJD(2) * t441;
t469 = qJD(5) * t430;
t464 = t431 * t481;
t460 = t438 * t472;
t459 = t442 * t472;
t458 = t433 * t471;
t457 = t433 * t470;
t456 = t500 * t430;
t404 = t417 * qJD(1) + t418 * qJD(2);
t452 = t431 * t480 - t404;
t450 = -t417 * t431 + t460;
t449 = t431 * t434 + t457;
t394 = t500 * t429 + t450 * t430;
t395 = -t417 * t482 + t429 * t460 - t456;
t447 = t499 * (-t417 * t429 - t430 * t477) + t491 * t395 + t451 * t394;
t392 = t452 * t429 - t430 * t459 + t431 * t483;
t393 = -t419 * t482 + t429 * t459 + t452 * t430;
t446 = t499 * (-t419 * t429 + t430 * t480) + t491 * t393 - t451 * t392;
t402 = -t429 * t464 + t449 * t430;
t445 = t499 * (-t429 * t481 + t430 * t434) + t491 * t402 + t451 * (-t449 * t429 - t430 * t464);
t444 = t494 * t469 - t496;
t443 = -pkin(10) - pkin(9);
t415 = t429 * t434 + t430 * t481;
t411 = t429 * t480 + t483;
t403 = -qJD(1) * t461 - t442 * t470 + t455 * t476;
t397 = -t450 * t429 + t456;
t385 = t393 * t439 - t403 * t435 + (-t411 * t435 + t418 * t439) * qJD(5);
t384 = -t393 * t435 - t403 * t439 + (-t411 * t439 - t418 * t435) * qJD(5);
t1 = [(t397 * t439 - t487) * r_i_i_C(1) + (-t397 * t435 - t486) * r_i_i_C(2) + t397 * pkin(4) - t406 * t428 + t405 * t443 + t491 * t394 + (t498 * r_i_i_C(1) - t497 * r_i_i_C(2)) * qJD(5) + (-pkin(1) * t442 - pkin(8) * t480) * qJD(1) + (-t436 * t460 + (t417 * t436 + t440 * t477) * qJD(3)) * pkin(3) (-t404 * t435 + t419 * t467) * r_i_i_C(1) + (-t404 * t439 - t419 * t468) * r_i_i_C(2) + t404 * t443 + t492 * t403 + t444 * t418 (t440 * t459 + t404 * t436 + (-t419 * t440 - t436 * t480) * qJD(3)) * pkin(3) + t446, t446, r_i_i_C(1) * t384 - t385 * r_i_i_C(2), 0; t393 * pkin(4) + t385 * r_i_i_C(1) + t384 * r_i_i_C(2) + t403 * t443 - t404 * t428 + t491 * t392 + (-pkin(1) * t438 + pkin(8) * t477) * qJD(1) + (t436 * t459 + (-t419 * t436 + t438 * t479) * qJD(3)) * pkin(3) (t406 * t435 + t417 * t467) * r_i_i_C(1) + (t406 * t439 - t417 * t468) * r_i_i_C(2) - t406 * t443 - t492 * t405 + t444 * t416 (t440 * t460 - t406 * t436 + (-t417 * t440 + t436 * t477) * qJD(3)) * pkin(3) + t447, t447 (-t395 * t435 + t486) * r_i_i_C(1) + (-t395 * t439 - t487) * r_i_i_C(2) + (t497 * r_i_i_C(1) + t498 * r_i_i_C(2)) * qJD(5), 0; 0 ((-qJD(2) * t492 + t453 * qJD(5)) * t437 + (-qJD(2) * t443 + t494 * (qJD(2) - t469) + t496) * t441) * t433 (-t436 * t457 + (-t434 * t436 - t437 * t479) * qJD(3)) * pkin(3) + t445, t445 (-t402 * t435 + t439 * t458) * r_i_i_C(1) + (-t402 * t439 - t435 * t458) * r_i_i_C(2) + ((-t415 * t439 + t435 * t478) * r_i_i_C(1) + (t415 * t435 + t439 * t478) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
