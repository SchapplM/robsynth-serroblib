% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:43
% EndTime: 2019-02-26 22:50:44
% DurationCPUTime: 0.63s
% Computational Cost: add. (1265->124), mult. (2263->185), div. (0->0), fcn. (2262->14), ass. (0->82)
t426 = sin(qJ(3));
t430 = cos(qJ(3));
t423 = qJ(4) + qJ(5);
t418 = cos(t423);
t429 = cos(qJ(4));
t405 = t429 * pkin(4) + pkin(5) * t418;
t403 = pkin(3) + t405;
t419 = qJ(6) + t423;
t414 = sin(t419);
t415 = cos(t419);
t446 = t415 * r_i_i_C(1) - t414 * r_i_i_C(2);
t442 = t403 + t446;
t473 = r_i_i_C(3) + pkin(12) + pkin(11) + pkin(10);
t417 = sin(t423);
t425 = sin(qJ(4));
t472 = pkin(4) * qJD(4);
t421 = qJD(4) + qJD(5);
t475 = pkin(5) * t421;
t397 = -t417 * t475 - t425 * t472;
t416 = qJD(6) + t421;
t445 = t414 * r_i_i_C(1) + t415 * r_i_i_C(2);
t479 = t445 * t416 - t397;
t432 = t430 * t479 + (t442 * t426 - t473 * t430) * qJD(3);
t428 = sin(qJ(1));
t431 = cos(qJ(2));
t471 = cos(pkin(6));
t476 = cos(qJ(1));
t447 = t471 * t476;
t427 = sin(qJ(2));
t456 = t428 * t471;
t448 = t427 * t456;
t457 = t476 * qJD(1);
t463 = qJD(2) * t427;
t387 = -qJD(1) * t448 - t428 * t463 + (qJD(2) * t447 + t457) * t431;
t400 = t427 * t447 + t428 * t431;
t424 = sin(pkin(6));
t462 = t424 * t476;
t450 = t430 * t462;
t470 = t424 * t428;
t460 = qJD(1) * t470;
t381 = (-qJD(3) * t400 + t460) * t426 - qJD(3) * t450 + t387 * t430;
t477 = t473 * t426 + t442 * t430 + pkin(2);
t404 = t425 * pkin(4) + pkin(5) * t417;
t474 = -pkin(9) - t404;
t469 = t424 * t430;
t468 = t424 * t431;
t467 = t428 * t427;
t443 = t431 * t447;
t461 = t476 * t431;
t384 = -qJD(1) * t443 - qJD(2) * t461 + (qJD(2) * t471 + qJD(1)) * t467;
t402 = t461 - t448;
t394 = t402 * t430 + t426 * t470;
t453 = -t394 * t416 - t384;
t401 = t476 * t427 + t431 * t456;
t385 = t400 * qJD(1) + t401 * qJD(2);
t441 = -t402 * t426 + t428 * t469;
t449 = t424 * t457;
t379 = t441 * qJD(3) - t385 * t430 + t426 * t449;
t455 = t401 * t416 + t379;
t374 = -t455 * t414 + t453 * t415;
t375 = t453 * t414 + t455 * t415;
t466 = t374 * r_i_i_C(1) - t375 * r_i_i_C(2);
t386 = t401 * qJD(1) + t400 * qJD(2);
t391 = t400 * t430 - t426 * t462;
t452 = t391 * t416 - t386;
t399 = -t443 + t467;
t454 = -t399 * t416 - t381;
t465 = (t454 * t414 - t452 * t415) * r_i_i_C(1) + (t452 * t414 + t454 * t415) * r_i_i_C(2);
t396 = t471 * t426 + t427 * t469;
t459 = t424 * t463;
t440 = -t396 * t416 + t459;
t437 = -t424 * t427 * t426 + t471 * t430;
t458 = qJD(2) * t468;
t389 = t437 * qJD(3) + t430 * t458;
t444 = t416 * t468 - t389;
t464 = (t444 * t414 + t440 * t415) * r_i_i_C(1) + (-t440 * t414 + t444 * t415) * r_i_i_C(2);
t439 = t445 - t474;
t398 = t418 * t475 + t429 * t472;
t435 = t446 * t416 + t398;
t433 = -t391 * qJD(3) - t387 * t426 + t430 * t460;
t378 = t394 * qJD(3) - t385 * t426 - t430 * t449;
t1 = [-t387 * pkin(2) - t391 * t397 - t399 * t398 - t442 * t381 + ((t391 * t414 - t399 * t415) * r_i_i_C(1) + (t391 * t415 + t399 * t414) * r_i_i_C(2)) * t416 - t439 * t386 + t473 * t433 + (-t476 * pkin(1) - pkin(8) * t470) * qJD(1), t384 * t477 - t439 * t385 + t432 * t401 + t435 * t402, -t442 * t378 + t473 * t379 - t441 * t479, -t379 * t404 - t384 * t405 - t394 * t398 + t401 * t397 + t466 ((-t394 * t421 - t384) * t418 + (-t401 * t421 - t379) * t417) * pkin(5) + t466, t466; -t385 * pkin(2) + t375 * r_i_i_C(1) + t374 * r_i_i_C(2) + t379 * t403 + t394 * t397 + t401 * t398 + t474 * t384 + t473 * t378 + (-pkin(1) * t428 + pkin(8) * t462) * qJD(1), -t386 * t477 + t439 * t387 + t432 * t399 + t435 * t400, t473 * t381 - t479 * (-t400 * t426 - t450) + t442 * t433, -t381 * t404 + t386 * t405 - t391 * t398 + t399 * t397 + t465 ((-t391 * t421 + t386) * t418 + (-t399 * t421 - t381) * t417) * pkin(5) + t465, t465; 0 ((-qJD(2) * t477 + t435) * t427 + (t439 * qJD(2) - t432) * t431) * t424, t473 * t389 - t479 * t437 + t442 * (-t396 * qJD(3) - t426 * t458) -t389 * t404 - t396 * t398 + (-t431 * t397 + t405 * t463) * t424 + t464 ((-t396 * t421 + t459) * t418 + (t421 * t468 - t389) * t417) * pkin(5) + t464, t464;];
JaD_transl  = t1;
