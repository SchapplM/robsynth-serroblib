% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP9
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
% Datum: 2019-02-26 22:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:44:27
% EndTime: 2019-02-26 22:44:28
% DurationCPUTime: 0.58s
% Computational Cost: add. (993->118), mult. (2126->178), div. (0->0), fcn. (2126->12), ass. (0->81)
t424 = sin(qJ(3));
t428 = cos(qJ(3));
t421 = qJ(4) + qJ(5);
t417 = cos(t421);
t427 = cos(qJ(4));
t404 = t427 * pkin(4) + pkin(5) * t417;
t402 = pkin(3) + t404;
t416 = sin(t421);
t445 = t417 * r_i_i_C(1) - t416 * r_i_i_C(2);
t441 = t402 + t445;
t473 = r_i_i_C(3) + qJ(6) + pkin(11) + pkin(10);
t423 = sin(qJ(4));
t472 = pkin(4) * qJD(4);
t420 = qJD(4) + qJD(5);
t475 = pkin(5) * t420;
t396 = -t416 * t475 - t423 * t472;
t444 = t416 * r_i_i_C(1) + t417 * r_i_i_C(2);
t479 = t444 * t420 - t396;
t430 = t479 * t428 + (t441 * t424 - t473 * t428) * qJD(3) - t424 * qJD(6);
t426 = sin(qJ(1));
t429 = cos(qJ(2));
t471 = cos(pkin(6));
t476 = cos(qJ(1));
t446 = t471 * t476;
t425 = sin(qJ(2));
t455 = t426 * t471;
t447 = t425 * t455;
t456 = t476 * qJD(1);
t463 = qJD(2) * t425;
t386 = -qJD(1) * t447 - t426 * t463 + (qJD(2) * t446 + t456) * t429;
t422 = sin(pkin(6));
t460 = t422 * t476;
t482 = -qJD(3) * t460 + t386;
t399 = t425 * t446 + t426 * t429;
t470 = t422 * t426;
t481 = qJD(1) * t470 - qJD(3) * t399;
t380 = t481 * t424 + t482 * t428;
t477 = t473 * t424 + t441 * t428 + pkin(2);
t403 = t423 * pkin(4) + pkin(5) * t416;
t474 = -pkin(9) - t403;
t469 = t422 * t428;
t468 = t422 * t429;
t467 = t426 * t425;
t442 = t429 * t446;
t459 = t476 * t429;
t383 = -qJD(1) * t442 - qJD(2) * t459 + (qJD(2) * t471 + qJD(1)) * t467;
t401 = t459 - t447;
t393 = t401 * t428 + t424 * t470;
t452 = -t393 * t420 - t383;
t400 = t476 * t425 + t429 * t455;
t384 = t399 * qJD(1) + t400 * qJD(2);
t392 = -t401 * t424 + t426 * t469;
t448 = t422 * t456;
t378 = t392 * qJD(3) - t384 * t428 + t424 * t448;
t454 = t400 * t420 + t378;
t373 = -t454 * t416 + t452 * t417;
t374 = t452 * t416 + t454 * t417;
t466 = t373 * r_i_i_C(1) - t374 * r_i_i_C(2);
t385 = t400 * qJD(1) + t399 * qJD(2);
t390 = t399 * t428 - t424 * t460;
t451 = t390 * t420 - t385;
t398 = -t442 + t467;
t453 = -t398 * t420 - t380;
t433 = t453 * t416 - t451 * t417;
t465 = t433 * r_i_i_C(1) + (t451 * t416 + t453 * t417) * r_i_i_C(2);
t395 = t471 * t424 + t425 * t469;
t440 = -t395 * t420 + t422 * t463;
t436 = -t422 * t425 * t424 + t471 * t428;
t457 = qJD(2) * t468;
t388 = t436 * qJD(3) + t428 * t457;
t443 = t420 * t468 - t388;
t431 = t443 * t416 + t440 * t417;
t464 = t431 * r_i_i_C(1) + (-t440 * t416 + t443 * t417) * r_i_i_C(2);
t439 = t444 - t474;
t437 = t399 * t424 + t428 * t460;
t397 = t417 * t475 + t427 * t472;
t434 = t445 * t420 + t397;
t379 = t482 * t424 - t481 * t428;
t387 = t395 * qJD(3) + t424 * t457;
t377 = t393 * qJD(3) - t384 * t424 - t428 * t448;
t1 = [-t390 * t396 - t437 * qJD(6) - t398 * t397 - t386 * pkin(2) - t441 * t380 + ((t390 * t416 - t398 * t417) * r_i_i_C(1) + (t390 * t417 + t398 * t416) * r_i_i_C(2)) * t420 - t439 * t385 - t473 * t379 + (-t476 * pkin(1) - pkin(8) * t470) * qJD(1), t477 * t383 - t439 * t384 + t430 * t400 + t434 * t401, t393 * qJD(6) - t441 * t377 + t473 * t378 - t392 * t479, -t378 * t403 - t383 * t404 - t393 * t397 + t400 * t396 + t466, pkin(5) * t373 + t466, t377; -t384 * pkin(2) + t374 * r_i_i_C(1) + t373 * r_i_i_C(2) - t392 * qJD(6) + t378 * t402 + t393 * t396 + t400 * t397 + t474 * t383 + t473 * t377 + (-pkin(1) * t426 + pkin(8) * t460) * qJD(1), -t385 * t477 + t386 * t439 + t398 * t430 + t399 * t434, t390 * qJD(6) - t441 * t379 + t473 * t380 + t437 * t479, -t380 * t403 + t385 * t404 - t390 * t397 + t398 * t396 + t465, pkin(5) * t433 + t465, t379; 0 ((-qJD(2) * t477 + t434) * t425 + (t439 * qJD(2) - t430) * t429) * t422, t395 * qJD(6) - t441 * t387 + t473 * t388 - t436 * t479, -t388 * t403 - t395 * t397 + (-t429 * t396 + t404 * t463) * t422 + t464, pkin(5) * t431 + t464, t387;];
JaD_transl  = t1;
