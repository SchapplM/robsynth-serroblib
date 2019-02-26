% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:23
% EndTime: 2019-02-26 22:36:24
% DurationCPUTime: 0.64s
% Computational Cost: add. (1021->114), mult. (2000->171), div. (0->0), fcn. (2003->14), ass. (0->76)
t426 = sin(qJ(3));
t430 = cos(qJ(3));
t423 = qJ(4) + pkin(12);
t403 = pkin(5) * cos(t423) + cos(qJ(4)) * pkin(4);
t401 = pkin(3) + t403;
t419 = qJ(6) + t423;
t415 = sin(t419);
t416 = cos(t419);
t445 = t416 * r_i_i_C(1) - t415 * r_i_i_C(2);
t441 = t401 + t445;
t472 = r_i_i_C(3) + pkin(11) + qJ(5) + pkin(10);
t402 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t423);
t399 = t402 * qJD(4);
t422 = qJD(4) + qJD(6);
t444 = t415 * r_i_i_C(1) + t416 * r_i_i_C(2);
t477 = t444 * t422 + t399;
t432 = t477 * t430 + (t441 * t426 - t472 * t430) * qJD(3) - t426 * qJD(5);
t428 = sin(qJ(1));
t431 = cos(qJ(2));
t471 = cos(pkin(6));
t474 = cos(qJ(1));
t446 = t471 * t474;
t427 = sin(qJ(2));
t455 = t428 * t471;
t447 = t427 * t455;
t456 = t474 * qJD(1);
t463 = qJD(2) * t427;
t385 = -qJD(1) * t447 - t428 * t463 + (qJD(2) * t446 + t456) * t431;
t424 = sin(pkin(6));
t460 = t424 * t474;
t480 = -qJD(3) * t460 + t385;
t396 = t427 * t446 + t428 * t431;
t470 = t424 * t428;
t479 = qJD(1) * t470 - qJD(3) * t396;
t379 = t479 * t426 + t480 * t430;
t475 = t472 * t426 + t441 * t430 + pkin(2);
t473 = -pkin(9) - t402;
t469 = t424 * t430;
t468 = t424 * t431;
t467 = t428 * t427;
t442 = t431 * t446;
t459 = t474 * t431;
t382 = -qJD(1) * t442 - qJD(2) * t459 + (qJD(2) * t471 + qJD(1)) * t467;
t398 = t459 - t447;
t392 = t398 * t430 + t426 * t470;
t452 = -t392 * t422 - t382;
t397 = t474 * t427 + t431 * t455;
t383 = t396 * qJD(1) + t397 * qJD(2);
t391 = -t398 * t426 + t428 * t469;
t448 = t424 * t456;
t377 = t391 * qJD(3) - t383 * t430 + t426 * t448;
t454 = t397 * t422 + t377;
t372 = -t454 * t415 + t452 * t416;
t373 = t452 * t415 + t454 * t416;
t466 = t372 * r_i_i_C(1) - t373 * r_i_i_C(2);
t384 = t397 * qJD(1) + t396 * qJD(2);
t389 = t396 * t430 - t426 * t460;
t451 = t389 * t422 - t384;
t395 = -t442 + t467;
t453 = -t395 * t422 - t379;
t465 = (t453 * t415 - t451 * t416) * r_i_i_C(1) + (t451 * t415 + t453 * t416) * r_i_i_C(2);
t394 = t471 * t426 + t427 * t469;
t440 = -t394 * t422 + t424 * t463;
t436 = -t424 * t427 * t426 + t471 * t430;
t457 = qJD(2) * t468;
t387 = t436 * qJD(3) + t430 * t457;
t443 = t422 * t468 - t387;
t464 = (t443 * t415 + t440 * t416) * r_i_i_C(1) + (-t440 * t415 + t443 * t416) * r_i_i_C(2);
t439 = t444 - t473;
t437 = t396 * t426 + t430 * t460;
t400 = t403 * qJD(4);
t434 = t445 * t422 + t400;
t378 = t480 * t426 - t479 * t430;
t386 = t394 * qJD(3) + t426 * t457;
t376 = t392 * qJD(3) - t383 * t426 - t430 * t448;
t1 = [t389 * t399 - t437 * qJD(5) - t395 * t400 - t385 * pkin(2) - t441 * t379 + ((t389 * t415 - t395 * t416) * r_i_i_C(1) + (t389 * t416 + t395 * t415) * r_i_i_C(2)) * t422 - t439 * t384 - t472 * t378 + (-t474 * pkin(1) - pkin(8) * t470) * qJD(1), t475 * t382 - t439 * t383 + t432 * t397 + t434 * t398, t392 * qJD(5) - t441 * t376 + t472 * t377 - t391 * t477, -t377 * t402 - t382 * t403 - t392 * t400 - t397 * t399 + t466, t376, t466; -t383 * pkin(2) + t373 * r_i_i_C(1) + t372 * r_i_i_C(2) - t391 * qJD(5) + t377 * t401 - t392 * t399 + t397 * t400 + t473 * t382 + t472 * t376 + (-pkin(1) * t428 + pkin(8) * t460) * qJD(1), -t384 * t475 + t439 * t385 + t432 * t395 + t434 * t396, t389 * qJD(5) - t441 * t378 + t472 * t379 + t437 * t477, -t379 * t402 + t384 * t403 - t389 * t400 - t395 * t399 + t465, t378, t465; 0 ((-qJD(2) * t475 + t434) * t427 + (t439 * qJD(2) - t432) * t431) * t424, t394 * qJD(5) - t441 * t386 + t472 * t387 - t436 * t477, -t387 * t402 - t394 * t400 + (t431 * t399 + t403 * t463) * t424 + t464, t386, t464;];
JaD_transl  = t1;
