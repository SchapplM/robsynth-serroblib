% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR9_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR9_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:30
% EndTime: 2019-02-26 21:05:31
% DurationCPUTime: 0.62s
% Computational Cost: add. (581->105), mult. (1716->177), div. (0->0), fcn. (1834->14), ass. (0->75)
t422 = cos(pkin(12));
t424 = cos(pkin(6));
t419 = sin(pkin(12));
t428 = sin(qJ(1));
t460 = t428 * t419;
t450 = t424 * t460;
t431 = cos(qJ(1));
t455 = qJD(1) * t431;
t400 = -qJD(1) * t450 + t422 * t455;
t420 = sin(pkin(7));
t421 = sin(pkin(6));
t463 = t421 * t431;
t449 = t420 * t463;
t486 = -qJD(3) * t449 + t400;
t423 = cos(pkin(7));
t457 = t431 * t422;
t439 = t424 * t457 - t460;
t480 = t439 * t423;
t485 = -t480 + t449;
t458 = t431 * t419;
t459 = t428 * t422;
t438 = t424 * t459 + t458;
t399 = t438 * qJD(1);
t403 = t424 * t458 + t459;
t427 = sin(qJ(3));
t430 = cos(qJ(3));
t464 = t421 * t428;
t447 = qJD(1) * t464;
t378 = (-qJD(3) * t403 - t399 * t423 + t420 * t447) * t427 + (qJD(3) * t480 + t486) * t430;
t470 = t399 * t420;
t390 = t423 * t447 + t470;
t418 = qJ(4) + pkin(13);
t416 = sin(t418);
t417 = cos(t418);
t484 = t378 * t416 - t390 * t417;
t483 = -t378 * t417 - t390 * t416;
t479 = t439 * t420 + t423 * t463;
t461 = t423 * t430;
t466 = t420 * t424;
t478 = (-t419 * t427 + t422 * t461) * t421 + t430 * t466;
t462 = t423 * t427;
t442 = -t403 * t430 - t439 * t462;
t383 = t427 * t449 + t442;
t477 = t383 * t417 + t416 * t479;
t476 = -t383 * t416 + t417 * t479;
t474 = t403 * t427 + t485 * t430;
t456 = qJD(1) * t430;
t465 = t421 * t420;
t446 = t456 * t465;
t473 = t442 * qJD(3) - t399 * t461 - t486 * t427 + t428 * t446;
t426 = sin(qJ(4));
t472 = t426 * pkin(4);
t471 = r_i_i_C(3) + qJ(5) + pkin(10);
t453 = t421 * qJD(2);
t429 = cos(qJ(4));
t415 = t429 * pkin(4) + pkin(3);
t443 = -t417 * r_i_i_C(1) + t416 * r_i_i_C(2) - t415;
t440 = t420 * t464 - t423 * t438;
t437 = t416 * r_i_i_C(1) + t417 * r_i_i_C(2) + t472;
t434 = qJD(4) * t437;
t405 = -t450 + t457;
t384 = -t405 * t427 + t440 * t430;
t385 = t405 * t430 + t440 * t427;
t392 = t427 * t466 + (t419 * t430 + t422 * t462) * t421;
t401 = -t422 * t465 + t424 * t423;
t398 = t403 * qJD(1);
t395 = t420 * t438 + t423 * t464;
t388 = t479 * qJD(1);
t387 = t392 * qJD(3);
t386 = t478 * qJD(3);
t376 = t485 * t427 * qJD(1) + t384 * qJD(3) - t398 * t430;
t375 = t385 * qJD(3) - t398 * t427 - t431 * t446 + t456 * t480;
t374 = t376 * t417 + t388 * t416 + (-t385 * t416 + t395 * t417) * qJD(4);
t373 = -t376 * t416 + t388 * t417 + (-t385 * t417 - t395 * t416) * qJD(4);
t1 = [t483 * r_i_i_C(1) + t484 * r_i_i_C(2) - t378 * t415 - t474 * qJD(5) - t390 * t472 - t400 * pkin(2) - pkin(9) * t470 + t431 * t453 + t471 * t473 + (-t431 * pkin(1) + (-pkin(9) * t423 - qJ(2)) * t464) * qJD(1) + (t476 * r_i_i_C(1) - t477 * r_i_i_C(2) + (-t383 * t426 + t429 * t479) * pkin(4)) * qJD(4), t421 * t455, t385 * qJD(5) + t443 * t375 + t471 * t376 - t384 * t434, t373 * r_i_i_C(1) - t374 * r_i_i_C(2) + (-t376 * t426 + t388 * t429 + (-t385 * t429 - t395 * t426) * qJD(4)) * pkin(4), t375, 0; t428 * t453 - t398 * pkin(2) + t374 * r_i_i_C(1) + t373 * r_i_i_C(2) - t384 * qJD(5) + t376 * t415 + t471 * t375 + (t388 * t426 + (-t385 * t426 + t395 * t429) * qJD(4)) * pkin(4) + (-t428 * pkin(1) + pkin(9) * t479 + qJ(2) * t463) * qJD(1), t447, -qJD(5) * t383 + t471 * t378 + t474 * t434 - t443 * t473, -t484 * r_i_i_C(1) + t483 * r_i_i_C(2) + (t477 * r_i_i_C(1) + t476 * r_i_i_C(2)) * qJD(4) + (-t378 * t426 + t390 * t429 + (t383 * t429 + t426 * t479) * qJD(4)) * pkin(4), -t473, 0; 0, 0, t392 * qJD(5) + t471 * t386 + t443 * t387 - t478 * t434, -t437 * t386 + ((-t392 * t417 - t401 * t416) * r_i_i_C(1) + (t392 * t416 - t401 * t417) * r_i_i_C(2) + (-t392 * t429 - t401 * t426) * pkin(4)) * qJD(4), t387, 0;];
JaD_transl  = t1;
