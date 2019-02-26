% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:36
% EndTime: 2019-02-26 22:29:37
% DurationCPUTime: 0.81s
% Computational Cost: add. (1012->129), mult. (2977->203), div. (0->0), fcn. (3075->10), ass. (0->78)
t416 = sin(qJ(2));
t417 = sin(qJ(1));
t420 = cos(qJ(2));
t421 = cos(qJ(1));
t471 = cos(pkin(6));
t440 = t421 * t471;
t399 = t416 * t440 + t417 * t420;
t441 = t417 * t471;
t400 = t421 * t416 + t420 * t441;
t385 = t400 * qJD(1) + t399 * qJD(2);
t418 = cos(qJ(4));
t413 = sin(pkin(6));
t419 = cos(qJ(3));
t415 = sin(qJ(3));
t461 = t415 * t421;
t391 = -t399 * t419 + t413 * t461;
t398 = t416 * t417 - t420 * t440;
t414 = sin(qJ(4));
t481 = t391 * t418 - t398 * t414;
t484 = -t481 * qJD(4) - t385 * t418;
t482 = t391 * t414 + t398 * t418;
t483 = t482 * qJD(4) + t385 * t414;
t451 = -r_i_i_C(3) - qJ(6) + pkin(10);
t454 = qJD(5) * t414;
t480 = (-pkin(3) * t415 + t451 * t419) * qJD(3) - t415 * qJD(6) + t419 * t454;
t452 = r_i_i_C(1) + pkin(5) + pkin(4);
t472 = r_i_i_C(2) + qJ(5);
t479 = -t454 + (t452 * t414 - t472 * t418) * qJD(4);
t464 = t413 * t419;
t397 = t471 * t415 + t416 * t464;
t463 = t413 * t420;
t478 = -t397 * t418 + t414 * t463;
t477 = t419 * pkin(3) + t451 * t415 + pkin(2);
t456 = qJD(3) * t415;
t476 = (qJD(2) * t419 - qJD(4)) * t416 + t420 * t456;
t426 = t472 * t414 + t452 * t418 + pkin(3);
t465 = t413 * t417;
t462 = t413 * t421;
t460 = t417 * t419;
t459 = t421 * t420;
t458 = qJD(1) * t417;
t457 = qJD(2) * t416;
t455 = qJD(4) * t419;
t453 = qJD(5) * t418;
t449 = t413 * t416 * t415;
t448 = t419 * t462;
t446 = t413 * t458;
t445 = qJD(1) * t464;
t444 = t413 * t457;
t443 = qJD(2) * t463;
t437 = t416 * t441;
t386 = -qJD(1) * t437 - t417 * t457 + (qJD(2) * t471 + qJD(1)) * t459;
t439 = qJD(3) * t448 - t386 * t419;
t384 = t399 * qJD(1) + t400 * qJD(2);
t436 = t400 * t455 - t384;
t435 = t398 * t455 + t386;
t401 = -t437 + t459;
t393 = t401 * t419 + t415 * t465;
t432 = t393 * t418 + t400 * t414;
t431 = -t393 * t414 + t400 * t418;
t430 = (qJD(2) - t455) * t420;
t429 = -t399 * t415 - t448;
t383 = qJD(2) * t437 + t416 * t458 + (-t471 * qJD(1) - qJD(2)) * t459;
t425 = qJD(4) * t401 + t383 * t419 + t400 * t456;
t424 = qJD(4) * t399 - t385 * t419 + t398 * t456;
t373 = t391 * qJD(3) - t386 * t415 + t417 * t445;
t392 = t401 * t415 - t413 * t460;
t388 = -qJD(3) * t449 + (t471 * qJD(3) + t443) * t419;
t387 = -t397 * qJD(3) - t415 * t443;
t377 = -t478 * qJD(4) + t388 * t414 - t418 * t444;
t376 = (qJD(3) * t399 - t446) * t415 + t439;
t374 = -t399 * t456 + t415 * t446 - t439;
t372 = -t384 * t419 - t401 * t456 + (qJD(1) * t461 + qJD(3) * t460) * t413;
t371 = t393 * qJD(3) - t384 * t415 - t421 * t445;
t363 = t374 * t414 + t484;
t362 = t431 * qJD(4) + t372 * t418 - t383 * t414;
t361 = t432 * qJD(4) + t372 * t414 + t383 * t418;
t1 = [-t429 * qJD(6) + t482 * qJD(5) + t376 * pkin(3) - t386 * pkin(2) - t385 * pkin(9) + t472 * (t376 * t414 - t484) + (-pkin(1) * t421 - pkin(8) * t465) * qJD(1) + t451 * t373 + t452 * (t376 * t418 - t483) -t401 * t453 - t384 * pkin(9) + t472 * (t425 * t414 - t436 * t418) + t452 * (t436 * t414 + t425 * t418) - t480 * t400 + t477 * t383, -t393 * qJD(6) - t426 * t371 + t451 * t372 + t479 * t392, t432 * qJD(5) - t452 * t361 + t472 * t362, t361, -t371; -t392 * qJD(6) - t431 * qJD(5) + t372 * pkin(3) - t384 * pkin(2) - t383 * pkin(9) + t472 * t361 + (-pkin(1) * t417 + pkin(8) * t462) * qJD(1) + t451 * t371 + t452 * t362, -t399 * t453 + t386 * pkin(9) + t472 * (t424 * t414 - t435 * t418) + t452 * (t435 * t414 + t424 * t418) - t480 * t398 - t477 * t385, qJD(6) * t391 + t426 * t373 + t451 * t374 - t429 * t479, -t481 * qJD(5) + t472 * (t374 * t418 + t483) - t452 * t363, t363, t373; 0 (-t472 * (t476 * t414 + t418 * t430) + t452 * (t414 * t430 - t476 * t418) - t416 * t453 + t480 * t420 + (t420 * pkin(9) - t416 * t477) * qJD(2)) * t413, -t397 * qJD(6) + t451 * t388 - t479 * (t471 * t419 - t449) + t426 * t387, -t478 * qJD(5) + t472 * (t414 * t444 + t388 * t418 + (-t397 * t414 - t418 * t463) * qJD(4)) - t452 * t377, t377, t387;];
JaD_transl  = t1;
