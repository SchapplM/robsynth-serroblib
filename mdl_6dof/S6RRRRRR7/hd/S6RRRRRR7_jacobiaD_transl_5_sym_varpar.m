% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_transl = S6RRRRRR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:53
% EndTime: 2019-02-26 22:50:54
% DurationCPUTime: 0.70s
% Computational Cost: add. (787->119), mult. (1843->197), div. (0->0), fcn. (1852->12), ass. (0->76)
t412 = sin(qJ(1));
t465 = cos(pkin(6));
t427 = qJD(2) * t465 + qJD(1);
t411 = sin(qJ(2));
t442 = t412 * t465;
t432 = t411 * t442;
t451 = qJD(2) * t411;
t415 = cos(qJ(2));
t416 = cos(qJ(1));
t456 = t416 * t415;
t382 = -qJD(1) * t432 - t412 * t451 + t427 * t456;
t441 = t416 * t465;
t393 = t411 * t441 + t412 * t415;
t410 = sin(qJ(3));
t414 = cos(qJ(3));
t408 = sin(pkin(6));
t452 = qJD(1) * t408;
t446 = t412 * t452;
t458 = t408 * t416;
t447 = t414 * t458;
t376 = (-qJD(3) * t393 + t446) * t410 - qJD(3) * t447 + t382 * t414;
t431 = t415 * t441;
t457 = t412 * t411;
t392 = -t431 + t457;
t406 = qJD(4) + qJD(5);
t438 = t392 * t406 + t376;
t413 = cos(qJ(4));
t403 = t413 * pkin(4) + pkin(3);
t407 = qJ(4) + qJ(5);
t404 = sin(t407);
t405 = cos(t407);
t430 = t405 * r_i_i_C(1) - t404 * r_i_i_C(2);
t426 = t403 + t430;
t466 = r_i_i_C(3) + pkin(11) + pkin(10);
t474 = (t426 * t410 - t466 * t414) * qJD(3);
t394 = t416 * t411 + t415 * t442;
t381 = t394 * qJD(1) + t393 * qJD(2);
t387 = -t393 * t414 + t410 * t458;
t473 = -t387 * t406 - t381;
t429 = t404 * r_i_i_C(1) + t405 * r_i_i_C(2);
t409 = sin(qJ(4));
t467 = t409 * pkin(4);
t472 = qJD(4) * t467 + t429 * t406;
t470 = t466 * t410 + t426 * t414 + pkin(2);
t463 = t404 * t406;
t462 = t405 * t406;
t461 = t408 * t412;
t460 = t408 * t414;
t459 = t408 * t415;
t450 = qJD(2) * t415;
t379 = -qJD(1) * t431 - t416 * t450 + t427 * t457;
t395 = -t432 + t456;
t389 = t395 * t414 + t410 * t461;
t437 = -t389 * t406 - t379;
t380 = t393 * qJD(1) + t394 * qJD(2);
t425 = -t395 * t410 + t412 * t460;
t445 = t416 * t452;
t374 = t425 * qJD(3) - t380 * t414 + t410 * t445;
t440 = t394 * t406 + t374;
t369 = -t440 * t404 + t437 * t405;
t370 = t437 * t404 + t440 * t405;
t455 = t369 * r_i_i_C(1) - t370 * r_i_i_C(2);
t454 = (-t404 * t438 - t405 * t473) * r_i_i_C(1) + (t404 * t473 - t405 * t438) * r_i_i_C(2);
t391 = t465 * t410 + t411 * t460;
t444 = t408 * t451;
t424 = -t391 * t406 + t444;
t422 = -t408 * t411 * t410 + t465 * t414;
t443 = t408 * t450;
t384 = t422 * qJD(3) + t414 * t443;
t428 = t406 * t459 - t384;
t453 = (t428 * t404 + t424 * t405) * r_i_i_C(1) + (-t424 * t404 + t428 * t405) * r_i_i_C(2);
t449 = qJD(4) * t413;
t419 = t387 * qJD(3) - t382 * t410 + t414 * t446;
t418 = t472 * t414 + t474;
t373 = t389 * qJD(3) - t380 * t410 - t414 * t445;
t1 = [-t382 * pkin(2) - t381 * pkin(9) - t376 * t403 + (-t438 * r_i_i_C(1) + r_i_i_C(2) * t473) * t405 + (r_i_i_C(1) * t473 + t438 * r_i_i_C(2)) * t404 + t466 * t419 + (-t416 * pkin(1) - pkin(8) * t461) * qJD(1) + (-t381 * t409 + (-t387 * t409 - t392 * t413) * qJD(4)) * pkin(4) (-t380 * t404 + t395 * t462) * r_i_i_C(1) + (-t380 * t405 - t395 * t463) * r_i_i_C(2) - t380 * pkin(9) + (-t380 * t409 + t395 * t449) * pkin(4) + t470 * t379 + t418 * t394, -t426 * t373 + t466 * t374 - t425 * t472 (-t374 * t409 - t379 * t413 + (-t389 * t413 - t394 * t409) * qJD(4)) * pkin(4) + t455, t455, 0; -t380 * pkin(2) - t379 * pkin(9) + t370 * r_i_i_C(1) + t369 * r_i_i_C(2) + t374 * t403 + t466 * t373 + (-pkin(1) * t412 + pkin(8) * t458) * qJD(1) + (-t379 * t409 + (-t389 * t409 + t394 * t413) * qJD(4)) * pkin(4) (t382 * t404 + t393 * t462) * r_i_i_C(1) + (t382 * t405 - t393 * t463) * r_i_i_C(2) + t382 * pkin(9) + (t382 * t409 + t393 * t449) * pkin(4) - t470 * t381 + t418 * t392, t466 * t376 - t472 * (-t393 * t410 - t447) + t426 * t419 (-t376 * t409 + t381 * t413 + (t387 * t413 - t392 * t409) * qJD(4)) * pkin(4) + t454, t454, 0; 0 ((pkin(4) * t449 - qJD(2) * t470 + t430 * t406) * t411 + (qJD(2) * pkin(9) + (-qJD(4) * t414 + qJD(2)) * t467 - t474 + t429 * (-t406 * t414 + qJD(2))) * t415) * t408, t466 * t384 - t472 * t422 + t426 * (-t391 * qJD(3) - t410 * t443) (t413 * t444 - t384 * t409 + (-t391 * t413 + t409 * t459) * qJD(4)) * pkin(4) + t453, t453, 0;];
JaD_transl  = t1;
