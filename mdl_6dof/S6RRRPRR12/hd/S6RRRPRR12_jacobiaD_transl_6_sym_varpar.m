% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR12_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR12_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:22:37
% EndTime: 2019-02-26 22:22:37
% DurationCPUTime: 0.60s
% Computational Cost: add. (996->114), mult. (1956->179), div. (0->0), fcn. (1978->14), ass. (0->76)
t419 = sin(qJ(3));
t422 = cos(qJ(3));
t416 = pkin(12) + qJ(5);
t413 = cos(t416);
t398 = pkin(5) * t413 + cos(pkin(12)) * pkin(4) + pkin(3);
t414 = qJ(6) + t416;
t410 = sin(t414);
t411 = cos(t414);
t437 = t411 * r_i_i_C(1) - t410 * r_i_i_C(2);
t433 = t398 + t437;
t467 = r_i_i_C(3) + pkin(11) + pkin(10) + qJ(4);
t412 = sin(t416);
t417 = qJD(5) + qJD(6);
t436 = t410 * r_i_i_C(1) + t411 * r_i_i_C(2);
t466 = pkin(5) * qJD(5);
t472 = t412 * t466 + t436 * t417;
t424 = t472 * t422 + (t433 * t419 - t467 * t422) * qJD(3) - t419 * qJD(4);
t421 = sin(qJ(1));
t423 = cos(qJ(2));
t465 = cos(pkin(6));
t469 = cos(qJ(1));
t438 = t465 * t469;
t420 = sin(qJ(2));
t447 = t421 * t465;
t439 = t420 * t447;
t448 = t469 * qJD(1);
t457 = qJD(2) * t420;
t384 = -qJD(1) * t439 - t421 * t457 + (qJD(2) * t438 + t448) * t423;
t418 = sin(pkin(6));
t453 = t418 * t469;
t475 = -qJD(3) * t453 + t384;
t395 = t420 * t438 + t421 * t423;
t464 = t418 * t421;
t474 = qJD(1) * t464 - qJD(3) * t395;
t378 = t474 * t419 + t475 * t422;
t470 = t467 * t419 + t433 * t422 + pkin(2);
t468 = -pkin(9) - pkin(5) * t412 - sin(pkin(12)) * pkin(4);
t463 = t418 * t422;
t462 = t418 * t423;
t461 = t421 * t420;
t434 = t423 * t438;
t452 = t469 * t423;
t381 = -qJD(1) * t434 - qJD(2) * t452 + (qJD(2) * t465 + qJD(1)) * t461;
t397 = t452 - t439;
t391 = t397 * t422 + t419 * t464;
t444 = -t391 * t417 - t381;
t396 = t469 * t420 + t423 * t447;
t382 = t395 * qJD(1) + t396 * qJD(2);
t390 = -t397 * t419 + t421 * t463;
t441 = t418 * t448;
t376 = t390 * qJD(3) - t382 * t422 + t419 * t441;
t446 = t396 * t417 + t376;
t371 = -t446 * t410 + t444 * t411;
t372 = t444 * t410 + t446 * t411;
t460 = t371 * r_i_i_C(1) - t372 * r_i_i_C(2);
t383 = t396 * qJD(1) + t395 * qJD(2);
t388 = t395 * t422 - t419 * t453;
t443 = t388 * t417 - t383;
t394 = -t434 + t461;
t445 = -t394 * t417 - t378;
t459 = (t445 * t410 - t443 * t411) * r_i_i_C(1) + (t443 * t410 + t445 * t411) * r_i_i_C(2);
t393 = t465 * t419 + t420 * t463;
t450 = t418 * t457;
t432 = -t393 * t417 + t450;
t428 = -t418 * t420 * t419 + t465 * t422;
t449 = qJD(2) * t462;
t386 = t428 * qJD(3) + t422 * t449;
t435 = t417 * t462 - t386;
t458 = (t435 * t410 + t432 * t411) * r_i_i_C(1) + (-t432 * t410 + t435 * t411) * r_i_i_C(2);
t431 = t436 - t468;
t429 = t395 * t419 + t422 * t453;
t377 = t475 * t419 - t474 * t422;
t426 = t413 * t466 + t437 * t417;
t385 = t393 * qJD(3) + t419 * t449;
t375 = t391 * qJD(3) - t382 * t419 - t422 * t441;
t1 = [-t429 * qJD(4) - t384 * pkin(2) - t433 * t378 + ((t388 * t410 - t394 * t411) * r_i_i_C(1) + (t388 * t411 + t394 * t410) * r_i_i_C(2)) * t417 - t431 * t383 - t467 * t377 + (-t469 * pkin(1) - pkin(8) * t464) * qJD(1) + (t388 * t412 - t394 * t413) * t466, t470 * t381 - t431 * t382 + t424 * t396 + t426 * t397, t391 * qJD(4) - t433 * t375 + t467 * t376 - t390 * t472, t375 (-t376 * t412 - t381 * t413 + (-t391 * t413 - t396 * t412) * qJD(5)) * pkin(5) + t460, t460; -t382 * pkin(2) + t372 * r_i_i_C(1) + t371 * r_i_i_C(2) - t390 * qJD(4) + t376 * t398 + t468 * t381 + t467 * t375 + (-pkin(1) * t421 + pkin(8) * t453) * qJD(1) + (-t391 * t412 + t396 * t413) * t466, -t383 * t470 + t431 * t384 + t424 * t394 + t426 * t395, t388 * qJD(4) - t433 * t377 + t467 * t378 + t429 * t472, t377 (-t378 * t412 + t383 * t413 + (-t388 * t413 - t394 * t412) * qJD(5)) * pkin(5) + t459, t459; 0 ((-qJD(2) * t470 + t426) * t420 + (t431 * qJD(2) - t424) * t423) * t418, t393 * qJD(4) - t433 * t385 + t467 * t386 - t428 * t472, t385 (t413 * t450 - t386 * t412 + (-t393 * t413 + t412 * t462) * qJD(5)) * pkin(5) + t458, t458;];
JaD_transl  = t1;
