% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR10_jacobiaD_transl_3_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_transl_3_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_transl_3_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobiaD_transl_3_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiaD_transl_3_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:19
% EndTime: 2018-11-23 11:27:19
% DurationCPUTime: 0.43s
% Computational Cost: add. (724->100), mult. (963->151), div. (0->0), fcn. (712->18), ass. (0->64)
t459 = pkin(7) - qJ(3);
t453 = cos(t459);
t427 = pkin(7) + qJ(3);
t457 = cos(t427) / 0.2e1;
t407 = t457 - t453 / 0.2e1;
t461 = pkin(6) - qJ(2);
t454 = cos(t461);
t460 = pkin(6) + qJ(2);
t470 = cos(t460) / 0.2e1;
t409 = t470 - t454 / 0.2e1;
t471 = pkin(11) + r_i_i_C(3);
t473 = t471 * cos(pkin(7)) - r_i_i_C(1) * t407 + pkin(10);
t447 = sin(t460) / 0.2e1;
t452 = sin(t461);
t406 = t447 - t452 / 0.2e1;
t434 = sin(qJ(1));
t436 = cos(qJ(2));
t437 = cos(qJ(1));
t388 = t406 * t437 + t434 * t436;
t458 = t454 / 0.2e1 + t470;
t400 = t458 * qJD(2);
t433 = sin(qJ(2));
t464 = qJD(2) * t437;
t472 = qJD(1) * t388 + t434 * t400 + t433 * t464;
t469 = sin(t427);
t451 = sin(t459);
t446 = t451 / 0.2e1;
t456 = t469 / 0.2e1;
t403 = t456 + t446;
t395 = t403 * qJD(3);
t468 = t395 * t434;
t467 = qJD(1) * t434;
t466 = qJD(1) * t437;
t465 = qJD(2) * t434;
t432 = sin(qJ(3));
t463 = qJD(3) * t432;
t435 = cos(qJ(3));
t462 = qJD(3) * t435;
t450 = t437 * t458;
t392 = t406 * t434 - t437 * t436;
t404 = t456 - t451 / 0.2e1;
t445 = t404 * r_i_i_C(1) - t471 * sin(pkin(7));
t443 = t437 * t433 + t434 * t458;
t408 = t457 + t453 / 0.2e1;
t442 = t408 * r_i_i_C(2) + t445;
t397 = t408 * qJD(3);
t441 = -t392 * t463 + t397 * t443 + t435 * t472;
t417 = t433 * t465;
t385 = -t406 * t467 - t417 + (qJD(1) * t436 + t400) * t437;
t387 = t434 * t433 - t450;
t396 = (t446 - t469 / 0.2e1) * qJD(3);
t440 = t385 * t432 + t387 * t396 + t388 * t462;
t439 = -t385 * t435 + t387 * t397 + t388 * t463;
t405 = t447 + t452 / 0.2e1;
t438 = t406 * qJD(2);
t431 = cos(pkin(6));
t429 = sin(pkin(6));
t401 = t409 * qJD(2);
t399 = t405 * qJD(2);
t398 = t407 * qJD(3);
t384 = qJD(1) * t443 + t436 * t465 + t437 * t438;
t381 = -qJD(1) * t450 + t433 * t467 + t434 * t438 - t436 * t464;
t380 = t392 * t462 + t381 * t408 + t472 * t432 - t443 * t396 + (t398 * t434 + t403 * t466) * t429;
t1 = [t439 * r_i_i_C(1) + t440 * r_i_i_C(2) - t385 * pkin(2) - pkin(1) * t466 + t442 * t384 + ((r_i_i_C(1) * t395 + r_i_i_C(2) * t398) * t437 + (-r_i_i_C(2) * t403 - t473) * t467) * t429 (t381 * t435 + t392 * t397 + t443 * t463) * r_i_i_C(1) + (-t381 * t432 + t392 * t396 + t443 * t462) * r_i_i_C(2) + t381 * pkin(2) + t442 * t472, t380 * r_i_i_C(1) + (-t381 * t404 + (t407 * t466 - t468) * t429 + t441) * r_i_i_C(2), 0, 0, 0; (t429 * t468 - t441) * r_i_i_C(1) + t380 * r_i_i_C(2) - t472 * pkin(2) + t445 * t381 + (t429 * t437 * t473 - t434 * pkin(1)) * qJD(1) (-t384 * t435 + t387 * t463 - t388 * t397) * r_i_i_C(1) + (t384 * t432 + t387 * t462 - t388 * t396) * r_i_i_C(2) - t384 * pkin(2) + t442 * (qJD(1) * t392 - t437 * t400 + t417) (-t384 * t408 - t440) * r_i_i_C(1) + (t384 * t404 + t439) * r_i_i_C(2) + ((-t398 * t437 + t403 * t467) * r_i_i_C(1) + (t395 * t437 + t407 * t467) * r_i_i_C(2)) * t429, 0, 0, 0; 0 (t397 * t409 + t401 * t435 - t405 * t463) * r_i_i_C(1) + (t396 * t409 - t401 * t432 - t405 * t462) * r_i_i_C(2) + t401 * pkin(2) - t442 * t399 (t396 * t405 + t398 * t431 - t399 * t432 + t401 * t408 + t409 * t462) * r_i_i_C(1) + (-t395 * t431 - t397 * t405 - t399 * t435 - t401 * t404 - t409 * t463) * r_i_i_C(2), 0, 0, 0;];
JaD_transl  = t1;
