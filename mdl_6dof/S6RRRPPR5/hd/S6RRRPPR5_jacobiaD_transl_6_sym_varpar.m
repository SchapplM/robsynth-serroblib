% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:52
% EndTime: 2019-02-26 22:05:53
% DurationCPUTime: 0.61s
% Computational Cost: add. (866->120), mult. (1643->191), div. (0->0), fcn. (1644->14), ass. (0->72)
t402 = sin(qJ(1));
t443 = cos(pkin(6));
t415 = qJD(2) * t443 + qJD(1);
t401 = sin(qJ(2));
t422 = t402 * t443;
t419 = t401 * t422;
t432 = qJD(2) * t401;
t404 = cos(qJ(2));
t405 = cos(qJ(1));
t434 = t405 * t404;
t364 = -qJD(1) * t419 - t402 * t432 + t415 * t434;
t397 = sin(pkin(6));
t436 = t397 * t405;
t452 = -qJD(3) * t436 + t364;
t421 = t405 * t443;
t375 = t401 * t421 + t402 * t404;
t433 = qJD(1) * t397;
t428 = t402 * t433;
t451 = -qJD(3) * t375 + t428;
t395 = qJ(3) + pkin(11);
t391 = sin(t395);
t393 = cos(t395);
t368 = t375 * t393 - t391 * t436;
t418 = t404 * t421;
t435 = t401 * t402;
t374 = -t418 + t435;
t394 = pkin(12) + qJ(6);
t390 = sin(t394);
t392 = cos(t394);
t450 = t368 * t390 - t374 * t392;
t449 = t368 * t392 + t374 * t390;
t400 = sin(qJ(3));
t416 = t390 * r_i_i_C(1) + t392 * r_i_i_C(2);
t412 = qJD(6) * t416;
t388 = cos(pkin(12)) * pkin(5) + pkin(4);
t417 = t392 * r_i_i_C(1) - t390 * r_i_i_C(2);
t414 = t388 + t417;
t444 = r_i_i_C(3) + pkin(10) + qJ(5);
t406 = (t400 * pkin(3) + t414 * t391 - t444 * t393) * qJD(3) - t391 * qJD(5) + t393 * t412;
t358 = t451 * t391 + t452 * t393;
t403 = cos(qJ(3));
t389 = pkin(3) * t403 + pkin(2);
t447 = t444 * t391 + t414 * t393 + t389;
t440 = t397 * t401;
t439 = t397 * t402;
t438 = t397 * t403;
t437 = t397 * t404;
t431 = qJD(2) * t404;
t427 = t405 * t433;
t426 = t397 * t431;
t424 = t397 * t432;
t423 = -sin(pkin(12)) * pkin(5) - qJ(4) - pkin(9);
t413 = t375 * t391 + t393 * t436;
t377 = -t419 + t434;
t370 = -t377 * t391 + t393 * t439;
t371 = t377 * t393 + t391 * t439;
t376 = t405 * t401 + t404 * t422;
t410 = -t391 * t440 + t443 * t393;
t373 = t443 * t391 + t393 * t440;
t409 = t416 - t423;
t408 = t417 * qJD(6) + qJD(4);
t357 = t452 * t391 - t451 * t393;
t366 = t410 * qJD(3) + t393 * t426;
t365 = t373 * qJD(3) + t391 * t426;
t363 = t376 * qJD(1) + t375 * qJD(2);
t362 = t375 * qJD(1) + t376 * qJD(2);
t361 = -qJD(1) * t418 - t405 * t431 + t415 * t435;
t356 = t370 * qJD(3) - t362 * t393 + t391 * t427;
t355 = t371 * qJD(3) - t362 * t391 - t393 * t427;
t354 = t356 * t392 - t361 * t390 + (-t371 * t390 + t376 * t392) * qJD(6);
t353 = -t356 * t390 - t361 * t392 + (-t371 * t392 - t376 * t390) * qJD(6);
t1 = [-t413 * qJD(5) - t364 * t389 - t374 * qJD(4) - t414 * t358 - t409 * t363 - t444 * t357 + (t450 * r_i_i_C(1) + t449 * r_i_i_C(2)) * qJD(6) + (-t405 * pkin(1) - pkin(8) * t439) * qJD(1) + (-t400 * t428 + (t375 * t400 + t403 * t436) * qJD(3)) * pkin(3), t447 * t361 - t409 * t362 + t406 * t376 + t408 * t377, t371 * qJD(5) + t444 * t356 - t370 * t412 - t414 * t355 + (t403 * t427 + t362 * t400 + (-t377 * t403 - t400 * t439) * qJD(3)) * pkin(3), -t361, t355, r_i_i_C(1) * t353 - t354 * r_i_i_C(2); t354 * r_i_i_C(1) + t353 * r_i_i_C(2) + t376 * qJD(4) - t370 * qJD(5) + t356 * t388 - t362 * t389 + t423 * t361 + t444 * t355 + (-pkin(1) * t402 + pkin(8) * t436) * qJD(1) + (t400 * t427 + (-t377 * t400 + t402 * t438) * qJD(3)) * pkin(3), -t363 * t447 + t409 * t364 + t406 * t374 + t408 * t375, t368 * qJD(5) + t444 * t358 + t413 * t412 - t414 * t357 + (t403 * t428 - t364 * t400 + (-t375 * t403 + t400 * t436) * qJD(3)) * pkin(3), t363, t357 (-t358 * t390 + t363 * t392) * r_i_i_C(1) + (-t358 * t392 - t363 * t390) * r_i_i_C(2) + (-t449 * r_i_i_C(1) + t450 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t447 + t408) * t401 + (t409 * qJD(2) - t406) * t404) * t397, t373 * qJD(5) + t444 * t366 - t410 * t412 - t414 * t365 + (-t400 * t426 + (-t443 * t400 - t401 * t438) * qJD(3)) * pkin(3), t424, t365 (-t366 * t390 + t392 * t424) * r_i_i_C(1) + (-t366 * t392 - t390 * t424) * r_i_i_C(2) + ((-t373 * t392 + t390 * t437) * r_i_i_C(1) + (t373 * t390 + t392 * t437) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
