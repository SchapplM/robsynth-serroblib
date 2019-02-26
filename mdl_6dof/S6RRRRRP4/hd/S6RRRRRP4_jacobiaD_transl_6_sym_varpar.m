% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:11
% EndTime: 2019-02-26 22:41:12
% DurationCPUTime: 0.53s
% Computational Cost: add. (959->97), mult. (946->131), div. (0->0), fcn. (771->10), ass. (0->75)
t357 = qJ(4) + qJ(5);
t351 = sin(t357);
t435 = r_i_i_C(3) + qJ(6);
t438 = t351 * t435;
t353 = cos(t357);
t355 = qJD(4) + qJD(5);
t364 = cos(qJ(1));
t414 = t355 * t364;
t397 = t351 * t414;
t361 = sin(qJ(1));
t412 = qJD(1) * t361;
t437 = t353 * t412 + t397;
t436 = pkin(5) + r_i_i_C(1);
t358 = qJ(2) + qJ(3);
t352 = sin(t358);
t365 = -pkin(10) - pkin(9);
t359 = sin(qJ(4));
t421 = pkin(4) * qJD(4);
t408 = t359 * t421;
t356 = qJD(2) + qJD(3);
t420 = t352 * t356;
t433 = t352 * t408 + t365 * t420;
t410 = qJD(6) * t351;
t432 = -t408 + t410;
t362 = cos(qJ(4));
t349 = pkin(4) * t362 + pkin(3);
t354 = cos(t358);
t360 = sin(qJ(2));
t422 = pkin(2) * qJD(2);
t407 = t360 * t422;
t424 = r_i_i_C(2) - t365;
t426 = pkin(4) * t359;
t431 = (t424 * t356 + t432) * t354 + (pkin(8) + pkin(7) + t426) * qJD(1) - t349 * t420 - t407;
t427 = pkin(2) * t360;
t425 = r_i_i_C(2) * t354;
t419 = t353 * t355;
t418 = t353 * t364;
t417 = t354 * t356;
t416 = t354 * t361;
t415 = t355 * t361;
t413 = t356 * t361;
t411 = qJD(1) * t364;
t409 = t353 * qJD(6);
t406 = t362 * t421;
t405 = t436 * t351;
t404 = t352 * t419;
t403 = t352 * t413;
t402 = t364 * t420;
t398 = t351 * t415;
t396 = t353 * t414;
t395 = t435 * t353 * t417 + t352 * t409;
t393 = t354 * t412;
t391 = t435 * t355;
t388 = qJD(1) * t354 - qJD(4);
t383 = t436 * t352 * t398 + t433 * t361 + t411 * t425;
t382 = (-qJD(4) * t354 + qJD(1)) * t362;
t380 = t351 * t361 + t354 * t418;
t378 = t351 * t411 + t353 * t415;
t377 = -t353 * t436 - t438;
t376 = -t353 * t391 - t410;
t311 = t351 * t402 - t354 * t396 - t398 + (t351 * t416 + t418) * qJD(1);
t312 = t354 * t397 + (t393 + t402) * t353 - t378;
t375 = qJD(6) * t380 + t436 * t311 - t435 * t312;
t313 = -t351 * t403 + t378 * t354 - t437;
t314 = t380 * qJD(1) - t353 * t403 - t354 * t398 - t396;
t374 = -qJD(6) * (t351 * t364 - t353 * t416) + t435 * t314 - t436 * t313;
t373 = -t349 + t377;
t372 = t433 * t364 + t365 * t393 + (t436 * t437 + (t349 + t438) * t412) * t352;
t371 = t373 * t352 - t354 * t365;
t363 = cos(qJ(2));
t370 = t406 - t409 + (-pkin(2) * t363 - t349 * t354 - t424 * t352 - pkin(1)) * qJD(1);
t369 = t376 * t352 + (-r_i_i_C(2) * t352 + t373 * t354) * t356;
t368 = -t363 * t422 + t369;
t367 = r_i_i_C(2) * t417 + t371 * t356 + (-t355 * t405 + t435 * t419 + t432) * t354;
t1 = [-t313 * t435 - t314 * t436 - t431 * t361 + t370 * t364 (-t425 + t427) * t412 + t368 * t364 + t372, -r_i_i_C(2) * t393 + t369 * t364 + t372 (t364 * t382 + (t388 * t361 + t402) * t359) * pkin(4) + t375, t375, -t311; -t311 * t435 - t312 * t436 + t370 * t361 + t431 * t364 (t371 - t427) * t411 + t368 * t361 + t383 (-t365 * t411 + t373 * t413) * t354 + ((-r_i_i_C(2) * t356 + t376) * t361 + t373 * t411) * t352 + t383 (t361 * t382 + (-t388 * t364 + t403) * t359) * pkin(4) + t374, t374, t313; 0, t367 - t407, t367 (-t405 - t426) * t417 + (t377 * t355 - t406) * t352 + t395, -t436 * t404 + (-t352 * t391 - t417 * t436) * t351 + t395, t351 * t417 + t404;];
JaD_transl  = t1;
