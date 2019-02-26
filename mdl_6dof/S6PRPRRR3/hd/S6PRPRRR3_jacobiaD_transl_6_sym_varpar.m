% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:50
% EndTime: 2019-02-26 19:54:50
% DurationCPUTime: 0.46s
% Computational Cost: add. (801->93), mult. (1063->161), div. (0->0), fcn. (1055->13), ass. (0->62)
t408 = pkin(10) + r_i_i_C(3);
t367 = sin(qJ(6));
t369 = cos(qJ(6));
t381 = t369 * r_i_i_C(1) - t367 * r_i_i_C(2);
t414 = pkin(5) + t381;
t394 = qJD(6) * t369;
t395 = qJD(6) * t367;
t413 = -r_i_i_C(1) * t395 - t394 * r_i_i_C(2);
t363 = pkin(12) + qJ(4);
t361 = qJ(5) + t363;
t357 = sin(t361);
t358 = cos(t361);
t359 = sin(t363);
t364 = qJD(4) + qJD(5);
t412 = -pkin(4) * qJD(4) * t359 - (t357 * t414 - t408 * t358) * t364;
t368 = sin(qJ(2));
t370 = cos(qJ(2));
t365 = sin(pkin(11));
t404 = cos(pkin(6));
t387 = t365 * t404;
t403 = cos(pkin(11));
t350 = -t368 * t387 + t403 * t370;
t410 = t367 * r_i_i_C(1) + t369 * r_i_i_C(2);
t402 = t357 * t364;
t401 = t358 * t364;
t366 = sin(pkin(6));
t400 = t365 * t366;
t399 = t366 * t368;
t398 = t366 * t370;
t397 = qJD(2) * t368;
t396 = qJD(6) * t358;
t391 = t357 * t399;
t390 = t358 * t399;
t389 = qJD(2) * t398;
t388 = t366 * t397;
t386 = t366 * t403;
t382 = t358 * t386;
t349 = t403 * t368 + t370 * t387;
t345 = t349 * qJD(2);
t380 = t364 * t400 - t345;
t379 = t404 * t403;
t377 = t370 * t379;
t376 = t404 * t364 + t389;
t348 = t365 * t370 + t368 * t379;
t360 = cos(t363);
t375 = -t408 * t357 - t414 * t358 - pkin(4) * t360 - cos(pkin(12)) * pkin(3) - pkin(2);
t329 = -t350 * t402 + t380 * t358;
t374 = t413 * (-t350 * t357 + t358 * t400) + t408 * t329 + t414 * (-t350 * t401 - t380 * t357);
t343 = -qJD(2) * t377 + t365 * t397;
t327 = -t343 * t358 - t348 * t402 - t364 * t382;
t373 = t413 * (-t348 * t357 - t382) + t408 * t327 + t414 * (-t348 * t401 + (t364 * t386 + t343) * t357);
t334 = t376 * t358 - t364 * t391;
t372 = t413 * (t404 * t358 - t391) + t408 * t334 + t414 * (-t376 * t357 - t364 * t390);
t371 = t410 * t396 - t412;
t362 = -pkin(9) - pkin(8) - qJ(3);
t347 = t365 * t368 - t377;
t346 = t350 * qJD(2);
t344 = t348 * qJD(2);
t342 = t404 * t357 + t390;
t338 = t350 * t358 + t357 * t400;
t336 = t348 * t358 - t357 * t386;
t1 = [0 (-t345 * t367 + t350 * t394) * r_i_i_C(1) + (-t345 * t369 - t350 * t395) * r_i_i_C(2) + t345 * t362 + t350 * qJD(3) + t375 * t346 + t371 * t349, t346 (t345 * t359 + (-t350 * t360 - t359 * t400) * qJD(4)) * pkin(4) + t374, t374 (-t329 * t367 + t346 * t369) * r_i_i_C(1) + (-t329 * t369 - t346 * t367) * r_i_i_C(2) + ((-t338 * t369 - t349 * t367) * r_i_i_C(1) + (t338 * t367 - t349 * t369) * r_i_i_C(2)) * qJD(6); 0 (-t343 * t367 + t348 * t394) * r_i_i_C(1) + (-t343 * t369 - t348 * t395) * r_i_i_C(2) + t343 * t362 + t348 * qJD(3) + t375 * t344 + t371 * t347, t344 (t343 * t359 + (-t348 * t360 + t359 * t386) * qJD(4)) * pkin(4) + t373, t373 (-t327 * t367 + t344 * t369) * r_i_i_C(1) + (-t327 * t369 - t344 * t367) * r_i_i_C(2) + ((-t336 * t369 - t347 * t367) * r_i_i_C(1) + (t336 * t367 - t347 * t369) * r_i_i_C(2)) * qJD(6); 0 ((t375 * qJD(2) + t381 * qJD(6) + qJD(3)) * t368 + (-qJD(2) * t362 + t410 * (qJD(2) - t396) + t412) * t370) * t366, t388 (-t359 * t389 + (-t404 * t359 - t360 * t399) * qJD(4)) * pkin(4) + t372, t372 (-t334 * t367 + t369 * t388) * r_i_i_C(1) + (-t334 * t369 - t367 * t388) * r_i_i_C(2) + ((-t342 * t369 + t367 * t398) * r_i_i_C(1) + (t342 * t367 + t369 * t398) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
