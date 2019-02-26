% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:42
% EndTime: 2019-02-26 22:07:42
% DurationCPUTime: 0.51s
% Computational Cost: add. (659->97), mult. (1950->157), div. (0->0), fcn. (1950->10), ass. (0->61)
t370 = sin(qJ(3));
t374 = cos(qJ(3));
t369 = sin(qJ(6));
t373 = cos(qJ(6));
t389 = -t373 * r_i_i_C(1) + t369 * r_i_i_C(2);
t413 = pkin(5) + qJ(4);
t385 = -t389 + t413;
t398 = r_i_i_C(3) + pkin(10) + pkin(4) + pkin(3);
t388 = t369 * r_i_i_C(1) + t373 * r_i_i_C(2);
t416 = t388 * qJD(6) - qJD(4);
t377 = t416 * t370 + (t398 * t370 - t385 * t374) * qJD(3);
t371 = sin(qJ(2));
t372 = sin(qJ(1));
t375 = cos(qJ(2));
t376 = cos(qJ(1));
t411 = cos(pkin(6));
t391 = t376 * t411;
t356 = t372 * t371 - t375 * t391;
t357 = t371 * t391 + t372 * t375;
t368 = sin(pkin(6));
t403 = t368 * t376;
t418 = t357 * t370 + t374 * t403;
t420 = t356 * t373 + t369 * t418;
t419 = -t356 * t369 + t373 * t418;
t392 = t372 * t411;
t390 = t371 * t392;
t400 = qJD(2) * t371;
t401 = t376 * t375;
t345 = -qJD(1) * t390 - t372 * t400 + (qJD(2) * t411 + qJD(1)) * t401;
t406 = t368 * t372;
t397 = t370 * t406;
t417 = -qJD(1) * t397 + qJD(3) * t418 - t345 * t374;
t414 = t385 * t370 + t398 * t374 + pkin(2);
t412 = -pkin(9) + qJ(5);
t381 = t390 - t401;
t407 = t381 * t370;
t405 = t368 * t374;
t404 = t368 * t375;
t402 = t370 * t376;
t399 = qJD(3) * t374;
t396 = t368 * t402;
t395 = qJD(1) * t405;
t394 = qJD(2) * t404;
t393 = t368 * t400;
t386 = -t374 * t381 + t397;
t384 = t388 + t412;
t358 = t376 * t371 + t375 * t392;
t354 = t368 * t371 * t370 - t411 * t374;
t382 = t411 * t370 + t371 * t405;
t379 = t389 * qJD(6) - qJD(5);
t338 = -qJD(3) * t396 + t345 * t370 + t357 * t399 - t372 * t395;
t351 = -t372 * t405 - t407;
t346 = t382 * qJD(3) + t370 * t394;
t344 = t358 * qJD(1) + t357 * qJD(2);
t343 = t357 * qJD(1) + t358 * qJD(2);
t342 = t356 * qJD(1) + t381 * qJD(2);
t337 = -t343 * t374 + qJD(3) * t407 + (qJD(1) * t402 + t372 * t399) * t368;
t336 = t386 * qJD(3) - t343 * t370 - t376 * t395;
t335 = t336 * t373 + t342 * t369 + (-t351 * t369 - t358 * t373) * qJD(6);
t334 = -t336 * t369 + t342 * t373 + (-t351 * t373 + t358 * t369) * qJD(6);
t1 = [-t345 * pkin(2) - t418 * qJD(4) + t356 * qJD(5) + t384 * t344 - t385 * t338 + (t420 * r_i_i_C(1) + t419 * r_i_i_C(2)) * qJD(6) + (-t376 * pkin(1) - pkin(8) * t406) * qJD(1) + t398 * t417, t414 * t342 + t384 * t343 + t377 * t358 - t379 * t381, -t398 * t336 + t385 * t337 - t386 * t416, t336, t342, t334 * r_i_i_C(1) - t335 * r_i_i_C(2); -t343 * pkin(2) + t335 * r_i_i_C(1) + t334 * r_i_i_C(2) + t351 * qJD(4) - t358 * qJD(5) + t412 * t342 + t413 * t336 + (-pkin(1) * t372 + pkin(8) * t403) * qJD(1) + t398 * t337, -t344 * t414 - t384 * t345 + t377 * t356 + t379 * t357, -t416 * (t357 * t374 - t396) - t385 * t417 - t398 * t338, t338, -t344 (-t338 * t369 - t344 * t373) * r_i_i_C(1) + (-t338 * t373 + t344 * t369) * r_i_i_C(2) + (-t419 * r_i_i_C(1) + t420 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t414 + t379) * t371 + (-t384 * qJD(2) - t377) * t375) * t368, -t416 * t382 + t385 * (-t354 * qJD(3) + t374 * t394) - t398 * t346, t346, -t393 (-t346 * t369 - t373 * t393) * r_i_i_C(1) + (-t346 * t373 + t369 * t393) * r_i_i_C(2) + ((-t354 * t373 - t369 * t404) * r_i_i_C(1) + (t354 * t369 - t373 * t404) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
