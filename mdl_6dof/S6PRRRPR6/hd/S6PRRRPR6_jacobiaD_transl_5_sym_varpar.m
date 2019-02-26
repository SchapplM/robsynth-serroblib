% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:13:27
% EndTime: 2019-02-26 20:13:28
% DurationCPUTime: 0.41s
% Computational Cost: add. (437->91), mult. (1367->161), div. (0->0), fcn. (1417->10), ass. (0->63)
t361 = sin(qJ(3));
t362 = sin(qJ(2));
t358 = sin(pkin(6));
t364 = cos(qJ(3));
t397 = t358 * t364;
t400 = cos(pkin(6));
t351 = t400 * t361 + t362 * t397;
t363 = cos(qJ(4));
t360 = sin(qJ(4));
t365 = cos(qJ(2));
t395 = t360 * t365;
t407 = -t351 * t363 + t358 * t395;
t391 = qJD(3) * t365;
t406 = (qJD(2) * t364 - qJD(4)) * t362 + t361 * t391;
t404 = pkin(9) + r_i_i_C(2);
t405 = -pkin(3) * t361 + t404 * t364;
t403 = -r_i_i_C(1) - pkin(4);
t401 = r_i_i_C(3) + qJ(5);
t398 = t358 * t361;
t396 = t360 * t364;
t394 = qJD(2) * t362;
t393 = qJD(2) * t365;
t392 = qJD(3) * t361;
t390 = qJD(4) * t364;
t388 = t358 * t394;
t387 = t358 * t393;
t385 = t362 * t400;
t384 = t365 * t400;
t357 = sin(pkin(11));
t382 = t357 * t385;
t359 = cos(pkin(11));
t381 = t359 * t384;
t342 = -qJD(2) * t381 + t357 * t394;
t346 = t357 * t362 - t381;
t380 = t346 * t390 - t342;
t348 = t357 * t384 + t359 * t362;
t344 = t348 * qJD(2);
t379 = t348 * t390 - t344;
t347 = t357 * t365 + t359 * t385;
t374 = -t347 * t364 + t359 * t398;
t378 = t346 * t360 - t363 * t374;
t349 = t359 * t365 - t382;
t337 = t349 * t364 + t357 * t398;
t377 = t337 * t363 + t348 * t360;
t376 = (qJD(2) - t390) * t365;
t375 = -t347 * t361 - t359 * t397;
t373 = -t349 * t361 + t357 * t397;
t372 = -t364 * pkin(3) - t404 * t361 - pkin(2);
t371 = qJD(3) * t405;
t370 = -t362 * t398 + t400 * t364;
t369 = t401 * t360 - t403 * t363 + pkin(3);
t343 = t347 * qJD(2);
t368 = qJD(4) * t347 - t343 * t364 + t346 * t392;
t345 = -qJD(2) * t382 + t359 * t393;
t367 = qJD(4) * t349 - t345 * t364 + t348 * t392;
t366 = qJD(5) * t360 + (t403 * t360 + t401 * t363) * qJD(4);
t339 = t370 * qJD(3) + t364 * t387;
t333 = t373 * qJD(3) - t344 * t364;
t331 = t375 * qJD(3) - t342 * t364;
t326 = -t407 * qJD(4) + t339 * t360 - t363 * t388;
t320 = t377 * qJD(4) + t333 * t360 - t345 * t363;
t318 = t378 * qJD(4) + t331 * t360 - t343 * t363;
t1 = [0 -(t348 * t396 + t349 * t363) * qJD(5) - t344 * pkin(8) - t403 * (t379 * t360 + t367 * t363) + t401 * (t367 * t360 - t379 * t363) - t348 * t371 + t372 * t345, t404 * t333 + t366 * t373 + t369 * (-t337 * qJD(3) + t344 * t361) t377 * qJD(5) + t401 * (t333 * t363 + t345 * t360 + (-t337 * t360 + t348 * t363) * qJD(4)) + t403 * t320, t320, 0; 0 -(t346 * t396 + t347 * t363) * qJD(5) - t342 * pkin(8) - t403 * (t380 * t360 + t368 * t363) + t401 * (t368 * t360 - t380 * t363) - t346 * t371 + t372 * t343, t404 * t331 + t366 * t375 + t369 * (t374 * qJD(3) + t342 * t361) t378 * qJD(5) + t401 * (t331 * t363 + t343 * t360 + (t346 * t363 + t360 * t374) * qJD(4)) + t403 * t318, t318, 0; 0 (-t403 * (t360 * t376 - t406 * t363) - t401 * (t406 * t360 + t363 * t376) - (t362 * t363 - t364 * t395) * qJD(5) + t405 * t391 + (t365 * pkin(8) + t372 * t362) * qJD(2)) * t358, t404 * t339 + t366 * t370 + t369 * (-t351 * qJD(3) - t361 * t387) -t407 * qJD(5) + t401 * (t360 * t388 + t339 * t363 + (-t358 * t363 * t365 - t351 * t360) * qJD(4)) + t403 * t326, t326, 0;];
JaD_transl  = t1;
