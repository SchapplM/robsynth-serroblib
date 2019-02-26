% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRP6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:17
% EndTime: 2019-02-26 19:53:17
% DurationCPUTime: 0.35s
% Computational Cost: add. (449->93), mult. (1404->158), div. (0->0), fcn. (1453->10), ass. (0->64)
t372 = sin(qJ(4));
t375 = cos(qJ(4));
t414 = pkin(9) + r_i_i_C(2);
t382 = -t372 * pkin(4) + t414 * t375 - qJ(3);
t415 = -pkin(8) - pkin(2);
t413 = -r_i_i_C(1) - pkin(5);
t412 = r_i_i_C(3) + qJ(6);
t368 = sin(pkin(6));
t411 = t368 * t372;
t410 = t368 * t375;
t376 = cos(qJ(2));
t409 = t368 * t376;
t370 = cos(pkin(6));
t373 = sin(qJ(2));
t408 = t370 * t373;
t407 = t370 * t376;
t374 = cos(qJ(5));
t406 = t373 * t374;
t405 = qJD(2) * t373;
t404 = qJD(2) * t376;
t403 = qJD(4) * t375;
t402 = qJD(5) * t372;
t371 = sin(qJ(5));
t401 = qJD(6) * t371;
t400 = qJD(6) * t374;
t367 = sin(pkin(10));
t399 = t367 * t411;
t398 = t367 * t405;
t397 = t368 * t404;
t369 = cos(pkin(10));
t396 = t369 * t404;
t395 = t368 * t405;
t393 = qJD(2) + t402;
t357 = t367 * t376 + t369 * t408;
t353 = t357 * qJD(2);
t392 = t357 * t402 + t353;
t355 = -t370 * t398 + t396;
t359 = -t367 * t408 + t369 * t376;
t391 = t359 * t402 + t355;
t358 = t367 * t407 + t369 * t373;
t343 = t358 * t372 + t367 * t410;
t390 = t343 * t374 + t359 * t371;
t356 = t367 * t373 - t369 * t407;
t387 = -t356 * t372 + t369 * t410;
t389 = t357 * t371 - t374 * t387;
t388 = (qJD(2) * t372 + qJD(5)) * t376;
t386 = t356 * t375 + t369 * t411;
t383 = -t370 * t375 + t372 * t409;
t385 = t368 * t371 * t373 - t374 * t383;
t384 = t370 * t372 + t375 * t409;
t381 = t412 * t371 - t413 * t374 + pkin(4);
t352 = -t370 * t396 + t398;
t380 = -qJD(5) * t356 - t352 * t372 + t357 * t403;
t354 = t358 * qJD(2);
t379 = -qJD(5) * t358 - t354 * t372 + t359 * t403;
t378 = t401 + (t413 * t371 + t412 * t374) * qJD(5);
t377 = t372 * t401 + qJD(3) + (pkin(4) * t375 + t414 * t372) * qJD(4);
t346 = t384 * qJD(4) - t372 * t395;
t340 = t386 * qJD(4) + t353 * t372;
t338 = qJD(4) * t399 - t355 * t372 - t358 * t403;
t334 = t385 * qJD(5) - t346 * t371 - t374 * t397;
t328 = t389 * qJD(5) + t340 * t371 + t352 * t374;
t326 = t390 * qJD(5) - t338 * t371 + t354 * t374;
t1 = [0, t358 * t400 + t415 * t355 - t413 * (-t391 * t371 + t379 * t374) + t412 * (t379 * t371 + t391 * t374) + t377 * t359 + t382 * t354, t355, -t414 * t338 + t378 * (t358 * t375 - t399) + t381 * (-t343 * qJD(4) + t355 * t375) t390 * qJD(6) + t412 * (-t338 * t374 - t354 * t371 + (-t343 * t371 + t359 * t374) * qJD(5)) + t413 * t326, t326; 0, t356 * t400 + t415 * t353 - t413 * (-t392 * t371 + t380 * t374) + t412 * (t380 * t371 + t392 * t374) + t377 * t357 + t382 * t352, t353, t414 * t340 + t378 * t386 + t381 * (t387 * qJD(4) + t353 * t375) t389 * qJD(6) + t412 * (t340 * t374 - t352 * t371 + (t357 * t374 + t371 * t387) * qJD(5)) + t413 * t328, t328; 0 (-t413 * (t374 * t388 + (-t393 * t371 + t374 * t403) * t373) + t412 * (t393 * t406 + (t373 * t403 + t388) * t371) - t376 * t400 + t377 * t373 + (t415 * t373 - t382 * t376) * qJD(2)) * t368, t395, -t414 * t346 - t378 * t384 + t381 * (t383 * qJD(4) + t375 * t395) t385 * qJD(6) + t412 * (t371 * t397 - t346 * t374 + (t368 * t406 + t371 * t383) * qJD(5)) + t413 * t334, t334;];
JaD_transl  = t1;
