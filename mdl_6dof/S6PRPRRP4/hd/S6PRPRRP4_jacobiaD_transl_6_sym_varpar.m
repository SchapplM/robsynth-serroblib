% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:03
% EndTime: 2019-02-26 19:52:03
% DurationCPUTime: 0.40s
% Computational Cost: add. (681->97), mult. (1392->166), div. (0->0), fcn. (1444->11), ass. (0->66)
t370 = pkin(11) + qJ(4);
t368 = sin(t370);
t369 = cos(t370);
t372 = sin(pkin(6));
t376 = sin(qJ(2));
t409 = t372 * t376;
t414 = cos(pkin(6));
t353 = t368 * t414 + t369 * t409;
t377 = cos(qJ(5));
t375 = sin(qJ(5));
t378 = cos(qJ(2));
t408 = t375 * t378;
t421 = -t353 * t377 + t372 * t408;
t418 = pkin(9) + r_i_i_C(2);
t420 = -pkin(4) * t368 + t369 * t418;
t404 = qJD(4) * t378;
t419 = (qJD(2) * t369 - qJD(5)) * t376 + t368 * t404;
t417 = -r_i_i_C(1) - pkin(5);
t415 = r_i_i_C(3) + qJ(6);
t412 = t369 * t375;
t371 = sin(pkin(10));
t411 = t371 * t372;
t373 = cos(pkin(10));
t410 = t372 * t373;
t407 = qJD(2) * t376;
t406 = qJD(2) * t378;
t405 = qJD(4) * t368;
t403 = qJD(5) * t369;
t401 = t372 * t406;
t400 = t372 * t407;
t398 = t376 * t414;
t397 = t378 * t414;
t395 = t371 * t398;
t394 = t373 * t397;
t354 = -qJD(2) * t394 + t371 * t407;
t358 = t371 * t376 - t394;
t393 = t358 * t403 - t354;
t360 = t371 * t397 + t373 * t376;
t356 = t360 * qJD(2);
t392 = t360 * t403 - t356;
t359 = t371 * t378 + t373 * t398;
t387 = -t359 * t369 + t368 * t410;
t391 = t358 * t375 - t377 * t387;
t361 = t373 * t378 - t395;
t349 = t361 * t369 + t368 * t411;
t390 = t349 * t377 + t360 * t375;
t389 = (qJD(2) - t403) * t378;
t388 = -t359 * t368 - t369 * t410;
t386 = -t361 * t368 + t369 * t411;
t385 = -t369 * pkin(4) - t418 * t368 - cos(pkin(11)) * pkin(3) - pkin(2);
t384 = qJD(4) * t420;
t383 = -t368 * t409 + t369 * t414;
t382 = t375 * t415 - t377 * t417 + pkin(4);
t355 = t359 * qJD(2);
t381 = qJD(5) * t359 - t355 * t369 + t358 * t405;
t357 = -qJD(2) * t395 + t373 * t406;
t380 = qJD(5) * t361 - t357 * t369 + t360 * t405;
t379 = qJD(6) * t375 + (t375 * t417 + t377 * t415) * qJD(5);
t374 = -pkin(8) - qJ(3);
t345 = qJD(4) * t383 + t369 * t401;
t343 = qJD(4) * t386 - t356 * t369;
t341 = qJD(4) * t388 - t354 * t369;
t336 = -t421 * qJD(5) + t345 * t375 - t377 * t400;
t330 = qJD(5) * t390 + t343 * t375 - t357 * t377;
t328 = qJD(5) * t391 + t341 * t375 - t355 * t377;
t1 = [0 -(t360 * t412 + t361 * t377) * qJD(6) + t356 * t374 + t361 * qJD(3) - t417 * (t375 * t392 + t377 * t380) + t415 * (t375 * t380 - t377 * t392) - t360 * t384 + t385 * t357, t357, t418 * t343 + t379 * t386 + t382 * (-qJD(4) * t349 + t356 * t368) t390 * qJD(6) + t415 * (t343 * t377 + t357 * t375 + (-t349 * t375 + t360 * t377) * qJD(5)) + t417 * t330, t330; 0 -(t358 * t412 + t359 * t377) * qJD(6) + t354 * t374 + t359 * qJD(3) - t417 * (t375 * t393 + t377 * t381) + t415 * (t375 * t381 - t377 * t393) - t358 * t384 + t385 * t355, t355, t418 * t341 + t379 * t388 + t382 * (qJD(4) * t387 + t354 * t368) t391 * qJD(6) + t415 * (t341 * t377 + t355 * t375 + (t358 * t377 + t375 * t387) * qJD(5)) + t417 * t328, t328; 0 (-t417 * (t375 * t389 - t419 * t377) - t415 * (t419 * t375 + t377 * t389) - (-t369 * t408 + t376 * t377) * qJD(6) + t376 * qJD(3) + t420 * t404 + (-t378 * t374 + t376 * t385) * qJD(2)) * t372, t400, t418 * t345 + t379 * t383 + t382 * (-qJD(4) * t353 - t368 * t401) -t421 * qJD(6) + t415 * (t375 * t400 + t345 * t377 + (-t372 * t377 * t378 - t353 * t375) * qJD(5)) + t417 * t336, t336;];
JaD_transl  = t1;
