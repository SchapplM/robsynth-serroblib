% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:47
% EndTime: 2019-02-26 20:08:47
% DurationCPUTime: 0.52s
% Computational Cost: add. (691->114), mult. (1620->184), div. (0->0), fcn. (1673->12), ass. (0->70)
t379 = cos(qJ(4));
t369 = t379 * pkin(4) + pkin(3);
t377 = sin(qJ(3));
t380 = cos(qJ(3));
t372 = qJ(4) + pkin(11);
t370 = sin(t372);
t411 = t370 * qJD(6);
t424 = r_i_i_C(2) + qJ(5) + pkin(9);
t376 = sin(qJ(4));
t425 = t376 * pkin(4);
t382 = -(t424 * qJD(3) - qJD(4) * t425 + t411) * t380 + (qJD(3) * t369 - qJD(5)) * t377;
t373 = sin(pkin(10));
t378 = sin(qJ(2));
t381 = cos(qJ(2));
t421 = cos(pkin(10));
t422 = cos(pkin(6));
t396 = t422 * t421;
t358 = t373 * t381 + t378 * t396;
t374 = sin(pkin(6));
t404 = t374 * t421;
t346 = t358 * t380 - t377 * t404;
t405 = t373 * t422;
t360 = -t378 * t405 + t421 * t381;
t417 = t374 * t380;
t362 = t422 * t377 + t378 * t417;
t371 = cos(t372);
t416 = t374 * t381;
t429 = -t362 * t371 + t370 * t416;
t414 = qJD(3) * t377;
t428 = (qJD(2) * t380 - qJD(4)) * t378 + t381 * t414;
t426 = -r_i_i_C(1) - pkin(5);
t423 = r_i_i_C(3) + qJ(6);
t418 = t374 * t377;
t415 = qJD(2) * t378;
t413 = qJD(4) * t379;
t412 = qJD(4) * t380;
t410 = t371 * qJD(6);
t408 = t374 * t415;
t407 = qJD(2) * t416;
t392 = t381 * t396;
t353 = -qJD(2) * t392 + t373 * t415;
t357 = t373 * t378 - t392;
t398 = t357 * t412 - t353;
t359 = t421 * t378 + t381 * t405;
t355 = t359 * qJD(2);
t397 = t359 * t412 - t355;
t395 = t346 * t371 + t357 * t370;
t348 = t360 * t380 + t373 * t418;
t394 = t348 * t371 + t359 * t370;
t393 = (qJD(2) - t412) * t381;
t390 = -t360 * t377 + t373 * t417;
t389 = -t380 * t369 - t424 * t377 - pkin(2);
t388 = -t358 * t377 - t380 * t404;
t387 = -t378 * t418 + t422 * t380;
t386 = -t423 * t370 + t426 * t371 - t369;
t354 = t358 * qJD(2);
t385 = qJD(4) * t358 - t354 * t380 + t357 * t414;
t356 = t360 * qJD(2);
t384 = qJD(4) * t360 - t356 * t380 + t359 * t414;
t383 = t411 + (t426 * t370 + t423 * t371 - t425) * qJD(4);
t350 = t387 * qJD(3) + t380 * t407;
t349 = t362 * qJD(3) + t377 * t407;
t344 = t390 * qJD(3) - t355 * t380;
t343 = t348 * qJD(3) - t355 * t377;
t342 = t388 * qJD(3) - t353 * t380;
t341 = t346 * qJD(3) - t353 * t377;
t337 = -t429 * qJD(4) + t350 * t370 - t371 * t408;
t331 = t394 * qJD(4) + t344 * t370 - t356 * t371;
t329 = t395 * qJD(4) + t342 * t370 - t354 * t371;
t1 = [0, -t360 * t410 - t355 * pkin(8) - t426 * (t397 * t370 + t384 * t371) + t423 * (t384 * t370 - t397 * t371) + (-t355 * t376 + t360 * t413) * pkin(4) + t389 * t356 + t382 * t359, t348 * qJD(5) + t386 * t343 + t424 * t344 + t383 * t390, t394 * qJD(6) + t423 * (t344 * t371 + t356 * t370 + (-t348 * t370 + t359 * t371) * qJD(4)) + t426 * t331 + (-t344 * t376 + t356 * t379 + (-t348 * t379 - t359 * t376) * qJD(4)) * pkin(4), t343, t331; 0, -t358 * t410 - t353 * pkin(8) - t426 * (t398 * t370 + t385 * t371) + t423 * (t385 * t370 - t398 * t371) + (-t353 * t376 + t358 * t413) * pkin(4) + t389 * t354 + t382 * t357, t346 * qJD(5) + t386 * t341 + t424 * t342 + t383 * t388, t395 * qJD(6) + t423 * (t342 * t371 + t354 * t370 + (-t346 * t370 + t357 * t371) * qJD(4)) + t426 * t329 + (-t342 * t376 + t354 * t379 + (-t346 * t379 - t357 * t376) * qJD(4)) * pkin(4), t341, t329; 0 (-t426 * (t370 * t393 - t428 * t371) - t423 * (t428 * t370 + t371 * t393) + (pkin(4) * t413 + t389 * qJD(2) - t410) * t378 + ((pkin(8) + t425) * qJD(2) - t382) * t381) * t374, t362 * qJD(5) + t386 * t349 + t424 * t350 + t383 * t387, -t429 * qJD(6) + t423 * (t370 * t408 + t350 * t371 + (-t362 * t370 - t371 * t416) * qJD(4)) + t426 * t337 + (t379 * t408 - t350 * t376 + (-t362 * t379 + t376 * t416) * qJD(4)) * pkin(4), t349, t337;];
JaD_transl  = t1;
