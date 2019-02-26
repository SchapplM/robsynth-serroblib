% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP9_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP9_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:08
% EndTime: 2019-02-26 22:30:08
% DurationCPUTime: 0.53s
% Computational Cost: add. (444->96), mult. (1331->167), div. (0->0), fcn. (1334->10), ass. (0->65)
t374 = sin(qJ(1));
t370 = cos(pkin(6));
t390 = qJD(2) * t370 + qJD(1);
t373 = sin(qJ(2));
t408 = t373 * t374;
t397 = t370 * t408;
t403 = qJD(2) * t373;
t377 = cos(qJ(2));
t378 = cos(qJ(1));
t405 = t377 * t378;
t348 = -qJD(1) * t397 - t374 * t403 + t390 * t405;
t406 = t374 * t377;
t407 = t373 * t378;
t359 = t370 * t407 + t406;
t372 = sin(qJ(3));
t376 = cos(qJ(3));
t369 = sin(pkin(6));
t404 = qJD(1) * t369;
t395 = t374 * t404;
t409 = t369 * t378;
t398 = t376 * t409;
t342 = (-qJD(3) * t359 + t395) * t372 - qJD(3) * t398 + t348 * t376;
t360 = t370 * t406 + t407;
t347 = t360 * qJD(1) + t359 * qJD(2);
t371 = sin(qJ(4));
t375 = cos(qJ(4));
t426 = t342 * t371 - t347 * t375;
t425 = -t342 * t375 - t347 * t371;
t388 = t375 * r_i_i_C(1) - t371 * r_i_i_C(2);
t386 = pkin(3) + t388;
t419 = r_i_i_C(3) + pkin(10);
t424 = (t386 * t372 - t419 * t376) * qJD(3);
t353 = -t359 * t376 + t372 * t409;
t396 = t370 * t405;
t358 = -t396 + t408;
t423 = -t353 * t371 - t358 * t375;
t422 = t353 * t375 - t358 * t371;
t387 = t371 * r_i_i_C(1) + t375 * r_i_i_C(2);
t420 = t419 * t372 + t386 * t376 + pkin(2);
t412 = t369 * t374;
t411 = t369 * t376;
t410 = t369 * t377;
t402 = qJD(2) * t377;
t401 = qJD(4) * t371;
t400 = qJD(4) * t375;
t399 = qJD(4) * t376;
t394 = t378 * t404;
t393 = t369 * t403;
t392 = t369 * t402;
t361 = -t397 + t405;
t385 = -t361 * t372 + t374 * t411;
t355 = t361 * t376 + t372 * t412;
t357 = t370 * t372 + t373 * t411;
t384 = -t369 * t372 * t373 + t370 * t376;
t383 = qJD(4) * t387;
t380 = t353 * qJD(3) - t348 * t372 + t376 * t395;
t379 = t387 * t399 + t424;
t350 = t384 * qJD(3) + t376 * t392;
t346 = t359 * qJD(1) + t360 * qJD(2);
t345 = -qJD(1) * t396 - t378 * t402 + t390 * t408;
t340 = t385 * qJD(3) - t346 * t376 + t372 * t394;
t339 = t355 * qJD(3) - t346 * t372 - t376 * t394;
t338 = t340 * t375 - t345 * t371 + (-t355 * t371 + t360 * t375) * qJD(4);
t337 = -t340 * t371 - t345 * t375 + (-t355 * t375 - t360 * t371) * qJD(4);
t1 = [t425 * r_i_i_C(1) + t426 * r_i_i_C(2) - t342 * pkin(3) - t348 * pkin(2) - t347 * pkin(9) + t419 * t380 + (t423 * r_i_i_C(1) - t422 * r_i_i_C(2)) * qJD(4) + (-pkin(1) * t378 - pkin(8) * t412) * qJD(1) (-t346 * t371 + t361 * t400) * r_i_i_C(1) + (-t346 * t375 - t361 * t401) * r_i_i_C(2) - t346 * pkin(9) + t420 * t345 + t379 * t360, -t386 * t339 + t419 * t340 - t385 * t383, r_i_i_C(1) * t337 - t338 * r_i_i_C(2), 0, 0; -t346 * pkin(2) + t340 * pkin(3) - t345 * pkin(9) + t338 * r_i_i_C(1) + t337 * r_i_i_C(2) + t419 * t339 + (-pkin(1) * t374 + pkin(8) * t409) * qJD(1) (t348 * t371 + t359 * t400) * r_i_i_C(1) + (t348 * t375 - t359 * t401) * r_i_i_C(2) + t348 * pkin(9) - t420 * t347 + t379 * t358, t419 * t342 - (-t359 * t372 - t398) * t383 + t386 * t380, -t426 * r_i_i_C(1) + t425 * r_i_i_C(2) + (t422 * r_i_i_C(1) + t423 * r_i_i_C(2)) * qJD(4), 0, 0; 0 ((-qJD(2) * t420 + t388 * qJD(4)) * t373 + (qJD(2) * pkin(9) - t424 + t387 * (qJD(2) - t399)) * t377) * t369, t419 * t350 - t384 * t383 + t386 * (-t357 * qJD(3) - t372 * t392) (-t350 * t371 + t375 * t393) * r_i_i_C(1) + (-t350 * t375 - t371 * t393) * r_i_i_C(2) + ((-t357 * t375 + t371 * t410) * r_i_i_C(1) + (t357 * t371 + t375 * t410) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
