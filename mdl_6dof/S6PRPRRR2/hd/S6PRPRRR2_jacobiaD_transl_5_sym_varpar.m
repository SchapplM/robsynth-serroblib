% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR2
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

function JaD_transl = S6PRPRRR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:13
% EndTime: 2019-02-26 19:54:14
% DurationCPUTime: 0.51s
% Computational Cost: add. (429->82), mult. (1352->159), div. (0->0), fcn. (1456->12), ass. (0->57)
t379 = sin(qJ(4));
t382 = cos(qJ(4));
t378 = sin(qJ(5));
t381 = cos(qJ(5));
t389 = (t378 * r_i_i_C(1) + t381 * r_i_i_C(2)) * qJD(5);
t393 = r_i_i_C(1) * t381 - r_i_i_C(2) * t378 + pkin(4);
t412 = pkin(9) + r_i_i_C(3);
t416 = (t393 * t379 - t412 * t382) * qJD(4) + t382 * t389;
t380 = sin(qJ(2));
t373 = sin(pkin(12));
t383 = cos(qJ(2));
t405 = t383 * t373;
t410 = cos(pkin(12));
t366 = -t380 * t410 - t405;
t377 = cos(pkin(6));
t362 = t366 * t377;
t399 = qJD(2) * t410;
t404 = qJD(2) * t380;
t414 = t373 * t404 - t383 * t399;
t388 = -t380 * t373 + t383 * t410;
t385 = t412 * t379 + t393 * t382 + pkin(3);
t411 = pkin(2) * qJD(2);
t375 = sin(pkin(6));
t409 = t375 * t379;
t408 = t375 * t382;
t407 = t377 * t380;
t403 = qJD(5) * t378;
t402 = qJD(5) * t381;
t359 = t414 * t377;
t364 = -qJD(2) * t405 - t380 * t399;
t374 = sin(pkin(11));
t376 = cos(pkin(11));
t342 = t359 * t376 - t364 * t374;
t344 = t374 * t359 + t364 * t376;
t361 = t366 * t375;
t353 = -t361 * t382 + t377 * t379;
t396 = t361 * t379 + t377 * t382;
t395 = -t362 * t376 + t374 * t388;
t394 = t362 * t374 + t376 * t388;
t392 = -t376 * t408 - t379 * t395;
t391 = t376 * t409 - t382 * t395;
t390 = t374 * t408 - t379 * t394;
t339 = t374 * t409 + t382 * t394;
t387 = t377 * t388;
t386 = qJD(2) * t362;
t363 = t388 * qJD(2);
t360 = t388 * t375;
t358 = qJD(2) * t361;
t357 = t414 * t375;
t350 = t366 * t376 - t374 * t387;
t347 = t374 * t366 + t376 * t387;
t343 = -t363 * t376 - t374 * t386;
t340 = -t374 * t363 + t376 * t386;
t335 = t396 * qJD(4) - t357 * t382;
t333 = t390 * qJD(4) + t344 * t382;
t331 = t392 * qJD(4) - t342 * t382;
t1 = [0 (t344 * t378 + t394 * t402) * r_i_i_C(1) + (t344 * t381 - t394 * t403) * r_i_i_C(2) + t344 * pkin(8) + (t374 * t407 - t376 * t383) * t411 + t385 * t343 - t416 * t350, 0, t412 * t333 - t390 * t389 + t393 * (-t339 * qJD(4) - t344 * t379) (-t333 * t378 - t343 * t381) * r_i_i_C(1) + (-t333 * t381 + t343 * t378) * r_i_i_C(2) + ((-t339 * t381 + t350 * t378) * r_i_i_C(1) + (t339 * t378 + t350 * t381) * r_i_i_C(2)) * qJD(5), 0; 0 (-t342 * t378 + t395 * t402) * r_i_i_C(1) + (-t342 * t381 - t395 * t403) * r_i_i_C(2) - t342 * pkin(8) + (-t374 * t383 - t376 * t407) * t411 + t385 * t340 - t416 * t347, 0, t412 * t331 - t392 * t389 + t393 * (t391 * qJD(4) + t342 * t379) (-t331 * t378 - t340 * t381) * r_i_i_C(1) + (-t331 * t381 + t340 * t378) * r_i_i_C(2) + ((t347 * t378 + t381 * t391) * r_i_i_C(1) + (t347 * t381 - t378 * t391) * r_i_i_C(2)) * qJD(5), 0; 0 (-t357 * t378 - t361 * t402) * r_i_i_C(1) + (-t357 * t381 + t361 * t403) * r_i_i_C(2) - t357 * pkin(8) - t375 * pkin(2) * t404 + t385 * t358 - t416 * t360, 0, t412 * t335 - t396 * t389 + t393 * (-t353 * qJD(4) + t357 * t379) (-t335 * t378 - t358 * t381) * r_i_i_C(1) + (-t335 * t381 + t358 * t378) * r_i_i_C(2) + ((-t353 * t381 + t360 * t378) * r_i_i_C(1) + (t353 * t378 + t360 * t381) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
