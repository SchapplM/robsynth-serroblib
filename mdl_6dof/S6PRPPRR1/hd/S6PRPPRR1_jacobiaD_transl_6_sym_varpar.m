% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:37
% EndTime: 2019-02-26 19:44:37
% DurationCPUTime: 0.49s
% Computational Cost: add. (575->89), mult. (1408->161), div. (0->0), fcn. (1519->13), ass. (0->58)
t386 = pkin(12) + qJ(5);
t384 = sin(t386);
t385 = cos(t386);
t393 = sin(qJ(6));
t395 = cos(qJ(6));
t402 = (t393 * r_i_i_C(1) + t395 * r_i_i_C(2)) * qJD(6);
t407 = t395 * r_i_i_C(1) - t393 * r_i_i_C(2) + pkin(5);
t424 = pkin(9) + r_i_i_C(3);
t428 = (t407 * t384 - t424 * t385) * qJD(5) + t385 * t402;
t391 = cos(pkin(6));
t387 = sin(pkin(11));
t394 = sin(qJ(2));
t421 = cos(pkin(11));
t423 = cos(qJ(2));
t401 = t423 * t387 + t394 * t421;
t372 = t401 * t391;
t410 = t423 * t421;
t416 = qJD(2) * t394;
t426 = -qJD(2) * t410 + t387 * t416;
t400 = -t394 * t387 + t410;
t397 = t424 * t384 + t407 * t385 + cos(pkin(12)) * pkin(4) + pkin(3);
t422 = pkin(2) * qJD(2);
t388 = sin(pkin(10));
t389 = sin(pkin(6));
t420 = t388 * t389;
t390 = cos(pkin(10));
t419 = t389 * t390;
t418 = t391 * t394;
t415 = qJD(6) * t393;
t414 = qJD(6) * t395;
t369 = t426 * t391;
t374 = t401 * qJD(2);
t351 = t390 * t369 + t388 * t374;
t353 = t388 * t369 - t390 * t374;
t371 = t401 * t389;
t362 = t371 * t385 + t391 * t384;
t408 = -t371 * t384 + t391 * t385;
t357 = t390 * t372 + t388 * t400;
t358 = t388 * t372 - t390 * t400;
t405 = -t357 * t384 - t385 * t419;
t404 = -t357 * t385 + t384 * t419;
t403 = t358 * t384 + t385 * t420;
t348 = -t358 * t385 + t384 * t420;
t399 = t400 * t391;
t398 = qJD(2) * t372;
t392 = -pkin(8) - qJ(4);
t373 = t400 * qJD(2);
t370 = t400 * t389;
t368 = qJD(2) * t371;
t367 = t426 * t389;
t359 = -t388 * t399 - t390 * t401;
t356 = -t388 * t401 + t390 * t399;
t352 = -t390 * t373 + t388 * t398;
t349 = -t388 * t373 - t390 * t398;
t344 = t408 * qJD(5) - t367 * t385;
t342 = t403 * qJD(5) + t353 * t385;
t340 = t405 * qJD(5) - t351 * t385;
t1 = [0 (t353 * t393 - t358 * t414) * r_i_i_C(1) + (t353 * t395 + t358 * t415) * r_i_i_C(2) - t353 * t392 - t358 * qJD(4) + (t388 * t418 - t423 * t390) * t422 + t397 * t352 - t428 * t359, 0, -t352, t424 * t342 - t403 * t402 + t407 * (-t348 * qJD(5) - t353 * t384) (-t342 * t393 - t352 * t395) * r_i_i_C(1) + (-t342 * t395 + t352 * t393) * r_i_i_C(2) + ((-t348 * t395 + t359 * t393) * r_i_i_C(1) + (t348 * t393 + t359 * t395) * r_i_i_C(2)) * qJD(6); 0 (-t351 * t393 + t357 * t414) * r_i_i_C(1) + (-t351 * t395 - t357 * t415) * r_i_i_C(2) + t351 * t392 + t357 * qJD(4) + (-t423 * t388 - t390 * t418) * t422 + t397 * t349 - t428 * t356, 0, -t349, t424 * t340 - t405 * t402 + t407 * (t404 * qJD(5) + t351 * t384) (-t340 * t393 - t349 * t395) * r_i_i_C(1) + (-t340 * t395 + t349 * t393) * r_i_i_C(2) + ((t356 * t393 + t395 * t404) * r_i_i_C(1) + (t356 * t395 - t393 * t404) * r_i_i_C(2)) * qJD(6); 0 (-t367 * t393 + t371 * t414) * r_i_i_C(1) + (-t367 * t395 - t371 * t415) * r_i_i_C(2) + t367 * t392 + t371 * qJD(4) - t389 * pkin(2) * t416 - t397 * t368 - t428 * t370, 0, t368, t424 * t344 - t408 * t402 + t407 * (-t362 * qJD(5) + t367 * t384) (-t344 * t393 + t368 * t395) * r_i_i_C(1) + (-t344 * t395 - t368 * t393) * r_i_i_C(2) + ((-t362 * t395 + t370 * t393) * r_i_i_C(1) + (t362 * t393 + t370 * t395) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
