% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:26
% EndTime: 2019-02-26 21:31:27
% DurationCPUTime: 0.46s
% Computational Cost: add. (661->80), mult. (1098->127), div. (0->0), fcn. (1018->10), ass. (0->54)
t389 = pkin(10) + qJ(5);
t353 = sin(t389);
t360 = cos(qJ(2));
t357 = sin(qJ(2));
t383 = cos(t389);
t381 = t357 * t383;
t409 = -t360 * t353 + t381;
t384 = pkin(4) * sin(pkin(10)) + qJ(3);
t401 = pkin(2) + cos(pkin(10)) * pkin(4) + pkin(3);
t404 = t401 * t357 - t384 * t360;
t415 = t404 * qJD(2) - t357 * qJD(3);
t379 = qJD(5) * t383;
t380 = qJD(2) * t383;
t392 = qJD(5) * t353;
t414 = t360 * t392 + (-t379 + t380) * t357;
t333 = t357 * t353 + t360 * t383;
t395 = qJD(2) * t357;
t327 = t333 * qJD(5) - t353 * t395 - t360 * t380;
t358 = sin(qJ(1));
t361 = cos(qJ(1));
t410 = t409 * qJD(1);
t325 = t327 * t358 - t410 * t361;
t394 = qJD(2) * t360;
t328 = t353 * t394 - t414;
t367 = qJD(1) * t333;
t326 = t328 * t358 + t361 * t367;
t356 = sin(qJ(6));
t359 = cos(qJ(6));
t376 = r_i_i_C(1) * t359 - r_i_i_C(2) * t356 + pkin(5);
t403 = -r_i_i_C(3) - pkin(9);
t402 = r_i_i_C(2) * t359;
t408 = qJD(6) * (r_i_i_C(1) * t356 + t402);
t413 = -t409 * t358 * t408 - t376 * t325 - t403 * t326;
t393 = qJD(2) * t361;
t386 = t360 * t393;
t323 = -t353 * t386 + t358 * t367 + t414 * t361;
t332 = t333 * t361;
t398 = t360 * t361;
t324 = -t361 * t357 * t392 + qJD(2) * t332 - t410 * t358 - t379 * t398;
t412 = (t353 * t398 - t361 * t381) * t408 + t403 * t323 + t376 * t324;
t411 = t403 * t327 - t376 * t328 + t333 * t408;
t400 = pkin(7) - pkin(8) - qJ(4);
t397 = qJD(1) * t358;
t396 = qJD(1) * t361;
t391 = qJD(6) * t409;
t330 = t333 * t358;
t378 = t330 * t359 + t356 * t361;
t377 = t330 * t356 - t359 * t361;
t369 = -t384 * t357 - t401 * t360;
t365 = -pkin(1) + t369;
t351 = t356 * t397;
t322 = -t356 * t396 - t323 * t359 + (-t332 * t356 - t358 * t359) * qJD(6);
t321 = -t359 * t396 + t323 * t356 + (-t332 * t359 + t356 * t358) * qJD(6);
t1 = [t351 * r_i_i_C(1) - t361 * qJD(4) - t376 * t326 + t403 * t325 + (t377 * r_i_i_C(1) + t378 * r_i_i_C(2)) * qJD(6) + t415 * t358 + ((-t400 + t402) * t358 + t365 * t361) * qJD(1) (-t384 * t393 + t401 * t397) * t357 + (-t384 * t397 + (-t401 * qJD(2) + qJD(3)) * t361) * t360 - t412, -t357 * t397 + t386, -t396, t412, r_i_i_C(1) * t321 - r_i_i_C(2) * t322; -t323 * pkin(5) + t322 * r_i_i_C(1) + t321 * r_i_i_C(2) - t358 * qJD(4) + t403 * t324 - t415 * t361 + (t365 * t358 + t400 * t361) * qJD(1), -t404 * t396 + (t369 * qJD(2) + qJD(3) * t360) * t358 - t413, t357 * t396 + t358 * t394, -t397, t413 (-t326 * t356 - t359 * t397) * r_i_i_C(1) + (-t326 * t359 + t351) * r_i_i_C(2) + (-t378 * r_i_i_C(1) + t377 * r_i_i_C(2)) * qJD(6); 0, -t415 - t411, t395, 0, t411 (t327 * t359 + t356 * t391) * r_i_i_C(2) + (t327 * t356 - t359 * t391) * r_i_i_C(1);];
JaD_transl  = t1;
