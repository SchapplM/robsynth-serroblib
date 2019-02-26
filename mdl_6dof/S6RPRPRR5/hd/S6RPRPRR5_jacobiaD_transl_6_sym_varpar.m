% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:17
% EndTime: 2019-02-26 20:51:17
% DurationCPUTime: 0.43s
% Computational Cost: add. (680->82), mult. (1056->133), div. (0->0), fcn. (985->9), ass. (0->56)
t354 = pkin(10) + qJ(3);
t353 = cos(t354);
t357 = sin(qJ(5));
t352 = sin(t354);
t397 = cos(qJ(5));
t383 = t352 * t397;
t405 = -t353 * t357 + t383;
t358 = sin(qJ(1));
t360 = cos(qJ(1));
t378 = qJD(3) * t397;
t375 = t358 * t378;
t406 = t405 * qJD(1);
t332 = t352 * t357 + t353 * t397;
t389 = qJD(3) * t357;
t407 = t332 * qJD(5) - t352 * t389;
t324 = -t353 * t375 + t407 * t358 - t406 * t360;
t377 = qJD(5) * t397;
t374 = t352 * t377;
t362 = t353 * t389 + t374;
t365 = qJD(1) * t332;
t386 = qJD(5) * t357;
t379 = t353 * t386;
t325 = -t352 * t375 + t360 * t365 + (t362 - t379) * t358;
t356 = sin(qJ(6));
t359 = cos(qJ(6));
t369 = r_i_i_C(1) * t359 - r_i_i_C(2) * t356 + pkin(5);
t398 = -r_i_i_C(3) - pkin(9);
t396 = r_i_i_C(2) * t359;
t404 = qJD(6) * (r_i_i_C(1) * t356 + t396);
t411 = -t405 * t358 * t404 - t369 * t324 - t398 * t325;
t363 = t352 * t378 + t379;
t387 = qJD(3) * t360;
t380 = t353 * t387;
t322 = -t357 * t380 + t358 * t365 + (t363 - t374) * t360;
t331 = t332 * t360;
t392 = t353 * t360;
t323 = -t360 * t352 * t386 + qJD(3) * t331 - t406 * t358 - t377 * t392;
t410 = (t357 * t392 - t360 * t383) * t404 + t398 * t322 + t369 * t323;
t326 = -t353 * t378 + t407;
t409 = t332 * t404 + t398 * t326 - t369 * (t362 - t363);
t399 = pkin(3) + pkin(4);
t408 = (-qJ(4) * t353 + t399 * t352) * qJD(3) - t352 * qJD(4);
t395 = -pkin(8) + pkin(7) + qJ(2);
t391 = qJD(1) * t358;
t390 = qJD(1) * t360;
t388 = qJD(3) * t358;
t385 = qJD(6) * t405;
t372 = -t399 * qJD(3) + qJD(4);
t329 = t332 * t358;
t371 = t329 * t359 + t356 * t360;
t370 = t329 * t356 - t359 * t360;
t366 = -qJ(4) * t352 - t399 * t353 - cos(pkin(10)) * pkin(2) - pkin(1);
t350 = t356 * t391;
t321 = -t356 * t390 - t322 * t359 + (-t331 * t356 - t358 * t359) * qJD(6);
t320 = -t359 * t390 + t322 * t356 + (-t331 * t359 + t356 * t358) * qJD(6);
t1 = [t350 * r_i_i_C(1) + t360 * qJD(2) - t369 * t325 + t398 * t324 + (t370 * r_i_i_C(1) + t371 * r_i_i_C(2)) * qJD(6) + t408 * t358 + ((-t395 + t396) * t358 + t366 * t360) * qJD(1), t390 (-qJ(4) * t387 + t399 * t391) * t352 + (-qJ(4) * t391 + t372 * t360) * t353 - t410, -t352 * t391 + t380, t410, r_i_i_C(1) * t320 - r_i_i_C(2) * t321; -t322 * pkin(5) + t321 * r_i_i_C(1) + t320 * r_i_i_C(2) + t358 * qJD(2) + t398 * t323 - t408 * t360 + (t366 * t358 + t395 * t360) * qJD(1), t391 (-qJ(4) * t388 - t399 * t390) * t352 + (qJ(4) * t390 + t372 * t358) * t353 - t411, t352 * t390 + t353 * t388, t411 (-t325 * t356 - t359 * t391) * r_i_i_C(1) + (-t325 * t359 + t350) * r_i_i_C(2) + (-t371 * r_i_i_C(1) + t370 * r_i_i_C(2)) * qJD(6); 0, 0, -t408 - t409, qJD(3) * t352, t409 (t326 * t359 + t356 * t385) * r_i_i_C(2) + (t326 * t356 - t359 * t385) * r_i_i_C(1);];
JaD_transl  = t1;
