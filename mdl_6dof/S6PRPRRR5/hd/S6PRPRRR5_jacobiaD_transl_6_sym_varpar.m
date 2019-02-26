% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:17
% EndTime: 2019-02-26 19:56:17
% DurationCPUTime: 0.55s
% Computational Cost: add. (607->80), mult. (1090->145), div. (0->0), fcn. (1076->12), ass. (0->60)
t358 = sin(qJ(6));
t361 = cos(qJ(6));
t406 = r_i_i_C(1) * t358 + r_i_i_C(2) * t361;
t403 = t406 * qJD(6);
t399 = pkin(10) + r_i_i_C(3);
t376 = r_i_i_C(1) * t361 - r_i_i_C(2) * t358;
t404 = pkin(5) + t376;
t357 = cos(pkin(6));
t356 = cos(pkin(11));
t363 = cos(qJ(2));
t386 = qJD(2) * t363;
t379 = t356 * t386;
t354 = sin(pkin(11));
t360 = sin(qJ(2));
t387 = qJD(2) * t360;
t381 = t354 * t387;
t341 = -t357 * t381 + t379;
t352 = qJD(4) + qJD(5);
t355 = sin(pkin(6));
t394 = t354 * t355;
t401 = -t352 * t394 + t341;
t371 = t376 * qJD(6);
t353 = qJ(4) + qJ(5);
t350 = sin(t353);
t351 = cos(t353);
t359 = sin(qJ(4));
t400 = t359 * pkin(4) + t350 * t404 - t399 * t351 + qJ(3);
t396 = t350 * t352;
t395 = t351 * t352;
t393 = t355 * t356;
t392 = t355 * t360;
t362 = cos(qJ(4));
t391 = t355 * t362;
t390 = t355 * t363;
t389 = t357 * t360;
t388 = t357 * t363;
t383 = t350 * t390;
t380 = t355 * t386;
t378 = t355 * t387;
t343 = t354 * t363 + t356 * t389;
t339 = t343 * qJD(2);
t374 = t352 * t393 + t339;
t372 = t350 * t357 + t351 * t390;
t344 = t354 * t388 + t356 * t360;
t370 = -pkin(2) - pkin(9) - pkin(8) - t406;
t321 = -t344 * t395 - t401 * t350;
t368 = -t403 * (t344 * t351 - t350 * t394) + t404 * (-t344 * t396 + t401 * t351) - t399 * t321;
t342 = t354 * t360 - t356 * t388;
t323 = t342 * t395 + t374 * t350;
t367 = -t403 * (t342 * t351 + t350 * t393) + t404 * (-t342 * t396 + t374 * t351) + t399 * t323;
t328 = -t350 * t378 + t372 * t352;
t366 = t403 * t372 + t404 * (t352 * t383 + (-t352 * t357 + t378) * t351) - t399 * t328;
t365 = qJD(4) * t362 * pkin(4) + qJD(3) - t350 * t403 + (t399 * t350 + t351 * t404) * t352;
t345 = -t354 * t389 + t356 * t363;
t340 = t344 * qJD(2);
t338 = -t357 * t379 + t381;
t337 = t351 * t357 - t383;
t333 = t342 * t350 - t351 * t393;
t331 = t344 * t350 + t351 * t394;
t1 = [0, -t340 * t400 + t370 * t341 - t344 * t371 + t365 * t345, t341 (t341 * t362 + (-t344 * t359 - t354 * t391) * qJD(4)) * pkin(4) + t368, t368 (t321 * t358 - t340 * t361) * r_i_i_C(1) + (t321 * t361 + t340 * t358) * r_i_i_C(2) + ((-t331 * t361 - t345 * t358) * r_i_i_C(1) + (t331 * t358 - t345 * t361) * r_i_i_C(2)) * qJD(6); 0, -t338 * t400 + t370 * t339 - t342 * t371 + t365 * t343, t339 (t339 * t362 + (-t342 * t359 + t356 * t391) * qJD(4)) * pkin(4) + t367, t367 (-t323 * t358 - t338 * t361) * r_i_i_C(1) + (-t323 * t361 + t338 * t358) * r_i_i_C(2) + ((-t333 * t361 - t343 * t358) * r_i_i_C(1) + (t333 * t358 - t343 * t361) * r_i_i_C(2)) * qJD(6); 0 ((t400 * qJD(2) + t371) * t363 + (t370 * qJD(2) + t365) * t360) * t355, t378 (t362 * t378 + (-t357 * t362 + t359 * t390) * qJD(4)) * pkin(4) + t366, t366 (t328 * t358 + t361 * t380) * r_i_i_C(1) + (t328 * t361 - t358 * t380) * r_i_i_C(2) + ((-t337 * t361 - t358 * t392) * r_i_i_C(1) + (t337 * t358 - t361 * t392) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
