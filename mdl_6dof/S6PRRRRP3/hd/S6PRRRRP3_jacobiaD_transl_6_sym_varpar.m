% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:22
% EndTime: 2019-02-26 20:16:22
% DurationCPUTime: 0.41s
% Computational Cost: add. (649->83), mult. (1308->134), div. (0->0), fcn. (1304->12), ass. (0->71)
t349 = sin(qJ(3));
t352 = cos(qJ(3));
t345 = qJ(4) + qJ(5);
t341 = cos(t345);
t351 = cos(qJ(4));
t334 = t351 * pkin(4) + pkin(5) * t341;
t340 = sin(t345);
t371 = t341 * r_i_i_C(1) - t340 * r_i_i_C(2);
t367 = pkin(3) + t334 + t371;
t394 = r_i_i_C(3) + qJ(6) + pkin(10) + pkin(9);
t348 = sin(qJ(4));
t393 = pkin(4) * qJD(4);
t344 = qJD(4) + qJD(5);
t395 = pkin(5) * t344;
t330 = -t340 * t395 - t348 * t393;
t370 = t340 * r_i_i_C(1) + t341 * r_i_i_C(2);
t397 = t370 * t344 - t330;
t354 = t397 * t352 + (t367 * t349 - t394 * t352) * qJD(3) - t349 * qJD(6);
t346 = sin(pkin(11));
t350 = sin(qJ(2));
t353 = cos(qJ(2));
t391 = cos(pkin(11));
t392 = cos(pkin(6));
t368 = t392 * t391;
t325 = t346 * t353 + t350 * t368;
t347 = sin(pkin(6));
t379 = t347 * t391;
t315 = t325 * t352 - t349 * t379;
t380 = t346 * t392;
t327 = -t350 * t380 + t391 * t353;
t389 = t347 * t349;
t388 = t347 * t352;
t387 = t347 * t353;
t321 = t325 * qJD(2);
t375 = t315 * t344 - t321;
t366 = t353 * t368;
t383 = qJD(2) * t350;
t320 = -qJD(2) * t366 + t346 * t383;
t362 = -t325 * t349 - t352 * t379;
t311 = t362 * qJD(3) - t320 * t352;
t324 = t346 * t350 - t366;
t377 = -t324 * t344 - t311;
t358 = t377 * t340 - t375 * t341;
t386 = t358 * r_i_i_C(1) + (t375 * t340 + t377 * t341) * r_i_i_C(2);
t317 = t327 * t352 + t346 * t389;
t323 = t327 * qJD(2);
t374 = t317 * t344 - t323;
t326 = t391 * t350 + t353 * t380;
t322 = t326 * qJD(2);
t365 = -t327 * t349 + t346 * t388;
t313 = t365 * qJD(3) - t322 * t352;
t376 = -t326 * t344 - t313;
t357 = t376 * t340 - t374 * t341;
t385 = t357 * r_i_i_C(1) + (t374 * t340 + t376 * t341) * r_i_i_C(2);
t329 = t392 * t349 + t350 * t388;
t364 = -t329 * t344 + t347 * t383;
t361 = -t350 * t389 + t392 * t352;
t381 = qJD(2) * t387;
t319 = t361 * qJD(3) + t352 * t381;
t369 = t344 * t387 - t319;
t355 = t369 * t340 + t364 * t341;
t384 = t355 * r_i_i_C(1) + (-t364 * t340 + t369 * t341) * r_i_i_C(2);
t333 = t348 * pkin(4) + pkin(5) * t340;
t363 = pkin(8) + t333 + t370;
t331 = t341 * t395 + t351 * t393;
t359 = t371 * t344 + t331;
t356 = -t394 * t349 - t367 * t352 - pkin(2);
t318 = t329 * qJD(3) + t349 * t381;
t312 = t317 * qJD(3) - t322 * t349;
t310 = t315 * qJD(3) - t320 * t349;
t1 = [0, -t363 * t322 + t356 * t323 + t354 * t326 + t359 * t327, t317 * qJD(6) - t367 * t312 + t394 * t313 - t365 * t397, -t313 * t333 - t317 * t331 + t323 * t334 + t326 * t330 + t385, t357 * pkin(5) + t385, t312; 0, -t363 * t320 + t356 * t321 + t354 * t324 + t359 * t325, t315 * qJD(6) - t367 * t310 + t394 * t311 - t362 * t397, -t311 * t333 - t315 * t331 + t321 * t334 + t324 * t330 + t386, t358 * pkin(5) + t386, t310; 0 ((t356 * qJD(2) + t359) * t350 + (t363 * qJD(2) - t354) * t353) * t347, t329 * qJD(6) - t367 * t318 + t394 * t319 - t361 * t397, -t319 * t333 - t329 * t331 + (-t353 * t330 + t334 * t383) * t347 + t384, t355 * pkin(5) + t384, t318;];
JaD_transl  = t1;
