% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:53
% EndTime: 2019-02-26 20:19:54
% DurationCPUTime: 0.45s
% Computational Cost: add. (889->90), mult. (1423->142), div. (0->0), fcn. (1422->14), ass. (0->70)
t350 = sin(qJ(3));
t353 = cos(qJ(3));
t346 = qJ(4) + qJ(5);
t341 = cos(t346);
t352 = cos(qJ(4));
t332 = t352 * pkin(4) + pkin(5) * t341;
t342 = qJ(6) + t346;
t337 = sin(t342);
t338 = cos(t342);
t370 = t338 * r_i_i_C(1) - t337 * r_i_i_C(2);
t366 = pkin(3) + t332 + t370;
t391 = r_i_i_C(3) + pkin(11) + pkin(10) + pkin(9);
t340 = sin(t346);
t349 = sin(qJ(4));
t390 = pkin(4) * qJD(4);
t344 = qJD(4) + qJD(5);
t392 = pkin(5) * t344;
t328 = -t340 * t392 - t349 * t390;
t339 = qJD(6) + t344;
t369 = t337 * r_i_i_C(1) + t338 * r_i_i_C(2);
t394 = t339 * t369 - t328;
t355 = t394 * t353 + (t350 * t366 - t353 * t391) * qJD(3);
t351 = sin(qJ(2));
t354 = cos(qJ(2));
t347 = sin(pkin(12));
t389 = cos(pkin(6));
t378 = t347 * t389;
t388 = cos(pkin(12));
t325 = -t351 * t378 + t388 * t354;
t348 = sin(pkin(6));
t387 = t348 * t350;
t386 = t348 * t353;
t385 = t348 * t354;
t367 = t389 * t388;
t323 = t347 * t354 + t351 * t367;
t319 = t323 * qJD(2);
t377 = t348 * t388;
t359 = -t323 * t353 + t350 * t377;
t373 = -t339 * t359 - t319;
t365 = t354 * t367;
t381 = qJD(2) * t351;
t318 = -qJD(2) * t365 + t347 * t381;
t361 = -t323 * t350 - t353 * t377;
t309 = qJD(3) * t361 - t318 * t353;
t322 = t347 * t351 - t365;
t375 = -t322 * t339 - t309;
t384 = (t337 * t375 - t338 * t373) * r_i_i_C(1) + (t337 * t373 + t338 * t375) * r_i_i_C(2);
t315 = t325 * t353 + t347 * t387;
t321 = t325 * qJD(2);
t372 = t315 * t339 - t321;
t324 = t351 * t388 + t354 * t378;
t320 = t324 * qJD(2);
t364 = -t325 * t350 + t347 * t386;
t311 = qJD(3) * t364 - t320 * t353;
t374 = -t324 * t339 - t311;
t383 = (t337 * t374 - t338 * t372) * r_i_i_C(1) + (t337 * t372 + t338 * t374) * r_i_i_C(2);
t327 = t350 * t389 + t351 * t386;
t380 = t348 * t381;
t363 = -t327 * t339 + t380;
t360 = -t351 * t387 + t353 * t389;
t379 = qJD(2) * t385;
t317 = qJD(3) * t360 + t353 * t379;
t368 = t339 * t385 - t317;
t382 = (t337 * t368 + t338 * t363) * r_i_i_C(1) + (-t337 * t363 + t338 * t368) * r_i_i_C(2);
t331 = t349 * pkin(4) + pkin(5) * t340;
t362 = pkin(8) + t331 + t369;
t329 = t341 * t392 + t352 * t390;
t357 = t339 * t370 + t329;
t356 = -t350 * t391 - t353 * t366 - pkin(2);
t1 = [0, -t320 * t362 + t321 * t356 + t324 * t355 + t325 * t357, t391 * t311 - t394 * t364 + t366 * (-qJD(3) * t315 + t320 * t350) -t311 * t331 - t315 * t329 + t321 * t332 + t324 * t328 + t383 ((-t315 * t344 + t321) * t341 + (-t324 * t344 - t311) * t340) * pkin(5) + t383, t383; 0, -t318 * t362 + t319 * t356 + t322 * t355 + t323 * t357, t391 * t309 - t394 * t361 + t366 * (qJD(3) * t359 + t318 * t350) -t309 * t331 + t319 * t332 + t322 * t328 + t329 * t359 + t384 ((t344 * t359 + t319) * t341 + (-t322 * t344 - t309) * t340) * pkin(5) + t384, t384; 0 ((qJD(2) * t356 + t357) * t351 + (t362 * qJD(2) - t355) * t354) * t348, t391 * t317 - t394 * t360 + t366 * (-qJD(3) * t327 - t350 * t379) -t317 * t331 - t327 * t329 + (-t354 * t328 + t332 * t381) * t348 + t382 ((-t327 * t344 + t380) * t341 + (t344 * t385 - t317) * t340) * pkin(5) + t382, t382;];
JaD_transl  = t1;
