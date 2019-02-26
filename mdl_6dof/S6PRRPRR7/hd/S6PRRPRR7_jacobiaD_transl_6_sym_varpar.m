% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:38
% EndTime: 2019-02-26 20:07:39
% DurationCPUTime: 0.44s
% Computational Cost: add. (550->76), mult. (1306->127), div. (0->0), fcn. (1310->12), ass. (0->61)
t343 = sin(qJ(3));
t346 = cos(qJ(3));
t338 = qJD(5) + qJD(6);
t345 = cos(qJ(5));
t339 = qJ(5) + qJ(6);
t336 = sin(t339);
t337 = cos(t339);
t364 = t337 * r_i_i_C(1) - t336 * r_i_i_C(2);
t388 = qJD(5) * pkin(5);
t351 = t338 * t364 + t345 * t388 + qJD(4);
t342 = sin(qJ(5));
t363 = -t336 * r_i_i_C(1) - t337 * r_i_i_C(2);
t353 = pkin(5) * t342 + qJ(4) - t363;
t377 = pkin(3) + r_i_i_C(3) + pkin(10) + pkin(9);
t349 = (t343 * t377 - t346 * t353) * qJD(3) - t351 * t343;
t340 = sin(pkin(11));
t344 = sin(qJ(2));
t347 = cos(qJ(2));
t386 = cos(pkin(11));
t387 = cos(pkin(6));
t361 = t387 * t386;
t326 = t340 * t347 + t344 * t361;
t341 = sin(pkin(6));
t373 = t341 * t386;
t390 = t326 * t346 - t343 * t373;
t374 = t340 * t387;
t328 = -t344 * t374 + t386 * t347;
t384 = t341 * t343;
t383 = t341 * t346;
t382 = t341 * t347;
t322 = t326 * qJD(2);
t354 = -t326 * t343 - t346 * t373;
t369 = t338 * t354 - t322;
t360 = t347 * t361;
t378 = qJD(2) * t344;
t321 = -qJD(2) * t360 + t340 * t378;
t311 = t390 * qJD(3) - t321 * t343;
t325 = t340 * t344 - t360;
t371 = t325 * t338 - t311;
t381 = (t336 * t369 - t337 * t371) * r_i_i_C(1) + (t336 * t371 + t337 * t369) * r_i_i_C(2);
t324 = t328 * qJD(2);
t359 = -t328 * t343 + t340 * t383;
t368 = t338 * t359 - t324;
t327 = t344 * t386 + t347 * t374;
t323 = t327 * qJD(2);
t358 = t328 * t346 + t340 * t384;
t313 = qJD(3) * t358 - t323 * t343;
t370 = t327 * t338 - t313;
t380 = (t336 * t368 - t337 * t370) * r_i_i_C(1) + (t336 * t370 + t337 * t368) * r_i_i_C(2);
t329 = t344 * t384 - t346 * t387;
t376 = t341 * t378;
t357 = -t329 * t338 - t376;
t355 = t343 * t387 + t344 * t383;
t375 = qJD(2) * t382;
t319 = qJD(3) * t355 + t343 * t375;
t362 = t338 * t382 + t319;
t379 = (t336 * t357 + t337 * t362) * r_i_i_C(1) + (-t336 * t362 + t337 * t357) * r_i_i_C(2);
t356 = t345 * pkin(5) + pkin(4) + pkin(8) + t364;
t352 = t338 * t363 - t342 * t388;
t350 = -t343 * t353 - t346 * t377 - pkin(2);
t1 = [0, -t323 * t356 + t324 * t350 + t327 * t349 + t328 * t352, -t377 * t313 + t351 * t358 + t353 * (qJD(3) * t359 - t323 * t346) t313 (t313 * t345 - t324 * t342 + (-t327 * t345 + t342 * t359) * qJD(5)) * pkin(5) + t380, t380; 0, -t321 * t356 + t322 * t350 + t325 * t349 + t326 * t352, -t377 * t311 + t351 * t390 + t353 * (qJD(3) * t354 - t321 * t346) t311 (t311 * t345 - t322 * t342 + (-t325 * t345 + t342 * t354) * qJD(5)) * pkin(5) + t381, t381; 0 ((qJD(2) * t350 + t352) * t344 + (t356 * qJD(2) - t349) * t347) * t341, -t377 * t319 + t351 * t355 + t353 * (-qJD(3) * t329 + t346 * t375) t319 (-t342 * t376 + t319 * t345 + (-t329 * t342 + t345 * t382) * qJD(5)) * pkin(5) + t379, t379;];
JaD_transl  = t1;
