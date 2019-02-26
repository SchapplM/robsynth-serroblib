% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRPRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRPRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:49
% EndTime: 2019-02-26 19:39:50
% DurationCPUTime: 0.28s
% Computational Cost: add. (265->66), mult. (885->136), div. (0->0), fcn. (967->14), ass. (0->55)
t337 = sin(pkin(13));
t342 = cos(pkin(13));
t350 = cos(qJ(3));
t362 = qJD(3) * t350;
t348 = sin(qJ(3));
t363 = qJD(3) * t348;
t374 = t337 * t363 - t342 * t362;
t340 = sin(pkin(7));
t316 = t374 * t340;
t345 = cos(pkin(7));
t318 = t374 * t345;
t330 = -t337 * t362 - t342 * t363;
t338 = sin(pkin(12));
t341 = sin(pkin(6));
t343 = cos(pkin(12));
t346 = cos(pkin(6));
t373 = t346 * t316 + (t318 * t343 - t330 * t338) * t341;
t372 = -pkin(9) - r_i_i_C(3);
t371 = pkin(3) * qJD(3);
t339 = sin(pkin(11));
t370 = t339 * t341;
t369 = t339 * t346;
t368 = t340 * t341;
t344 = cos(pkin(11));
t367 = t341 * t344;
t366 = t341 * t345;
t365 = t344 * t346;
t347 = sin(qJ(5));
t349 = cos(qJ(5));
t359 = -t347 * r_i_i_C(1) - t349 * r_i_i_C(2);
t357 = t350 * t337 + t348 * t342;
t356 = t348 * t337 - t350 * t342;
t355 = t349 * r_i_i_C(1) - t347 * r_i_i_C(2) + pkin(4);
t354 = qJD(5) * t359;
t353 = qJD(3) * t357;
t325 = -t339 * t338 + t343 * t365;
t326 = t338 * t365 + t339 * t343;
t352 = t316 * t367 - t325 * t318 + t326 * t330;
t327 = -t344 * t338 - t343 * t369;
t328 = -t338 * t369 + t344 * t343;
t351 = t316 * t370 + t327 * t318 - t328 * t330;
t329 = t356 * qJD(3);
t324 = -t343 * t368 + t346 * t345;
t323 = t357 * t345;
t322 = t356 * t345;
t321 = t357 * t340;
t320 = t356 * t340;
t319 = t345 * t353;
t317 = t340 * t353;
t315 = -t327 * t340 + t339 * t366;
t314 = -t325 * t340 - t344 * t366;
t313 = t346 * t321 + (t323 * t343 - t338 * t356) * t341;
t308 = t321 * t370 + t327 * t323 - t328 * t356;
t306 = -t321 * t367 + t325 * t323 - t326 * t356;
t1 = [0, 0, t372 * t351 + (-t320 * t370 - t327 * t322 - t328 * t357) * t354 + t355 * (-t317 * t370 - t327 * t319 + t328 * t329) + (-t328 * t350 + (-t327 * t345 - t339 * t368) * t348) * t371, 0, -t359 * t351 + ((-t308 * t349 - t315 * t347) * r_i_i_C(1) + (t308 * t347 - t315 * t349) * r_i_i_C(2)) * qJD(5), 0; 0, 0, -t372 * t352 + (t320 * t367 - t325 * t322 - t326 * t357) * t354 + t355 * (t317 * t367 - t325 * t319 + t326 * t329) + (-t326 * t350 + (-t325 * t345 + t340 * t367) * t348) * t371, 0, t359 * t352 + ((-t306 * t349 - t314 * t347) * r_i_i_C(1) + (t306 * t347 - t314 * t349) * r_i_i_C(2)) * qJD(5), 0; 0, 0, t372 * t373 + (-t346 * t320 + (-t322 * t343 - t338 * t357) * t341) * t354 + t355 * (-t346 * t317 + (-t319 * t343 + t329 * t338) * t341) + (-t340 * t346 * t348 + (-t343 * t345 * t348 - t338 * t350) * t341) * t371, 0, -t359 * t373 + ((-t313 * t349 - t324 * t347) * r_i_i_C(1) + (t313 * t347 - t324 * t349) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
