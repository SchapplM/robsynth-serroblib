% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:18
% EndTime: 2019-02-26 22:18:18
% DurationCPUTime: 0.45s
% Computational Cost: add. (610->80), mult. (663->111), div. (0->0), fcn. (506->10), ass. (0->68)
t315 = qJ(2) + qJ(3);
t311 = cos(t315);
t314 = qJ(5) + qJ(6);
t308 = sin(t314);
t310 = cos(t314);
t316 = sin(qJ(5));
t342 = pkin(5) * t316 + qJ(4);
t379 = r_i_i_C(1) * t308 + r_i_i_C(2) * t310 + t342;
t386 = t311 * t379;
t319 = cos(qJ(5));
t371 = pkin(5) * t319;
t382 = qJD(5) * t371 + qJD(4);
t385 = pkin(3) + r_i_i_C(3);
t309 = sin(t315);
t321 = cos(qJ(1));
t356 = qJD(1) * t321;
t344 = t309 * t356;
t318 = sin(qJ(1));
t313 = qJD(2) + qJD(3);
t365 = t311 * t313;
t347 = t318 * t365;
t384 = t344 + t347;
t312 = qJD(5) + qJD(6);
t366 = t311 * t312;
t348 = t310 * t366;
t381 = r_i_i_C(1) * t348 + t382 * t311;
t317 = sin(qJ(2));
t367 = pkin(2) * qJD(2);
t350 = t317 * t367;
t322 = -pkin(10) - pkin(9);
t353 = -t322 + t385;
t380 = (t353 * t313 - t382) * t309 - (pkin(4) + t371 + pkin(8) + pkin(7)) * qJD(1) - t342 * t365 + t350;
t340 = t309 * t312 + qJD(1);
t378 = t318 * t340;
t377 = t321 * t340;
t373 = pkin(2) * t317;
t369 = r_i_i_C(2) * t308;
t364 = t313 * t309;
t363 = t313 * t319;
t362 = t313 * t321;
t358 = qJD(1) * t309;
t339 = -t312 - t358;
t329 = t339 * t321 - t347;
t276 = t308 * t378 + t329 * t310;
t277 = t329 * t308 - t310 * t378;
t361 = -t276 * r_i_i_C(1) + t277 * r_i_i_C(2);
t346 = t311 * t362;
t328 = t339 * t318 + t346;
t278 = -t308 * t377 + t328 * t310;
t279 = t328 * t308 + t310 * t377;
t360 = t278 * r_i_i_C(1) - t279 * r_i_i_C(2);
t357 = qJD(1) * t318;
t354 = qJD(5) * t316;
t352 = r_i_i_C(1) * t309 * t310;
t349 = t308 * t366;
t345 = t309 * t357;
t337 = qJD(5) + t358;
t335 = (-qJD(5) * t309 - qJD(1)) * t316;
t334 = t381 * t321 + t322 * t346 + t385 * t345;
t333 = -t312 * t369 - t313 * t385;
t332 = r_i_i_C(1) * t349 + r_i_i_C(2) * t348 + t313 * t352 - t364 * t369;
t330 = t381 * t318 + t384 * t322 + t356 * t386;
t320 = cos(qJ(2));
t327 = -pkin(5) * t354 + (-pkin(2) * t320 - t342 * t309 - t353 * t311 - pkin(1)) * qJD(1);
t326 = t312 * t352 + t322 * t364 + (t333 + t382) * t309 + t379 * t365;
t325 = -r_i_i_C(2) * t349 + (-t309 * t379 - t311 * t385) * t313;
t324 = -t320 * t367 + t325;
t1 = [t277 * r_i_i_C(1) + t276 * r_i_i_C(2) + t380 * t318 + t327 * t321 (-t309 * t322 + t373 - t386) * t357 + t324 * t321 + t334 (-t322 * t357 - t362 * t379) * t309 + (t333 * t321 - t357 * t379) * t311 + t334, -t345 + t346 (t321 * t335 + (-t337 * t318 + t346) * t319) * pkin(5) + t360, t360; t279 * r_i_i_C(1) + t278 * r_i_i_C(2) + t327 * t318 - t380 * t321 (-t309 * t385 - t373) * t356 + t324 * t318 + t330, t325 * t318 - t344 * t385 + t330, t384 (t337 * t321 * t319 + (t311 * t363 + t335) * t318) * pkin(5) + t361, t361; 0, t326 - t350, t326, t364 (t309 * t363 + t311 * t354) * pkin(5) + t332, t332;];
JaD_transl  = t1;
