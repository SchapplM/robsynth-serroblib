% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR10_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR10_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:50
% EndTime: 2019-02-26 22:35:50
% DurationCPUTime: 0.37s
% Computational Cost: add. (657->83), mult. (1178->130), div. (0->0), fcn. (1139->10), ass. (0->56)
t383 = pkin(4) - r_i_i_C(2);
t381 = r_i_i_C(3) + qJ(5);
t342 = cos(pkin(6));
t345 = sin(qJ(1));
t344 = sin(qJ(2));
t374 = t345 * t344;
t367 = t342 * t374;
t369 = qJD(2) * t344;
t347 = cos(qJ(2));
t348 = cos(qJ(1));
t371 = t348 * t347;
t321 = -qJD(1) * t367 - t345 * t369 + (qJD(2) * t342 + qJD(1)) * t371;
t340 = qJ(3) + qJ(4);
t337 = sin(t340);
t338 = cos(t340);
t339 = qJD(3) + qJD(4);
t372 = t348 * t344;
t373 = t345 * t347;
t328 = t342 * t372 + t373;
t341 = sin(pkin(6));
t370 = qJD(1) * t341;
t365 = t345 * t370;
t356 = -t328 * t339 + t365;
t375 = t341 * t348;
t385 = (t339 * t375 - t321) * t337 + t356 * t338;
t346 = cos(qJ(3));
t336 = pkin(3) * t346 + pkin(2);
t384 = t337 * t381 + t338 * t383 + t336;
t382 = r_i_i_C(1) + pkin(10) + pkin(9);
t357 = t367 - t371;
t380 = t357 * t338;
t379 = t337 * t339;
t378 = t341 * t344;
t377 = t341 * t345;
t376 = t341 * t346;
t368 = t338 * t378;
t366 = t338 * t375;
t364 = t348 * t370;
t363 = qJD(2) * t341 * t347;
t362 = -t321 * t338 + t339 * t366;
t358 = t342 * t373 + t372;
t319 = qJD(1) * t328 + qJD(2) * t358;
t361 = t339 * t377 - t319;
t359 = t342 * t371 - t374;
t355 = t339 * t342 + t363;
t307 = t337 * t361 - t338 * t364 - t339 * t380;
t308 = t337 * t364 + t338 * t361 + t357 * t379;
t354 = -(-t337 * t377 + t380) * qJD(5) + t381 * t308 - t383 * t307;
t353 = -(-t328 * t338 + t337 * t375) * qJD(5) + t381 * (-t328 * t379 + t337 * t365 - t362) + t383 * t385;
t316 = t337 * t355 + t339 * t368;
t352 = -(-t337 * t342 - t368) * qJD(5) + t381 * (t338 * t355 - t378 * t379) - t383 * t316;
t343 = sin(qJ(3));
t350 = -qJD(3) * t343 * pkin(3) + t337 * qJD(5) + (-t337 * t383 + t338 * t381) * t339;
t320 = qJD(1) * t358 + qJD(2) * t328;
t318 = -qJD(1) * t359 + qJD(2) * t357;
t1 = [-(t328 * t337 + t366) * qJD(5) - t321 * t336 - t382 * t320 + t383 * (-t337 * t356 + t362) + t381 * t385 + (-pkin(1) * t348 - pkin(8) * t377) * qJD(1) + (-t343 * t365 + (t328 * t343 + t346 * t375) * qJD(3)) * pkin(3), t318 * t384 - t319 * t382 - t350 * t358 (t346 * t364 + t319 * t343 + (-t343 * t377 + t346 * t357) * qJD(3)) * pkin(3) + t354, t354, t307, 0; -(t337 * t357 + t338 * t377) * qJD(5) - t319 * t336 - t382 * t318 + t383 * t308 + t381 * t307 + (-pkin(1) * t345 + pkin(8) * t375) * qJD(1) + (t343 * t364 + (t343 * t357 + t345 * t376) * qJD(3)) * pkin(3), -t320 * t384 + t321 * t382 + t350 * t359 (t346 * t365 - t321 * t343 + (-t328 * t346 + t343 * t375) * qJD(3)) * pkin(3) + t353, t353, -t385, 0; 0 (-t384 * t369 + (qJD(2) * t382 + t350) * t347) * t341 (-t343 * t363 + (-t342 * t343 - t344 * t376) * qJD(3)) * pkin(3) + t352, t352, t316, 0;];
JaD_transl  = t1;
