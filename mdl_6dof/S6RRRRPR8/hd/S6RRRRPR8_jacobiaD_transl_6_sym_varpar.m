% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:44
% EndTime: 2019-02-26 22:34:45
% DurationCPUTime: 0.70s
% Computational Cost: add. (1124->113), mult. (1546->166), div. (0->0), fcn. (1441->10), ass. (0->75)
t339 = qJ(3) + qJ(4);
t336 = sin(t339);
t343 = sin(qJ(1));
t346 = cos(qJ(2));
t388 = t343 * t346;
t337 = cos(t339);
t347 = cos(qJ(1));
t393 = t337 * t347;
t317 = t336 * t388 + t393;
t342 = sin(qJ(2));
t383 = qJD(2) * t342;
t371 = t347 * t383;
t338 = qJD(3) + qJD(4);
t374 = t338 * t393;
t394 = t336 * t343;
t376 = t338 * t394;
t313 = t317 * qJD(1) + t336 * t371 - t346 * t374 - t376;
t384 = qJD(1) * t347;
t389 = t343 * t337;
t355 = t336 * t384 + t338 * t389;
t387 = t347 * t336;
t375 = t338 * t387;
t385 = qJD(1) * t346;
t314 = t346 * t375 + (t343 * t385 + t371) * t337 - t355;
t319 = t346 * t387 - t389;
t320 = t346 * t393 + t394;
t340 = sin(qJ(6));
t344 = cos(qJ(6));
t303 = -t313 * t340 - t314 * t344 + (t319 * t344 - t320 * t340) * qJD(6);
t403 = t313 * t344 - t314 * t340 + (t319 * t340 + t320 * t344) * qJD(6);
t414 = -t403 * r_i_i_C(1) - r_i_i_C(2) * t303;
t366 = r_i_i_C(1) * t344 - r_i_i_C(2) * t340;
t413 = pkin(4) + pkin(5);
t356 = t366 + t413;
t341 = sin(qJ(3));
t396 = pkin(3) * qJD(3);
t379 = t341 * t396;
t380 = -r_i_i_C(3) - pkin(10) + pkin(9) + pkin(8);
t381 = qJD(6) - t338;
t397 = r_i_i_C(2) * t344;
t398 = r_i_i_C(1) * t340;
t349 = (-t338 * qJ(5) + (t397 + t398) * t381) * t337 - (t366 * qJD(6) - t356 * t338 + qJD(5)) * t336 - t380 * qJD(2) + t379;
t345 = cos(qJ(3));
t335 = pkin(3) * t345 + pkin(2);
t410 = (-t342 * t335 + t380 * t346) * qJD(2) - t346 * t379;
t359 = t336 * t340 + t337 * t344;
t360 = t336 * t344 - t337 * t340;
t382 = qJD(2) * t346;
t405 = t381 * t342;
t409 = -(t359 * t382 + t360 * t405) * r_i_i_C(2) - (t359 * t405 - t360 * t382) * r_i_i_C(1);
t372 = t343 * t383;
t386 = qJD(1) * t343;
t315 = -t336 * t372 - t337 * t386 + t355 * t346 - t375;
t308 = t315 * t344;
t316 = t320 * qJD(1) - t337 * t372 - t346 * t376 - t374;
t408 = -t316 * t340 + t308;
t370 = -qJ(5) - t398;
t402 = (-t370 + t397) * t336 + t356 * t337 + t335;
t399 = pkin(3) * t341;
t391 = t338 * t342;
t378 = t345 * t396;
t377 = t337 * t391;
t373 = pkin(7) + t399;
t368 = -qJD(3) + t385;
t367 = -t409 + (qJ(5) * t382 + qJD(5) * t342) * t337;
t364 = -t315 * t340 - t316 * t344;
t318 = t337 * t388 - t387;
t363 = t317 * t344 - t318 * t340;
t362 = t317 * t340 + t318 * t344;
t358 = (-qJD(3) * t346 + qJD(1)) * t345;
t354 = -t335 * t346 - t380 * t342 - pkin(1);
t353 = -t314 * qJ(5) + t320 * qJD(5) + t413 * t313 - t414;
t352 = t316 * qJ(5) + t318 * qJD(5) + (t363 * qJD(6) - t364) * r_i_i_C(2) + (t362 * qJD(6) - t408) * r_i_i_C(1) - t413 * t315;
t350 = qJD(2) * t402;
t1 = [t347 * t378 - t308 * r_i_i_C(2) - t317 * qJD(5) + t370 * t315 - t356 * t316 + (-t363 * r_i_i_C(1) + t362 * r_i_i_C(2)) * qJD(6) - t410 * t343 + (-t373 * t343 + t354 * t347) * qJD(1) (-t347 * t350 - t380 * t386) * t346 + (t349 * t347 + t402 * t386) * t342 (t347 * t358 + (t368 * t343 + t371) * t341) * pkin(3) + t353, t353, -t313, t414; t343 * t378 + t303 * r_i_i_C(1) - t403 * r_i_i_C(2) - t313 * qJ(5) + t319 * qJD(5) - t413 * t314 + t410 * t347 + (t354 * t343 + t373 * t347) * qJD(1) (-t343 * t350 + t380 * t384) * t346 + (t349 * t343 - t384 * t402) * t342 (t343 * t358 + (-t368 * t347 + t372) * t341) * pkin(3) + t352, t352, t315, t408 * r_i_i_C(1) + t364 * r_i_i_C(2) + (-t362 * r_i_i_C(1) - t363 * r_i_i_C(2)) * qJD(6); 0, -t342 * t350 - t349 * t346 (-t336 * t413 - t399) * t382 + (-t378 + (-t336 * qJ(5) - t337 * t413) * t338) * t342 + t367, -t413 * t377 + (-qJ(5) * t391 - t382 * t413) * t336 + t367, t336 * t382 + t377, t409;];
JaD_transl  = t1;
