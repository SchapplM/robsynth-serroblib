% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:40:11
% EndTime: 2019-02-26 22:40:12
% DurationCPUTime: 0.44s
% Computational Cost: add. (953->84), mult. (900->108), div. (0->0), fcn. (706->10), ass. (0->66)
t330 = sin(qJ(5));
t333 = cos(qJ(5));
t335 = cos(qJ(1));
t375 = qJD(5) * t335;
t332 = sin(qJ(1));
t378 = qJD(1) * t332;
t346 = t330 * t375 + t333 * t378;
t389 = pkin(5) + r_i_i_C(1);
t385 = r_i_i_C(3) + qJ(6);
t393 = t385 * t330;
t329 = qJ(2) + qJ(3);
t326 = qJ(4) + t329;
t322 = cos(t326);
t388 = pkin(10) + r_i_i_C(2);
t367 = t388 * t322;
t327 = qJD(2) + qJD(3);
t323 = qJD(4) + t327;
t366 = t388 * t323;
t331 = sin(qJ(2));
t384 = pkin(2) * qJD(2);
t370 = t331 * t384;
t324 = sin(t329);
t386 = pkin(3) * t327;
t372 = t324 * t386;
t374 = qJD(6) * t330;
t321 = sin(t326);
t383 = t321 * t323;
t392 = -pkin(4) * t383 + (t366 + t374) * t322 - t370 - t372;
t376 = qJD(5) * t333;
t344 = -t385 * t376 - t374;
t325 = cos(t329);
t347 = -t389 * t333 - t393;
t343 = -pkin(4) + t347;
t368 = t388 * t321;
t338 = t344 * t321 + (t343 * t322 - t368) * t323;
t336 = -t325 * t386 + t338;
t387 = pkin(3) * t324;
t382 = t322 * t323;
t381 = t322 * t332;
t380 = t332 * t330;
t379 = t333 * t335;
t377 = qJD(1) * t335;
t373 = t333 * qJD(6);
t369 = t389 * t330;
t365 = t332 * t383;
t364 = t335 * t383;
t359 = qJD(5) * t380;
t357 = t333 * t375;
t356 = t389 * t321 * t359 + t377 * t367;
t350 = t322 * t379 + t380;
t334 = cos(qJ(2));
t349 = -pkin(2) * t334 - pkin(3) * t325 - pkin(4) * t322 - pkin(1) - t368;
t348 = (t389 * t346 + (pkin(4) + t393) * t378) * t321;
t345 = t330 * t377 + t332 * t376;
t342 = t343 * t321;
t341 = t343 * t323;
t340 = t321 * t341 + t388 * t382 + (-qJD(5) * t369 - t344) * t322;
t339 = t340 - t372;
t337 = -t334 * t384 + t336;
t328 = -pkin(9) - pkin(8) - pkin(7);
t313 = -pkin(2) * t331 - t387;
t290 = t350 * qJD(1) - t322 * t359 - t333 * t365 - t357;
t289 = t345 * t322 - t330 * t365 - t346;
t288 = t346 * t322 + t333 * t364 - t345;
t287 = t330 * t364 - t322 * t357 - t359 + (t322 * t380 + t379) * qJD(1);
t1 = [-t335 * t373 - t389 * t290 - t385 * t289 - t392 * t332 + (t332 * t328 + t349 * t335) * qJD(1) (-t313 - t367) * t378 + t337 * t335 + t348 (-t367 + t387) * t378 + t336 * t335 + t348, t338 * t335 - t367 * t378 + t348, t350 * qJD(6) + t389 * t287 - t385 * t288, -t287; -t332 * t373 - t389 * t288 - t385 * t287 + t392 * t335 + (-t335 * t328 + t349 * t332) * qJD(1) (t313 + t342) * t377 + t337 * t332 + t356 (t342 - t387) * t377 + t336 * t332 + t356, t341 * t381 + ((-t366 + t344) * t332 + t343 * t377) * t321 + t356 -(t330 * t335 - t333 * t381) * qJD(6) + t385 * t290 - t389 * t289, t289; 0, t339 - t370, t339, t340 (t385 * t333 - t369) * t382 + (t347 * qJD(5) + t373) * t321, t321 * t376 + t330 * t382;];
JaD_transl  = t1;
