% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR6
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
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:56
% EndTime: 2019-02-26 19:56:56
% DurationCPUTime: 0.37s
% Computational Cost: add. (499->77), mult. (1130->129), div. (0->0), fcn. (1134->12), ass. (0->63)
t335 = sin(pkin(11));
t337 = cos(pkin(11));
t340 = sin(qJ(2));
t343 = cos(qJ(2));
t381 = cos(pkin(6));
t365 = t343 * t381;
t322 = t335 * t365 + t337 * t340;
t342 = cos(qJ(4));
t336 = sin(pkin(6));
t339 = sin(qJ(4));
t379 = t336 * t339;
t385 = -t322 * t342 + t335 * t379;
t333 = qJD(5) + qJD(6);
t341 = cos(qJ(5));
t334 = qJ(5) + qJ(6);
t331 = sin(t334);
t332 = cos(t334);
t358 = t332 * r_i_i_C(1) - t331 * r_i_i_C(2);
t382 = pkin(5) * qJD(5);
t347 = -t358 * t333 - t341 * t382;
t355 = t341 * pkin(5) + pkin(4) + t358;
t383 = r_i_i_C(3) + pkin(10) + pkin(9);
t384 = t355 * t339 - t383 * t342 + qJ(3);
t378 = t336 * t340;
t377 = t336 * t342;
t376 = t336 * t343;
t309 = t322 * t339 + t335 * t377;
t318 = t322 * qJD(2);
t362 = t309 * t333 + t318;
t366 = t340 * t381;
t360 = t335 * t366;
t371 = qJD(2) * t343;
t319 = -qJD(2) * t360 + t337 * t371;
t304 = t385 * qJD(4) - t319 * t339;
t323 = t337 * t343 - t360;
t364 = -t323 * t333 + t304;
t375 = (t364 * t331 - t362 * t332) * r_i_i_C(1) + (t362 * t331 + t364 * t332) * r_i_i_C(2);
t359 = t337 * t365;
t372 = qJD(2) * t340;
t316 = -qJD(2) * t359 + t335 * t372;
t320 = t335 * t340 - t359;
t354 = -t320 * t339 + t337 * t377;
t361 = -t333 * t354 + t316;
t321 = t335 * t343 + t337 * t366;
t317 = t321 * qJD(2);
t353 = t320 * t342 + t337 * t379;
t306 = t353 * qJD(4) + t317 * t339;
t363 = -t321 * t333 - t306;
t374 = (t363 * t331 - t361 * t332) * r_i_i_C(1) + (t361 * t331 + t363 * t332) * r_i_i_C(2);
t351 = t339 * t376 - t381 * t342;
t368 = t336 * t371;
t352 = t333 * t351 + t368;
t350 = t381 * t339 + t342 * t376;
t367 = t336 * t372;
t312 = t350 * qJD(4) - t339 * t367;
t356 = -t333 * t378 + t312;
t373 = (t356 * t331 + t352 * t332) * r_i_i_C(1) + (-t352 * t331 + t356 * t332) * r_i_i_C(2);
t357 = -t331 * r_i_i_C(1) - t332 * r_i_i_C(2);
t338 = sin(qJ(5));
t349 = -t338 * pkin(5) - pkin(2) - pkin(8) + t357;
t348 = t357 * t333 - t338 * t382;
t345 = qJD(3) + t348 * t339 + (t383 * t339 + t355 * t342) * qJD(4);
t1 = [0, -t318 * t384 + t349 * t319 + t347 * t322 + t345 * t323, t319, -t383 * t304 - t348 * t385 + t355 * (-t309 * qJD(4) + t319 * t342) (t304 * t338 - t318 * t341 + (-t309 * t341 - t323 * t338) * qJD(5)) * pkin(5) + t375, t375; 0, -t316 * t384 + t349 * t317 + t347 * t320 + t345 * t321, t317, t383 * t306 + t348 * t353 + t355 * (t354 * qJD(4) + t317 * t342) (-t306 * t338 - t316 * t341 + (-t321 * t338 + t341 * t354) * qJD(5)) * pkin(5) + t374, t374; 0 ((t384 * qJD(2) - t347) * t343 + (t349 * qJD(2) + t345) * t340) * t336, t367, -t383 * t312 - t348 * t350 + t355 * (t351 * qJD(4) + t342 * t367) (t341 * t368 + t312 * t338 + (-t338 * t378 + t341 * t351) * qJD(5)) * pkin(5) + t373, t373;];
JaD_transl  = t1;
