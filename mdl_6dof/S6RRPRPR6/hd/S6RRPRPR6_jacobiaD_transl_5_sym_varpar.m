% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:39
% EndTime: 2019-02-26 21:40:40
% DurationCPUTime: 0.42s
% Computational Cost: add. (538->86), mult. (1584->146), div. (0->0), fcn. (1696->10), ass. (0->64)
t353 = sin(pkin(11));
t360 = cos(qJ(2));
t393 = cos(pkin(11));
t376 = qJD(2) * t393;
t357 = sin(qJ(2));
t382 = qJD(2) * t357;
t400 = t353 * t382 - t360 * t376;
t355 = cos(pkin(6));
t334 = t400 * t355;
t345 = -t360 * t353 - t357 * t393;
t340 = t345 * t355;
t361 = cos(qJ(1));
t368 = -t357 * t353 + t360 * t393;
t358 = sin(qJ(1));
t384 = qJD(1) * t358;
t381 = qJD(2) * t360;
t342 = -t353 * t381 - t357 * t376;
t388 = t358 * t342;
t323 = t388 + t340 * t384 + (qJD(1) * t368 - t334) * t361;
t356 = sin(qJ(4));
t359 = cos(qJ(4));
t328 = -t361 * t340 + t358 * t368;
t354 = sin(pkin(6));
t391 = t354 * t361;
t371 = t328 * t356 + t359 * t391;
t379 = t354 * t384;
t399 = t371 * qJD(4) - t323 * t359 - t356 * t379;
t370 = -t328 * t359 + t356 * t391;
t398 = t370 * qJD(4) - t323 * t356 + t359 * t379;
t394 = r_i_i_C(3) + qJ(5);
t397 = pkin(4) - r_i_i_C(2);
t365 = t394 * t356 + t397 * t359 + pkin(3);
t396 = pkin(9) + r_i_i_C(1);
t395 = pkin(2) * t355;
t372 = -t358 * t340 - t361 * t368;
t392 = t372 * t356;
t390 = t357 * t358;
t389 = t357 * t361;
t387 = t358 * t359;
t386 = t358 * t360;
t385 = t360 * t361;
t383 = qJD(1) * t361;
t380 = pkin(2) * t382;
t378 = t354 * t383;
t338 = t345 * t354;
t374 = -t338 * t359 + t355 * t356;
t339 = t368 * t355;
t373 = t361 * t339 + t358 * t345;
t369 = t358 * t354 * t356 - t359 * t372;
t364 = qJD(2) * t345;
t363 = t368 * qJD(2);
t362 = qJD(5) * t356 + (-t397 * t356 + t394 * t359) * qJD(4);
t320 = -t328 * qJD(1) + t358 * t334 + t361 * t342;
t352 = t360 * pkin(2) + pkin(1);
t343 = -t354 * qJD(3) + t381 * t395;
t341 = t357 * t395 + (-pkin(8) - qJ(3)) * t354;
t335 = t355 * t364;
t332 = t400 * t354;
t324 = t374 * qJD(4) - t332 * t356;
t322 = t361 * t335 - t339 * t384 + t345 * t383 - t358 * t363;
t319 = t373 * qJD(1) + t358 * t335 + t361 * t363;
t313 = t320 * t359 + qJD(4) * t392 + (qJD(4) * t387 + t356 * t383) * t354;
t312 = t369 * qJD(4) + t320 * t356 - t359 * t378;
t1 = [-t371 * qJD(5) - t323 * pkin(3) + t358 * t380 - t361 * t343 + t396 * t322 + t397 * t399 + t394 * t398 + (t358 * t341 - t361 * t352) * qJD(1), t396 * t320 + ((t355 * t390 - t385) * qJD(2) + (-t355 * t385 + t390) * qJD(1)) * pkin(2) + t362 * (-t358 * t339 + t361 * t345) - t365 * t319, t378, t369 * qJD(5) - t397 * t312 + t394 * t313, t312, 0; -(t354 * t387 + t392) * qJD(5) + t320 * pkin(3) - t361 * t380 - t358 * t343 + t396 * t319 + t397 * t313 + t394 * t312 + (-t361 * t341 - t358 * t352) * qJD(1), -t396 * (t372 * qJD(1) + t361 * t334 - t388) + ((-t355 * t389 - t386) * qJD(2) + (-t355 * t386 - t389) * qJD(1)) * pkin(2) + t362 * t373 + t365 * t322, t379, -t370 * qJD(5) - t394 * t399 + t397 * t398, -t398, 0; 0, -t396 * t332 + (t362 * t368 + t364 * t365 - t380) * t354, 0, t374 * qJD(5) + t394 * (-t332 * t359 + (t338 * t356 + t355 * t359) * qJD(4)) - t397 * t324, t324, 0;];
JaD_transl  = t1;
