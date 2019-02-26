% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR11_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR11_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiaD_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:33
% EndTime: 2019-02-26 20:54:33
% DurationCPUTime: 0.32s
% Computational Cost: add. (281->64), mult. (951->112), div. (0->0), fcn. (978->12), ass. (0->55)
t355 = sin(qJ(3));
t349 = sin(pkin(7));
t350 = sin(pkin(6));
t358 = cos(qJ(1));
t385 = t350 * t358;
t375 = t349 * t385;
t353 = cos(pkin(7));
t352 = cos(pkin(12));
t354 = cos(pkin(6));
t382 = t354 * t358;
t348 = sin(pkin(12));
t356 = sin(qJ(1));
t388 = t348 * t356;
t365 = t352 * t382 - t388;
t398 = t365 * t353;
t368 = -t398 + t375;
t400 = t368 * t355;
t376 = t354 * t388;
t379 = qJD(1) * t358;
t336 = -qJD(1) * t376 + t352 * t379;
t399 = -qJD(3) * t375 + t336;
t357 = cos(qJ(3));
t384 = t353 * t355;
t387 = t349 * t354;
t397 = (t348 * t357 + t352 * t384) * t350 + t355 * t387;
t340 = t352 * t358 - t376;
t381 = t356 * t352;
t364 = t348 * t358 + t354 * t381;
t386 = t350 * t356;
t367 = t349 * t386 - t364 * t353;
t396 = t340 * t357 + t367 * t355;
t335 = t364 * qJD(1);
t380 = qJD(1) * t357;
t373 = t349 * t350 * t380;
t383 = t353 * t357;
t338 = t348 * t382 + t381;
t390 = t338 * t357;
t395 = (-t365 * t384 - t390) * qJD(3) - t335 * t383 + t356 * t373 - t399 * t355;
t374 = qJD(1) * t386;
t394 = (-qJD(3) * t338 - t335 * t353 + t349 * t374) * t355 + (qJD(3) * t398 + t399) * t357;
t393 = r_i_i_C(3) + qJ(4);
t392 = t335 * t349;
t378 = t350 * qJD(2);
t347 = sin(pkin(13));
t351 = cos(pkin(13));
t370 = -t351 * r_i_i_C(1) + t347 * r_i_i_C(2) - pkin(3);
t360 = -t340 * t355 + t367 * t357;
t359 = t365 * t349 + t353 * t385;
t334 = t338 * qJD(1);
t331 = -t353 * t374 - t392;
t330 = t359 * qJD(1);
t329 = t397 * qJD(3);
t323 = qJD(1) * t400 + t360 * qJD(3) - t334 * t357;
t322 = t396 * qJD(3) - t334 * t355 - t358 * t373 + t380 * t398;
t1 = [(t331 * t347 - t351 * t394) * r_i_i_C(1) + (t331 * t351 + t347 * t394) * r_i_i_C(2) - t394 * pkin(3) - (t338 * t355 + t368 * t357) * qJD(4) - t336 * pkin(2) - pkin(9) * t392 + t358 * t378 + t393 * t395 + (-t358 * pkin(1) + (-pkin(9) * t353 - qJ(2)) * t386) * qJD(1), t350 * t379, t396 * qJD(4) + t370 * t322 + t393 * t323, t322, 0, 0; (t323 * t351 + t330 * t347) * r_i_i_C(1) + (-t323 * t347 + t330 * t351) * r_i_i_C(2) + t323 * pkin(3) - t360 * qJD(4) - t334 * pkin(2) + t356 * t378 + t393 * t322 + (-t356 * pkin(1) + t359 * pkin(9) + qJ(2) * t385) * qJD(1), t374 -(-t390 + t400) * qJD(4) + t393 * t394 - t370 * t395, -t395, 0, 0; 0, 0, t397 * qJD(4) - t393 * (-t357 * t387 + (t348 * t355 - t352 * t383) * t350) * qJD(3) + t370 * t329, t329, 0, 0;];
JaD_transl  = t1;
