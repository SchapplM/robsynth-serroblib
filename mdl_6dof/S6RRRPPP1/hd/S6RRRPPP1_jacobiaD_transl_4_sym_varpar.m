% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPP1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPP1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:50
% EndTime: 2019-02-26 22:02:50
% DurationCPUTime: 0.73s
% Computational Cost: add. (291->88), mult. (998->154), div. (0->0), fcn. (918->10), ass. (0->56)
t325 = sin(qJ(2));
t320 = sin(pkin(10));
t322 = cos(pkin(10));
t324 = sin(qJ(3));
t327 = cos(qJ(3));
t323 = cos(pkin(6));
t370 = t323 * t324;
t321 = sin(pkin(6));
t375 = r_i_i_C(3) + qJ(4);
t382 = t321 * t375;
t333 = -t324 * t382 + (t320 * t370 - t322 * t327) * r_i_i_C(1) + (t320 * t327 + t322 * t370) * r_i_i_C(2) - t327 * pkin(3);
t332 = -pkin(2) + t333;
t328 = cos(qJ(2));
t340 = pkin(9) + (r_i_i_C(1) * t320 + r_i_i_C(2) * t322) * t321;
t385 = t375 * t323 + t340;
t335 = t385 * t328;
t331 = t332 * t325 + t335;
t386 = t323 * qJ(4) + t340;
t369 = t323 * t327;
t334 = t327 * t382 - (t320 * t369 + t322 * t324) * r_i_i_C(1) - (-t320 * t324 + t322 * t369) * r_i_i_C(2) - t324 * pkin(3);
t356 = t321 * qJD(4);
t384 = t334 * qJD(3) + t324 * t356;
t377 = t325 * pkin(2);
t326 = sin(qJ(1));
t362 = qJD(1) * t328;
t346 = -qJD(3) + t362;
t347 = qJD(3) * t328 - qJD(1);
t329 = cos(qJ(1));
t358 = qJD(2) * t329;
t351 = t325 * t358;
t365 = t327 * t329;
t311 = -t347 * t365 + (t346 * t326 + t351) * t324;
t374 = t311 * t321;
t373 = t320 * t323;
t372 = t322 * t323;
t368 = t324 * t329;
t367 = t325 * t323;
t366 = t326 * t328;
t364 = t328 * t329;
t363 = qJD(1) * t326;
t361 = qJD(1) * t329;
t360 = qJD(2) * t326;
t359 = qJD(2) * t328;
t357 = qJD(3) * t327;
t355 = t323 * qJD(4);
t354 = -pkin(2) * t328 - pkin(1);
t352 = t325 * t355;
t350 = qJD(3) * t368;
t338 = t324 * t361 + t326 * t357;
t336 = t385 * t325;
t330 = t328 * t355 - t384 * t325 + (t332 * t328 - t336) * qJD(2);
t314 = -t346 * t365 + (qJD(2) * t325 * t327 + t347 * t324) * t326;
t313 = -t325 * t324 * t360 - t327 * t363 + t338 * t328 - t350;
t312 = t328 * t350 + (t326 * t362 + t351) * t327 - t338;
t310 = -t374 + (-t325 * t363 + t328 * t358) * t323;
t1 = [(t313 * t373 + t314 * t322) * r_i_i_C(1) + (t313 * t372 - t314 * t320) * r_i_i_C(2) + t314 * pkin(3) - t326 * t352 + (-(t324 * t366 + t365) * qJD(4) - t375 * t313) * t321 + (-t335 + t377) * t360 + (-t326 * pkin(8) + (-t336 + t354) * t329) * qJD(1), t330 * t329 - t331 * t363 (t311 * t322 + t312 * t373) * r_i_i_C(1) + (-t311 * t320 + t312 * t372) * r_i_i_C(2) + t311 * pkin(3) + (-(-t326 * t324 - t327 * t364) * qJD(4) - t375 * t312) * t321, t310, 0, 0; (t311 * t373 - t312 * t322) * r_i_i_C(1) + (t311 * t372 + t312 * t320) * r_i_i_C(2) + t310 * r_i_i_C(3) - t312 * pkin(3) - qJ(4) * t374 + (-(-t324 * t364 + t326 * t327) * t321 + t329 * t367) * qJD(4) + (t386 * t328 - t377) * t358 + (t329 * pkin(8) + (-t386 * t325 + t354) * t326) * qJD(1), t330 * t326 + t331 * t361 (-t313 * t322 + t314 * t373) * r_i_i_C(1) + (t313 * t320 + t314 * t372) * r_i_i_C(2) - t313 * pkin(3) + (-(-t327 * t366 + t368) * qJD(4) - t375 * t314) * t321, t313 * t321 + (t325 * t361 + t326 * t359) * t323, 0, 0; 0, t331 * qJD(2) + t384 * t328 + t352, t334 * t359 + (t333 * qJD(3) + t327 * t356) * t325, t321 * t325 * t357 + (t328 * t324 * t321 + t367) * qJD(2), 0, 0;];
JaD_transl  = t1;
