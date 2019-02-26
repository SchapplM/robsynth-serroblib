% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR10_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR10_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:48
% EndTime: 2019-02-26 22:08:49
% DurationCPUTime: 0.32s
% Computational Cost: add. (445->70), mult. (1325->111), div. (0->0), fcn. (1308->10), ass. (0->54)
t323 = cos(pkin(6));
t325 = sin(qJ(2));
t329 = cos(qJ(1));
t354 = t329 * t325;
t326 = sin(qJ(1));
t328 = cos(qJ(2));
t355 = t326 * t328;
t309 = t323 * t354 + t355;
t324 = sin(qJ(3));
t327 = cos(qJ(3));
t321 = sin(pkin(6));
t358 = t321 * t329;
t365 = t309 * t324 + t327 * t358;
t341 = qJD(2) * t323 + qJD(1);
t356 = t326 * t325;
t347 = t323 * t356;
t352 = qJD(2) * t325;
t353 = t329 * t328;
t302 = -qJD(1) * t347 - t326 * t352 + t341 * t353;
t360 = t321 * t326;
t348 = t324 * t360;
t364 = -qJD(1) * t348 + qJD(3) * t365 - t302 * t327;
t320 = sin(pkin(11));
t322 = cos(pkin(11));
t340 = t320 * r_i_i_C(1) + t322 * r_i_i_C(2) + qJ(4);
t349 = r_i_i_C(3) + qJ(5) + pkin(3);
t363 = t340 * t324 + t349 * t327 + pkin(2);
t311 = -t347 + t353;
t361 = t311 * t324;
t359 = t321 * t327;
t357 = t324 * t329;
t351 = qJD(2) * t328;
t350 = qJD(3) * t327;
t346 = t321 * t357;
t345 = t323 * t353;
t344 = qJD(1) * t359;
t343 = t321 * t351;
t339 = t322 * r_i_i_C(1) - t320 * r_i_i_C(2) + pkin(4) + pkin(9);
t337 = -t309 * t327 + t346;
t336 = t311 * t327 + t348;
t335 = t323 * t324 + t325 * t359;
t334 = t321 * t325 * t324 - t323 * t327;
t333 = t323 * t355 + t354;
t295 = -qJD(3) * t346 + t302 * t324 + t309 * t350 - t326 * t344;
t330 = t324 * qJD(4) + t327 * qJD(5) + (-t349 * t324 + t340 * t327) * qJD(3);
t305 = -t326 * t359 + t361;
t304 = -qJD(3) * t334 + t327 * t343;
t303 = qJD(3) * t335 + t324 * t343;
t301 = qJD(1) * t333 + qJD(2) * t309;
t300 = qJD(1) * t309 + qJD(2) * t333;
t299 = -qJD(1) * t345 - t329 * t351 + t341 * t356;
t294 = -t300 * t327 - qJD(3) * t361 + (qJD(1) * t357 + t326 * t350) * t321;
t293 = qJD(3) * t336 - t300 * t324 - t329 * t344;
t1 = [t337 * qJD(5) - t365 * qJD(4) - t302 * pkin(2) - t340 * t295 - t339 * t301 + (-t329 * pkin(1) - pkin(8) * t360) * qJD(1) + t349 * t364, t363 * t299 - t339 * t300 - t330 * t333, qJD(4) * t336 - t305 * qJD(5) - t349 * t293 + t340 * t294, t293, t294, 0; -t300 * pkin(2) + t305 * qJD(4) + t336 * qJD(5) + t340 * t293 - t339 * t299 + (-t326 * pkin(1) + pkin(8) * t358) * qJD(1) + t349 * t294, t339 * t302 - t363 * t301 + t330 * (t345 - t356) -t337 * qJD(4) - qJD(5) * t365 - t349 * t295 - t340 * t364, t295, -t364, 0; 0 (-t363 * t352 + (qJD(2) * t339 + t330) * t328) * t321, t335 * qJD(4) - t334 * qJD(5) - t349 * t303 + t340 * t304, t303, t304, 0;];
JaD_transl  = t1;
