% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:47
% EndTime: 2019-02-26 19:58:47
% DurationCPUTime: 0.40s
% Computational Cost: add. (462->76), mult. (966->136), div. (0->0), fcn. (960->12), ass. (0->50)
t316 = qJ(3) + pkin(11);
t314 = sin(t316);
t315 = cos(t316);
t321 = sin(qJ(3));
t320 = sin(qJ(6));
t323 = cos(qJ(6));
t339 = t323 * r_i_i_C(1) - t320 * r_i_i_C(2);
t328 = t339 * qJD(6) + qJD(5);
t338 = -t320 * r_i_i_C(1) - t323 * r_i_i_C(2);
t336 = qJ(5) - t338;
t347 = pkin(4) + pkin(9) + r_i_i_C(3);
t326 = (t321 * pkin(3) + t347 * t314 - t336 * t315) * qJD(3) - t328 * t314;
t317 = sin(pkin(10));
t322 = sin(qJ(2));
t325 = cos(qJ(2));
t353 = cos(pkin(10));
t354 = cos(pkin(6));
t337 = t354 * t353;
t306 = t317 * t325 + t322 * t337;
t318 = sin(pkin(6));
t343 = t318 * t353;
t357 = t306 * t315 - t314 * t343;
t344 = t317 * t354;
t308 = -t322 * t344 + t353 * t325;
t351 = t317 * t318;
t350 = t318 * t322;
t349 = t318 * t325;
t348 = qJD(2) * t322;
t346 = qJD(2) * t349;
t345 = t318 * t348;
t335 = t325 * t337;
t334 = -t308 * t314 + t315 * t351;
t333 = t308 * t315 + t314 * t351;
t332 = pkin(5) + qJ(4) + pkin(8) + t339;
t299 = t314 * t350 - t354 * t315;
t331 = t354 * t314 + t315 * t350;
t330 = -t306 * t314 - t315 * t343;
t329 = t338 * qJD(6) + qJD(4);
t307 = t353 * t322 + t325 * t344;
t324 = cos(qJ(3));
t327 = -t324 * pkin(3) - t336 * t314 - t347 * t315 - pkin(2);
t305 = t317 * t322 - t335;
t304 = t308 * qJD(2);
t303 = t307 * qJD(2);
t302 = t306 * qJD(2);
t301 = -qJD(2) * t335 + t317 * t348;
t293 = t331 * qJD(3) + t314 * t346;
t291 = t333 * qJD(3) - t303 * t314;
t289 = t357 * qJD(3) - t301 * t314;
t1 = [0, -t332 * t303 + t327 * t304 + t326 * t307 + t329 * t308, t328 * t333 + t336 * (t334 * qJD(3) - t303 * t315) - t347 * t291 + (t303 * t321 + (-t308 * t324 - t321 * t351) * qJD(3)) * pkin(3), t304, t291 (t291 * t323 - t304 * t320) * r_i_i_C(1) + (-t291 * t320 - t304 * t323) * r_i_i_C(2) + ((-t307 * t323 + t320 * t334) * r_i_i_C(1) + (t307 * t320 + t323 * t334) * r_i_i_C(2)) * qJD(6); 0, -t332 * t301 + t327 * t302 + t326 * t305 + t329 * t306, t328 * t357 + t336 * (t330 * qJD(3) - t301 * t315) - t347 * t289 + (t301 * t321 + (-t306 * t324 + t321 * t343) * qJD(3)) * pkin(3), t302, t289 (t289 * t323 - t302 * t320) * r_i_i_C(1) + (-t289 * t320 - t302 * t323) * r_i_i_C(2) + ((-t305 * t323 + t320 * t330) * r_i_i_C(1) + (t305 * t320 + t323 * t330) * r_i_i_C(2)) * qJD(6); 0 ((t327 * qJD(2) + t329) * t322 + (t332 * qJD(2) - t326) * t325) * t318, t328 * t331 + t336 * (-t299 * qJD(3) + t315 * t346) - t347 * t293 + (-t321 * t346 + (-t354 * t321 - t324 * t350) * qJD(3)) * pkin(3), t345, t293 (t293 * t323 - t320 * t345) * r_i_i_C(1) + (-t293 * t320 - t323 * t345) * r_i_i_C(2) + ((-t299 * t320 + t323 * t349) * r_i_i_C(1) + (-t299 * t323 - t320 * t349) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
