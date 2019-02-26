% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:35
% EndTime: 2019-02-26 19:51:35
% DurationCPUTime: 0.46s
% Computational Cost: add. (497->86), mult. (1003->146), div. (0->0), fcn. (1007->11), ass. (0->61)
t318 = sin(qJ(5));
t320 = cos(qJ(5));
t361 = pkin(5) + r_i_i_C(1);
t365 = -t318 * r_i_i_C(2) + t361 * t320;
t363 = t320 * r_i_i_C(2) + t361 * t318;
t313 = pkin(11) + qJ(4);
t311 = sin(t313);
t312 = cos(t313);
t333 = pkin(4) + t365;
t358 = r_i_i_C(3) + qJ(6) + pkin(9);
t364 = -(t333 * t311 - t358 * t312) * qJD(4) + t311 * qJD(6);
t314 = sin(pkin(10));
t319 = sin(qJ(2));
t321 = cos(qJ(2));
t356 = cos(pkin(10));
t357 = cos(pkin(6));
t338 = t357 * t356;
t302 = t314 * t321 + t319 * t338;
t315 = sin(pkin(6));
t343 = t315 * t356;
t292 = t302 * t312 - t311 * t343;
t344 = t314 * t357;
t304 = -t319 * t344 + t356 * t321;
t354 = t314 * t315;
t353 = t315 * t319;
t352 = t315 * t321;
t351 = qJD(2) * t319;
t350 = qJD(5) * t312;
t349 = qJD(5) * t318;
t348 = qJD(5) * t320;
t346 = qJD(2) * t352;
t345 = t315 * t351;
t332 = t321 * t338;
t297 = -qJD(2) * t332 + t314 * t351;
t325 = -t302 * t311 - t312 * t343;
t286 = t325 * qJD(4) - t297 * t312;
t298 = t302 * qJD(2);
t337 = -t286 * t318 + t298 * t320;
t303 = t356 * t319 + t321 * t344;
t299 = t303 * qJD(2);
t330 = -t304 * t311 + t312 * t354;
t288 = t330 * qJD(4) - t299 * t312;
t300 = t304 * qJD(2);
t336 = -t288 * t318 + t300 * t320;
t301 = t314 * t319 - t332;
t335 = -t292 * t320 - t301 * t318;
t294 = t304 * t312 + t311 * t354;
t334 = -t294 * t320 - t303 * t318;
t296 = t357 * t311 + t312 * t353;
t331 = -t296 * t320 + t318 * t352;
t324 = -t311 * t353 + t357 * t312;
t290 = t324 * qJD(4) + t312 * t346;
t327 = -t290 * t318 + t320 * t345;
t326 = qJD(5) * t363;
t323 = -t358 * t311 - t333 * t312 - cos(pkin(11)) * pkin(3) - pkin(2);
t322 = t363 * t350 - t364;
t317 = -pkin(8) - qJ(3);
t289 = t296 * qJD(4) + t311 * t346;
t287 = t294 * qJD(4) - t299 * t311;
t285 = t292 * qJD(4) - t297 * t311;
t1 = [0 (-t299 * t320 - t304 * t349) * r_i_i_C(2) + t299 * t317 + t304 * qJD(3) + t323 * t300 + t322 * t303 + t361 * (-t299 * t318 + t304 * t348) t300, t294 * qJD(6) - t333 * t287 + t358 * t288 - t330 * t326, t336 * r_i_i_C(1) + (-t288 * t320 - t300 * t318) * r_i_i_C(2) + (t334 * r_i_i_C(1) + (t294 * t318 - t303 * t320) * r_i_i_C(2)) * qJD(5) + (t334 * qJD(5) + t336) * pkin(5), t287; 0 (-t297 * t320 - t302 * t349) * r_i_i_C(2) + t297 * t317 + t302 * qJD(3) + t323 * t298 + t322 * t301 + t361 * (-t297 * t318 + t302 * t348) t298, t292 * qJD(6) - t333 * t285 + t358 * t286 - t325 * t326, t337 * r_i_i_C(1) + (-t286 * t320 - t298 * t318) * r_i_i_C(2) + (t335 * r_i_i_C(1) + (t292 * t318 - t301 * t320) * r_i_i_C(2)) * qJD(5) + (t335 * qJD(5) + t337) * pkin(5), t285; 0 ((t323 * qJD(2) + t365 * qJD(5) + qJD(3)) * t319 + (-qJD(2) * t317 + t363 * (qJD(2) - t350) + t364) * t321) * t315, t345, t296 * qJD(6) - t333 * t289 + t358 * t290 - t324 * t326, t327 * r_i_i_C(1) + (-t290 * t320 - t318 * t345) * r_i_i_C(2) + (t331 * r_i_i_C(1) + (t296 * t318 + t320 * t352) * r_i_i_C(2)) * qJD(5) + (t331 * qJD(5) + t327) * pkin(5), t289;];
JaD_transl  = t1;
