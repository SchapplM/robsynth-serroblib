% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:14
% EndTime: 2019-02-26 19:58:15
% DurationCPUTime: 0.40s
% Computational Cost: add. (508->83), mult. (905->142), div. (0->0), fcn. (904->14), ass. (0->54)
t318 = qJ(3) + pkin(11);
t314 = sin(t318);
t316 = cos(t318);
t324 = sin(qJ(3));
t317 = pkin(12) + qJ(6);
t313 = sin(t317);
t315 = cos(t317);
t339 = t313 * r_i_i_C(1) + t315 * r_i_i_C(2);
t334 = qJD(6) * t339;
t340 = t315 * r_i_i_C(1) - t313 * r_i_i_C(2);
t337 = cos(pkin(12)) * pkin(5) + pkin(4) + t340;
t356 = r_i_i_C(3) + pkin(9) + qJ(5);
t328 = (t324 * pkin(3) + t337 * t314 - t356 * t316) * qJD(3) - t314 * qJD(5) + t316 * t334;
t320 = sin(pkin(10));
t325 = sin(qJ(2));
t327 = cos(qJ(2));
t354 = cos(pkin(10));
t355 = cos(pkin(6));
t338 = t355 * t354;
t304 = t320 * t327 + t325 * t338;
t321 = sin(pkin(6));
t344 = t321 * t354;
t294 = t304 * t316 - t314 * t344;
t345 = t320 * t355;
t306 = -t325 * t345 + t354 * t327;
t352 = t320 * t321;
t351 = t321 * t325;
t350 = t321 * t327;
t349 = qJD(2) * t325;
t347 = qJD(2) * t350;
t346 = t321 * t349;
t336 = t327 * t338;
t335 = -t306 * t314 + t316 * t352;
t296 = t306 * t316 + t314 * t352;
t333 = -t304 * t314 - t316 * t344;
t332 = -t314 * t351 + t355 * t316;
t298 = t355 * t314 + t316 * t351;
t331 = sin(pkin(12)) * pkin(5) + qJ(4) + pkin(8) + t339;
t330 = t340 * qJD(6) + qJD(4);
t305 = t354 * t325 + t327 * t345;
t326 = cos(qJ(3));
t329 = -t326 * pkin(3) - t356 * t314 - t337 * t316 - pkin(2);
t303 = t320 * t325 - t336;
t302 = t306 * qJD(2);
t301 = t305 * qJD(2);
t300 = t304 * qJD(2);
t299 = -qJD(2) * t336 + t320 * t349;
t292 = t332 * qJD(3) + t316 * t347;
t291 = t298 * qJD(3) + t314 * t347;
t290 = t335 * qJD(3) - t301 * t316;
t289 = t296 * qJD(3) - t301 * t314;
t288 = t333 * qJD(3) - t299 * t316;
t287 = t294 * qJD(3) - t299 * t314;
t1 = [0, -t331 * t301 + t329 * t302 + t328 * t305 + t330 * t306, t296 * qJD(5) + t356 * t290 - t335 * t334 - t337 * t289 + (t301 * t324 + (-t306 * t326 - t324 * t352) * qJD(3)) * pkin(3), t302, t289 (-t290 * t313 + t302 * t315) * r_i_i_C(1) + (-t290 * t315 - t302 * t313) * r_i_i_C(2) + ((-t296 * t315 - t305 * t313) * r_i_i_C(1) + (t296 * t313 - t305 * t315) * r_i_i_C(2)) * qJD(6); 0, -t331 * t299 + t329 * t300 + t328 * t303 + t330 * t304, t294 * qJD(5) + t356 * t288 - t333 * t334 - t337 * t287 + (t299 * t324 + (-t304 * t326 + t324 * t344) * qJD(3)) * pkin(3), t300, t287 (-t288 * t313 + t300 * t315) * r_i_i_C(1) + (-t288 * t315 - t300 * t313) * r_i_i_C(2) + ((-t294 * t315 - t303 * t313) * r_i_i_C(1) + (t294 * t313 - t303 * t315) * r_i_i_C(2)) * qJD(6); 0 ((t329 * qJD(2) + t330) * t325 + (t331 * qJD(2) - t328) * t327) * t321, t298 * qJD(5) + t356 * t292 - t332 * t334 - t337 * t291 + (-t324 * t347 + (-t355 * t324 - t326 * t351) * qJD(3)) * pkin(3), t346, t291 (-t292 * t313 + t315 * t346) * r_i_i_C(1) + (-t292 * t315 - t313 * t346) * r_i_i_C(2) + ((-t298 * t315 + t313 * t350) * r_i_i_C(1) + (t298 * t313 + t315 * t350) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
