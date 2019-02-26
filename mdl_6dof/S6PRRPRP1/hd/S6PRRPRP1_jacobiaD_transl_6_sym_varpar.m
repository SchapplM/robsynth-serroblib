% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:12
% EndTime: 2019-02-26 20:01:12
% DurationCPUTime: 0.49s
% Computational Cost: add. (515->96), mult. (1064->162), div. (0->0), fcn. (1060->12), ass. (0->63)
t321 = sin(qJ(5));
t324 = cos(qJ(5));
t367 = pkin(5) + r_i_i_C(1);
t371 = -t321 * r_i_i_C(2) + t367 * t324;
t369 = t324 * r_i_i_C(2) + t367 * t321;
t316 = qJ(3) + pkin(11);
t314 = sin(t316);
t315 = cos(t316);
t322 = sin(qJ(3));
t338 = pkin(4) + t371;
t363 = r_i_i_C(3) + qJ(6) + pkin(9);
t370 = -(t322 * pkin(3) + t338 * t314 - t363 * t315) * qJD(3) + t314 * qJD(6);
t317 = sin(pkin(10));
t323 = sin(qJ(2));
t326 = cos(qJ(2));
t361 = cos(pkin(10));
t362 = cos(pkin(6));
t343 = t362 * t361;
t305 = t317 * t326 + t323 * t343;
t318 = sin(pkin(6));
t348 = t318 * t361;
t295 = t305 * t315 - t314 * t348;
t349 = t317 * t362;
t307 = -t323 * t349 + t361 * t326;
t359 = t317 * t318;
t358 = t318 * t323;
t357 = t318 * t326;
t356 = qJD(2) * t323;
t355 = qJD(5) * t315;
t354 = qJD(5) * t321;
t353 = qJD(5) * t324;
t351 = qJD(2) * t357;
t350 = t318 * t356;
t337 = t326 * t343;
t300 = -qJD(2) * t337 + t317 * t356;
t330 = -t305 * t314 - t315 * t348;
t289 = t330 * qJD(3) - t300 * t315;
t301 = t305 * qJD(2);
t342 = -t289 * t321 + t301 * t324;
t306 = t361 * t323 + t326 * t349;
t302 = t306 * qJD(2);
t335 = -t307 * t314 + t315 * t359;
t291 = t335 * qJD(3) - t302 * t315;
t303 = t307 * qJD(2);
t341 = -t291 * t321 + t303 * t324;
t304 = t317 * t323 - t337;
t340 = -t295 * t324 - t304 * t321;
t297 = t307 * t315 + t314 * t359;
t339 = -t297 * t324 - t306 * t321;
t299 = t362 * t314 + t315 * t358;
t336 = -t299 * t324 + t321 * t357;
t329 = -t314 * t358 + t362 * t315;
t293 = t329 * qJD(3) + t315 * t351;
t332 = -t293 * t321 + t324 * t350;
t331 = qJD(5) * t369;
t325 = cos(qJ(3));
t328 = -t325 * pkin(3) - t363 * t314 - t338 * t315 - pkin(2);
t327 = t369 * t355 - t370;
t320 = -qJ(4) - pkin(8);
t292 = t299 * qJD(3) + t314 * t351;
t290 = t297 * qJD(3) - t302 * t314;
t288 = t295 * qJD(3) - t300 * t314;
t1 = [0 (-t302 * t324 - t307 * t354) * r_i_i_C(2) + t302 * t320 + t307 * qJD(4) + t328 * t303 + t327 * t306 + t367 * (-t302 * t321 + t307 * t353) t297 * qJD(6) + t363 * t291 - t338 * t290 - t335 * t331 + (t302 * t322 + (-t307 * t325 - t322 * t359) * qJD(3)) * pkin(3), t303, t341 * r_i_i_C(1) + (-t291 * t324 - t303 * t321) * r_i_i_C(2) + (t339 * r_i_i_C(1) + (t297 * t321 - t306 * t324) * r_i_i_C(2)) * qJD(5) + (t339 * qJD(5) + t341) * pkin(5), t290; 0 (-t300 * t324 - t305 * t354) * r_i_i_C(2) + t300 * t320 + t305 * qJD(4) + t328 * t301 + t327 * t304 + t367 * (-t300 * t321 + t305 * t353) t295 * qJD(6) + t363 * t289 - t338 * t288 - t330 * t331 + (t300 * t322 + (-t305 * t325 + t322 * t348) * qJD(3)) * pkin(3), t301, t342 * r_i_i_C(1) + (-t289 * t324 - t301 * t321) * r_i_i_C(2) + (t340 * r_i_i_C(1) + (t295 * t321 - t304 * t324) * r_i_i_C(2)) * qJD(5) + (t340 * qJD(5) + t342) * pkin(5), t288; 0 ((t328 * qJD(2) + t371 * qJD(5) + qJD(4)) * t323 + (-qJD(2) * t320 + t369 * (qJD(2) - t355) + t370) * t326) * t318, t299 * qJD(6) + t363 * t293 - t338 * t292 - t329 * t331 + (-t322 * t351 + (-t362 * t322 - t325 * t358) * qJD(3)) * pkin(3), t350, t332 * r_i_i_C(1) + (-t293 * t324 - t321 * t350) * r_i_i_C(2) + (t336 * r_i_i_C(1) + (t299 * t321 + t324 * t357) * r_i_i_C(2)) * qJD(5) + (t336 * qJD(5) + t332) * pkin(5), t292;];
JaD_transl  = t1;
