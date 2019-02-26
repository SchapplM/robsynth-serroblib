% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:16
% EndTime: 2019-02-26 19:48:16
% DurationCPUTime: 0.36s
% Computational Cost: add. (490->73), mult. (844->126), div. (0->0), fcn. (851->13), ass. (0->52)
t315 = pkin(11) + qJ(4);
t311 = sin(t315);
t313 = cos(t315);
t314 = pkin(12) + qJ(6);
t310 = sin(t314);
t312 = cos(t314);
t334 = t310 * r_i_i_C(1) + t312 * r_i_i_C(2);
t329 = qJD(6) * t334;
t335 = t312 * r_i_i_C(1) - t310 * r_i_i_C(2);
t332 = cos(pkin(12)) * pkin(5) + pkin(4) + t335;
t351 = r_i_i_C(3) + pkin(9) + qJ(5);
t323 = (t332 * t311 - t351 * t313) * qJD(4) - t311 * qJD(5) + t313 * t329;
t317 = sin(pkin(10));
t321 = sin(qJ(2));
t322 = cos(qJ(2));
t349 = cos(pkin(10));
t350 = cos(pkin(6));
t333 = t350 * t349;
t301 = t317 * t322 + t321 * t333;
t318 = sin(pkin(6));
t339 = t318 * t349;
t291 = t301 * t313 - t311 * t339;
t340 = t317 * t350;
t303 = -t321 * t340 + t349 * t322;
t347 = t317 * t318;
t346 = t318 * t321;
t345 = t318 * t322;
t344 = qJD(2) * t321;
t342 = qJD(2) * t345;
t341 = t318 * t344;
t331 = t322 * t333;
t330 = -t303 * t311 + t313 * t347;
t293 = t303 * t313 + t311 * t347;
t328 = -t301 * t311 - t313 * t339;
t327 = -t311 * t346 + t350 * t313;
t295 = t350 * t311 + t313 * t346;
t326 = sin(pkin(12)) * pkin(5) + pkin(8) + qJ(3) + t334;
t325 = t335 * qJD(6) + qJD(3);
t302 = t349 * t321 + t322 * t340;
t324 = -t351 * t311 - t332 * t313 - cos(pkin(11)) * pkin(3) - pkin(2);
t300 = t317 * t321 - t331;
t299 = t303 * qJD(2);
t298 = t302 * qJD(2);
t297 = t301 * qJD(2);
t296 = -qJD(2) * t331 + t317 * t344;
t289 = t327 * qJD(4) + t313 * t342;
t288 = t295 * qJD(4) + t311 * t342;
t287 = t330 * qJD(4) - t298 * t313;
t286 = t293 * qJD(4) - t298 * t311;
t285 = t328 * qJD(4) - t296 * t313;
t284 = t291 * qJD(4) - t296 * t311;
t1 = [0, -t326 * t298 + t324 * t299 + t323 * t302 + t325 * t303, t299, t293 * qJD(5) - t332 * t286 + t351 * t287 - t330 * t329, t286 (-t287 * t310 + t299 * t312) * r_i_i_C(1) + (-t287 * t312 - t299 * t310) * r_i_i_C(2) + ((-t293 * t312 - t302 * t310) * r_i_i_C(1) + (t293 * t310 - t302 * t312) * r_i_i_C(2)) * qJD(6); 0, -t326 * t296 + t324 * t297 + t323 * t300 + t325 * t301, t297, t291 * qJD(5) - t332 * t284 + t351 * t285 - t328 * t329, t284 (-t285 * t310 + t297 * t312) * r_i_i_C(1) + (-t285 * t312 - t297 * t310) * r_i_i_C(2) + ((-t291 * t312 - t300 * t310) * r_i_i_C(1) + (t291 * t310 - t300 * t312) * r_i_i_C(2)) * qJD(6); 0 ((t324 * qJD(2) + t325) * t321 + (t326 * qJD(2) - t323) * t322) * t318, t341, t295 * qJD(5) - t332 * t288 + t351 * t289 - t327 * t329, t288 (-t289 * t310 + t312 * t341) * r_i_i_C(1) + (-t289 * t312 - t310 * t341) * r_i_i_C(2) + ((-t295 * t312 + t310 * t345) * r_i_i_C(1) + (t295 * t310 + t312 * t345) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
