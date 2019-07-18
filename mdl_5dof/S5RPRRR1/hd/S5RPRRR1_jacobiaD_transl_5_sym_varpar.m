% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:27
% EndTime: 2019-07-18 13:26:27
% DurationCPUTime: 0.35s
% Computational Cost: add. (185->67), mult. (600->128), div. (0->0), fcn. (566->8), ass. (0->57)
t311 = cos(qJ(4));
t358 = -qJD(5) * t311 + qJD(3);
t308 = sin(qJ(3));
t309 = sin(qJ(1));
t312 = cos(qJ(3));
t329 = qJD(1) * t312 - qJD(4);
t313 = cos(qJ(1));
t342 = qJD(3) * t313;
t357 = t308 * t342 + t329 * t309;
t343 = qJD(3) * t312;
t344 = qJD(1) * t313;
t307 = sin(qJ(4));
t347 = t313 * t307;
t348 = t309 * t312;
t318 = -qJD(5) * (t311 * t348 - t347) + t308 * t344 + t309 * t343;
t346 = t313 * t311;
t301 = t309 * t307 + t312 * t346;
t332 = qJD(3) * t308 * t309;
t333 = t307 * t348;
t340 = qJD(4) * t311;
t297 = t301 * qJD(1) - qJD(4) * t333 - t311 * t332 - t313 * t340;
t338 = qJD(5) * t308;
t325 = t309 * t338 + t297;
t356 = t318 * r_i_i_C(1) - t325 * r_i_i_C(2);
t306 = sin(qJ(5));
t310 = cos(qJ(5));
t341 = qJD(4) * t307;
t319 = t306 * t341 + t358 * t310;
t320 = -t358 * t306 + t310 * t341;
t355 = t320 * r_i_i_C(1) - t319 * r_i_i_C(2) - r_i_i_C(3) * t340;
t350 = t308 * t311;
t351 = t307 * r_i_i_C(3);
t354 = (-t306 * t312 + t310 * t350) * r_i_i_C(1) - (t306 * t350 + t310 * t312) * r_i_i_C(2) + t308 * t351;
t353 = t306 * r_i_i_C(1);
t352 = t306 * r_i_i_C(2);
t349 = t309 * t311;
t345 = qJD(1) * t309;
t339 = qJD(5) * t306;
t337 = qJD(5) * t310;
t330 = -qJD(4) * t312 + qJD(1);
t328 = qJD(3) * t311 - qJD(5);
t327 = -t310 * r_i_i_C(1) + t352;
t324 = t330 * t313;
t295 = t307 * t324 - t357 * t311;
t326 = t313 * t338 + t295;
t323 = t328 * t310;
t317 = -qJD(5) * t301 - t308 * t345 + t312 * t342;
t316 = -r_i_i_C(1) * t323 - qJD(3) * t351 + t328 * t352;
t315 = -t325 * r_i_i_C(1) - t318 * r_i_i_C(2);
t314 = t355 * t308 + t316 * t312;
t300 = -t312 * t347 + t349;
t298 = -t333 - t346;
t296 = t330 * t349 + (-t329 * t313 + t332) * t307;
t294 = t357 * t307 + t311 * t324;
t293 = t317 * t306 + t326 * t310;
t292 = -t326 * t306 + t317 * t310;
t1 = [t296 * r_i_i_C(3) - qJ(2) * t345 + t313 * qJD(2) - t356 * t306 + t315 * t310, t344, t314 * t313 + t354 * t345, t295 * r_i_i_C(3) + (-t294 * t306 - t300 * t337) * r_i_i_C(2) + (t294 * t310 - t300 * t339) * r_i_i_C(1), t292 * r_i_i_C(1) - t293 * r_i_i_C(2); t293 * r_i_i_C(1) + t292 * r_i_i_C(2) - t294 * r_i_i_C(3) + qJ(2) * t344 + t309 * qJD(2), t345, t314 * t309 - t354 * t344, t297 * r_i_i_C(3) + (-t296 * t306 - t298 * t337) * r_i_i_C(2) + (t296 * t310 - t298 * t339) * r_i_i_C(1), t315 * t306 + t356 * t310; 0, 0, t316 * t308 - t355 * t312, (t327 * t308 * qJD(4) + r_i_i_C(3) * t343) * t311 + (t327 * t343 + (-qJD(4) * r_i_i_C(3) + (t310 * r_i_i_C(2) + t353) * qJD(5)) * t308) * t307, (-r_i_i_C(2) * t323 - t328 * t353) * t312 + (t319 * r_i_i_C(1) + t320 * r_i_i_C(2)) * t308;];
JaD_transl  = t1;
