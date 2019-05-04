% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14V3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14V3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:06
% EndTime: 2019-04-12 15:12:06
% DurationCPUTime: 0.32s
% Computational Cost: add. (197->79), mult. (636->146), div. (0->0), fcn. (595->8), ass. (0->61)
t311 = cos(qJ(4));
t361 = -qJD(5) * t311 + qJD(2);
t308 = sin(qJ(2));
t309 = sin(qJ(1));
t312 = cos(qJ(2));
t332 = qJD(1) * t312 - qJD(4);
t313 = cos(qJ(1));
t344 = qJD(2) * t313;
t360 = t308 * t344 + t332 * t309;
t306 = sin(qJ(5));
t310 = cos(qJ(5));
t307 = sin(qJ(4));
t343 = qJD(4) * t307;
t318 = t306 * t343 + t361 * t310;
t319 = -t361 * t306 + t310 * t343;
t342 = qJD(4) * t311;
t359 = t319 * r_i_i_C(1) - t318 * r_i_i_C(2) - r_i_i_C(3) * t342 - qJD(2) * qJ(3);
t354 = t308 * t311;
t355 = t307 * r_i_i_C(3);
t358 = (-t306 * t312 + t310 * t354) * r_i_i_C(1) - (t306 * t354 + t310 * t312) * r_i_i_C(2) - t312 * qJ(3) + t308 * t355;
t357 = t306 * r_i_i_C(1);
t356 = t306 * r_i_i_C(2);
t353 = t308 * t313;
t352 = t309 * t311;
t351 = t309 * t312;
t350 = t313 * t307;
t349 = t313 * t311;
t348 = qJD(1) * t309;
t347 = qJD(1) * t313;
t346 = qJD(2) * t308;
t345 = qJD(2) * t312;
t341 = qJD(5) * t306;
t340 = qJD(5) * t309;
t339 = qJD(5) * t310;
t337 = t307 * t351;
t336 = t309 * t346;
t335 = t309 * t345;
t333 = -qJD(4) * t312 + qJD(1);
t331 = qJD(2) * t311 - qJD(5);
t330 = -t310 * r_i_i_C(1) + t356;
t327 = t333 * t313;
t295 = t307 * t327 - t360 * t311;
t329 = qJD(5) * t353 + t295;
t301 = t309 * t307 + t312 * t349;
t297 = t301 * qJD(1) - qJD(4) * t337 - t311 * t336 - t313 * t342;
t328 = -t308 * t340 - t297;
t326 = t331 * t310;
t321 = t308 * t347 + t335;
t320 = -t308 * t348 + t312 * t344;
t299 = t311 * t351 - t350;
t317 = -qJD(5) * t299 + t321;
t316 = -qJD(5) * t301 + t320;
t315 = -r_i_i_C(1) * t326 - qJD(2) * t355 + t331 * t356 + qJD(3);
t314 = t359 * t308 + t315 * t312;
t300 = -t312 * t350 + t352;
t298 = -t337 - t349;
t296 = t333 * t352 + (-t332 * t313 + t336) * t307;
t294 = t360 * t307 + t311 * t327;
t293 = t316 * t306 + t329 * t310;
t292 = -t329 * t306 + t316 * t310;
t1 = [(-t297 * t310 + t299 * t341 - t306 * t335) * r_i_i_C(1) + (t297 * t306 + t299 * t339 - t310 * t335) * r_i_i_C(2) + t296 * r_i_i_C(3) - qJ(3) * t335 + ((-t306 * t347 - t309 * t339) * r_i_i_C(1) + (t306 * t340 - t310 * t347) * r_i_i_C(2) - qJ(3) * t347 - t309 * qJD(3)) * t308, t314 * t313 + t358 * t348, t320, t295 * r_i_i_C(3) + (-t294 * t306 - t300 * t339) * r_i_i_C(2) + (t294 * t310 - t300 * t341) * r_i_i_C(1), t292 * r_i_i_C(1) - t293 * r_i_i_C(2), 0; t293 * r_i_i_C(1) + t292 * r_i_i_C(2) - t294 * r_i_i_C(3) + t320 * qJ(3) + qJD(3) * t353, t314 * t309 - t358 * t347, t321, t297 * r_i_i_C(3) + (-t296 * t306 - t298 * t339) * r_i_i_C(2) + (t296 * t310 - t298 * t341) * r_i_i_C(1) (t317 * r_i_i_C(1) + t328 * r_i_i_C(2)) * t310 + (t328 * r_i_i_C(1) - t317 * r_i_i_C(2)) * t306, 0; 0, t315 * t308 - t359 * t312, t346 (t330 * t308 * qJD(4) + r_i_i_C(3) * t345) * t311 + (t330 * t345 + (-qJD(4) * r_i_i_C(3) + (t310 * r_i_i_C(2) + t357) * qJD(5)) * t308) * t307 (-r_i_i_C(2) * t326 - t331 * t357) * t312 + (t318 * r_i_i_C(1) + t319 * r_i_i_C(2)) * t308, 0;];
JaD_transl  = t1;
