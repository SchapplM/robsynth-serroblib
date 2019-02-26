% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:05
% EndTime: 2019-02-26 20:03:05
% DurationCPUTime: 0.48s
% Computational Cost: add. (382->82), mult. (1191->137), div. (0->0), fcn. (1192->10), ass. (0->60)
t313 = sin(qJ(3));
t316 = cos(qJ(3));
t312 = sin(qJ(5));
t315 = cos(qJ(5));
t358 = pkin(5) + r_i_i_C(1);
t331 = -t315 * r_i_i_C(2) - t312 * t358;
t326 = -qJ(4) + t331;
t346 = pkin(3) + r_i_i_C(3) + qJ(6) + pkin(9);
t363 = -(t313 * t346 + t326 * t316) * qJD(3) + t316 * qJD(6);
t357 = t312 * r_i_i_C(2);
t320 = qJD(4) + (t315 * t358 - t357) * qJD(5);
t309 = sin(pkin(10));
t314 = sin(qJ(2));
t317 = cos(qJ(2));
t354 = cos(pkin(10));
t355 = cos(pkin(6));
t337 = t355 * t354;
t299 = t309 * t317 + t314 * t337;
t310 = sin(pkin(6));
t342 = t310 * t354;
t360 = t299 * t316 - t313 * t342;
t343 = t309 * t355;
t301 = -t314 * t343 + t354 * t317;
t352 = t310 * t313;
t351 = t310 * t316;
t350 = t310 * t317;
t349 = qJD(2) * t314;
t348 = qJD(5) * t313;
t345 = t310 * t349;
t344 = qJD(2) * t350;
t340 = t315 * pkin(5) + pkin(4) + pkin(8) - t357;
t332 = t317 * t337;
t294 = -qJD(2) * t332 + t309 * t349;
t284 = t360 * qJD(3) - t294 * t313;
t295 = t299 * qJD(2);
t336 = t284 * t315 - t295 * t312;
t300 = t314 * t354 + t317 * t343;
t296 = t300 * qJD(2);
t329 = t301 * t316 + t309 * t352;
t286 = qJD(3) * t329 - t296 * t313;
t297 = t301 * qJD(2);
t335 = t286 * t315 - t297 * t312;
t298 = t309 * t314 - t332;
t321 = -t299 * t313 - t316 * t342;
t334 = -t298 * t315 + t312 * t321;
t330 = -t301 * t313 + t309 * t351;
t333 = -t300 * t315 + t312 * t330;
t302 = t314 * t352 - t316 * t355;
t328 = -t302 * t312 + t315 * t350;
t327 = -t315 * r_i_i_C(1) - t340;
t322 = t313 * t355 + t314 * t351;
t292 = qJD(3) * t322 + t313 * t344;
t324 = t292 * t315 - t312 * t345;
t323 = t331 * qJD(5);
t319 = t313 * t326 - t316 * t346 - pkin(2);
t318 = -t320 * t313 - t363;
t293 = -qJD(3) * t302 + t316 * t344;
t287 = qJD(3) * t330 - t296 * t316;
t285 = qJD(3) * t321 - t294 * t316;
t1 = [0, t296 * t327 + t297 * t319 + t300 * t318 + t301 * t323, qJD(6) * t330 - t286 * t346 - t287 * t326 + t320 * t329, t286, t335 * r_i_i_C(1) + (-t286 * t312 - t297 * t315) * r_i_i_C(2) + (t333 * r_i_i_C(1) + (t300 * t312 + t315 * t330) * r_i_i_C(2)) * qJD(5) + (qJD(5) * t333 + t335) * pkin(5), t287; 0, t294 * t327 + t295 * t319 + t298 * t318 + t299 * t323, qJD(6) * t321 - t346 * t284 - t285 * t326 + t320 * t360, t284, t336 * r_i_i_C(1) + (-t284 * t312 - t295 * t315) * r_i_i_C(2) + (t334 * r_i_i_C(1) + (t298 * t312 + t315 * t321) * r_i_i_C(2)) * qJD(5) + (qJD(5) * t334 + t336) * pkin(5), t285; 0 ((qJD(2) * t319 + t323) * t314 + ((-qJD(5) * t357 + qJD(4)) * t313 + t340 * qJD(2) + ((qJD(2) + t348) * r_i_i_C(1) + pkin(5) * t348) * t315 + t363) * t317) * t310, -t302 * qJD(6) - t292 * t346 - t293 * t326 + t320 * t322, t292, t324 * r_i_i_C(1) + (-t292 * t312 - t315 * t345) * r_i_i_C(2) + (t328 * r_i_i_C(1) + (-t302 * t315 - t312 * t350) * r_i_i_C(2)) * qJD(5) + (qJD(5) * t328 + t324) * pkin(5), t293;];
JaD_transl  = t1;
