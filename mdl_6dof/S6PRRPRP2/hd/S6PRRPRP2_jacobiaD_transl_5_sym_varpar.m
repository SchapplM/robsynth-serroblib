% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP2
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

function JaD_transl = S6PRRPRP2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:49
% EndTime: 2019-02-26 20:01:49
% DurationCPUTime: 0.40s
% Computational Cost: add. (378->86), mult. (811->155), div. (0->0), fcn. (804->12), ass. (0->56)
t307 = qJ(3) + pkin(11);
t305 = sin(t307);
t306 = cos(t307);
t313 = sin(qJ(3));
t312 = sin(qJ(5));
t315 = cos(qJ(5));
t327 = t315 * r_i_i_C(1) - t312 * r_i_i_C(2);
t325 = pkin(4) + t327;
t349 = pkin(9) + r_i_i_C(3);
t351 = (t313 * pkin(3) + t325 * t305 - t349 * t306) * qJD(3);
t326 = t312 * r_i_i_C(1) + t315 * r_i_i_C(2);
t345 = cos(pkin(6));
t308 = sin(pkin(10));
t309 = sin(pkin(6));
t344 = t308 * t309;
t310 = cos(pkin(10));
t343 = t309 * t310;
t342 = t309 * t313;
t314 = sin(qJ(2));
t341 = t309 * t314;
t317 = cos(qJ(2));
t340 = t309 * t317;
t339 = qJD(2) * t314;
t338 = qJD(2) * t317;
t337 = qJD(5) * t306;
t336 = qJD(5) * t312;
t335 = qJD(5) * t315;
t334 = t309 * t338;
t333 = t309 * t339;
t332 = t314 * t345;
t331 = t317 * t345;
t329 = t308 * t332;
t328 = t310 * t331;
t298 = t308 * t317 + t310 * t332;
t324 = -t298 * t305 - t306 * t343;
t323 = -t298 * t306 + t305 * t343;
t300 = t310 * t317 - t329;
t322 = -t300 * t305 + t306 * t344;
t290 = t300 * t306 + t305 * t344;
t321 = qJD(5) * t326;
t299 = t308 * t331 + t310 * t314;
t320 = -t305 * t341 + t345 * t306;
t292 = t345 * t305 + t306 * t341;
t316 = cos(qJ(3));
t319 = -t316 * pkin(3) - t349 * t305 - t325 * t306 - pkin(2);
t318 = t326 * t337 + t351;
t311 = -qJ(4) - pkin(8);
t297 = t308 * t314 - t328;
t296 = -qJD(2) * t329 + t310 * t338;
t295 = t299 * qJD(2);
t294 = t298 * qJD(2);
t293 = -qJD(2) * t328 + t308 * t339;
t286 = t320 * qJD(3) + t306 * t334;
t284 = t322 * qJD(3) - t295 * t306;
t282 = t324 * qJD(3) - t293 * t306;
t1 = [0 (-t295 * t312 + t300 * t335) * r_i_i_C(1) + (-t295 * t315 - t300 * t336) * r_i_i_C(2) + t295 * t311 + t300 * qJD(4) + t319 * t296 + t318 * t299, t349 * t284 - t322 * t321 + t325 * (-t290 * qJD(3) + t295 * t305) + (t295 * t313 + (-t300 * t316 - t308 * t342) * qJD(3)) * pkin(3), t296 (-t284 * t312 + t296 * t315) * r_i_i_C(1) + (-t284 * t315 - t296 * t312) * r_i_i_C(2) + ((-t290 * t315 - t299 * t312) * r_i_i_C(1) + (t290 * t312 - t299 * t315) * r_i_i_C(2)) * qJD(5), 0; 0 (-t293 * t312 + t298 * t335) * r_i_i_C(1) + (-t293 * t315 - t298 * t336) * r_i_i_C(2) + t293 * t311 + t298 * qJD(4) + t319 * t294 + t318 * t297, t349 * t282 - t324 * t321 + t325 * (t323 * qJD(3) + t293 * t305) + (t293 * t313 + (-t298 * t316 + t310 * t342) * qJD(3)) * pkin(3), t294 (-t282 * t312 + t294 * t315) * r_i_i_C(1) + (-t282 * t315 - t294 * t312) * r_i_i_C(2) + ((-t297 * t312 + t315 * t323) * r_i_i_C(1) + (-t297 * t315 - t312 * t323) * r_i_i_C(2)) * qJD(5), 0; 0 ((t319 * qJD(2) + t327 * qJD(5) + qJD(4)) * t314 + (-qJD(2) * t311 - t351 + t326 * (qJD(2) - t337)) * t317) * t309, t349 * t286 - t320 * t321 + t325 * (-t292 * qJD(3) - t305 * t334) + (-t313 * t334 + (-t345 * t313 - t316 * t341) * qJD(3)) * pkin(3), t333 (-t286 * t312 + t315 * t333) * r_i_i_C(1) + (-t286 * t315 - t312 * t333) * r_i_i_C(2) + ((-t292 * t315 + t312 * t340) * r_i_i_C(1) + (t292 * t312 + t315 * t340) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
