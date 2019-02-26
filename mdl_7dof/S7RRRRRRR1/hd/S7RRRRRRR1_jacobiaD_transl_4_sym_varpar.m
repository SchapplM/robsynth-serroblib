% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JaD_transl [3x7]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S7RRRRRRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:22
% EndTime: 2019-02-26 22:54:22
% DurationCPUTime: 0.36s
% Computational Cost: add. (190->67), mult. (618->126), div. (0->0), fcn. (577->8), ass. (0->57)
t301 = cos(qJ(3));
t348 = -qJD(4) * t301 + qJD(2);
t297 = sin(qJ(3));
t303 = cos(qJ(1));
t338 = t303 * t297;
t299 = sin(qJ(1));
t302 = cos(qJ(2));
t339 = t299 * t302;
t289 = t301 * t339 + t338;
t298 = sin(qJ(2));
t334 = qJD(2) * t302;
t335 = qJD(1) * t303;
t313 = t298 * t335 + t299 * t334;
t308 = -qJD(4) * t289 + t313;
t324 = qJD(1) * t302 + qJD(3);
t326 = t299 * qJD(2) * t298;
t327 = t297 * t339;
t336 = qJD(1) * t299;
t337 = t303 * t301;
t287 = -qJD(3) * t327 - t297 * t336 - t301 * t326 + t324 * t337;
t330 = qJD(4) * t298;
t320 = t299 * t330 + t287;
t347 = t308 * r_i_i_C(1) - t320 * r_i_i_C(2);
t296 = sin(qJ(4));
t300 = cos(qJ(4));
t332 = qJD(3) * t297;
t309 = t296 * t332 + t348 * t300;
t310 = -t348 * t296 + t300 * t332;
t341 = qJD(3) * r_i_i_C(3);
t346 = qJD(2) * pkin(2) + t310 * r_i_i_C(1) - t309 * r_i_i_C(2) + t301 * t341;
t340 = t298 * t301;
t342 = t297 * r_i_i_C(3);
t345 = (-t296 * t302 + t300 * t340) * r_i_i_C(1) - (t296 * t340 + t300 * t302) * r_i_i_C(2) + t302 * pkin(2) - t298 * t342;
t344 = t296 * r_i_i_C(1);
t343 = t296 * r_i_i_C(2);
t333 = qJD(2) * t303;
t331 = qJD(4) * t296;
t329 = qJD(4) * t300;
t325 = qJD(3) * t302 + qJD(1);
t323 = qJD(2) * t301 - qJD(4);
t322 = -t300 * r_i_i_C(1) + t343;
t311 = t298 * t333 + t324 * t299;
t285 = t311 * t301 + t325 * t338;
t321 = t303 * t330 - t285;
t319 = t323 * t300;
t314 = t299 * t301 + t302 * t338;
t312 = t298 * t336 - t302 * t333;
t307 = -qJD(4) * (-t299 * t297 + t302 * t337) - t312;
t306 = -r_i_i_C(1) * t319 + qJD(2) * t342 + t323 * t343;
t305 = -t320 * r_i_i_C(1) - t308 * r_i_i_C(2);
t304 = t346 * t298 + t306 * t302;
t288 = -t327 + t337;
t286 = t314 * qJD(1) + t289 * qJD(3) - t297 * t326;
t284 = t311 * t297 - t325 * t337;
t283 = t307 * t296 + t321 * t300;
t282 = -t321 * t296 + t307 * t300;
t1 = [t313 * pkin(2) + t286 * r_i_i_C(3) - t347 * t296 + t305 * t300, t304 * t303 + t345 * t336, t285 * r_i_i_C(3) + (-t284 * t296 + t314 * t329) * r_i_i_C(2) + (t284 * t300 + t314 * t331) * r_i_i_C(1), t282 * r_i_i_C(1) - t283 * r_i_i_C(2), 0, 0, 0; t312 * pkin(2) + t283 * r_i_i_C(1) + t282 * r_i_i_C(2) + t284 * r_i_i_C(3), t304 * t299 - t345 * t335, -t287 * r_i_i_C(3) + (t286 * t296 - t288 * t329) * r_i_i_C(2) + (-t286 * t300 - t288 * t331) * r_i_i_C(1), t305 * t296 + t347 * t300, 0, 0, 0; 0, t306 * t298 - t346 * t302 (t322 * t298 * qJD(3) - r_i_i_C(3) * t334) * t301 + (t322 * t334 + (t341 + (t300 * r_i_i_C(2) + t344) * qJD(4)) * t298) * t297 (-r_i_i_C(2) * t319 - t323 * t344) * t302 + (t309 * r_i_i_C(1) + t310 * r_i_i_C(2)) * t298, 0, 0, 0;];
JaD_transl  = t1;
