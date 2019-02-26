% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR10_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR10_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:14
% EndTime: 2019-02-26 22:21:14
% DurationCPUTime: 0.49s
% Computational Cost: add. (376->90), mult. (1158->144), div. (0->0), fcn. (1082->8), ass. (0->57)
t285 = sin(qJ(2));
t335 = (qJD(3) - qJD(5)) * t285;
t284 = sin(qJ(3));
t288 = cos(qJ(3));
t290 = cos(qJ(1));
t322 = t288 * t290;
t286 = sin(qJ(1));
t289 = cos(qJ(2));
t323 = t286 * t289;
t267 = t284 * t323 + t322;
t321 = t290 * t284;
t268 = t288 * t323 - t321;
t283 = sin(qJ(5));
t287 = cos(qJ(5));
t302 = t267 * t283 + t268 * t287;
t303 = t267 * t287 - t268 * t283;
t334 = (t302 * r_i_i_C(1) + t303 * r_i_i_C(2)) * qJD(5);
t326 = -pkin(3) - pkin(4);
t305 = t283 * r_i_i_C(2) + t326;
t296 = t287 * r_i_i_C(1) - t305;
t306 = -t283 * r_i_i_C(1) - qJ(4);
t297 = t287 * r_i_i_C(2) - t306;
t298 = t283 * t284 + t287 * t288;
t299 = t283 * t288 - t284 * t287;
t312 = -r_i_i_C(3) - pkin(9) + pkin(8);
t291 = -t284 * qJD(4) - t312 * qJD(2) + (t299 * r_i_i_C(1) + t298 * r_i_i_C(2)) * qJD(5) + (t296 * t284 - t297 * t288) * qJD(3);
t332 = -t285 * pkin(2) + t312 * t289;
t327 = t297 * t284 + t296 * t288 + pkin(2);
t315 = qJD(3) * t288;
t319 = qJD(1) * t290;
t294 = t284 * t319 + t286 * t315;
t314 = qJD(3) * t290;
t308 = t284 * t314;
t318 = qJD(2) * t286;
t311 = t285 * t318;
t320 = qJD(1) * t286;
t265 = -t284 * t311 - t288 * t320 + t294 * t289 - t308;
t262 = t265 * t287;
t324 = t286 * t284;
t317 = qJD(2) * t289;
t316 = qJD(2) * t290;
t310 = qJD(3) * t324;
t309 = t285 * t316;
t307 = t288 * t314;
t304 = -r_i_i_C(1) * (-t298 * t335 + t299 * t317) - (t298 * t317 + t299 * t335) * r_i_i_C(2);
t269 = -t286 * t288 + t289 * t321;
t270 = t289 * t322 + t324;
t301 = t269 * t287 - t270 * t283;
t300 = t269 * t283 + t270 * t287;
t295 = -pkin(2) * t289 - t312 * t285 - pkin(1);
t292 = qJD(2) * t327;
t266 = t270 * qJD(1) - t288 * t311 - t289 * t310 - t307;
t264 = t289 * t308 + (t289 * t320 + t309) * t288 - t294;
t263 = t267 * qJD(1) + t284 * t309 - t289 * t307 - t310;
t259 = t301 * qJD(5) - t263 * t283 - t264 * t287;
t258 = -t300 * qJD(5) - t263 * t287 + t264 * t283;
t1 = [-t262 * r_i_i_C(2) - t267 * qJD(4) + t306 * t265 - t296 * t266 + (-t303 * r_i_i_C(1) + t302 * r_i_i_C(2)) * qJD(5) - t332 * t318 + (-t286 * pkin(7) + t295 * t290) * qJD(1) (-t290 * t292 - t312 * t320) * t289 + (t291 * t290 + t327 * t320) * t285, t270 * qJD(4) - t297 * t264 + t296 * t263 + (t300 * r_i_i_C(1) + t301 * r_i_i_C(2)) * qJD(5), -t263, r_i_i_C(1) * t258 - r_i_i_C(2) * t259, 0; t259 * r_i_i_C(1) + t258 * r_i_i_C(2) - t263 * qJ(4) + t269 * qJD(4) + t326 * t264 + t332 * t316 + (pkin(7) * t290 + t295 * t286) * qJD(1) (-t286 * t292 + t312 * t319) * t289 + (t291 * t286 - t319 * t327) * t285, -t262 * r_i_i_C(1) + t268 * qJD(4) + t305 * t265 + t297 * t266 + t334, t265 (-t266 * t283 + t262) * r_i_i_C(1) + (-t265 * t283 - t266 * t287) * r_i_i_C(2) - t334, 0; 0, -t285 * t292 - t291 * t289 (t288 * qJ(4) + t326 * t284) * t317 + (qJD(4) * t288 + (-t284 * qJ(4) + t326 * t288) * qJD(3)) * t285 - t304, t284 * t317 + t285 * t315, t304, 0;];
JaD_transl  = t1;
