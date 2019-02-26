% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:05
% EndTime: 2019-02-26 21:06:05
% DurationCPUTime: 0.42s
% Computational Cost: add. (380->89), mult. (1166->140), div. (0->0), fcn. (1088->8), ass. (0->57)
t286 = sin(qJ(3));
t285 = sin(qJ(4));
t291 = cos(qJ(1));
t323 = t291 * t285;
t287 = sin(qJ(1));
t289 = cos(qJ(4));
t324 = t287 * t289;
t276 = t286 * t323 + t324;
t322 = t291 * t289;
t313 = t286 * t322;
t325 = t287 * t285;
t277 = t313 - t325;
t284 = sin(qJ(6));
t288 = cos(qJ(6));
t301 = t276 * t284 + t277 * t288;
t302 = t276 * t288 - t277 * t284;
t334 = (t301 * r_i_i_C(1) + t302 * r_i_i_C(2)) * qJD(6);
t326 = pkin(4) + pkin(5);
t308 = t284 * r_i_i_C(2) - t326;
t296 = t288 * r_i_i_C(1) - t308;
t310 = -t284 * r_i_i_C(1) - qJ(5);
t297 = t288 * r_i_i_C(2) - t310;
t299 = t284 * t285 + t288 * t289;
t300 = t284 * t289 - t285 * t288;
t314 = r_i_i_C(3) + pkin(9) - pkin(8);
t292 = -t285 * qJD(5) + t314 * qJD(3) + (t300 * r_i_i_C(1) + t299 * r_i_i_C(2)) * qJD(6) + (t296 * t285 - t297 * t289) * qJD(4);
t290 = cos(qJ(3));
t333 = pkin(3) * t286 + t314 * t290 + qJ(2);
t331 = (-qJD(4) + qJD(6)) * t290;
t294 = t297 * t285 + t296 * t289 + pkin(3);
t327 = -pkin(1) - pkin(7);
t306 = qJD(1) * t286 + qJD(4);
t317 = qJD(3) * t291;
t312 = t290 * t317;
t320 = qJD(1) * t291;
t270 = -qJD(4) * t313 - t285 * t312 - t289 * t320 + t306 * t325;
t269 = t270 * t288;
t321 = qJD(1) * t287;
t319 = qJD(3) * t286;
t318 = qJD(3) * t290;
t316 = qJD(4) * t290;
t309 = t314 * t286;
t307 = qJD(4) * t286 + qJD(1);
t305 = -(t299 * t331 - t300 * t319) * r_i_i_C(1) + (t299 * t319 + t300 * t331) * r_i_i_C(2);
t274 = t286 * t325 - t322;
t275 = t286 * t324 + t323;
t304 = t274 * t288 - t275 * t284;
t303 = t274 * t284 + t275 * t288;
t298 = t306 * t291;
t295 = qJD(2) + (t290 * pkin(3) - t309) * qJD(3);
t293 = qJD(3) * t294;
t273 = t289 * t298 + (-t307 * t285 + t289 * t318) * t287;
t272 = t307 * t324 + (t287 * t318 + t298) * t285;
t271 = t275 * qJD(1) + t276 * qJD(4) - t289 * t312;
t266 = t304 * qJD(6) + t272 * t284 + t273 * t288;
t265 = -t303 * qJD(6) + t272 * t288 - t273 * t284;
t1 = [-t269 * r_i_i_C(2) + t276 * qJD(5) + t310 * t270 - t296 * t271 + (t302 * r_i_i_C(1) - t301 * r_i_i_C(2)) * qJD(6) + t295 * t291 + (-t333 * t287 + t327 * t291) * qJD(1), t320 (t294 * t290 - t309) * t320 + (-t286 * t293 - t292 * t290) * t287, t275 * qJD(5) + t297 * t273 - t296 * t272 + (t303 * r_i_i_C(1) + t304 * r_i_i_C(2)) * qJD(6), t272, t265 * r_i_i_C(1) - t266 * r_i_i_C(2); t266 * r_i_i_C(1) + t265 * r_i_i_C(2) + t272 * qJ(5) + t274 * qJD(5) + t326 * t273 + t295 * t287 + (t327 * t287 + t333 * t291) * qJD(1), t321 (t294 * t317 - t314 * t321) * t286 + (t292 * t291 + t294 * t321) * t290, -t269 * r_i_i_C(1) - t277 * qJD(5) + t308 * t270 + t297 * t271 - t334, t270 (-t271 * t284 + t269) * r_i_i_C(1) + (-t270 * t284 - t271 * t288) * r_i_i_C(2) + t334; 0, 0, t292 * t286 - t290 * t293 (-qJ(5) * t316 + t326 * t319) * t285 + (-qJ(5) * t319 + (-t326 * qJD(4) + qJD(5)) * t290) * t289 - t305, -t285 * t319 + t289 * t316, t305;];
JaD_transl  = t1;
