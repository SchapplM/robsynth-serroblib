% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:14
% EndTime: 2019-02-26 20:49:14
% DurationCPUTime: 0.30s
% Computational Cost: add. (510->64), mult. (416->95), div. (0->0), fcn. (312->12), ass. (0->52)
t277 = qJ(3) + pkin(11);
t274 = qJ(5) + t277;
t268 = sin(t274);
t281 = cos(qJ(6));
t325 = r_i_i_C(1) * t281 + pkin(5);
t330 = t268 * t325;
t269 = cos(t274);
t309 = qJD(6) * t281;
t276 = qJD(3) + qJD(5);
t279 = sin(qJ(6));
t314 = t276 * t279;
t329 = t268 * t309 + t269 * t314;
t319 = pkin(9) + r_i_i_C(3);
t327 = t319 * t269;
t328 = (-pkin(5) * t268 + t327) * t276;
t265 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t277);
t258 = t265 * qJD(3);
t326 = t258 + t328;
t310 = qJD(6) * t279;
t298 = t268 * t310;
t322 = r_i_i_C(1) * t298 + t329 * r_i_i_C(2);
t292 = qJD(1) * t269 - qJD(6);
t321 = t281 * t292;
t293 = qJD(6) * t269 - qJD(1);
t304 = t268 * t314;
t320 = t293 * t281 - t304;
t316 = r_i_i_C(2) * t279;
t313 = t276 * t281;
t278 = qJ(1) + pkin(10);
t271 = sin(t278);
t312 = qJD(1) * t271;
t273 = cos(t278);
t311 = qJD(1) * t273;
t308 = qJD(1) * t316;
t307 = t319 * t268;
t305 = t319 * t276;
t303 = t268 * t313;
t291 = t325 * t276;
t290 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t277);
t289 = t322 * t273 + t312 * t330;
t288 = t292 * t279;
t287 = t273 * t268 * t308 + t322 * t271 + t311 * t327;
t286 = -pkin(5) * t269 - pkin(2) + t290 - t307;
t285 = t293 * t279 + t303;
t284 = t290 * qJD(3) + (-t269 * t325 - t307) * t276;
t283 = (-t269 * t310 - t303) * r_i_i_C(1) + (-t269 * t309 + t304) * r_i_i_C(2) + t328;
t275 = -pkin(8) - qJ(4) - pkin(7);
t249 = t285 * t271 - t273 * t321;
t248 = t320 * t271 + t273 * t288;
t247 = t271 * t321 + t285 * t273;
t246 = t271 * t288 - t320 * t273;
t1 = [t249 * r_i_i_C(1) + t248 * r_i_i_C(2) + t273 * qJD(4) - t326 * t271 + (-cos(qJ(1)) * pkin(1) + t271 * t275 + t286 * t273) * qJD(1), 0 (-t268 * t316 - t265 - t327) * t312 + t284 * t273 + t289, t311 (-t271 * t308 - t273 * t305) * t268 + (-t273 * t291 - t319 * t312) * t269 + t289, t246 * r_i_i_C(1) + t247 * r_i_i_C(2); -t247 * r_i_i_C(1) + t246 * r_i_i_C(2) + t271 * qJD(4) + t326 * t273 + (-sin(qJ(1)) * pkin(1) - t273 * t275 + t286 * t271) * qJD(1), 0 (t265 - t330) * t311 + t284 * t271 + t287, t312, -t271 * t269 * t291 + (-t271 * t305 - t311 * t325) * t268 + t287, -t248 * r_i_i_C(1) + t249 * r_i_i_C(2); 0, 0, t258 + t283, 0, t283 (-t269 * t313 + t298) * r_i_i_C(2) - t329 * r_i_i_C(1);];
JaD_transl  = t1;
