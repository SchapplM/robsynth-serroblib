% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:02
% EndTime: 2019-02-26 21:54:02
% DurationCPUTime: 0.32s
% Computational Cost: add. (745->76), mult. (542->107), div. (0->0), fcn. (399->12), ass. (0->64)
t278 = qJ(2) + pkin(11);
t275 = qJ(4) + t278;
t271 = qJ(5) + t275;
t267 = sin(t271);
t282 = cos(qJ(6));
t338 = r_i_i_C(1) * t282 + pkin(5);
t301 = t338 * t267;
t279 = sin(qJ(6));
t321 = qJD(6) * t282;
t268 = cos(t271);
t277 = qJD(2) + qJD(4);
t274 = qJD(5) + t277;
t329 = t268 * t274;
t340 = t267 * t321 + t279 * t329;
t335 = pkin(10) + r_i_i_C(3);
t315 = t335 * t268;
t298 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t278);
t293 = t298 * qJD(2);
t269 = sin(t275);
t333 = pkin(4) * t277;
t320 = t269 * t333;
t332 = pkin(5) * t267;
t339 = (t315 - t332) * t274 + t293 - t320;
t322 = qJD(6) * t279;
t308 = t267 * t322;
t336 = r_i_i_C(1) * t308 + t340 * r_i_i_C(2);
t270 = cos(t275);
t300 = t338 * t268;
t316 = t335 * t267;
t287 = (-t300 - t316) * t274 - t270 * t333;
t334 = pkin(4) * t269;
t330 = r_i_i_C(2) * t279;
t281 = sin(qJ(1));
t328 = t274 * t281;
t327 = t274 * t282;
t284 = cos(qJ(1));
t326 = t274 * t284;
t325 = t282 * t284;
t324 = qJD(1) * t281;
t323 = qJD(1) * t284;
t318 = t267 * t330;
t317 = qJD(1) * t330;
t314 = t335 * t281;
t313 = t267 * t327;
t303 = qJD(6) * t268 - qJD(1);
t302 = qJD(1) * t268 - qJD(6);
t299 = t338 * t284;
t297 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t278);
t296 = t336 * t284 + t324 * t301;
t295 = t303 * t279;
t294 = t284 * t267 * t317 + t336 * t281 + t323 * t315;
t292 = -t315 - t318;
t291 = -pkin(4) * t270 - pkin(5) * t268 - pkin(1) + t297 - t316;
t290 = t267 * t326 + t302 * t281;
t288 = t297 * qJD(2) + t287;
t286 = -t268 * r_i_i_C(2) * t321 + (-t268 * t322 - t313) * r_i_i_C(1) + t335 * t329 + (-t332 + t318) * t274;
t285 = t286 - t320;
t276 = -pkin(9) - pkin(8) - qJ(3) - pkin(7);
t252 = t298 - t334;
t248 = -t302 * t325 + (t295 + t313) * t281;
t247 = t303 * t282 * t281 + (-t267 * t328 + t302 * t284) * t279;
t246 = t290 * t282 + t284 * t295;
t245 = t290 * t279 - t303 * t325;
t1 = [t248 * r_i_i_C(1) + t247 * r_i_i_C(2) + t284 * qJD(3) - t339 * t281 + (t276 * t281 + t291 * t284) * qJD(1) (-t252 + t292) * t324 + t288 * t284 + t296, t323 (t292 + t334) * t324 + t287 * t284 + t296 (-t281 * t317 - t335 * t326) * t267 + (-qJD(1) * t314 - t274 * t299) * t268 + t296, t245 * r_i_i_C(1) + t246 * r_i_i_C(2); -t246 * r_i_i_C(1) + t245 * r_i_i_C(2) + t281 * qJD(3) + t339 * t284 + (-t276 * t284 + t291 * t281) * qJD(1) (t252 - t301) * t323 + t288 * t281 + t294, t324 (-t301 - t334) * t323 + t287 * t281 + t294, -t300 * t328 + (-qJD(1) * t299 - t274 * t314) * t267 + t294, -t247 * r_i_i_C(1) + t248 * r_i_i_C(2); 0, t293 + t285, 0, t285, t286 (-t268 * t327 + t308) * r_i_i_C(2) - t340 * r_i_i_C(1);];
JaD_transl  = t1;
