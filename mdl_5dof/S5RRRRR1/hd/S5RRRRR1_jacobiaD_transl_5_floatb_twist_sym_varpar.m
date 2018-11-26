% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR1_jacobiaD_transl_5_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_jacobiaD_transl_5_floatb_twist_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_jacobiaD_transl_5_floatb_twist_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR1_jacobiaD_transl_5_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_jacobiaD_transl_5_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:52:31
% EndTime: 2018-11-16 14:52:31
% DurationCPUTime: 0.36s
% Computational Cost: add. (570->67), mult. (520->100), div. (0->0), fcn. (384->10), ass. (0->61)
t269 = qJ(2) + qJ(3);
t267 = qJ(4) + t269;
t262 = sin(t267);
t273 = cos(qJ(5));
t325 = r_i_i_C(1) * t273 + pkin(4);
t289 = t325 * t262;
t263 = cos(t267);
t268 = qJD(2) + qJD(3);
t264 = qJD(4) + t268;
t270 = sin(qJ(5));
t308 = qJD(5) * t273;
t327 = t263 * t264 * t270 + t262 * t308;
t322 = pkin(6) + r_i_i_C(3);
t303 = t322 * t263;
t265 = sin(t269);
t320 = pkin(3) * t268;
t261 = t265 * t320;
t271 = sin(qJ(2));
t316 = pkin(2) * qJD(2);
t305 = t271 * t316;
t319 = pkin(4) * t262;
t326 = (t303 - t319) * t264 - t261 - t305;
t309 = qJD(5) * t270;
t296 = t262 * t309;
t323 = r_i_i_C(1) * t296 + t327 * r_i_i_C(2);
t266 = cos(t269);
t288 = t325 * t263;
t304 = t322 * t262;
t276 = (-t288 - t304) * t264 - t266 * t320;
t321 = pkin(3) * t265;
t317 = r_i_i_C(2) * t270;
t272 = sin(qJ(1));
t315 = t264 * t272;
t314 = t264 * t273;
t275 = cos(qJ(1));
t313 = t264 * t275;
t312 = t273 * t275;
t311 = qJD(1) * t272;
t310 = qJD(1) * t275;
t306 = qJD(1) * t317;
t302 = t322 * t272;
t301 = t262 * t314;
t291 = qJD(5) * t263 + qJD(1);
t290 = qJD(1) * t263 + qJD(5);
t287 = t325 * t275;
t286 = t323 * t275 + t311 * t289;
t285 = t291 * t270;
t284 = t275 * t262 * t306 + t323 * t272 + t310 * t303;
t283 = -t262 * t317 - t303;
t274 = cos(qJ(2));
t282 = qJD(1) * (-t274 * pkin(2) - pkin(3) * t266 - pkin(4) * t263 - pkin(1) - t304);
t281 = t262 * t313 + t290 * t272;
t279 = -t274 * t316 + t276;
t278 = t263 * r_i_i_C(2) * t308 + (t283 + t319) * t264 + (t263 * t309 + t301) * r_i_i_C(1);
t277 = t261 + t278;
t260 = -t271 * pkin(2) - t321;
t241 = -t290 * t312 + (t285 + t301) * t272;
t240 = t291 * t273 * t272 + (-t262 * t315 + t290 * t275) * t270;
t239 = t281 * t273 + t275 * t285;
t238 = t281 * t270 - t291 * t312;
t1 = [t241 * r_i_i_C(1) + t240 * r_i_i_C(2) - t326 * t272 + t275 * t282 (-t260 + t283) * t311 + t279 * t275 + t286 (t283 + t321) * t311 + t276 * t275 + t286 (-t272 * t306 - t322 * t313) * t262 + (-qJD(1) * t302 - t264 * t287) * t263 + t286, t238 * r_i_i_C(1) + t239 * r_i_i_C(2); -t239 * r_i_i_C(1) + t238 * r_i_i_C(2) + t272 * t282 + t326 * t275 (t260 - t289) * t310 + t279 * t272 + t284 (-t289 - t321) * t310 + t276 * t272 + t284, -t288 * t315 + (-qJD(1) * t287 - t264 * t302) * t262 + t284, -t240 * r_i_i_C(1) + t241 * r_i_i_C(2); 0, t277 + t305, t277, t278 (t263 * t314 - t296) * r_i_i_C(2) + t327 * r_i_i_C(1);];
JaD_transl  = t1;
