% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:49
% EndTime: 2019-02-26 22:41:49
% DurationCPUTime: 0.35s
% Computational Cost: add. (800->81), mult. (691->108), div. (0->0), fcn. (542->10), ass. (0->68)
t256 = qJ(3) + qJ(4);
t253 = qJ(5) + t256;
t247 = sin(t253);
t259 = sin(qJ(1));
t262 = cos(qJ(1));
t255 = qJD(3) + qJD(4);
t249 = qJD(5) + t255;
t261 = cos(qJ(2));
t289 = qJD(1) * t261;
t277 = -t249 + t289;
t258 = sin(qJ(2));
t287 = qJD(2) * t258;
t305 = -t259 * t287 + t277 * t262;
t311 = t305 * t247;
t250 = sin(t256);
t301 = pkin(5) * t247;
t302 = pkin(4) * t255;
t237 = -t249 * t301 - t250 * t302;
t257 = sin(qJ(3));
t296 = pkin(3) * qJD(3);
t235 = -t257 * t296 + t237;
t248 = cos(t253);
t251 = cos(t256);
t244 = pkin(4) * t251 + pkin(5) * t248;
t260 = cos(qJ(3));
t241 = t260 * pkin(3) + t244;
t239 = pkin(2) + t241;
t243 = -pkin(4) * t250 - t301;
t240 = t257 * pkin(3) - t243;
t284 = t258 * qJD(6);
t297 = r_i_i_C(3) + qJ(6) + pkin(10) + pkin(9) + pkin(8);
t307 = t297 * t261;
t310 = (-t239 * t258 + t307) * qJD(2) + (pkin(7) + t240) * qJD(1) + t261 * t235 + t284;
t300 = r_i_i_C(2) * t247;
t270 = r_i_i_C(1) * t248 + t239 - t300;
t264 = -t270 * t258 + t307;
t299 = r_i_i_C(2) * t248;
t276 = r_i_i_C(1) * t247 + t299;
t306 = t276 * t249 - t235;
t303 = -pkin(5) - r_i_i_C(1);
t294 = t248 * t249;
t293 = t249 * t258;
t285 = qJD(2) * t262;
t266 = t258 * t285 + t277 * t259;
t278 = t249 * t261 - qJD(1);
t273 = t248 * t278;
t231 = t266 * t247 - t262 * t273;
t232 = t278 * t262 * t247 + t266 * t248;
t292 = t231 * r_i_i_C(1) + t232 * r_i_i_C(2);
t272 = t278 * t259;
t233 = t248 * t272 + t311;
t234 = t247 * t272 - t248 * t305;
t291 = -t233 * r_i_i_C(1) + t234 * r_i_i_C(2);
t290 = qJD(1) * t259;
t288 = qJD(1) * t262;
t286 = qJD(2) * t261;
t283 = r_i_i_C(1) * t294;
t281 = t297 * t258;
t275 = t240 * t289 + t235;
t274 = t243 * t289 - t237;
t238 = -pkin(5) * t294 - t251 * t302;
t236 = t260 * t296 - t238;
t269 = qJD(1) * t241 - t236 * t261 + t240 * t287;
t268 = qJD(1) * t244 + t238 * t261 - t243 * t287;
t265 = t236 + (-t239 * t261 - pkin(1) - t281) * qJD(1);
t263 = qJD(6) * t261 + t306 * t258 + (-t270 * t261 - t281) * qJD(2);
t242 = t293 * t300;
t1 = [t234 * r_i_i_C(1) + t233 * r_i_i_C(2) - t310 * t259 + t265 * t262, t263 * t262 - t264 * t290, t275 * t259 + t269 * t262 + t292, -t274 * t259 + t268 * t262 + t292, pkin(5) * t231 + t292, -t258 * t290 + t261 * t285; -t232 * r_i_i_C(1) + t231 * r_i_i_C(2) + t265 * t259 + t310 * t262, t259 * t263 + t264 * t288, t269 * t259 - t275 * t262 + t291, t268 * t259 + t274 * t262 + t291 (-t259 * t273 - t311) * pkin(5) + t291, t258 * t288 + t259 * t286; 0, t264 * qJD(2) - t306 * t261 + t284, t242 + (-t236 - t283) * t258 + (-t240 - t276) * t286, t242 + (t238 - t283) * t258 + (t243 - t276) * t286, t242 + t303 * t248 * t293 + (t303 * t247 - t299) * t286, t287;];
JaD_transl  = t1;
