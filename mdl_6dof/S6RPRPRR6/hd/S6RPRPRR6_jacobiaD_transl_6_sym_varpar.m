% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR6
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
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:43
% EndTime: 2019-02-26 20:51:43
% DurationCPUTime: 0.30s
% Computational Cost: add. (518->63), mult. (455->90), div. (0->0), fcn. (367->11), ass. (0->56)
t252 = pkin(11) + qJ(5);
t248 = cos(t252);
t241 = pkin(5) * t248 + cos(pkin(11)) * pkin(4) + pkin(3);
t253 = pkin(10) + qJ(3);
t247 = sin(t253);
t249 = cos(t253);
t246 = sin(t252);
t289 = pkin(5) * qJD(5);
t276 = t246 * t289;
t277 = t247 * qJD(4);
t294 = pkin(5) * t246;
t290 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
t298 = t290 * t249;
t302 = (-t241 * t247 + t298) * qJD(3) + (t294 + sin(pkin(11)) * pkin(4) + pkin(7) + qJ(2)) * qJD(1) - t249 * t276 + t277;
t250 = qJ(6) + t252;
t243 = sin(t250);
t292 = r_i_i_C(2) * t243;
t244 = cos(t250);
t293 = r_i_i_C(1) * t244;
t263 = t241 - t292 + t293;
t260 = -t247 * t263 + t298;
t257 = cos(qJ(1));
t254 = qJD(5) + qJD(6);
t269 = t249 * t254 - qJD(1);
t300 = t257 * t269;
t291 = r_i_i_C(2) * t244;
t266 = r_i_i_C(1) * t243 + t291;
t297 = t254 * t266 + t276;
t283 = qJD(1) * t249;
t268 = -t254 + t283;
t256 = sin(qJ(1));
t279 = qJD(3) * t256;
t274 = t247 * t279;
t296 = t257 * t268 - t274;
t287 = t247 * t254;
t278 = qJD(3) * t257;
t273 = t247 * t278;
t261 = t256 * t268 + t273;
t236 = t243 * t261 - t244 * t300;
t237 = t243 * t300 + t244 * t261;
t286 = r_i_i_C(1) * t236 + r_i_i_C(2) * t237;
t265 = t269 * t256;
t238 = t243 * t296 + t244 * t265;
t239 = t243 * t265 - t244 * t296;
t285 = -r_i_i_C(1) * t238 + r_i_i_C(2) * t239;
t282 = qJD(1) * t256;
t281 = qJD(1) * t257;
t280 = qJD(3) * t249;
t275 = t248 * t289;
t272 = t290 * t247;
t267 = -qJD(5) + t283;
t264 = t248 * (-qJD(5) * t249 + qJD(1));
t259 = t275 + qJD(2) + (-t241 * t249 - cos(pkin(10)) * pkin(2) - pkin(1) - t272) * qJD(1);
t258 = qJD(4) * t249 + t297 * t247 + (-t249 * t263 - t272) * qJD(3);
t240 = t287 * t292;
t1 = [t239 * r_i_i_C(1) + t238 * r_i_i_C(2) - t256 * t302 + t257 * t259, t281, t257 * t258 - t260 * t282, -t247 * t282 + t249 * t278 (t257 * t264 + (t256 * t267 + t273) * t246) * pkin(5) + t286, t286; -t237 * r_i_i_C(1) + t236 * r_i_i_C(2) + t256 * t259 + t257 * t302, t282, t256 * t258 + t260 * t281, t247 * t281 + t249 * t279 (t256 * t264 + (-t257 * t267 + t274) * t246) * pkin(5) + t285, t285; 0, 0, qJD(3) * t260 - t249 * t297 + t277, qJD(3) * t247, t240 + (-t254 * t293 - t275) * t247 + (-t266 - t294) * t280, -t280 * t291 + t240 + (-t243 * t280 - t244 * t287) * r_i_i_C(1);];
JaD_transl  = t1;
