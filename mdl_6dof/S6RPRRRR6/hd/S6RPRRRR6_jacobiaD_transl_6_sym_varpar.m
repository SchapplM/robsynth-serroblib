% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:17:37
% EndTime: 2019-02-26 21:17:38
% DurationCPUTime: 0.31s
% Computational Cost: add. (708->76), mult. (578->106), div. (0->0), fcn. (458->11), ass. (0->65)
t263 = qJ(4) + qJ(5);
t256 = sin(t263);
t265 = sin(qJ(4));
t299 = pkin(4) * qJD(4);
t261 = qJD(4) + qJD(5);
t304 = pkin(5) * t261;
t244 = -t256 * t304 - t265 * t299;
t257 = cos(t263);
t267 = cos(qJ(4));
t248 = t267 * pkin(4) + pkin(5) * t257;
t246 = pkin(3) + t248;
t260 = pkin(11) + qJ(3);
t253 = sin(t260);
t254 = cos(t260);
t300 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t282 = t300 * t254;
t310 = t254 * t244 + (-t246 * t253 + t282) * qJD(3);
t268 = cos(qJ(1));
t255 = qJD(6) + t261;
t281 = t254 * t255 - qJD(1);
t308 = t268 * t281;
t258 = qJ(6) + t263;
t251 = sin(t258);
t252 = cos(t258);
t301 = r_i_i_C(2) * t252;
t278 = r_i_i_C(1) * t251 + t301;
t307 = -t278 * t255 + t244;
t292 = qJD(1) * t254;
t280 = -t255 + t292;
t266 = sin(qJ(1));
t289 = qJD(3) * t253;
t284 = t266 * t289;
t306 = t280 * t268 - t284;
t305 = pkin(5) * t256;
t303 = r_i_i_C(1) * t252;
t302 = r_i_i_C(2) * t251;
t297 = t253 * t255;
t283 = t268 * t289;
t270 = t280 * t266 + t283;
t239 = t270 * t251 - t252 * t308;
t240 = t251 * t308 + t270 * t252;
t295 = t239 * r_i_i_C(1) + t240 * r_i_i_C(2);
t275 = t281 * t266;
t241 = t306 * t251 + t252 * t275;
t242 = t251 * t275 - t306 * t252;
t294 = -t241 * r_i_i_C(1) + t242 * r_i_i_C(2);
t247 = pkin(4) * t265 + t305;
t293 = t247 + pkin(7) + qJ(2);
t291 = qJD(1) * t266;
t290 = qJD(1) * t268;
t288 = qJD(3) * t254;
t286 = t257 * t304;
t245 = t267 * t299 + t286;
t287 = qJD(2) + t245;
t285 = t255 * t303;
t279 = -t261 + t292;
t277 = t247 * t292 + t244;
t276 = t257 * (-t254 * t261 + qJD(1));
t274 = t246 - t302 + t303;
t273 = qJD(3) * t274;
t272 = -t246 * t254 - t300 * t253 - cos(pkin(11)) * pkin(2) - pkin(1);
t271 = qJD(1) * t248 - t245 * t254 + t247 * t289;
t269 = -t300 * qJD(3) - t307;
t243 = t297 * t302;
t1 = [t242 * r_i_i_C(1) + t241 * r_i_i_C(2) + t287 * t268 - t310 * t266 + (-t293 * t266 + t272 * t268) * qJD(1), t290 (-t268 * t273 - t300 * t291) * t254 + (t269 * t268 + t274 * t291) * t253, t277 * t266 + t271 * t268 + t295 (t268 * t276 + (t279 * t266 + t283) * t256) * pkin(5) + t295, t295; -t240 * r_i_i_C(1) + t239 * r_i_i_C(2) + t287 * t266 + t310 * t268 + (t272 * t266 + t293 * t268) * qJD(1), t291 (-t266 * t273 + t300 * t290) * t254 + (t269 * t266 - t274 * t290) * t253, t271 * t266 - t277 * t268 + t294 (t266 * t276 + (-t279 * t268 + t284) * t256) * pkin(5) + t294, t294; 0, 0, t307 * t254 + (-t274 * t253 + t282) * qJD(3), t243 + (-t245 - t285) * t253 + (-t247 - t278) * t288, t243 + (-t285 - t286) * t253 + (-t278 - t305) * t288, -t288 * t301 + t243 + (-t251 * t288 - t252 * t297) * r_i_i_C(1);];
JaD_transl  = t1;
