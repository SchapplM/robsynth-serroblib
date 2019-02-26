% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:48
% EndTime: 2019-02-26 22:48:48
% DurationCPUTime: 0.43s
% Computational Cost: add. (1068->90), mult. (741->120), div. (0->0), fcn. (584->12), ass. (0->74)
t264 = qJ(3) + qJ(4);
t260 = qJ(5) + t264;
t254 = sin(t260);
t258 = sin(t264);
t263 = qJD(3) + qJD(4);
t257 = qJD(5) + t263;
t307 = pkin(5) * t257;
t309 = pkin(4) * t263;
t241 = -t254 * t307 - t258 * t309;
t265 = sin(qJ(3));
t301 = pkin(3) * qJD(3);
t239 = -t265 * t301 + t241;
t255 = cos(t260);
t259 = cos(t264);
t248 = pkin(4) * t259 + pkin(5) * t255;
t268 = cos(qJ(3));
t245 = t268 * pkin(3) + t248;
t243 = pkin(2) + t245;
t266 = sin(qJ(2));
t269 = cos(qJ(2));
t302 = r_i_i_C(3) + pkin(11) + pkin(10) + pkin(9) + pkin(8);
t286 = t302 * t269;
t314 = (-t243 * t266 + t286) * qJD(2) + t269 * t239;
t270 = cos(qJ(1));
t252 = qJD(6) + t257;
t285 = t252 * t269 - qJD(1);
t312 = t270 * t285;
t256 = qJ(6) + t260;
t250 = sin(t256);
t251 = cos(t256);
t304 = r_i_i_C(2) * t251;
t282 = r_i_i_C(1) * t250 + t304;
t311 = -t282 * t252 + t239;
t294 = qJD(1) * t269;
t284 = -t252 + t294;
t267 = sin(qJ(1));
t292 = qJD(2) * t266;
t288 = t267 * t292;
t310 = t284 * t270 - t288;
t308 = pkin(5) * t254;
t306 = r_i_i_C(1) * t251;
t305 = r_i_i_C(2) * t250;
t247 = -pkin(4) * t258 - t308;
t244 = t265 * pkin(3) - t247;
t303 = pkin(7) + t244;
t299 = t252 * t266;
t287 = t270 * t292;
t272 = t284 * t267 + t287;
t235 = t272 * t250 - t251 * t312;
t236 = t250 * t312 + t272 * t251;
t297 = t235 * r_i_i_C(1) + t236 * r_i_i_C(2);
t278 = t285 * t267;
t237 = t310 * t250 + t251 * t278;
t238 = t250 * t278 - t310 * t251;
t296 = -t237 * r_i_i_C(1) + t238 * r_i_i_C(2);
t295 = qJD(1) * t267;
t293 = qJD(1) * t270;
t291 = qJD(2) * t269;
t290 = t255 * t307;
t289 = t252 * t306;
t283 = -t257 + t294;
t281 = t244 * t294 + t239;
t280 = t247 * t294 - t241;
t279 = t255 * (-t257 * t269 + qJD(1));
t277 = t243 - t305 + t306;
t242 = -t259 * t309 - t290;
t276 = -t243 * t269 - t302 * t266 - pkin(1);
t275 = qJD(2) * t277;
t240 = t268 * t301 - t242;
t274 = qJD(1) * t245 - t240 * t269 + t244 * t292;
t273 = qJD(1) * t248 + t242 * t269 - t247 * t292;
t271 = -t302 * qJD(2) - t311;
t246 = t299 * t305;
t1 = [t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t270 * t240 - t314 * t267 + (-t303 * t267 + t276 * t270) * qJD(1) (-t270 * t275 - t302 * t295) * t269 + (t271 * t270 + t277 * t295) * t266, t281 * t267 + t274 * t270 + t297, -t280 * t267 + t273 * t270 + t297 (t270 * t279 + (t283 * t267 + t287) * t254) * pkin(5) + t297, t297; -t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t267 * t240 + t314 * t270 + (t276 * t267 + t303 * t270) * qJD(1) (-t267 * t275 + t302 * t293) * t269 + (t271 * t267 - t277 * t293) * t266, t274 * t267 - t281 * t270 + t296, t273 * t267 + t280 * t270 + t296 (t267 * t279 + (-t283 * t270 + t288) * t254) * pkin(5) + t296, t296; 0, t311 * t269 + (-t277 * t266 + t286) * qJD(2), t246 + (-t240 - t289) * t266 + (-t244 - t282) * t291, t246 + (t242 - t289) * t266 + (t247 - t282) * t291, t246 + (-t289 - t290) * t266 + (-t282 - t308) * t291, -t291 * t304 + t246 + (-t250 * t291 - t251 * t299) * r_i_i_C(1);];
JaD_transl  = t1;
