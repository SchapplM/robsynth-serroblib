% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR3
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
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:01
% EndTime: 2019-02-26 21:16:01
% DurationCPUTime: 0.29s
% Computational Cost: add. (728->72), mult. (576->99), div. (0->0), fcn. (454->12), ass. (0->62)
t255 = sin(qJ(3));
t253 = qJ(4) + qJ(5);
t247 = cos(t253);
t256 = cos(qJ(4));
t239 = t256 * pkin(4) + pkin(5) * t247;
t237 = pkin(3) + t239;
t248 = qJ(6) + t253;
t241 = sin(t248);
t289 = r_i_i_C(2) * t241;
t242 = cos(t248);
t290 = r_i_i_C(1) * t242;
t264 = t237 - t289 + t290;
t257 = cos(qJ(3));
t286 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t295 = t286 * t257;
t259 = -t264 * t255 + t295;
t301 = qJD(1) * t259;
t246 = sin(t253);
t254 = sin(qJ(4));
t285 = pkin(4) * qJD(4);
t250 = qJD(4) + qJD(5);
t291 = pkin(5) * t250;
t234 = -t246 * t291 - t254 * t285;
t300 = (-t237 * t255 + t295) * qJD(3) + t257 * t234;
t251 = qJ(1) + pkin(11);
t244 = cos(t251);
t245 = qJD(6) + t250;
t270 = t245 * t257 - qJD(1);
t298 = t244 * t270;
t279 = qJD(1) * t257;
t297 = t246 * (-t250 + t279);
t288 = r_i_i_C(2) * t242;
t267 = r_i_i_C(1) * t241 + t288;
t294 = t267 * t245 - t234;
t243 = sin(t251);
t269 = -t245 + t279;
t278 = qJD(3) * t255;
t293 = -t243 * t278 + t269 * t244;
t292 = pkin(5) * t246;
t238 = pkin(4) * t254 + t292;
t287 = pkin(7) + t238;
t283 = t245 * t255;
t261 = t269 * t243 + t244 * t278;
t230 = t261 * t241 - t242 * t298;
t231 = t241 * t298 + t261 * t242;
t281 = t230 * r_i_i_C(1) + t231 * r_i_i_C(2);
t265 = t270 * t243;
t232 = t293 * t241 + t242 * t265;
t233 = t241 * t265 - t293 * t242;
t280 = -t232 * r_i_i_C(1) + t233 * r_i_i_C(2);
t277 = qJD(3) * t257;
t276 = t247 * t291;
t275 = t245 * t290;
t273 = t286 * t255;
t266 = t238 * t279 + t234;
t263 = -t237 * t257 - pkin(2) - t273;
t235 = t256 * t285 + t276;
t262 = qJD(1) * t239 - t235 * t257 + t238 * t278;
t260 = t246 * t278 + (-t250 * t257 + qJD(1)) * t247;
t258 = t294 * t255 + (-t264 * t257 - t273) * qJD(3);
t236 = t283 * t289;
t1 = [t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t244 * t235 - t300 * t243 + (-cos(qJ(1)) * pkin(1) - t287 * t243 + t263 * t244) * qJD(1), 0, -t243 * t301 + t258 * t244, t266 * t243 + t262 * t244 + t281 (t243 * t297 + t260 * t244) * pkin(5) + t281, t281; -t231 * r_i_i_C(1) + t230 * r_i_i_C(2) + t243 * t235 + t300 * t244 + (-sin(qJ(1)) * pkin(1) + t287 * t244 + t263 * t243) * qJD(1), 0, t258 * t243 + t244 * t301, t262 * t243 - t266 * t244 + t280 (t260 * t243 - t244 * t297) * pkin(5) + t280, t280; 0, 0, t259 * qJD(3) - t294 * t257, t236 + (-t235 - t275) * t255 + (-t238 - t267) * t277, t236 + (-t275 - t276) * t255 + (-t267 - t292) * t277, -t277 * t288 + t236 + (-t241 * t277 - t242 * t283) * r_i_i_C(1);];
JaD_transl  = t1;
