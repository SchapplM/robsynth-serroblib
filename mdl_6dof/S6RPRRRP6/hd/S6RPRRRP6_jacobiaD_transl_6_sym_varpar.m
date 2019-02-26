% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:51
% EndTime: 2019-02-26 21:10:51
% DurationCPUTime: 0.29s
% Computational Cost: add. (502->68), mult. (528->95), div. (0->0), fcn. (416->9), ass. (0->56)
t255 = qJ(4) + qJ(5);
t249 = sin(t255);
t253 = pkin(10) + qJ(3);
t247 = sin(t253);
t260 = cos(qJ(1));
t254 = qJD(4) + qJD(5);
t248 = cos(t253);
t286 = qJD(1) * t248;
t273 = -t254 + t286;
t258 = sin(qJ(1));
t281 = qJD(3) * t258;
t298 = -t247 * t281 + t273 * t260;
t304 = t249 * t298;
t257 = sin(qJ(4));
t291 = t249 * t254;
t293 = pkin(4) * qJD(4);
t239 = -pkin(5) * t291 - t257 * t293;
t250 = cos(t255);
t259 = cos(qJ(4));
t244 = t259 * pkin(4) + pkin(5) * t250;
t242 = pkin(3) + t244;
t243 = t257 * pkin(4) + pkin(5) * t249;
t279 = t247 * qJD(6);
t294 = r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
t300 = t294 * t248;
t303 = (-t242 * t247 + t300) * qJD(3) + (t243 + pkin(7) + qJ(2)) * qJD(1) + t248 * t239 + t279;
t267 = r_i_i_C(1) * t250 - r_i_i_C(2) * t249 + t242;
t263 = -t267 * t247 + t300;
t295 = r_i_i_C(2) * t250;
t272 = r_i_i_C(1) * t249 + t295;
t299 = t272 * t254 - t239;
t296 = -pkin(5) - r_i_i_C(1);
t290 = t250 * t254;
t280 = qJD(3) * t260;
t264 = t247 * t280 + t273 * t258;
t274 = t248 * t254 - qJD(1);
t270 = t250 * t274;
t235 = t264 * t249 - t260 * t270;
t236 = t274 * t260 * t249 + t264 * t250;
t289 = t235 * r_i_i_C(1) + t236 * r_i_i_C(2);
t269 = t274 * t258;
t237 = t250 * t269 + t304;
t238 = t249 * t269 - t250 * t298;
t288 = -t237 * r_i_i_C(1) + t238 * r_i_i_C(2);
t285 = qJD(1) * t258;
t284 = qJD(1) * t260;
t283 = qJD(3) * t247;
t282 = qJD(3) * t248;
t277 = t294 * t247;
t271 = t243 * t286 + t239;
t240 = pkin(5) * t290 + t259 * t293;
t266 = qJD(1) * t244 - t240 * t248 + t243 * t283;
t262 = qJD(2) + t240 + (-t242 * t248 - cos(pkin(10)) * pkin(2) - pkin(1) - t277) * qJD(1);
t261 = qJD(6) * t248 + t299 * t247 + (-t267 * t248 - t277) * qJD(3);
t241 = t247 * r_i_i_C(2) * t291;
t1 = [t238 * r_i_i_C(1) + t237 * r_i_i_C(2) - t303 * t258 + t262 * t260, t284, t261 * t260 - t263 * t285, t271 * t258 + t266 * t260 + t289, t235 * pkin(5) + t289, -t247 * t285 + t248 * t280; -t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t262 * t258 + t303 * t260, t285, t261 * t258 + t263 * t284, t266 * t258 - t271 * t260 + t288 (-t258 * t270 - t304) * pkin(5) + t288, t247 * t284 + t248 * t281; 0, 0, t263 * qJD(3) - t248 * t299 + t279, t241 + (-r_i_i_C(1) * t290 - t240) * t247 + (-t243 - t272) * t282, t241 + t296 * t247 * t290 + (t296 * t249 - t295) * t282, t283;];
JaD_transl  = t1;
