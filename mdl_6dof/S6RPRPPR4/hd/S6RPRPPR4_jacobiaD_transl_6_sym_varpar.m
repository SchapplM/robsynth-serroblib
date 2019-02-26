% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:55
% EndTime: 2019-02-26 20:40:55
% DurationCPUTime: 0.44s
% Computational Cost: add. (380->71), mult. (665->113), div. (0->0), fcn. (606->9), ass. (0->51)
t258 = pkin(9) + qJ(3);
t256 = sin(t258);
t259 = sin(pkin(10));
t260 = cos(pkin(10));
t262 = sin(qJ(6));
t264 = cos(qJ(6));
t296 = pkin(4) + pkin(5);
t271 = t264 * r_i_i_C(1) - t262 * r_i_i_C(2) + t296;
t272 = t262 * r_i_i_C(1) + t264 * r_i_i_C(2) + qJ(5);
t297 = t272 * t259 + t271 * t260 + pkin(3);
t257 = cos(t258);
t283 = -r_i_i_C(3) - pkin(8) + qJ(4);
t299 = t283 * t257;
t302 = t297 * t256 - t299;
t285 = t256 * qJD(4);
t301 = (-t256 * pkin(3) + t299) * qJD(3) + t285;
t273 = t259 * t262 + t260 * t264;
t274 = t259 * t264 - t260 * t262;
t269 = t274 * r_i_i_C(1) - t273 * r_i_i_C(2);
t298 = t259 * qJD(5) + t269 * qJD(6);
t265 = cos(qJ(1));
t294 = t260 * t265;
t263 = sin(qJ(1));
t293 = t263 * t259;
t292 = t263 * t260;
t291 = t265 * t259;
t290 = qJD(1) * t263;
t289 = qJD(1) * t265;
t288 = qJD(3) * t257;
t287 = qJD(3) * t263;
t286 = qJD(3) * t265;
t282 = t257 * t291;
t281 = t256 * t287;
t280 = t256 * t286;
t279 = t283 * t256;
t247 = t257 * t293 + t294;
t248 = t257 * t292 - t291;
t276 = -t247 * t264 + t248 * t262;
t275 = t247 * t262 + t248 * t264;
t250 = t257 * t294 + t293;
t270 = -pkin(3) * t257 - cos(pkin(9)) * pkin(2) - pkin(1) - t279;
t266 = t257 * qJD(4) - t298 * t256 + (-t257 * t297 - t279) * qJD(3);
t261 = -pkin(7) - qJ(2);
t249 = t282 - t292;
t246 = t250 * qJD(1) - t260 * t281;
t245 = qJD(1) * t282 - t259 * t281 - t260 * t290;
t244 = -t248 * qJD(1) - t260 * t280;
t243 = t247 * qJD(1) + t259 * t280;
t242 = -t243 * t262 + t244 * t264 + (t249 * t264 - t250 * t262) * qJD(6);
t241 = -t243 * t264 - t244 * t262 + (-t249 * t262 - t250 * t264) * qJD(6);
t1 = [t265 * qJD(2) - t247 * qJD(5) - t272 * t245 - t271 * t246 + (t276 * r_i_i_C(1) + t275 * r_i_i_C(2)) * qJD(6) - t301 * t263 + (t263 * t261 + t270 * t265) * qJD(1), t289, t266 * t265 + t302 * t290, -t256 * t290 + t257 * t286, -t243, r_i_i_C(1) * t241 - r_i_i_C(2) * t242; t242 * r_i_i_C(1) + t241 * r_i_i_C(2) - t243 * qJ(5) + t263 * qJD(2) + t249 * qJD(5) + t296 * t244 + t301 * t265 + (-t261 * t265 + t270 * t263) * qJD(1), t290, t266 * t263 - t289 * t302, t256 * t289 + t257 * t287, t245 (t245 * t264 - t246 * t262) * r_i_i_C(1) + (-t245 * t262 - t246 * t264) * r_i_i_C(2) + (-t275 * r_i_i_C(1) + t276 * r_i_i_C(2)) * qJD(6); 0, 0, -qJD(3) * t302 + t298 * t257 + t285, qJD(3) * t256, t259 * t288 (-t273 * r_i_i_C(1) - t274 * r_i_i_C(2)) * t256 * qJD(6) + t269 * t288;];
JaD_transl  = t1;
