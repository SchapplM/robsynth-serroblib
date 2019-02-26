% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:00
% EndTime: 2019-02-26 20:54:00
% DurationCPUTime: 0.27s
% Computational Cost: add. (406->58), mult. (457->86), div. (0->0), fcn. (367->10), ass. (0->52)
t238 = pkin(10) + qJ(5);
t235 = cos(t238);
t228 = pkin(5) * t235 + cos(pkin(10)) * pkin(4) + pkin(3);
t240 = sin(qJ(3));
t265 = qJD(5) * t235;
t242 = cos(qJ(3));
t276 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
t284 = t276 * t242;
t296 = -pkin(5) * t265 + (-t228 * t240 - qJ(2) + t284) * qJD(1);
t239 = qJD(5) + qJD(6);
t234 = sin(t238);
t281 = pkin(5) * t234;
t263 = qJD(5) * t281;
t236 = qJ(6) + t238;
t232 = sin(t236);
t233 = cos(t236);
t293 = r_i_i_C(1) * t232 + r_i_i_C(2) * t233;
t246 = t239 * t293 + t263;
t292 = -r_i_i_C(1) * t233 + r_i_i_C(2) * t232;
t248 = t228 - t292;
t295 = (-t240 * t248 + t284) * qJD(3) + qJD(4) * t240 - t246 * t242;
t258 = t276 * t240;
t294 = t242 * t248 + t258;
t288 = t235 * (qJD(5) * t240 + qJD(1));
t241 = sin(qJ(1));
t272 = qJD(1) * t240;
t254 = t239 + t272;
t243 = cos(qJ(1));
t267 = qJD(3) * t243;
t260 = t242 * t267;
t283 = t241 * t254 - t260;
t268 = qJD(3) * t242;
t261 = t241 * t268;
t282 = t243 * t254 + t261;
t255 = -t239 * t240 - qJD(1);
t249 = t255 * t243;
t223 = t232 * t283 + t233 * t249;
t224 = t232 * t249 - t233 * t283;
t274 = -r_i_i_C(1) * t223 + r_i_i_C(2) * t224;
t250 = t255 * t241;
t225 = -t232 * t282 + t233 * t250;
t226 = t232 * t250 + t233 * t282;
t273 = r_i_i_C(1) * t225 - r_i_i_C(2) * t226;
t271 = qJD(1) * t241;
t270 = qJD(1) * t243;
t269 = qJD(3) * t240;
t264 = t242 * qJD(4);
t252 = -qJD(5) - t272;
t247 = t239 * t242 * t292 + t269 * t293;
t245 = qJD(1) * t294;
t244 = -t240 * t263 - t264 + qJD(2) + (t228 * t242 + t258) * qJD(3) + (-pkin(1) - pkin(7) - t281 - sin(pkin(10)) * pkin(4)) * qJD(1);
t1 = [t224 * r_i_i_C(1) + t223 * r_i_i_C(2) + t241 * t296 + t243 * t244, t270, t241 * t295 + t243 * t245, t241 * t269 - t242 * t270 (-t241 * t288 + (t243 * t252 - t261) * t234) * pkin(5) + t273, t273; t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t241 * t244 - t243 * t296, t271, t241 * t245 - t243 * t295, -t240 * t267 - t242 * t271 (t243 * t288 + (t241 * t252 + t260) * t234) * pkin(5) + t274, t274; 0, 0, -qJD(3) * t294 + t246 * t240 + t264, t268 (t234 * t269 - t242 * t265) * pkin(5) + t247, t247;];
JaD_transl  = t1;
