% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:38:12
% EndTime: 2019-02-26 20:38:13
% DurationCPUTime: 0.26s
% Computational Cost: add. (392->67), mult. (412->99), div. (0->0), fcn. (311->9), ass. (0->57)
t306 = pkin(9) + r_i_i_C(3);
t265 = cos(qJ(6));
t303 = r_i_i_C(1) * t265;
t312 = pkin(5) + t303;
t261 = pkin(10) + qJ(4);
t256 = sin(t261);
t258 = qJ(5) + t261;
t255 = cos(t258);
t290 = t306 * t255;
t254 = sin(t258);
t305 = pkin(5) * t254;
t310 = t290 - qJ(2) - pkin(4) * t256 - sin(pkin(10)) * pkin(3) - t305;
t266 = cos(qJ(1));
t276 = qJD(1) * t254 + qJD(6);
t309 = t266 * t276;
t308 = t306 * t254;
t264 = sin(qJ(1));
t262 = qJD(4) + qJD(5);
t297 = t262 * t266;
t307 = -t255 * t297 + t264 * t276;
t304 = pkin(5) * t255;
t263 = sin(qJ(6));
t302 = r_i_i_C(2) * t263;
t301 = -pkin(1) - pkin(8) - pkin(7) - qJ(3);
t300 = pkin(4) * qJD(4);
t299 = t254 * t263;
t298 = t262 * t265;
t296 = qJD(1) * t264;
t259 = qJD(1) * t266;
t295 = qJD(6) * t263;
t294 = qJD(6) * t265;
t293 = t255 * t302;
t292 = t256 * t300;
t257 = cos(t261);
t291 = t257 * t300;
t289 = t262 * t299;
t287 = t255 * t262 * t264;
t283 = t255 * t259;
t282 = t255 * t295;
t281 = t255 * t294;
t280 = r_i_i_C(2) * t289;
t279 = qJD(1) * t255 * t303;
t278 = r_i_i_C(2) * t281;
t277 = -qJD(6) * t254 - qJD(1);
t274 = t277 * t266;
t273 = pkin(5) * t283 + t264 * t280 + t266 * t279 + t306 * (t254 * t259 + t287);
t272 = qJD(1) * (pkin(4) * t257 - t293);
t271 = t254 * t298 + t282;
t270 = t264 * t279 + t296 * t304 + (r_i_i_C(1) * t282 + t278) * t266 + (t306 * t296 + t312 * t297) * t254;
t269 = t291 + qJD(2) + (t304 + t308) * t262;
t268 = (-t312 * t255 + t293 - t308) * t262 + (r_i_i_C(1) * t295 + r_i_i_C(2) * t294) * t254;
t267 = -r_i_i_C(1) * t271 - t262 * t305 - t278;
t234 = t265 * t309 + (t255 * t298 + t263 * t277) * t264;
t233 = t277 * t265 * t264 + (-t287 - t309) * t263;
t232 = t263 * t274 - t307 * t265;
t231 = t307 * t263 + t265 * t274;
t1 = [t232 * r_i_i_C(1) + t231 * r_i_i_C(2) - t264 * qJD(3) + t269 * t266 + (t310 * t264 + t301 * t266) * qJD(1), t259, -t296, t266 * t272 + (t267 - t292) * t264 + t273, t264 * t267 - t283 * t302 + t273, t233 * r_i_i_C(1) - t234 * r_i_i_C(2); t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t266 * qJD(3) + t269 * t264 + (t301 * t264 - t310 * t266) * qJD(1), t296, t259, t264 * t272 + (t292 + (-r_i_i_C(2) * t299 - t290) * t262) * t266 + t270, -t266 * t280 + (-t296 * t302 - t297 * t306) * t255 + t270, -t231 * r_i_i_C(1) + t232 * r_i_i_C(2); 0, 0, 0, t268 - t291, t268, t271 * r_i_i_C(2) + (-t281 + t289) * r_i_i_C(1);];
JaD_transl  = t1;
