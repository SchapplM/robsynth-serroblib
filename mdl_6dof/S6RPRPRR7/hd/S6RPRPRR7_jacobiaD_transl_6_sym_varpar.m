% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR7
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
% Datum: 2019-02-26 20:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:52:13
% EndTime: 2019-02-26 20:52:13
% DurationCPUTime: 0.27s
% Computational Cost: add. (399->66), mult. (424->99), div. (0->0), fcn. (318->10), ass. (0->56)
t307 = pkin(9) + r_i_i_C(3);
t268 = cos(qJ(6));
t304 = r_i_i_C(1) * t268;
t314 = pkin(5) + t304;
t264 = qJ(3) + pkin(10);
t254 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t264);
t260 = qJ(5) + t264;
t257 = cos(t260);
t294 = t307 * t257;
t256 = sin(t260);
t306 = pkin(5) * t256;
t312 = t294 - qJ(2) - t254 - t306;
t270 = cos(qJ(1));
t280 = qJD(1) * t256 + qJD(6);
t311 = t280 * t270;
t310 = t307 * t256;
t255 = cos(qJ(3)) * pkin(3) + pkin(4) * cos(t264);
t309 = t255 * qJD(3);
t267 = sin(qJ(1));
t263 = qJD(3) + qJD(5);
t299 = t263 * t270;
t308 = -t257 * t299 + t280 * t267;
t305 = pkin(5) * t257;
t265 = sin(qJ(6));
t303 = r_i_i_C(2) * t265;
t302 = -pkin(1) - pkin(8) - qJ(4) - pkin(7);
t301 = t256 * t265;
t300 = t263 * t268;
t298 = qJD(1) * t267;
t261 = qJD(1) * t270;
t297 = qJD(6) * t265;
t296 = qJD(6) * t268;
t295 = t257 * t303;
t293 = t263 * t301;
t291 = t257 * t263 * t267;
t287 = t257 * t261;
t286 = t257 * t297;
t285 = t257 * t296;
t284 = r_i_i_C(2) * t293;
t283 = qJD(1) * t257 * t304;
t282 = r_i_i_C(2) * t285;
t281 = -qJD(6) * t256 - qJD(1);
t278 = t281 * t270;
t277 = qJD(1) * (t255 - t295);
t276 = pkin(5) * t287 + t267 * t284 + t270 * t283 + t307 * (t256 * t261 + t291);
t275 = t256 * t300 + t286;
t274 = t267 * t283 + t298 * t305 + (r_i_i_C(1) * t286 + t282) * t270 + (t307 * t298 + t314 * t299) * t256;
t273 = qJD(2) + t309 + (t305 + t310) * t263;
t272 = (-t314 * t257 + t295 - t310) * t263 + (r_i_i_C(1) * t297 + r_i_i_C(2) * t296) * t256;
t271 = -t275 * r_i_i_C(1) - t263 * t306 - t282;
t240 = t254 * qJD(3);
t233 = t268 * t311 + (t257 * t300 + t281 * t265) * t267;
t232 = t281 * t268 * t267 + (-t291 - t311) * t265;
t231 = t265 * t278 - t308 * t268;
t230 = t308 * t265 + t268 * t278;
t1 = [t231 * r_i_i_C(1) + t230 * r_i_i_C(2) - t267 * qJD(4) + t273 * t270 + (t312 * t267 + t302 * t270) * qJD(1), t261, t270 * t277 + (-t240 + t271) * t267 + t276, -t298, t271 * t267 - t287 * t303 + t276, t232 * r_i_i_C(1) - t233 * r_i_i_C(2); t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t270 * qJD(4) + t273 * t267 + (t302 * t267 - t312 * t270) * qJD(1), t298, t267 * t277 + (t240 + (-r_i_i_C(2) * t301 - t294) * t263) * t270 + t274, t261, -t270 * t284 + (-t298 * t303 - t307 * t299) * t257 + t274, -t230 * r_i_i_C(1) + t231 * r_i_i_C(2); 0, 0, t272 - t309, 0, t272, t275 * r_i_i_C(2) + (-t285 + t293) * r_i_i_C(1);];
JaD_transl  = t1;
