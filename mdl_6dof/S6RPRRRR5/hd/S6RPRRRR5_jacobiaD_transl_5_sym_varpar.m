% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR5
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

function JaD_transl = S6RPRRRR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:17:01
% EndTime: 2019-02-26 21:17:01
% DurationCPUTime: 0.29s
% Computational Cost: add. (387->62), mult. (400->98), div. (0->0), fcn. (303->9), ass. (0->55)
t262 = pkin(11) + qJ(3);
t260 = qJ(4) + t262;
t256 = sin(t260);
t266 = cos(qJ(5));
t312 = r_i_i_C(1) * t266 + pkin(4);
t315 = t256 * t312;
t264 = sin(qJ(5));
t295 = qJD(5) * t266;
t257 = cos(t260);
t263 = qJD(3) + qJD(4);
t303 = t257 * t263;
t314 = t256 * t295 + t264 * t303;
t309 = pkin(9) + r_i_i_C(3);
t290 = t309 * t257;
t258 = sin(t262);
t304 = pkin(3) * qJD(3);
t293 = t258 * t304;
t307 = pkin(4) * t256;
t313 = -t293 + (t290 - t307) * t263;
t296 = qJD(5) * t264;
t283 = t256 * t296;
t310 = r_i_i_C(1) * t283 + t314 * r_i_i_C(2);
t308 = pkin(3) * t258;
t305 = r_i_i_C(2) * t264;
t265 = sin(qJ(1));
t302 = t263 * t265;
t301 = t263 * t266;
t267 = cos(qJ(1));
t300 = t263 * t267;
t299 = t266 * t267;
t298 = qJD(1) * t265;
t297 = qJD(1) * t267;
t294 = t256 * t305;
t292 = qJD(1) * t305;
t291 = t309 * t256;
t289 = t309 * t265;
t288 = t256 * t301;
t278 = qJD(5) * t257 - qJD(1);
t277 = qJD(1) * t257 - qJD(5);
t276 = t312 * t257;
t275 = t312 * t267;
t274 = t310 * t267 + t298 * t315;
t273 = t278 * t264;
t272 = t267 * t256 * t292 + t310 * t265 + t297 * t290;
t259 = cos(t262);
t271 = -pkin(4) * t257 - pkin(3) * t259 - cos(pkin(11)) * pkin(2) - pkin(1) - t291;
t270 = t256 * t300 + t277 * t265;
t269 = -t259 * t304 + (-t276 - t291) * t263;
t268 = -t257 * r_i_i_C(2) * t295 + (-t257 * t296 - t288) * r_i_i_C(1) + t309 * t303 + (-t307 + t294) * t263;
t261 = -pkin(8) - pkin(7) - qJ(2);
t240 = -t277 * t299 + (t273 + t288) * t265;
t239 = t278 * t266 * t265 + (-t256 * t302 + t277 * t267) * t264;
t238 = t270 * t266 + t267 * t273;
t237 = t270 * t264 - t278 * t299;
t1 = [t240 * r_i_i_C(1) + t239 * r_i_i_C(2) + t267 * qJD(2) - t313 * t265 + (t261 * t265 + t271 * t267) * qJD(1), t297 (-t290 - t294 + t308) * t298 + t269 * t267 + t274 (-t265 * t292 - t309 * t300) * t256 + (-qJD(1) * t289 - t263 * t275) * t257 + t274, t237 * r_i_i_C(1) + t238 * r_i_i_C(2), 0; -t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t265 * qJD(2) + t313 * t267 + (-t261 * t267 + t265 * t271) * qJD(1), t298 (-t308 - t315) * t297 + t269 * t265 + t272, -t276 * t302 + (-qJD(1) * t275 - t263 * t289) * t256 + t272, -t239 * r_i_i_C(1) + t240 * r_i_i_C(2), 0; 0, 0, t268 - t293, t268 (-t257 * t301 + t283) * r_i_i_C(2) - t314 * r_i_i_C(1), 0;];
JaD_transl  = t1;
