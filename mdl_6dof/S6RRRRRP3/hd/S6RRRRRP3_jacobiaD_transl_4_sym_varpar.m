% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:40:42
% EndTime: 2019-02-26 22:40:43
% DurationCPUTime: 0.29s
% Computational Cost: add. (270->58), mult. (394->95), div. (0->0), fcn. (297->8), ass. (0->54)
t258 = qJ(2) + qJ(3);
t255 = sin(t258);
t262 = cos(qJ(4));
t310 = r_i_i_C(1) * t262 + pkin(3);
t313 = t255 * t310;
t259 = sin(qJ(4));
t293 = qJD(4) * t262;
t256 = cos(t258);
t257 = qJD(2) + qJD(3);
t301 = t256 * t257;
t312 = t255 * t293 + t259 * t301;
t307 = pkin(9) + r_i_i_C(3);
t288 = t307 * t256;
t260 = sin(qJ(2));
t302 = pkin(2) * qJD(2);
t290 = t260 * t302;
t305 = pkin(3) * t255;
t311 = -t290 + (t288 - t305) * t257;
t294 = qJD(4) * t259;
t281 = t255 * t294;
t308 = r_i_i_C(1) * t281 + t312 * r_i_i_C(2);
t306 = pkin(2) * t260;
t303 = r_i_i_C(2) * t259;
t261 = sin(qJ(1));
t300 = t257 * t261;
t299 = t257 * t262;
t264 = cos(qJ(1));
t298 = t257 * t264;
t297 = t262 * t264;
t296 = qJD(1) * t261;
t295 = qJD(1) * t264;
t292 = t255 * t303;
t291 = qJD(1) * t303;
t289 = t307 * t255;
t287 = t307 * t261;
t286 = t255 * t299;
t276 = qJD(4) * t256 - qJD(1);
t275 = qJD(1) * t256 - qJD(4);
t274 = t310 * t256;
t273 = t310 * t264;
t272 = t308 * t264 + t296 * t313;
t271 = t276 * t259;
t270 = t264 * t255 * t291 + t308 * t261 + t295 * t288;
t263 = cos(qJ(2));
t269 = -t263 * pkin(2) - pkin(3) * t256 - pkin(1) - t289;
t268 = t255 * t298 + t275 * t261;
t267 = -t263 * t302 + (-t274 - t289) * t257;
t266 = -t256 * r_i_i_C(2) * t293 + (-t256 * t294 - t286) * r_i_i_C(1) + t307 * t301 + (-t305 + t292) * t257;
t265 = -pkin(8) - pkin(7);
t239 = -t275 * t297 + (t271 + t286) * t261;
t238 = t276 * t262 * t261 + (-t255 * t300 + t275 * t264) * t259;
t237 = t268 * t262 + t264 * t271;
t236 = t268 * t259 - t276 * t297;
t1 = [t239 * r_i_i_C(1) + t238 * r_i_i_C(2) - t311 * t261 + (t261 * t265 + t269 * t264) * qJD(1) (-t288 - t292 + t306) * t296 + t267 * t264 + t272 (-t261 * t291 - t307 * t298) * t255 + (-qJD(1) * t287 - t257 * t273) * t256 + t272, t236 * r_i_i_C(1) + t237 * r_i_i_C(2), 0, 0; -t237 * r_i_i_C(1) + t236 * r_i_i_C(2) + t311 * t264 + (t269 * t261 - t264 * t265) * qJD(1) (-t306 - t313) * t295 + t267 * t261 + t270, -t274 * t300 + (-qJD(1) * t273 - t257 * t287) * t255 + t270, -t238 * r_i_i_C(1) + t239 * r_i_i_C(2), 0, 0; 0, t266 - t290, t266 (-t256 * t299 + t281) * r_i_i_C(2) - t312 * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
