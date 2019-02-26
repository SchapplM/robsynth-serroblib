% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:06
% EndTime: 2019-02-26 22:11:06
% DurationCPUTime: 0.30s
% Computational Cost: add. (334->60), mult. (461->87), div. (0->0), fcn. (348->8), ass. (0->55)
t258 = qJ(2) + qJ(3);
t256 = cos(t258);
t259 = sin(qJ(5));
t262 = cos(qJ(5));
t310 = r_i_i_C(1) * t259 + r_i_i_C(2) * t262 + qJ(4);
t314 = t256 * t310;
t257 = qJD(2) + qJD(3);
t301 = t256 * t257;
t250 = qJ(4) * t301;
t255 = sin(t258);
t292 = pkin(3) + pkin(9) + r_i_i_C(3);
t281 = t292 * t257;
t260 = sin(qJ(2));
t302 = pkin(2) * qJD(2);
t291 = t260 * t302;
t313 = (qJD(4) - t281) * t255 + (pkin(4) + pkin(8) + pkin(7)) * qJD(1) + t250 - t291;
t293 = qJD(5) * t262;
t284 = t256 * t293;
t312 = r_i_i_C(1) * t284 + qJD(4) * t256;
t279 = qJD(5) * t255 + qJD(1);
t309 = t262 * t279;
t308 = t279 * t259;
t306 = pkin(2) * t260;
t300 = t257 * t259;
t299 = t257 * t262;
t264 = cos(qJ(1));
t298 = t257 * t264;
t261 = sin(qJ(1));
t297 = qJD(1) * t261;
t296 = qJD(1) * t264;
t294 = qJD(5) * t259;
t290 = t256 * t299;
t289 = t261 * t301;
t288 = t256 * t298;
t286 = t255 * t297;
t285 = t256 * t294;
t283 = t292 * t255;
t282 = t292 * t256;
t278 = -qJD(1) * t255 - qJD(5);
t277 = t261 * t312 + t296 * t314;
t276 = t264 * t312 + t292 * t286;
t275 = t278 * t264;
t272 = t310 * t255;
t271 = -r_i_i_C(2) * t294 - t281;
t263 = cos(qJ(2));
t270 = qJD(1) * (-t263 * pkin(2) - qJ(4) * t255 - pkin(1) - t282);
t269 = t278 * t261 + t288;
t268 = t256 * r_i_i_C(1) * t300 + r_i_i_C(2) * t290 + t250 + (r_i_i_C(1) * t293 + qJD(4) + t271) * t255;
t267 = -r_i_i_C(2) * t285 + (-t282 - t272) * t257;
t266 = -t263 * t302 + t267;
t238 = t269 * t259 + t264 * t309;
t237 = t269 * t262 - t264 * t308;
t236 = -t261 * t309 + (t275 - t289) * t259;
t235 = t262 * t275 + (-t290 + t308) * t261;
t1 = [t236 * r_i_i_C(1) + t235 * r_i_i_C(2) - t313 * t261 + t264 * t270 (t306 - t314) * t297 + t266 * t264 + t276, -t272 * t298 + (t271 * t264 - t297 * t310) * t256 + t276, -t286 + t288, t237 * r_i_i_C(1) - t238 * r_i_i_C(2), 0; t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t261 * t270 + t313 * t264 (-t283 - t306) * t296 + t266 * t261 + t277, t267 * t261 - t283 * t296 + t277, t255 * t296 + t289, -t235 * r_i_i_C(1) + t236 * r_i_i_C(2), 0; 0, t268 - t291, t268, t257 * t255 (-t255 * t300 + t284) * r_i_i_C(2) + (t255 * t299 + t285) * r_i_i_C(1), 0;];
JaD_transl  = t1;
