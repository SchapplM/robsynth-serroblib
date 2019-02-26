% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:37:05
% EndTime: 2019-02-26 20:37:05
% DurationCPUTime: 0.24s
% Computational Cost: add. (279->60), mult. (414->86), div. (0->0), fcn. (311->8), ass. (0->50)
t307 = pkin(9) + r_i_i_C(3);
t263 = qJ(4) + qJ(5);
t261 = cos(t263);
t267 = cos(qJ(6));
t312 = (r_i_i_C(1) * t267 + pkin(5)) * t261;
t262 = qJD(4) + qJD(5);
t260 = sin(t263);
t292 = t307 * t260;
t268 = cos(qJ(4));
t303 = pkin(4) * qJD(4);
t293 = t268 * t303;
t311 = (pkin(5) * t261 + t292) * t262 + qJD(3) + t293;
t306 = pkin(4) * t268;
t302 = t260 * t262;
t301 = t261 * t262;
t300 = t262 * t267;
t269 = cos(qJ(1));
t299 = t267 * t269;
t298 = qJ(2) - pkin(8) - pkin(7);
t266 = sin(qJ(1));
t297 = qJD(1) * t266;
t259 = qJD(1) * t269;
t264 = sin(qJ(6));
t296 = qJD(6) * t264;
t295 = qJD(6) * t267;
t294 = r_i_i_C(2) * t261 * t264;
t291 = t264 * t302;
t290 = t266 * t301;
t289 = t269 * t301;
t286 = t261 * t295;
t285 = r_i_i_C(2) * t291;
t284 = qJD(1) * t294;
t283 = qJD(6) * t260 + qJD(1);
t282 = qJD(1) * t260 + qJD(6);
t281 = t266 * t284 + t269 * t285 + t307 * t289;
t279 = t283 * t264;
t278 = t259 * t312 + t266 * t285 + t307 * (t260 * t259 + t290);
t277 = t260 * t300 + t261 * t296;
t276 = -t292 - t312;
t275 = t282 * t266 - t289;
t265 = sin(qJ(4));
t274 = -pkin(4) * t265 - pkin(5) * t260 + t307 * t261 - pkin(1) - qJ(3);
t273 = (t276 + t294) * t262 + (r_i_i_C(1) * t296 + r_i_i_C(2) * t295) * t260;
t272 = -pkin(5) * t302 - t277 * r_i_i_C(1) - r_i_i_C(2) * t286;
t271 = -t265 * t303 + t272;
t244 = -t282 * t299 + (-t261 * t300 + t279) * t266;
t243 = t283 * t267 * t266 + (t282 * t269 + t290) * t264;
t242 = t275 * t267 + t269 * t279;
t241 = t275 * t264 - t283 * t299;
t1 = [t244 * r_i_i_C(1) + t243 * r_i_i_C(2) + t269 * qJD(2) - t311 * t266 + (-t266 * t298 + t269 * t274) * qJD(1), t259, -t297, t271 * t269 + (t276 - t306) * t297 + t281, t269 * t272 + t276 * t297 + t281, t241 * r_i_i_C(1) + t242 * r_i_i_C(2); -t242 * r_i_i_C(1) + t241 * r_i_i_C(2) + t266 * qJD(2) + t311 * t269 + (t266 * t274 + t269 * t298) * qJD(1), t297, t259 (-t294 + t306) * t259 + t271 * t266 + t278, t266 * t272 - t269 * t284 + t278, -t243 * r_i_i_C(1) + t244 * r_i_i_C(2); 0, 0, 0, t273 - t293, t273, t277 * r_i_i_C(2) + (-t286 + t291) * r_i_i_C(1);];
JaD_transl  = t1;
