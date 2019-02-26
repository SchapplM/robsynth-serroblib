% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:53
% EndTime: 2019-02-26 20:05:54
% DurationCPUTime: 0.36s
% Computational Cost: add. (331->67), mult. (1044->117), div. (0->0), fcn. (1080->10), ass. (0->44)
t260 = sin(pkin(11));
t262 = cos(pkin(11));
t269 = cos(qJ(2));
t263 = cos(pkin(6));
t266 = sin(qJ(2));
t289 = t263 * t266;
t251 = t260 * t269 + t262 * t289;
t265 = sin(qJ(3));
t261 = sin(pkin(6));
t268 = cos(qJ(3));
t290 = t261 * t268;
t240 = t251 * t265 + t262 * t290;
t291 = t261 * t265;
t241 = t251 * t268 - t262 * t291;
t264 = sin(qJ(5));
t267 = cos(qJ(5));
t299 = ((t240 * t264 + t241 * t267) * r_i_i_C(1) + (t240 * t267 - t241 * t264) * r_i_i_C(2)) * qJD(5);
t272 = t260 * t289 - t262 * t269;
t243 = t260 * t291 - t268 * t272;
t274 = t260 * t290 + t265 * t272;
t298 = ((t243 * t267 - t264 * t274) * r_i_i_C(1) + (-t243 * t264 - t267 * t274) * r_i_i_C(2)) * qJD(5);
t254 = -t263 * t268 + t266 * t291;
t255 = t263 * t265 + t266 * t290;
t297 = ((t254 * t264 + t255 * t267) * r_i_i_C(1) + (t254 * t267 - t255 * t264) * r_i_i_C(2)) * qJD(5);
t275 = t267 * r_i_i_C(1) - t264 * r_i_i_C(2) + pkin(3) + pkin(4);
t276 = t264 * r_i_i_C(1) + t267 * r_i_i_C(2) + qJ(4);
t293 = t276 * t265 + t275 * t268 + pkin(2);
t288 = t263 * t269;
t287 = qJD(2) * t266;
t286 = r_i_i_C(3) + pkin(9) - pkin(8);
t284 = t262 * t288;
t283 = qJD(2) * t261 * t269;
t273 = t260 * t288 + t262 * t266;
t271 = qJD(2) * t293;
t270 = t265 * qJD(4) + ((-t264 * t268 + t265 * t267) * r_i_i_C(1) + (-t264 * t265 - t267 * t268) * r_i_i_C(2)) * qJD(5) + (-t275 * t265 + t276 * t268) * qJD(3);
t248 = t273 * qJD(2);
t246 = -qJD(2) * t284 + t260 * t287;
t245 = -t254 * qJD(3) + t268 * t283;
t244 = t255 * qJD(3) + t265 * t283;
t239 = t274 * qJD(3) - t248 * t268;
t238 = t243 * qJD(3) - t248 * t265;
t237 = -t240 * qJD(3) - t246 * t268;
t236 = t241 * qJD(3) - t246 * t265;
t1 = [0, t286 * t248 - t270 * t273 + t272 * t271, t243 * qJD(4) - t275 * t238 + t276 * t239 + t298, t238 (t238 * t267 - t239 * t264) * r_i_i_C(1) + (-t238 * t264 - t239 * t267) * r_i_i_C(2) - t298, 0; 0, t286 * t246 - t251 * t271 + t270 * (-t260 * t266 + t284) t241 * qJD(4) - t275 * t236 + t276 * t237 + t299, t236 (t236 * t267 - t237 * t264) * r_i_i_C(1) + (-t236 * t264 - t237 * t267) * r_i_i_C(2) - t299, 0; 0 (-t293 * t287 + (-t286 * qJD(2) + t270) * t269) * t261, t255 * qJD(4) - t275 * t244 + t276 * t245 + t297, t244 (t244 * t267 - t245 * t264) * r_i_i_C(1) + (-t244 * t264 - t245 * t267) * r_i_i_C(2) - t297, 0;];
JaD_transl  = t1;
