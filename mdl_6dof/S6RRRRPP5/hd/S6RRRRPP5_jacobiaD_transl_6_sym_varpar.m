% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:27:40
% EndTime: 2019-02-26 22:27:40
% DurationCPUTime: 0.36s
% Computational Cost: add. (641->73), mult. (908->105), div. (0->0), fcn. (763->8), ass. (0->58)
t268 = qJ(3) + qJ(4);
t266 = cos(t268);
t317 = r_i_i_C(2) + qJ(5);
t328 = t317 * t266;
t270 = sin(qJ(2));
t273 = cos(qJ(2));
t272 = cos(qJ(3));
t264 = pkin(3) * t272 + pkin(2);
t265 = sin(t268);
t303 = -r_i_i_C(1) - pkin(5) - pkin(4);
t291 = t303 * t266;
t284 = -t317 * t265 + t291;
t283 = -t264 + t284;
t302 = r_i_i_C(3) + qJ(6) - pkin(9) - pkin(8);
t327 = -t283 * t270 + t302 * t273;
t267 = qJD(3) + qJD(4);
t269 = sin(qJ(3));
t316 = pkin(3) * qJD(3);
t279 = t302 * qJD(2) - qJD(5) * t265 + t269 * t316;
t292 = t303 * t265;
t326 = (-t292 - t328) * t267 + t279;
t318 = pkin(3) * t269;
t325 = t279 * t273 - (pkin(7) + t318) * qJD(1) + (qJD(2) * t264 + qJD(6)) * t270;
t274 = cos(qJ(1));
t315 = t266 * t274;
t314 = t267 * t270;
t271 = sin(qJ(1));
t313 = t267 * t271;
t312 = t267 * t274;
t311 = t271 * t273;
t310 = qJD(1) * t271;
t309 = qJD(1) * t273;
t308 = qJD(1) * t274;
t307 = qJD(2) * t270;
t306 = qJD(2) * t273;
t305 = qJD(2) * t274;
t304 = t266 * qJD(5);
t301 = t272 * t316;
t300 = t265 * t313;
t299 = t265 * t312;
t298 = t266 * t312;
t297 = t270 * t304 + t306 * t328;
t294 = t271 * t307;
t293 = t270 * t305;
t290 = -qJD(3) + t309;
t288 = (-qJD(3) * t273 + qJD(1)) * t272;
t286 = t265 * t271 + t273 * t315;
t285 = t265 * t308 + t266 * t313;
t245 = t265 * t293 - t273 * t298 - t300 + (t265 * t311 + t315) * qJD(1);
t246 = t273 * t299 + (t271 * t309 + t293) * t266 - t285;
t282 = t286 * qJD(5) - t303 * t245 - t317 * t246;
t247 = -t265 * t294 - t266 * t310 + t285 * t273 - t299;
t248 = t286 * qJD(1) - t266 * t294 - t273 * t300 - t298;
t281 = -(t265 * t274 - t266 * t311) * qJD(5) + t317 * t248 + t303 * t247;
t278 = t283 * qJD(2) - qJD(6);
t277 = t301 - t304 + (-t264 * t273 + t302 * t270 - pkin(1)) * qJD(1);
t276 = t326 * t270 + t278 * t273;
t1 = [-t317 * t247 + t303 * t248 + t325 * t271 + t277 * t274, t276 * t274 + t327 * t310 (t274 * t288 + (t290 * t271 + t293) * t269) * pkin(3) + t282, t282, -t245, t270 * t310 - t273 * t305; -t317 * t245 + t303 * t246 + t277 * t271 - t325 * t274, t276 * t271 - t327 * t308 (t271 * t288 + (-t290 * t274 + t294) * t269) * pkin(3) + t281, t281, t247, -t270 * t308 - t271 * t306; 0, t278 * t270 - t326 * t273 (t292 - t318) * t306 + (t284 * t267 - t301) * t270 + t297, t291 * t314 + (t303 * t306 - t317 * t314) * t265 + t297, t265 * t306 + t266 * t314, -t307;];
JaD_transl  = t1;
