% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:43
% EndTime: 2019-02-26 22:31:43
% DurationCPUTime: 0.35s
% Computational Cost: add. (708->74), mult. (610->96), div. (0->0), fcn. (451->10), ass. (0->64)
t273 = qJ(2) + qJ(3);
t270 = qJ(4) + t273;
t266 = cos(t270);
t274 = sin(qJ(6));
t277 = cos(qJ(6));
t330 = r_i_i_C(1) * t274 + r_i_i_C(2) * t277 + qJ(5);
t288 = t330 * t266;
t271 = qJD(2) + qJD(3);
t267 = qJD(4) + t271;
t320 = t266 * t267;
t254 = qJ(5) * t320;
t265 = sin(t270);
t311 = pkin(4) + pkin(10) + r_i_i_C(3);
t298 = t311 * t267;
t275 = sin(qJ(2));
t321 = pkin(2) * qJD(2);
t308 = t275 * t321;
t268 = sin(t273);
t325 = pkin(3) * t271;
t310 = t268 * t325;
t333 = (qJD(5) - t298) * t265 + (pkin(5) + pkin(9) + pkin(8) + pkin(7)) * qJD(1) + t254 - t308 - t310;
t312 = qJD(6) * t277;
t301 = t266 * t312;
t332 = r_i_i_C(1) * t301 + qJD(5) * t266;
t296 = qJD(6) * t265 + qJD(1);
t329 = t277 * t296;
t328 = t296 * t274;
t269 = cos(t273);
t289 = t330 * t265;
t299 = t311 * t266;
t313 = qJD(6) * t274;
t302 = t266 * t313;
t282 = -r_i_i_C(2) * t302 + (-t299 - t289) * t267;
t280 = -t269 * t325 + t282;
t326 = pkin(3) * t268;
t319 = t267 * t274;
t318 = t267 * t277;
t279 = cos(qJ(1));
t317 = t267 * t279;
t276 = sin(qJ(1));
t316 = qJD(1) * t276;
t315 = qJD(1) * t279;
t307 = t266 * t318;
t306 = t276 * t320;
t305 = t266 * t317;
t303 = t265 * t316;
t300 = t311 * t265;
t295 = -qJD(1) * t265 - qJD(6);
t293 = t332 * t276 + t315 * t288;
t292 = t332 * t279 + t311 * t303;
t291 = t295 * t279;
t287 = -r_i_i_C(2) * t313 - t298;
t278 = cos(qJ(2));
t286 = qJD(1) * (-t278 * pkin(2) - pkin(3) * t269 - qJ(5) * t265 - pkin(1) - t299);
t285 = t295 * t276 + t305;
t284 = t266 * r_i_i_C(1) * t319 + r_i_i_C(2) * t307 + t254 + (r_i_i_C(1) * t312 + qJD(5) + t287) * t265;
t283 = t284 - t310;
t281 = -t278 * t321 + t280;
t261 = -t275 * pkin(2) - t326;
t245 = t285 * t274 + t279 * t329;
t244 = t285 * t277 - t279 * t328;
t243 = -t276 * t329 + (t291 - t306) * t274;
t242 = t277 * t291 + (-t307 + t328) * t276;
t1 = [t243 * r_i_i_C(1) + t242 * r_i_i_C(2) - t333 * t276 + t279 * t286 (-t261 - t288) * t316 + t281 * t279 + t292 (-t288 + t326) * t316 + t280 * t279 + t292, -t289 * t317 + (t287 * t279 - t316 * t330) * t266 + t292, -t303 + t305, r_i_i_C(1) * t244 - r_i_i_C(2) * t245; t245 * r_i_i_C(1) + t244 * r_i_i_C(2) + t276 * t286 + t333 * t279 (t261 - t300) * t315 + t281 * t276 + t293 (-t300 - t326) * t315 + t280 * t276 + t293, t282 * t276 - t300 * t315 + t293, t265 * t315 + t306, -r_i_i_C(1) * t242 + r_i_i_C(2) * t243; 0, t283 - t308, t283, t284, t267 * t265 (-t265 * t319 + t301) * r_i_i_C(2) + (t265 * t318 + t302) * r_i_i_C(1);];
JaD_transl  = t1;
