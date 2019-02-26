% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:34:44
% EndTime: 2019-02-26 20:34:44
% DurationCPUTime: 0.30s
% Computational Cost: add. (503->63), mult. (404->95), div. (0->0), fcn. (305->11), ass. (0->54)
t271 = pkin(11) + qJ(4);
t269 = qJ(5) + t271;
t263 = sin(t269);
t275 = cos(qJ(6));
t320 = r_i_i_C(1) * t275 + pkin(5);
t325 = t263 * t320;
t264 = cos(t269);
t302 = qJD(6) * t275;
t272 = qJD(4) + qJD(5);
t274 = sin(qJ(6));
t307 = t272 * t274;
t324 = t263 * t302 + t264 * t307;
t314 = pkin(9) + r_i_i_C(3);
t322 = t314 * t264;
t323 = (-pkin(5) * t263 + t322) * t272;
t265 = sin(t271);
t309 = pkin(4) * qJD(4);
t301 = t265 * t309;
t321 = -t301 + t323;
t303 = qJD(6) * t274;
t290 = t263 * t303;
t317 = r_i_i_C(1) * t290 + t324 * r_i_i_C(2);
t284 = qJD(1) * t264 - qJD(6);
t316 = t275 * t284;
t285 = qJD(6) * t264 - qJD(1);
t296 = t263 * t307;
t315 = t285 * t275 - t296;
t313 = pkin(4) * t265;
t310 = r_i_i_C(2) * t274;
t306 = t272 * t275;
t273 = qJ(1) + pkin(10);
t266 = sin(t273);
t305 = qJD(1) * t266;
t268 = cos(t273);
t304 = qJD(1) * t268;
t300 = qJD(1) * t310;
t299 = t314 * t263;
t297 = t314 * t272;
t295 = t263 * t306;
t283 = t320 * t272;
t282 = t317 * t268 + t305 * t325;
t281 = t284 * t274;
t280 = t268 * t263 * t300 + t317 * t266 + t304 * t322;
t267 = cos(t271);
t279 = -pkin(5) * t264 - pkin(4) * t267 - cos(pkin(11)) * pkin(3) - pkin(2) - t299;
t278 = t285 * t274 + t295;
t277 = -t267 * t309 + (-t264 * t320 - t299) * t272;
t276 = (-t264 * t303 - t295) * r_i_i_C(1) + (-t264 * t302 + t296) * r_i_i_C(2) + t323;
t270 = -pkin(8) - pkin(7) - qJ(3);
t247 = t278 * t266 - t268 * t316;
t246 = t315 * t266 + t268 * t281;
t245 = t266 * t316 + t278 * t268;
t244 = t266 * t281 - t315 * t268;
t1 = [t247 * r_i_i_C(1) + t246 * r_i_i_C(2) + t268 * qJD(3) - t321 * t266 + (-cos(qJ(1)) * pkin(1) + t266 * t270 + t279 * t268) * qJD(1), 0, t304 (-t263 * t310 + t313 - t322) * t305 + t277 * t268 + t282 (-t266 * t300 - t268 * t297) * t263 + (-t268 * t283 - t314 * t305) * t264 + t282, t244 * r_i_i_C(1) + t245 * r_i_i_C(2); -t245 * r_i_i_C(1) + t244 * r_i_i_C(2) + t266 * qJD(3) + t321 * t268 + (-sin(qJ(1)) * pkin(1) - t268 * t270 + t279 * t266) * qJD(1), 0, t305 (-t313 - t325) * t304 + t277 * t266 + t280, -t266 * t264 * t283 + (-t266 * t297 - t304 * t320) * t263 + t280, -t246 * r_i_i_C(1) + t247 * r_i_i_C(2); 0, 0, 0, t276 - t301, t276 (-t264 * t306 + t290) * r_i_i_C(2) - t324 * r_i_i_C(1);];
JaD_transl  = t1;
