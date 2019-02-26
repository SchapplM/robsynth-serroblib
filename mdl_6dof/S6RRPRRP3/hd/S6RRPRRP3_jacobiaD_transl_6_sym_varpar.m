% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:10
% EndTime: 2019-02-26 21:47:11
% DurationCPUTime: 0.30s
% Computational Cost: add. (509->69), mult. (548->96), div. (0->0), fcn. (429->10), ass. (0->56)
t260 = sin(qJ(4));
t258 = qJ(4) + qJ(5);
t252 = sin(t258);
t256 = qJD(4) + qJD(5);
t296 = t252 * t256;
t298 = pkin(4) * qJD(4);
t242 = -pkin(5) * t296 - t260 * t298;
t253 = cos(t258);
t263 = cos(qJ(4));
t247 = t263 * pkin(4) + pkin(5) * t253;
t245 = pkin(3) + t247;
t246 = t260 * pkin(4) + pkin(5) * t252;
t257 = qJ(2) + pkin(10);
t250 = sin(t257);
t251 = cos(t257);
t299 = r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
t273 = t251 * t299 - sin(qJ(2)) * pkin(2);
t284 = t250 * qJD(6);
t311 = (-t245 * t250 + t273) * qJD(2) + (t246 + qJ(3) + pkin(7)) * qJD(1) + t251 * t242 + t284;
t265 = cos(qJ(1));
t291 = qJD(1) * t251;
t280 = -t256 + t291;
t262 = sin(qJ(1));
t286 = qJD(2) * t262;
t305 = -t250 * t286 + t265 * t280;
t310 = t305 * t252;
t274 = r_i_i_C(1) * t253 - r_i_i_C(2) * t252 + t245;
t267 = -t250 * t274 + t273;
t307 = -t299 * t250 - cos(qJ(2)) * pkin(2);
t301 = r_i_i_C(2) * t253;
t279 = r_i_i_C(1) * t252 + t301;
t306 = t256 * t279 - t242;
t303 = -pkin(5) - r_i_i_C(1);
t295 = t253 * t256;
t285 = qJD(2) * t265;
t269 = t250 * t285 + t262 * t280;
t281 = t251 * t256 - qJD(1);
t277 = t253 * t281;
t238 = t252 * t269 - t265 * t277;
t239 = t252 * t265 * t281 + t253 * t269;
t294 = t238 * r_i_i_C(1) + t239 * r_i_i_C(2);
t276 = t281 * t262;
t240 = t253 * t276 + t310;
t241 = t252 * t276 - t253 * t305;
t293 = -t240 * r_i_i_C(1) + t241 * r_i_i_C(2);
t290 = qJD(1) * t262;
t289 = qJD(1) * t265;
t288 = qJD(2) * t250;
t287 = qJD(2) * t251;
t278 = t246 * t291 + t242;
t243 = pkin(5) * t295 + t263 * t298;
t271 = qJD(1) * t247 - t243 * t251 + t246 * t288;
t268 = qJD(3) + t243 + (-t245 * t251 - pkin(1) + t307) * qJD(1);
t266 = qJD(6) * t251 + t306 * t250 + (-t251 * t274 + t307) * qJD(2);
t244 = t250 * r_i_i_C(2) * t296;
t1 = [t241 * r_i_i_C(1) + t240 * r_i_i_C(2) - t311 * t262 + t268 * t265, t266 * t265 - t267 * t290, t289, t262 * t278 + t265 * t271 + t294, pkin(5) * t238 + t294, -t250 * t290 + t251 * t285; -t239 * r_i_i_C(1) + t238 * r_i_i_C(2) + t268 * t262 + t311 * t265, t262 * t266 + t267 * t289, t290, t262 * t271 - t265 * t278 + t293 (-t262 * t277 - t310) * pkin(5) + t293, t250 * t289 + t251 * t286; 0, t267 * qJD(2) - t306 * t251 + t284, 0, t244 + (-r_i_i_C(1) * t295 - t243) * t250 + (-t246 - t279) * t287, t244 + t303 * t250 * t295 + (t252 * t303 - t301) * t287, t288;];
JaD_transl  = t1;
