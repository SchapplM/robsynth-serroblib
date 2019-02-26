% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:50
% EndTime: 2019-02-26 21:21:50
% DurationCPUTime: 0.43s
% Computational Cost: add. (387->70), mult. (685->110), div. (0->0), fcn. (619->10), ass. (0->50)
t261 = qJ(2) + pkin(9);
t259 = sin(t261);
t260 = cos(t261);
t289 = -r_i_i_C(3) - pkin(8) + qJ(4);
t277 = t289 * t260 - sin(qJ(2)) * pkin(2);
t291 = t259 * qJD(4);
t310 = (-t259 * pkin(3) + t277) * qJD(2) - qJD(1) * (-qJ(3) - pkin(7)) + t291;
t262 = sin(pkin(10));
t263 = cos(pkin(10));
t265 = sin(qJ(6));
t268 = cos(qJ(6));
t304 = pkin(4) + pkin(5);
t278 = t268 * r_i_i_C(1) - t265 * r_i_i_C(2) + t304;
t279 = t265 * r_i_i_C(1) + t268 * r_i_i_C(2) + qJ(5);
t305 = t279 * t262 + t278 * t263 + pkin(3);
t308 = t305 * t259 - t277;
t307 = -t289 * t259 - cos(qJ(2)) * pkin(2);
t281 = t262 * t265 + t263 * t268;
t282 = t262 * t268 - t263 * t265;
t274 = t282 * r_i_i_C(1) - t281 * r_i_i_C(2);
t306 = t262 * qJD(5) + t274 * qJD(6);
t267 = sin(qJ(1));
t300 = t267 * t262;
t299 = t267 * t263;
t270 = cos(qJ(1));
t298 = t270 * t262;
t297 = t270 * t263;
t296 = qJD(1) * t267;
t295 = qJD(1) * t270;
t294 = qJD(2) * t260;
t293 = qJD(2) * t267;
t292 = qJD(2) * t270;
t288 = t260 * t298;
t287 = t259 * t293;
t286 = t259 * t292;
t250 = t260 * t300 + t297;
t251 = t260 * t299 - t298;
t284 = -t250 * t268 + t251 * t265;
t283 = t250 * t265 + t251 * t268;
t253 = t260 * t297 + t300;
t275 = qJD(3) + (-pkin(3) * t260 - pkin(1) + t307) * qJD(1);
t271 = t260 * qJD(4) - t306 * t259 + (-t260 * t305 + t307) * qJD(2);
t252 = t288 - t299;
t249 = t253 * qJD(1) - t263 * t287;
t248 = qJD(1) * t288 - t262 * t287 - t263 * t296;
t247 = -t251 * qJD(1) - t263 * t286;
t246 = t250 * qJD(1) + t262 * t286;
t245 = -t246 * t265 + t247 * t268 + (t252 * t268 - t253 * t265) * qJD(6);
t244 = -t246 * t268 - t247 * t265 + (-t252 * t265 - t253 * t268) * qJD(6);
t1 = [-t250 * qJD(5) - t279 * t248 - t278 * t249 + (t284 * r_i_i_C(1) + t283 * r_i_i_C(2)) * qJD(6) + t275 * t270 - t310 * t267, t271 * t270 + t308 * t296, t295, -t259 * t296 + t260 * t292, -t246, t244 * r_i_i_C(1) - t245 * r_i_i_C(2); t245 * r_i_i_C(1) + t244 * r_i_i_C(2) - t246 * qJ(5) + t252 * qJD(5) + t304 * t247 + t275 * t267 + t310 * t270, t271 * t267 - t295 * t308, t296, t259 * t295 + t260 * t293, t248 (t248 * t268 - t249 * t265) * r_i_i_C(1) + (-t248 * t265 - t249 * t268) * r_i_i_C(2) + (-t283 * r_i_i_C(1) + t284 * r_i_i_C(2)) * qJD(6); 0, -qJD(2) * t308 + t306 * t260 + t291, 0, qJD(2) * t259, t262 * t294 (-t281 * r_i_i_C(1) - t282 * r_i_i_C(2)) * t259 * qJD(6) + t274 * t294;];
JaD_transl  = t1;
