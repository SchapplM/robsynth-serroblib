% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP2
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
% Datum: 2019-02-26 21:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:46:36
% EndTime: 2019-02-26 21:46:36
% DurationCPUTime: 0.30s
% Computational Cost: add. (394->63), mult. (412->98), div. (0->0), fcn. (310->10), ass. (0->53)
t268 = qJ(2) + pkin(10);
t265 = qJ(4) + t268;
t261 = sin(t265);
t272 = cos(qJ(5));
t317 = r_i_i_C(1) * t272 + pkin(4);
t320 = t261 * t317;
t269 = sin(qJ(5));
t302 = qJD(5) * t272;
t262 = cos(t265);
t267 = qJD(2) + qJD(4);
t310 = t262 * t267;
t319 = t261 * t302 + t269 * t310;
t314 = pkin(9) + r_i_i_C(3);
t298 = t314 * t262;
t258 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t268);
t251 = t258 * qJD(2);
t313 = pkin(4) * t261;
t318 = t251 + (t298 - t313) * t267;
t303 = qJD(5) * t269;
t291 = t261 * t303;
t315 = r_i_i_C(1) * t291 + t319 * r_i_i_C(2);
t311 = r_i_i_C(2) * t269;
t271 = sin(qJ(1));
t309 = t267 * t271;
t308 = t267 * t272;
t274 = cos(qJ(1));
t307 = t267 * t274;
t306 = t272 * t274;
t305 = qJD(1) * t271;
t304 = qJD(1) * t274;
t301 = t261 * t311;
t300 = qJD(1) * t311;
t299 = t314 * t261;
t297 = t314 * t271;
t296 = t261 * t308;
t286 = qJD(5) * t262 - qJD(1);
t285 = qJD(1) * t262 - qJD(5);
t284 = t317 * t262;
t283 = t317 * t274;
t282 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t268);
t281 = t315 * t274 + t305 * t320;
t280 = t286 * t269;
t279 = t274 * t261 * t300 + t315 * t271 + t304 * t298;
t278 = -pkin(4) * t262 - pkin(1) + t282 - t299;
t277 = t261 * t307 + t285 * t271;
t276 = t282 * qJD(2) + (-t284 - t299) * t267;
t275 = -t262 * r_i_i_C(2) * t302 + (-t262 * t303 - t296) * r_i_i_C(1) + t314 * t310 + (-t313 + t301) * t267;
t266 = -pkin(8) - qJ(3) - pkin(7);
t242 = -t285 * t306 + (t280 + t296) * t271;
t241 = t286 * t272 * t271 + (-t261 * t309 + t285 * t274) * t269;
t240 = t277 * t272 + t274 * t280;
t239 = t277 * t269 - t286 * t306;
t1 = [t242 * r_i_i_C(1) + t241 * r_i_i_C(2) + t274 * qJD(3) - t318 * t271 + (t266 * t271 + t278 * t274) * qJD(1) (-t258 - t298 - t301) * t305 + t276 * t274 + t281, t304 (-t271 * t300 - t314 * t307) * t261 + (-qJD(1) * t297 - t267 * t283) * t262 + t281, t239 * r_i_i_C(1) + t240 * r_i_i_C(2), 0; -t240 * r_i_i_C(1) + t239 * r_i_i_C(2) + t271 * qJD(3) + t318 * t274 + (-t266 * t274 + t278 * t271) * qJD(1) (t258 - t320) * t304 + t276 * t271 + t279, t305, -t284 * t309 + (-qJD(1) * t283 - t267 * t297) * t261 + t279, -t241 * r_i_i_C(1) + t242 * r_i_i_C(2), 0; 0, t251 + t275, 0, t275 (-t262 * t308 + t291) * r_i_i_C(2) - t319 * r_i_i_C(1), 0;];
JaD_transl  = t1;
