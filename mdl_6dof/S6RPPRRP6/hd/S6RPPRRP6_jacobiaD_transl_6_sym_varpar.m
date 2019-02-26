% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_jacobiaD_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:33:09
% EndTime: 2019-02-26 20:33:10
% DurationCPUTime: 0.22s
% Computational Cost: add. (173->59), mult. (526->86), div. (0->0), fcn. (439->6), ass. (0->38)
t248 = cos(qJ(4));
t244 = sin(qJ(5));
t247 = cos(qJ(5));
t272 = r_i_i_C(3) + qJ(6);
t275 = -r_i_i_C(1) - pkin(5);
t277 = t272 * t244 - t275 * t247 + pkin(4);
t245 = sin(qJ(4));
t276 = pkin(8) + r_i_i_C(2);
t279 = t276 * t245;
t282 = t277 * t248 + t279;
t263 = qJD(6) * t244;
t281 = (pkin(4) * t248 + t279) * qJD(4) + t245 * t263 + qJD(3);
t250 = t263 + (t275 * t244 + t272 * t247) * qJD(5);
t273 = pkin(7) - qJ(2);
t246 = sin(qJ(1));
t271 = t246 * t244;
t270 = t246 * t247;
t249 = cos(qJ(1));
t269 = t249 * t244;
t268 = qJD(1) * t246;
t243 = qJD(1) * t249;
t267 = qJD(4) * t245;
t266 = qJD(4) * t249;
t265 = qJD(5) * t245;
t264 = qJD(5) * t248;
t260 = t276 * t248;
t259 = t249 * t245 * t247;
t258 = t248 * t266;
t257 = qJD(1) * t245 + qJD(5);
t256 = t247 * qJD(6) + qJD(2);
t254 = t245 * t270 + t269;
t253 = -pkin(4) * t245 - pkin(1) - qJ(3) + t260;
t251 = t246 * qJD(4) * t248 + t257 * t249;
t236 = -t244 * t268 + t251 * t247 - t265 * t271;
t235 = (qJD(1) + t265) * t270 + t251 * t244;
t234 = -t247 * t258 + (t245 * t269 + t270) * qJD(5) + t254 * qJD(1);
t233 = -qJD(5) * t259 - t247 * t243 - t244 * t258 + t257 * t271;
t1 = [t256 * t249 + t275 * t236 - t272 * t235 - t281 * t246 + (t273 * t246 + t253 * t249) * qJD(1), t243, -t268 (-t266 * t277 - t276 * t268) * t245 + (-t277 * t268 + (t276 * qJD(4) + t250) * t249) * t248 -(-t259 + t271) * qJD(6) - t272 * t234 - t275 * t233, -t233; t256 * t246 + t275 * t234 - t272 * t233 + t281 * t249 + (t253 * t246 - t273 * t249) * qJD(1), t268, t243, t282 * t243 + (t250 * t248 + (-t245 * t277 + t260) * qJD(4)) * t246, t254 * qJD(6) + t275 * t235 + t272 * t236, t235; 0, 0, 0, -t282 * qJD(4) - t250 * t245 (-t272 * t264 - t275 * t267) * t244 + (-t272 * t267 + (t275 * qJD(5) + qJD(6)) * t248) * t247, -t244 * t267 + t247 * t264;];
JaD_transl  = t1;
