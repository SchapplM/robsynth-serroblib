% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:27:28
% EndTime: 2019-02-26 21:27:28
% DurationCPUTime: 0.30s
% Computational Cost: add. (353->59), mult. (613->88), div. (0->0), fcn. (510->8), ass. (0->42)
t247 = sin(qJ(2));
t249 = cos(qJ(2));
t244 = pkin(9) + qJ(5);
t243 = cos(t244);
t261 = qJD(6) * t243 - qJD(3);
t264 = pkin(4) * sin(pkin(9)) + qJ(3);
t269 = t249 * qJD(4);
t268 = pkin(2) + r_i_i_C(2) + pkin(8) + qJ(4);
t284 = t268 * t247;
t289 = (-t264 * t249 + t284) * qJD(2) - (pkin(7) + cos(pkin(9)) * pkin(4) + pkin(3)) * qJD(1) + t261 * t247 - t269;
t242 = sin(t244);
t279 = r_i_i_C(3) + qJ(6);
t281 = -r_i_i_C(1) - pkin(5);
t255 = t281 * t242 + t279 * t243 - t264;
t287 = t255 * t249 + t284;
t248 = sin(qJ(1));
t272 = qJD(2) * t249;
t266 = t248 * t272;
t250 = cos(qJ(1));
t274 = qJD(1) * t250;
t283 = t247 * t274 + t266;
t278 = t248 * t243;
t277 = t250 * t242;
t276 = t250 * t243;
t275 = qJD(1) * t248;
t273 = qJD(2) * t247;
t271 = qJD(2) * t250;
t270 = qJD(5) * t249;
t265 = t249 * t271;
t263 = qJD(5) * t247 + qJD(1);
t262 = qJD(1) * t247 + qJD(5);
t259 = t268 * t249;
t257 = t263 * t248;
t256 = t247 * t277 + t278;
t254 = (t279 * t242 - t281 * t243) * qJD(5) - t261;
t253 = t242 * qJD(6) + (-t264 * t247 - pkin(1) - t259) * qJD(1);
t251 = -qJD(4) * t247 + t254 * t249 + (t255 * t247 - t259) * qJD(2);
t236 = t263 * t276 + (-t262 * t248 + t265) * t242;
t235 = -t243 * t265 + t256 * qJD(5) + (t247 * t278 + t277) * qJD(1);
t234 = t243 * t257 + (t262 * t250 + t266) * t242;
t233 = -qJD(5) * t276 + t242 * t257 - t283 * t243;
t1 = [-t279 * t233 + t281 * t234 + t289 * t248 + t253 * t250, t251 * t250 + t287 * t275, -t247 * t275 + t265, -t247 * t271 - t249 * t275, t256 * qJD(6) + t281 * t235 + t279 * t236, t235; t279 * t235 - t281 * t236 + t253 * t248 - t289 * t250, t251 * t248 - t274 * t287, t283, -t248 * t273 + t249 * t274 -(-t248 * t247 * t242 + t276) * qJD(6) + t279 * t234 + t281 * t233, t233; 0, -qJD(2) * t287 + t254 * t247 + t269, t273, t272 (-t279 * t270 - t281 * t273) * t243 + (t279 * t273 + (-t281 * qJD(5) - qJD(6)) * t249) * t242, -t242 * t270 - t243 * t273;];
JaD_transl  = t1;
