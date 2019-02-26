% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP8
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
% Datum: 2019-02-26 21:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:54
% EndTime: 2019-02-26 21:49:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (402->59), mult. (449->88), div. (0->0), fcn. (361->10), ass. (0->55)
t237 = pkin(10) + qJ(4);
t234 = cos(t237);
t229 = pkin(4) * t234 + cos(pkin(10)) * pkin(3) + pkin(2);
t239 = sin(qJ(2));
t241 = cos(qJ(2));
t233 = sin(t237);
t273 = pkin(4) * qJD(4);
t261 = t233 * t273;
t262 = t239 * qJD(3);
t279 = pkin(4) * t233;
t274 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3);
t283 = t274 * t241;
t287 = (-t229 * t239 + t283) * qJD(2) + (pkin(7) + t279 + sin(pkin(10)) * pkin(3)) * qJD(1) - t241 * t261 + t262;
t235 = qJ(5) + t237;
t231 = sin(t235);
t277 = r_i_i_C(2) * t231;
t232 = cos(t235);
t278 = r_i_i_C(1) * t232;
t248 = t229 - t277 + t278;
t245 = -t248 * t239 + t283;
t242 = cos(qJ(1));
t238 = qJD(4) + qJD(5);
t254 = t238 * t241 - qJD(1);
t285 = t242 * t254;
t276 = r_i_i_C(2) * t232;
t251 = r_i_i_C(1) * t231 + t276;
t282 = t251 * t238 + t261;
t267 = qJD(1) * t241;
t253 = -t238 + t267;
t240 = sin(qJ(1));
t265 = qJD(2) * t239;
t259 = t240 * t265;
t281 = t253 * t242 - t259;
t271 = t238 * t239;
t263 = qJD(2) * t242;
t258 = t239 * t263;
t246 = t253 * t240 + t258;
t224 = t246 * t231 - t232 * t285;
t225 = t231 * t285 + t246 * t232;
t270 = t224 * r_i_i_C(1) + t225 * r_i_i_C(2);
t250 = t254 * t240;
t226 = t281 * t231 + t232 * t250;
t227 = t231 * t250 - t281 * t232;
t269 = -t226 * r_i_i_C(1) + t227 * r_i_i_C(2);
t268 = qJD(1) * t240;
t266 = qJD(1) * t242;
t264 = qJD(2) * t241;
t260 = t234 * t273;
t257 = t274 * t239;
t252 = -qJD(4) + t267;
t249 = t234 * (-qJD(4) * t241 + qJD(1));
t244 = t260 + (-t229 * t241 - pkin(1) - t257) * qJD(1);
t243 = qJD(3) * t241 + t282 * t239 + (-t248 * t241 - t257) * qJD(2);
t228 = t271 * t277;
t1 = [t227 * r_i_i_C(1) + t226 * r_i_i_C(2) - t287 * t240 + t244 * t242, t243 * t242 - t245 * t268, -t239 * t268 + t241 * t263 (t242 * t249 + (t252 * t240 + t258) * t233) * pkin(4) + t270, t270, 0; -t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t244 * t240 + t287 * t242, t243 * t240 + t245 * t266, t239 * t266 + t240 * t264 (t240 * t249 + (-t252 * t242 + t259) * t233) * pkin(4) + t269, t269, 0; 0, t245 * qJD(2) - t282 * t241 + t262, t265, t228 + (-t238 * t278 - t260) * t239 + (-t251 - t279) * t264, -t264 * t276 + t228 + (-t231 * t264 - t232 * t271) * r_i_i_C(1), 0;];
JaD_transl  = t1;
