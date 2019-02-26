% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR12_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR12_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:14
% EndTime: 2019-02-26 21:44:14
% DurationCPUTime: 0.21s
% Computational Cost: add. (298->82), mult. (727->130), div. (0->0), fcn. (682->10), ass. (0->50)
t261 = cos(qJ(4));
t291 = t261 * pkin(4);
t290 = pkin(8) + pkin(3) + t291;
t256 = cos(pkin(6));
t259 = sin(qJ(2));
t263 = cos(qJ(1));
t282 = t263 * t259;
t260 = sin(qJ(1));
t262 = cos(qJ(2));
t283 = t260 * t262;
t242 = t256 * t282 + t283;
t243 = t256 * t283 + t282;
t239 = t243 * qJD(1) + t242 * qJD(2);
t254 = qJ(4) + pkin(11);
t252 = sin(t254);
t289 = t239 * t252;
t253 = cos(t254);
t288 = t239 * t253;
t255 = sin(pkin(6));
t287 = t255 * t260;
t286 = t255 * t262;
t285 = t255 * t263;
t284 = t260 * t259;
t281 = t263 * t262;
t280 = qJD(1) * t260;
t279 = qJD(1) * t263;
t278 = qJD(2) * t259;
t277 = qJD(2) * t262;
t273 = t256 * t281;
t241 = -t273 + t284;
t276 = qJD(4) * t241;
t275 = -r_i_i_C(3) - qJ(5) - pkin(9) - pkin(2);
t274 = t256 * t284;
t272 = t255 * t280;
t271 = t255 * t279;
t270 = t255 * t278;
t269 = qJD(2) * t256 + qJD(1);
t268 = r_i_i_C(1) * t253 - r_i_i_C(2) * t252;
t258 = sin(qJ(4));
t267 = -t258 * pkin(4) - t252 * r_i_i_C(1) - t253 * r_i_i_C(2);
t266 = t268 + t291;
t265 = qJ(3) - t267;
t264 = t266 * qJD(4) + qJD(3);
t244 = -t274 + t281;
t240 = -qJD(1) * t274 - t260 * t278 + t269 * t281;
t238 = t242 * qJD(1) + t243 * qJD(2);
t237 = -qJD(1) * t273 - t263 * t277 + t269 * t284;
t236 = t253 * t271 - t237 * t252 + (t243 * t253 - t252 * t287) * qJD(4);
t235 = -t252 * t271 - t237 * t253 + (-t243 * t252 - t253 * t287) * qJD(4);
t1 = [(-t253 * t276 - t289) * r_i_i_C(1) + (t252 * t276 - t288) * r_i_i_C(2) - t242 * qJD(5) - t239 * qJ(3) - t241 * qJD(3) - pkin(1) * t279 + (-t239 * t258 - t261 * t276) * pkin(4) + t275 * t240 + (t267 * t263 * qJD(4) + (-t268 - t290) * t280) * t255, -t243 * qJD(5) - t275 * t237 - t265 * t238 + t264 * t244, -t237, t235 * r_i_i_C(1) - t236 * r_i_i_C(2) + (-t258 * t271 - t237 * t261 + (-t243 * t258 - t261 * t287) * qJD(4)) * pkin(4), -t238, 0; t236 * r_i_i_C(1) + t235 * r_i_i_C(2) - t237 * qJ(3) + t243 * qJD(3) + t244 * qJD(5) + t275 * t238 + (-pkin(1) * t260 + t290 * t285) * qJD(1) + (-t237 * t258 + (t243 * t261 - t258 * t287) * qJD(4)) * pkin(4), -t241 * qJD(5) + t275 * t239 + t265 * t240 + t264 * t242, t239 (-t252 * t272 + t288) * r_i_i_C(1) + (-t253 * t272 - t289) * r_i_i_C(2) + ((-t241 * t252 + t253 * t285) * r_i_i_C(1) + (-t241 * t253 - t252 * t285) * r_i_i_C(2)) * qJD(4) + (-t258 * t272 + t239 * t261 + (-t241 * t258 + t261 * t285) * qJD(4)) * pkin(4), t240, 0; 0 (qJD(5) * t262 + t264 * t259 + (t275 * t259 + t265 * t262) * qJD(2)) * t255, t270, t266 * t270 + ((t252 * t286 - t253 * t256) * r_i_i_C(1) + (t252 * t256 + t253 * t286) * r_i_i_C(2) + (-t256 * t261 + t258 * t286) * pkin(4)) * qJD(4), t255 * t277, 0;];
JaD_transl  = t1;
