% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:49
% EndTime: 2019-02-26 21:59:49
% DurationCPUTime: 0.32s
% Computational Cost: add. (596->72), mult. (642->100), div. (0->0), fcn. (504->10), ass. (0->61)
t249 = qJ(4) + qJ(5);
t244 = cos(t249);
t239 = pkin(5) * t244;
t253 = cos(qJ(4));
t238 = t253 * pkin(4) + t239;
t251 = sin(qJ(2));
t254 = cos(qJ(2));
t274 = pkin(2) + r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t267 = t274 * t251;
t247 = qJD(4) + qJD(5);
t287 = pkin(4) * qJD(4);
t232 = t247 * t239 + t253 * t287;
t275 = qJD(3) + t232;
t243 = sin(t249);
t250 = sin(qJ(4));
t237 = t250 * pkin(4) + pkin(5) * t243;
t284 = qJ(3) + t237;
t299 = (-t284 * t254 + t267) * qJD(2) - (pkin(7) + pkin(3) + t238) * qJD(1) - t275 * t251;
t245 = qJ(6) + t249;
t240 = sin(t245);
t241 = cos(t245);
t298 = r_i_i_C(1) * t240 + r_i_i_C(2) * t241;
t252 = sin(qJ(1));
t242 = qJD(6) + t247;
t270 = t242 * t251 + qJD(1);
t296 = t252 * t270;
t255 = cos(qJ(1));
t295 = t255 * t270;
t291 = r_i_i_C(1) * t241;
t290 = r_i_i_C(2) * t240;
t285 = t243 * t247;
t281 = qJD(1) * t251;
t269 = -t242 - t281;
t277 = qJD(2) * t254;
t272 = t252 * t277;
t259 = t255 * t269 - t272;
t227 = t240 * t296 + t241 * t259;
t228 = t240 * t259 - t241 * t296;
t283 = -t227 * r_i_i_C(1) + t228 * r_i_i_C(2);
t276 = qJD(2) * t255;
t271 = t254 * t276;
t258 = t252 * t269 + t271;
t229 = -t240 * t295 + t241 * t258;
t230 = t240 * t258 + t241 * t295;
t282 = t229 * r_i_i_C(1) - t230 * r_i_i_C(2);
t280 = qJD(1) * t252;
t279 = qJD(1) * t255;
t278 = qJD(2) * t251;
t273 = t298 * t242 * t254 + t278 * t291;
t268 = t247 + t281;
t266 = t274 * t254;
t265 = t238 * t281 + t232;
t264 = t243 * (-t247 * t251 - qJD(1));
t263 = t284 + t298;
t262 = -t278 * t290 + t273;
t261 = (-t290 + t291) * t242 + t275;
t231 = -pkin(5) * t285 - t250 * t287;
t260 = -qJD(1) * t237 + t231 * t251 + t238 * t277;
t257 = t231 + (-t284 * t251 - pkin(1) - t266) * qJD(1);
t256 = t254 * t263 - t267;
t1 = [t228 * r_i_i_C(1) + t227 * r_i_i_C(2) + t299 * t252 + t257 * t255 (-t263 * t276 + t274 * t280) * t251 + (-t263 * t280 + (-qJD(2) * t274 + t261) * t255) * t254, -t251 * t280 + t271, -t252 * t265 + t255 * t260 + t282 (t255 * t264 + (-t252 * t268 + t271) * t244) * pkin(5) + t282, t282; t230 * r_i_i_C(1) + t229 * r_i_i_C(2) + t257 * t252 - t299 * t255, t256 * t279 + (t261 * t254 + (-t251 * t263 - t266) * qJD(2)) * t252, t251 * t279 + t272, t252 * t260 + t255 * t265 + t283 (t252 * t264 + (t255 * t268 + t272) * t244) * pkin(5) + t283, t283; 0, qJD(2) * t256 + t251 * t261, t278, -t254 * t231 + (t238 - t290) * t278 + t273 (t244 * t278 + t254 * t285) * pkin(5) + t262, t262;];
JaD_transl  = t1;
