% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:51:39
% EndTime: 2019-02-26 21:51:39
% DurationCPUTime: 0.32s
% Computational Cost: add. (396->66), mult. (592->93), div. (0->0), fcn. (462->8), ass. (0->56)
t241 = qJ(4) + qJ(5);
t236 = sin(t241);
t237 = cos(t241);
t244 = sin(qJ(1));
t240 = qJD(4) + qJD(5);
t243 = sin(qJ(2));
t262 = t240 * t243 + qJD(1);
t256 = t262 * t244;
t247 = cos(qJ(1));
t276 = qJD(1) * t243;
t261 = t240 + t276;
t246 = cos(qJ(2));
t272 = qJD(2) * t246;
t264 = t244 * t272;
t287 = t261 * t247 + t264;
t297 = -t236 * t256 + t287 * t237;
t245 = cos(qJ(4));
t233 = t245 * pkin(4) + pkin(5) * t237;
t281 = pkin(4) * qJD(4);
t285 = pkin(5) * t240;
t228 = t237 * t285 + t245 * t281;
t269 = qJD(3) + t228;
t270 = t246 * qJD(6);
t242 = sin(qJ(4));
t232 = t242 * pkin(4) + pkin(5) * t236;
t279 = qJ(3) + t232;
t268 = pkin(2) + r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
t289 = t268 * t243;
t296 = (-t279 * t246 + t289) * qJD(2) - (pkin(7) + pkin(3) + t233) * qJD(1) - t269 * t243 - t270;
t273 = qJD(2) * t243;
t280 = t240 * t246;
t295 = t236 * t280 + t237 * t273;
t283 = r_i_i_C(2) * t237;
t255 = r_i_i_C(1) * t236 + t279 + t283;
t293 = -t255 * t246 + t289;
t291 = t247 * t262;
t284 = r_i_i_C(2) * t236;
t224 = -t236 * t287 - t237 * t256;
t278 = r_i_i_C(1) * t297 + t224 * r_i_i_C(2);
t271 = qJD(2) * t247;
t263 = t246 * t271;
t251 = -t261 * t244 + t263;
t225 = -t236 * t291 + t251 * t237;
t226 = t251 * t236 + t237 * t291;
t277 = t225 * r_i_i_C(1) - t226 * r_i_i_C(2);
t275 = qJD(1) * t244;
t274 = qJD(1) * t247;
t266 = t295 * r_i_i_C(1) + t280 * t283;
t259 = t268 * t246;
t257 = t233 * t276 + t228;
t254 = (r_i_i_C(1) * t237 - t284) * t240 + t269;
t227 = -t236 * t285 - t242 * t281;
t253 = -qJD(1) * t232 + t227 * t243 + t233 * t272;
t250 = t227 + (-t279 * t243 - pkin(1) - t259) * qJD(1);
t248 = -qJD(6) * t243 + t254 * t246 + (-t255 * t243 - t259) * qJD(2);
t1 = [t224 * r_i_i_C(1) - r_i_i_C(2) * t297 + t296 * t244 + t250 * t247, t248 * t247 + t293 * t275, -t243 * t275 + t263, -t257 * t244 + t253 * t247 + t277, t225 * pkin(5) + t277, -t243 * t271 - t246 * t275; t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t250 * t244 - t296 * t247, t248 * t244 - t274 * t293, t243 * t274 + t264, t253 * t244 + t257 * t247 + t278, t297 * pkin(5) + t278, -t244 * t273 + t246 * t274; 0, -qJD(2) * t293 + t254 * t243 + t270, t273, -t246 * t227 + (t233 - t284) * t273 + t266, t295 * pkin(5) - t273 * t284 + t266, t272;];
JaD_transl  = t1;
