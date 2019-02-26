% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_transl = S6RRPRRR11_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR11_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:59
% EndTime: 2019-02-26 21:59:59
% DurationCPUTime: 0.25s
% Computational Cost: add. (279->55), mult. (490->85), div. (0->0), fcn. (384->8), ass. (0->48)
t236 = sin(qJ(2));
t239 = cos(qJ(2));
t238 = cos(qJ(4));
t271 = t238 * pkin(4);
t250 = qJD(4) * t271 + qJD(3);
t259 = pkin(2) + r_i_i_C(3) + pkin(9) + pkin(8);
t252 = t259 * t236;
t235 = sin(qJ(4));
t256 = pkin(4) * t235 + qJ(3);
t283 = qJD(2) * (-t256 * t239 + t252) - (pkin(7) + pkin(3) + t271) * qJD(1) - t250 * t236;
t234 = qJ(4) + qJ(5);
t231 = sin(t234);
t232 = cos(t234);
t281 = r_i_i_C(1) * t231 + r_i_i_C(2) * t232;
t280 = r_i_i_C(1) * t232 - r_i_i_C(2) * t231;
t237 = sin(qJ(1));
t233 = qJD(4) + qJD(5);
t255 = t233 * t236 + qJD(1);
t279 = t237 * t255;
t240 = cos(qJ(1));
t278 = t240 * t255;
t266 = qJD(1) * t236;
t254 = -t233 - t266;
t262 = qJD(2) * t239;
t258 = t237 * t262;
t246 = t254 * t240 - t258;
t223 = t231 * t279 + t246 * t232;
t224 = t246 * t231 - t232 * t279;
t268 = -t223 * r_i_i_C(1) + t224 * r_i_i_C(2);
t261 = qJD(2) * t240;
t257 = t239 * t261;
t245 = t254 * t237 + t257;
t225 = -t231 * t278 + t245 * t232;
t226 = t245 * t231 + t232 * t278;
t267 = t225 * r_i_i_C(1) - t226 * r_i_i_C(2);
t265 = qJD(1) * t237;
t264 = qJD(1) * t240;
t263 = qJD(2) * t236;
t260 = qJD(4) * t235;
t253 = qJD(4) + t266;
t251 = t259 * t239;
t249 = (-qJD(4) * t236 - qJD(1)) * t235;
t248 = t233 * t239 * t281 + t263 * t280;
t247 = t256 + t281;
t244 = t233 * t280 + t250;
t243 = t247 * t239 - t252;
t242 = -pkin(4) * t260 + (-t256 * t236 - pkin(1) - t251) * qJD(1);
t1 = [t224 * r_i_i_C(1) + t223 * r_i_i_C(2) + t283 * t237 + t242 * t240 (-t247 * t261 + t259 * t265) * t236 + (-t247 * t265 + (-t259 * qJD(2) + t244) * t240) * t239, -t236 * t265 + t257 (t240 * t249 + (-t253 * t237 + t257) * t238) * pkin(4) + t267, t267, 0; t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t242 * t237 - t283 * t240, t243 * t264 + (t244 * t239 + (-t247 * t236 - t251) * qJD(2)) * t237, t236 * t264 + t258 (t253 * t240 * t238 + (t238 * t262 + t249) * t237) * pkin(4) + t268, t268, 0; 0, t243 * qJD(2) + t244 * t236, t263 (t238 * t263 + t239 * t260) * pkin(4) + t248, t248, 0;];
JaD_transl  = t1;
