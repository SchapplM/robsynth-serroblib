% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:12:31
% EndTime: 2019-02-26 21:12:31
% DurationCPUTime: 0.30s
% Computational Cost: add. (382->65), mult. (530->94), div. (0->0), fcn. (416->8), ass. (0->56)
t240 = qJ(4) + qJ(5);
t235 = sin(t240);
t236 = cos(t240);
t246 = cos(qJ(1));
t239 = qJD(4) + qJD(5);
t242 = sin(qJ(3));
t258 = t239 * t242 + qJD(1);
t253 = t258 * t246;
t243 = sin(qJ(1));
t245 = cos(qJ(3));
t275 = qJD(1) * t242;
t257 = t239 + t275;
t270 = qJD(3) * t246;
t286 = t257 * t243 - t245 * t270;
t222 = t286 * t235 - t236 * t253;
t244 = cos(qJ(4));
t279 = t236 * t239;
t280 = pkin(4) * qJD(4);
t227 = pkin(5) * t279 + t244 * t280;
t231 = t244 * pkin(4) + pkin(5) * t236;
t229 = pkin(3) + t231;
t281 = r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
t287 = t281 * t245;
t297 = (-t229 * t242 - qJ(2) + t287) * qJD(1) - t227;
t241 = sin(qJ(4));
t284 = pkin(5) * t235;
t226 = -t239 * t284 - t241 * t280;
t282 = r_i_i_C(2) * t236;
t251 = (r_i_i_C(1) * t235 + t282) * t239 - t226;
t283 = r_i_i_C(2) * t235;
t252 = r_i_i_C(1) * t236 + t229 - t283;
t296 = (-t252 * t242 + t287) * qJD(3) + qJD(6) * t242 - t251 * t245;
t261 = t281 * t242;
t295 = t252 * t245 + t261;
t271 = qJD(3) * t245;
t285 = t243 * t271 + t257 * t246;
t278 = t239 * t245;
t223 = -t235 * t253 - t236 * t286;
t277 = -t222 * r_i_i_C(1) + t223 * r_i_i_C(2);
t254 = t258 * t243;
t224 = -t285 * t235 - t236 * t254;
t225 = -t235 * t254 + t285 * t236;
t276 = t224 * r_i_i_C(1) - t225 * r_i_i_C(2);
t274 = qJD(1) * t243;
t273 = qJD(1) * t246;
t272 = qJD(3) * t242;
t268 = t245 * qJD(6);
t267 = t236 * t278;
t265 = t235 * t272;
t266 = r_i_i_C(1) * t265 + t272 * t282 + t278 * t283;
t230 = t241 * pkin(4) + t284;
t255 = -t230 * t275 + t226;
t250 = qJD(1) * t231 + t227 * t242 + t230 * t271;
t248 = qJD(1) * t295;
t247 = -t268 + t242 * t226 + qJD(2) + (t229 * t245 + t261) * qJD(3) + (-pkin(1) - pkin(7) - t230) * qJD(1);
t1 = [t223 * r_i_i_C(1) + t222 * r_i_i_C(2) + t297 * t243 + t247 * t246, t273, t296 * t243 + t246 * t248, -t250 * t243 + t255 * t246 + t276, t224 * pkin(5) + t276, t243 * t272 - t245 * t273; t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t247 * t243 - t297 * t246, t274, t243 * t248 - t296 * t246, t255 * t243 + t250 * t246 + t277, -t222 * pkin(5) + t277, -t242 * t270 - t245 * t274; 0, 0, -t295 * qJD(3) + t251 * t242 + t268, t230 * t272 + (-r_i_i_C(1) * t279 - t227) * t245 + t266, -r_i_i_C(1) * t267 + (t265 - t267) * pkin(5) + t266, t271;];
JaD_transl  = t1;
