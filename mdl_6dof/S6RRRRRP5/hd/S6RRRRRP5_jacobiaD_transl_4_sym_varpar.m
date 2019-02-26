% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP5_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP5_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:54
% EndTime: 2019-02-26 22:41:54
% DurationCPUTime: 0.26s
% Computational Cost: add. (261->56), mult. (420->90), div. (0->0), fcn. (332->8), ass. (0->52)
t235 = cos(qJ(3));
t227 = pkin(3) * t235 + pkin(2);
t233 = sin(qJ(2));
t236 = cos(qJ(2));
t267 = r_i_i_C(3) + pkin(9) + pkin(8);
t251 = t267 * t236;
t232 = sin(qJ(3));
t266 = pkin(3) * qJD(3);
t256 = t232 * t266;
t276 = (-t227 * t233 + t251) * qJD(2) - t236 * t256;
t237 = cos(qJ(1));
t230 = qJD(3) + qJD(4);
t250 = t230 * t236 - qJD(1);
t274 = t237 * t250;
t231 = qJ(3) + qJ(4);
t228 = sin(t231);
t229 = cos(t231);
t268 = r_i_i_C(2) * t229;
t246 = r_i_i_C(1) * t228 + t268;
t273 = -t230 * t246 - t256;
t260 = qJD(1) * t236;
t249 = -t230 + t260;
t234 = sin(qJ(1));
t258 = qJD(2) * t233;
t253 = t234 * t258;
t272 = t237 * t249 - t253;
t271 = pkin(3) * t232;
t270 = r_i_i_C(1) * t229;
t269 = r_i_i_C(2) * t228;
t264 = t230 * t233;
t252 = t237 * t258;
t240 = t234 * t249 + t252;
t222 = t228 * t240 - t229 * t274;
t223 = t228 * t274 + t229 * t240;
t263 = r_i_i_C(1) * t222 + r_i_i_C(2) * t223;
t245 = t250 * t234;
t224 = t228 * t272 + t229 * t245;
t225 = t228 * t245 - t229 * t272;
t262 = -r_i_i_C(1) * t224 + r_i_i_C(2) * t225;
t261 = qJD(1) * t234;
t259 = qJD(1) * t237;
t257 = qJD(2) * t236;
t255 = t235 * t266;
t254 = pkin(7) + t271;
t247 = -qJD(3) + t260;
t244 = (-qJD(3) * t236 + qJD(1)) * t235;
t243 = t227 - t269 + t270;
t242 = -t227 * t236 - t233 * t267 - pkin(1);
t241 = qJD(2) * t243;
t239 = -qJD(2) * t267 - t273;
t226 = t264 * t269;
t1 = [t237 * t255 + t225 * r_i_i_C(1) + t224 * r_i_i_C(2) - t276 * t234 + (-t234 * t254 + t237 * t242) * qJD(1) (-t237 * t241 - t261 * t267) * t236 + (t237 * t239 + t243 * t261) * t233 (t237 * t244 + (t234 * t247 + t252) * t232) * pkin(3) + t263, t263, 0, 0; t234 * t255 - t223 * r_i_i_C(1) + t222 * r_i_i_C(2) + t276 * t237 + (t234 * t242 + t237 * t254) * qJD(1) (-t234 * t241 + t259 * t267) * t236 + (t234 * t239 - t243 * t259) * t233 (t234 * t244 + (-t237 * t247 + t253) * t232) * pkin(3) + t262, t262, 0, 0; 0, t273 * t236 + (-t233 * t243 + t251) * qJD(2), t226 + (-t230 * t270 - t255) * t233 + (-t246 - t271) * t257, -t257 * t268 + t226 + (-t228 * t257 - t229 * t264) * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
