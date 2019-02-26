% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:49
% EndTime: 2019-02-26 21:30:49
% DurationCPUTime: 0.21s
% Computational Cost: add. (217->64), mult. (636->104), div. (0->0), fcn. (594->8), ass. (0->45)
t244 = sin(qJ(5));
t247 = cos(qJ(5));
t275 = r_i_i_C(2) * t247;
t253 = r_i_i_C(1) * t244 + t275;
t251 = qJD(5) * t253;
t274 = pkin(8) - qJ(4);
t242 = sin(pkin(6));
t246 = sin(qJ(1));
t273 = t242 * t246;
t272 = t242 * t247;
t249 = cos(qJ(1));
t271 = t242 * t249;
t245 = sin(qJ(2));
t270 = t245 * t246;
t269 = t245 * t249;
t248 = cos(qJ(2));
t268 = t246 * t248;
t267 = t248 * t249;
t266 = qJD(1) * t246;
t265 = qJD(1) * t249;
t264 = qJD(2) * t245;
t263 = qJD(2) * t248;
t243 = cos(pkin(6));
t233 = t243 * t269 + t268;
t262 = qJD(5) * t233;
t261 = -pkin(2) - pkin(3) - pkin(4);
t260 = r_i_i_C(3) + pkin(9) - qJ(3);
t259 = t243 * t270;
t258 = t243 * t267;
t257 = t242 * t266;
t256 = t242 * t265;
t255 = qJD(2) * t243 + qJD(1);
t254 = -r_i_i_C(1) * t247 + r_i_i_C(2) * t244;
t252 = t243 * t268 + t269;
t250 = -t254 - t261;
t236 = t244 * t257;
t235 = -t259 + t267;
t232 = -t258 + t270;
t231 = -qJD(1) * t259 - t246 * t264 + t255 * t267;
t230 = t252 * qJD(1) + t233 * qJD(2);
t229 = t233 * qJD(1) + t252 * qJD(2);
t228 = -qJD(1) * t258 - t249 * t263 + t255 * t270;
t227 = -t244 * t256 - t229 * t247 + (-t235 * t244 - t246 * t272) * qJD(5);
t226 = -t247 * t256 + t229 * t244 + (-t235 * t247 + t244 * t273) * qJD(5);
t1 = [(t244 * t262 + t236) * r_i_i_C(1) + t262 * t275 - t232 * qJD(3) - pkin(1) * t265 - t250 * t231 + t260 * t230 + ((t254 * qJD(5) - qJD(4)) * t249 + (-t274 + t275) * t266) * t242, qJD(3) * t235 + t250 * t228 + t260 * t229 + t251 * t252, -t228, -t256, r_i_i_C(1) * t226 - t227 * r_i_i_C(2), 0; -qJD(4) * t273 + t227 * r_i_i_C(1) + t226 * r_i_i_C(2) + t252 * qJD(3) + t261 * t229 + t260 * t228 + (-pkin(1) * t246 + t274 * t271) * qJD(1), qJD(3) * t233 - t250 * t230 - t260 * t231 + t232 * t251, t230, -t257 (-t231 * t244 - t247 * t257) * r_i_i_C(1) + (-t231 * t247 + t236) * r_i_i_C(2) + ((-t233 * t247 - t244 * t271) * r_i_i_C(1) + (t233 * t244 - t247 * t271) * r_i_i_C(2)) * qJD(5), 0; 0 (qJD(3) * t245 - t248 * t251 + (-t250 * t245 - t260 * t248) * qJD(2)) * t242, t242 * t264, 0, -t253 * t242 * t263 + ((t243 * t244 - t245 * t272) * r_i_i_C(1) + (t242 * t244 * t245 + t243 * t247) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
