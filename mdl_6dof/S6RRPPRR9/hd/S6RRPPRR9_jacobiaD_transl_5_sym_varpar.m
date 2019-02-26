% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR9
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
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR9_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR9_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:14
% EndTime: 2019-02-26 21:33:14
% DurationCPUTime: 0.23s
% Computational Cost: add. (213->61), mult. (627->100), div. (0->0), fcn. (592->8), ass. (0->43)
t248 = sin(qJ(5));
t251 = cos(qJ(5));
t257 = r_i_i_C(1) * t251 - r_i_i_C(2) * t248;
t278 = t257 * qJD(5) + qJD(4);
t277 = -pkin(2) - qJ(4);
t246 = sin(pkin(6));
t276 = t246 * t248;
t275 = t246 * t251;
t253 = cos(qJ(1));
t274 = t246 * t253;
t249 = sin(qJ(2));
t250 = sin(qJ(1));
t273 = t249 * t250;
t272 = t249 * t253;
t252 = cos(qJ(2));
t271 = t250 * t252;
t270 = t252 * t253;
t269 = qJD(1) * t250;
t268 = qJD(1) * t253;
t267 = qJD(2) * t249;
t266 = qJD(2) * t252;
t265 = pkin(3) + pkin(4) + pkin(8);
t264 = r_i_i_C(3) + pkin(9) - qJ(3);
t247 = cos(pkin(6));
t263 = t247 * t273;
t262 = t247 * t270;
t261 = t246 * t269;
t260 = t246 * t268;
t259 = t246 * t266;
t258 = qJD(2) * t247 + qJD(1);
t256 = -r_i_i_C(1) * t248 - r_i_i_C(2) * t251;
t238 = t247 * t271 + t272;
t237 = t247 * t272 + t271;
t255 = -t256 - t277;
t239 = -t263 + t270;
t236 = -t262 + t273;
t235 = -qJD(1) * t263 - t250 * t267 + t258 * t270;
t234 = t238 * qJD(1) + t237 * qJD(2);
t233 = t237 * qJD(1) + t238 * qJD(2);
t232 = -qJD(1) * t262 - t253 * t266 + t258 * t273;
t231 = t251 * t260 - t233 * t248 + (t239 * t251 - t250 * t276) * qJD(5);
t230 = -t248 * t260 - t233 * t251 + (-t239 * t248 - t250 * t275) * qJD(5);
t1 = [-pkin(1) * t268 - t236 * qJD(3) - t278 * t237 - t255 * t235 + t264 * t234 + (t256 * t253 * qJD(5) + (-t257 - t265) * t269) * t246, t239 * qJD(3) + t255 * t232 + t264 * t233 - t238 * t278, -t232, -t233, r_i_i_C(1) * t230 - t231 * r_i_i_C(2), 0; t231 * r_i_i_C(1) + t230 * r_i_i_C(2) + t238 * qJD(3) + t239 * qJD(4) + t277 * t233 + t264 * t232 + (-pkin(1) * t250 + t265 * t274) * qJD(1), qJD(3) * t237 - t255 * t234 - t264 * t235 - t236 * t278, t234, t235 (t235 * t251 - t248 * t261) * r_i_i_C(1) + (-t235 * t248 - t251 * t261) * r_i_i_C(2) + ((-t237 * t248 + t251 * t274) * r_i_i_C(1) + (-t237 * t251 - t248 * t274) * r_i_i_C(2)) * qJD(5), 0; 0 (qJD(3) * t249 + t278 * t252 + (-t255 * t249 - t264 * t252) * qJD(2)) * t246, t246 * t267, t259, t257 * t259 + ((-t247 * t251 - t249 * t276) * r_i_i_C(1) + (t247 * t248 - t249 * t275) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
