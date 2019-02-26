% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RPRRRP12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP12_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP12_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobiaD_transl_3_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:13
% EndTime: 2019-02-26 21:14:13
% DurationCPUTime: 0.16s
% Computational Cost: add. (103->49), mult. (358->92), div. (0->0), fcn. (352->10), ass. (0->47)
t248 = cos(pkin(7));
t278 = r_i_i_C(3) + pkin(9);
t279 = t278 * t248 + qJ(2);
t250 = sin(qJ(3));
t277 = t250 * r_i_i_C(2);
t246 = sin(pkin(6));
t251 = sin(qJ(1));
t276 = t246 * t251;
t253 = cos(qJ(1));
t275 = t246 * t253;
t274 = t248 * t250;
t252 = cos(qJ(3));
t273 = t248 * t252;
t244 = sin(pkin(12));
t272 = t251 * t244;
t247 = cos(pkin(12));
t271 = t251 * t247;
t270 = t253 * t244;
t269 = t253 * t247;
t268 = qJD(1) * t251;
t267 = qJD(1) * t253;
t245 = sin(pkin(7));
t266 = qJD(3) * t245;
t265 = t278 * t245;
t249 = cos(pkin(6));
t264 = t249 * t272;
t263 = t246 * t268;
t262 = t246 * t267;
t261 = t266 * t275;
t260 = r_i_i_C(1) * t250 + r_i_i_C(2) * t252;
t236 = -qJD(1) * t264 + t247 * t267;
t237 = -t249 * t269 + t272;
t259 = qJD(3) * t237 * t248 - t236;
t258 = t260 * t245;
t256 = t249 * t271 + t270;
t257 = t245 * t276 - t256 * t248;
t238 = t249 * t270 + t271;
t233 = t237 * qJD(1);
t255 = t233 * t248 + t245 * t262;
t235 = t256 * qJD(1);
t254 = -qJD(3) * t238 - t235 * t248 + t245 * t263;
t241 = t252 * t261;
t240 = -t264 + t269;
t234 = t238 * qJD(1);
t232 = -t234 * t252 + t255 * t250 + (-t240 * t250 + t257 * t252) * qJD(3);
t231 = t234 * t250 + t255 * t252 + (-t240 * t252 - t257 * t250) * qJD(3);
t1 = [-pkin(1) * t267 + t241 * r_i_i_C(1) + (-t252 * r_i_i_C(1) - pkin(2) + t277) * t236 + (t260 * t248 - t265) * t235 + ((t237 * t273 + t238 * t250) * r_i_i_C(1) + (-t237 * t274 + t238 * t252) * r_i_i_C(2)) * qJD(3) + ((-t266 * t277 + qJD(2)) * t253 + (-t258 - t279) * t268) * t246, t262, t231 * r_i_i_C(1) - t232 * r_i_i_C(2), 0, 0, 0; qJD(2) * t276 - t234 * pkin(2) + t232 * r_i_i_C(1) + t231 * r_i_i_C(2) - t233 * t265 + (-t251 * pkin(1) + t279 * t275) * qJD(1), t263, t241 * r_i_i_C(2) + (t254 * r_i_i_C(1) + t259 * r_i_i_C(2)) * t252 + ((t259 + t261) * r_i_i_C(1) - t254 * r_i_i_C(2)) * t250, 0, 0, 0; 0, 0 (-t249 * t258 + ((-t244 * t252 - t247 * t274) * r_i_i_C(1) + (t244 * t250 - t247 * t273) * r_i_i_C(2)) * t246) * qJD(3), 0, 0, 0;];
JaD_transl  = t1;
