% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:15
% EndTime: 2019-02-26 20:50:15
% DurationCPUTime: 0.28s
% Computational Cost: add. (530->64), mult. (453->93), div. (0->0), fcn. (363->12), ass. (0->56)
t244 = sin(qJ(3));
t241 = pkin(11) + qJ(5);
t237 = cos(t241);
t231 = pkin(5) * t237 + cos(pkin(11)) * pkin(4) + pkin(3);
t239 = qJ(6) + t241;
t233 = sin(t239);
t278 = r_i_i_C(2) * t233;
t234 = cos(t239);
t279 = r_i_i_C(1) * t234;
t251 = t231 - t278 + t279;
t245 = cos(qJ(3));
t275 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
t283 = t275 * t245;
t247 = -t251 * t244 + t283;
t288 = qJD(1) * t247;
t235 = sin(t241);
t274 = pkin(5) * qJD(5);
t264 = t235 * t274;
t265 = t244 * qJD(4);
t287 = (-t231 * t244 + t283) * qJD(3) - t245 * t264 + t265;
t243 = qJ(1) + pkin(10);
t238 = cos(t243);
t242 = qJD(5) + qJD(6);
t257 = t242 * t245 - qJD(1);
t285 = t238 * t257;
t268 = qJD(1) * t245;
t256 = -t242 + t268;
t236 = sin(t243);
t267 = qJD(3) * t244;
t262 = t236 * t267;
t282 = t256 * t238 - t262;
t277 = r_i_i_C(2) * t234;
t254 = r_i_i_C(1) * t233 + t277;
t281 = t254 * t242 + t264;
t280 = pkin(5) * t235;
t276 = pkin(7) + t280 + sin(pkin(11)) * pkin(4);
t272 = t242 * t244;
t261 = t238 * t267;
t248 = t256 * t236 + t261;
t226 = t248 * t233 - t234 * t285;
t227 = t233 * t285 + t248 * t234;
t271 = t226 * r_i_i_C(1) + t227 * r_i_i_C(2);
t253 = t257 * t236;
t228 = t282 * t233 + t234 * t253;
t229 = t233 * t253 - t282 * t234;
t270 = -t228 * r_i_i_C(1) + t229 * r_i_i_C(2);
t269 = qJD(1) * t244;
t266 = qJD(3) * t245;
t263 = t237 * t274;
t260 = t275 * t244;
t255 = -qJD(5) + t268;
t252 = (-qJD(5) * t245 + qJD(1)) * t237;
t250 = -t231 * t245 - pkin(2) - t260;
t246 = qJD(4) * t245 + t281 * t244 + (-t251 * t245 - t260) * qJD(3);
t230 = t272 * t278;
t1 = [t238 * t263 + t229 * r_i_i_C(1) + t228 * r_i_i_C(2) - t287 * t236 + (-cos(qJ(1)) * pkin(1) - t276 * t236 + t250 * t238) * qJD(1), 0, -t236 * t288 + t246 * t238, -t236 * t269 + t238 * t266 (t238 * t252 + (t255 * t236 + t261) * t235) * pkin(5) + t271, t271; t236 * t263 - t227 * r_i_i_C(1) + t226 * r_i_i_C(2) + t287 * t238 + (-sin(qJ(1)) * pkin(1) + t276 * t238 + t250 * t236) * qJD(1), 0, t246 * t236 + t238 * t288, t236 * t266 + t238 * t269 (t236 * t252 + (-t255 * t238 + t262) * t235) * pkin(5) + t270, t270; 0, 0, t247 * qJD(3) - t281 * t245 + t265, t267, t230 + (-t242 * t279 - t263) * t244 + (-t254 - t280) * t266, -t266 * t277 + t230 + (-t233 * t266 - t234 * t272) * r_i_i_C(1);];
JaD_transl  = t1;
